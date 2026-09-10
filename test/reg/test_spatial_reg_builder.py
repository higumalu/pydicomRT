"""
Spatial Registration builder.

A registration object that is merely *readable* is not enough -- a matrix in the wrong
direction, at the wrong nesting level, or missing its identity item still loads without
complaint and quietly misplaces an image. The structure assertions here mirror what a real
TPS emits; see test/real_data/test_clinical_reg_builder.py for the comparison against an
actual one.
"""

import numpy as np
import pytest
from pydicom import dcmread
from pydicom.uid import ImplicitVRLittleEndian

from synthetic import make_ct_series
from pydicomrt.reg import (
    SpatialRegistrationBuilder,
    check_spatial_reg_iod,
    get_spatial_registrations,
)
from pydicomrt.reg.builder import DEFAULT_UID_PREFIX, validate_transformation_matrix

SPATIAL_REGISTRATION_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.66.1"


def translation_matrix(offset):
    matrix = np.identity(4)
    matrix[:3, 3] = offset
    return matrix.ravel().tolist()


@pytest.fixture
def fixed_series():
    return make_ct_series(shape=(4, 8, 8))


@pytest.fixture
def moving_series():
    return make_ct_series(shape=(4, 8, 8))


@pytest.fixture
def reg_ds(fixed_series, moving_series):
    builder = SpatialRegistrationBuilder(fixed_series)
    builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    builder.add_registration(moving_series, translation_matrix([5.0, -3.0, 2.0]))
    return builder.build()


class TestIodConformance:
    def test_passes_the_iod_check(self, reg_ds):
        result = check_spatial_reg_iod(reg_ds)

        assert result["result"] is True, result["content"]

    def test_carries_the_spatial_registration_sop_class(self, reg_ds):
        assert reg_ds.SOPClassUID == SPATIAL_REGISTRATION_SOP_CLASS_UID
        assert reg_ds.file_meta.MediaStorageSOPClassUID == SPATIAL_REGISTRATION_SOP_CLASS_UID
        assert reg_ds.Modality == "REG"

    def test_content_identification_is_present(self, reg_ds):
        """Instance Number, Content Date and Content Time are Type 1 and were all missing."""
        assert reg_ds.InstanceNumber == 1
        assert reg_ds.ContentDate
        assert reg_ds.ContentTime
        assert reg_ds.ContentLabel
        # ContentLabel is VR CS: 16 characters maximum
        assert len(reg_ds.ContentLabel) <= 16

    def test_implementation_class_uid_is_a_real_uid(self, fixed_series, moving_series):
        """Without a prefix this used to be built as the string "1"."""
        builder = SpatialRegistrationBuilder(fixed_series)
        builder.add_registration(moving_series, translation_matrix([0.0, 0.0, 0.0]))
        ds = builder.build()

        uid = ds.file_meta.ImplementationClassUID
        assert uid.startswith(DEFAULT_UID_PREFIX)
        assert len(uid) > 1 and all(part.isdigit() for part in uid.split("."))

    def test_inherits_patient_and_study_context(self, reg_ds, fixed_series):
        assert reg_ds.PatientID == fixed_series[0].PatientID
        assert reg_ds.StudyInstanceUID == fixed_series[0].StudyInstanceUID
        assert reg_ds.FrameOfReferenceUID == fixed_series[0].FrameOfReferenceUID


class TestRegistrationSequenceStructure:
    def test_includes_an_identity_item_for_the_fixed_frame(self, reg_ds, fixed_series):
        """
        Real systems emit this, and it is what marks which Frame of Reference the other
        matrices are relative to. The builder used to write only the moving item.
        """
        assert len(reg_ds.RegistrationSequence) == 2

        by_frame = {item.FrameOfReferenceUID: item for item in reg_ds.RegistrationSequence}
        identity_item = by_frame[fixed_series[0].FrameOfReferenceUID]
        matrix = (identity_item.MatrixRegistrationSequence[0]
                  .MatrixSequence[0].FrameOfReferenceTransformationMatrix)

        assert np.allclose(np.array(matrix, dtype=float).reshape(4, 4), np.identity(4))

    def test_identity_item_can_be_suppressed(self, fixed_series, moving_series):
        builder = SpatialRegistrationBuilder(fixed_series)
        builder.add_registration(moving_series, translation_matrix([1.0, 0.0, 0.0]))

        ds = builder.build(include_identity=False)

        assert len(ds.RegistrationSequence) == 1

    def test_registration_type_code_sits_inside_matrix_registration_sequence(self, reg_ds):
        """
        PS3.3 C.20.2 puts Registration Type Code Sequence alongside Matrix Sequence, not in
        the Registration Sequence item. It used to be written one level too high, where a
        conformant reader would never look for it.
        """
        for item in reg_ds.RegistrationSequence:
            assert not hasattr(item, "RegistrationTypeCodeSequence")

            matrix_registration = item.MatrixRegistrationSequence[0]
            assert hasattr(matrix_registration, "RegistrationTypeCodeSequence")
            assert hasattr(matrix_registration, "MatrixSequence")

    def test_identity_and_moving_items_use_different_codes(self, reg_ds, fixed_series):
        codes = {}
        for item in reg_ds.RegistrationSequence:
            code = item.MatrixRegistrationSequence[0].RegistrationTypeCodeSequence[0]
            codes[item.FrameOfReferenceUID] = code.CodeValue

        assert codes[fixed_series[0].FrameOfReferenceUID] == "125021"  # FoR Identity
        moving_codes = [v for k, v in codes.items() if k != fixed_series[0].FrameOfReferenceUID]
        assert moving_codes == ["125024"]  # Image Content-based Alignment

    def test_references_every_moving_instance(self, reg_ds, moving_series):
        by_frame = {item.FrameOfReferenceUID: item for item in reg_ds.RegistrationSequence}
        moving_item = by_frame[moving_series[0].FrameOfReferenceUID]

        referenced = [ref.ReferencedSOPInstanceUID
                      for ref in moving_item.ReferencedImageSequence]
        assert referenced == [ds.SOPInstanceUID for ds in moving_series]


class TestCommonInstanceReference:
    def test_lists_both_series(self, reg_ds, fixed_series, moving_series):
        series_uids = {item.SeriesInstanceUID for item in reg_ds.ReferencedSeriesSequence}

        assert series_uids == {
            fixed_series[0].SeriesInstanceUID,
            moving_series[0].SeriesInstanceUID,
        }

    def test_lists_every_instance_of_each_series(self, reg_ds, fixed_series):
        by_series = {item.SeriesInstanceUID: item for item in reg_ds.ReferencedSeriesSequence}
        fixed_item = by_series[fixed_series[0].SeriesInstanceUID]

        assert len(fixed_item.ReferencedInstanceSequence) == len(fixed_series)


class TestRoundTrip:
    def test_survives_write_and_read(self, reg_ds, tmp_path, fixed_series, moving_series):
        path = tmp_path / "registration.dcm"
        reg_ds.save_as(path, enforce_file_format=True)

        reloaded = dcmread(path)

        assert reloaded.file_meta.TransferSyntaxUID == ImplicitVRLittleEndian
        assert check_spatial_reg_iod(reloaded)["result"] is True
        assert len(reloaded.RegistrationSequence) == 2

    def test_matrix_reads_back_unchanged(self, reg_ds, tmp_path, fixed_series, moving_series):
        path = tmp_path / "registration.dcm"
        reg_ds.save_as(path, enforce_file_format=True)

        parsed = get_spatial_registrations(dcmread(path))
        matrix = parsed[fixed_series[0].FrameOfReferenceUID][
            moving_series[0].FrameOfReferenceUID
        ]

        assert np.array(matrix, dtype=float).reshape(4, 4)[:3, 3] == pytest.approx(
            [5.0, -3.0, 2.0]
        )


class TestMatrixValidation:
    """
    A malformed matrix produces a file that loads and misplaces the image, so it has to be
    rejected at build time rather than discovered later.
    """

    def test_accepts_a_valid_rigid_matrix(self):
        result = validate_transformation_matrix(translation_matrix([1.0, 2.0, 3.0]))

        assert len(result) == 16
        assert result[3] == 1.0

    def test_accepts_a_4x4_array(self):
        result = validate_transformation_matrix(np.identity(4))

        assert result == np.identity(4).ravel().tolist()

    @pytest.mark.parametrize(
        "matrix, message",
        [
            ([1, 2, 3], "16 elements"),
            (np.identity(4).ravel().tolist()[:15] + [2.0], r"\[0, 0, 0, 1\]"),
            (np.diag([2.0, 1.0, 1.0, 1.0]).ravel().tolist(), "orthonormal"),
            (np.diag([-1.0, 1.0, 1.0, 1.0]).ravel().tolist(), "determinant"),
        ],
        ids=["too-short", "bad-bottom-row", "scaled", "reflection"],
    )
    def test_rejects_malformed_matrices(self, matrix, message):
        with pytest.raises(ValueError, match=message):
            validate_transformation_matrix(matrix)

    def test_rejects_non_finite_values(self):
        matrix = np.identity(4)
        matrix[0, 3] = np.nan

        with pytest.raises(ValueError, match="NaN or infinity"):
            validate_transformation_matrix(matrix)

    def test_scaling_is_allowed_when_declared(self):
        matrix = np.diag([2.0, 1.0, 1.0, 1.0]).ravel().tolist()

        assert validate_transformation_matrix(matrix, "RIGID_SCALE")[0] == 2.0

    def test_rejects_an_unknown_matrix_type(self):
        with pytest.raises(ValueError, match="matrix_type"):
            validate_transformation_matrix(np.identity(4), "ELASTIC")

    def test_builder_rejects_a_bad_matrix(self, fixed_series, moving_series):
        builder = SpatialRegistrationBuilder(fixed_series)

        with pytest.raises(ValueError, match="16 elements"):
            builder.add_registration(moving_series, [1, 2, 3])


class TestBuilderValidation:
    def test_empty_fixed_series_is_rejected(self):
        with pytest.raises(ValueError, match="fixed_series"):
            SpatialRegistrationBuilder([])

    def test_empty_moving_series_is_rejected(self, fixed_series):
        builder = SpatialRegistrationBuilder(fixed_series)

        with pytest.raises(ValueError, match="moving_series"):
            builder.add_registration([], translation_matrix([0.0, 0.0, 0.0]))

    def test_building_without_a_registration_is_rejected(self, fixed_series):
        """Used to emit a REG with an empty Registration Sequence, which is Type 1."""
        builder = SpatialRegistrationBuilder(fixed_series)

        with pytest.raises(ValueError, match="no registration added"):
            builder.build()

    def test_multiple_moving_series_each_get_an_item(self, fixed_series):
        second_moving = make_ct_series(shape=(4, 8, 8))
        third_moving = make_ct_series(shape=(4, 8, 8))

        builder = SpatialRegistrationBuilder(fixed_series)
        builder.add_registration(second_moving, translation_matrix([1.0, 0.0, 0.0]))
        builder.add_registration(third_moving, translation_matrix([0.0, 1.0, 0.0]))
        ds = builder.build()

        assert len(ds.RegistrationSequence) == 3  # two moving plus identity
        assert len(ds.ReferencedSeriesSequence) == 3
