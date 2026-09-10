"""
RT Dose builder.

Dose is stored as scaled integers, so the failure modes are quiet arithmetic ones: a wrong
signedness declaration turns high dose negative, and a scaling factor that does not fit
wraps around. Both produce a file that opens.
"""

import numpy as np
import pytest
import SimpleITK as sitk
from pydicom import dcmread
from pydicom.dataset import Dataset
from pydicom.uid import generate_uid

from synthetic import make_ct_series
from pydicomrt.dose import (
    RTDoseBuilder,
    add_dose_grid_to_ds,
    check_rtdose_iod,
    cp_information_from_ds,
    generate_base_dataset,
    get_dose_image,
)
from pydicomrt.dose.builder import (
    add_rf_rt_plan_seq_from_plan_ds,
    check_dose_plan_reference,
)

RT_PLAN_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.5"


def make_dose_image(peak=60.0, size=(12, 10, 8), spacing=(2.5, 3.5, 4.0), origin=(-30.0, -20.0, -10.0)):
    """A smooth dose blob. Anisotropic spacing so axis mix-ups are visible."""
    z, y, x = np.mgrid[0:size[2], 0:size[1], 0:size[0]]
    array = peak * np.exp(
        -(((x - size[0] / 2) ** 2) / 20.0
          + ((y - size[1] / 2) ** 2) / 20.0
          + ((z - size[2] / 2) ** 2) / 10.0)
    )
    image = sitk.GetImageFromArray(array.astype(np.float32))
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    return image


def make_plan_ds():
    plan = Dataset()
    plan.SOPClassUID = RT_PLAN_SOP_CLASS_UID
    plan.SOPInstanceUID = generate_uid()
    plan.file_meta = Dataset()
    plan.file_meta.MediaStorageSOPClassUID = RT_PLAN_SOP_CLASS_UID
    return plan


@pytest.fixture
def reference_ds():
    ds = make_ct_series(shape=(4, 8, 8))[0]
    # Distinctive, so "the dose kept its own series number" is a real assertion rather than
    # a coincidence between two defaults.
    ds.SeriesNumber = 7
    return ds


@pytest.fixture
def dose_ds(reference_ds):
    ds = generate_base_dataset()
    ds = cp_information_from_ds(ds, reference_ds)
    ds = add_dose_grid_to_ds(ds, make_dose_image())
    add_rf_rt_plan_seq_from_plan_ds(ds, make_plan_ds())
    return ds


class TestGeometry:
    def test_grid_round_trips(self, dose_ds):
        source = make_dose_image()

        recovered = get_dose_image(dose_ds)

        assert recovered.GetSize() == source.GetSize()
        assert recovered.GetSpacing() == pytest.approx(source.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(source.GetOrigin())

    def test_pixel_spacing_is_row_then_column(self, dose_ds):
        # source spacing is (x, y, z) = (2.5, 3.5, 4.0)
        assert list(dose_ds.PixelSpacing) == pytest.approx([3.5, 2.5])
        assert dose_ds.Rows == 10
        assert dose_ds.Columns == 12

    def test_orientation_is_vm_six(self, dose_ds):
        """It used to be written with all nine direction elements."""
        assert len(dose_ds.ImageOrientationPatient) == 6

    def test_frame_offsets_match_the_slice_spacing(self, dose_ds):
        offsets = np.array(dose_ds.GridFrameOffsetVector, dtype=float)

        assert len(offsets) == dose_ds.NumberOfFrames
        assert np.diff(offsets) == pytest.approx(4.0)


class TestDoseValues:
    def test_dose_round_trips_within_the_scaling(self, dose_ds):
        source = sitk.GetArrayFromImage(make_dose_image())

        recovered = sitk.GetArrayFromImage(get_dose_image(dose_ds))

        assert np.abs(recovered - source).max() < dose_ds.DoseGridScaling * 2

    def test_stored_values_are_declared_unsigned(self, dose_ds):
        """
        The data is written as uint32. Declaring it signed made every dose above about
        215 Gy read back negative.
        """
        assert dose_ds.PixelRepresentation == 0
        assert dose_ds.BitsAllocated == 32

    def test_high_dose_does_not_wrap(self, reference_ds):
        """300 Gy is past the signed 32-bit boundary at the default scaling."""
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        ds = add_dose_grid_to_ds(ds, make_dose_image(peak=300.0))

        recovered = sitk.GetArrayFromImage(get_dose_image(ds))

        assert recovered.max() == pytest.approx(300.0, rel=1e-4)
        assert recovered.min() >= 0.0

    def test_dose_beyond_the_default_scaling_rescales_instead_of_overflowing(self, reference_ds):
        """The default 1e-7 tops out near 429 Gy; the builder raises the factor rather than wrap."""
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        ds = add_dose_grid_to_ds(ds, make_dose_image(peak=5000.0))

        assert ds.DoseGridScaling > 1e-7
        recovered = sitk.GetArrayFromImage(get_dose_image(ds))
        assert recovered.max() == pytest.approx(5000.0, rel=1e-4)

    def test_explicit_scaling_that_does_not_fit_is_rejected(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)

        with pytest.raises(ValueError, match="does not fit"):
            add_dose_grid_to_ds(ds, make_dose_image(peak=5000.0), dose_grid_scaling=1e-9)

    def test_negative_dose_is_rejected(self, reference_ds):
        """Unsigned storage would turn a negative dose into an enormous positive one."""
        image = make_dose_image()
        array = sitk.GetArrayFromImage(image)
        array[0, 0, 0] = -1.0
        negative = sitk.GetImageFromArray(array)
        negative.CopyInformation(image)

        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        with pytest.raises(ValueError, match="negative"):
            add_dose_grid_to_ds(ds, negative)

    def test_non_positive_scaling_is_rejected(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)

        with pytest.raises(ValueError, match="must be positive"):
            add_dose_grid_to_ds(ds, make_dose_image(), dose_grid_scaling=0.0)


class TestPlanReference:
    def test_plan_requiring_summation_needs_a_plan(self, reference_ds):
        """
        DoseSummationType PLAN makes Referenced RT Plan Sequence Type 1C. The builder leaves
        the sequence empty, so a dose written without adding one claimed a summation it
        could not substantiate.
        """
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        ds = add_dose_grid_to_ds(ds, make_dose_image(), dose_summation_type="PLAN")

        problems = check_dose_plan_reference(ds)

        assert problems
        assert "ReferencedRTPlanSequence" in problems[0]

    def test_satisfied_once_a_plan_is_referenced(self, dose_ds):
        assert check_dose_plan_reference(dose_ds) == []

    def test_summation_type_without_the_requirement(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        ds = add_dose_grid_to_ds(ds, make_dose_image(), dose_summation_type="RECORD")

        assert check_dose_plan_reference(ds) == []

    def test_non_plan_reference_is_rejected(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        not_a_plan = make_plan_ds()
        not_a_plan.file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"

        with pytest.raises(ValueError, match="not an RT Plan"):
            add_rf_rt_plan_seq_from_plan_ds(ds, not_a_plan)


class TestDatasetIdentity:
    def test_passes_the_iod_check(self, dose_ds):
        result = check_rtdose_iod(dose_ds)

        assert result["result"] is True, result["content"]

    def test_iod_check_catches_the_missing_plan_reference(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds)
        ds = add_dose_grid_to_ds(ds, make_dose_image())

        assert check_rtdose_iod(ds)["result"] is False

    def test_instance_number_is_present(self, dose_ds):
        """General Image module, Type 2 -- was missing entirely."""
        assert dose_ds.InstanceNumber == 1

    def test_gets_its_own_series(self, dose_ds, reference_ds):
        """
        The dose is a new series. It used to copy the reference image's SeriesNumber, which
        made it claim membership of the series it was merely derived from.
        """
        assert dose_ds.SeriesInstanceUID != reference_ds.SeriesInstanceUID
        assert int(dose_ds.SeriesNumber) == 1
        assert int(reference_ds.SeriesNumber) == 7
        assert dose_ds.Modality == "RTDOSE"

    def test_series_number_is_settable(self, reference_ds):
        ds = cp_information_from_ds(generate_base_dataset(), reference_ds, series_number=42)

        assert int(ds.SeriesNumber) == 42

    def test_inherits_patient_study_and_frame_of_reference(self, dose_ds, reference_ds):
        assert dose_ds.PatientID == reference_ds.PatientID
        assert dose_ds.StudyInstanceUID == reference_ds.StudyInstanceUID
        assert dose_ds.FrameOfReferenceUID == reference_ds.FrameOfReferenceUID

    def test_default_dose_type_is_physical(self, dose_ds):
        """EFFECTIVE means biologically weighted, which a plain dose grid is not."""
        assert dose_ds.DoseType == "PHYSICAL"

    def test_survives_write_and_read(self, dose_ds, tmp_path):
        path = tmp_path / "rtdose.dcm"
        dose_ds.save_as(path, enforce_file_format=True)

        reloaded = dcmread(path)

        assert check_rtdose_iod(reloaded)["result"] is True
        recovered = get_dose_image(reloaded)
        assert recovered.GetSpacing() == pytest.approx(make_dose_image().GetSpacing())


class TestRTDoseBuilder:
    """The class wrapper; the underlying encoding is covered above."""

    def test_calls_chain(self, reference_ds):
        builder = RTDoseBuilder(reference_ds)

        assert builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.") is builder
        assert builder.set_series_number(3) is builder
        assert builder.set_dose_grid(make_dose_image()) is builder
        assert builder.add_referenced_plan(make_plan_ds()) is builder

    def test_builds_in_one_expression(self, reference_ds):
        dose_ds = (
            RTDoseBuilder(reference_ds)
            .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
            .set_dose_grid(make_dose_image())
            .add_referenced_plan(make_plan_ds())
            .build()
        )

        assert check_rtdose_iod(dose_ds)["result"] is True
        assert dose_ds.SOPInstanceUID.startswith("1.2.826.0.1.3680043.2.1125.")
        assert dose_ds.SeriesInstanceUID.startswith("1.2.826.0.1.3680043.2.1125.")

    def test_uid_prefix_is_not_a_no_op(self, reference_ds):
        """set_uid_prefix() used to be stored and then ignored by the underlying functions."""
        default_ds = RTDoseBuilder(reference_ds).set_dose_grid(make_dose_image()).build()
        prefixed_ds = (
            RTDoseBuilder(reference_ds)
            .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
            .set_dose_grid(make_dose_image())
            .build()
        )

        assert not default_ds.SOPInstanceUID.startswith("1.2.826.0.1.3680043.2.1125.")
        assert prefixed_ds.SOPInstanceUID.startswith("1.2.826.0.1.3680043.2.1125.")

    def test_series_number_is_settable(self, reference_ds):
        dose_ds = (
            RTDoseBuilder(reference_ds)
            .set_series_number(42)
            .set_dose_grid(make_dose_image())
            .build()
        )

        assert int(dose_ds.SeriesNumber) == 42

    def test_grid_options_pass_through(self, reference_ds):
        dose_ds = (
            RTDoseBuilder(reference_ds)
            .set_dose_grid(
                make_dose_image(), dose_type="EFFECTIVE",
                dose_summation_type="BEAM", dose_units="RELATIVE",
            )
            .add_referenced_plan(make_plan_ds())
            .build()
        )

        assert dose_ds.DoseType == "EFFECTIVE"
        assert dose_ds.DoseSummationType == "BEAM"
        assert dose_ds.DoseUnits == "RELATIVE"

    def test_building_without_a_grid_is_rejected(self, reference_ds):
        with pytest.raises(ValueError, match="no dose grid set"):
            RTDoseBuilder(reference_ds).build()

    def test_a_non_plan_reference_is_rejected_at_build(self, reference_ds):
        not_a_plan = make_plan_ds()
        not_a_plan.file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"

        builder = (
            RTDoseBuilder(reference_ds)
            .set_dose_grid(make_dose_image())
            .add_referenced_plan(not_a_plan)
        )
        with pytest.raises(ValueError, match="not an RT Plan"):
            builder.build()
