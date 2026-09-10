"""
Compare what SpatialRegistrationBuilder emits against a registration a TPS actually wrote.

The synthetic builder tests check the object against my reading of PS3.3. This checks it
against a file that real clinical software produced and real clinical software consumes,
which is the only thing that catches a misreading of the standard rather than a regression
against it.

Structure only. The matrices differ -- one is a real registration, the other is built from a
synthetic series -- so nothing here compares values.
"""

import numpy as np
import pytest

from pydicomrt.reg import SpatialRegistrationBuilder, check_spatial_reg_iod


@pytest.fixture(scope="module")
def built_reg(planning_ct_series, cbct_series):
    """
    A registration built by this library from the same two real series.

    The matrix is deliberately not the identity, so the moving item stays distinguishable
    from the identity item the builder appends for the fixed frame.
    """
    moving_to_fixed = np.identity(4)
    moving_to_fixed[:3, 3] = [-15.0, -178.0, -31.0]

    builder = SpatialRegistrationBuilder(planning_ct_series)
    builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    builder.add_registration(cbct_series, moving_to_fixed.ravel().tolist())
    return builder.build()


def top_level_keywords(ds):
    return {element.keyword for element in ds if element.keyword}


class TestAgainstTheClinicalObject:
    def test_both_pass_the_same_iod_check(self, built_reg, clinical_reg_ds):
        """
        If the checker only passed our own output it would be describing the builder rather
        than the standard.
        """
        assert check_spatial_reg_iod(clinical_reg_ds)["result"] is True
        assert check_spatial_reg_iod(built_reg)["result"] is True

    def test_no_required_top_level_element_is_missing(self, built_reg, clinical_reg_ds):
        """
        Anything the TPS writes and we do not is worth knowing about. A short allow-list
        covers what is genuinely vendor-specific or not ours to invent.
        """
        # All Type 3, and all describing the acquisition device rather than the software
        # that wrote the object -- not ours to invent.
        vendor_specific = {
            "DeviceSerialNumber",
            "StationName",
            "InstitutionalDepartmentName",
        }

        missing = top_level_keywords(clinical_reg_ds) - top_level_keywords(built_reg)

        assert not (missing - vendor_specific), f"not written by the builder: {missing}"

    def test_registration_sequence_has_the_same_shape(self, built_reg, clinical_reg_ds):
        """One item per frame of reference, including the identity item for the RCS."""
        assert len(clinical_reg_ds.RegistrationSequence) == 2
        assert len(built_reg.RegistrationSequence) == 2

    def test_registration_items_carry_the_same_elements(self, built_reg, clinical_reg_ds):
        clinical_keywords = {
            element.keyword
            for item in clinical_reg_ds.RegistrationSequence
            for element in item
            if element.keyword
        }
        built_keywords = {
            element.keyword
            for item in built_reg.RegistrationSequence
            for element in item
            if element.keyword
        }

        assert clinical_keywords <= built_keywords, (
            f"registration item elements the TPS writes and we do not: "
            f"{clinical_keywords - built_keywords}"
        )

    def test_registration_type_code_is_nested_the_same_way(self, built_reg, clinical_reg_ds):
        """
        The TPS puts it inside Matrix Registration Sequence. The builder used to put it in
        the Registration Sequence item, one level too high.
        """
        for ds in (clinical_reg_ds, built_reg):
            for item in ds.RegistrationSequence:
                assert not hasattr(item, "RegistrationTypeCodeSequence")
                assert hasattr(
                    item.MatrixRegistrationSequence[0], "RegistrationTypeCodeSequence"
                )

    def test_identity_item_uses_the_same_code(self, built_reg, clinical_reg_ds):
        def identity_codes(ds):
            codes = set()
            for item in ds.RegistrationSequence:
                matrix_registration = item.MatrixRegistrationSequence[0]
                matrix = np.array(
                    matrix_registration.MatrixSequence[0]
                    .FrameOfReferenceTransformationMatrix,
                    dtype=float,
                ).reshape(4, 4)
                if np.allclose(matrix, np.identity(4)):
                    codes.add(matrix_registration.RegistrationTypeCodeSequence[0].CodeValue)
            return codes

        assert identity_codes(clinical_reg_ds) == {"125021"}
        assert identity_codes(built_reg) == {"125021"}

    def test_frame_of_reference_is_the_planning_ct(
        self, built_reg, clinical_reg_ds, planning_ct_series
    ):
        """The registered RCS is the fixed series' frame in both."""
        expected = planning_ct_series[0].FrameOfReferenceUID

        assert clinical_reg_ds.FrameOfReferenceUID == expected
        assert built_reg.FrameOfReferenceUID == expected

    def test_referenced_series_sequence_covers_both_series(
        self, built_reg, clinical_reg_ds, planning_ct_series, cbct_series
    ):
        expected = {
            planning_ct_series[0].SeriesInstanceUID,
            cbct_series[0].SeriesInstanceUID,
        }

        clinical = {item.SeriesInstanceUID for item in clinical_reg_ds.ReferencedSeriesSequence}
        built = {item.SeriesInstanceUID for item in built_reg.ReferencedSeriesSequence}

        assert clinical == expected
        assert built == expected

    def test_every_instance_of_each_series_is_referenced(
        self, built_reg, planning_ct_series, cbct_series
    ):
        counts = {
            item.SeriesInstanceUID: len(item.ReferencedInstanceSequence)
            for item in built_reg.ReferencedSeriesSequence
        }

        assert counts[planning_ct_series[0].SeriesInstanceUID] == len(planning_ct_series)
        assert counts[cbct_series[0].SeriesInstanceUID] == len(cbct_series)
