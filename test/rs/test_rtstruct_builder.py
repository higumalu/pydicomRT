"""
RTStructBuilder.

Covers the builder wrapper specifically -- the mask/contour round trip itself is pinned in
test_rs_roundtrip.py.
"""

import numpy as np
import pytest

from synthetic import make_ct_series
from pydicomrt.rs import RTStructBuilder, check_rtstruct_iod, get_roi_names, rtstruct_to_masks
from pydicomrt.rs import calc_image_series_affine_mapping

SHAPE = (12, 64, 64)


def make_mask(z_range=(3, 9), radius=15, shape=SHAPE):
    mask = np.zeros(shape, dtype=np.uint8)
    rows, cols = np.mgrid[0:shape[1], 0:shape[2]]
    disc = ((cols - shape[2] // 2) ** 2 + (rows - shape[1] // 2) ** 2) < radius ** 2
    for z in range(*z_range):
        mask[z][disc] = 1
    return mask


@pytest.fixture
def image_series():
    return make_ct_series(shape=SHAPE, pixel_spacing=(2.0, 2.0), slice_spacing=3.0)


class TestBuilderShape:
    def test_calls_chain(self, image_series):
        """Every builder in the package returns self from set_* and add_*."""
        builder = RTStructBuilder(image_series)

        assert builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.") is builder
        assert builder.set_structure_set_label("PLAN_A") is builder
        assert builder.add_roi(mask=make_mask(), name="CTV") is builder

    def test_builds_in_one_expression(self, image_series):
        rs_ds = (
            RTStructBuilder(image_series)
            .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
            .add_roi(mask=make_mask(), name="CTV", color=[0, 255, 0])
            .build()
        )

        assert rs_ds.Modality == "RTSTRUCT"
        assert check_rtstruct_iod(rs_ds)["result"] is True
        assert rs_ds.SOPInstanceUID.startswith("1.2.826.0.1.3680043.2.1125.")

    def test_build_is_repeatable(self, image_series):
        """build() assembles from the accumulated spec, so it can be called twice."""
        builder = RTStructBuilder(image_series).add_roi(mask=make_mask(), name="CTV")

        first, second = builder.build(), builder.build()

        assert len(first.StructureSetROISequence) == len(second.StructureSetROISequence) == 1
        assert first.SOPInstanceUID != second.SOPInstanceUID


class TestRois:
    def test_roi_numbers_are_assigned_in_sequence(self, image_series):
        rs_ds = (
            RTStructBuilder(image_series)
            .add_roi(mask=make_mask(z_range=(2, 5)), name="First")
            .add_roi(mask=make_mask(z_range=(6, 9)), name="Second")
            .build()
        )

        assert [int(r.ROINumber) for r in rs_ds.StructureSetROISequence] == [1, 2]
        assert {int(k): v for k, v in get_roi_names(rs_ds).items()} == {1: "First", 2: "Second"}

    def test_explicit_roi_number_is_honoured(self, image_series):
        rs_ds = (
            RTStructBuilder(image_series)
            .add_roi(mask=make_mask(), name="CTV", number=7)
            .build()
        )

        assert int(rs_ds.StructureSetROISequence[0].ROINumber) == 7

    def test_duplicate_roi_number_is_rejected(self, image_series):
        builder = RTStructBuilder(image_series).add_roi(mask=make_mask(), name="A", number=3)

        with pytest.raises(ValueError, match="already added"):
            builder.add_roi(mask=make_mask(), name="B", number=3)

    def test_roi_metadata_reaches_the_dataset(self, image_series):
        rs_ds = (
            RTStructBuilder(image_series)
            .add_roi(
                mask=make_mask(), name="Cord", color=[0, 0, 255],
                description="spinal cord", interpreted_type="ORGAN",
            )
            .build()
        )

        assert rs_ds.StructureSetROISequence[0].ROIName == "Cord"
        assert rs_ds.StructureSetROISequence[0].ROIDescription == "spinal cord"
        assert list(rs_ds.ROIContourSequence[0].ROIDisplayColor) == [0, 0, 255]
        assert rs_ds.RTROIObservationsSequence[0].RTROIInterpretedType == "ORGAN"

    def test_mask_round_trips(self, image_series):
        mask = make_mask()
        rs_ds = RTStructBuilder(image_series).add_roi(mask=mask, name="CTV").build()

        affine, shape = calc_image_series_affine_mapping(image_series)
        recovered = np.asarray(rtstruct_to_masks(rs_ds, affine, shape)["CTV"]["mask_volume"])

        overlap = 2 * np.logical_and(mask > 0, recovered > 0).sum() / (mask.sum() + recovered.sum())
        assert overlap > 0.95

    def test_structure_set_label_override(self, image_series):
        rs_ds = (
            RTStructBuilder(image_series)
            .set_structure_set_label("PLAN_A")
            .add_roi(mask=make_mask(), name="CTV")
            .build()
        )

        assert rs_ds.StructureSetLabel == "PLAN_A"


class TestValidation:
    def test_empty_series_is_rejected(self):
        with pytest.raises(ValueError, match="image_series"):
            RTStructBuilder([])

    def test_building_without_an_roi_is_rejected(self, image_series):
        with pytest.raises(ValueError, match="no ROI added"):
            RTStructBuilder(image_series).build()

    def test_unnamed_roi_is_rejected(self, image_series):
        with pytest.raises(ValueError, match="needs a name"):
            RTStructBuilder(image_series).add_roi(mask=make_mask(), name="")

    def test_mask_shape_mismatch_is_rejected(self, image_series):
        """
        A mask of the wrong shape would otherwise be sliced against the wrong images and
        produce contours in the wrong place rather than an error.
        """
        wrong = make_mask(shape=(12, 32, 32))

        with pytest.raises(ValueError, match="does not match the image series"):
            RTStructBuilder(image_series).add_roi(mask=wrong, name="CTV")
