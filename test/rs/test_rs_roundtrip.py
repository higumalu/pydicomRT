"""
RTSTRUCT mask <-> contour round trip.

This is the library's most-used path and the one where a geometry mistake is hardest to
spot: a contour that is offset by half a voxel, mirrored, or scaled still loads and still
draws, it just sits in the wrong place. Comparing a mask against itself after a full round
trip is what makes that visible.
"""

import numpy as np
import pytest

from synthetic import make_ct_series
from pydicomrt.rs import (
    add_contour_sequence_from_mask3d,
    calc_image_series_affine_mapping,
    check_rtstruct_iod,
    create_roi_into_rs_ds,
    create_rtstruct_dataset,
    get_roi_names,
    is_rtstruct_matching_series,
    rtstruct_to_masks,
)

SERIES_SHAPE = (12, 64, 64)


def make_cylinder_mask(shape=SERIES_SHAPE, radius=15, z_range=(3, 9)):
    """A blob large enough to survive the contour noise filtering."""
    mask = np.zeros(shape, dtype=np.uint8)
    rows, cols = np.mgrid[0:shape[1], 0:shape[2]]
    disc = ((cols - shape[2] // 2) ** 2 + (rows - shape[1] // 2) ** 2) < radius ** 2
    for z in range(*z_range):
        mask[z][disc] = 1
    return mask


def dice(a, b):
    a, b = a > 0, b > 0
    total = a.sum() + b.sum()
    if total == 0:
        return 1.0
    return 2.0 * np.logical_and(a, b).sum() / total


@pytest.fixture
def series():
    return make_ct_series(shape=SERIES_SHAPE, pixel_spacing=(2.0, 2.0), slice_spacing=3.0)


@pytest.fixture
def rtstruct_with_cylinder(series):
    mask = make_cylinder_mask()
    rs_ds = create_rtstruct_dataset(series)
    rs_ds = create_roi_into_rs_ds(rs_ds, [255, 0, 0], 1, "TestROI", "round trip")
    rs_ds = add_contour_sequence_from_mask3d(rs_ds, series, 1, mask)
    return rs_ds, mask


class TestRoundTrip:
    def test_mask_survives_the_round_trip(self, series, rtstruct_with_cylinder):
        rs_ds, mask = rtstruct_with_cylinder

        affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
        recovered = np.asarray(
            rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)["TestROI"]["mask_volume"]
        )

        assert recovered.shape == mask.shape
        # Contour extraction resamples through polygon space, so exact equality is not the
        # bar; a geometry error would drop this far below 0.95, not shave a percent off it.
        assert dice(recovered, mask) > 0.95

    def test_round_trip_is_not_merely_returning_a_full_volume(
        self, series, rtstruct_with_cylinder
    ):
        """Guards the guard: a mask of everything would also score well against a big blob."""
        rs_ds, mask = rtstruct_with_cylinder

        affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
        recovered = np.asarray(
            rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)["TestROI"]["mask_volume"]
        )

        assert recovered.sum() < 0.5 * recovered.size
        assert dice(recovered, np.ones_like(mask)) < 0.5

    def test_contours_land_on_the_slices_that_held_the_mask(self, rtstruct_with_cylinder):
        rs_ds, mask = rtstruct_with_cylinder

        occupied_slices = int((mask.sum(axis=(1, 2)) > 0).sum())
        assert len(rs_ds.ROIContourSequence[0].ContourSequence) == occupied_slices

    @pytest.mark.parametrize(
        "pixel_spacing", [(2.0, 2.0), (3.0, 1.5), (1.0, 2.5)], ids=["square", "wide", "tall"]
    )
    def test_round_trip_holds_for_anisotropic_pixels(self, pixel_spacing):
        """
        Contour writing and mask reading derive their affines separately, so a row/column
        spacing mix-up in either one shows up as a stretched mask -- but only when the
        pixels are not square, which is exactly when square-pixel tests stay green.
        """
        series = make_ct_series(
            shape=SERIES_SHAPE, pixel_spacing=pixel_spacing, slice_spacing=3.0
        )
        mask = make_cylinder_mask()

        rs_ds = create_rtstruct_dataset(series)
        rs_ds = create_roi_into_rs_ds(rs_ds, [255, 0, 0], 1, "TestROI", "")
        rs_ds = add_contour_sequence_from_mask3d(rs_ds, series, 1, mask)

        affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
        recovered = np.asarray(
            rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)["TestROI"]["mask_volume"]
        )

        assert dice(recovered, mask) > 0.95

    def test_contour_points_are_triplets_in_patient_space(self, series, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder
        contour = rs_ds.ROIContourSequence[0].ContourSequence[0]

        assert contour.ContourGeometricType == "CLOSED_PLANAR"
        assert len(contour.ContourData) % 3 == 0
        assert int(contour.NumberOfContourPoints) == len(contour.ContourData) // 3

        # Every point must sit on one of the series slice planes
        z_values = {round(float(ds.ImagePositionPatient[2]), 3) for ds in series}
        contour_z = {round(float(z), 3) for z in contour.ContourData[2::3]}
        assert contour_z <= z_values


class TestRoiMetadata:
    def test_roi_name_and_number_survive(self, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder

        roi_map = get_roi_names(rs_ds)
        assert {int(number): name for number, name in roi_map.items()} == {1: "TestROI"}
        assert rs_ds.StructureSetROISequence[0].ROIName == "TestROI"
        assert rs_ds.ROIContourSequence[0].ReferencedROINumber == 1
        assert list(rs_ds.ROIContourSequence[0].ROIDisplayColor) == [255, 0, 0]

    def test_roi_numbers_are_returned_as_pydicom_is_not_str(self, rtstruct_with_cylinder):
        """
        The keys print like strings but are ``pydicom.valuerep.IS`` (an int subclass), so
        ``roi_map[1]`` finds the ROI and ``roi_map["1"]`` raises KeyError. Pinned here
        because the repr makes the opposite look true.
        """
        rs_ds, _ = rtstruct_with_cylinder
        roi_map = get_roi_names(rs_ds)

        assert roi_map[1] == "TestROI"
        with pytest.raises(KeyError):
            roi_map["1"]

    def test_multiple_rois_stay_independent(self, series):
        first = make_cylinder_mask(radius=15, z_range=(3, 6))
        second = make_cylinder_mask(radius=10, z_range=(7, 10))

        rs_ds = create_rtstruct_dataset(series)
        rs_ds = create_roi_into_rs_ds(rs_ds, [255, 0, 0], 1, "First", "")
        rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 2, "Second", "")
        rs_ds = add_contour_sequence_from_mask3d(rs_ds, series, 1, first)
        rs_ds = add_contour_sequence_from_mask3d(rs_ds, series, 2, second)

        affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
        masks = rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)

        assert set(masks) == {"First", "Second"}
        assert dice(np.asarray(masks["First"]["mask_volume"]), first) > 0.95
        assert dice(np.asarray(masks["Second"]["mask_volume"]), second) > 0.95

    def test_unknown_roi_number_is_rejected(self, series):
        rs_ds = create_rtstruct_dataset(series)
        rs_ds = create_roi_into_rs_ds(rs_ds, [255, 0, 0], 1, "First", "")

        with pytest.raises(ValueError, match="roi_number not found"):
            add_contour_sequence_from_mask3d(rs_ds, series, 99, make_cylinder_mask())


class TestStructureSetLinkage:
    def test_rtstruct_matches_the_series_it_was_built_from(self, series, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder

        assert is_rtstruct_matching_series(rs_ds, series)

    def test_rtstruct_does_not_match_an_unrelated_series(self, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder
        other_series = make_ct_series(shape=SERIES_SHAPE)

        assert not is_rtstruct_matching_series(rs_ds, other_series)

    def test_inherits_patient_and_study_context(self, series, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder

        assert rs_ds.PatientID == series[0].PatientID
        assert rs_ds.StudyInstanceUID == series[0].StudyInstanceUID
        assert rs_ds.Modality == "RTSTRUCT"

    def test_passes_its_own_iod_check(self, rtstruct_with_cylinder):
        rs_ds, _ = rtstruct_with_cylinder

        result = check_rtstruct_iod(rs_ds)

        assert result["result"] is True, result.get("content")
