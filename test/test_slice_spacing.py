"""Series whose slices are not evenly spaced: the warning, and slice_spacing="resample"."""
import warnings

import numpy as np
import pytest
import SimpleITK as sitk
from pydicom.uid import generate_uid

from pydicomrt.rs import RTStructBuilder, calc_image_series_affine_mapping
from pydicomrt.utils import (
    NonUniformSliceSpacingWarning,
    get_pixel_to_patient_transformation_matrix,
    get_spacing_between_slices,
    image_series_to_sitk_image,
    parse_image_series,
)
from synthetic import make_ct_slice

ORIGIN = np.array([-50.0, -60.0, -20.0])
AXIAL = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
# Rows run along x, columns along a y tilted out of the axial plane: the slice normal is
# not the z axis, so only a projection onto the normal measures the gaps correctly.
TILT = np.deg2rad(25)
TILTED = (1.0, 0.0, 0.0, 0.0, np.cos(TILT), np.sin(TILT))


def series_at(offsets, orientation=AXIAL, shear=0.0):
    """
    A CT series with slices ``offsets`` mm along the slice normal; HU = 10 * offset + 2 * row.

    ``shear`` slides each slice along its columns by that many mm per mm of offset, as a
    tilted gantry does: the slice origins then do not lie along the normal.
    """
    normal = np.cross(orientation[:3], orientation[3:]) + shear * np.asarray(orientation[3:])
    uids = generate_uid(), generate_uid(), generate_uid()
    rows = np.arange(6).reshape(6, 1) * np.ones((1, 5))
    return [
        make_ct_slice(
            pixel_array=np.round(1024 + 10 * offset + 2 * rows).astype(np.uint16),
            image_position=tuple(ORIGIN + offset * normal),
            pixel_spacing=(2.0, 1.5),
            slice_thickness=3.0,
            study_uid=uids[0], series_uid=uids[1], frame_of_reference_uid=uids[2],
            instance_number=index + 1,
            orientation=orientation,
        )
        for index, offset in enumerate(offsets)
    ]


def assert_same_image(image, expected):
    assert image.GetSize() == expected.GetSize()
    assert image.GetSpacing() == pytest.approx(expected.GetSpacing(), abs=1e-9)
    assert image.GetOrigin() == pytest.approx(expected.GetOrigin(), abs=1e-9)
    assert image.GetDirection() == pytest.approx(expected.GetDirection(), abs=1e-9)
    assert sitk.GetArrayFromImage(image) == pytest.approx(sitk.GetArrayFromImage(expected), abs=1e-6)


EVEN = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0]


@pytest.mark.parametrize("offsets", [
    EVEN,
    [0.0, 3.004, 5.997, 9.0, 12.003, 15.0, 18.0],   # the rounding a scanner writes
    [0.0],
], ids=["even", "rounding", "single slice"])
def test_evenly_spaced_series_do_not_warn_and_read_the_same_either_way(offsets):
    series = series_at(offsets)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean = image_series_to_sitk_image(series)
        resampled = image_series_to_sitk_image(series, slice_spacing="resample")
        get_pixel_to_patient_transformation_matrix(series)
    assert_same_image(resampled, mean)
    assert sitk.GetArrayFromImage(resampled).tolist() == sitk.GetArrayFromImage(mean).tolist()


@pytest.mark.parametrize("offsets, description", [
    ([0.0, 3.0, 6.0, 12.0, 15.0, 18.0], "slices 2 and 3 are 6.000 mm apart (1 slice missing)"),
    ([0.0, 3.0, 6.0, 15.0, 18.0], "slices 2 and 3 are 9.000 mm apart (2 slices missing)"),
    ([0.0, 3.0, 6.0, 6.0, 9.0, 12.0], "slices 2 and 3 are 0.000 mm apart (duplicated position)"),
    ([0.0, 3.0, 6.0, 9.0, 10.0, 13.0, 16.0], "slices 3 and 4 are 1.000 mm apart (different spacing)"),
], ids=["missing", "two missing", "duplicate", "odd gap"])
def test_uneven_spacing_warns_and_says_where(offsets, description):
    series = series_at(offsets)
    with pytest.warns(NonUniformSliceSpacingWarning) as caught:
        spacing = get_spacing_between_slices(series)
    # The value itself is what it always was: existing code keeps its result.
    assert spacing == pytest.approx((offsets[-1] - offsets[0]) / (len(offsets) - 1))
    assert description in str(caught[0].message)
    assert 'slice_spacing="resample"' in str(caught[0].message)


def test_every_route_from_series_to_geometry_warns():
    series = series_at([0.0, 3.0, 6.0, 12.0, 15.0])
    for read in (
        image_series_to_sitk_image,
        parse_image_series,
        get_pixel_to_patient_transformation_matrix,
        calc_image_series_affine_mapping,
    ):
        with pytest.warns(NonUniformSliceSpacingWarning):
            read(series)


def test_the_warning_can_be_made_an_error():
    series = series_at([0.0, 3.0, 6.0, 12.0, 15.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error", NonUniformSliceSpacingWarning)
        with pytest.raises(NonUniformSliceSpacingWarning):
            image_series_to_sitk_image(series)


@pytest.mark.parametrize("orientation", [AXIAL, TILTED], ids=["axial", "tilted"])
@pytest.mark.parametrize("offsets", [
    [0.0, 3.0, 6.0, 12.0, 15.0, 18.0],
    [0.0, 3.0, 6.0, 6.0, 9.0, 12.0, 12.0, 15.0, 18.0],
    [0.0, 3.0, 9.0, 9.0, 12.0, 15.0, 18.0],
], ids=["missing", "duplicates", "missing and duplicate"])
def test_resample_recovers_the_evenly_spaced_image(offsets, orientation):
    # HU is linear along the normal, so linear interpolation must rebuild a lost slice exactly.
    expected = image_series_to_sitk_image(series_at(EVEN, orientation))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        image = image_series_to_sitk_image(series_at(offsets, orientation), slice_spacing="resample")
    assert_same_image(image, expected)


def test_mean_spacing_misplaces_what_resample_places():
    offsets = [0.0, 3.0, 6.0, 12.0, 15.0, 18.0]
    with pytest.warns(NonUniformSliceSpacingWarning):
        mean = image_series_to_sitk_image(series_at(offsets))
    resampled = image_series_to_sitk_image(series_at(offsets), slice_spacing="resample")
    # 9 mm up is where the lost slice was; HU there is 10 * 9 at row 0. Spaced by the mean
    # (3.6 mm) the stack has no slice there, only its neighbours from 6 and 12 mm.
    point = tuple(ORIGIN + (0.0, 0.0, 9.0))
    assert resampled[resampled.TransformPhysicalPointToIndex(point)] == pytest.approx(90.0)
    assert mean[mean.TransformPhysicalPointToIndex(point)] in (60.0, 120.0)


def test_resample_uses_the_finest_spacing_of_a_mixed_series():
    offsets = [0.0, 5.0, 10.0, 15.0, 17.5, 20.0, 22.5, 25.0]
    image = image_series_to_sitk_image(series_at(offsets), slice_spacing="resample")
    assert image.GetSpacing()[2] == pytest.approx(2.5)
    assert image.GetSize()[2] == 11
    column = sitk.GetArrayFromImage(image)[:, 0, 0]
    assert column == pytest.approx(10 * 2.5 * np.arange(11), abs=1e-6)


def test_resample_stops_at_the_last_slice_when_the_span_is_not_a_multiple():
    image = image_series_to_sitk_image(series_at([0.0, 2.0, 4.0, 6.0, 9.0]), slice_spacing="resample")
    assert image.GetSpacing()[2] == pytest.approx(2.0)
    assert image.GetSize()[2] == 5            # 0, 2, 4, 6, 8 -- not past 9
    assert sitk.GetArrayFromImage(image)[:, 0, 0] == pytest.approx([0, 20, 40, 60, 80], abs=1e-6)


def test_unknown_mode_is_rejected():
    with pytest.raises(ValueError, match="slice_spacing must be one of"):
        image_series_to_sitk_image(series_at(EVEN), slice_spacing="median")


# ---------------------------------------------------------------------------------------- #
# Contours are placed by the slice they reference
# ---------------------------------------------------------------------------------------- #

def contours_by_slice(series, mask):
    """``{slice index: (n, 3) points}`` for a one-ROI RTSTRUCT built from ``mask``."""
    index_of = {ds.SOPInstanceUID: index for index, ds in enumerate(series)}
    dataset = RTStructBuilder(series).add_roi(mask, name="roi").build()
    found = {}
    for contour in dataset.ROIContourSequence[0].ContourSequence:
        index = index_of[contour.ContourImageSequence[0].ReferencedSOPInstanceUID]
        found[index] = np.array(contour.ContourData, dtype=float).reshape(-1, 3)
    return found


@pytest.mark.parametrize("orientation", [AXIAL, TILTED], ids=["axial", "tilted"])
@pytest.mark.parametrize("offsets, shear", [
    (EVEN, 0.0),
    ([0.0, 3.0, 6.0, 12.0, 15.0, 18.0], 0.0),
    ([0.0, 3.0, 6.0, 6.0, 9.0, 12.0], 0.0),
    ([0.0, 5.0, 10.0, 12.5, 15.0, 17.5], 0.0),
    (EVEN, 0.3),
    ([0.0, 3.0, 6.0, 12.0, 15.0, 18.0], 0.3),
], ids=["even", "missing", "duplicate", "mixed", "gantry tilt", "gantry tilt, missing"])
def test_contours_lie_in_the_plane_of_the_slice_they_reference(offsets, shear, orientation):
    # A TPS drops a contour that is not in the plane of the image it references. Placing
    # contours by the series-wide matrix (mean spacing, origins along the normal) fails that
    # on every case here but the first.
    series = series_at(offsets, orientation, shear)
    mask = np.zeros((len(series), 6, 5), dtype=bool)
    mask[:, 1:5, 1:4] = True
    contours = contours_by_slice(series, mask)
    assert sorted(contours) == list(range(len(series)))

    for index, points in contours.items():
        image_slice = series[index]
        position = np.array(image_slice.ImagePositionPatient, dtype=float)
        row_direction = np.array(image_slice.ImageOrientationPatient[:3], dtype=float)
        column_direction = np.array(image_slice.ImageOrientationPatient[3:], dtype=float)
        relative = points - position
        assert relative @ np.cross(row_direction, column_direction) == pytest.approx(0.0, abs=1e-6)
        # ... and around the pixels the mask marks, columns 1-3 and rows 1-4 (the outline is
        # smoothed, so not exactly on them). PixelSpacing is (2.0, 1.5).
        columns = relative @ row_direction / 1.5
        rows = relative @ column_direction / 2.0
        assert 0.5 < columns.min() < 1.5 and 2.5 < columns.max() < 3.5
        assert 0.5 < rows.min() < 1.5 and 3.5 < rows.max() < 4.5


def test_contours_on_an_even_series_are_where_the_series_matrix_puts_them():
    # The per-slice placement is the same mapping as the series-wide one wherever that one holds.
    series = series_at(EVEN, TILTED)
    mask = np.zeros((len(series), 6, 5), dtype=bool)
    mask[:, 1:5, 1:4] = True
    matrix = get_pixel_to_patient_transformation_matrix(series).astype(float)
    inverse = np.linalg.inv(matrix)
    for index, points in contours_by_slice(series, mask).items():
        pixels = np.c_[points, np.ones(len(points))] @ inverse.T
        assert pixels[:, 2] == pytest.approx(index, abs=1e-4)
        assert 0.5 < pixels[:, 0].min() < 1.5 and 2.5 < pixels[:, 0].max() < 3.5
        assert 0.5 < pixels[:, 1].min() < 1.5 and 3.5 < pixels[:, 1].max() < 4.5
