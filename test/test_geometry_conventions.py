"""
DICOM <-> SimpleITK geometry conventions.

These are the conversions that silently produce a plausible-but-wrong image: nothing raises,
the volume looks fine, and only the physical coordinates are off. Square in-plane pixels and
axis-aligned axial series hide every mistake here, so each test deliberately uses
anisotropic spacing or a rotated orientation.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_ct_series
from pydicomrt.utils.coordinate_transform import (
    get_patient_to_pixel_transformation_matrix,
    get_pixel_to_patient_transformation_matrix,
)
from pydicomrt.utils.sitk_transform import (
    SimpleITKImageBuilder,
    parse_image_series,
    sitk_direction_to_image_orientation_patient,
    sitk_spacing_to_pixel_spacing,
)

# Row direction +y, column direction -x: a 90 degree in-plane rotation. The resulting
# direction matrix is not symmetric, so a transposed matrix gives different coordinates.
ROTATED_ORIENTATION = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0]


class TestPixelSpacingOrder:
    """DICOM PixelSpacing is [row spacing, column spacing]; SimpleITK spacing is (x, y, z)."""

    def test_row_and_column_spacing_map_to_the_right_axes(self, anisotropic_ct_series):
        # fixture uses PixelSpacing = [3.0 (between rows), 1.5 (between columns)]
        _, _, spacing, _ = parse_image_series(anisotropic_ct_series)

        # x follows the column index, so it takes the column spacing
        assert spacing[0] == pytest.approx(1.5)
        assert spacing[1] == pytest.approx(3.0)
        assert spacing[2] == pytest.approx(5.0)

    def test_stepping_one_column_moves_by_the_column_spacing(self, anisotropic_ct_series):
        image = SimpleITKImageBuilder().from_image_series(anisotropic_ct_series)

        origin = np.array(image.TransformIndexToPhysicalPoint((0, 0, 0)))
        one_column = np.array(image.TransformIndexToPhysicalPoint((1, 0, 0)))
        one_row = np.array(image.TransformIndexToPhysicalPoint((0, 1, 0)))

        assert np.linalg.norm(one_column - origin) == pytest.approx(1.5)
        assert np.linalg.norm(one_row - origin) == pytest.approx(3.0)

    def test_pixel_spacing_round_trips_through_the_sitk_helper(self):
        assert sitk_spacing_to_pixel_spacing((1.5, 3.0, 5.0)) == [3.0, 1.5]


class TestDirectionMatrix:
    """A SimpleITK direction matrix holds the axis directions in its columns."""

    def test_direction_columns_are_the_dicom_direction_vectors(self):
        ds_list = make_ct_series(shape=(3, 8, 8), orientation=ROTATED_ORIENTATION)
        _, _, _, direction = parse_image_series(ds_list)

        assert direction[:, 0] == pytest.approx(ROTATED_ORIENTATION[:3])
        assert direction[:, 1] == pytest.approx(ROTATED_ORIENTATION[3:])
        assert direction[:, 2] == pytest.approx(
            np.cross(ROTATED_ORIENTATION[:3], ROTATED_ORIENTATION[3:])
        )

    def test_rotated_series_places_voxels_where_dicom_says(self):
        origin = (10.0, 20.0, 30.0)
        column_spacing = 2.0
        ds_list = make_ct_series(
            shape=(3, 8, 8),
            origin=origin,
            pixel_spacing=(column_spacing, column_spacing),
            orientation=ROTATED_ORIENTATION,
        )
        image = SimpleITKImageBuilder().from_image_series(ds_list)

        # Advancing one column travels along the row direction by the column spacing
        expected = np.array(origin) + column_spacing * np.array(ROTATED_ORIENTATION[:3])
        assert image.TransformIndexToPhysicalPoint((1, 0, 0)) == pytest.approx(expected)

        # Advancing one row travels along the column direction
        expected = np.array(origin) + column_spacing * np.array(ROTATED_ORIENTATION[3:])
        assert image.TransformIndexToPhysicalPoint((0, 1, 0)) == pytest.approx(expected)

    def test_transposed_direction_would_be_caught(self):
        """Guards the guard: confirm the rotated fixture actually distinguishes a transpose."""
        ds_list = make_ct_series(shape=(3, 8, 8), orientation=ROTATED_ORIENTATION)
        _, _, _, direction = parse_image_series(ds_list)

        assert not np.allclose(direction, direction.T)

    def test_orientation_helper_takes_columns_not_the_flattened_prefix(self):
        matrix = np.array(
            [[0.0, -1.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 0.0, 1.0]]
        )
        orientation = sitk_direction_to_image_orientation_patient(matrix.ravel())

        assert orientation == pytest.approx([0.0, 1.0, 0.0, -1.0, 0.0, 0.0])
        assert orientation != pytest.approx(matrix.ravel()[:6].tolist())


class TestPixelToPatientMatrix:
    """The affine used to place contour points into patient space."""

    def test_matches_sitk_for_an_anisotropic_series(self, anisotropic_ct_series):
        matrix = get_pixel_to_patient_transformation_matrix(anisotropic_ct_series)
        image = SimpleITKImageBuilder().from_image_series(anisotropic_ct_series)

        for index in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (5, 3, 2)]:
            homogeneous = np.array([*index, 1.0])
            from_matrix = (matrix @ homogeneous)[:3]
            from_sitk = image.TransformIndexToPhysicalPoint(index)
            assert from_matrix == pytest.approx(from_sitk, abs=1e-4)

    def test_matches_sitk_for_a_rotated_series(self):
        ds_list = make_ct_series(
            shape=(4, 8, 8), pixel_spacing=(3.0, 1.5), orientation=ROTATED_ORIENTATION
        )
        matrix = get_pixel_to_patient_transformation_matrix(ds_list)
        image = SimpleITKImageBuilder().from_image_series(ds_list)

        for index in [(0, 0, 0), (2, 1, 3), (7, 7, 3)]:
            homogeneous = np.array([*index, 1.0])
            from_matrix = (matrix @ homogeneous)[:3]
            assert from_matrix == pytest.approx(
                image.TransformIndexToPhysicalPoint(index), abs=1e-4
            )

    def test_patient_to_pixel_inverts_pixel_to_patient(self, anisotropic_ct_series):
        forward = get_pixel_to_patient_transformation_matrix(anisotropic_ct_series)
        backward = get_patient_to_pixel_transformation_matrix(anisotropic_ct_series)

        assert (backward @ forward) == pytest.approx(np.identity(4), abs=1e-4)


class TestVolumeParsing:
    def test_rescale_slope_and_intercept_are_applied(self):
        stored = np.full((2, 4, 4), 1200, dtype=np.uint16)
        ds_list = make_ct_series(volume=stored, rescale_intercept=-1024.0)

        volume, _, _, _ = parse_image_series(ds_list)

        assert np.allclose(volume, 1200 - 1024.0)

    def test_origin_is_the_first_slice_image_position(self):
        ds_list = make_ct_series(origin=(-11.0, -22.0, -33.0), slice_spacing=2.5)
        _, origin, spacing, _ = parse_image_series(ds_list)

        assert origin == pytest.approx([-11.0, -22.0, -33.0])
        assert spacing[2] == pytest.approx(2.5)

    def test_volume_axes_are_z_row_column(self):
        ds_list = make_ct_series(shape=(5, 8, 12))
        volume, _, _, _ = parse_image_series(ds_list)

        assert volume.shape == (5, 8, 12)

        image = SimpleITKImageBuilder().from_image_series(ds_list)
        # SimpleITK reports size as (x, y, z), the reverse of the numpy axis order
        assert image.GetSize() == (12, 8, 5)
