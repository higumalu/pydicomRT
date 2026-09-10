"""
DICOM series loading and slice ordering.

Slice order is load-bearing everywhere downstream: a series read in the wrong order still
produces a valid-looking volume, just mirrored along z, and every contour and registration
built on it inherits the flip.
"""

import numpy as np
import pytest

from synthetic import make_ct_series
from pydicomrt.utils.coordinate_transform import InvalidImageOrientationError
from pydicomrt.utils.image_series_loader import (
    get_slice_directions,
    get_slice_position,
    load_sorted_image_series,
    sort_image_series,
)


@pytest.fixture
def series_on_disk(tmp_path):
    """Write a series to disk in deliberately shuffled filename order."""
    series = make_ct_series(shape=(6, 8, 8), origin=(0.0, 0.0, -10.0), slice_spacing=2.5)
    shuffled = [series[i] for i in [3, 0, 5, 1, 4, 2]]

    for index, ds in enumerate(shuffled):
        ds.save_as(tmp_path / f"slice_{index}.dcm", enforce_file_format=True)

    return tmp_path, series


class TestSorting:
    def test_sorts_by_position_along_the_slice_normal(self):
        series = make_ct_series(shape=(5, 8, 8), origin=(0.0, 0.0, -10.0), slice_spacing=2.0)
        shuffled = [series[i] for i in [2, 4, 0, 3, 1]]

        ordered = sort_image_series(shuffled)

        positions = [get_slice_position(ds) for ds in ordered]
        assert positions == sorted(positions)
        assert [ds.SOPInstanceUID for ds in ordered] == [ds.SOPInstanceUID for ds in series]

    def test_slice_directions_are_orthonormal(self):
        series = make_ct_series(shape=(2, 4, 4))

        row, column, normal = get_slice_directions(series[0])

        assert np.dot(row, column) == pytest.approx(0.0, abs=1e-6)
        assert np.linalg.norm(normal) == pytest.approx(1.0)
        assert normal == pytest.approx(np.cross(row, column))

    def test_non_orthogonal_orientation_is_rejected(self):
        """
        A typed error, not a bare Exception: callers need to be able to distinguish bad
        geometry from a missing file without matching on message text.
        """
        series = make_ct_series(shape=(2, 4, 4))
        series[0].ImageOrientationPatient = [1.0, 0.0, 0.0, 0.9, 0.1, 0.0]

        with pytest.raises(InvalidImageOrientationError, match="orthogonal unit vectors"):
            get_slice_directions(series[0])

    def test_the_orientation_error_is_a_value_error(self):
        assert issubclass(InvalidImageOrientationError, ValueError)


class TestLoading:
    def test_loads_a_directory_in_slice_order(self, series_on_disk):
        directory, expected = series_on_disk

        loaded = load_sorted_image_series(str(directory))

        assert len(loaded) == len(expected)
        assert [ds.SOPInstanceUID for ds in loaded] == [ds.SOPInstanceUID for ds in expected]

    def test_loads_an_explicit_file_list_in_slice_order(self, series_on_disk):
        directory, expected = series_on_disk
        paths = sorted(str(p) for p in directory.glob("*.dcm"))

        loaded = load_sorted_image_series(paths)

        assert [ds.SOPInstanceUID for ds in loaded] == [ds.SOPInstanceUID for ds in expected]

    def test_empty_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No DICOM files"):
            load_sorted_image_series(str(tmp_path))

    def test_missing_file_in_a_list_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist"):
            load_sorted_image_series([str(tmp_path / "nope.dcm")])
