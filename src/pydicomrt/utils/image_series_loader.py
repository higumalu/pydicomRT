"""
Load DICOM image series from disk, sorted along the slice normal.

Slice order is load-bearing: an unsorted series still produces a valid-looking volume, just
flipped in z, and every contour and registration built on it inherits the flip.
"""

import glob
import os
from pathlib import Path
from typing import List, Union

from pydicom import dcmread
from pydicom.dataset import Dataset

from pydicomrt.utils.coordinate_transform import (
    get_slice_directions,
    get_slice_position,
)

__all__ = [
    "sort_image_series",
    "load_sorted_image_series",
    "get_slice_directions",
    "get_slice_position",
]


def sort_image_series(image_series: List[Dataset]) -> List[Dataset]:
    """
    Sort slices by their position along the slice normal, ascending.

    Use this when the slices came from somewhere other than
    :func:`load_sorted_image_series` -- a PACS query, a database, a list you filtered --
    since every other function in the library assumes this ordering.

    Parameters
    ----------
    image_series : list of Dataset
        Slices of one series, in any order.

    Returns
    -------
    list of Dataset
        A new list, ascending along the slice normal. The input is not modified and the
        ``Dataset`` objects are shared, not copied.

    Raises
    ------
    InvalidImageOrientationError
        If any slice's ImageOrientationPatient is not two orthogonal unit vectors.

    Notes
    -----
    Sorting is by projection onto the slice normal, not by ``InstanceNumber`` and not by
    ``ImagePositionPatient[2]``. Scanners assign instance numbers in acquisition order,
    which need not be spatial order, and the z component of the position is only the
    slice-normal projection for axis-aligned axial series.
    """
    return sorted(image_series, key=get_slice_position)


def load_sorted_image_series(
    image_series_path: Union[str, List[str]],
    ) -> List[Dataset]:
    """
    Read a DICOM image series and return it sorted along the slice normal.

    Parameters
    ----------
    image_series_path : str or list of str
        A directory containing ``*.dcm`` files, or an explicit list of file paths. The
        directory form globs ``*.dcm`` only and does not recurse, so extensionless DICOM
        (common on PACS exports) must be passed as a list.

    Returns
    -------
    list of Dataset
        Slices in ascending slice-normal order, fully read into memory.

    Raises
    ------
    FileNotFoundError
        If the directory holds no ``*.dcm`` files, or a listed path does not exist. Also
        raised when a *file* path is passed where a directory was expected, since that
        matches neither branch and leaves the path list empty.
    ValueError
        If a file cannot be parsed as DICOM, or the series geometry is unusable. The
        message names the offending file.

    See Also
    --------
    sort_image_series : Sorting alone, for slices you already have.
    image_series_to_sitk_image : The usual next call.

    Examples
    --------
    >>> series = load_sorted_image_series("/path/to/CT")            # doctest: +SKIP
    >>> len(series)                                                 # doctest: +SKIP
    120

    Extensionless files, or a subset:

    >>> from pathlib import Path                                    # doctest: +SKIP
    >>> paths = sorted(str(p) for p in Path(folder).iterdir())      # doctest: +SKIP
    >>> series = load_sorted_image_series(paths)                    # doctest: +SKIP

    Notes
    -----
    No filtering by SeriesInstanceUID is performed. A directory holding two series
    produces one interleaved list that will sort into a physically meaningless order.
    """
    dcm_path_list: List[str] = []

    if isinstance(image_series_path, str) and Path(image_series_path).is_dir():
        dcm_path_list = glob.glob(os.path.join(image_series_path, "*.dcm"))
    elif isinstance(image_series_path, (list, tuple)):
        for path in image_series_path:
            if not Path(path).exists():
                raise FileNotFoundError(f"File {path} does not exist")
        dcm_path_list = list(image_series_path)

    if not dcm_path_list:
        raise FileNotFoundError(f"No DICOM files found in {image_series_path}")

    image_series = []
    for dcm_path in dcm_path_list:
        try:
            image_series.append(dcmread(dcm_path))
        except Exception as error:
            raise ValueError(f"Error reading DICOM file {dcm_path}: {error}") from error

    try:
        return sort_image_series(image_series)
    except ValueError as error:
        raise ValueError(
            f"Error sorting DICOM files: {error}. Check ImageOrientationPatient."
        ) from error
