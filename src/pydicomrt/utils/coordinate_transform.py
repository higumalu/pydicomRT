"""
Slice geometry and pixel <-> patient coordinate transforms.

This module is the single home for reading geometry out of a DICOM image series.
``image_series_loader`` and ``sitk_transform`` import from here rather than carrying their
own copies, which is how the three used to drift apart.
"""

import warnings
from typing import List, Tuple

import numpy as np
from pydicom.dataset import Dataset


class InvalidImageOrientationError(ValueError):
    """
    ImageOrientationPatient (0020,0037) is not two orthogonal unit vectors.

    A subclass of ``ValueError``, so ``except ValueError`` still catches it. Raised by
    :func:`get_slice_directions` and therefore by everything that reads series geometry:
    loading, sorting, contour placement and SimpleITK conversion.

    In practice this means the file is corrupt or non-conformant. The tolerance is 1e-3 on
    both the dot product and the normal's length, which is loose enough for the rounding
    real scanners write and tight enough to catch a genuinely bad orientation.
    """


def get_slice_directions(image_slice: Dataset) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Split ImageOrientationPatient into its three unit direction vectors.

    Parameters
    ----------
    image_slice : Dataset
        Any slice of the series; orientation is shared across a series.

    Returns
    -------
    row_direction : np.ndarray
        Travelled when the *column* index increases -- IOP[0:3].
    column_direction : np.ndarray
        Travelled when the *row* index increases -- IOP[3:6].
    slice_direction : np.ndarray
        Their cross product, i.e. the slice normal.

    Raises
    ------
    InvalidImageOrientationError
        If the two stored vectors are not orthogonal unit vectors.
    """
    orientation = image_slice.ImageOrientationPatient
    row_direction = np.array(orientation[:3], dtype=float)
    column_direction = np.array(orientation[3:], dtype=float)
    slice_direction = np.cross(row_direction, column_direction)

    if not np.allclose(np.dot(row_direction, column_direction), 0.0, atol=1e-3) or \
            not np.allclose(np.linalg.norm(slice_direction), 1.0, atol=1e-3):
        raise InvalidImageOrientationError(
            f"ImageOrientationPatient {list(orientation)} is not two orthogonal unit vectors"
        )

    return row_direction, column_direction, slice_direction


def get_slice_position(image_slice: Dataset) -> float:
    """
    Position of a slice along the slice normal.

    This is the value :func:`sort_image_series` orders by. It is not the same as
    ``ImagePositionPatient[2]`` unless the series is axial and axis-aligned, and it is not
    ``InstanceNumber``, which scanners are free to assign in any order.

    Parameters
    ----------
    image_slice : Dataset
        Slice carrying ImagePositionPatient and ImageOrientationPatient.

    Returns
    -------
    float
        Signed distance in mm from the patient origin, projected onto the slice normal.

    Raises
    ------
    InvalidImageOrientationError
        If ImageOrientationPatient is not two orthogonal unit vectors.
    """
    _, _, slice_direction = get_slice_directions(image_slice)
    return float(np.dot(slice_direction, image_slice.ImagePositionPatient))


class NonUniformSliceSpacingWarning(UserWarning):
    """
    The slices of a series are not evenly spaced along the slice normal.

    Issued by :func:`get_spacing_between_slices`, and therefore by everything that turns a
    series into geometry: the pixel/patient matrices and SimpleITK conversion. The library
    places slice ``i`` at ``first + i * mean spacing``, so on such a series -- a missing or
    duplicated slice, or an acquisition whose spacing changes part-way -- slices are put
    where they are not, by up to half a slice for one missing slice and by centimetres for
    mixed spacings. Nothing else fails; only physical positions are wrong.

    A gap counts as irregular when it differs from the series' median gap by more than
    1 % or 0.01 mm, whichever is larger, which passes the rounding scanners write.

    To get a correctly placed image, read with
    ``image_series_to_sitk_image(series, slice_spacing="resample")``. To refuse such series
    outright, turn the warning into an error:
    ``warnings.simplefilter("error", NonUniformSliceSpacingWarning)``.
    """


def _slice_positions(image_series: List[Dataset]) -> np.ndarray:
    """Position of every slice along the slice normal, in series order."""
    return np.array([get_slice_position(image_slice) for image_slice in image_series])


def _gap_tolerance(nominal: float) -> float:
    """How far a gap may be from the nominal one before it counts as irregular (mm)."""
    return max(0.01, 0.01 * abs(nominal))


def _irregular_gaps(gaps: np.ndarray) -> Tuple[float, np.ndarray]:
    """The nominal (median) gap, and the indices of the gaps that are not it."""
    nominal = float(np.median(gaps))
    return nominal, np.flatnonzero(np.abs(gaps - nominal) > _gap_tolerance(nominal))


def _describe_irregular_gaps(gaps: np.ndarray, nominal: float, irregular: np.ndarray) -> str:
    tolerance = _gap_tolerance(nominal)
    described = []
    for index in irregular[:3]:
        gap = float(gaps[index])
        multiple = gap / nominal if nominal else 0.0
        if abs(gap) <= tolerance:
            kind = "duplicated position"
        elif round(multiple) >= 2 and abs(gap - round(multiple) * nominal) <= tolerance * round(multiple):
            missing = round(multiple) - 1
            kind = f"{missing} slice{'s' if missing > 1 else ''} missing"
        else:
            kind = "different spacing"
        described.append(f"slices {index} and {index + 1} are {gap:.3f} mm apart ({kind})")
    if len(irregular) > 3:
        described.append(f"and {len(irregular) - 3} more")
    return (
        f"{len(irregular)} of {len(gaps)} slice gaps differ from the usual {nominal:.3f} mm: "
        + "; ".join(described)
    )


def _warn_if_not_uniform(positions: np.ndarray, assumed: str, stacklevel: int = 2) -> None:
    """Issue NonUniformSliceSpacingWarning if ``positions`` are unevenly spaced; ``assumed`` names the spacing used."""
    gaps = np.diff(positions)
    if not len(gaps):
        return
    nominal, irregular = _irregular_gaps(gaps)
    if len(irregular):
        warnings.warn(
            f"slice spacing is not uniform, and slices are placed as if it were {assumed}: "
            f"{_describe_irregular_gaps(gaps, nominal, irregular)}. "
            "Read the image with slice_spacing=\"resample\" to place it correctly",
            NonUniformSliceSpacingWarning,
            stacklevel=stacklevel,
        )


def _uniform_slice_plan(positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    How to turn slices at ``positions`` (ascending) into an evenly spaced stack.

    Returns the indices of the slices to keep -- one per position, so duplicates are
    dropped -- the positions of the output slices, and their spacing. The spacing is the
    smallest gap in the series, so a missing slice is filled in and a series that mixes
    spacings keeps the resolution of its finest part. The output starts on the first slice
    and stops at or before the last. A series that is already even, duplicates aside, comes
    back as its own positions.
    """
    if len(positions) < 2:
        return np.arange(len(positions)), positions, 1.0
    nominal, _ = _irregular_gaps(np.diff(positions))
    tolerance = _gap_tolerance(nominal)
    keep = [0]
    for index in range(1, len(positions)):
        if positions[index] - positions[keep[-1]] > tolerance:
            keep.append(index)
    keep = np.array(keep)
    kept = positions[keep]
    if len(kept) < 2:
        return keep, kept, 1.0
    gaps = np.diff(kept)
    nominal, irregular = _irregular_gaps(gaps)
    span = float(kept[-1] - kept[0])
    if not len(irregular):
        return keep, kept, span / (len(kept) - 1)
    # Gaps within tolerance of the usual one are the usual one, not a finer spacing.
    spacing = nominal if gaps.min() >= nominal - _gap_tolerance(nominal) else float(gaps.min())
    count = int(np.floor(span / spacing + 1e-6)) + 1
    return keep, kept[0] + spacing * np.arange(count), spacing


def get_spacing_between_slices(image_series: List[Dataset]) -> float:
    """
    Mean distance between slice centres, along the slice normal.

    Measured from the first and last slice rather than read from SliceThickness or
    SpacingBetweenSlices, both of which are routinely absent or wrong.

    Parameters
    ----------
    image_series : list of Dataset
        Slices **sorted along the slice normal**. On an unsorted series this returns a
        value that is too small, or negative.

    Returns
    -------
    float
        Mean centre-to-centre distance in mm. Negative if the series runs against the
        slice normal. Exactly ``1.0`` for a single-slice series, which has no spacing to
        measure -- the fallback keeps the transformation matrices invertible.

    Warns
    -----
    NonUniformSliceSpacingWarning
        If the gaps between neighbouring slices are not all the same, beyond rounding. The
        mean is still returned.

    Notes
    -----
    Averaging over the whole series rather than differencing neighbours means a series
    with one duplicated or missing slice still yields a usable spacing, at the cost of
    placing its slices slightly wrong. Uniform spacing is assumed throughout the library;
    the warning says when that assumption does not hold.
    """
    if len(image_series) > 1:
        positions = _slice_positions(image_series)
        mean = float(positions[-1] - positions[0]) / (len(positions) - 1)
        _warn_if_not_uniform(positions, f"the mean, {mean:.3f} mm", stacklevel=3)
        return mean

    return 1.0


def get_pixel_to_patient_transformation_matrix(image_series: List[Dataset]) -> np.ndarray:
    """
    4x4 affine mapping pixel indices to patient coordinates (mm).

    Points are ordered ``(column index, row index, slice index)``, which is the reverse of
    DICOM's PixelSpacing ordering -- see the spacing pairing below.

    Parameters
    ----------
    image_series : list[Dataset]
        Slices sorted along the slice normal.

    Returns
    -------
    np.ndarray
        4x4 matrix; multiply column vectors ``[x_index, y_index, z_index, 1]``.
    """
    first_slice = image_series[0]

    offset = np.array(first_slice.ImagePositionPatient, dtype=float)
    # DICOM PixelSpacing is [row spacing, column spacing]: row_spacing is the step taken
    # when the *row* index advances, which travels along column_direction. Pixel points are
    # ordered (column index, row index, slice index), so the spacings cross over.
    row_spacing, column_spacing = first_slice.PixelSpacing
    slice_spacing = get_spacing_between_slices(image_series)
    row_direction, column_direction, slice_direction = get_slice_directions(first_slice)

    matrix = np.identity(4, dtype=np.float32)
    matrix[:3, 0] = row_direction * float(column_spacing)
    matrix[:3, 1] = column_direction * float(row_spacing)
    matrix[:3, 2] = slice_direction * slice_spacing
    matrix[:3, 3] = offset

    return matrix


def _slice_pixel_to_patient_matrix(image_slice: Dataset) -> np.ndarray:
    """
    4x4 affine mapping ``(column index, row index, 0)`` of one slice to patient mm.

    DICOM's own per-image equation (PS3.3 C.7.6.2.1.1), from this slice's position,
    orientation and pixel spacing alone -- nothing is assumed about the rest of the series,
    so the result lies in the slice's plane whatever the spacing, tilt or shear between
    slices. :func:`get_pixel_to_patient_transformation_matrix` is the series-wide form, which
    places slice ``i`` by the mean spacing along the normal instead.
    """
    row_direction, column_direction, _ = get_slice_directions(image_slice)
    row_spacing, column_spacing = image_slice.PixelSpacing
    matrix = np.identity(4)
    matrix[:3, 0] = row_direction * float(column_spacing)
    matrix[:3, 1] = column_direction * float(row_spacing)
    matrix[:3, 3] = [float(value) for value in image_slice.ImagePositionPatient]
    return matrix


def get_patient_to_pixel_transformation_matrix(image_series: List[Dataset]) -> np.ndarray:
    """
    4x4 affine mapping patient coordinates (mm) to pixel indices.

    The inverse of :func:`get_pixel_to_patient_transformation_matrix`, computed in closed
    form rather than by inverting that matrix numerically.

    Parameters
    ----------
    image_series : list of Dataset
        Slices sorted along the slice normal.

    Returns
    -------
    np.ndarray
        4x4. Multiply column vectors ``[x_mm, y_mm, z_mm, 1]``; the result is ordered
        ``(column index, row index, slice index)`` and is generally not integral -- round
        or floor according to what you are doing.

    See Also
    --------
    apply_transformation_to_3d_points : Apply this to an ``(n, 3)`` point array.
    calc_image_series_affine_mapping : What RTSTRUCT reading needs -- this matrix plus the
        volume shape.
    """
    first_slice = image_series[0]

    offset = np.array(first_slice.ImagePositionPatient, dtype=float)
    row_spacing, column_spacing = first_slice.PixelSpacing
    slice_spacing = get_spacing_between_slices(image_series)
    row_direction, column_direction, slice_direction = get_slice_directions(first_slice)

    # See get_pixel_to_patient_transformation_matrix for the row/column spacing pairing.
    linear = np.identity(3, dtype=np.float32)
    linear[0, :3] = row_direction / float(column_spacing)
    linear[1, :3] = column_direction / float(row_spacing)
    linear[2, :3] = slice_direction / slice_spacing

    matrix = np.identity(4, dtype=np.float32)
    matrix[:3, :3] = linear
    matrix[:3, 3] = offset.dot(-linear.T)

    return matrix


def apply_transformation_to_3d_points(
    points: np.ndarray,
    transformation_matrix: np.ndarray,
    ) -> np.ndarray:
    """
    Apply a 4x4 affine to an ``(n, 3)`` array of points.

    Each point is augmented with a 1 so the translation column applies, then the extra
    coordinate is dropped again.

    Parameters
    ----------
    points : np.ndarray
        ``(n, 3)``. A single point must still be shaped ``(1, 3)``.
    transformation_matrix : np.ndarray
        4x4 affine, such as either of the transformation matrices in this module.

    Returns
    -------
    np.ndarray
        ``(n, 3)``, transformed.

    Examples
    --------
    >>> import numpy as np
    >>> shift = np.identity(4); shift[:3, 3] = [10.0, 0.0, 0.0]
    >>> apply_transformation_to_3d_points(np.array([[1.0, 2.0, 3.0]]), shift)
    array([[11.,  2.,  3.]])
    """
    vec = np.concatenate((points, np.ones((points.shape[0], 1))), axis=1)
    return vec.dot(transformation_matrix.T)[:, :3]
