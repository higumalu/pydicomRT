"""
Conversion between DICOM image series and SimpleITK images.

This is where the two libraries' geometry conventions meet, and where a mistake produces a
plausible image in the wrong place rather than an error. The four helpers at the top own
that translation; nothing else in the package should index geometry tags by hand.

Axis ordering, for reference:

===========================  ===========================
numpy volume                 ``(slice, row, column)``
SimpleITK size/spacing       ``(x, y, z)``
DICOM ``PixelSpacing``       ``[row spacing, column spacing]``
===========================  ===========================
"""

from typing import List, Sequence, Tuple

import numpy as np
import SimpleITK as sitk
from pydicom.dataset import Dataset

from pydicomrt.utils.coordinate_transform import (
    get_slice_directions,
    get_spacing_between_slices,
)
from pydicomrt.utils.image_series_loader import sort_image_series

__all__ = [
    "SimpleITKImageBuilder",
    "get_sitk_spacing",
    "get_sitk_direction",
    "sitk_direction_to_image_orientation_patient",
    "sitk_spacing_to_pixel_spacing",
    "parse_image_series",
    "image_series_to_sitk_image",
    "resample_to_reference_image",
]


# --------------------------------------------------------------------------------------- #
# Geometry conventions
# --------------------------------------------------------------------------------------- #

def get_sitk_spacing(image_slice: Dataset, slice_spacing: float) -> np.ndarray:
    """
    Build a SimpleITK spacing triple from a DICOM slice.

    DICOM PixelSpacing (0028,0030) is ordered ``[row spacing, column spacing]``: the first
    value is the distance between adjacent *rows*, i.e. the step taken when the row index
    (SimpleITK's y axis) increases. SimpleITK spacing is ordered ``(x, y, z)``, so the two
    in-plane values must be swapped rather than copied across.

    Parameters
    ----------
    image_slice : Dataset
        Slice carrying PixelSpacing.
    slice_spacing : float
        Distance between slice centres.

    Returns
    -------
    np.ndarray
        Spacing as ``(x, y, z)``.
    """
    row_spacing = float(image_slice.PixelSpacing[0])
    column_spacing = float(image_slice.PixelSpacing[1])
    return np.array([column_spacing, row_spacing, float(slice_spacing)])


def get_sitk_direction(image_slice: Dataset) -> np.ndarray:
    """
    Build a SimpleITK direction matrix from a DICOM slice.

    A SimpleITK direction matrix holds the axis directions in its *columns*: column j is the
    unit vector travelled when index axis j increases. ``get_slice_directions`` returns the
    row, column and slice directions as three separate vectors, so stacking them produces
    those vectors as rows -- the transpose of what SimpleITK expects. Stacking and
    transposing keeps oblique and rotated acquisitions correct; for axis-aligned axial
    series the matrix is symmetric and the distinction is invisible.

    Parameters
    ----------
    image_slice : Dataset
        Slice carrying ImageOrientationPatient.

    Returns
    -------
    np.ndarray
        3x3 direction matrix in SimpleITK (column-vector) convention.
    """
    row_direction, column_direction, slice_direction = get_slice_directions(image_slice)
    return np.array([row_direction, column_direction, slice_direction]).T


def sitk_direction_to_image_orientation_patient(direction: Sequence[float]) -> List[float]:
    """
    Convert a SimpleITK direction matrix into DICOM ImageOrientationPatient (0020,0037).

    ImageOrientationPatient has VM 6: the row direction followed by the column direction.
    Those are the first two *columns* of a SimpleITK direction matrix, which is not the same
    as the first six entries of its row-major flattening -- the two only coincide when the
    matrix is symmetric, as it is for axis-aligned axial series.

    Parameters
    ----------
    direction : sequence
        SimpleITK direction as a flat 9-sequence or a 3x3 array.

    Returns
    -------
    list[float]
        Six floats: row direction then column direction.
    """
    matrix = np.asarray(direction, dtype=float).reshape(3, 3)
    return matrix[:, 0].tolist() + matrix[:, 1].tolist()


def sitk_spacing_to_pixel_spacing(spacing: Sequence[float]) -> List[float]:
    """
    Convert SimpleITK ``(x, y, z)`` spacing into DICOM PixelSpacing (0028,0030).

    Parameters
    ----------
    spacing : sequence of float
        SimpleITK spacing, ``(x, y, z)``. Only the first two entries are used.

    Returns
    -------
    list of float
        ``[row spacing, column spacing]``, i.e. ``(y, x)`` -- the DICOM order.

    See Also
    --------
    get_sitk_spacing : The inverse direction, DICOM to SimpleITK.

    Examples
    --------
    >>> sitk_spacing_to_pixel_spacing((0.5, 1.5, 3.0))
    [1.5, 0.5]
    """
    return [float(spacing[1]), float(spacing[0])]


# --------------------------------------------------------------------------------------- #
# Series -> SimpleITK
# --------------------------------------------------------------------------------------- #

def parse_image_series(
    image_series: List[Dataset],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read volume and geometry out of a DICOM image series.

    Rescale slope and intercept are applied, so the volume is in the modality's real units
    (HU for CT). The series is sorted along the slice normal first.

    Parameters
    ----------
    image_series : list[Dataset]
        Slices of one series.

    Returns
    -------
    volume : np.ndarray
        ``(slice, row, column)``, rescaled.
    origin : np.ndarray
        ImagePositionPatient of the first slice.
    spacing : np.ndarray
        ``(x, y, z)``.
    direction : np.ndarray
        3x3, axis directions in columns.

    Notes
    -----
    Assumes a single, uniformly spaced series; spacing and orientation are taken from the
    first slice.
    """
    image_series = sort_image_series(image_series)
    first_slice = image_series[0]

    volume = np.stack([ds.pixel_array for ds in image_series], axis=0)
    slope = getattr(first_slice, "RescaleSlope", 1)
    intercept = getattr(first_slice, "RescaleIntercept", 0)
    volume = volume * slope + intercept

    slice_spacing = get_spacing_between_slices(image_series)
    origin = np.array([float(v) for v in first_slice.ImagePositionPatient])
    spacing = get_sitk_spacing(first_slice, slice_spacing)
    direction = get_sitk_direction(first_slice)

    return volume, origin, spacing, direction


def image_series_to_sitk_image(image_series: List[Dataset]) -> sitk.Image:
    """
    Convert a DICOM image series straight to a ``sitk.Image``.

    The one-call form of :meth:`SimpleITKImageBuilder.from_image_series`, and the usual
    entry point for reading images.

    Parameters
    ----------
    image_series : list of Dataset
        Slices of one series, in any order. They are sorted along the slice normal.

    Returns
    -------
    sitk.Image
        The volume with origin, spacing and direction set, in the modality's real units
        (HU for CT -- rescale slope and intercept are applied).

    See Also
    --------
    parse_image_series : Same reading, returning the four pieces separately.
    SimpleITKImageBuilder : Assemble an image from a volume of your own.

    Examples
    --------
    >>> series = load_sorted_image_series("/path/to/CT")   # doctest: +SKIP
    >>> image = image_series_to_sitk_image(series)         # doctest: +SKIP
    >>> image.GetSize()                                    # doctest: +SKIP
    (512, 512, 120)
    """
    return SimpleITKImageBuilder().from_image_series(image_series)


class SimpleITKImageBuilder:
    """
    Assemble a ``sitk.Image`` from a volume plus geometry.

    Either set the four pieces individually and call :meth:`build`, or use one of the
    ``from_*`` shortcuts, which build and return a ``sitk.Image`` in one step.

    The setters return ``self`` so calls chain; the ``from_*`` shortcuts do not, because
    they return the finished image.

    Methods
    -------
    set_volume(volume)
        The voxels, ``(slice, row, column)``.
    set_origin(origin)
        Patient coordinates of voxel ``[0, 0, 0]``.
    set_spacing(spacing)
        Voxel size as ``(x, y, z)``.
    set_direction(direction)
        3x3 with axis directions in *columns*.
    build()
        Assemble and return the image.
    from_image_series(image_series)
        Shortcut: build from DICOM slices.
    from_reference_image(volume, reference_image)
        Shortcut: build from a volume, borrowing geometry.
    from_dicom_directory(directory)
        Shortcut: build via SimpleITK's GDCM reader.

    See Also
    --------
    image_series_to_sitk_image : The common case, as a single function call.

    Examples
    --------
    Put a volume of your own onto an existing series' geometry:

    >>> volume, origin, spacing, direction = parse_image_series(series)   # doctest: +SKIP
    >>> image = (                                                        # doctest: +SKIP
    ...     SimpleITKImageBuilder()
    ...     .set_volume(my_segmentation)
    ...     .set_origin(origin)
    ...     .set_spacing(spacing)
    ...     .set_direction(direction)
    ...     .build()
    ... )
    """

    def __init__(self):
        self.volume = None
        self.origin = None
        self.spacing = None
        self.direction = None
        self.image = None

    def set_volume(self, volume: np.ndarray) -> "SimpleITKImageBuilder":
        """
        Set the voxel data.

        Parameters
        ----------
        volume : np.ndarray
            3D, ordered ``(slice, row, column)`` -- the order a series stacks to, and the
            reverse of what ``sitk.Image.GetSize()`` reports.

        Returns
        -------
        SimpleITKImageBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``volume`` is not 3D.
        """
        if np.ndim(volume) != 3:
            raise ValueError(f"volume must be 3D, got shape {np.shape(volume)}")
        self.volume = volume
        return self

    def set_origin(self, origin: Sequence[float]) -> "SimpleITKImageBuilder":
        """
        Set the patient coordinates of voxel ``[0, 0, 0]``.

        Parameters
        ----------
        origin : sequence of float
            Three values in mm. For a DICOM series this is ImagePositionPatient of the
            first slice *after sorting along the slice normal*.

        Returns
        -------
        SimpleITKImageBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``origin`` is not a 3-element vector.
        """
        origin = np.asarray(origin, dtype=float)
        if origin.shape != (3,):
            raise ValueError(f"origin must be a 3-element vector, got shape {origin.shape}")
        self.origin = origin.tolist()
        return self

    def set_spacing(self, spacing: Sequence[float]) -> "SimpleITKImageBuilder":
        """
        Set the voxel size.

        Parameters
        ----------
        spacing : sequence of float
            ``(x, y, z)`` in mm -- SimpleITK's order. DICOM PixelSpacing is
            ``[row, column]``, i.e. ``(y, x)``, so the in-plane pair must be swapped, not
            copied across. Use :func:`get_sitk_spacing` rather than doing it by hand.

        Returns
        -------
        SimpleITKImageBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``spacing`` is not a 3-element vector.
        """
        spacing = np.asarray(spacing, dtype=float)
        if spacing.shape != (3,):
            raise ValueError(f"spacing must be a 3-element vector, got shape {spacing.shape}")
        self.spacing = spacing.tolist()
        return self

    def set_direction(self, direction: Sequence[float]) -> "SimpleITKImageBuilder":
        """
        Set the axis directions.

        Parameters
        ----------
        direction : sequence of float
            3x3 array or flat 9-sequence. The axis vectors go in the **columns**: column
            ``j`` is the unit vector travelled when index axis ``j`` increases. Building it
            with the vectors as rows gives the transpose, which is identical for
            axis-aligned axial series and mirrored for oblique ones. Use
            :func:`get_sitk_direction`.

        Returns
        -------
        SimpleITKImageBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``direction`` is not 3x3 or a flat 9-sequence.
        """
        direction = np.asarray(direction, dtype=float)
        if direction.shape not in {(3, 3), (9,)}:
            raise ValueError(
                f"direction must be 3x3 or a flat 9-sequence, got shape {direction.shape}"
            )
        self.direction = direction.reshape(3, 3)
        return self

    def from_image_series(self, image_series: List[Dataset]) -> sitk.Image:
        """
        Build directly from a DICOM image series.

        Parameters
        ----------
        image_series : list of Dataset
            Slices of one series, in any order.

        Returns
        -------
        sitk.Image
            The finished image -- this shortcut calls :meth:`build` for you, so it does
            not return ``self``.
        """
        volume, origin, spacing, direction = parse_image_series(image_series)
        return (
            self.set_volume(volume)
            .set_origin(origin)
            .set_spacing(spacing)
            .set_direction(direction)
            .build()
        )

    def from_reference_image(self, volume: np.ndarray, reference_image: sitk.Image) -> sitk.Image:
        """
        Build from a volume, borrowing geometry from an existing image.

        The usual way to wrap a numpy result -- a segmentation, a dose grid, a filtered
        volume -- so it sits on the same grid as the image it was computed from.

        Parameters
        ----------
        volume : np.ndarray
            3D, ``(slice, row, column)``. Must match ``reference_image`` in shape, since
            it inherits that image's geometry without being resampled.
        reference_image : sitk.Image
            Supplies origin, spacing and direction.

        Returns
        -------
        sitk.Image
            The finished image; this shortcut does not return ``self``.

        See Also
        --------
        resample_to_reference_image : When the volume is on a *different* grid and needs
            resampling rather than relabelling.
        """
        return (
            self.set_volume(volume)
            .set_origin(reference_image.GetOrigin())
            .set_spacing(reference_image.GetSpacing())
            .set_direction(reference_image.GetDirection())
            .build()
        )

    def from_dicom_directory(self, directory: str) -> sitk.Image:
        """
        Build via SimpleITK's own GDCM series reader.

        Reads the directory with ITK's reader rather than this library's series loader.
        Useful as an independent cross-check, and when the directory holds a series this
        library's loader rejects.

        Parameters
        ----------
        directory : str
            Directory holding one DICOM series. If it holds several, GDCM picks one.

        Returns
        -------
        sitk.Image
            The finished image; this shortcut does not return ``self``.

        Notes
        -----
        Unlike :meth:`from_image_series` this path gives you no access to the underlying
        ``Dataset`` objects, so ROI and registration builders -- which need patient and
        frame-of-reference context -- cannot be driven from it.
        """
        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(directory))
        image = reader.Execute()
        return self.from_reference_image(sitk.GetArrayFromImage(image), image)

    def build(self) -> sitk.Image:
        """
        Assemble the image.

        Returns
        -------
        sitk.Image
            Volume with origin, spacing and direction applied.

        Raises
        ------
        ValueError
            If volume, origin, spacing or direction has not been set. The message names
            the missing pieces.
        """
        missing = [
            name for name in ("volume", "origin", "spacing", "direction")
            if getattr(self, name) is None
        ]
        if missing:
            raise ValueError(f"cannot build: {', '.join(missing)} not set")

        self.image = sitk.GetImageFromArray(self.volume)
        self.image.SetOrigin(self.origin)
        self.image.SetSpacing(self.spacing)
        self.image.SetDirection(self.direction.ravel().tolist())
        return self.image


# --------------------------------------------------------------------------------------- #
# Resampling
# --------------------------------------------------------------------------------------- #

def resample_to_reference_image(
    reference_image: sitk.Image,
    source_image: sitk.Image,
    default_value: float = -1000.0,
    interpolator: int = sitk.sitkLinear,
    ) -> sitk.Image:
    """
    Resample the source image onto the reference image grid (identity transform).

    Parameters
    ----------
    reference_image : sitk.Image
        Defines the output grid and pixel type.
    source_image : sitk.Image
        Image to resample.
    default_value : float
        Padding for samples outside the source extent. -1000 is CT air; pass 0 for MR/PET
        or already-clipped images, otherwise the padding invents a soft-tissue shell.
    interpolator : int
        SimpleITK interpolator. Use ``sitkNearestNeighbor`` for masks and labels.

    Returns
    -------
    sitk.Image
        Source image on the reference grid.
    """
    return sitk.Resample(
        source_image,
        reference_image,
        sitk.Transform(),
        interpolator,
        default_value,
        reference_image.GetPixelID(),
    )
