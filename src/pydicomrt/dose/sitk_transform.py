
import logging
import numpy as np
import SimpleITK as sitk

logger = logging.getLogger(__name__)

def get_dose_spacing(dose_ds):
    """
    Voxel size of an RT Dose grid, in SimpleITK order.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object carrying PixelSpacing and GridFrameOffsetVector.

    Returns
    -------
    list of float
        ``(x, y, z)`` in mm. The in-plane pair is swapped out of DICOM's
        ``[row, column]`` order; ``z`` is the mean step in GridFrameOffsetVector.

    Warns
    -----
    Logs a warning on the ``pydicomrt.dose.sitk_transform`` logger when the frame offsets
    are not uniformly spaced (standard deviation above 0.001 mm). A dose grid is allowed
    to be non-uniform in z; SimpleITK images are not, so the mean is used and the
    departure is worth knowing about.

    See Also
    --------
    get_dose_image : Assembles this with the array, origin and direction.
    """
    # DICOM PixelSpacing is [row spacing, column spacing]; SimpleITK wants (x, y, z), and
    # the column spacing is the one that applies to x.
    row_spacing, column_spacing = dose_ds.PixelSpacing
    x_spacing, y_spacing = column_spacing, row_spacing
    frame_offset_vector_list = dose_ds.GridFrameOffsetVector
    frame_distance_list = [np.linalg.norm(frame_offset_vector_list[i + 1] - frame_offset_vector_list[i])
                           for i in range(len(frame_offset_vector_list) - 1)]
    z_spacing = np.mean(frame_distance_list)
    std_distance = np.std(frame_distance_list)
    if std_distance > 0.001:
        logger.warning("Dose grid frame offsets are not uniformly spaced (std %.4f mm)", std_distance)
    return [x_spacing, y_spacing, z_spacing]

def get_dose_direction(dose_ds):
    """
    Axis directions of an RT Dose grid, in SimpleITK convention.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object carrying ImageOrientationPatient.

    Returns
    -------
    np.ndarray
        3x3 with the axis vectors in **columns**. The slice direction is the cross product
        of the two stored vectors.

    Notes
    -----
    No orthonormality check is made here, unlike :func:`get_slice_directions` for images.
    """
    ori = dose_ds.ImageOrientationPatient
    row = ori[0:3]
    col = ori[3:6]
    cross = np.cross(row, col)
    dose_direction = np.asarray([row, col, cross]).T
    return dose_direction

def get_dose_origin(dose_ds):
    """
    Patient coordinates of the dose grid's first voxel.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object.

    Returns
    -------
    list of float
        ImagePositionPatient, ``[x, y, z]`` in mm. Returned as pydicom's ``DSfloat``
        values; wrap in ``[float(v) for v in ...]`` if you need plain floats.
    """
    return dose_ds.ImagePositionPatient

def get_dose_array(dose_ds):
    """
    Dose values as a numpy array, scaling applied.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object.

    Returns
    -------
    np.ndarray
        ``(frame, row, column)`` in DoseUnits -- Gy for a conventional plan dose. Stored
        dose is scaled unsigned integers; multiplying by DoseGridScaling is what makes the
        numbers physical, so never read ``dose_ds.pixel_array`` directly.

    See Also
    --------
    get_dose_image : The same values with geometry attached.
    """
    return dose_ds.pixel_array * dose_ds.DoseGridScaling

def get_dose_image(dose_ds):
    """
    Read an RT Dose object into a ``sitk.Image``.

    The read path matching :class:`RTDoseBuilder`: scaling applied, geometry attached.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object.

    Returns
    -------
    sitk.Image
        Dose in DoseUnits (Gy for a conventional plan dose), on the dose grid -- which is
        typically coarser than the CT and need not cover the same extent. Use
        :func:`resample_to_reference_image` with ``default_value=0.0`` to put it on the CT
        grid; the -1000 default would invent negative dose.

    See Also
    --------
    RTDoseBuilder : The write path.
    get_dose_array : Values without geometry.

    Examples
    --------
    >>> dose_image = get_dose_image(dcmread("rtdose.dcm"))          # doctest: +SKIP
    >>> on_ct = resample_to_reference_image(ct_image, dose_image, default_value=0.0)  # doctest: +SKIP
    """
    dose_array = get_dose_array(dose_ds)
    spacing = get_dose_spacing(dose_ds)
    origin = get_dose_origin(dose_ds)
    direction = get_dose_direction(dose_ds)
    sitk_image = sitk.GetImageFromArray(dose_array)
    sitk_image.SetSpacing(spacing)
    sitk_image.SetOrigin(origin)
    sitk_image.SetDirection(direction.flatten())
    return sitk_image
