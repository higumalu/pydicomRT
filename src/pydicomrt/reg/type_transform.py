import numpy as np
import SimpleITK as sitk
from pydicom import Dataset

from pydicomrt.utils.sitk_transform import sitk_direction_to_image_orientation_patient

def affine_to_homogeneous_matrix(transform: sitk.AffineTransform) -> np.ndarray:
    """
    Convert a SimpleITK affine transform to a 4x4 homogeneous matrix.

    The bridge from a registration result to something DICOM can store.

    Parameters
    ----------
    transform : sitk.AffineTransform
        A 3D linear transform exposing ``GetMatrix()`` and ``GetTranslation()``.
        A non-zero rotation centre is folded into the matrix offset automatically.

    Returns
    -------
    np.ndarray
        4x4, ``float64``, row-major: rotation in ``[:3, :3]``, translation in ``[:3, 3]``,
        bottom row ``[0, 0, 0, 1]``.

    See Also
    --------
    to_centre_free_affine : The equivalent conversion when a SimpleITK affine is needed.
    SpatialRegistrationBuilder.add_registration : Consumes the flattened result.

    Notes
    -----
    Direction is preserved, not flipped. Registration functions return *fixed to moving*
    and Spatial REG stores *moving to fixed*, so for that IOD the usual call is
    ``affine_to_homogeneous_matrix(transform.GetInverse())``.

    Flatten the matrix into 16 row-major values for the builder, for example:
    ``matrix.astype("float32").ravel().tolist()``.
    """
    matrix3x3 = np.array(transform.GetMatrix()).reshape((3, 3))
    translation = np.array(transform.GetTranslation())
    centre = np.array(transform.GetCenter()) if hasattr(transform, "GetCenter") else np.zeros(3)
    hom_mat = np.eye(4)
    hom_mat[:3, :3] = matrix3x3
    hom_mat[:3, 3] = centre + translation - matrix3x3 @ centre
    return hom_mat


def sitk_displacement_field_to_deformable_registration_grid(transform: sitk.DisplacementFieldTransform) -> Dataset:
    """
    Convert a SimpleITK displacement field into a Deformable Registration Grid item.

    Parameters
    ----------
    transform : sitk.DisplacementFieldTransform
        The displacement field. Its geometry -- origin, spacing, direction and size --
        becomes the grid's, so the field defines the region the deformation is described
        over.

    Returns
    -------
    Dataset
        One Deformable Registration Grid Sequence item, carrying GridDimensions,
        GridResolution, ImagePositionPatient, ImageOrientationPatient and VectorGridData.

    See Also
    --------
    DeformableSpatialRegistrationBuilder.add_registration : Uses this internally.
    get_deformable_registrations : Read the result back.

    Notes
    -----
    VectorGridData is written as ``float32`` in DICOM's x-fastest order, the reverse of
    the numpy ``(z, y, x, 3)`` layout the field is held in.
    """
    # 先從 DisplacementFieldTransform 取出對應的影像，再從影像取得空間資訊
    displacement_image = transform.GetDisplacementField()
    if displacement_image.GetDimension() != 3 or displacement_image.GetNumberOfComponentsPerPixel() != 3:
        raise ValueError("DICOM deformation grids require a 3D field with three components")
    axes = np.asarray(displacement_image.GetDirection()).reshape(3, 3)
    if not np.allclose(axes.T @ axes, np.eye(3), atol=1e-4) or not np.isclose(np.linalg.det(axes), 1.0, atol=1e-4):
        raise ValueError("DICOM deformation grid direction must be right-handed and orthonormal")
    displacement_field = sitk.GetArrayFromImage(displacement_image)
    displacement_field = displacement_field.astype("<f4", copy=False)
    vector_grid_data = displacement_field.ravel(order="C").tobytes()

    origin = list(map(float, displacement_image.GetOrigin()))
    spacing = list(map(float, displacement_image.GetSpacing()))
    size = list(map(int, displacement_image.GetSize()))
    # ImageOrientationPatient is the row and column direction, i.e. the first two columns of
    # the direction matrix -- not the first six entries of its row-major flattening.
    orientation = sitk_direction_to_image_orientation_patient(displacement_image.GetDirection())

    deformable_registration_grid = Dataset()
    deformable_registration_grid.ImagePositionPatient = origin
    deformable_registration_grid.ImageOrientationPatient = orientation
    deformable_registration_grid.GridDimensions = size
    deformable_registration_grid.GridResolution = spacing
    deformable_registration_grid.VectorGridData = vector_grid_data

    return deformable_registration_grid
