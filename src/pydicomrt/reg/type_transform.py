import numpy as np
import SimpleITK as sitk
from pydicom import Dataset

def affine_to_homogeneous_matrix(transform: sitk.AffineTransform) -> np.ndarray:
    """
    Convert a SimpleITK AffineTransform to a 4x4 homogeneous transformation matrix.

    Args:
        transform (sitk.AffineTransform): The affine transform to convert.

    Returns:
        np.ndarray: 4x4 homogeneous matrix (dtype float64).
    """
    matrix3x3 = np.array(transform.GetMatrix()).reshape((3, 3))
    translation = np.array(transform.GetTranslation())
    hom_mat = np.eye(4)
    hom_mat[:3, :3] = matrix3x3
    hom_mat[:3, 3] = translation
    return hom_mat


def sitk_displacement_field_to_deformable_registration_grid(transform: sitk.DisplacementFieldTransform) -> Dataset:
    """
    Convert a SimpleITK DisplacementFieldTransform to a dictionary.
    Args:
        transform (sitk.DisplacementFieldTransform): The displacement field transform to convert.
    Returns:
        Dataset: Dataset containing deformable registration grid information.
    """
    # 先從 DisplacementFieldTransform 取出對應的影像，再從影像取得空間資訊
    displacement_image = transform.GetDisplacementField()
    displacement_field = sitk.GetArrayFromImage(displacement_image)
    displacement_field = displacement_field.astype("<f4", copy=False)
    vector_grid_data = displacement_field.ravel(order="C").tobytes()

    origin = list(map(float, displacement_image.GetOrigin()))
    spacing = list(map(float, displacement_image.GetSpacing()))
    size = list(map(int, displacement_image.GetSize()))
    direction = displacement_image.GetDirection()
    orientation = list(map(float, direction[:6]))

    deformable_registration_grid = Dataset()
    deformable_registration_grid.ImagePositionPatient = origin
    deformable_registration_grid.ImageOrientationPatient = orientation
    deformable_registration_grid.GridDimensions = size
    deformable_registration_grid.GridResolution = spacing
    deformable_registration_grid.VectorGridData = vector_grid_data

    return deformable_registration_grid
