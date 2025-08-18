import numpy as np
import SimpleITK as sitk

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


def displacement_field_to_dict(transform: sitk.DisplacementFieldTransform) -> dict:
    """
    Convert a SimpleITK DisplacementFieldTransform to a dictionary.
    Args:
        transform (sitk.DisplacementFieldTransform): The displacement field transform to convert.
    Returns:
        dict: Dictionary containing displacement field information.
    """
    displacement_field = sitk.GetArrayFromImage(transform.GetDisplacementField())
    origin = transform.GetOrigin()
    spacing = transform.GetSpacing()
    direction = transform.GetDirection()
    size = transform.GetSize()
    return {
        "vector_grid": displacement_field,
        "origin": origin,
        "spacing": spacing,
        "direction": direction,
        "size": size,
    }
