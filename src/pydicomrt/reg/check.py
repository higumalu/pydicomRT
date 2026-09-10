"""IOD validation for the two registration objects."""

from pydicom.dataset import Dataset
import numpy as np

from pydicomrt.utils.validate_dcm_info import build_check_result, check_iod

from .iod import DEFORMABLE_SPATIAL_REGISTRATION_IOD, SPATIAL_REGISTRATION_IOD


def _matrix_errors(items):
    from .builder import validate_transformation_matrix

    errors = []
    for item in items:
        try:
            validate_transformation_matrix(item.FrameOfReferenceTransformationMatrix,
                                           item.FrameOfReferenceTransformationMatrixType)
        except (AttributeError, TypeError, ValueError) as exc:
            errors.append(f"Invalid transformation matrix: {exc}")
    return errors


def check_spatial_reg_iod(reg_ds: Dataset) -> dict:
    """
    Validate a Spatial Registration object against its CIOD.

    Parameters
    ----------
    reg_ds : Dataset
        Spatial Registration object.

    Returns
    -------
    dict
        ``{"result": bool, "content": list of str}``. ``result`` is True when conformant
        and ``content`` empty; otherwise one message per missing or empty required
        element.

    See Also
    --------
    SpatialRegistrationBuilder : Produces conformant objects.
    check_deformable_reg_iod : The deformable equivalent.

    Notes
    -----
    Structural conformance only. A registration can pass here and still point the wrong
    way -- the matrix direction is not something a field-presence check can see.
    """
    errors = check_iod(ds=reg_ds, config_map=SPATIAL_REGISTRATION_IOD)
    if str(getattr(reg_ds, "SOPClassUID", "")) != "1.2.840.10008.5.1.4.1.1.66.1":
        errors.append("SOPClassUID must identify Spatial Registration")
    for reg in getattr(reg_ds, "RegistrationSequence", []):
        if not getattr(reg, "FrameOfReferenceUID", None) and not getattr(reg, "ReferencedImageSequence", None):
            errors.append("Registration item requires FrameOfReferenceUID or ReferencedImageSequence")
        for matrix_reg in getattr(reg, "MatrixRegistrationSequence", []):
            errors.extend(_matrix_errors(getattr(matrix_reg, "MatrixSequence", [])))
    return build_check_result(errors)


def check_deformable_reg_iod(reg_ds: Dataset) -> dict:
    """
    Validate a Deformable Spatial Registration object against its CIOD.

    Parameters
    ----------
    reg_ds : Dataset
        Deformable Spatial Registration object.

    Returns
    -------
    dict
        ``{"result": bool, "content": list of str}``, as for
        :func:`check_spatial_reg_iod`.

    See Also
    --------
    DeformableSpatialRegistrationBuilder : Produces conformant objects.
    check_spatial_reg_iod : The rigid/affine equivalent.

    Notes
    -----
    Structural conformance only; the displacement field's values are not examined.
    """
    errors = check_iod(ds=reg_ds, config_map=DEFORMABLE_SPATIAL_REGISTRATION_IOD)
    if str(getattr(reg_ds, "SOPClassUID", "")) != "1.2.840.10008.5.1.4.1.1.66.3":
        errors.append("SOPClassUID must identify Deformable Spatial Registration")
    grids = []
    for reg in getattr(reg_ds, "DeformableRegistrationSequence", []):
        for key in ("PreDeformationMatrixRegistrationSequence", "PostDeformationMatrixRegistrationSequence"):
            errors.extend(_matrix_errors(getattr(reg, key, [])))
        grids.extend(getattr(reg, "DeformableRegistrationGridSequence", []))
    if not grids:
        errors.append("At least one DeformableRegistrationGridSequence item is required")
    for grid in grids:
        try:
            dims = np.asarray(grid.GridDimensions, dtype=float)
            spacing = np.asarray(grid.GridResolution, dtype=float)
            origin = np.asarray(grid.ImagePositionPatient, dtype=float)
            orientation = np.asarray(grid.ImageOrientationPatient, dtype=float).reshape(2, 3)
            if dims.shape != (3,) or not np.isfinite(dims).all() or np.any(dims <= 0) or np.any(dims != np.floor(dims)):
                raise ValueError("GridDimensions must be three positive integers")
            if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0):
                raise ValueError("GridResolution must be three positive finite values")
            if origin.shape != (3,) or not np.isfinite(origin).all():
                raise ValueError("ImagePositionPatient must contain three finite values")
            if not np.allclose(orientation @ orientation.T, np.eye(2), atol=1e-4):
                raise ValueError("ImageOrientationPatient must contain orthonormal directions")
            if len(grid.VectorGridData) != int(np.prod(dims)) * 12:
                raise ValueError("VectorGridData byte length does not match GridDimensions")
        except (AttributeError, TypeError, ValueError) as exc:
            errors.append(f"Invalid deformation grid: {exc}")
    return build_check_result(errors)
