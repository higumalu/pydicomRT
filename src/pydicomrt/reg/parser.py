
import logging
import numpy as np
import SimpleITK as sitk

from pydicom.dataset import Dataset

logger = logging.getLogger(__name__)

def get_spatial_registrations(reg_ds: Dataset, calc_inverse_matrix: bool = False) -> dict:
    """
    Read the matrices out of a Spatial Registration object.

    Parameters
    ----------
    reg_ds : Dataset
        Spatial Registration object.
    calc_inverse_matrix : bool, optional
        Also populate the reverse direction, by numerically inverting each matrix, so the
        result can be indexed either way round. Default False.

    Returns
    -------
    dict
        Nested by frame of reference, outer key fixed and inner key moving::

            {fixed_frame_uid: {moving_frame_uid: [16 floats]}}

        Matrices within each MatrixSequence are composed in DICOM order (last @ ... @ first).
        Each returned matrix is a flat row-major 4x4 running **moving to fixed** -- DICOM's
        direction, and the inverse of what ``sitk.Resample`` wants. Empty if the object
        carries no FrameOfReferenceUID.

        Registrations whose matrix is missing or unreadable are skipped with a warning on
        the ``pydicomrt.reg.parser`` logger, not raised.

    See Also
    --------
    SpatialRegistrationBuilder : Write this object.
    get_deformable_registrations : The deformable equivalent.

    Notes
    -----
    An object built by this library includes an identity item for the fixed frame itself,
    so the fixed UID normally appears as its own inner key mapping to the identity matrix.

    To *apply* a matrix, invert it and load it into a transform::

        inverse = np.linalg.inv(np.array(matrix, dtype=float).reshape(4, 4))
        transform = sitk.AffineTransform(3)
        transform.SetMatrix(inverse[:3, :3].ravel().tolist())
        transform.SetTranslation(inverse[:3, 3].tolist())

    Examples
    --------
    >>> reg = get_spatial_registrations(dcmread("registration.dcm"))   # doctest: +SKIP
    >>> matrix = reg[fixed_frame_uid][moving_frame_uid]                # doctest: +SKIP
    >>> len(matrix)                                                    # doctest: +SKIP
    16
    """
    fixed_frame_of_reference_uid = getattr(reg_ds, "FrameOfReferenceUID", None)
    if fixed_frame_of_reference_uid is None:
        return {}

    reg_dict = {}
    reg_dict[fixed_frame_of_reference_uid] = {}

    for reg_index, reg in enumerate(getattr(reg_ds, "RegistrationSequence", [])):
        moving_frame_of_reference_uid = getattr(reg, "FrameOfReferenceUID", None)
        matrix_registration_sequence = getattr(reg, "MatrixRegistrationSequence", None)
        if moving_frame_of_reference_uid is None:
            continue
        if matrix_registration_sequence is None or len(matrix_registration_sequence) == 0:
            logger.warning("Matrix registration sequence is empty for %s", moving_frame_of_reference_uid)
            continue

        try:
            if len(matrix_registration_sequence) != 1:
                raise ValueError("expected one MatrixRegistrationSequence item")
            matrices = matrix_registration_sequence[0].MatrixSequence
            if not matrices:
                raise ValueError("empty MatrixSequence")
            matrix = np.eye(4)
            for item in matrices:
                matrix = np.asarray(item.FrameOfReferenceTransformationMatrix, dtype=float).reshape(4, 4) @ matrix
            matrix_matrix = matrix.ravel().tolist()
        except (AttributeError, IndexError, TypeError, ValueError):
            logger.warning("No transformation matrix found for %s", moving_frame_of_reference_uid)
            continue
        reg_dict[fixed_frame_of_reference_uid][moving_frame_of_reference_uid] = matrix_matrix

    if calc_inverse_matrix:
        reversed_reg_dict = {}
        for fixed_frame_of_reference_uid, moving_frame_of_reference_uid_dict in reg_dict.items():
            for moving_frame_of_reference_uid, matrix_matrix in moving_frame_of_reference_uid_dict.items():
                matrix_matrix = np.array(matrix_matrix).reshape((4, 4))
                inverse_matrix = np.linalg.inv(matrix_matrix)
                inverse_matrix = inverse_matrix.ravel().tolist()
                if moving_frame_of_reference_uid not in reversed_reg_dict:
                    reversed_reg_dict[moving_frame_of_reference_uid] = {}
                reversed_reg_dict[moving_frame_of_reference_uid][fixed_frame_of_reference_uid] = inverse_matrix
        for uid, entries in reversed_reg_dict.items():
            reg_dict.setdefault(uid, {}).update(entries)

        all_items = list(reg_dict.items())
        for src_uid, dst_dict in all_items:
            for dst_uid, mat_list in list(dst_dict.items()):
                if dst_uid in reg_dict and src_uid in reg_dict[dst_uid]:
                    continue
                mat = np.array(mat_list, dtype=float).reshape((4, 4))
                inv_mat = np.linalg.inv(mat).ravel().tolist()
                if dst_uid not in reg_dict:
                    reg_dict[dst_uid] = {}
                reg_dict[dst_uid][src_uid] = inv_mat
    return reg_dict


def get_deformable_registrations(reg_ds: Dataset) -> list:
    """
    Read the deformation fields out of a Deformable Spatial Registration object.

    Parameters
    ----------
    reg_ds : Dataset
        Deformable Spatial Registration object.

    Returns
    -------
    list of dict
        One entry per registered frame of reference::

            [{"SourceFrameOfReferenceUID": str,
              "PreDeformationMatrixRegistration":  [16 floats],
              "PostDeformationMatrixRegistration": [16 floats],
              "DeformableRegistrationGrid": {
                  "GridDimensions":          [nx, ny, nz],
                  "GridResolution":          [dx, dy, dz],
                  "ImagePositionPatient":    [x, y, z],
                  "ImageOrientationPatient": [6 floats],
                  "VectorGridData":          np.ndarray,   # (z, y, x, 3)
              }}]

        ``VectorGridData`` is reshaped into ``(z, y, x, 3)``; the last axis holds the
        ``(dx, dy, dz)`` displacement in mm. Note the axis order is numpy's, the reverse
        of the ``GridDimensions`` it was reshaped from.

        The transform composes as pre-matrix, then grid, then post-matrix, and runs
        **fixed to moving** (Registered RCS to Source RCS). Absent pre/post matrices
        are returned as identity matrices. Only grid-bearing entries are returned;
        matrix-only entries and entries with malformed or missing grid elements are skipped.

    See Also
    --------
    DeformableSpatialRegistrationBuilder : Write this object.
    get_spatial_registrations : The rigid/affine equivalent.

    Examples
    --------
    >>> regs = get_deformable_registrations(dcmread("deformable.dcm"))  # doctest: +SKIP
    >>> field = regs[0]["DeformableRegistrationGrid"]["VectorGridData"] # doctest: +SKIP
    >>> field.shape                                                     # doctest: +SKIP
    (120, 512, 512, 3)
    """
    reg_dict_list = []
    for reg in getattr(reg_ds, "DeformableRegistrationSequence", []):
        try:
            source_for_uid = getattr(reg, "SourceFrameOfReferenceUID", None)
            pre_deform_seq = getattr(reg, "PreDeformationMatrixRegistrationSequence", None)
            post_deform_seq = getattr(reg, "PostDeformationMatrixRegistrationSequence", None)
            grid_seq = getattr(reg, "DeformableRegistrationGridSequence", None)
            if source_for_uid is None or not grid_seq or len(grid_seq) != 1:
                continue

            matrices = []
            for seq in (pre_deform_seq, post_deform_seq):
                if seq is None:
                    matrices.append(np.eye(4).ravel().tolist())
                elif len(seq) == 1:
                    matrices.append(np.asarray(seq[0].FrameOfReferenceTransformationMatrix, dtype=float).reshape(16).tolist())
                else:
                    raise ValueError("expected one pre/post matrix item")
            pre_matrix, post_matrix = matrices
            deformed_grid_data = grid_seq[0]
            grid_dim = getattr(deformed_grid_data, "GridDimensions", None)
            vector_grid_data = getattr(deformed_grid_data, "VectorGridData", None)
            if pre_matrix is None or post_matrix is None or grid_dim is None or vector_grid_data is None:
                continue
            if any(not hasattr(deformed_grid_data, key) for key in
                   ("GridResolution", "ImagePositionPatient", "ImageOrientationPatient")):
                continue

            reg_dict = {
                'SourceFrameOfReferenceUID': source_for_uid,
                'PreDeformationMatrixRegistration': pre_matrix,
                'PostDeformationMatrixRegistration': post_matrix,
                'DeformableRegistrationGrid': {
                    'GridDimensions': grid_dim,
                    'GridResolution': getattr(deformed_grid_data, "GridResolution", None),
                    'ImagePositionPatient': getattr(deformed_grid_data, "ImagePositionPatient", None),
                    'ImageOrientationPatient': getattr(deformed_grid_data, "ImageOrientationPatient", None),
                },
            }

            vector_grid_unpack_data = np.frombuffer(vector_grid_data, dtype='<f4')
            deformed_array = vector_grid_unpack_data.reshape((grid_dim[2], grid_dim[1], grid_dim[0], 3))
            reg_dict['DeformableRegistrationGrid']['VectorGridData'] = deformed_array
            reg_dict_list.append(reg_dict)
        except (AttributeError, IndexError, TypeError, ValueError):
            continue

    return reg_dict_list


if __name__ == '__main__':
    import pydicom
    deformable_dcm_path = 'example/data/DF_001/REG/DR.dcm'
    reg_ds = pydicom.dcmread(deformable_dcm_path)
    reg_dict_list = get_deformable_registrations(reg_ds)
    # print(reg_dict_list)
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['VectorGridData'].shape)
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['GridDimensions'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['GridResolution'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['ImagePositionPatient'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['ImageOrientationPatient'])

    spatial_reg_dcm_path = 'example/data/DF_001/REG/RE.dcm'
    spatial_reg_ds = pydicom.dcmread(spatial_reg_dcm_path)
    spatial_reg_dict = get_spatial_registrations(spatial_reg_ds, calc_inverse_matrix=True)
    print(spatial_reg_dict)
