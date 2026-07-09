
import numpy as np
import SimpleITK as sitk

from pydicom.dataset import Dataset

def get_spatial_reg_dict(reg_ds: Dataset, calc_inverse_matrix: bool = False) -> dict:
    """
    Get the spatial registration dictionary from the registration dataset.
    Args:
        reg_ds (Dataset): The registration dataset.
        calc_inverse_matrix (bool): Whether to calculate the inverse matrix.
    Returns:
        dict: The spatial registration dictionary.

        dict: {
            a_foruid: {
                b_foruid: matrix_matrix,
                c_foruid: matrix_matrix,
                ...
            },
            b_foruid: {
                a_foruid: matrix_matrix,
                ...
            },
            c_foruid: {
                a_foruid: matrix_matrix,
                ...
            },

        }
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
            print(f"Matrix registration sequence is empty for {moving_frame_of_reference_uid}.")
            continue

        try:
            matrix_matrix = matrix_registration_sequence[-1].MatrixSequence[-1].FrameOfReferenceTransformationMatrix
        except (AttributeError, IndexError, TypeError):
            print(f"Matrix is not found for {moving_frame_of_reference_uid}.")
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
        reg_dict.update(reversed_reg_dict)

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


def get_deformable_reg_list(reg_ds: Dataset) -> list:
    reg_dict_list = []
    for reg in getattr(reg_ds, "DeformableRegistrationSequence", []):
        try:
            source_for_uid = getattr(reg, "SourceFrameOfReferenceUID", None)
            pre_deform_seq = getattr(reg, "PreDeformationMatrixRegistrationSequence", None)
            post_deform_seq = getattr(reg, "PostDeformationMatrixRegistrationSequence", None)
            grid_seq = getattr(reg, "DeformableRegistrationGridSequence", None)
            if source_for_uid is None or not pre_deform_seq or not post_deform_seq or not grid_seq:
                continue

            pre_matrix = getattr(pre_deform_seq[-1], "FrameOfReferenceTransformationMatrix", None)
            post_matrix = getattr(post_deform_seq[-1], "FrameOfReferenceTransformationMatrix", None)
            deformed_grid_data = grid_seq[-1]
            grid_dim = getattr(deformed_grid_data, "GridDimensions", None)
            vector_grid_data = getattr(deformed_grid_data, "VectorGridData", None)
            if pre_matrix is None or post_matrix is None or grid_dim is None or vector_grid_data is None:
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
    reg_dict_list = get_deformable_reg_list(reg_ds)
    # print(reg_dict_list)
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['VectorGridData'].shape)
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['GridDimensions'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['GridResolution'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['ImagePositionPatient'])
    # print(reg_dict_list[0]['DeformableRegistrationGrid']['ImageOrientationPatient'])

    spatial_reg_dcm_path = 'example/data/DF_001/REG/RE.dcm'
    spatial_reg_ds = pydicom.dcmread(spatial_reg_dcm_path)
    spatial_reg_dict = get_spatial_reg_dict(spatial_reg_ds, calc_inverse_matrix=True)
    print(spatial_reg_dict)