
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
    reg_dict = {}
    fixed_frame_of_reference_uid = reg_ds.FrameOfReferenceUID
    reg_dict[fixed_frame_of_reference_uid] = {}

    for reg_index, reg in enumerate(reg_ds.RegistrationSequence):
        moving_frame_of_reference_uid = getattr(reg, "FrameOfReferenceUID", None)
        matrix_registration_sequence = getattr(reg, "MatrixRegistrationSequence", None)
        if moving_frame_of_reference_uid is None:
            raise ValueError("Moving frame of reference UID is not found.")
        if matrix_registration_sequence is None or len(matrix_registration_sequence) == 0:
            print(f"Matrix registration sequence is empty for {moving_frame_of_reference_uid}.")
            continue

        try:
            matrix_matrix = matrix_registration_sequence[-1].MatrixSequence[-1].FrameOfReferenceTransformationMatrix
        except:
            print(f"Matrix is not found for {moving_frame_of_reference_uid}.")
            continue
        reg_dict[fixed_frame_of_reference_uid][moving_frame_of_reference_uid] = matrix_matrix

        # for matrix_reg_index, matrix_registration in enumerate(matrix_registration_sequence):
        #     matrix_sequence = matrix_registration.MatrixSequence
        #     if matrix_sequence is None or len(matrix_sequence) == 0:
        #         print(f"Matrix sequence is empty for {moving_frame_of_reference_uid}.")
        #         continue
        #     for matrix_index, matrix in enumerate(matrix_sequence):
        #         matrix_type = getattr(matrix, "FrameOfReferenceTransformationMatrixType", None)
        #         matrix_matrix = getattr(matrix, "FrameOfReferenceTransformationMatrix", None)
        #         if matrix_type is None or matrix_matrix is None:
        #             raise ValueError("Matrix type or matrix is not found.")
        #         if matrix_matrix is None:
        #             raise ValueError("Matrix is not found.")
        #         reg_dict[fixed_frame_of_reference_uid][moving_frame_of_reference_uid] = matrix_matrix

    if calc_inverse_matrix:
        reversed_reg_dict = {}
        for fixed_frame_of_reference_uid, moving_frame_of_reference_uid_dict in reg_dict.items():
            for moving_frame_of_reference_uid, matrix_matrix in moving_frame_of_reference_uid_dict.items():
                matrix_matrix = np.array(matrix_matrix).reshape((4, 4))
                inverse_matrix = np.linalg.inv(matrix_matrix)
                inverse_matrix = inverse_matrix.ravel().tolist()
                reversed_reg_dict[moving_frame_of_reference_uid] = {}
                reversed_reg_dict[moving_frame_of_reference_uid][fixed_frame_of_reference_uid] = inverse_matrix
        reg_dict.update(reversed_reg_dict)
    return reg_dict


def get_deformable_reg_list(reg_ds: Dataset) -> list:
    reg_dict_list = []
    for reg in reg_ds.DeformableRegistrationSequence:
        reg_dict = {}
        reg_dict['SourceFrameOfReferenceUID'] = reg.SourceFrameOfReferenceUID
        reg_dict['PreDeformationMatrixRegistration'] = reg.PreDeformationMatrixRegistrationSequence[-1].FrameOfReferenceTransformationMatrix
        reg_dict['PostDeformationMatrixRegistration'] = reg.PostDeformationMatrixRegistrationSequence[-1].FrameOfReferenceTransformationMatrix
        reg_dict['DeformableRegistrationGrid'] = {}
        deformed_grid_data = reg.DeformableRegistrationGridSequence[-1]

        reg_dict['DeformableRegistrationGrid']['GridDimensions'] = deformed_grid_data.GridDimensions
        reg_dict['DeformableRegistrationGrid']['GridResolution'] = deformed_grid_data.GridResolution
        reg_dict['DeformableRegistrationGrid']['ImagePositionPatient'] = deformed_grid_data.ImagePositionPatient
        reg_dict['DeformableRegistrationGrid']['ImageOrientationPatient'] = deformed_grid_data.ImageOrientationPatient

        grid_dim = deformed_grid_data.GridDimensions
        vector_grid_data = deformed_grid_data.VectorGridData
        # vector_grid_data = unpack(f"<{len(vector_grid_data) // 4}f", vector_grid_data)
        # deformed_array = np.reshape(vector_grid_data, grid_dim[::-1] + [3,])

        vector_grid_unpack_data = np.frombuffer(vector_grid_data, dtype='<f4')
        deformed_array = vector_grid_unpack_data.reshape((grid_dim[2], grid_dim[1], grid_dim[0], 3))

        reg_dict['DeformableRegistrationGrid']['VectorGridData'] = deformed_array
        reg_dict_list.append(reg_dict)

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