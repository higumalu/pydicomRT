
import numpy as np

from pydicom.dataset import Dataset

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
    print(reg_dict_list[0]['DeformableRegistrationGrid']['VectorGridData'].shape)
    print(reg_dict_list[0]['DeformableRegistrationGrid']['GridDimensions'])
    print(reg_dict_list[0]['DeformableRegistrationGrid']['GridResolution'])
    print(reg_dict_list[0]['DeformableRegistrationGrid']['ImagePositionPatient'])
    print(reg_dict_list[0]['DeformableRegistrationGrid']['ImageOrientationPatient'])
