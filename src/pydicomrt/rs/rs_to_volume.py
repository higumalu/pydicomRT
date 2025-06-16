'''
Author:   BearLu
Status:   Final
Created:  2024-10-11
Python-Version: 3.10.0
Coding: UTF-8
'''

import numpy as np
import cv2


def self_round(x: float, decimal: int) -> float:
    x = x * 10**decimal
    if x >= 0:
        return float(int(x + 0.5)) / 10**decimal
    else:
        return float(int(x - 0.5)) / 10**decimal


def get_normal_vector(contour_data):
    ctrs = np.reshape(contour_data, (-1, 3))
    normal_list = [np.cross(ctrs[i + 2] - ctrs[0], ctrs[i + 1] - ctrs[0]) for i in range(0, len(ctrs) - 2, 2)]
    normal_list = [normal / np.linalg.norm(normal) for normal in normal_list]
    normal_vector = np.mean(normal_list, axis=0)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)
    return normal_vector


def get_normal_vector_bear(contour_data):
    ctrs = np.reshape(contour_data, (-1, 3))
    min_index = np.unravel_index(np.argmin(ctrs[:, 0] + ctrs[:, 1]), ctrs[:, 0].shape)
    x_max_index = np.unravel_index(np.argmax(ctrs[:, 0]), ctrs[:, 0].shape)
    y_max_index = np.unravel_index(np.argmax(ctrs[:, 1]), ctrs[:, 1].shape)

    v1 = ctrs[min_index] - ctrs[x_max_index]
    v2 = ctrs[min_index] - ctrs[y_max_index]
    normal_vector = np.cross(v1, v2)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)
    return normal_vector


def get_slice_index(contour_data, normal_vector):
    ctrs = np.reshape(contour_data, (-1, 3))
    index = np.dot(ctrs[0], normal_vector)
    return index


def project_point_onto_plane(point, plane_normal, d):
    '''
    a, b, c = plane_normal
    x0, y0, z0 = origin
    d = a * x0 + b * y0 + c * z0
    '''
    x_p, y_p, z_p = point
    a, b, c = plane_normal
    D = (a * x_p + b * y_p + c * z_p + d) / np.sqrt(a**2 + b**2 + c**2)
    proj_vector = np.array([a * D, b * D, c * D]) / np.sqrt(a**2 + b**2 + c**2)
    projected_point = np.array([x_p, y_p, z_p]) - proj_vector
    return projected_point


def calc_normal_projection(points: np.ndarray, normal_vector: np.ndarray) -> np.ndarray:
    """
    Calculate normal projection for a set of 3D points.
    Args:
        points: 3D points in homogeneous coordinates. (n*3 array)
        normal_vector: Normal vector of the plane.  (1*3 array)
    Returns:
        Normal projection of the points in the plane.    (n*1 array)
    """
    projection_list = [np.dot(point, normal_vector) for point in points]
    return projection_list


def calc_obb(points: np.ndarray, normal_vector: np.ndarray) -> np.ndarray:
    """
    Calculate oriented bounding box for a set of 3D points.
    Args:
        points: 3D points in homogeneous coordinates. (n*3 array)
        normal_vector: Normal vector of the plane.  (1*3 array)
    Returns:
        Oriented bounding box in homogeneous coordinates.    (8*3 array)
    """
    v1 = np.cross([0, 1, 0], normal_vector)
    if np.linalg.norm(v1) == 0:
        v1 = np.cross([1, 0, 0], normal_vector)
    v1 = v1 / np.linalg.norm(v1)
    v2 = np.cross(normal_vector, v1)
    v2 = v2 / np.linalg.norm(v2)
    R = np.array([v1, v2, normal_vector]).T
    points_rotated = points @ R
    min_bound = np.min(points_rotated, axis=0)
    max_bound = np.max(points_rotated, axis=0)
    obb_corners_local = np.array([[min_bound[0], min_bound[1], min_bound[2]],   # origin
                                  [min_bound[0], min_bound[1], max_bound[2]],   # z
                                  [min_bound[0], max_bound[1], min_bound[2]],   # y
                                  [max_bound[0], min_bound[1], min_bound[2]],   # x
                                  [min_bound[0], max_bound[1], max_bound[2]],   # yz
                                  [max_bound[0], min_bound[1], max_bound[2]],   # zx
                                  [max_bound[0], max_bound[1], min_bound[2]],   # xy
                                  [max_bound[0], max_bound[1], max_bound[2]]])  # xyz
    obb_corners_original = obb_corners_local @ R.T
    obb_corners_original
    return obb_corners_original


def apply_transformation_to_3d_points(points: np.ndarray, transformation_matrix: np.ndarray):
    vec = np.concatenate((points, np.ones((points.shape[0], 1))), axis=1)
    return vec.dot(transformation_matrix.T)[:, :3]


def round_3d_points(points):
    return (points + 0.5).astype(int)


def get_roi_number_name_map(rs_ds):
    roi_number_name_map = {}
    if hasattr(rs_ds, 'StructureSetROISequence'):
        for ssroi in rs_ds.StructureSetROISequence:
            roi_name = ssroi.ROIName
            roi_number = ssroi.ROINumber
            roi_number_name_map[roi_number] = roi_name
    return roi_number_name_map


def get_series_z_spacing(series_ds_list):
    slice_position_list = [ds.ImagePositionPatient for ds in series_ds_list]
    slice_distance_list = [np.linalg.norm(np.array(slice_position_list[index + 1]) - np.array(slice_position_list[index])) for index in range(len(slice_position_list) - 1)]
    slice_distance_list = [dist for dist in slice_distance_list if dist > 0.1]
    if not slice_distance_list:
        raise Exception("No valid slice distance found")
    return min(slice_distance_list)


def calc_image_series_affine_mapping(series_ds_list):
    first_slice = series_ds_list[0]
    origin = np.array(first_slice.ImagePositionPatient)
    x_spacing, y_spacing = first_slice.PixelSpacing
    z_spacing = get_series_z_spacing(series_ds_list)
    v_x = np.array(first_slice.ImageOrientationPatient[:3])
    v_y = np.array(first_slice.ImageOrientationPatient[3:])
    v_z = np.cross(v_x, v_y)
    # check image orientation is legal
    if not np.allclose(np.dot(v_x, v_y), 0.0, atol=1e-3) or not np.allclose(np.linalg.norm(v_z), 1.0, atol=1e-3):
        raise Exception("Invalid Image Orientation (Patient) attribute")
    linear = np.identity(3, dtype=np.float32)
    linear[0, :3] = v_x / x_spacing
    linear[1, :3] = v_y / y_spacing
    linear[2, :3] = v_z / z_spacing
    affine_mapping = np.identity(4, dtype=np.float32)
    affine_mapping[:3, :3] = linear
    affine_mapping[:3, 3] = origin.dot(-linear.T)
    mask_volume_shape = (len(series_ds_list), first_slice.Rows, first_slice.Columns)
    return affine_mapping, mask_volume_shape


def calc_rs_affine_mapping(rs_ds):
    if hasattr(rs_ds, 'ROIContourSequence'):
        all_slice_ctr = {}
        all_ctr = []
        all_normal_vector = []
        for roi_contour in rs_ds.ROIContourSequence:
            if hasattr(roi_contour, 'ContourSequence'):
                for rtroi in roi_contour.ContourSequence:
                    poly = rtroi.ContourData
                    all_ctr.extend(poly)
                    if len(poly) < 9: continue
                    # v = get_normal_vector_bear(poly)
                    normal = get_normal_vector(poly)
                    all_normal_vector.append(normal)
                    sop_instance_uid = rtroi.ContourImageSequence[0].ReferencedSOPInstanceUID   # careful use
                    if sop_instance_uid not in all_slice_ctr:
                        all_slice_ctr[sop_instance_uid] = [poly]
                    else:
                        all_slice_ctr[sop_instance_uid].append(poly)

        ##### calc normal vector #####
        all_normal_vector = np.array(all_normal_vector)
        mean_normal = np.mean(all_normal_vector, axis=0)
        normal = mean_normal / np.linalg.norm(mean_normal)
        # print(normal)
        all_ctr = np.reshape(all_ctr, (-1, 3))
        # print(all_ctr.shape)
        # print(all_ctr[:, 0].min(), all_ctr[:, 0].max())
        # print(all_ctr[:, 1].min(), all_ctr[:, 1].max())
        # print(all_ctr[:, 2].min(), all_ctr[:, 2].max())
        # print(all_ctr)
        ##### calc orented bounding box #####
        obb = calc_obb(all_ctr, normal)
        print(obb, obb.shape)
        origin = obb[0]
        # xyz = obb[7]
        v_z = obb[1] - obb[0]
        v_y = obb[2] - obb[0]
        v_x = obb[3] - obb[0]
        v_z = v_z / np.linalg.norm(v_z)
        v_y = v_y / np.linalg.norm(v_y)
        v_x = v_x / np.linalg.norm(v_x)
        # print(v_z, v_y, v_x)
        ##### calc z_spacing #####
        projection_list = calc_normal_projection(all_ctr, normal)
        # print(projection_list)
        projection_list, index_list = np.unique(projection_list, return_index=True)
        sorted_index = np.argsort(projection_list)
        projection_list = projection_list[sorted_index]
        # print(projection_list)
        slice_distance_list = [self_round(abs(projection_list[i + 1] - projection_list[i]), 2) for i in range(len(projection_list) - 1)]
        # print(slice_distance_list)
        slice_distance_list = [x for x in slice_distance_list if x != 0]
        if len(slice_distance_list) == 0:
            raise ValueError('no slice distance list')
        # print(slice_distance_list)
        slice_distance_list = list(dict.fromkeys(slice_distance_list))
        # print(slice_distance_list)
        # print(np.min(slice_distance_list), np.max(slice_distance_list), np.median(slice_distance_list))
        ##### calc affine matrix #####
        x_spacing = 0.5
        y_spacing = 0.5
        z_spacing = min(slice_distance_list)
        linear = np.identity(3, dtype=np.float32)
        linear[0, :3] = v_x / x_spacing
        linear[1, :3] = v_y / y_spacing
        linear[2, :3] = v_z / z_spacing
        affine_mapping = np.identity(4, dtype=np.float32)
        affine_mapping[:3, :3] = linear
        affine_mapping[:3, 3] = origin.dot(-linear.T)
        # print(affine_mapping)
        new_points = apply_transformation_to_3d_points(all_ctr, affine_mapping)
        new_points = round_3d_points(new_points)
        # print(new_points[:, 0].min(), new_points[:, 0].max())
        # print(new_points[:, 1].min(), new_points[:, 1].max())
        # print(new_points[:, 2].min(), new_points[:, 2].max())
        # for point in new_points:
        #     print(point)
        mask_volume_shape = (new_points[:, 2].max() + 1 - new_points[:, 2].min(),
                             new_points[:, 1].max() + 1 - new_points[:, 1].min(),
                             new_points[:, 0].max() + 1 - new_points[:, 0].min())
    return affine_mapping, mask_volume_shape


def inverse_affine_mapping(matrix):
    # Extract the linear part (top-left 3x3 submatrix)
    linear_part = matrix[:3, :3]
    # Extract the translation vector (top-right 3x1 vector)
    translation = matrix[:3, 3]

    # Compute the inverse of the linear part
    linear_part_inv = np.linalg.inv(linear_part)

    # Compute the new translation vector
    translation_inv = -np.dot(linear_part_inv, translation)

    # Create the inverse affine matrix
    inverse_matrix = np.eye(4)  # Start with identity matrix
    inverse_matrix[:3, :3] = linear_part_inv  # Set the inverse linear part
    inverse_matrix[:3, 3] = translation_inv  # Set the inverse translation

    return inverse_matrix


def rtstruct_to_mask_dict(rs_ds, affine_mapping, mask_volume_shape, roi_list=["all"], packbits=False):
    """
    Convert a DICOM RT Structure Set to a dictionary of masks.

    Parameters:
    - rs_ds: DICOM RT Structure Set dataset.
    - affine_mapping: Affine mapping matrix.
    - mask_volume_shape: Shape of the mask volume.
    - roi_list: List of ROI names to include in the dictionary. If "all", include all ROIs.
    - packbits: Whether to pack the masks using packbits.

    Returns:
    - roi_dict: Dictionary of masks, where keys are ROI names and values are masks.
    """
    roi_number_name_map = get_roi_number_name_map(rs_ds)
    if hasattr(rs_ds, 'ROIContourSequence'):
        roi_dict = {}
        for roi_contour in rs_ds.ROIContourSequence:
            roi_number = roi_contour.ReferencedROINumber
            roi_name = roi_number_name_map[roi_number]
            if roi_list != ["all"] and roi_name not in roi_list:
                continue
            roi_dict[roi_name] = {}
            if hasattr(roi_contour, 'ContourSequence'):
                mask_volume = np.zeros(mask_volume_shape, dtype=np.int32)
                for rtroi in roi_contour.ContourSequence:
                    # sop_instance_uid = rtroi.ContourImageSequence[0].ReferencedSOPInstanceUID
                    poly = rtroi.ContourData
                    if len(poly) < 9: continue
                    ctrs = np.reshape(poly, (-1, 3))
                    transform_ctrs = apply_transformation_to_3d_points(ctrs, affine_mapping)
                    mask_index = round_3d_points(transform_ctrs)
                    mask_2d_index = mask_index[:, :2].astype(np.int32)
                    mask_slice_index = mask_index[0, 2]     # TODO check index are same
                    # print(mask_slice_index)
                    slice_mask = mask_volume[mask_slice_index, :, :].astype(np.int32)
                    tmp_mask = np.zeros_like(slice_mask, dtype=np.int32)
                    # print(slice_mask.shape)
                    # print(mask_2d_index)
                    tmp_mask = cv2.fillPoly(img=tmp_mask, pts=[mask_2d_index], color=1)
                    # print(tmp_mask.shape, tmp_mask.dtype, tmp_mask.sum(), tmp_mask.max(), tmp_mask.min())
                    slice_mask = slice_mask + tmp_mask
                    mask_volume[mask_slice_index, :, :] = slice_mask
                    # print(slice_mask.shape, slice_mask.dtype, slice_mask.sum())
                mask_volume[mask_volume % 2 == 0] = 0
                mask_volume[mask_volume % 2 == 1] = 1
                # print(mask_volume.shape, mask_volume.dtype, mask_volume.sum())
                if packbits:
                    mask_volume = np.packbits(mask_volume, axis=-1)
                roi_dict[roi_name]['mask_volume'] = mask_volume
                roi_dict[roi_name]['affine_mapping'] = affine_mapping

    return roi_dict


# if __name__ == '__main__':
    # start = time.time()

    # from rs_to_volume import load_sorted_image_series

    # mr_path = 'data/mr_data/MR'
    # # rs_path = 'data/ct_data/RS/107.dcm'
    # rs_path = 'data/mr_data/RS/920.dcm'

    # ds_list = load_sorted_image_series(mr_path)
    # rs_ds = pydicom.dcmread(rs_path)
    # # print(rs_ds[0x30060050])

    # # af_map, mask_shape = calc_image_series_affine_mapping(ds_list)
    # af_map, mask_shape = calc_rs_affine_mapping(rs_ds)
    # roi_dict = rtstruct_to_mask_dict(rs_ds, af_map, mask_shape)
    # inv_af_map = inverse_affine_mapping(af_map)
    # print(af_map, inv_af_map)

    # import random
    # from rs_builder_v2 import create_rtstruct_dataset, edit_required_elements
    # from rs_builder_v2 import add_mask3d_into_rsds
    # rs_ds = create_rtstruct_dataset(ds_list)
    # rs_ds = edit_required_elements(rs_ds, structure_set_label="RSlabel")
    # for index, roi_name in enumerate(roi_dict.keys()):
    #     roi_number = str(index + 1)
    #     roi_color = [random.randint(0, 255) for _ in range(3)]
    #     mask_volume = roi_dict[roi_name]['mask_volume']
    #     roi_description = roi_name
    #     rs_ds = add_mask3d_into_rsds(rs_ds, mask_volume, inv_af_map, roi_color, roi_number, roi_name, roi_description)

    # rs_ds.save_as('data/mr_data/RS/rs.dcm')

    # from matplotlib import pyplot as plt

    # for roi_name in roi_dict.keys():
    #     print(roi_name)
    #     if roi_name != 'MR_OpticNrv_R': continue
    #     mask_volume = roi_dict[roi_name]['mask_volume']
    #     affine_mapping = roi_dict[roi_name]['affine_mapping']
    #     for index in range(0, mask_volume.shape[0], 1):
    #         # plt.imshow(mask_volume[index, :, :], cmap='gray')
    #         # plt.show()
    #         pass

    # print(time.time() - start)
