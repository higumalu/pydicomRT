import cv2
import numpy as np

from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from pydicomrt.utils.coordinate_transform import apply_transformation_to_3d_points, get_pixel_to_patient_transformation_matrix
from .contour_process_method import contour_process

# ---------------------------------- interface ----------------------------------------- #

def add_contour_sequence_from_dcm_ctr_dict(
    rs_ds: Dataset,
    image_ds_list: list[Dataset],
    roi_number: int,
    dcm_ctr_dict: dict
    ) -> Sequence:
    if not hasattr(rs_ds, 'ROIContourSequence'):
        raise ValueError("ROIContourSequence does not exist")

    for index, roi_contour_sequence in enumerate(rs_ds.ROIContourSequence):
        roi_contour_sequence.get("ReferencedROINumber", None)
        if roi_contour_sequence.ReferencedROINumber == roi_number:
            roi_contour_sequence.ContourSequence = create_contour_sequence_from_dcm_ctr_dict(
                image_ds_list,
                dcm_ctr_dict,
            )
            rs_ds.ROIContourSequence[index] = roi_contour_sequence
            return rs_ds
        else:
            continue
    raise ValueError("roi_number not found")

def add_contour_sequence_from_mask3d(
    rs_ds: Dataset,
    image_ds_list: list[Dataset],
    roi_number: int,
    mask_volume: np.ndarray,
    ctr_config: dict = {
        "ex_noise_size": 10,
        "in_noise_size": 10,
        "lowpass_ratio": 10,
        "ctr_precision": 8
        }
    ) -> Sequence:
    if not hasattr(rs_ds, 'ROIContourSequence'):
        raise ValueError("ROIContourSequence does not exist")

    for index, roi_contour_sequence in enumerate(rs_ds.ROIContourSequence):
        roi_contour_sequence.get("ReferencedROINumber", None)
        if roi_contour_sequence.ReferencedROINumber == roi_number:
            roi_contour_sequence.ContourSequence = create_contour_sequence_from_mask3d(
                image_ds_list,
                mask_volume,
                ctr_config,
            )
            rs_ds.ROIContourSequence[index] = roi_contour_sequence
            return rs_ds
        else:
            continue
    raise ValueError("roi_number not found")


# ------------------------------ from_dcm_ctr_dict ------------------------------------------------- #

def create_contour_sequence_from_dcm_ctr_dict(
    image_ds_list: list[Dataset],
    dcm_ctr_dict: dict
    ) -> Sequence:
    """
    dcm_ctr_dict: {
        "sop_instance_uid": {
            "contours": list,
            "sop_class_uid": str,
        }
    }
    """
    contour_sequence = Sequence()
    for sop_instance_uid, ctr_info in dcm_ctr_dict.items():
        sop_class_uid = ctr_info["sop_class_uid"]
        ctr_list = ctr_info["contours"]
        for contour_data in ctr_list:
            contour_sequence.append(create_contour_sequence_block(contour_data, sop_class_uid, sop_instance_uid))
    return contour_sequence

def create_contour_sequence_block(contour_data, sop_class_uid, sop_instance_uid) -> Dataset:
    contour_image = Dataset()
    contour_image.ReferencedSOPClassUID = sop_class_uid
    contour_image.ReferencedSOPInstanceUID = sop_instance_uid

    contour_image_sequence = Sequence()
    contour_image_sequence.append(contour_image)

    contour = Dataset()
    contour.ContourGeometricType = "CLOSED_PLANAR"
    contour.ContourImageSequence = contour_image_sequence
    contour.NumberOfContourPoints = (len(contour_data) / 3)
    contour.ContourData = list(contour_data)

    return contour


# ------------------------------------ from_mask3d ------------------------------------------------- #

def create_contour_sequence_from_mask3d(
    image_ds_list: list[Dataset],
    mask_volume: np.ndarray,
    ctr_config: dict
    ) -> Sequence:
    contour_sequence = Sequence()
    contour_data_seq = volume_to_contour_list(mask_volume, image_ds_list, ctr_config)
    for index, [contour_data, image_slice] in enumerate(contour_data_seq):
        sop_class_uid = image_slice.SOPClassUID
        sop_instance_uid = image_slice.SOPInstanceUID
        contour = create_contour_sequence_block(contour_data, sop_class_uid, sop_instance_uid)
        contour_sequence.append(contour)

    return contour_sequence

def volume_to_contour_list(mask_volume, image_series, ctr_config):
    formatted_contours = []
    for index, series_slice in enumerate(image_series):
        mask_slice = mask_volume[index, :, :]           # mask_volume.shape = (z, y, x)     ex: (144, 512, 512)
        if np.sum(mask_slice) == 0:
            continue
        cv2_contours, hierarchy = cv2.findContours(
            mask_slice.astype(np.uint8),
            cv2.RETR_CCOMP,
            cv2.CHAIN_APPROX_NONE
            )
        if cv2_contours is None or hierarchy is None: continue
        for hier_index, contour in enumerate(cv2_contours):
            poly_hierarchy = hierarchy[0][hier_index][3]
            x_points, y_points = contour.T[:, 0]
            x_points, y_points = contour_process(x_points, y_points, poly_hierarchy,
                                                 external_noise_size=ctr_config["ex_noise_size"],
                                                 internal_noise_size=ctr_config["in_noise_size"],
                                                 low_pass_ratio=ctr_config["lowpass_ratio"])
            if len(x_points) < 3 or len(y_points) < 3: continue
            contour = np.stack([x_points, y_points]).T
            contour = np.concatenate((np.array(contour), np.full((len(contour), 1), index)), axis=1)    # add index into contour-z-axis
            transformation_matrix = get_pixel_to_patient_transformation_matrix(image_series)
            transformed_contour = apply_transformation_to_3d_points(contour, transformation_matrix)
            dicom_formatted_contour = np.ravel(transformed_contour).tolist()
            dicom_formatted_contour = [round(crood, ctr_config["ctr_precision"]) for crood in dicom_formatted_contour]
            formatted_contours.append([dicom_formatted_contour, series_slice])
    return formatted_contours
