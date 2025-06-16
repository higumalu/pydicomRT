import numpy as np
import cv2

from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from utils.coordinate_transform import apply_transformation_to_3d_points
from .contour_process_method import contour_process


def add_mask3d_into_rsds(
    rs_ds,
    mask_volume,
    affine_map,
    roi_color,
    roi_number,
    roi_name,
    roi_description,
    roi_interpreted_type="ORGAN",
    contour_config: dict = {
        "ex_noise_size": 0,
        "in_noise_size": 0,
        "lowpass_ratio": 10,
        "ctr_precision": 8
        }
    ) -> Dataset:
    rs_ds.StructureSetROISequence.append(create_structure_set_roi(roi_number, roi_name, roi_description))
    rs_ds.ROIContourSequence.append(create_roi_contour_sequence(roi_number, roi_color, mask_volume, affine_map, contour_config))
    rs_ds.RTROIObservationsSequence.append(create_rtroi_observation(roi_number, roi_interpreted_type))
    return rs_ds

# StructureSetROISequence  --> [ROINumber, ReferencedFrameOfReferenceUID, ROIName, ROIGenerationAlgorithm=AUTOMATIC]
################################ StructureSetROISequence ####################################################################
def create_structure_set_roi(roi_number,
                             roi_name, roi_description,
                             roi_generation_algorithm="AUTOMATIC",
                             roi_generation_description="KumaGenerated") -> Dataset:
    # Structure Set ROI Sequence: Structure Set ROI 1
    structure_set_roi = Dataset()
    structure_set_roi.ROINumber = roi_number
    structure_set_roi.ROIName = roi_name
    structure_set_roi.ROIDescription = roi_description
    structure_set_roi.ROIGenerationAlgorithm = roi_generation_algorithm
    structure_set_roi.ROIGenerationDescription = roi_generation_description
    return structure_set_roi


# create ROIContourSequence -> ContourSequence -> [ContourImageSequence / ContourData / ContourNumber / ...]
################################ ROIContourSequence ####################################################################
def create_roi_contour_sequence(roi_number, roi_color, mask_volume, affine_map, ctr_config) -> Dataset:     # == create_roi_contour() in rt-utils
    roi_contour_sequence = Dataset()
    roi_contour_sequence.ReferencedROINumber = roi_number    # also in create_structure_set_roi
    roi_contour_sequence.ROIDisplayColor = roi_color
    roi_contour_sequence.ContourSequence = create_contour_sequence(mask_volume, affine_map, ctr_config)
    return roi_contour_sequence


def create_roi_contour_sequence_from_ctr_dict(roi_number, roi_color, ctr_dict, affine_map) -> Dataset:     # == create_roi_contour() in rt-utils
    roi_contour_sequence = Dataset()
    roi_contour_sequence.ReferencedROINumber = roi_number    # also in create_structure_set_roi
    roi_contour_sequence.ROIDisplayColor = roi_color
    roi_contour_sequence.ContourSequence = create_contour_sequence_from_ctr_dict(ctr_dict, affine_map)
    return roi_contour_sequence


# 3Dmask -[for loop]-> 2Dmask -[cv2.findContours]-> poly_list -[for loop]-> 2Dcontour
def create_contour_sequence(mask_volume, affine_map, ctr_config) -> Sequence:
    contour_sequence = Sequence()
    contour_data_seq = volume_to_contour_list(mask_volume, affine_map, ctr_config)
    for index, contour_data in enumerate(contour_data_seq):
        contour = create_contour_sequence_block(contour_data, precision=ctr_config['ctr_precision'])
        contour_sequence.append(contour)
    return contour_sequence


def create_contour_sequence_from_ctr_dict(ctr_dict, affine_map) -> Sequence:
    contour_sequence = Sequence()
    for z_key, polys in ctr_dict.items():
        for poly in polys:
            contour = Dataset()
            z_col = np.full((poly.shape[0], 1), z_key)
            contour_data = np.hstack((poly, z_col))
            contour_data = apply_transformation_to_3d_points(contour_data, affine_map)
            contour_data = list(contour_data.reshape(-1))
            contour.ContourData = contour_data
            contour.NumberOfContourPoints = (len(contour_data) / 3)
            contour.ContourGeometricType = "CLOSED_PLANAR"
            contour_sequence.append(contour)
    return contour_sequence


# ContourSequence -> [ContourImageSequence, ContourGeometricType, NumberOfContourPoints, ContourNumber, ContourData]
# -> ContourImageSequence -> [ReferencedSOPClassUID, ReferencedSOPInstanceUID]
def create_contour_sequence_block(contour_data, precision=4) -> Dataset:
    # contour_image = Dataset()
    # contour_image.ReferencedSOPClassUID = image_slice.SOPClassUID
    # contour_image.ReferencedSOPInstanceUID = image_slice.SOPInstanceUID
    # contour_image_sequence = Sequence()
    # contour_image_sequence.append(contour_image)

    contour = Dataset()
    # contour.ContourImageSequence = contour_image_sequence
    contour.ContourGeometricType = "CLOSED_PLANAR"
    contour.NumberOfContourPoints = (len(contour_data) / 3)
    contour_data = [round(coord, precision) for coord in contour_data]
    contour.ContourData = contour_data

    return contour

# 3Dmask
# -[for loop with ds_list]-> 2Dmask with series slice
# -[cv2.findContours]-> cv2_contours
# -[for loop]-> contour(polygon)
# -[smooth_ctr (maybe)]-> contour
# -[coordination transform]-> dicom fromat contour (ContourData)
def volume_to_contour_list(mask_volume, affine_map, ctr_config):
    formatted_contours = []
    for index, mask_slice in enumerate(mask_volume):
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
            # interface for contour smoothing and other process
            x_points, y_points = contour_process(
                x_points, y_points, poly_hierarchy,
                external_noise_size=ctr_config["ex_noise_size"],
                internal_noise_size=ctr_config["in_noise_size"],
                low_pass_ratio=ctr_config["lowpass_ratio"]
                )
            # check if contour is not empty
            if len(x_points) < 3 or len(y_points) < 3: continue

            contour = np.stack([x_points, y_points]).T
            contour = np.concatenate((np.array(contour), np.full((len(contour), 1), index)), axis=1)
            # transformation_matrix = get_pixel_to_patient_transformation_matrix(image_series)
            transformed_contour = apply_transformation_to_3d_points(contour, affine_map)

            dicom_formatted_contour = np.ravel(transformed_contour).tolist()
            formatted_contours.append(dicom_formatted_contour)

    return formatted_contours

# ------------------------------- RTROIObservationsSequence --------------------------------------- #
def create_rtroi_observation(roi_number, roi_interpreted_type='ORGAN') -> Dataset:
    rtroi_observation = Dataset()
    rtroi_observation.ObservationNumber = roi_number
    rtroi_observation.ReferencedROINumber = roi_number

    rtroi_observation.ROIObservationDescription = "Type:Soft,Range:*/*,Fill:0,Opacity:0.0,Thickness:1,LineThickness:2,read-only:false"
    rtroi_observation.private_creators = "higumalu"
    rtroi_observation.RTROIInterpretedType = roi_interpreted_type
    rtroi_observation.ROIInterpreter = ""
    return rtroi_observation
