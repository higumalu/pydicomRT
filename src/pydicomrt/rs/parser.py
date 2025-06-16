from pydicom.dataset import Dataset


def get_contour_dict(
    rs_ds: Dataset,
    ) -> dict:
    """
    Get contour dict from rs_ds.
    Args:
        rs_ds (Dataset): rs_ds
    Returns:
        dict: contour_dict
    Limitations:
        1. Only one image is referenced.
        2. Only one contour is referenced.
    Structure of contour_dict:
    {
        roi_number: {
            'color': [R, G, B],  # ROI display color
            'name': str,         # ROI name
            'dcm_contour': {     # Contours organized by image
                sop_instance_uid: {
                    'sop_class_uid': str,  # SOP Class UID of referenced image
                    'contours': [          # List of contour data arrays
                        [x1, y1, z1, x2, y2, z2, ...],  # First contour
                        [x1, y1, z1, x2, y2, z2, ...],  # Second contour
                        ...
                    ]
                }
            }
        }
    }
    """
    contour_dict = {}
    number_name_map = get_roi_number_to_name(rs_ds)
    for roi_contour_sequence in rs_ds.ROIContourSequence:
        roi_number = roi_contour_sequence.ReferencedROINumber
        roi_color = roi_contour_sequence.ROIDisplayColor
        roi_name = number_name_map[roi_number]
        contour_dict[roi_number] = {}
        contour_dict[roi_number]['color'] = roi_color
        contour_dict[roi_number]['name'] = roi_name
        contour_dict[roi_number]['dcm_contour'] = {}
        for contour_sequence in roi_contour_sequence.ContourSequence:
            sop_instance_uid = contour_sequence.ContourImageSequence[0].ReferencedSOPInstanceUID    # TODO ensure only one image is referenced
            sop_class_uid = contour_sequence.ContourImageSequence[0].ReferencedSOPClassUID
            contour_data = contour_sequence.ContourData

            if sop_instance_uid not in contour_dict[roi_number]['dcm_contour']:
                contour_dict[roi_number]['dcm_contour'][sop_instance_uid] = {}
                contour_dict[roi_number]['dcm_contour'][sop_instance_uid]['sop_class_uid'] = sop_class_uid
                contour_dict[roi_number]['dcm_contour'][sop_instance_uid]['contours'] = []

            contour_dict[roi_number]['dcm_contour'][sop_instance_uid]['contours'].append(contour_data)

    for roi_observations_sequence in rs_ds.RTROIObservationsSequence:
        roi_number = roi_observations_sequence.ReferencedROINumber
        contour_dict[roi_number]['interpreted_type'] = getattr(roi_observations_sequence, "RTROIInterpretedType", "NOTAG")

    return contour_dict


def get_roi_number_to_name(
    rs_ds: Dataset,
    ) -> dict:
    """
    Get roi number to name map from rs_ds.
    Args:
        rs_ds (Dataset): rs_ds
    Returns:
        dict: roi_number_to_name
    """
    roi_number_to_name = {}
    for ssroi in rs_ds.StructureSetROISequence:
        roi_number_to_name[ssroi.ROINumber] = ssroi.ROIName
    return roi_number_to_name
