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
    for roi_contour_sequence in getattr(rs_ds, "ROIContourSequence", []):
        roi_number = getattr(roi_contour_sequence, "ReferencedROINumber", None)
        if roi_number is None:
            continue

        roi_color = getattr(roi_contour_sequence, "ROIDisplayColor", [255, 255, 255])
        roi_name = number_name_map.get(roi_number, f"ROI_{roi_number}")
        contour_dict[roi_number] = {
            'color': roi_color,
            'name': roi_name,
            'dcm_contour': {},
        }

        for contour_sequence in getattr(roi_contour_sequence, "ContourSequence", []):
            try:
                contour_image_sequence = getattr(contour_sequence, "ContourImageSequence", None)
                if not contour_image_sequence:
                    continue
                contour_image = contour_image_sequence[0]
                sop_instance_uid = getattr(contour_image, "ReferencedSOPInstanceUID", None)
                sop_class_uid = getattr(contour_image, "ReferencedSOPClassUID", None)
                contour_data = getattr(contour_sequence, "ContourData", None)
                if sop_instance_uid is None or contour_data is None:
                    continue

                if sop_instance_uid not in contour_dict[roi_number]['dcm_contour']:
                    contour_dict[roi_number]['dcm_contour'][sop_instance_uid] = {
                        'sop_class_uid': sop_class_uid,
                        'contours': [],
                    }

                contour_dict[roi_number]['dcm_contour'][sop_instance_uid]['contours'].append(contour_data)
            except (AttributeError, IndexError, TypeError):
                continue

    for roi_observations_sequence in getattr(rs_ds, "RTROIObservationsSequence", []):
        roi_number = getattr(roi_observations_sequence, "ReferencedROINumber", None)
        if roi_number is None or roi_number not in contour_dict:
            continue
        contour_dict[roi_number]['interpreted_type'] = getattr(
            roi_observations_sequence, "RTROIInterpretedType", "NOTAG"
        )

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
    for ssroi in getattr(rs_ds, "StructureSetROISequence", []):
        roi_number = getattr(ssroi, "ROINumber", None)
        if roi_number is None:
            continue
        roi_name = getattr(ssroi, "ROIName", None)
        roi_number_to_name[roi_number] = roi_name if roi_name is not None else f"ROI_{roi_number}"
    return roi_number_to_name
