"""
Check RTSTRUCT dataset Class Information Object Definition (IOD)

Author: Higumalu
Date: 2025-06-13
"""
from pydicom.dataset import Dataset

from pydicomrt.utils.validate_dcm_info import check_iod, ValidationError
from .rs_ds_iod import RT_STRUCTURE_SET_IOD


def contour_data_validator(value):
    if hasattr(value, '__len__'):
        if len(value) == 0:
            raise ValidationError("Contour data is empty")
        elif len(value) < 9:
            raise ValidationError("Contour data less then 3 points")
        elif len(value) % 3 != 0:
            raise ValidationError("Contour data length is not multiple of 3")
    else:
        raise ValidationError("Contour data must be a sequence type that supports len()")


RS_VALIDATORS_MAP = {
    "ContourDataValidator": contour_data_validator,
}


def check_rs_iod(rs_ds: Dataset) -> dict:
    result_dict = {
        "result": True,
        "content": []
    }
    failed_item_list = []
    failed_item_list = check_iod(ds=rs_ds, config_map=RT_STRUCTURE_SET_IOD, validators=RS_VALIDATORS_MAP, path="")

    if len(failed_item_list) > 0:
        result_dict["result"] = False
        result_dict["content"] = failed_item_list

    return result_dict

# --------------------------------------------------------------------------------------------------------------------- #

def is_rtstruct_matching_series(
    rs_ds: Dataset,
    series_ds_list: list) -> bool:
    """
    Check if the RTSTRUCT dataset is matching the series dataset
    :param rs_ds: RTSTRUCT dataset
    :param series_ds_list: list of series dataset
    :return: True if matching, False otherwise
    """
    if not series_ds_list:
        return False

    status = 0
    first_slice = series_ds_list[0]

    rs_for_uid = getattr(rs_ds, "FrameOfReferenceUID", None)
    series_for_uid = getattr(first_slice, "FrameOfReferenceUID", None)
    if rs_for_uid is not None and series_for_uid is not None and rs_for_uid == series_for_uid:
        status += 1

    referenced_for_seq = getattr(rs_ds, "ReferencedFrameOfReferenceSequence", None)
    if referenced_for_seq:
        try:
            rt_referenced_series_seq = getattr(referenced_for_seq[0], "RTReferencedSeriesSequence", None)
            if rt_referenced_series_seq:
                rs_series_uid = getattr(rt_referenced_series_seq[0], "SeriesInstanceUID", None)
                series_instance_uid = getattr(first_slice, "SeriesInstanceUID", None)
                if rs_series_uid is not None and series_instance_uid is not None and rs_series_uid == series_instance_uid:
                    status += 1
        except (AttributeError, IndexError, TypeError):
            pass

    sop_instance_uid_list = [
        uid for uid in (getattr(ds, "SOPInstanceUID", None) for ds in series_ds_list)
        if uid is not None
    ]
    for roi_contour_sequence in getattr(rs_ds, "ROIContourSequence", []):
        for contour in getattr(roi_contour_sequence, "ContourSequence", []):
            for contour_image in getattr(contour, "ContourImageSequence", []):
                referenced_sop_uid = getattr(contour_image, "ReferencedSOPInstanceUID", None)
                if referenced_sop_uid is not None and referenced_sop_uid not in sop_instance_uid_list:
                    status = 0
                    break

    return status > 0
