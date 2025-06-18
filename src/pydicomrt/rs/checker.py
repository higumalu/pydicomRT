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
    status = 0
    first_slice = series_ds_list[0]
    if rs_ds.FrameOfReferenceUID == first_slice.FrameOfReferenceUID:
        status += 1

    try:
        rs_ds_series_uid_list = rs_ds.ReferencedFrameOfReferenceSequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID
        if rs_ds_series_uid_list == first_slice.SeriesInstanceUID:
            status += 1
    except Exception:
        pass

    try:
        sop_instance_uid_list = [ds.SOPInstanceUID for ds in series_ds_list]
        for roi_contour_sequence in rs_ds.ROIContourSequence:
            contour_sequence = roi_contour_sequence.ContourSequence
            for contour in contour_sequence:
                for contour_image_sequence in contour.ContourImageSequence:
                    if contour_image_sequence.ReferencedSOPInstanceUID not in sop_instance_uid_list:
                        status = 0
                        break
    except Exception:
        pass

    return status > 0
