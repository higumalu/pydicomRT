"""
Check RTSTRUCT dataset Class Information Object Definition (IOD)

Author: Higumalu
Date: 2025-06-13
"""

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


def check_rs_iod(rs_ds):
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

