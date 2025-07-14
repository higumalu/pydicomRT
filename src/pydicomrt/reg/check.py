from pydicom.dataset import Dataset
from .ds_reg_ds_iod import DEFORMABLE_SPATIAL_REGSITRATION_IOD
from .s_reg_ds_iod import SPATIAL_REGSITRATION_IOD
from pydicomrt.utils.validate_dcm_info import check_iod


def check_ds_reg_iod(rs_ds: Dataset) -> dict:
    result_dict = {
        "result": True,
        "content": []
    }
    failed_item_list = []
    failed_item_list = check_iod(ds=rs_ds, config_map=DEFORMABLE_SPATIAL_REGSITRATION_IOD)
    if len(failed_item_list) > 0:
        result_dict["result"] = False
        result_dict["content"] = failed_item_list
    return result_dict


def check_s_reg_iod(rs_ds: Dataset) -> dict:
    result_dict = {
        "result": True,
        "content": []
    }
    failed_item_list = []
    failed_item_list = check_iod(ds=rs_ds, config_map=SPATIAL_REGSITRATION_IOD)
    if len(failed_item_list) > 0:
        result_dict["result"] = False
        result_dict["content"] = failed_item_list
    return result_dict
