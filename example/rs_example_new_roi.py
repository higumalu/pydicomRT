import numpy as np

from pydicomrt.rs.make_contour_sequence import add_contour_sequence_from_dcm_ctr_dict, add_contour_sequence_from_mask3d
from pydicomrt.rs.add_new_roi import create_roi_into_rs_ds
from pydicomrt.rs.builder import create_rtstruct_dataset
from pydicomrt.utils.image_series_loader import load_sorted_image_series

from pydicomrt.rs.parser import get_roi_number_to_name, get_contour_dict

ds_list = load_sorted_image_series("./example/data/RV_002/CT")     # load image series
rs_ds = create_rtstruct_dataset(ds_list)                        # create empty RTSTRUCT dataset
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "CTV")  # create ROI into RTSTRUCT dataset

mask = np.zeros((len(ds_list), 512, 512))
mask[100:200, 100:400, 100:400] = 1
mask[120:180, 200:300, 200:300] = 0

rs_ds = add_contour_sequence_from_mask3d(rs_ds, ds_list, 1, mask)   # add contour sequence into RTSTRUCT dataset

ctr_dict = get_contour_dict(rs_ds)  # get contour dict from RTSTRUCT dataset
roi_map = get_roi_number_to_name(rs_ds)  # get roi map from RTSTRUCT dataset
print(roi_map)

dcm_ctr_dict = ctr_dict[1]['dcm_contour']
rs_ds = create_roi_into_rs_ds(rs_ds, [255, 255, 0], 2, "PTV", "PTV")      # create ROI into RTSTRUCT dataset
rs_ds = add_contour_sequence_from_dcm_ctr_dict(rs_ds, ds_list, 2, dcm_ctr_dict)   # add contour sequence into RTSTRUCT dataset

rs_ds.save_as("./test/data/RV_002/RS/test.dcm", write_like_original=False)  # save RTSTRUCT dataset
