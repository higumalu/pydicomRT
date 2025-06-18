
import numpy as np
import pydicom
from pydicom.dataset import Dataset

from pydicomrt.rs.checker import check_rs_iod, is_rtstruct_matching_series
from pydicomrt.rs.make_contour_sequence import add_contour_sequence_from_mask3d
from pydicomrt.rs.add_new_roi import create_roi_into_rs_ds
from pydicomrt.rs.builder import create_rtstruct_dataset
from pydicomrt.utils.image_series_loader import load_sorted_image_series

ds_list = load_sorted_image_series("example/data/RV_002/CT")     # load image series
rs_ds = create_rtstruct_dataset(ds_list)
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "CTV")  # create ROI into RTSTRUCT dataset

mask = np.zeros((len(ds_list), 512, 512))
mask[100:200, 100:400, 100:400] = 1
mask[120:180, 200:300, 200:300] = 0

rs_ds = add_contour_sequence_from_mask3d(rs_ds, ds_list, 1, mask)   # add contour sequence into RTSTRUCT dataset
rs_ds.ROIContourSequence[0].ContourSequence[0] = Dataset()
ctr_data = rs_ds.ROIContourSequence[0].ContourSequence[1].ContourData
rs_ds.ROIContourSequence[0].ContourSequence[1].ContourData = [1, 1, 1]
series_uid = rs_ds.ROIContourSequence[0].ContourSequence[1].ContourImageSequence[0].ReferencedSOPInstanceUID
roi_inter_type = rs_ds.RTROIObservationsSequence[0].RTROIInterpretedType
print(type(roi_inter_type))
print(len(roi_inter_type))

rs_ds = pydicom.dcmread("example/data/RV_002/RS/4021.dcm")
print(check_rs_iod(rs_ds))
print(is_rtstruct_matching_series(rs_ds, ds_list))
