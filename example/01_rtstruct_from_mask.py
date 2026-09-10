"""
Build an RTSTRUCT from a 3D mask, then read the mask back out.

Run:  python example/01_rtstruct_from_mask.py
"""

import numpy as np

from _demo_data import make_demo_series
from pydicomrt.rs import (
    RTStructBuilder,
    calc_image_series_affine_mapping,
    check_rtstruct_iod,
    rtstruct_to_masks,
)

# 1. Load the reference image series. Slices must be sorted along the slice normal --
#    load_sorted_image_series() does that for you when reading from disk.
series = make_demo_series()

# 2. Build the mask. Shape is (slices, rows, columns) -- the same axis order as
#    numpy arrays stacked from the series, i.e. (z, y, x).
mask = np.zeros((len(series), series[0].Rows, series[0].Columns), dtype=np.uint8)
rows, cols = np.mgrid[0:series[0].Rows, 0:series[0].Columns]
disc = ((cols - 48) ** 2 + (rows - 48) ** 2) < 20 ** 2
mask[8:16][:, disc] = 1

# 3. Every builder in the package reads the same way: set_* / add_* / build().
#    ROI numbers are assigned in sequence unless you pass one.
rs_ds = (
    RTStructBuilder(series)
    .add_roi(mask=mask, name="CTV", color=[0, 255, 0], description="demo target")
    .build()
)

# 4. Save. pydicom 3 uses enforce_file_format; write_like_original was removed.
# rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)

# 5. Read it back. The affine and volume shape come from the image series, not the RTSTRUCT.
affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
masks = rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)
recovered = np.asarray(masks["CTV"]["mask_volume"])

overlap = 2 * np.logical_and(mask > 0, recovered > 0).sum() / ((mask > 0).sum() + (recovered > 0).sum())
print(f"contours written : {len(rs_ds.ROIContourSequence[0].ContourSequence)}")
print(f"IOD check        : {check_rtstruct_iod(rs_ds)['result']}")
print(f"round-trip Dice  : {overlap:.4f}")
