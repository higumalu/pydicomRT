"""
Write an RTDOSE from a SimpleITK dose grid, then read it back.

Run:  python example/04_rtdose.py
"""

import numpy as np
import SimpleITK as sitk
from pydicom.dataset import Dataset
from pydicom.uid import RTPlanStorage, generate_uid

from _demo_data import make_demo_series
from pydicomrt.dose import RTDoseBuilder, check_rtdose_iod, get_dose_image

series = make_demo_series()
reference_ds = series[0]

# A dose grid does not have to match the CT grid; it carries its own geometry.
dose_array = np.zeros((20, 40, 40), dtype=np.float32)
z, y, x = np.mgrid[0:20, 0:40, 0:40]
dose_array = 60.0 * np.exp(-(((x - 20) ** 2 + (y - 20) ** 2) / 200.0 + ((z - 10) ** 2) / 60.0))

dose_image = sitk.GetImageFromArray(dose_array)
dose_image.SetSpacing((2.5, 2.5, 3.0))
dose_image.SetOrigin((-50.0, -50.0, -25.0))

# DoseSummationType "PLAN" (the default) makes Referenced RT Plan Sequence Type 1C, so the
# dose has to say which plan it belongs to. Stand-in here; use your real plan dataset.
plan_ds = Dataset()
plan_ds.SOPClassUID = RTPlanStorage
plan_ds.SOPInstanceUID = generate_uid()
plan_ds.file_meta = Dataset()
plan_ds.file_meta.MediaStorageSOPClassUID = RTPlanStorage

dose_ds = (
    RTDoseBuilder(reference_ds)
    .set_dose_grid(dose_image)
    .add_referenced_plan(plan_ds)
    .build()
)

# dose_ds.save_as("rtdose.dcm", enforce_file_format=True)

# Read back. get_dose_image applies DoseGridScaling and rebuilds the geometry.
recovered = get_dose_image(dose_ds)
recovered_array = sitk.GetArrayFromImage(recovered)

print(f"modality            : {dose_ds.Modality}")
print(f"dose units          : {dose_ds.DoseUnits}")
print(f"IOD check           : {check_rtdose_iod(dose_ds)['result']}")
print(f"dose grid scaling   : {dose_ds.DoseGridScaling:g}")
print(f"grid size           : {recovered.GetSize()}  (SimpleITK reports x, y, z)")
print(f"spacing round trip  : {tuple(round(v, 4) for v in recovered.GetSpacing())}")
print(f"origin round trip   : {tuple(round(v, 4) for v in recovered.GetOrigin())}")
print(f"max dose written    : {dose_array.max():.4f} Gy")
print(f"max dose read back  : {recovered_array.max():.4f} Gy")
