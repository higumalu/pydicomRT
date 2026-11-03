"""
MR to planning CT: rigid with mutual information, then demons with metric="mind", exported
as a DICOM Deformable Spatial Registration.

Demons compares intensities by default, which only works when the same tissue has the same
intensity in both images. MR and CT do not: here the bone-like insert is bright on CT and
nearly black on MR, the lung-like one the reverse. metric="mind" compares where the edges are
instead (MIND-SSC descriptors). It runs on PyTorch -- a GPU when there is one, the CPU
otherwise -- after `pip install "pydicomrt[gpu]"`.

Run:  python example/05_mr_to_ct_deformable.py
"""

import numpy as np
import SimpleITK as sitk

from _demo_data import make_demo_mr_series, make_demo_series, true_mr_point
from pydicomrt.reg import DeformableSpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import demons_registration, rigid_registration
from pydicomrt.reg.pipeline import compose_transforms
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

try:
    import torch  # noqa: F401  -- metric="mind" needs it
except ImportError:
    print('metric="mind" runs on PyTorch, which is optional: pip install "pydicomrt[gpu]"')
    raise SystemExit(0)

# 1. A planning CT and an MR of the same phantom. The MR has its own contrast, sits 4 mm and
#    3 mm off, and is locally deformed by up to 4 mm -- the answer is known exactly.
GEOMETRY = dict(shape=(24, 96, 96), origin=(-60.0, -60.0, -30.0), pixel_spacing=(1.5, 1.5), slice_spacing=2.5)
SHIFT_MM, BULGE_MM = (4.0, 3.0, 0.0), 4.0
ct_series = make_demo_series(**GEOMETRY)
mr_series = make_demo_mr_series(shift_mm=SHIFT_MM, bulge_mm=BULGE_MM, **GEOMETRY)

ct_image = SimpleITKImageBuilder().from_image_series(ct_series)
mr_image = SimpleITKImageBuilder().from_image_series(mr_series)

# 2. Rigid stage. Mattes mutual information does not need the intensities to match, so this
#    is the same call as for CT to CT. It returns a resampling transform (CT -> MR).
rigid_transform = rigid_registration(ct_image, mr_image)

# 3. Bring the MR onto the CT grid. Pad with 0, not -1000: MR background is 0.
mr_on_ct = sitk.Resample(mr_image, ct_image, rigid_transform, sitk.sitkLinear, 0.0, mr_image.GetPixelID())

# 4. Deformable stage with metric="mind". The schedule, regularisation and smoothing defaults
#    change with the metric; the raw CT numbers and MR signal go in as they are -- MIND does
#    not care about gain or offset, so no normalisation or window clipping is needed.
_, deformable_transform, displacement_field = demons_registration(
    ct_image, mr_on_ct, backend="torch", metric="mind", default_value=0.0,
)

# For comparison: the default intensity metric on the same pair.
_, intensity_transform, _ = demons_registration(
    ct_image, mr_on_ct, backend="torch", default_value=0.0,
)

# 5. Check against the known answer at points around the deformation. The full mapping is
#    rigid(deform(point)): compose_transforms takes the stages in order.
centre = np.array(ct_image.TransformContinuousIndexToPhysicalPoint([(n - 1) / 2 for n in ct_image.GetSize()]))
probes = [centre + np.array([dx, dy, dz]) for dx in (-15, 0, 15) for dy in (-10, 0, 10) for dz in (-10, 0, 10)]
truth = [true_mr_point(p, shift_mm=SHIFT_MM, bulge_mm=BULGE_MM, **GEOMETRY) for p in probes]


def residual_mm(*stages):
    mapping = compose_transforms(list(stages))
    return np.mean([np.linalg.norm(np.subtract(mapping.TransformPoint(p.tolist()), t)) for p, t in zip(probes, truth)])


print(f"distance from the known answer (mean of {len(probes)} probes):")
print(f"  rigid only                 : {residual_mm(rigid_transform):.2f} mm")
print(f"  rigid + demons, intensity  : {residual_mm(rigid_transform, intensity_transform):.2f} mm")
print(f"  rigid + demons, mind       : {residual_mm(rigid_transform, deformable_transform):.2f} mm")

# 6. Export exactly as for CT to CT (example 03): Deformable REG maps fixed -> moving, the
#    residual field first, the forward rigid matrix after it.
post_transform = affine_to_homogeneous_matrix(rigid_transform).astype(np.float32).ravel().tolist()
identity = np.eye(4, dtype=np.float32).ravel().tolist()
reg_ds = (
    DeformableSpatialRegistrationBuilder(ct_series)
    .add_registration(
        moving_series=mr_series,
        vectorial_field_transform=deformable_transform,
        pre_transform=identity,
        post_transform=post_transform,
    )
    .build()
)
# reg_ds.save_as("mr_to_ct_registration.dcm", enforce_file_format=True)

grid = reg_ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0]
print(f"REG grid dimensions          : {list(grid.GridDimensions)}")
print(f"moving series modality       : {mr_series[0].Modality}")
