"""
Rigid pre-alignment followed by demons deformable registration, exported as a DICOM
Deformable Spatial Registration.

Run:  python example/03_deformable_registration_to_reg.py
"""

import numpy as np
import SimpleITK as sitk

from _demo_data import make_demo_series
from pydicomrt.reg import DeformableSpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import demons_registration, rigid_registration
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

GEOMETRY = dict(shape=(24, 96, 96), pixel_spacing=(1.5, 1.5), slice_spacing=2.5)
fixed_series = make_demo_series(**GEOMETRY)
moving_series = make_demo_series(shift_mm=(4.0, 3.0, 0.0), **GEOMETRY)

fixed_image = SimpleITKImageBuilder().from_image_series(fixed_series)
moving_image = SimpleITKImageBuilder().from_image_series(moving_series)

# 1. Rigid stage. Both transforms below are resampling transforms (fixed -> moving).
rigid_transform = rigid_registration(fixed_image, moving_image)

# 2. Bring the moving image into the fixed frame before the deformable stage. -1000 is the
#    CT air value; use 0 for MR or PET, otherwise the padding invents a soft-tissue shell.
moving_rigid = sitk.Resample(
    moving_image, fixed_image, rigid_transform, sitk.sitkLinear, -1000.0, moving_image.GetPixelID()
)

# 3. Demons. Returns the resampled image, a displacement-field transform, and the field.
_, deformable_transform, displacement_field = demons_registration(
    fixed_image,
    moving_rigid,
    resolution_staging=(4, 2),
    iteration_staging=(10, 10),
)

# 4. Deformable REG maps fixed -> moving, unlike Spatial REG. The pipeline maps
#    rigid(deform(point)), so put the forward rigid transform AFTER the residual field.
post_transform = affine_to_homogeneous_matrix(
    rigid_transform
).astype(np.float32).ravel().tolist()
identity = np.eye(4, dtype=np.float32).ravel().tolist()

reg_ds = (
    DeformableSpatialRegistrationBuilder(fixed_series)
    .add_registration(
        moving_series=moving_series,
        vectorial_field_transform=deformable_transform,
        pre_transform=identity,
        post_transform=post_transform,
    )
    .build()
)
# reg_ds.save_as("deformable_registration.dcm", enforce_file_format=True)

grid = reg_ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0]
vectors = np.frombuffer(grid.VectorGridData, dtype="<f4").reshape(-1, 3)

print(f"rigid stage (fixed->moving): {tuple(round(v, 2) for v in rigid_transform.GetTranslation())}")
print(f"grid dimensions            : {list(grid.GridDimensions)}")
print(f"grid resolution (mm)       : {[round(float(v), 3) for v in grid.GridResolution]}")
print(f"displacement vectors       : {vectors.shape[0]} x 3 float32")
print(f"max displacement (mm)      : {np.abs(vectors).max():.2f}")
