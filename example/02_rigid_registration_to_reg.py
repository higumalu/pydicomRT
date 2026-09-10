"""
Rigid registration between two series, exported as a DICOM Spatial Registration (REG).

The important detail is the inversion in step 4. Getting it wrong produces a REG file that
loads fine and shifts the image the wrong way.

Run:  python example/02_rigid_registration_to_reg.py
"""

import numpy as np

from _demo_data import make_demo_series
from pydicomrt.reg import (
    SpatialRegistrationBuilder,
    affine_to_homogeneous_matrix,
    get_spatial_registrations,
)
from pydicomrt.reg.method import rigid_registration
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

# 1. Load both series. The moving series here sits 6 mm off in x and 4 mm in y.
KNOWN_SHIFT = (6.0, 4.0, 0.0)
GEOMETRY = dict(shape=(32, 128, 128), pixel_spacing=(1.0, 1.0), slice_spacing=2.0)
fixed_series = make_demo_series(**GEOMETRY)
moving_series = make_demo_series(shift_mm=KNOWN_SHIFT, **GEOMETRY)

# 2. Convert to SimpleITK. Geometry (origin, spacing, direction) comes from the DICOM tags,
#    so the two images stay in their true relative positions in patient space.
fixed_image = SimpleITKImageBuilder().from_image_series(fixed_series)
moving_image = SimpleITKImageBuilder().from_image_series(moving_series)

# 3. Register. The result is a RESAMPLING transform: it maps a point on the fixed grid
#    back into the moving image, which is the direction sitk.Resample expects.
resampling_transform = rigid_registration(fixed_image, moving_image)
print(f"known shift            : {KNOWN_SHIFT}")
print(f"recovered (fixed->moving): {tuple(round(v, 2) for v in resampling_transform.GetTranslation())}")

# 4. Invert before export. A DICOM REG matrix maps the referenced (moving) Frame of
#    Reference INTO the registered RCS, i.e. moving -> fixed -- the opposite direction.
reg_transform = resampling_transform.GetInverse()
matrix = affine_to_homogeneous_matrix(reg_transform).astype(np.float32).ravel().tolist()

# 5. Build the REG. The builder is constructed with the FIXED series (it supplies patient,
#    study and frame-of-reference context); the moving series is passed per registration.
reg_ds = (
    SpatialRegistrationBuilder(fixed_series)
    .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    .add_registration(moving_series, matrix)
    .build()
)
# reg_ds.save_as("registration.dcm", enforce_file_format=True)

# 6. Read it back and confirm the stored direction
parsed = get_spatial_registrations(reg_ds)
fixed_frame = fixed_series[0].FrameOfReferenceUID
moving_frame = moving_series[0].FrameOfReferenceUID
stored = np.array(parsed[fixed_frame][moving_frame], dtype=float).reshape(4, 4)

print(f"stored in REG (moving->fixed): {tuple(round(v, 2) for v in stored[:3, 3])}")

# The stored matrix must be the inverse of the resampling transform. Comparing the whole
# matrix rather than just negating the translation, since any rotation makes the inverse
# translation differ from the plain negation.
expected = np.linalg.inv(affine_to_homogeneous_matrix(resampling_transform))
print("stored matrix is the inverse of the resampling transform:",
      np.allclose(stored, expected, atol=1e-4))
