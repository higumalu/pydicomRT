# pydicomRT

Read and write radiotherapy DICOM objects — RTSTRUCT, Spatial and Deformable Registration,
RTDOSE, CT Image — using the array types Python already speaks: `numpy` and `SimpleITK`.

The goal is to keep the DICOM standard out of your way for the common tasks (turn a mask
into contours, turn contours back into a mask, export a registration) while being explicit
about the coordinate conventions that quietly corrupt results when you get them wrong.

**中文說明：[README_zh.md](README_zh.md) · Working with an AI coding agent? [AGENTS.md](AGENTS.md)**

This page is the tour. For the task-by-task manual — every parameter, every error message,
and the 0.8 → 0.9 migration table — see **[docs/user-guide.md](docs/user-guide.md)**. For
the per-function reference, generated from the docstrings, see
**[docs/api-reference.md](docs/api-reference.md)**.

---

## Install

```bash
pip install pydicomrt
```

From source:

```bash
git clone https://github.com/higumalu/pydicomRT.git
cd pydicomRT
pip install -e ".[dev]"
```

Python ≥ 3.10, with `pydicom`, `numpy`, `SimpleITK`, `opencv-python`, `scipy`. The library
is written against pydicom 3.

---

## What you can do with it

| Task | Entry point |
|---|---|
| Load and sort a DICOM image series | `utils.load_sorted_image_series` |
| DICOM series → `sitk.Image` | `utils.sitk_transform.SimpleITKImageBuilder` |
| Build an RTSTRUCT from masks | `rs.RTStructBuilder` |
| RTSTRUCT → 3D masks | `rs.rtstruct_to_masks` |
| Validate an RTSTRUCT / REG / dose / CT | `rs.check_rtstruct_iod`, `reg.check_spatial_reg_iod`, `dose.check_rtdose_iod`, `ct.check_ct_iod` |
| Rigid / deformable registration | `reg.method.rigid_registration`, `reg.method.demons_registration` |
| Export a registration to DICOM REG | `reg.SpatialRegistrationBuilder`, `reg.DeformableSpatialRegistrationBuilder` |
| Read a DICOM REG | `reg.get_spatial_registrations`, `reg.get_deformable_registrations` |
| Build an RTDOSE / read one back | `dose.RTDoseBuilder`, `dose.get_dose_image` |
| Emit CT slices from a volume | `ct.CTBuilder` |
| Denoise, N4 bias correction | `utils.sitk_image_process` |

---

## Conventions you need to know

These four are where things go wrong silently — no exception, a perfectly normal-looking
image, and only the physical coordinates are off. Read this section before the examples.

### 0. Every builder reads the same way

The five object builders share one vocabulary, so learning any of them teaches the rest:

```python
Builder(reference)      # patient/study context comes from the images
    .set_something(...)  # single-valued configuration
    .add_something(...)  # repeatable content
    .build()             # -> FileDataset, or list[FileDataset] for CT
```

| Builder | Constructed with | Content added by |
|---|---|---|
| `rs.RTStructBuilder` | image series | `.add_roi(mask, name, ...)` |
| `dose.RTDoseBuilder` | reference dataset | `.set_dose_grid(...)`, `.add_referenced_plan(...)` |
| `ct.CTBuilder` | image series | `.set_volume(...)`, `.set_plane(...)` |
| `reg.SpatialRegistrationBuilder` | fixed series | `.add_registration(moving, matrix)` |
| `reg.DeformableSpatialRegistrationBuilder` | fixed series | `.add_registration(moving, dvf, pre, post)` |

`set_*` and `add_*` return the builder, so calls chain. `build()` validates and raises rather
than emitting an object that is missing a Type 1 element.

Validation follows the same shape: `check_rtstruct_iod`, `check_spatial_reg_iod`,
`check_deformable_reg_iod`, `check_rtdose_iod` and `check_ct_iod` all return
`{"result": bool, "content": [...]}`.

### 1. Axis order flips three times

| Thing | Order |
|---|---|
| numpy volume stacked from a series | `(slice, row, column)` = `(z, y, x)` |
| `sitk.Image` size / spacing / origin | `(x, y, z)` |
| DICOM `PixelSpacing` | `[row spacing, column spacing]` = `(y, x)` |

`PixelSpacing[0]` is the gap between **rows**, so it belongs on **y**, not **x**. Masks you
hand to `RTStructBuilder.add_roi()` are `(slice, row, column)`.

### 2. Registration transforms point from fixed to moving

Everything this library returns from a registration is a **resampling transform**: it takes
a point on the *fixed* grid and maps it back into the *moving* image. That is the direction
`sitk.Resample` expects, so this works directly:

```python
transform = rigid_registration(fixed_image, moving_image)
registered = sitk.Resample(moving_image, fixed_image, transform, sitk.sitkLinear, -1000.0)
```

A DICOM Spatial REG matrix points the **other way** — moving → fixed. Invert before exporting:

```python
matrix = affine_to_homogeneous_matrix(transform.GetInverse())
```

Both directions produce a file that loads without complaint, so the mistake surfaces as a
shift applied backwards, often only noticed on the treatment machine. See
`example/02_rigid_registration_to_reg.py`.

### 3. Compose transforms, then resample once

`sitk.CompositeTransform([a, b])` applies **b first**. Stage-ordered lists therefore need no
reversing — use `compose_transforms([rigid, deformable])`. Apply the result in a single
`Resample`. Resampling once per stage crops to the fixed grid after each stage, so anything
one stage pushes out of view is replaced by padding and can never be brought back.

### 4. Padding is not free

Resampling outside a source image invents voxels. For CT that value is **-1000** (air).
Padding with 0 wraps the patient in a shell of soft tissue that a similarity metric will
happily try to match. `infer_default_pixel_value()` guesses from the intensity minimum,
which only works on raw data — a window-clipped CT no longer looks like a CT. Pass
`default_value` explicitly for MR, PET, or preprocessed images.

---

## Examples

Every script in `example/` runs standalone on synthetic data — no patient files needed — and
the snippets below are taken from them.

```bash
python example/01_rtstruct_from_mask.py
python example/02_rigid_registration_to_reg.py
python example/03_deformable_registration_to_reg.py
python example/04_rtdose.py
```

### RTSTRUCT from a 3D mask, and back

```python
import numpy as np
from pydicomrt.rs import RTStructBuilder, calc_image_series_affine_mapping, rtstruct_to_masks
from pydicomrt.utils import load_sorted_image_series

series = load_sorted_image_series("path/to/ct")     # sorted along the slice normal

mask = np.zeros((len(series), series[0].Rows, series[0].Columns), dtype=np.uint8)
mask[8:16, 30:60, 30:60] = 1                        # (slice, row, column)

rs_ds = (
    RTStructBuilder(series)
    .add_roi(mask=mask, name="CTV", color=[0, 255, 0])
    .build()
)

rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)

# read it back — the affine and shape come from the image series, not the RTSTRUCT
affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
masks = rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)
recovered = np.asarray(masks["CTV"]["mask_volume"])
```

### Rigid registration exported as DICOM REG

Note the `.GetInverse()` — see convention 2 above.

```python
import numpy as np
from pydicomrt.reg import SpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import rigid_registration
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

fixed_image = SimpleITKImageBuilder().from_image_series(fixed_series)
moving_image = SimpleITKImageBuilder().from_image_series(moving_series)

resampling_transform = rigid_registration(fixed_image, moving_image)   # fixed -> moving

reg_matrix = affine_to_homogeneous_matrix(resampling_transform.GetInverse())  # moving -> fixed
reg_ds = (
    SpatialRegistrationBuilder(fixed_series)           # built with the FIXED series
    .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    .add_registration(moving_series, reg_matrix.astype(np.float32).ravel().tolist())
    .build()
)
reg_ds.save_as("registration.dcm", enforce_file_format=True)
```

The builder validates the matrix before writing it — 16 elements, a `[0, 0, 0, 1]` bottom
row, and a genuine rotation when the type is `RIGID`. Pass `matrix_type="RIGID_SCALE"` or
`"AFFINE"` if scaling or shear is intended. It also appends an identity item for the fixed
Frame of Reference, which is what marks the frame the other matrices are relative to; pass
`build(include_identity=False)` if a downstream system objects.

Reading one back:

```python
from pydicom import dcmread
from pydicomrt.reg import get_spatial_registrations

reg = get_spatial_registrations(dcmread("registration.dcm"))
matrix = reg[fixed_frame_of_reference_uid][moving_frame_of_reference_uid]   # moving -> fixed
```

### Deformable registration exported as DICOM REG

Deformable REG maps fixed → moving. For rigid followed by deformable registration,
export the residual field with an identity pre-matrix and the forward rigid post-matrix.
Do not invert the rigid transform for this object.

```python
import SimpleITK as sitk
from pydicomrt.reg import DeformableSpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import demons_registration, rigid_registration

rigid_transform = rigid_registration(fixed_image, moving_image)
moving_rigid = sitk.Resample(
    moving_image, fixed_image, rigid_transform, sitk.sitkLinear, -1000.0, moving_image.GetPixelID()
)

_, deformable_transform, displacement_field = demons_registration(fixed_image, moving_rigid)

reg_ds = (
    DeformableSpatialRegistrationBuilder(fixed_series)
    .add_registration(
        moving_series=moving_series,
        vectorial_field_transform=deformable_transform,
        pre_transform=np.eye(4).ravel().tolist(),
        post_transform=affine_to_homogeneous_matrix(rigid_transform).ravel().tolist(),
    )
    .build()
)
reg_ds.save_as("deformable_registration.dcm", enforce_file_format=True)
```

### RTDOSE

```python
from pydicomrt.dose import RTDoseBuilder, get_dose_image

dose_ds = (
    RTDoseBuilder(reference_ds)          # patient/study context
    .set_dose_grid(dose_sitk_image)      # geometry + scaled pixel data
    .add_referenced_plan(plan_ds)        # DoseSummationType PLAN makes this Type 1C
    .build()
)
dose_ds.save_as("rtdose.dcm", enforce_file_format=True)

dose_image = get_dose_image(dose_ds)     # back to sitk, DoseGridScaling applied
```

---

## Registration: what is and is not validated

Registration quality is the part most likely to disappoint, so here is the measured state
rather than a claim.

**Verified against a real planning CT / CBCT pair and the TPS registration relating them**
(one pelvis case, single-threaded):

- Series loading matches SimpleITK's own GDCM reader exactly — geometry and voxel values.
- `get_spatial_registrations` reads the TPS registration correctly, in the documented direction.
- `rigid_registration` reproduces the clinical registration to 7.5 mm overall, with the
  superior-inferior axis — the one that used to be worst — under 1 mm.

Measured offsets from the TPS transform, in mm, as `(x, y, z)`, single-threaded:

| `optimizer=` | offset | distance | time |
|---|---|---|---|
| `"regular_step"` (default) | `(+2.8, −7.0, −0.9)` | **7.5** | ~9 min |
| `"gradient_descent"` | `(+1.4, −8.0, −36.7)` | 37.5 | ~1 min |

`"gradient_descent"` can take an estimated first step, overshoot into a flat region, and
have its convergence window then report success — returning an alignment *worse than its own
starting point* without raising anything. `"regular_step"` shrinks its step whenever the
gradient reverses, so it cannot run away. The default is the slow, reliable one; the fast
path stays selectable for previews and pipelines that re-run often.

**Known gaps:**

- **The remaining ~7 mm is almost entirely anterior-posterior, and is not obviously the
  library's error.** Sweeping y around the clinical answer shows full-field-of-view mutual
  information peaking about 10 mm away from it — on this pair, the metric and the TPS
  genuinely disagree. Narrowing it would mean restricting the metric to an ROI, which is not
  implemented. One case is not a validation: check the residual against your own tolerance.
- **Results are only reproducible single-threaded.** ITK reduces metric values in thread
  completion order; thread count moved the answer by more than 15 mm on the real pair. Call
  `sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)` when you need determinism.
- **Only axial acquisitions have been run end to end.** Oblique and anisotropic geometry is
  implemented and unit-tested, but not exercised on a real study.

Two gaps that used to sit here have since closed, and both had been written up with a
plausible explanation first. The superior-inferior error was blamed on the CBCT's limited
field of view; it was the optimizer walking away from a correct position. The pipeline was
blamed on field-of-view mismatch; it was resampling both images onto one grid before
registering, which left the centred initializer nothing to estimate. Worth remembering
before accepting the next "the data is just hard" explanation.

---

## Module layout

| Module | Contents |
|---|---|
| `rs` | `builder`, `add_new_roi`, `make_contour_sequence`, `parser`, `check`, `iod`, `rs_to_volume`, `contour_process_method` |
| `reg` | `builder`, `parser`, `check`, `iod`, `type_transform` |
| `reg.method` | `rigid`, `demons`, `soft_demons`, `bspline`, `common` |
| `reg.pipeline` | `registration_pipeline`, `compose_transforms`, `preprocessing` |
| `dose` | `builder`, `check`, `iod`, `sitk_transform` |
| `ct` | `builder`, `check`, `iod` |
| `utils` | `image_series_loader`, `coordinate_transform`, `sitk_transform`, `validate_dcm_info`, `sitk_image_process` |

Each package re-exports its public API, so `from pydicomrt.rs import ...` is enough — you
should not need to reach into submodules.

See [docs/user-guide.md](docs/user-guide.md) for the full manual,
[docs/api-reference.md](docs/api-reference.md) for every function's parameters, and
[docs/architecture.md](docs/architecture.md) for data flow and extension points.

UID roots: RTSTRUCT and RTDOSE honour the `DICOM_UID_PREFIX` environment variable; the
registration and CT builders take `set_uid_prefix()`.

---

## Development

```bash
pip install -e ".[dev]"
pytest                  # ~295 tests, ~34 s
pytest -m slow          # real-data registrations, ~22 min
```

Tests build their own DICOM series and phantoms in memory — no patient data required.
Builders live in `test/synthetic.py`, fixtures in `test/conftest.py`.

### Testing against real clinical data

`test/real_data/` is an optional layer that runs against a real planning CT, an
on-treatment CBCT, and the TPS registration relating them. It skips cleanly when the data is
absent. Lay the directory out as:

```
<testdata>/CTSIM/*.dcm          planning CT series
<testdata>/CBCT/*.dcm           on-treatment CBCT series
<testdata>/registration.dcm     Spatial Registration produced by the TPS
```

```bash
PYDICOMRT_TESTDATA=/path/to/testdata pytest test/real_data
PYDICOMRT_TESTDATA=/path/to/testdata pytest test/real_data -m slow
```

Keep the data out of the repository — `testdata/` is gitignored.

### If you are adding tests

- Assert on the **mapping** (where a transform sends probe points), not on raw transform
  parameters. One mapping has many parameter representations, and a centred transform's
  `GetTranslation()` is not its offset.
- Geometry changes need a case with **anisotropic** spacing or a **rotated** orientation.
  Square-pixel axial fixtures pass whether or not the code is correct.
- `conftest.py` pins SimpleITK to one thread, for the reason given above.

---

## Reference

- [SimpleITK](https://simpleitk.org/) · [pydicom](https://pydicom.github.io/)
- [RT-Utils](https://github.com/qurit/rt-utils) · [PlatiPy](https://github.com/pyplati/platipy)

## License

MIT — see [LICENSE](LICENSE).

## Author

Higumalu (higuma.lu@gmail.com)
