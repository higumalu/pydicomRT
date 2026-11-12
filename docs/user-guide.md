# pydicomRT User Guide

Task-oriented reference for pydicomRT 0.9. [README.md](../README.md) is the tour; this is
the manual you come back to. For how the library is put together internally, see
[architecture.md](architecture.md).

**中文版：[user-guide_zh.md](user-guide_zh.md)** · **Per-function reference:
[api-reference.md](api-reference.md)**

**0.9 renames most of the public API.** If you are coming from 0.8, start at
[Migrating from 0.8](#migrating-from-08).

- [Install](#install)
- [The four conventions](#the-four-conventions)
- [Loading images](#loading-images)
- [RT Structure Sets](#rt-structure-sets)
- [Registration](#registration)
- [RT Dose](#rt-dose)
- [CT Image](#ct-image)
- [Validation](#validation)
- [Choosing and tuning a registration](#choosing-and-tuning-a-registration)
- [Troubleshooting](#troubleshooting)
- [Migrating from 0.8](#migrating-from-08)

---

## Install

```bash
pip install pydicomrt
```

Python ≥ 3.10 and **pydicom ≥ 3.0**. Pydicom 3 is a hard requirement, not a preference: the
builders take their encoding from the transfer syntax in `file_meta`, where pydicom 2
required the `is_little_endian` / `is_implicit_VR` attributes that pydicom 4 removes.

GPU-accelerated demons registration is an optional extra. PyTorch is a multi-gigabyte
install, so it is left out by default and nothing imports it unless you ask for
`backend="torch"`:

```bash
pip install "pydicomrt[gpu]"
```

If pip's default PyTorch wheel is CPU-only on your platform, install a CUDA build of PyTorch
first (see pytorch.org); the extra only adds `torch` as a requirement.

The wheels include a small compiled C++ library, the inner loops of `rigid_registration`'s
default algorithm. It has no dependencies of its own, so it is not an extra: every install
gets it. Installing from source (`pip install -e .`, or a platform without a wheel) builds it
and needs a C++17 compiler; CMake is fetched automatically. Without a compiler the install
still succeeds, and `rigid_registration` runs the same algorithm on numpy, about ten times
slower.

Saving uses the pydicom 3 spelling throughout:

```python
ds.save_as("out.dcm", enforce_file_format=True)   # write_like_original is gone
```

---

## The four conventions

These are where results go wrong without an exception — the file loads, the image looks
normal, and only the physical coordinates are off. The README covers them in full; this is
the short form.

**1. Axis order flips three times.** A numpy volume stacked from a series is
`(slice, row, column)`. A `sitk.Image` reports `(x, y, z)`. DICOM `PixelSpacing` is
`[row spacing, column spacing]`, i.e. `(y, x)`. Masks handed to `RTStructBuilder.add_roi()`
are `(slice, row, column)`.

**2. Registration transforms point fixed → moving.** Everything the library returns from a
registration is a *resampling* transform, which is the direction `sitk.Resample` wants. A
DICOM Spatial REG matrix points the other way, so invert before exporting to that IOD.
Deformable REG keeps the fixed → moving direction.

**3. Compose transforms, then resample once.** `sitk.CompositeTransform([a, b])` applies
`b` first, so stage-ordered lists need no reversing. Use `compose_transforms()`.

**4. Padding is not free.** Resampling outside a source image invents voxels; for CT that
value is `-1000`. Pass `default_value` explicitly for MR, PET or window-clipped images,
because the inference reads the intensity minimum and a clipped CT no longer looks like one.

---

## Loading images

```python
from pydicomrt.utils import load_sorted_image_series, image_series_to_sitk_image

series = load_sorted_image_series("/path/to/CT")        # a directory, or a list of paths
image = image_series_to_sitk_image(series)              # -> sitk.Image, HU applied
```

`load_sorted_image_series` sorts along the slice normal. Do not assume filename order or
`InstanceNumber`: an unsorted series still produces a valid-looking volume, just flipped in
z, and every contour and registration built on it inherits the flip.

Rescale slope and intercept are applied on conversion, so the `sitk.Image` is in real units.

Need the pieces separately, or a volume of your own on a series' geometry:

```python
from pydicomrt.utils import SimpleITKImageBuilder, parse_image_series

volume, origin, spacing, direction = parse_image_series(series)

image = (
    SimpleITKImageBuilder()
    .set_volume(my_array)          # (slice, row, column)
    .set_origin(origin)
    .set_spacing(spacing)
    .set_direction(direction)
    .build()
)

# or borrow geometry from an image you already have
image = SimpleITKImageBuilder().from_reference_image(my_array, reference_image)
```

---

## RT Structure Sets

### Mask to RTSTRUCT

```python
import numpy as np
from pydicomrt.rs import RTStructBuilder

mask = np.zeros((len(series), series[0].Rows, series[0].Columns), dtype=np.uint8)
mask[8:16, 30:60, 30:60] = 1                    # (slice, row, column)

rs_ds = (
    RTStructBuilder(series)
    .add_roi(mask=mask, name="CTV", color=[0, 255, 0], interpreted_type="CTV")
    .add_roi(mask=cord_mask, name="SpinalCord")
    .build()
)
rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)
```

ROI numbers are assigned in sequence unless you pass `number=`. The mask shape is checked
against the series, because a mismatched mask would otherwise be sliced against the wrong
images and produce contours in the wrong place rather than an error.

`add_roi()` also takes `description` and `contour_config` (the contour extraction tuning —
noise size and low-pass ratio; see `rs.builder.DEFAULT_CONTOUR_CONFIG`).

### RTSTRUCT to masks

```python
from pydicomrt.rs import calc_image_series_affine_mapping, rtstruct_to_masks

affine, shape = calc_image_series_affine_mapping(series)   # from the images, not the RTSTRUCT
masks = rtstruct_to_masks(rs_ds, affine, shape)

ctv = np.asarray(masks["CTV"]["mask_volume"])              # (slice, row, column)
```

The affine and the volume shape come from the image series. An RTSTRUCT stores contour
points in patient coordinates and does not know the grid you want them on.

### Reading structure metadata

```python
from pydicomrt.rs import get_roi_names, get_contours, is_rtstruct_matching_series

get_roi_names(rs_ds)                    # {roi_number: roi_name}
get_contours(rs_ds)                     # contours grouped by ROI and referenced image
is_rtstruct_matching_series(rs_ds, series)
```

`get_roi_names` keys are `pydicom.valuerep.IS`, an int subclass. `names[1]` works and
`names["1"]` raises `KeyError`, despite the repr showing `{'1': ...}`.

### Contours you already have

If your contours are already points in patient coordinates rather than a mask:

```python
RTStructBuilder(series).add_roi_from_contours(
    contours={sop_instance_uid: {"sop_class_uid": ..., "contours": [[x, y, z, ...], ...]}},
    name="CTV",
).build()
```

---

## Registration

### Rigid, exported as DICOM REG

```python
import numpy as np
from pydicomrt.reg import SpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import rigid_registration
from pydicomrt.utils import image_series_to_sitk_image

fixed_image = image_series_to_sitk_image(fixed_series)
moving_image = image_series_to_sitk_image(moving_series)

transform = rigid_registration(fixed_image, moving_image)      # fixed -> moving

registered = sitk.Resample(
    moving_image, fixed_image, transform, sitk.sitkLinear, -1000.0
)

# DICOM wants moving -> fixed, so invert
matrix = affine_to_homogeneous_matrix(transform.GetInverse())
reg_ds = (
    SpatialRegistrationBuilder(fixed_series)                   # built with the FIXED series
    .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    .add_registration(moving_series, matrix.ravel().tolist())
    .build()
)
reg_ds.save_as("registration.dcm", enforce_file_format=True)
```

The builder validates the matrix before writing: 16 elements, a `[0, 0, 0, 1]` bottom row,
finite values, and — when the type is `RIGID` — a genuine rotation. `RIGID_SCALE` is checked
too, and rejects shear. Declare `matrix_type="AFFINE"` if shear is intended rather than
letting it through under a narrower label.

`build()` appends an identity item for the fixed Frame of Reference, which is what marks
which frame the other matrices are relative to. Pass `build(include_identity=False)` if a
downstream system objects. It then runs the object's own IOD checker and raises
`Invalid Spatial Registration: ...` (or `Invalid Deformable Registration: ...`) rather than
handing back a dataset that fails conformance. The constructor, for its part, requires every
reference image to share one Frame of Reference, so a series list assembled from two studies
is rejected before any registration is added.

Both REG objects carry four Type 1 equipment attributes — Manufacturer, Manufacturer's
Model Name, Device Serial Number and Software Versions — and a TPS refuses the import if
any is missing or empty. They default to `pydicomRT` (Software Versions also records the
installed version). If the receiving system expects your site's own values, set them with
`set_equipment()`; arguments you leave out keep their default:

```python
builder.set_equipment(manufacturer="My Hospital", device_serial_number="SN-001")
```

### Deformable, exported as DICOM REG

Deformable REG maps fixed → moving. For rigid followed by deformable registration,
export the residual field with the forward rigid transform as the post-matrix. Do not
invert it for this object. `pre_transform` and `post_transform` both default to `None`,
which omits the sequence entirely — there is no need to pass an identity matrix.

```python
import SimpleITK as sitk
from pydicomrt.reg import DeformableSpatialRegistrationBuilder, affine_to_homogeneous_matrix
from pydicomrt.reg.method import demons_registration, rigid_registration

rigid = rigid_registration(fixed_image, moving_image)
moving_rigid = sitk.Resample(
    moving_image, fixed_image, rigid, sitk.sitkLinear, -1000.0, moving_image.GetPixelID()
)

registered, deform_transform, field = demons_registration(fixed_image, moving_rigid)

reg_ds = (
    DeformableSpatialRegistrationBuilder(fixed_series)
    .add_registration(
        moving_series=moving_series,
        vectorial_field_transform=deform_transform,
        post_transform=affine_to_homogeneous_matrix(rigid).ravel().tolist(),
    )
    .build()
)
```

`demons_registration` and `bspline_registration` both return
`(registered_image, transform, deformation_field)`, so they are interchangeable.

### Reading an existing REG

```python
from pydicom import dcmread
from pydicomrt.reg import get_spatial_registrations, get_deformable_registrations

reg = get_spatial_registrations(dcmread("registration.dcm"))
matrix = reg[fixed_frame_uid][moving_frame_uid]        # moving -> fixed, 16 floats

deformable = get_deformable_registrations(dcmread("deformable.dcm"))
field = deformable[0]["DeformableRegistrationGrid"]["VectorGridData"]   # (z, y, x, 3)
```

**This applies to Spatial REG matrices only.** To use one as a resampling transform, invert
it — Spatial REG stores moving → fixed and `sitk.Resample` wants fixed → moving:

```python
matrix_4x4 = np.array(matrix, dtype=float).reshape(4, 4)
inverse = np.linalg.inv(matrix_4x4)
transform = sitk.AffineTransform(3)
transform.SetMatrix(inverse[:3, :3].ravel().tolist())
transform.SetTranslation(inverse[:3, 3].tolist())
```

A deformable object's pre/post matrices already run fixed → moving, so use them as they
come — inverting them is the mistake this direction split exists to prevent. Absent pre/post
sequences are returned as identity matrices, and only grid-bearing entries are returned.

### The pipeline

```python
from pydicomrt.reg.pipeline import registration_pipeline

registered, rigid, deformable, field = registration_pipeline(
    fixed_image, moving_image,
    perform_rigid=True,
    perform_deformable=True,
    preprocess_config={"rigid": {"window_clip": [-200, 800]}},
)
```

The returned image is always on the `fixed_image` grid, carries the moving image's pixel
type, and is a single resample of the untouched moving image through the composed transform
— so preprocessing never leaks into the output, and the image agrees exactly with the
transforms it hands back.

With no preprocessing configured the pipeline does exactly what calling
`rigid_registration` directly does.

---

## RT Dose

```python
from pydicomrt.dose import RTDoseBuilder, get_dose_image

dose_ds = (
    RTDoseBuilder(reference_ds)              # a planning CT slice, or the plan
    .set_dose_grid(dose_sitk_image)          # values in Gy
    .add_referenced_plan(plan_ds)
    .build()
)
dose_ds.save_as("rtdose.dcm", enforce_file_format=True)

dose_image = get_dose_image(dose_ds)         # back to sitk, DoseGridScaling applied
```

`DoseSummationType` defaults to `"PLAN"`, which makes Referenced RT Plan Sequence Type 1C —
so a plan reference is required. `check_rtdose_iod()` reports it if missing.
`set_dose_grid()` also takes `dose_type` (default `"PHYSICAL"`), `dose_summation_type`,
`dose_units` and `dose_grid_scaling`.

Dose is stored as scaled unsigned 32-bit integers. The default scaling of 1e-7 covers about
429 Gy; beyond that the builder raises the scaling automatically rather than wrapping.
Negative dose is rejected, since unsigned storage would turn it into a very large positive.

---

## CT Image

```python
from pydicomrt.ct import CTBuilder

slices = (
    CTBuilder(reference_series)              # patient/study context
    .set_volume(sitk_image)                  # values treated as HU
    .set_plane("AXIAL")                      # or CORONAL, SAGITTAL
    .build()
)
for index, slice_ds in enumerate(slices):
    slice_ds.save_as(f"ct_{index:04d}.dcm", enforce_file_format=True)
```

Planes are named for the plane each slice lies in, so the axis the slices stack along is the
one *not* in the name: axial stacks along z, coronal along y, sagittal along x.

Output is marked `DERIVED\SECONDARY` and gets its own Series Instance UID, so it cannot be
mistaken for the acquired series. Values are clipped to what unsigned 16-bit can hold with
the -1024 intercept, i.e. -1024 to 64511 HU.

`set_volume_from_array(volume, origin, spacing, direction)` takes a numpy volume plus
geometry instead. `set_copy_all_attributes(True)` inherits every attribute from the
reference series rather than building a fresh dataset — useful for acquisition parameters,
but it inherits everything else the source carried too.

---

## Validation

Every modality has a checker with the same shape:

```python
from pydicomrt.rs import check_rtstruct_iod
from pydicomrt.reg import check_spatial_reg_iod, check_deformable_reg_iod
from pydicomrt.dose import check_rtdose_iod
from pydicomrt.ct import check_ct_iod

result = check_rtstruct_iod(rs_ds)
# {"result": True, "content": []}          -- conformant
# {"result": False, "content": ["Missing InstanceNumber in root", ...]}
```

These check required fields against the object's CIOD, and `check_rtdose_iod` additionally
checks the conditional requirement most Dose Summation Types place on Referenced RT Plan
Sequence.

The two REG checkers go a little further than field presence: they verify the SOP Class,
that transformation matrices are orthonormal for the type they declare, and that a
deformation grid is internally consistent — dimensions, spacing, orthonormal orientation,
and a `VectorGridData` byte count matching `GridDimensions`. What none of them can see is
whether the registration points the right *way*: a backwards transform is structurally
perfect.

CT is a single-frame IOD, so `check_ct_iod` takes one slice at a time.

---

## Choosing and tuning a registration

### Which method

| | Use when | Cost |
|---|---|---|
| `rigid_registration` | Position and orientation only | ~1 s for a brain MR/CT with 8 threads |
| `demons_registration` | Deformable, the default choice — same modality (CT↔CT, CT↔CBCT) | fast, forgiving of its parameters |
| `demons_registration(metric="mind")` | Deformable, different modalities (MR→CT) | a few seconds on a GPU; needs `[gpu]` |
| `bspline_registration` | Deformable, when you need a parametric transform | ~40x demons for the same accuracy |
| `demons_with_soft_mask` | Deformable, restricted to a structure | as demons |

**Prefer demons for deformable work.** On a synthetic phantom recovering a 2.11 mm
deformation, demons scores 0.80 mm in 0.4 s against B-spline's 0.82 mm in 17 s. Reach for
B-spline when you specifically need a parametric transform — control points you can store or
manipulate — rather than a displacement field. Demons can also run on a GPU; see below.

**Demons compares intensities**, so by default it only works when the same tissue has the
same intensity in both images. For MR to CT, use `metric="mind"`; see
[MR to CT](#mr-to-ct-metricmind).

### Rigid

```python
transform = rigid_registration(fixed_image, moving_image)
```

The library's own implementation of elastix's strategy: Mattes mutual information on fresh
random points every iteration, adaptive stochastic gradient ascent over a four-level
smoothing pyramid, then a deterministic finish on one fixed set of 500 000 points. It uses
the exact derivative of the linear interpolant, and its schedule is fixed.

Measured against elastix's default rigid (mean error, mm; failure = more than 5 mm):

| | `rigid_registration` | elastix |
|---|---:|---:|
| RIRE brain MR → CT, fiducial gold standard (6 MR sequences) | 1.38 | 1.32 |
| RIRE, agreement of T1/T2/PD from one session (16 patients) | 1.33 | 1.31 |
| RIRE, MR header moved by up to ±10°, ±20 mm (10 starts) | 1.24 | 1.35 |
| SynthRAD2023 brain MR → CT, header moved (9 pairs) | 0.53 | 1.00 |
| Pelvis CT/CBCT against the TPS registration (1 pair) | 5.1 | 16.7 |
| Time, brain MR/CT, 8 threads | ~1.2 s | ~1.2 s |

Both fail on one RIRE patient (about 8 mm between that session's sequences), probably real
motion between them. On the pelvis pair it takes about 15 s single-threaded.

- `backend`: `"auto"` (default) uses the compiled core, or numpy — about ten times slower —
  in an installation built from source without a C++ compiler. Both draw the same random
  points; rounding still tips the odd point into a neighbouring histogram bin and the
  stochastic steps carry that forward, so the two backends' answers differ by up to about
  0.2 mm.
- `seed`: the random points. Seeds 0–2 on the pelvis pair gave 5.09–5.16 mm.
- `threads`: the compiled core's threads, SimpleITK's global default when `None`. The answer
  does not depend on it.

`algorithm="itk"` keeps the previous default, SimpleITK's `ImageRegistrationMethod`, with
its parameters (`histogram_bins`, `optimizer`, `shrink_factors`, …). It is only
bit-reproducible single-threaded, since ITK reduces the metric in thread completion order.

### Demons on a GPU

`backend="torch"` runs the same multiscale demons on PyTorch, on a GPU when one is available.
It needs the `[gpu]` extra (see [Install](#install)) and returns the same three objects:

```python
registered, transform, field = demons_registration(
    fixed_image, moving_rigid, backend="torch"       # device="cuda:1" picks a card
)

# or through the pipeline
registration_pipeline(fixed_image, moving_image,
                      deformable_kwargs={"backend": "torch"})
```

It is a port of the ITK filter, not a reimplementation of the idea — boundary handling,
smoothing kernels and the early stop follow the ITK sources — and the tests hold it to the
SimpleITK backend on anisotropic, rotated and mismatched grids, where the two agree to
thousandths of a millimetre. On a real 418 × 417 × 126 planning CT (22 M voxels; RTX 5090,
24-core CPU) with default settings:

| backend | time |
|---|---|
| `"simpleitk"`, one thread (the `ncores` default) | ~63 s |
| `"simpleitk"`, `ncores=24` | ~8 s |
| `"torch"`, CUDA | ~0.7 s |

Both scale linearly with iterations, so the ratio holds as you raise them.

- **Memory** follows the finest pyramid level, at roughly 180 bytes per voxel — about 4 GB
  for that 22 M-voxel CT. A 1024 × 1024 × 200 volume would need ~39 GB at full resolution;
  with `isotropic_resample=True, resolution_staging=(8, 4, 1.5)` it ran in 5.2 GB and 3.7 s
  (36.7 s on 24 threads) at the same accuracy, and the field still comes back on the fixed
  grid. Running out raises `torch.OutOfMemoryError` with an estimate.
- **RAM with `backend="simpleitk"`** follows the finest level too, at roughly 150 bytes per
  voxel — about 6 GB for a 512 × 512 × 160 CT, most of it inside ITK's demons filter. On a
  machine with less, coarsen the finest level the same way: at 2 mm
  (`isotropic_resample=True, resolution_staging=(8, 4, 2)`) that CT needs 2.5 GB, of which
  2 GB is the returned transform and field, which are full-resolution `float64`.
- **Spend extra iterations at the finest level.** On known synthetic deformations of a real
  CT, raising every level's count made the result worse — 1.15 mm at the default
  `(10, 10, 10)`, 4.45 mm at `(100, 100, 100)` — while raising only the last improved it:
  0.55 mm at `(10, 10, 100)`. That holds for both backends; the GPU just makes the long
  schedules affordable.
- **Image edges are unstable on either backend.** Where anatomy is cut off by the volume
  edge, demons is chaotic: perturbing the input by a thousandth of a HU moves the SimpleITK
  result there by millimetres. Expect the two backends to differ there by as much; inside the
  volume they agree.
- Repeated runs on one device give identical results. The torch backend is 3D only,
  interpolates linearly, and ignores `ncores` and `verbose`. With no CUDA device it runs on
  the CPU and logs a warning.

### MR to CT: `metric="mind"`

Demons' force is the intensity difference between the images, which means nothing between MR
and CT — bone is bright on one and dark on the other. On undeformed brain MR/CT pairs the
default metric moves the brain by 9 mm on average, and so does `registration_pipeline`, whose
deformable stage is demons. `metric="mind"` replaces the intensities with MIND-SSC
descriptors (Heinrich et al., 2013): at every voxel, how alike its neighbouring patches are.
That depends on where the edges are, not on which side is brighter, so it is the same in MR
and CT wherever both show the same boundary.

```python
rigid = rigid_registration(ct_image, mr_image)            # mutual information: MR is fine
mr_on_ct = sitk.Resample(mr_image, ct_image, rigid, sitk.sitkLinear, 0.0)  # pad MR with 0

registered, transform, field = demons_registration(
    ct_image, mr_on_ct, backend="torch", metric="mind", default_value=0.0,
)

# or through the pipeline
registration_pipeline(ct_image, mr_image, default_value=0.0,
                      deformable_kwargs={"backend": "torch", "metric": "mind"})
```

`example/05_mr_to_ct_deformable.py` runs this end to end on a synthetic pair with a known
answer and writes the Deformable REG; the export is the same as for CT to CT.

- **It needs `backend="torch"`** (the `[gpu]` extra; a CUDA card is optional — the CPU works,
  slowly).
- **Feed it raw values.** MIND is unchanged by gain and offset, so CT numbers and MR signal go
  in as they are. Do not set `window_clip` in the pipeline for an MR pair: it applies one
  window to both images, and a CT window flattens the MR.
- **Its defaults differ**: `resolution_staging=(4, 2, 1)`, `iteration_staging=(100, 100, 200)`,
  `regularization_kernel_mm=3.0`, `smoothing_sigma_factor=0.5`. The intensity defaults
  ((8, 4, 1), 1.5 mm) each raised the invented motion to 2–4 mm with this metric. Anything you
  pass explicitly is used as given.
- **Memory**: about 360 bytes per voxel of the finest level, twice the intensity metric. For a
  full-resolution planning CT, register at 1–1.5 mm: `isotropic_resample=True,
  resolution_staging=(6, 3, 1.5)`. The field still comes back on the fixed grid.
- It does not reproduce SimpleITK: it carries one field through all pyramid levels instead
  of composing a refinement per level, which is what keeps a finer level able to correct a
  coarser one.

Measured on nine public brain MR/CT pairs (SynthRAD2023: radiotherapy planning CT and 3D
T1, rigidly aligned by the organisers), with known deformations applied to the CT:

| | brain moved on an undeformed pair | error, known 2.7 mm deformation | time |
|---|---|---|---|
| no deformable stage | 0 | 2.74 mm | — |
| `metric="intensity"` | 9.1 mm | 9.3 mm | ~1 s |
| `bspline_registration(metric="mutual_information")` (one case) | 6.6 mm | 6.6 mm | 13 min |
| `metric="mind"` | 0.62 mm | 0.86 mm | ~9 s |

Mutual information works for rigid registration — six parameters — but a dense or B-spline
deformation can raise it by distorting what lies outside the brain, and does. Part of the
0.62 mm on undeformed pairs is real (MR distortion, posture; the pairs were only rigidly
aligned): on a synthetic MR made from the CT itself, with no residual at all, `"mind"` moves
the brain by about 0.1 mm.

**For a brain target, rigid is usually the answer.** In those measurements the deformable
stage moved the brain by well under a millimetre; the real deformation is in the neck and
skull base, where an 8° flexion was only a third recovered. And a deformable result moves
real anatomical change — tumour growth, oedema, a resection cavity — out of sight. Review it
before contouring on it. MIND sees edges, so 2D MR series with 3–5 mm slices give it little
to go on along the slice direction; register the 3D T1 to the CT and bring other series of
the same MR study across through their rigid alignment to that T1.

### B-spline

The parameter that governs whether a B-spline fit means anything is the ratio of metric
samples to transform coefficients. Too few samples and the optimizer drives the metric to
near zero by contorting the grid, reports convergence, and returns a deformation that is
unconstrained wherever the metric did not look. Nothing raises.

`describe_sampling()` reports that ratio per level, and `bspline_registration` warns when a
level is underdetermined.

**The sampling rate is not a speed dial.** It looks like one:

| sampling | error | time |
|---|---|---|
| `None` (dense, the default) | 0.89 mm | 37 s |
| `REGULAR` at 100% | 1.50 mm | 47 s |
| `REGULAR` at 10% | 2.04 mm | 4 s |
| `RANDOM` at 1% | 1.64 mm | 1 s |

`REGULAR` at 100% samples the same voxels as `None` and still scores 1.50 mm, so the gap is
the derivative path ITK takes once a sampled point set is in play — not the sample count.
Every sampled configuration lands between 1.5 and 2.0 mm whatever the rate. Leave it at
`None`.

To trade accuracy for time, use the dials that vary smoothly:

| dial | effect |
|---|---|
| `ncores` | near-linear; 111 s on one thread to 19 s on sixteen |
| `number_of_iterations` | mean error by cap: 5 → 1.09, 10 → 1.01, **20 → 0.82**, 30 → 0.88, 50 → 1.07 mm |
| `initial_grid_spacing` | coarser is faster; 96 mm → 1.00 mm at 11 s, 64 mm → 0.89 mm at 37 s |
| `resolution_staging` | fewer levels is faster; `(2, 1)` → 0.95 mm at 11 s |

Note that the iteration cap has an optimum rather than a ceiling: past roughly twenty the
fit starts chasing noise, so raising it costs time *and* accuracy.

`optimizer="LBFGSB"` requires `grid_scale_factors=None` — ITK sizes the optimizer scales
once, and a scaled grid changes the coefficient count between levels.

### What is and is not validated

Registration accuracy has been measured against a TPS registration on **one** real
planning CT / CBCT pelvis pair, and against known deformations on synthetic phantoms.

- Series loading matches SimpleITK's own GDCM reader exactly, geometry and voxel values.
- `rigid_registration` reproduces the clinical registration to 5.1 mm overall,
  `(+1.6, −4.7, −1.1)` mm. Most of the residual is anterior-posterior, where full-field-of-view
  mutual information peaks away from the clinical answer — so on that pair the metric and the
  TPS genuinely disagree. It was also measured on public brain MR/CT against a fiducial gold
  standard (RIRE: one patient, 1990s 4 mm slices).
- **Deformable registration has never been run against a clinical ground truth.** The
  B-spline measurements above recover a B-spline deformation with a B-spline, which is the
  most favourable possible test. Real anatomical deformation is not B-spline shaped.
- `metric="mind"` was measured on nine public brain MR/CT pairs against *known, synthetic*
  deformations applied to real images — 3D T1 only, no T2/FLAIR, no tumour change.

One case is not a validation. Check the residual against your own tolerance.

---

## Troubleshooting

### Loading

| Message | Cause |
|---|---|
| `No DICOM files found in ...` | The directory has no `*.dcm`, or you passed a file path where a directory was expected |
| `Error reading DICOM file ...` | Not DICOM, or truncated |
| `ImageOrientationPatient [...] is not two orthogonal unit vectors` | Corrupt or non-conformant orientation. Raised as `InvalidImageOrientationError`, a `ValueError` subclass |

### Builders

| Message | Fix |
|---|---|
| `no ROI added; call add_roi() before build()` | The builder had nothing to write |
| `no registration added; call add_registration() before build()` | Same, for REG |
| `no dose grid set; call set_dose_grid() before build()` | Same, for dose |
| `no volume set; call set_volume() before build()` | Same, for CT |
| `mask shape (...) does not match the image series (...)` | The mask is `(slice, row, column)` and must match the series |
| `roi_number N was already added` | Two ROIs claimed the same number |
| `the bottom row of the matrix must be [0, 0, 0, 1]` | The matrix is row-major; translation goes in the last column |
| `matrix_type is RIGID but the upper 3x3 is not orthonormal` | Declare `RIGID_SCALE` or `AFFINE` |
| `RIGID_SCALE requires orthogonal, non-zero axes; use AFFINE for shear` | The matrix has shear, which `RIGID_SCALE` does not cover |
| `Invalid Spatial Registration: ...` / `Invalid Deformable Registration: ...` | `build()`'s own conformance check failed; the message lists what is missing |
| `reference images must share one non-empty FrameOfReferenceUID` | The series handed to the REG builder mixes Frames of Reference |
| `DICOM deformation grids require a 3D field with three components` | The displacement field is not 3D-with-3-components |
| `DICOM deformation grid direction must be right-handed and orthonormal` | The field's direction matrix is a reflection or non-orthonormal |
| `reference_ds is not an RT Plan` | `add_referenced_plan()` needs an RT Plan (SOP class `...481.5`) |
| `dose grid holds negative values` | Stored dose is unsigned; clip or offset first |

### Registration

| Symptom | Likely cause |
|---|---|
| Warning: "the fit is underdetermined" | B-spline has more coefficients than metric samples. Raise `sampling_rate`, coarsen `initial_grid_spacing`, or drop the coarsest level |
| `optimizer must be one of (...)` | Check spelling; the list is `bspline.OPTIMIZERS` |
| `optimizer='LBFGSB' cannot be combined with grid_scale_factors` | Pass `grid_scale_factors=None` |
| Result differs between runs | With `algorithm="itk"`, ITK reduces the metric in thread completion order: pin to one thread. The default does not depend on it; check `seed` |
| Registration returns something worse than doing nothing | If using `algorithm="itk", optimizer="gradient_descent"`, switch to the default |
| `ImportError: backend='torch' needs PyTorch` | `pip install "pydicomrt[gpu]"` |
| `ImportError: backend='cpp' needs the compiled core` | pydicomrt was built from source without a C++ compiler. Install a wheel, or reinstall with a C++17 compiler; `backend="auto"` runs on numpy meanwhile |
| `torch.OutOfMemoryError` from `backend="torch"` | The finest pyramid level is too big for the card. Coarsen it with `isotropic_resample=True, resolution_staging=(8, 4, 1.5)` (`(6, 3, 1.5)` for `metric="mind"`) |
| `MemoryError`, or the process is killed, during `demons_registration` | The finest pyramid level does not fit in RAM (~150 bytes per voxel). Coarsen it with `isotropic_resample=True, resolution_staging=(8, 4, 2)`; the field still comes back on the fixed grid |
| `metric='mind' runs on backend='torch' only` | Add `backend="torch"`; MIND has no SimpleITK implementation |
| MR registered to CT is distorted, several mm of motion where none should be | The default intensity metric on a multimodal pair. Use `metric="mind"` |

### Output looks shifted or stretched

Almost always one of the four conventions:

- **Shifted the wrong way** — a transform used in the wrong direction. Registrations return
  fixed → moving; DICOM Spatial REG stores moving → fixed; Deformable REG maps fixed → moving.
- **Stretched or squashed in-plane** — row/column spacing swapped somewhere. `PixelSpacing`
  is `[row, column]`, SimpleITK spacing is `(x, y, z)`.
- **Mirrored in z** — the series was not sorted along the slice normal.
- **A soft-tissue shell around the patient** — padded with 0 instead of -1000.

---

## Migrating from 0.8

0.9 renames the public API for consistency and keeps no aliases. Behaviour changed too;
read [What changed behaviourally](#what-changed-behaviourally) even if your imports still
resolve.

### Object construction

All five builders now share one shape: `Builder(reference).set_*(...).add_*(...).build()`.

```python
# 0.8
rs_ds = create_rtstruct_dataset(series)
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "target")
rs_ds = add_contour_sequence_from_mask3d(rs_ds, series, 1, mask)

# 0.9
rs_ds = RTStructBuilder(series).add_roi(mask=mask, name="CTV", color=[0, 255, 0]).build()
```

```python
# 0.8
dose_ds = generate_base_dataset()
dose_ds = cp_information_from_ds(dose_ds, reference_ds)
dose_ds = add_dose_grid_to_ds(dose_ds, dose_image)

# 0.9
dose_ds = RTDoseBuilder(reference_ds).set_dose_grid(dose_image).add_referenced_plan(plan).build()
```

```python
# 0.8
slices = CTBuilder(series).build_from_sitk_image(image, "AXIAL")

# 0.9
slices = CTBuilder(series).set_volume(image).set_plane("AXIAL").build()
```

The 0.8 functions still exist in their submodules and still work; they are the building
blocks the builders are assembled from. They are no longer the documented entry point.

### Renames

| 0.8 | 0.9 |
|---|---|
| `check_rs_iod` | `check_rtstruct_iod` |
| `check_s_reg_iod` | `check_spatial_reg_iod` |
| `check_ds_reg_iod` | `check_deformable_reg_iod` |
| `get_contour_dict` | `get_contours` |
| `get_roi_number_to_name` | `get_roi_names` |
| `rtstruct_to_mask_dict` | `rtstruct_to_masks` |
| `get_spatial_reg_dict` | `get_spatial_registrations` |
| `get_deformable_reg_list` | `get_deformable_registrations` |
| `get_dose_sitk_image` | `get_dose_image` |
| `get_dose_array_spacing` | `get_dose_spacing` |
| `parse_ds_list` | `parse_image_series` |
| `ds_list_to_sitk_image` | `image_series_to_sitk_image` |
| `sort_ds_list` | `sort_image_series` |
| `SimpleITKImageBuilder.from_ds_list` | `.from_image_series` |
| `SimpleITKImageBuilder.from_ref_sitk_image` | `.from_reference_image` |
| `SimpleITKImageBuilder.from_dcms_dir` | `.from_dicom_directory` |
| `add_rigid_registration` / `add_deformable_registration` | `add_registration` |
| `n4bfc` | `n4_bias_field_correction` |

### Moved and removed modules

| 0.8 | 0.9 |
|---|---|
| `rs.checker` | `rs.check` |
| `rs.rs_ds_iod` (`RT_STRUCTURE_SET_IOD`) | `rs.iod` (`RTSTRUCT_IOD`) |
| `reg.s_reg_ds_iod` (`SPATIAL_REGSITRATION_IOD`) | `reg.iod` (`SPATIAL_REGISTRATION_IOD`) |
| `reg.ds_reg_ds_iod` | `reg.iod` (`DEFORMABLE_SPATIAL_REGISTRATION_IOD`) |
| `dose.dose_ds_iod` (`RT_DOSE_IOD`) | `dose.iod` (`RTDOSE_IOD`) |
| `ct.ct_ds_iod` | `ct.iod` |
| `rs.packer` | removed — an unused older copy of `make_contour_sequence` + `add_new_roi` |
| `utils.rs_from_altas` | removed — the file was empty |

Series parameters are named `image_series`, or `fixed_series` / `moving_series` where the
role matters. `optimiser` is spelled `optimizer`.

### What changed behaviourally

Your imports may still resolve while the output differs.

**Geometry fixes.** These were wrong and are now right, which means output changes for
affected data. Square-pixel axial CT is unaffected.

- `PixelSpacing` was copied straight into SimpleITK's `(x, y)`, swapping the in-plane
  spacings for anisotropic pixels. A series written from a 1.5 × 3.0 mm volume read back as
  3.0 × 1.5.
- Direction matrices were built with the axis vectors as rows rather than columns, so
  oblique and rotated series landed mirrored.
- `ImageOrientationPatient` was written from the first six entries of a row-major direction
  flattening rather than the first two columns.

**Pixel encoding.** `CTBuilder` and the dose builder wrote unsigned data while declaring
`PixelRepresentation = 1` (signed), so values past the signed midpoint read back negative —
dense implants on CT, dose above about 215 Gy. Both now declare unsigned.

**Registration.** `rigid_registration` defaults to a regular-step optimizer; the previous
unbounded gradient descent could return a result worse than its own initialisation. On a
real CT/CBCT pair this took the error from 37.5 mm to 7.5 mm. `registration_pipeline` no
longer resamples both images onto a shared union-extent grid before the rigid stage, which
had reduced the centred initializer to a no-op: 194 mm to 7.5 mm on the same pair.
`bspline_registration` returns three values rather than two, matching demons.

**Deformable REG direction — check any file you have already written.** Deformable Spatial
Registration maps fixed → moving (PS3.3 C.20.3.1.1), the *same* direction registration
functions return and the opposite of a Spatial REG matrix. Earlier versions documented and
exported it as moving → fixed, inverting the rigid stage into the pre-matrix. Every
deformable object written before 0.9.0 therefore points the wrong way and needs re-exporting;
the correct form is the residual field in the grid with the forward rigid matrix as
`post_transform`. Spatial REG is unaffected — it really is moving → fixed.

Two supporting fixes landed with it. `affine_to_homogeneous_matrix` now folds a non-zero
rotation centre into the matrix offset instead of silently dropping it, so a centred rigid
transform no longer needs `to_centre_free_affine` first. And `get_spatial_registrations`
composes every item in a Matrix Sequence in DICOM order, where it previously read only the
last one.

**DICOM conformance.** The REG builders were missing `InstanceNumber`, `ContentDate` and
`ContentTime` (all Type 1), wrote Registration Type Code Sequence at the wrong nesting
level, and omitted the identity item for the fixed Frame of Reference. Files written by 0.8
are readable but not conformant on those points. Both REG builders now run their checker in
`build()`, so a malformed object fails where it is constructed rather than downstream.

**pydicom floor raised to 3.0.** See [Install](#install).
