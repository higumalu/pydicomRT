# AGENTS.md

Working notes for coding agents. Humans: [README.md](README.md) is the friendlier entry
point; this file is the dense version, and everything in it is true for you too.

## What this library is

`pydicomrt` converts between DICOM radiotherapy objects and the array/image types Python
already speaks — `numpy` and `SimpleITK`. It builds and parses RTSTRUCT, Spatial and
Deformable Registration, RTDOSE, and CT Image objects.

It is a **conversion and construction** library. It does not manage a database, talk to a
PACS, or make clinical decisions.

## Repository map

```
src/pydicomrt/
  <package>/  each modality follows the same layout:
              builder.py  construct · parser.py  read · check.py + iod.py  validate
  rs/       RTSTRUCT: RTStructBuilder, mask <-> contour, ROI parsing
  reg/      Registration: two REG builders/parsers, algorithms, pipeline
    type_transform.py  sitk transform <-> DICOM matrix / deformation grid
    method/     rigid, demons, bspline, soft_demons, common  (SimpleITK algorithms)
    pipeline/   registration_pipeline + preprocessing helpers
  dose/     RTDOSE: RTDoseBuilder, sitk_transform.py (get_dose_image)
  ct/       CT Image: CTBuilder
  utils/    series loading, coordinate transforms, SimpleITK conversion, check_iod
    sitk_image_process/  bilateral/median denoise, N4 bias correction — importable
                         only from this subpackage, not re-exported by utils
test/
  synthetic.py        in-memory DICOM/phantom builders (no patient data)
  conftest.py         fixtures; pins SimpleITK to one thread
  test_*_conventions.py  guardrails, not feature tests — see "Conventions the suite enforces"
  real_data/          opt-in tests against real clinical DICOM, skipped when absent
example/              01–04 are verified and standalone; the README snippets come from these.
                      rs_example_*.py are 0.8-style and read `example/data/`, which is not
                      in the repo — they fail on a clean checkout. Do not cite them.
docs/user-guide.md    task-by-task manual: parameters, errors, 0.8 -> 0.9 migration
docs/api-reference.md generated from docstrings — regenerate, never hand-edit
docs/architecture.md  module responsibilities and data flow
README.md, README_zh.md, docs/user-guide_zh.md  the zh files mirror the en ones; a user-
                      facing change belongs in both
```

Nothing is exported from the top-level package — `pydicomrt/__init__.py` holds only a
placeholder `main()`. Import from the subpackages: `from pydicomrt.rs import ...`,
`pydicomrt.reg`, `pydicomrt.reg.method`, `pydicomrt.reg.pipeline`, `pydicomrt.dose`,
`pydicomrt.ct`, `pydicomrt.utils`.

## Commands

```bash
pip install -e ".[dev]"
pytest                                    # ~1000 tests, ~35 s
pytest -m slow                            # 9 full-size registrations, tens of minutes
python example/01_rtstruct_from_mask.py   # 01–04 run standalone
python docs/generate_api_reference.py     # after any public docstring change
```

Counts are parametrised, so the number is per test *case*: 1001 collected by default,
992 passing here plus 9 skips from `test_docstring_conventions.py` opting out of parameter
checks for exception classes and no-argument callables — those skips are expected. Without
`testdata/` another 33 (all of `test/real_data`) skip instead of running. `-m slow` is 8
real-data tests plus one synthetic B-spline accuracy test.

Real-data tests need `PYDICOMRT_TESTDATA` (or a `testdata/` directory). They skip cleanly
when it is absent — a skip there is not a failure.

## Builder vocabulary

Every object builder is the same shape, so one example teaches the rest:

```python
Builder(reference).set_something(...).add_something(...).build()
```

`set_*` takes single-valued configuration, `add_*` takes repeatable content, both return the
builder, and `build()` validates and raises rather than emitting an object missing a Type 1
element. `RTStructBuilder`, `RTDoseBuilder`, `CTBuilder`, `SpatialRegistrationBuilder` and
`DeformableSpatialRegistrationBuilder` all follow it; `build()` returns a `FileDataset`
except for `CTBuilder`, which returns one per slice.

Validation matches: `check_rtstruct_iod`, `check_spatial_reg_iod`, `check_deformable_reg_iod`,
`check_rtdose_iod` and `check_ct_iod` all return `{"result": bool, "content": [...]}`.

Within a package the layout is fixed too: `builder.py` constructs, `parser.py` reads,
`check.py` validates against `iod.py`.

The 0.8 functional API is still exported and still works — `create_rtstruct_dataset`,
`create_roi_into_rs_ds`, `add_contour_sequence_from_mask3d` in `rs`, and
`generate_base_dataset`, `cp_information_from_ds`, `add_dose_grid_to_ds` in `dose`. It is
not deprecated in code, but the builders are what the docs teach; write new code and
examples with the builders. `docs/user-guide.md` has the 0.8 → 0.9 mapping.

## Conventions the suite enforces

Three test files assert on the shape of the code rather than its behaviour, so a change
that reads fine can fail them:

- `test_module_conventions.py` — every builder exposes only `set_*` / `add_*` / `build`,
  every `set_*` / `add_*` is *annotated* to return the builder (the class or a `BuilderT`
  TypeVar), every checker takes exactly one `Dataset` and returns
  `{"result", "content"}`, every modality package has `builder`/`check`/`iod`, and every
  package declares `__all__` with no missing names.
- `test_docstring_conventions.py` — everything reachable through a package's `__all__`,
  public methods included, needs a numpydoc docstring whose `Parameters` match the
  signature. `docs/api-reference.md` is generated from these, so a drifted docstring
  publishes something false. Adding a name to `__all__` therefore also commits you to its
  docstring.
- `test_geometry_conventions.py` — the axis-order and direction-matrix rules below,
  asserted against anisotropic and rotated fixtures.

## Conventions that will bite you

These are the mistakes that do not raise. The code looks right, the image looks right, and
only the physical coordinates are wrong.

### Axis order

| Thing | Order |
|---|---|
| numpy volume from a series | `(slice, row, column)` = `(z, y, x)` |
| `sitk.Image.GetSize()` / `GetSpacing()` / `GetOrigin()` | `(x, y, z)` — reversed |
| DICOM `PixelSpacing` | `[row spacing, column spacing]` = `(y, x)` — reversed again |
| Mask passed to `RTStructBuilder.add_roi()` | `(slice, row, column)` |

`PixelSpacing[0]` is the distance between **rows**, so it belongs on the **y** axis. Use
`get_sitk_spacing()` / `sitk_spacing_to_pixel_spacing()` in `utils.sitk_transform` rather
than indexing by hand.

### Direction matrices

A SimpleITK direction matrix holds each axis direction in a **column**. DICOM
`ImageOrientationPatient` is the row direction followed by the column direction, i.e. the
first two **columns** — not the first six entries of the row-major flattening. Those two
readings coincide only for axis-aligned axial series, which is why mistakes here survive
until someone loads an oblique or gantry-tilted study. Use
`sitk_direction_to_image_orientation_patient()`.

### Transform direction — the one that costs the most

Two opposite conventions are in play:

| Object | Maps |
|---|---|
| Anything returned by `rigid_registration`, `demons_registration`, `registration_pipeline` | **fixed → moving** (a *resampling* transform, what `sitk.Resample` wants) |
| A DICOM Spatial REG matrix | **moving → fixed** |
| A DICOM Deformable REG (grid plus pre/post matrices) | **fixed → moving** |

So **invert before writing to Spatial REG**:

```python
resampling = rigid_registration(fixed_image, moving_image)   # fixed -> moving
matrix = affine_to_homogeneous_matrix(resampling.GetInverse())  # moving -> fixed
builder.add_registration(moving_series, matrix.astype(np.float32).ravel().tolist())
```

For Deformable REG, keep the resampling direction. A rigid-then-deformable pipeline
maps `rigid(deform(point))`: use identity pre-matrix, the residual deformation grid,
and the forward rigid matrix as post-matrix. Do not invert it.

Both directions produce a file that loads without complaint. The only way to tell them
apart is to check that anatomy actually lines up — see
`test/real_data/test_clinical_series.py::TestClinicalTransformDirection`.

### Writing a DICOM REG

Nesting differs between the two registration IODs, and getting it wrong produces a file that
loads while a conformant reader finds nothing:

| | Registration Type Code Sequence lives |
|---|---|
| Spatial Registration (rigid) | inside `MatrixRegistrationSequence`, beside `MatrixSequence` |
| Deformable Spatial Registration | directly in the `DeformableRegistrationSequence` item |

`SpatialRegistrationBuilder.build()` also appends an identity item for the fixed Frame of
Reference (code 125021, "Frame of Reference Identity"; pass `include_identity=False` to
suppress it) — that item is what marks which frame the other
matrices are relative to. Moving items use 125024, "Image Content-based Alignment".

Matrices are validated at build time: 16 elements, `[0, 0, 0, 1]` bottom row, and a genuine
rotation for `RIGID`. Declare `RIGID_SCALE` or `AFFINE` rather than sneaking scale past a
`RIGID` label.

### Composing transforms

`sitk.CompositeTransform([a, b])` evaluates **back to front**: `b` applies first. Stage-order
lists (rigid first, then deformable) therefore pass through unchanged — use
`compose_transforms()` from `pydicomrt.reg.pipeline` rather than hand-rolling it, and apply
the composed transform in **one** `Resample` call. Resampling once per stage crops to the
fixed grid after every stage, permanently discarding anything a later stage would have
brought back into view.

### Padding values

Resampling outside a source image invents voxels. For CT that value is `-1000` (air); pad
with `0` and you have wrapped the patient in a shell of soft tissue that the similarity
metric will try to match. `infer_default_pixel_value()` picks it from the intensity
minimum, which only works on **raw** intensities — after window clipping a CT no longer
reads as a CT. Pass `default_value` explicitly for MR, PET, or preprocessed images.

### Slice ordering

`load_sorted_image_series()` sorts along the slice normal. Do not assume filename order or
`InstanceNumber`. An unsorted series still produces a valid-looking volume — just flipped
in z, with every contour and registration built on it inheriting the flip.

### pydicom version

The code requires pydicom 3. `save_as(..., write_like_original=False)` is gone; use
`save_as(path, enforce_file_format=True)`. The builders take the encoding from the transfer
syntax in `file_meta` rather than setting the `is_little_endian` / `is_implicit_VR`
attributes pydicom 4 removes — which is what puts the floor at 3.0, since pydicom 2 needed
those attributes.

### Writing pixel data

Pixel Representation must match what you actually wrote. Both image builders got this wrong
in the same way: `uint16`/`uint32` data declared signed, so every value past the signed
midpoint read back negative. Nothing raises — it only shows up on dense implants (CT) or
very high dose.

| Builder | Storage | Pixel Representation |
|---|---|---|
| `CTBuilder` | `uint16`, `RescaleIntercept = -1024` | `0` |
| `RTDoseBuilder` | `uint32`, scaled by `DoseGridScaling` | `0` |

Clip *after* applying the intercept offset, not before, or the clip bound itself overflows.

## Known limitations

Be straight with users about these rather than presenting registration output as final.

- **`rigid_registration` on CT ↔ CBCT is not clinically validated.** Measured against a TPS
  registration on one real pelvis pair it is `(+2.8, −7.0, −0.9)` mm out, 7.5 mm total.
  Almost all of that is anterior-posterior, and full-FOV mutual information peaks about
  10 mm from the clinical answer on this pair, so the metric and the TPS genuinely disagree
  there. Narrowing it means restricting the metric to an ROI, which is not implemented.
  One case is not a validation.
- **`optimizer="gradient_descent"` is a fast path, not an equivalent.** It can return an
  alignment worse than its own initialisation and still report convergence — 37.5 mm out on
  the pair above, against 7.5 mm for the `"regular_step"` default. It costs about a minute
  instead of nine on full-size volumes. Do not switch to it to make tests faster.
- **Registration is only reproducible single-threaded.** ITK reduces metric values in thread
  completion order; on the real CT/CBCT pair the thread count moved the answer by more than
  15 mm. Set `sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)` when you need
  determinism.
- **Only axial acquisitions are exercised end to end.** The geometry conventions above are
  implemented and unit-tested for oblique and anisotropic data, but no real oblique study
  has been run through the library.

## If you change something

- Registration behaviour: assert on the **mapping** (where the transform sends probe
  points), never on raw transform parameters — one mapping has many parameter
  representations, and a centred transform's `GetTranslation()` is not its offset.
- Geometry code: add a case with **anisotropic** spacing or a **rotated** orientation.
  Square-pixel axial fixtures pass whether or not the code is correct.
- A public docstring or a new name in `__all__`: re-run
  `python docs/generate_api_reference.py` in the same change. The reference is generated,
  so editing `docs/api-reference.md` by hand loses the edit at the next run.
- Anything user-facing: `README.md` / `docs/user-guide.md` have `_zh` mirrors. Update both,
  or they drift apart silently.
- Anything in `example/`: `01`–`04` are the documentation. Run them.
- Do not commit patient data. `testdata/` is gitignored; keep it that way.
