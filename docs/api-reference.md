# pydicomRT API Reference

Every public function and class, with parameters, return values and the behaviour
that is not visible from the signature.

This page is **generated from the docstrings** by `docs/generate_api_reference.py`.
Edit the docstrings and re-run it; edits made here are lost on the next run. The
same text is what `help(...)` prints and what an IDE shows on hover.

For task-oriented instructions see [user-guide.md](user-guide.md); for how the
library is put together see [architecture.md](architecture.md).

---

## Contents

**[`utils`](#utils)** — Series loading, geometry, and SimpleITK conversion.

- [`load_sorted_image_series`](#utilsload-sorted-image-series) <sub>func</sub> — Read a DICOM image series and return it sorted along the slice normal.
- [`sort_image_series`](#utilssort-image-series) <sub>func</sub> — Sort slices by their position along the slice normal, ascending.
- [`InvalidImageOrientationError`](#utilsinvalidimageorientationerror) <sub>class</sub> — ImageOrientationPatient (0020,0037) is not two orthogonal unit vectors.
- [`get_slice_directions`](#utilsget-slice-directions) <sub>func</sub> — Split ImageOrientationPatient into its three unit direction vectors.
- [`get_slice_position`](#utilsget-slice-position) <sub>func</sub> — Position of a slice along the slice normal.
- [`get_spacing_between_slices`](#utilsget-spacing-between-slices) <sub>func</sub> — Mean distance between slice centres, along the slice normal.
- [`get_pixel_to_patient_transformation_matrix`](#utilsget-pixel-to-patient-transformation-matrix) <sub>func</sub> — 4x4 affine mapping pixel indices to patient coordinates (mm).
- [`get_patient_to_pixel_transformation_matrix`](#utilsget-patient-to-pixel-transformation-matrix) <sub>func</sub> — 4x4 affine mapping patient coordinates (mm) to pixel indices.
- [`apply_transformation_to_3d_points`](#utilsapply-transformation-to-3d-points) <sub>func</sub> — Apply a 4x4 affine to an `(n, 3)` array of points.
- [`SimpleITKImageBuilder`](#utilssimpleitkimagebuilder) <sub>class</sub> — Assemble a `sitk.Image` from a volume plus geometry.
- [`image_series_to_sitk_image`](#utilsimage-series-to-sitk-image) <sub>func</sub> — Convert a DICOM image series straight to a `sitk.Image`.
- [`parse_image_series`](#utilsparse-image-series) <sub>func</sub> — Read volume and geometry out of a DICOM image series.
- [`resample_to_reference_image`](#utilsresample-to-reference-image) <sub>func</sub> — Resample the source image onto the reference image grid (identity transform).
- [`sitk_direction_to_image_orientation_patient`](#utilssitk-direction-to-image-orientation-patient) <sub>func</sub> — Convert a SimpleITK direction matrix into DICOM ImageOrientationPatient (0020,0037).
- [`sitk_spacing_to_pixel_spacing`](#utilssitk-spacing-to-pixel-spacing) <sub>func</sub> — Convert SimpleITK `(x, y, z)` spacing into DICOM PixelSpacing (0028,0030).
- [`check_iod`](#utilscheck-iod) <sub>func</sub> — Check a DICOM dataset against an IOD description.

**[`rs`](#rs)** — RT Structure Sets: masks to contours and back.

- [`RTStructBuilder`](#rsrtstructbuilder) <sub>class</sub> — Build an RT Structure Set from masks or contours.
- [`create_rtstruct_dataset`](#rscreate-rtstruct-dataset) <sub>func</sub> — Create an empty RT Structure Set bound to an image series.
- [`create_roi_into_rs_ds`](#rscreate-roi-into-rs-ds) <sub>func</sub> — Declare an ROI in all three sequences an RTSTRUCT needs it in.
- [`add_contour_sequence_from_mask3d`](#rsadd-contour-sequence-from-mask3d) <sub>func</sub> — Extract contours from a 3D mask and attach them to an existing ROI.
- [`add_contour_sequence_from_dcm_ctr_dict`](#rsadd-contour-sequence-from-dcm-ctr-dict) <sub>func</sub> — Attach contour points to an ROI that already exists in the dataset.
- [`get_contours`](#rsget-contours) <sub>func</sub> — Read every contour out of an RT Structure Set, grouped by ROI and by image.
- [`get_roi_names`](#rsget-roi-names) <sub>func</sub> — Map ROI numbers to ROI names.
- [`is_rtstruct_matching_series`](#rsis-rtstruct-matching-series) <sub>func</sub> — Check whether an RT Structure Set belongs to a given image series.
- [`check_rtstruct_iod`](#rscheck-rtstruct-iod) <sub>func</sub> — Validate an RT Structure Set against its CIOD.
- [`rtstruct_to_masks`](#rsrtstruct-to-masks) <sub>func</sub> — Rasterise an RT Structure Set into 3D binary masks.
- [`calc_image_series_affine_mapping`](#rscalc-image-series-affine-mapping) <sub>func</sub> — Derive the patient-to-pixel affine and volume shape from an image series.
- [`calc_rs_affine_mapping`](#rscalc-rs-affine-mapping) <sub>func</sub> — Infer an affine and volume shape from the contours themselves.

**[`reg`](#reg)** — DICOM registration objects: build, read, validate.

- [`SpatialRegistrationBuilder`](#regspatialregistrationbuilder) <sub>class</sub> — Build a Spatial Registration (rigid/affine) object.
- [`DeformableSpatialRegistrationBuilder`](#regdeformablespatialregistrationbuilder) <sub>class</sub> — Build a Deformable Spatial Registration object.
- [`get_spatial_registrations`](#regget-spatial-registrations) <sub>func</sub> — Read the matrices out of a Spatial Registration object.
- [`get_deformable_registrations`](#regget-deformable-registrations) <sub>func</sub> — Read the deformation fields out of a Deformable Spatial Registration object.
- [`check_spatial_reg_iod`](#regcheck-spatial-reg-iod) <sub>func</sub> — Validate a Spatial Registration object against its CIOD.
- [`check_deformable_reg_iod`](#regcheck-deformable-reg-iod) <sub>func</sub> — Validate a Deformable Spatial Registration object against its CIOD.
- [`affine_to_homogeneous_matrix`](#regaffine-to-homogeneous-matrix) <sub>func</sub> — Convert a SimpleITK affine transform to a 4x4 homogeneous matrix.
- [`sitk_displacement_field_to_deformable_registration_grid`](#regsitk-displacement-field-to-deformable-registration-grid) <sub>func</sub> — Convert a SimpleITK displacement field into a Deformable Registration Grid item.

**[`reg.method`](#regmethod)** — Registration algorithms.

- [`demons_registration`](#regmethoddemons-registration) <sub>func</sub> — Deformable registration via Fast Symmetric-Forces Demons (SimpleITK) with multi-resolution strategy.
- [`rigid_registration`](#regmethodrigid-registration) <sub>func</sub> — Rigidly align a moving image to a fixed image.
- [`bspline_registration`](#regmethodbspline-registration) <sub>func</sub> — Deformable registration with a B-spline transform.
- [`demons_with_soft_mask`](#regmethoddemons-with-soft-mask) <sub>func</sub> — Run Demons registration weighted towards a region of interest.

**[`reg.pipeline`](#regpipeline)** — Staged registration and transform composition.

- [`registration_pipeline`](#regpipelineregistration-pipeline) <sub>func</sub> — Registration pipeline integrating rigid and deformable registration.
- [`compose_transforms`](#regpipelinecompose-transforms) <sub>func</sub> — Compose a stage-ordered list of transforms into a single resampling transform.
- [`preprocess_image`](#regpipelinepreprocess-image) <sub>func</sub> — Preprocess the image by applying multiple preprocessing steps according to the config.
- [`infer_default_pixel_value`](#regpipelineinfer-default-pixel-value) <sub>func</sub> — Resolve the padding value to use when resampling outside an image's extent.
- [`window_clip`](#regpipelinewindow-clip) <sub>func</sub> — Apply window clipping (clamping) to the image.
- [`align_image_extents`](#regpipelinealign-image-extents) <sub>func</sub> — Align two images' physical extents by cropping them to their intersection.
- [`get_image_physical_extent`](#regpipelineget-image-physical-extent) <sub>func</sub> — Compute the physical extent (bounding box) of the image in physical space.
- [`get_intersection_extent`](#regpipelineget-intersection-extent) <sub>func</sub> — Compute the intersection of two images' physical extents.
- [`crop_image_to_extent`](#regpipelinecrop-image-to-extent) <sub>func</sub> — Crop the image to the given physical extent.
- [`get_initial_rigid_transform`](#regpipelineget-initial-rigid-transform) <sub>func</sub> — Get initial rigid transform using SimpleITK's CenteredTransformInitializer.
- [`get_image_center`](#regpipelineget-image-center) <sub>func</sub> — Compute the center point of the image in physical space.
- [`get_images_distance`](#regpipelineget-images-distance) <sub>func</sub> — Compute the Euclidean distance between the centers of two images.
- [`create_reference_image_from_extent`](#regpipelinecreate-reference-image-from-extent) <sub>func</sub> — Create a reference image from the given physical extent.

**[`dose`](#dose)** — RT Dose: build and read dose grids.

- [`RTDoseBuilder`](#dosertdosebuilder) <sub>class</sub> — Build an RT Dose object from a dose grid.
- [`generate_base_dataset`](#dosegenerate-base-dataset) <sub>func</sub> — Create an empty RT Dose dataset with its required elements.
- [`cp_information_from_ds`](#dosecp-information-from-ds) <sub>func</sub> — Copy patient, study and frame-of-reference context into a dose dataset.
- [`add_dose_grid_to_ds`](#doseadd-dose-grid-to-ds) <sub>func</sub> — Attach a dose grid to an RT Dose dataset.
- [`check_rtdose_iod`](#dosecheck-rtdose-iod) <sub>func</sub> — Validate an RT Dose dataset against its CIOD.
- [`get_dose_image`](#doseget-dose-image) <sub>func</sub> — Read an RT Dose object into a `sitk.Image`.
- [`get_dose_array`](#doseget-dose-array) <sub>func</sub> — Dose values as a numpy array, scaling applied.
- [`get_dose_spacing`](#doseget-dose-spacing) <sub>func</sub> — Voxel size of an RT Dose grid, in SimpleITK order.
- [`get_dose_origin`](#doseget-dose-origin) <sub>func</sub> — Patient coordinates of the dose grid's first voxel.
- [`get_dose_direction`](#doseget-dose-direction) <sub>func</sub> — Axis directions of an RT Dose grid, in SimpleITK convention.

**[`ct`](#ct)** — CT Image: emit a series from a volume.

- [`CTBuilder`](#ctctbuilder) <sub>class</sub> — Emit a CT Image series from a volume, borrowing patient and study context from a
- [`check_ct_iod`](#ctcheck-ct-iod) <sub>func</sub> — Validate a CT Image dataset against its CIOD.

---

<a id="utils"></a>

## utils

Series loading, geometry, and SimpleITK conversion.

<a id="utilsload-sorted-image-series"></a>

### `utils.load_sorted_image_series`

```python
load_sorted_image_series(image_series_path: Union[str, List[str]]) -> List[pydicom.dataset.Dataset]
```

Read a DICOM image series and return it sorted along the slice normal.

**Parameters**

- **`image_series_path`** — *str or list of str*
  A directory containing `*.dcm` files, or an explicit list of file paths. The
  directory form globs `*.dcm` only and does not recurse, so extensionless DICOM
  (common on PACS exports) must be passed as a list.

**Returns**

- **`list of Dataset`**
  Slices in ascending slice-normal order, fully read into memory.

**Raises**

- **`FileNotFoundError`**
  If the directory holds no `*.dcm` files, or a listed path does not exist. Also
  raised when a *file* path is passed where a directory was expected, since that
  matches neither branch and leaves the path list empty.

- **`ValueError`**
  If a file cannot be parsed as DICOM, or the series geometry is unusable. The
  message names the offending file.

**See Also**

- [`sort_image_series`](#utilssort-image-series) — Sorting alone, for slices you already have.
- [`image_series_to_sitk_image`](#utilsimage-series-to-sitk-image) — The usual next call.

**Examples**

```python
>>> series = load_sorted_image_series("/path/to/CT")
>>> len(series)
120
```

Extensionless files, or a subset:

```python
>>> from pathlib import Path
>>> paths = sorted(str(p) for p in Path(folder).iterdir())
>>> series = load_sorted_image_series(paths)
```

**Notes**

No filtering by SeriesInstanceUID is performed. A directory holding two series
produces one interleaved list that will sort into a physically meaningless order.

<a id="utilssort-image-series"></a>

### `utils.sort_image_series`

```python
sort_image_series(image_series: List[pydicom.dataset.Dataset]) -> List[pydicom.dataset.Dataset]
```

Sort slices by their position along the slice normal, ascending.

Use this when the slices came from somewhere other than
`load_sorted_image_series` -- a PACS query, a database, a list you filtered --
since every other function in the library assumes this ordering.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices of one series, in any order.

**Returns**

- **`list of Dataset`**
  A new list, ascending along the slice normal. The input is not modified and the
  `Dataset` objects are shared, not copied.

**Raises**

- **`InvalidImageOrientationError`**
  If any slice's ImageOrientationPatient is not two orthogonal unit vectors.

**Notes**

Sorting is by projection onto the slice normal, not by `InstanceNumber` and not by
`ImagePositionPatient[2]`. Scanners assign instance numbers in acquisition order,
which need not be spatial order, and the z component of the position is only the
slice-normal projection for axis-aligned axial series.

<a id="utilsinvalidimageorientationerror"></a>

### `utils.InvalidImageOrientationError`

```python
InvalidImageOrientationError
```

ImageOrientationPatient (0020,0037) is not two orthogonal unit vectors.

A subclass of `ValueError`, so `except ValueError` still catches it. Raised by
`get_slice_directions` and therefore by everything that reads series geometry:
loading, sorting, contour placement and SimpleITK conversion.

In practice this means the file is corrupt or non-conformant. The tolerance is 1e-3 on
both the dot product and the normal's length, which is loose enough for the rounding
real scanners write and tight enough to catch a genuinely bad orientation.

<a id="utilsget-slice-directions"></a>

### `utils.get_slice_directions`

```python
get_slice_directions(image_slice: pydicom.dataset.Dataset) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
```

Split ImageOrientationPatient into its three unit direction vectors.

**Parameters**

- **`image_slice`** — *Dataset*
  Any slice of the series; orientation is shared across a series.

**Returns**

- **`row_direction`** — *np.ndarray*
  Travelled when the *column* index increases -- IOP[0:3].

- **`column_direction`** — *np.ndarray*
  Travelled when the *row* index increases -- IOP[3:6].

- **`slice_direction`** — *np.ndarray*
  Their cross product, i.e. the slice normal.

**Raises**

- **`InvalidImageOrientationError`**
  If the two stored vectors are not orthogonal unit vectors.

<a id="utilsget-slice-position"></a>

### `utils.get_slice_position`

```python
get_slice_position(image_slice: pydicom.dataset.Dataset) -> float
```

Position of a slice along the slice normal.

This is the value `sort_image_series` orders by. It is not the same as
`ImagePositionPatient[2]` unless the series is axial and axis-aligned, and it is not
`InstanceNumber`, which scanners are free to assign in any order.

**Parameters**

- **`image_slice`** — *Dataset*
  Slice carrying ImagePositionPatient and ImageOrientationPatient.

**Returns**

- **`float`**
  Signed distance in mm from the patient origin, projected onto the slice normal.

**Raises**

- **`InvalidImageOrientationError`**
  If ImageOrientationPatient is not two orthogonal unit vectors.

<a id="utilsget-spacing-between-slices"></a>

### `utils.get_spacing_between_slices`

```python
get_spacing_between_slices(image_series: List[pydicom.dataset.Dataset]) -> float
```

Mean distance between slice centres, along the slice normal.

Measured from the first and last slice rather than read from SliceThickness or
SpacingBetweenSlices, both of which are routinely absent or wrong.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices **sorted along the slice normal**. On an unsorted series this returns a
  value that is too small, or negative.

**Returns**

- **`float`**
  Mean centre-to-centre distance in mm. Negative if the series runs against the
  slice normal. Exactly `1.0` for a single-slice series, which has no spacing to
  measure -- the fallback keeps the transformation matrices invertible.

**Notes**

Averaging over the whole series rather than differencing neighbours means a series
with one duplicated or missing slice still yields a usable spacing, at the cost of
hiding the irregularity. Uniform spacing is assumed throughout the library.

<a id="utilsget-pixel-to-patient-transformation-matrix"></a>

### `utils.get_pixel_to_patient_transformation_matrix`

```python
get_pixel_to_patient_transformation_matrix(image_series: List[pydicom.dataset.Dataset]) -> numpy.ndarray
```

4x4 affine mapping pixel indices to patient coordinates (mm).

Points are ordered `(column index, row index, slice index)`, which is the reverse of
DICOM's PixelSpacing ordering -- see the spacing pairing below.

**Parameters**

- **`image_series`** — *list[Dataset]*
  Slices sorted along the slice normal.

**Returns**

- **`np.ndarray`**
  4x4 matrix; multiply column vectors `[x_index, y_index, z_index, 1]`.

<a id="utilsget-patient-to-pixel-transformation-matrix"></a>

### `utils.get_patient_to_pixel_transformation_matrix`

```python
get_patient_to_pixel_transformation_matrix(image_series: List[pydicom.dataset.Dataset]) -> numpy.ndarray
```

4x4 affine mapping patient coordinates (mm) to pixel indices.

The inverse of `get_pixel_to_patient_transformation_matrix`, computed in closed
form rather than by inverting that matrix numerically.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices sorted along the slice normal.

**Returns**

- **`np.ndarray`**
  4x4. Multiply column vectors `[x_mm, y_mm, z_mm, 1]`; the result is ordered
  `(column index, row index, slice index)` and is generally not integral -- round
  or floor according to what you are doing.

**See Also**

- [`apply_transformation_to_3d_points`](#utilsapply-transformation-to-3d-points) — Apply this to an `(n, 3)` point array.
- [`calc_image_series_affine_mapping`](#rscalc-image-series-affine-mapping) — What RTSTRUCT reading needs -- this matrix plus the volume shape.

<a id="utilsapply-transformation-to-3d-points"></a>

### `utils.apply_transformation_to_3d_points`

```python
apply_transformation_to_3d_points(points: numpy.ndarray, transformation_matrix: numpy.ndarray) -> numpy.ndarray
```

Apply a 4x4 affine to an `(n, 3)` array of points.

Each point is augmented with a 1 so the translation column applies, then the extra
coordinate is dropped again.

**Parameters**

- **`points`** — *np.ndarray*
  `(n, 3)`. A single point must still be shaped `(1, 3)`.

- **`transformation_matrix`** — *np.ndarray*
  4x4 affine, such as either of the transformation matrices in this module.

**Returns**

- **`np.ndarray`**
  `(n, 3)`, transformed.

**Examples**

```python
>>> import numpy as np
>>> shift = np.identity(4); shift[:3, 3] = [10.0, 0.0, 0.0]
>>> apply_transformation_to_3d_points(np.array([[1.0, 2.0, 3.0]]), shift)
array([[11.,  2.,  3.]])
```

<a id="utilssimpleitkimagebuilder"></a>

### `utils.SimpleITKImageBuilder`

```python
SimpleITKImageBuilder()
```

Assemble a `sitk.Image` from a volume plus geometry.

Either set the four pieces individually and call `build`, or use one of the
`from_*` shortcuts, which build and return a `sitk.Image` in one step.

The setters return `self` so calls chain; the `from_*` shortcuts do not, because
they return the finished image.

**Methods**

- **`set_volume(volume)`**
  The voxels, `(slice, row, column)`.

- **`set_origin(origin)`**
  Patient coordinates of voxel `[0, 0, 0]`.

- **`set_spacing(spacing)`**
  Voxel size as `(x, y, z)`.

- **`set_direction(direction)`**
  3x3 with axis directions in *columns*.

- **`build()`**
  Assemble and return the image.

- **`from_image_series(image_series)`**
  Shortcut: build from DICOM slices.

- **`from_reference_image(volume, reference_image)`**
  Shortcut: build from a volume, borrowing geometry.

- **`from_dicom_directory(directory)`**
  Shortcut: build via SimpleITK's GDCM reader.

**See Also**

- [`image_series_to_sitk_image`](#utilsimage-series-to-sitk-image) — The common case, as a single function call.

**Examples**

Put a volume of your own onto an existing series' geometry:

```python
>>> volume, origin, spacing, direction = parse_image_series(series)
>>> image = (
...     SimpleITKImageBuilder()
...     .set_volume(my_segmentation)
...     .set_origin(origin)
...     .set_spacing(spacing)
...     .set_direction(direction)
...     .build()
... )
```

#### `SimpleITKImageBuilder.build`

```python
build(self) -> SimpleITK.SimpleITK.Image
```

Assemble the image.

**Returns**

- **`sitk.Image`**
  Volume with origin, spacing and direction applied.

**Raises**

- **`ValueError`**
  If volume, origin, spacing or direction has not been set. The message names
  the missing pieces.

#### `SimpleITKImageBuilder.from_dicom_directory`

```python
from_dicom_directory(self, directory: str) -> SimpleITK.SimpleITK.Image
```

Build via SimpleITK's own GDCM series reader.

Reads the directory with ITK's reader rather than this library's series loader.
Useful as an independent cross-check, and when the directory holds a series this
library's loader rejects.

**Parameters**

- **`directory`** — *str*
  Directory holding one DICOM series. If it holds several, GDCM picks one.

**Returns**

- **`sitk.Image`**
  The finished image; this shortcut does not return `self`.

**Notes**

Unlike `from_image_series` this path gives you no access to the underlying
`Dataset` objects, so ROI and registration builders -- which need patient and
frame-of-reference context -- cannot be driven from it.

#### `SimpleITKImageBuilder.from_image_series`

```python
from_image_series(self, image_series: List[pydicom.dataset.Dataset]) -> SimpleITK.SimpleITK.Image
```

Build directly from a DICOM image series.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices of one series, in any order.

**Returns**

- **`sitk.Image`**
  The finished image -- this shortcut calls `build` for you, so it does
  not return `self`.

#### `SimpleITKImageBuilder.from_reference_image`

```python
from_reference_image(self, volume: numpy.ndarray, reference_image: SimpleITK.SimpleITK.Image) -> SimpleITK.SimpleITK.Image
```

Build from a volume, borrowing geometry from an existing image.

The usual way to wrap a numpy result -- a segmentation, a dose grid, a filtered
volume -- so it sits on the same grid as the image it was computed from.

**Parameters**

- **`volume`** — *np.ndarray*
  3D, `(slice, row, column)`. Must match `reference_image` in shape, since
  it inherits that image's geometry without being resampled.

- **`reference_image`** — *sitk.Image*
  Supplies origin, spacing and direction.

**Returns**

- **`sitk.Image`**
  The finished image; this shortcut does not return `self`.

**See Also**

- [`resample_to_reference_image`](#utilsresample-to-reference-image) — When the volume is on a *different* grid and needs resampling rather than relabelling.

#### `SimpleITKImageBuilder.set_direction`

```python
set_direction(self, direction: Sequence[float]) -> 'SimpleITKImageBuilder'
```

Set the axis directions.

**Parameters**

- **`direction`** — *sequence of float*
  3x3 array or flat 9-sequence. The axis vectors go in the **columns**: column
  `j` is the unit vector travelled when index axis `j` increases. Building it
  with the vectors as rows gives the transpose, which is identical for
  axis-aligned axial series and mirrored for oblique ones. Use
  `get_sitk_direction`.

**Returns**

- **`SimpleITKImageBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `direction` is not 3x3 or a flat 9-sequence.

#### `SimpleITKImageBuilder.set_origin`

```python
set_origin(self, origin: Sequence[float]) -> 'SimpleITKImageBuilder'
```

Set the patient coordinates of voxel `[0, 0, 0]`.

**Parameters**

- **`origin`** — *sequence of float*
  Three values in mm. For a DICOM series this is ImagePositionPatient of the
  first slice *after sorting along the slice normal*.

**Returns**

- **`SimpleITKImageBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `origin` is not a 3-element vector.

#### `SimpleITKImageBuilder.set_spacing`

```python
set_spacing(self, spacing: Sequence[float]) -> 'SimpleITKImageBuilder'
```

Set the voxel size.

**Parameters**

- **`spacing`** — *sequence of float*
  `(x, y, z)` in mm -- SimpleITK's order. DICOM PixelSpacing is
  `[row, column]`, i.e. `(y, x)`, so the in-plane pair must be swapped, not
  copied across. Use `get_sitk_spacing` rather than doing it by hand.

**Returns**

- **`SimpleITKImageBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `spacing` is not a 3-element vector.

#### `SimpleITKImageBuilder.set_volume`

```python
set_volume(self, volume: numpy.ndarray) -> 'SimpleITKImageBuilder'
```

Set the voxel data.

**Parameters**

- **`volume`** — *np.ndarray*
  3D, ordered `(slice, row, column)` -- the order a series stacks to, and the
  reverse of what `sitk.Image.GetSize()` reports.

**Returns**

- **`SimpleITKImageBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `volume` is not 3D.

<a id="utilsimage-series-to-sitk-image"></a>

### `utils.image_series_to_sitk_image`

```python
image_series_to_sitk_image(image_series: List[pydicom.dataset.Dataset]) -> SimpleITK.SimpleITK.Image
```

Convert a DICOM image series straight to a `sitk.Image`.

The one-call form of `SimpleITKImageBuilder.from_image_series`, and the usual
entry point for reading images.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices of one series, in any order. They are sorted along the slice normal.

**Returns**

- **`sitk.Image`**
  The volume with origin, spacing and direction set, in the modality's real units
  (HU for CT -- rescale slope and intercept are applied).

**See Also**

- [`parse_image_series`](#utilsparse-image-series) — Same reading, returning the four pieces separately.
- [`SimpleITKImageBuilder`](#utilssimpleitkimagebuilder) — Assemble an image from a volume of your own.

**Examples**

```python
>>> series = load_sorted_image_series("/path/to/CT")
>>> image = image_series_to_sitk_image(series)
>>> image.GetSize()
(512, 512, 120)
```

<a id="utilsparse-image-series"></a>

### `utils.parse_image_series`

```python
parse_image_series(image_series: List[pydicom.dataset.Dataset]) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
```

Read volume and geometry out of a DICOM image series.

Rescale slope and intercept are applied, so the volume is in the modality's real units
(HU for CT). The series is sorted along the slice normal first.

**Parameters**

- **`image_series`** — *list[Dataset]*
  Slices of one series.

**Returns**

- **`volume`** — *np.ndarray*
  `(slice, row, column)`, rescaled.

- **`origin`** — *np.ndarray*
  ImagePositionPatient of the first slice.

- **`spacing`** — *np.ndarray*
  `(x, y, z)`.

- **`direction`** — *np.ndarray*
  3x3, axis directions in columns.

**Notes**

Assumes a single, uniformly spaced series; spacing and orientation are taken from the
first slice.

<a id="utilsresample-to-reference-image"></a>

### `utils.resample_to_reference_image`

```python
resample_to_reference_image(reference_image: SimpleITK.SimpleITK.Image, source_image: SimpleITK.SimpleITK.Image, default_value: float = -1000.0, interpolator: int = 2) -> SimpleITK.SimpleITK.Image
```

Resample the source image onto the reference image grid (identity transform).

**Parameters**

- **`reference_image`** — *sitk.Image*
  Defines the output grid and pixel type.

- **`source_image`** — *sitk.Image*
  Image to resample.

- **`default_value`** — *float*
  Padding for samples outside the source extent. -1000 is CT air; pass 0 for MR/PET
  or already-clipped images, otherwise the padding invents a soft-tissue shell.

- **`interpolator`** — *int*
  SimpleITK interpolator. Use `sitkNearestNeighbor` for masks and labels.

**Returns**

- **`sitk.Image`**
  Source image on the reference grid.

<a id="utilssitk-direction-to-image-orientation-patient"></a>

### `utils.sitk_direction_to_image_orientation_patient`

```python
sitk_direction_to_image_orientation_patient(direction: Sequence[float]) -> List[float]
```

Convert a SimpleITK direction matrix into DICOM ImageOrientationPatient (0020,0037).

ImageOrientationPatient has VM 6: the row direction followed by the column direction.
Those are the first two *columns* of a SimpleITK direction matrix, which is not the same
as the first six entries of its row-major flattening -- the two only coincide when the
matrix is symmetric, as it is for axis-aligned axial series.

**Parameters**

- **`direction`** — *sequence*
  SimpleITK direction as a flat 9-sequence or a 3x3 array.

**Returns**

- **`list[float]`**
  Six floats: row direction then column direction.

<a id="utilssitk-spacing-to-pixel-spacing"></a>

### `utils.sitk_spacing_to_pixel_spacing`

```python
sitk_spacing_to_pixel_spacing(spacing: Sequence[float]) -> List[float]
```

Convert SimpleITK `(x, y, z)` spacing into DICOM PixelSpacing (0028,0030).

**Parameters**

- **`spacing`** — *sequence of float*
  SimpleITK spacing, `(x, y, z)`. Only the first two entries are used.

**Returns**

- **`list of float`**
  `[row spacing, column spacing]`, i.e. `(y, x)` -- the DICOM order.

**See Also**

- `get_sitk_spacing` — The inverse direction, DICOM to SimpleITK.

**Examples**

```python
>>> sitk_spacing_to_pixel_spacing((0.5, 1.5, 3.0))
[1.5, 0.5]
```

<a id="utilscheck-iod"></a>

### `utils.check_iod`

```python
check_iod(ds, config_map: Dict[str, dict], validators: Dict[str, Callable] = None, path: str = '') -> List[str]
```

Check a DICOM dataset against an IOD description.

The engine behind every `check_*_iod` function. Call those instead unless you are
describing a new object type; this one takes the IOD as data so a new modality needs
an `iod.py` rather than new validation code.

**Parameters**

- **`ds`** — *Dataset*
  The dataset to check. Nested sequences are recursed into via `submap`.

- **`config_map`** — *dict*
  Keyword to requirement, nestable:

  ```python
  {
      "PatientID":     {},
      "ROIContourSequence": {
          "min_items": 1,
          "submap": {
              "ReferencedROINumber": {"nonempty": True},
          },
      },
      "DoseSummationType": {"nonempty": True, "validator": ["dose_summation"]},
  }
  ```

  Elements are required unless `optional=True`. `nonempty=True` rejects empty
  values; `min_items` and `max_items` constrain sequence cardinality. `type`
  is a Python value class, not the DICOM requirement number. `validator` names
  entries in `validators`. `submap` describes sequence items.

- **`validators`** — *dict of str to callable, optional*
  Named checks beyond presence, as `{name: func}`. Each is called with the element
  value and raises `ValidationError` on failure. Used
  for the conditional (Type 1C) requirements presence alone cannot express.

- **`path`** — *str, optional*
  Prefix used when reporting nested elements; set by the recursion. Leave it alone
  at the top level, where it produces messages ending "in root".

**Returns**

- **`list of str`**
  One message per problem, empty when the dataset conforms. Pass it through
  `build_check_result` to get the `{"result", "content"}` shape the public
  checkers return.

**See Also**

- `build_check_result` — Wrap the returned list in the public result shape.

**Notes**

This checks that required elements are *present*, not that their values are correct.
A conformant object can still describe the wrong geometry.

---

<a id="rs"></a>

## rs

RT Structure Sets: masks to contours and back.

<a id="rsrtstructbuilder"></a>

### `rs.RTStructBuilder`

```python
RTStructBuilder(image_series: List[pydicom.dataset.Dataset])
```

Build an RT Structure Set from masks or contours.

Patient, study and frame-of-reference context is taken from the image series the
structures were drawn on, so the output files into the same study as those images.

**Parameters**

- **`image_series`** — *list of Dataset*
  The slices the structures were drawn on, sorted along the slice normal. Every mask
  added must match this series in shape, and the contours are placed using its
  geometry.

**Methods**

- **`add_roi(mask, name, ...)`**
  Add an ROI from a 3D binary mask.

- **`add_roi_from_contours(contours, name, ...)`**
  Add an ROI from points already in patient coordinates.

- **`set_structure_set_label(label)`**
  Override the generated StructureSetLabel.

- **`set_uid_prefix(uid_prefix)`**
  Override the UID root.

- **`build()`**
  Assemble the `FileDataset`.

**Raises**

- **`ValueError`**
  If `image_series` is empty.

**See Also**

- [`rtstruct_to_masks`](#rsrtstruct-to-masks) — The reverse direction, RTSTRUCT to voxels.
- [`check_rtstruct_iod`](#rscheck-rtstruct-iod) — Validate the result before writing it.

**Notes**

The UID root defaults to the `DICOM_UID_PREFIX` environment variable. Set it, or
call `set_uid_prefix`, before writing files that leave your machine.

**Examples**

```python
>>> rs_ds = (
...     RTStructBuilder(image_series)
...     .add_roi(mask=ctv_mask, name="CTV", color=[0, 255, 0], interpreted_type="CTV")
...     .add_roi(mask=cord_mask, name="SpinalCord", interpreted_type="ORGAN")
...     .build()
... )
>>> rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)
```

#### `RTStructBuilder.add_roi`

```python
add_roi(self, mask: numpy.ndarray, name: str, number: Optional[int] = None, color: Sequence[int] = (255, 0, 0), description: str = '', interpreted_type: str = 'ORGAN', contour_config: Optional[Dict] = None) -> 'RTStructBuilder'
```

Add an ROI from a 3D mask.

**Parameters**

- **`mask`** — *np.ndarray*
  `(slice, row, column)` -- the same axis order as the image series stacks to.
  Non-zero voxels are inside the ROI.

- **`name`** — *str*
  ROIName, as it appears in the TPS.

- **`number`** — *int, optional*
  ROINumber. Assigned in sequence when omitted.

- **`color`** — *sequence[int]*
  ROIDisplayColor, RGB.

- **`description`** — *str*
  ROIDescription.

- **`interpreted_type`** — *str*
  RTROIInterpretedType, e.g. `"ORGAN"`, `"PTV"`, `"EXTERNAL"`.

- **`contour_config`** — *dict, optional*
  Contour extraction tuning; see `DEFAULT_CONTOUR_CONFIG`.

**Returns**

- **`RTStructBuilder`**
  self, so calls chain.

#### `RTStructBuilder.add_roi_from_contours`

```python
add_roi_from_contours(self, contours: Dict, name: str, number: Optional[int] = None, color: Sequence[int] = (255, 0, 0), description: str = '', interpreted_type: str = 'ORGAN') -> 'RTStructBuilder'
```

Add an ROI from contour points already in patient coordinates.

**Parameters**

- **`contours`** — *dict*
  Keyed by the SOP Instance UID of the image each contour lies on:

  ```python
  {sop_instance_uid: {"sop_class_uid": str,
                      "contours": [[x1, y1, z1, x2, y2, z2, ...], ...]}}
  ```

  Points are patient coordinates in mm, interleaved into one flat list per
  contour -- the DICOM ContourData layout, and what `get_contours` returns
  under `"dcm_contour"`. The UIDs must belong to the builder's image series.

- **`name`** — *str*
  ROIName, as it appears in the TPS.

- **`number`** — *int, optional*
  ROINumber. Assigned in sequence when omitted.

- **`color`** — *sequence of int*
  ROIDisplayColor, RGB 0-255. Default `(255, 0, 0)`.

- **`description`** — *str*
  ROIDescription. Default empty.

- **`interpreted_type`** — *str*
  RTROIInterpretedType, e.g. `"ORGAN"`, `"PTV"`, `"EXTERNAL"`. Default
  `"ORGAN"`.

**Returns**

- **`RTStructBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `name` is empty, or `number` was already added.

**See Also**

- `add_roi` — When you have a mask rather than points.
- [`get_contours`](#rsget-contours) — Produces this structure from an existing RTSTRUCT.

**Notes**

The points are written as given -- no contour filtering, no resampling, no
planarity check. This is the path to use when contours must round-trip unchanged.

#### `RTStructBuilder.build`

```python
build(self) -> pydicom.dataset.FileDataset
```

Assemble the RT Structure Set.

Contours are extracted here, not when `add_roi` is called, so the cost of
`cv2.findContours` over every mask lands on this call.

**Returns**

- **`FileDataset`**
  Ready to `save_as(path, enforce_file_format=True)`.

**Raises**

- **`ValueError`**
  If no ROI has been added.

**Notes**

The builder is not consumed. Calling `build()` twice produces two datasets with
*different* SOP Instance UIDs and the same content, which is usually not what you
want -- keep the first result rather than rebuilding.

#### `RTStructBuilder.set_structure_set_label`

```python
set_structure_set_label(self, label: str) -> 'RTStructBuilder'
```

Override the generated StructureSetLabel.

**Parameters**

- **`label`** — *str*
  VR SH, so 16 characters. Longer labels are what TPSs display, and a
  non-conformant length is a common reason for an import to be refused.

**Returns**

- **`RTStructBuilder`**
  self, so calls chain.

#### `RTStructBuilder.set_uid_prefix`

```python
set_uid_prefix(self, uid_prefix: str) -> 'RTStructBuilder'
```

Set the root under which generated UIDs are minted.

**Parameters**

- **`uid_prefix`** — *str*
  Organisation UID root, dot-terminated, e.g. `"1.2.826.0.1.3680043.2.1125."`.
  Defaults to the `DICOM_UID_PREFIX` environment variable.

**Returns**

- **`RTStructBuilder`**
  self, so calls chain.

**Notes**

UIDs generated under a root you do not own can collide with another
organisation's. That matters for anything leaving your machine and not at all for
scratch files.

<a id="rscreate-rtstruct-dataset"></a>

### `rs.create_rtstruct_dataset`

```python
create_rtstruct_dataset(image_series, uid_prefix: str = None) -> pydicom.dataset.FileDataset
```

Create an empty RT Structure Set bound to an image series.

Everything except the ROIs: file meta, patient and study context copied from the
images, the referenced frame-of-reference sequence listing every slice, and the three
empty ROI sequences ready to be filled.

Called by `RTStructBuilder`. Use it directly when you need to interleave your
own dataset edits with ROI creation.

**Parameters**

- **`image_series`** — *list of Dataset*
  Slices sorted along the slice normal. Patient identity, study, frame of reference
  and the per-slice references all come from here.

- **`uid_prefix`** — *str, optional*
  UID root. Defaults to the `DICOM_UID_PREFIX` environment variable.

**Returns**

- **`FileDataset`**
  A conformant but ROI-less RT Structure Set. It will not pass
  `check_rtstruct_iod` until at least one ROI is added.

**See Also**

- [`RTStructBuilder`](#rsrtstructbuilder) — The supported entry point.
- [`create_roi_into_rs_ds`](#rscreate-roi-into-rs-ds) — Add an ROI to the result.

<a id="rscreate-roi-into-rs-ds"></a>

### `rs.create_roi_into_rs_ds`

```python
create_roi_into_rs_ds(rs_ds: pydicom.dataset.Dataset, roi_color: list, roi_number: int, roi_name: str, roi_description: str, roi_interpreted_type: str = 'ORGAN') -> pydicom.dataset.Dataset
```

Declare an ROI in all three sequences an RTSTRUCT needs it in.

An ROI is not one element but three parallel entries, cross-referenced by ROINumber:
StructureSetROISequence (name and frame of reference), ROIContourSequence (colour and,
later, the contours), and RTROIObservationsSequence (interpreted type). Adding it to
only some of them produces a file most TPSs silently ignore.

This creates all three, with an empty ContourSequence for the contour extraction to
fill. Called by `RTStructBuilder.add_roi`, which is the supported entry point.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set with the three sequences present, from
  `create_rtstruct_dataset`. Modified in place.

- **`roi_color`** — *list of int*
  ROIDisplayColor as RGB, 0-255.

- **`roi_number`** — *int*
  ROINumber. Must be unique within the dataset; nothing here enforces that.

- **`roi_name`** — *str*
  ROIName, as the TPS displays it.

- **`roi_description`** — *str*
  ROIDescription. May be empty.

- **`roi_interpreted_type`** — *str, optional*
  RTROIInterpretedType, e.g. `"ORGAN"`, `"PTV"`, `"CTV"`, `"EXTERNAL"`.
  Default `"ORGAN"`.

**Returns**

- **`Dataset`**
  `rs_ds`, modified in place.

**See Also**

- [`RTStructBuilder.add_roi`](#rsrtstructbuilder) — Adds the ROI *and* its contours, and rejects duplicate numbers.
- [`add_contour_sequence_from_mask3d`](#rsadd-contour-sequence-from-mask3d) — Fills in the contours afterwards.

<a id="rsadd-contour-sequence-from-mask3d"></a>

### `rs.add_contour_sequence_from_mask3d`

```python
add_contour_sequence_from_mask3d(rs_ds: pydicom.dataset.Dataset, image_ds_list: list[pydicom.dataset.Dataset], roi_number: int, mask_volume: numpy.ndarray, ctr_config: dict = {'ex_noise_size': 10, 'in_noise_size': 10, 'lowpass_ratio': 10, 'ctr_precision': 8}) -> pydicom.sequence.Sequence
```

Extract contours from a 3D mask and attach them to an existing ROI.

The building block behind `RTStructBuilder.add_roi`, and where the mask actually
becomes contours: `cv2.findContours` per slice, noise removal and low-pass filtering
per contour, then pixel-to-patient conversion.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set with an ROIContourSequence containing `roi_number`. Modified
  in place.

- **`image_ds_list`** — *list of Dataset*
  The referenced image series, sorted along the slice normal. Supplies the geometry
  that turns pixel indices into patient coordinates.

- **`roi_number`** — *int*
  ROINumber to attach to. Must already exist.

- **`mask_volume`** — *np.ndarray*
  `(slice, row, column)`, matching the series. Non-zero is inside the ROI.

- **`ctr_config`** — *dict, optional*
  Contour extraction tuning:

  `ex_noise_size` : int
  Drop exterior contours enclosing fewer than this many pixels.
  `in_noise_size` : int
  Same for interior contours -- holes.
  `lowpass_ratio` : int
  Fourier low-pass strength; higher keeps more detail and more staircasing.
  `ctr_precision` : int
  Decimal places kept in the written coordinates.

  Defaults to `DEFAULT_CONTOUR_CONFIG`. Note this is a mutable default shared
  between calls -- pass a fresh dict rather than mutating it.

**Returns**

- **`Dataset`**
  `rs_ds`, modified in place. Despite the annotation this is the dataset, not a
  `Sequence`.

**Raises**

- **`ValueError`**
  If `rs_ds` has no ROIContourSequence, or no item matches `roi_number`.

**See Also**

- [`RTStructBuilder.add_roi`](#rsrtstructbuilder) — The supported entry point, which also checks the mask shape.
- [`rtstruct_to_masks`](#rsrtstruct-to-masks) — The reverse direction.

**Notes**

Filtering is lossy in both directions: it removes single-voxel speckle that would
otherwise become degenerate contours, and it rounds off genuine fine detail. For a
mask that must round-trip exactly, raise `lowpass_ratio` and set both noise sizes
to 0.

<a id="rsadd-contour-sequence-from-dcm-ctr-dict"></a>

### `rs.add_contour_sequence_from_dcm_ctr_dict`

```python
add_contour_sequence_from_dcm_ctr_dict(rs_ds: pydicom.dataset.Dataset, image_ds_list: list[pydicom.dataset.Dataset], roi_number: int, dcm_ctr_dict: dict) -> pydicom.sequence.Sequence
```

Attach contour points to an ROI that already exists in the dataset.

A building block of `RTStructBuilder.add_roi_from_contours`; call the builder
unless you are assembling a dataset element by element.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set with an ROIContourSequence containing `roi_number`. Modified
  in place.

- **`image_ds_list`** — *list of Dataset*
  The referenced image series, used to match points to slices.

- **`roi_number`** — *int*
  ROINumber to attach to. Must already exist -- create it with
  `create_roi_into_rs_ds` first.

- **`dcm_ctr_dict`** — *dict*
  Points in patient coordinates, keyed by referenced SOP Instance UID:

  ```python
  {sop_instance_uid: {"sop_class_uid": str,
                      "contours": [[x1, y1, z1, x2, y2, z2, ...], ...]}}
  ```

  The layout `get_contours` returns under `"dcm_contour"`.

**Returns**

- **`Dataset`**
  `rs_ds`, modified in place. Despite the annotation this is the dataset, not a
  `Sequence`.

**Raises**

- **`ValueError`**
  If `rs_ds` has no ROIContourSequence, or no item matches `roi_number`.

**See Also**

- [`RTStructBuilder.add_roi_from_contours`](#rsrtstructbuilder) — The supported entry point.
- [`add_contour_sequence_from_mask3d`](#rsadd-contour-sequence-from-mask3d) — Same, starting from a mask.

<a id="rsget-contours"></a>

### `rs.get_contours`

```python
get_contours(rs_ds: pydicom.dataset.Dataset) -> dict
```

Read every contour out of an RT Structure Set, grouped by ROI and by image.

Returns the points as stored -- patient coordinates in mm, no rasterisation and no
image series required. Use this to inspect, edit or copy contours; use
`rtstruct_to_masks` when you want voxels.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set. Missing sequences are tolerated and yield an empty result rather
  than raising.

**Returns**

- **`dict`**
  Keyed by ROI *number*:

  ```python
  {roi_number: {
      "color": [R, G, B],
      "name": str,
      "dcm_contour": {
          sop_instance_uid: {
              "sop_class_uid": str,
              "contours": [[x1, y1, z1, x2, y2, z2, ...], ...],
          },
      },
  }}
  ```

  Each contour is one flat list of interleaved xyz triples, the DICOM ContourData
  layout. Reshape with `np.reshape(contour, (-1, 3))`.

  Keys are `pydicom.valuerep.IS`, an `int` subclass: `contours[1]` works,
  `contours["1"]` raises `KeyError` despite the repr showing `{'1': ...}`.

**See Also**

- [`rtstruct_to_masks`](#rsrtstruct-to-masks) — Rasterise onto an image grid instead.
- [`get_roi_names`](#rsget-roi-names) — Just the number-to-name mapping.
- [`RTStructBuilder.add_roi_from_contours`](#rsrtstructbuilder) — Feed this structure back into a new RTSTRUCT.

**Notes**

Contours are grouped under the SOP Instance UID of the image they reference. A contour
whose ContourImageSequence is absent or empty is skipped, since it cannot be attributed
to an image.

**Examples**

```python
>>> contours = get_contours(rs_ds)
>>> {n: d["name"] for n, d in contours.items()}
{1: 'CTV', 2: 'SpinalCord'}
```

<a id="rsget-roi-names"></a>

### `rs.get_roi_names`

```python
get_roi_names(rs_ds: pydicom.dataset.Dataset) -> dict
```

Map ROI numbers to ROI names.

The cheapest way to see what is in a structure set: it reads StructureSetROISequence
only and never touches contour data.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set.

**Returns**

- **`dict`**
  `{roi_number: roi_name}`. An ROI with no ROIName is reported as
  `"ROI_<number>"` rather than `None`, so the value is always usable as a label.

  Keys are `pydicom.valuerep.IS`, an `int` subclass. `names[1]` works;
  `names["1"]` raises `KeyError` even though the repr prints `{'1': 'CTV'}`.

**See Also**

- [`get_contours`](#rsget-contours) — The same ROIs with their contour points.

**Notes**

DICOM does not require ROI names to be unique. Compare `len(names)` against
`len(set(names.values()))` before keying anything by name --
`rtstruct_to_masks` does key by name, and collapses duplicates.

**Examples**

```python
>>> get_roi_names(rs_ds)
{1: 'CTV', 2: 'SpinalCord', 3: 'BODY'}
```

<a id="rsis-rtstruct-matching-series"></a>

### `rs.is_rtstruct_matching_series`

```python
is_rtstruct_matching_series(rs_ds: pydicom.dataset.Dataset, image_series: list) -> bool
```

Check whether an RT Structure Set belongs to a given image series.

Worth calling before `rtstruct_to_masks`, since rasterising a structure set
against the wrong series produces masks in the wrong place rather than an error.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set.

- **`image_series`** — *list of Dataset*
  Candidate image series.

**Returns**

- **`bool`**
  True if the two are linked and no contour references an image outside the series.

**Notes**

The test is deliberately lenient in one direction and strict in the other. Evidence
*for* a match is either a shared FrameOfReferenceUID or a matching SeriesInstanceUID
in ReferencedFrameOfReferenceSequence -- one is enough, since real exports frequently
carry only one. Evidence *against* is decisive: a contour referencing a SOP Instance
UID not present in the series returns False regardless.

A structure set covering a *subset* of the series still matches; only references
pointing outside it fail. An empty `image_series` is False.

<a id="rscheck-rtstruct-iod"></a>

### `rs.check_rtstruct_iod`

```python
check_rtstruct_iod(rs_ds: pydicom.dataset.Dataset) -> dict
```

Validate an RT Structure Set against its CIOD.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set, as built or as read from disk.

**Returns**

- **`dict`**
  `{"result": bool, "content": list of str}`. `result` is True when the object
  is conformant and `content` is empty; otherwise `content` holds one message
  per missing or empty required element, each naming its path in the dataset.

**See Also**

- [`is_rtstruct_matching_series`](#rsis-rtstruct-matching-series) — Whether the structure set belongs to a given series.

**Notes**

Structural conformance only: this reports elements that are absent or empty, not
contours in the wrong place. An RTSTRUCT can pass here and still be geometrically
wrong.

**Examples**

```python
>>> result = check_rtstruct_iod(rs_ds)
>>> if not result["result"]:
...     print("\n".join(result["content"]))
```

<a id="rsrtstruct-to-masks"></a>

### `rs.rtstruct_to_masks`

```python
rtstruct_to_masks(rs_ds, affine_mapping, mask_volume_shape, roi_list=['all'], packbits=False)
```

Rasterise an RT Structure Set into 3D binary masks.

Each contour is transformed into pixel indices and filled with `cv2.fillPoly`,
slice by slice. Overlapping contours on the same slice cancel, which is how nested
contours become holes -- a ring drawn inside another ring on the same slice reads as
exterior, matching the even-odd rule TPSs use.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set.

- **`affine_mapping`** — *np.ndarray*
  4x4 patient-to-pixel affine, from `calc_image_series_affine_mapping`.

- **`mask_volume_shape`** — *tuple of int*
  `(slice, row, column)`, from the same call.

- **`roi_list`** — *list of str, optional*
  ROI names to rasterise. The default `["all"]` does every ROI. Filtering here is
  much cheaper than rasterising everything and discarding.

- **`packbits`** — *bool, optional*
  Pack each mask along the last axis with `np.packbits`, 8x smaller in memory.
  Unpack with `np.unpackbits(mask, axis=-1)` before use. Default False.

**Returns**

- **`dict`**
  Keyed by ROI *name*, not number:

  ```python
  {"CTV": {"mask_volume": np.ndarray,      # (slice, row, column), 0 or 1
           "affine_mapping": np.ndarray}}  # the affine passed in
  ```

  An ROI declared without contours yields an empty inner dict, so read
  `masks[name]["mask_volume"]` defensively when the source is unknown.

**See Also**

- [`calc_image_series_affine_mapping`](#rscalc-image-series-affine-mapping) — Produces both geometry arguments.
- [`RTStructBuilder.add_roi`](#rsrtstructbuilder) — The reverse direction, mask to contours.

**Notes**

ROI names are not required to be unique in DICOM. Two ROIs sharing a name collapse to
one entry, the later overwriting the earlier. Use `get_roi_names` to detect this
before trusting the mapping.

Each contour is assigned to the slice of its **first** point. Contours that are not
planar, or that sit exactly between two slices, land on one slice rather than being
split.

**Examples**

```python
>>> affine, shape = calc_image_series_affine_mapping(series)
>>> masks = rtstruct_to_masks(rs_ds, affine, shape, roi_list=["CTV"])
>>> masks["CTV"]["mask_volume"].sum()
18422
```

<a id="rscalc-image-series-affine-mapping"></a>

### `rs.calc_image_series_affine_mapping`

```python
calc_image_series_affine_mapping(series_ds_list)
```

Derive the patient-to-pixel affine and volume shape from an image series.

The companion call to `rtstruct_to_masks`, which needs to know the grid the
contours should be rasterised onto. That grid comes from the *images*, not from the
RTSTRUCT: an RTSTRUCT stores points in patient coordinates and carries no grid of its
own.

**Parameters**

- **`series_ds_list`** — *list of Dataset*
  Slices sorted along the slice normal. Sorting matters -- the origin is taken from
  `series_ds_list[0]`, so an unsorted series produces masks flipped in z.

**Returns**

- **`affine_mapping`** — *np.ndarray*
  4x4 mapping patient mm to pixel indices.

- **`mask_volume_shape`** — *tuple of int*
  `(slice, row, column)` -- number of slices, then Rows, then Columns.

**Raises**

- **`Exception`**
  If ImageOrientationPatient is not two orthogonal unit vectors.

**See Also**

- [`rtstruct_to_masks`](#rsrtstruct-to-masks) — Consumes both return values.
- [`calc_rs_affine_mapping`](#rscalc-rs-affine-mapping) — The fallback when the image series is unavailable.

**Examples**

```python
>>> affine, shape = calc_image_series_affine_mapping(series)
>>> masks = rtstruct_to_masks(rs_ds, affine, shape)
```

<a id="rscalc-rs-affine-mapping"></a>

### `rs.calc_rs_affine_mapping`

```python
calc_rs_affine_mapping(rs_ds)
```

Infer an affine and volume shape from the contours themselves.

A fallback for when the referenced image series is not available. The grid is
reconstructed from the contour points: the in-plane axes come from an oriented
bounding box, the slice axis from the contour normals, and the extent from the point
cloud's bounds.

Prefer `calc_image_series_affine_mapping` whenever you have the images. The
inferred grid is *not* the acquisition grid -- it is bounded by the contours, so it
is generally smaller, differently placed, and differently sampled. Masks produced
through it will not line up voxel-for-voxel with the CT.

**Parameters**

- **`rs_ds`** — *Dataset*
  RT Structure Set with a populated ROIContourSequence.

**Returns**

- **`affine_mapping`** — *np.ndarray*
  4x4 mapping patient mm to indices of the inferred grid.

- **`mask_volume_shape`** — *tuple of int*
  `(slice, row, column)` of the inferred grid.

**See Also**

- [`calc_image_series_affine_mapping`](#rscalc-image-series-affine-mapping) — The accurate path, when the images are on hand.

**Notes**

Contours with fewer than three points contribute no normal and are skipped when
estimating orientation.

---

<a id="reg"></a>

## reg

DICOM registration objects: build, read, validate.

<a id="regspatialregistrationbuilder"></a>

### `reg.SpatialRegistrationBuilder`

```python
SpatialRegistrationBuilder(fixed_series)
```

Build a Spatial Registration (rigid/affine) object.

Serialises a 4x4 matrix as DICOM Spatial Registration (SOP class
`1.2.840.10008.5.1.4.1.1.66.1`). Several moving series can be registered to one
fixed series in a single object.

**Parameters**

- **`fixed_series`** — *list of Dataset*
  The series everything is registered *to*. Patient and study context is copied
  from it, and its Frame of Reference becomes the registered RCS -- the frame every
  stored matrix is expressed relative to.

**Methods**

- **`add_registration(moving_series, matrix, matrix_type)`**
  Register one moving series. Call repeatedly for several.

- **`set_instance_number(instance_number)`**
  Set InstanceNumber.

- **`set_uid_prefix(uid_prefix)`**
  Set the UID root.

- **`build(include_identity=True)`**
  Assemble the `FileDataset`.

**Raises**

- **`ValueError`**
  If `fixed_series` is empty.

**See Also**

- [`get_spatial_registrations`](#regget-spatial-registrations) — Read the result back.
- [`DeformableSpatialRegistrationBuilder`](#regdeformablespatialregistrationbuilder) — When the transform is a displacement field.
- [`check_spatial_reg_iod`](#regcheck-spatial-reg-iod) — Validate before writing.

**Notes**

**The matrix runs moving to fixed**, which is the opposite of what the registration
functions in `pydicomrt.reg.method` return. Those return *resampling* transforms
(fixed to moving, the direction `sitk.Resample` wants), so invert before passing one
here. Getting this backwards produces a valid file that misplaces the image by twice
the offset.

**Examples**

```python
>>> transform = rigid_registration(fixed_image, moving_image)
>>> matrix = affine_to_homogeneous_matrix(transform.GetInverse())
>>> reg_ds = (
...     SpatialRegistrationBuilder(fixed_series)
...     .add_registration(moving_series, matrix.astype("float32").ravel().tolist())
...     .build()
... )
```

#### `SpatialRegistrationBuilder.add_registration`

```python
add_registration(self, moving_series, rigid_transform_matrix, matrix_type: str = 'RIGID') -> 'SpatialRegistrationBuilder'
```

Register a moving series to the fixed series' Frame of Reference.

**Parameters**

- **`moving_series`** — *list[Dataset]*
  The series being registered.

- **`rigid_transform_matrix`** — *sequence*
  16 numbers, row-major 4x4, mapping **moving -> fixed**. A SimpleITK resampling
  transform runs fixed -> moving, so invert it first.

- **`matrix_type`** — *str*
  `"RIGID"`, `"RIGID_SCALE"` or `"AFFINE"`. RIGID is validated as a rotation.

**Returns**

- **`SpatialRegistrationBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If the moving series is empty or the matrix does not validate.

#### `SpatialRegistrationBuilder.build`

```python
build(self, include_identity: bool = True) -> pydicom.dataset.FileDataset
```

Assemble the dataset.

**Parameters**

- **`include_identity`** — *bool*
  Append the identity item for the fixed Frame of Reference. Leave it on unless a
  downstream system objects to it.

**Raises**

- **`ValueError`**
  If no registration has been added or the assembled dataset fails validation.

#### `SpatialRegistrationBuilder.set_instance_number`

```python
set_instance_number(self: ~BuilderT, instance_number: int) -> ~BuilderT
```

Set InstanceNumber (0020,0013).

**Parameters**

- **`instance_number`** — *int*
  Position within the series. Default 1. Give each object a distinct number when
  writing several registrations into one series -- InstanceNumber is Type 1 here,
  and some archives reject duplicates within a series.

**Returns**

- **`SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder`**
  self, so calls chain.

#### `SpatialRegistrationBuilder.set_uid_prefix`

```python
set_uid_prefix(self: ~BuilderT, uid_prefix: str) -> ~BuilderT
```

Set the root under which generated UIDs are minted.

**Parameters**

- **`uid_prefix`** — *str*
  Organisation UID root, dot-terminated, e.g. `"1.2.826.0.1.3680043.2.1125."`.
  Unlike the RTSTRUCT and dose builders, the registration builders do **not**
  read `DICOM_UID_PREFIX`; unset, they fall back to pydicom's default root.

**Returns**

- **`SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder`**
  self, so calls chain.

<a id="regdeformablespatialregistrationbuilder"></a>

### `reg.DeformableSpatialRegistrationBuilder`

```python
DeformableSpatialRegistrationBuilder(fixed_series)
```

Build a Deformable Spatial Registration object.

Serialises a displacement field, optionally bracketed by two affine matrices, as
DICOM Deformable Spatial Registration (SOP class `1.2.840.10008.5.1.4.1.1.66.3`).

The three parts compose in a fixed order: **pre-matrix, then displacement grid, then
post-matrix**. That is what makes the usual rigid-then-deformable pipeline
expressible: for a resampling pipeline `rigid(deform(point))`, use identity
for `pre_transform`, the residual field for the grid, and the forward rigid
resampling matrix for `post_transform`.

**Parameters**

- **`fixed_series`** — *list of Dataset*
  The series everything is registered to. Supplies patient and study context, and
  its Frame of Reference becomes the registered RCS.

**Methods**

- **`add_registration(moving_series, vectorial_field_transform, pre_transform, post_transform)`**
  Register one moving series.

- **`set_instance_number(instance_number)`**
  Set InstanceNumber.

- **`set_uid_prefix(uid_prefix)`**
  Set the UID root.

- **`build()`**
  Assemble the `FileDataset`.

**Raises**

- **`ValueError`**
  If `fixed_series` is empty.

**See Also**

- [`get_deformable_registrations`](#regget-deformable-registrations) — Read the result back.
- [`SpatialRegistrationBuilder`](#regspatialregistrationbuilder) — When the transform is a single matrix.
- [`demons_registration`](#regmethoddemons-registration), [`bspline_registration`](#regmethodbspline-registration) — Produce the displacement transform.

**Notes**

Unlike Spatial Registration, this object maps **fixed to moving** (Registered RCS
to Source RCS), as specified by DICOM PS3.3 C.20.3.1.1. Do not invert the
resampling transform. The displacement grid is located in the fixed frame.
At a fixed grid point `x`, the DICOM equation is
`source = post(pre(x) + displacement(x))`: the vector is indexed on the original
registered grid, even when a pre-matrix is present.

A displacement grid is bulky -- three `float32` per voxel, so roughly 12 bytes where
the CT stores 2. Full-resolution fields on a whole CT run to hundreds of megabytes.

**Examples**

```python
>>> rigid = rigid_registration(fixed_image, moving_image)
>>> moving_rigid = sitk.Resample(
...     moving_image, fixed_image, rigid, sitk.sitkLinear, -1000.0,
...     moving_image.GetPixelID())
>>> _, deform, _ = demons_registration(fixed_image, moving_rigid)
>>> reg_ds = (
...     DeformableSpatialRegistrationBuilder(fixed_series)
...     .add_registration(
...         moving_series=moving_series,
...         vectorial_field_transform=deform,
...         pre_transform=np.eye(4).ravel().tolist(),
...         post_transform=affine_to_homogeneous_matrix(rigid).ravel().tolist(),
...     )
...     .build()
... )
```

#### `DeformableSpatialRegistrationBuilder.add_registration`

```python
add_registration(self, moving_series, vectorial_field_transform, pre_transform=None, post_transform=None, pre_transform_type: str = 'RIGID', post_transform_type: str = 'RIGID') -> 'DeformableSpatialRegistrationBuilder'
```

Register a moving series to the fixed series' Frame of Reference.

**Parameters**

- **`moving_series`** — *list[Dataset]*
  The series being registered.

- **`vectorial_field_transform`** — *sitk.DisplacementFieldTransform*
  The displacement field, converted to a Deformable Registration Grid.

- **`pre_transform, post_transform`** — *sequence or None, optional*
  16 numbers each, row-major 4x4, applied before and after the grid. Pass an
  identity matrix or None for either if unused (None omits the sequence).

- **`pre_transform_type, post_transform_type`** — *str*
  `"RIGID"`, `"RIGID_SCALE"` or `"AFFINE"`.

**Returns**

- **`DeformableSpatialRegistrationBuilder`**
  self, so calls chain.

#### `DeformableSpatialRegistrationBuilder.build`

```python
build(self) -> pydicom.dataset.FileDataset
```

Assemble the Deformable Spatial Registration object.

**Returns**

- **`FileDataset`**
  Ready to `save_as(path, enforce_file_format=True)`.

**Raises**

- **`ValueError`**
  If no registration has been added or the assembled dataset fails validation.

**Notes**

Unlike `SpatialRegistrationBuilder.build` there is no `include_identity`
parameter. The deformable IOD carries the fixed frame in
the top-level FrameOfReferenceUID rather than needing an identity item to mark it.

#### `DeformableSpatialRegistrationBuilder.set_instance_number`

```python
set_instance_number(self: ~BuilderT, instance_number: int) -> ~BuilderT
```

Set InstanceNumber (0020,0013).

**Parameters**

- **`instance_number`** — *int*
  Position within the series. Default 1. Give each object a distinct number when
  writing several registrations into one series -- InstanceNumber is Type 1 here,
  and some archives reject duplicates within a series.

**Returns**

- **`SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder`**
  self, so calls chain.

#### `DeformableSpatialRegistrationBuilder.set_uid_prefix`

```python
set_uid_prefix(self: ~BuilderT, uid_prefix: str) -> ~BuilderT
```

Set the root under which generated UIDs are minted.

**Parameters**

- **`uid_prefix`** — *str*
  Organisation UID root, dot-terminated, e.g. `"1.2.826.0.1.3680043.2.1125."`.
  Unlike the RTSTRUCT and dose builders, the registration builders do **not**
  read `DICOM_UID_PREFIX`; unset, they fall back to pydicom's default root.

**Returns**

- **`SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder`**
  self, so calls chain.

<a id="regget-spatial-registrations"></a>

### `reg.get_spatial_registrations`

```python
get_spatial_registrations(reg_ds: pydicom.dataset.Dataset, calc_inverse_matrix: bool = False) -> dict
```

Read the matrices out of a Spatial Registration object.

**Parameters**

- **`reg_ds`** — *Dataset*
  Spatial Registration object.

- **`calc_inverse_matrix`** — *bool, optional*
  Also populate the reverse direction, by numerically inverting each matrix, so the
  result can be indexed either way round. Default False.

**Returns**

- **`dict`**
  Nested by frame of reference, outer key fixed and inner key moving:

  ```python
  {fixed_frame_uid: {moving_frame_uid: [16 floats]}}
  ```

  Matrices within each MatrixSequence are composed in DICOM order (last @ ... @ first).
  Each returned matrix is a flat row-major 4x4 running **moving to fixed** -- DICOM's
  direction, and the inverse of what `sitk.Resample` wants. Empty if the object
  carries no FrameOfReferenceUID.

  Registrations whose matrix is missing or unreadable are skipped with a warning on
  the `pydicomrt.reg.parser` logger, not raised.

**See Also**

- [`SpatialRegistrationBuilder`](#regspatialregistrationbuilder) — Write this object.
- [`get_deformable_registrations`](#regget-deformable-registrations) — The deformable equivalent.

**Notes**

An object built by this library includes an identity item for the fixed frame itself,
so the fixed UID normally appears as its own inner key mapping to the identity matrix.

To *apply* a matrix, invert it and load it into a transform:

```python
inverse = np.linalg.inv(np.array(matrix, dtype=float).reshape(4, 4))
transform = sitk.AffineTransform(3)
transform.SetMatrix(inverse[:3, :3].ravel().tolist())
transform.SetTranslation(inverse[:3, 3].tolist())
```

**Examples**

```python
>>> reg = get_spatial_registrations(dcmread("registration.dcm"))
>>> matrix = reg[fixed_frame_uid][moving_frame_uid]
>>> len(matrix)
16
```

<a id="regget-deformable-registrations"></a>

### `reg.get_deformable_registrations`

```python
get_deformable_registrations(reg_ds: pydicom.dataset.Dataset) -> list
```

Read the deformation fields out of a Deformable Spatial Registration object.

**Parameters**

- **`reg_ds`** — *Dataset*
  Deformable Spatial Registration object.

**Returns**

- **`list of dict`**
  One entry per registered frame of reference:

  ```python
  [{"SourceFrameOfReferenceUID": str,
    "PreDeformationMatrixRegistration":  [16 floats],
    "PostDeformationMatrixRegistration": [16 floats],
    "DeformableRegistrationGrid": {
        "GridDimensions":          [nx, ny, nz],
        "GridResolution":          [dx, dy, dz],
        "ImagePositionPatient":    [x, y, z],
        "ImageOrientationPatient": [6 floats],
        "VectorGridData":          np.ndarray,   # (z, y, x, 3)
    }}]
  ```

  `VectorGridData` is reshaped into `(z, y, x, 3)`; the last axis holds the
  `(dx, dy, dz)` displacement in mm. Note the axis order is numpy's, the reverse
  of the `GridDimensions` it was reshaped from.

  The transform composes as pre-matrix, then grid, then post-matrix, and runs
  **fixed to moving** (Registered RCS to Source RCS). Absent pre/post matrices
  are returned as identity matrices. Only grid-bearing entries are returned;
  matrix-only entries and entries with malformed or missing grid elements are skipped.

**See Also**

- [`DeformableSpatialRegistrationBuilder`](#regdeformablespatialregistrationbuilder) — Write this object.
- [`get_spatial_registrations`](#regget-spatial-registrations) — The rigid/affine equivalent.

**Examples**

```python
>>> regs = get_deformable_registrations(dcmread("deformable.dcm"))
>>> field = regs[0]["DeformableRegistrationGrid"]["VectorGridData"]
>>> field.shape
(120, 512, 512, 3)
```

<a id="regcheck-spatial-reg-iod"></a>

### `reg.check_spatial_reg_iod`

```python
check_spatial_reg_iod(reg_ds: pydicom.dataset.Dataset) -> dict
```

Validate a Spatial Registration object against its CIOD.

**Parameters**

- **`reg_ds`** — *Dataset*
  Spatial Registration object.

**Returns**

- **`dict`**
  `{"result": bool, "content": list of str}`. `result` is True when conformant
  and `content` empty; otherwise one message per missing or empty required
  element.

**See Also**

- [`SpatialRegistrationBuilder`](#regspatialregistrationbuilder) — Produces conformant objects.
- [`check_deformable_reg_iod`](#regcheck-deformable-reg-iod) — The deformable equivalent.

**Notes**

Structural conformance only. A registration can pass here and still point the wrong
way -- the matrix direction is not something a field-presence check can see.

<a id="regcheck-deformable-reg-iod"></a>

### `reg.check_deformable_reg_iod`

```python
check_deformable_reg_iod(reg_ds: pydicom.dataset.Dataset) -> dict
```

Validate a Deformable Spatial Registration object against its CIOD.

**Parameters**

- **`reg_ds`** — *Dataset*
  Deformable Spatial Registration object.

**Returns**

- **`dict`**
  `{"result": bool, "content": list of str}`, as for
  `check_spatial_reg_iod`.

**See Also**

- [`DeformableSpatialRegistrationBuilder`](#regdeformablespatialregistrationbuilder) — Produces conformant objects.
- [`check_spatial_reg_iod`](#regcheck-spatial-reg-iod) — The rigid/affine equivalent.

**Notes**

Structural conformance only; the displacement field's values are not examined.

<a id="regaffine-to-homogeneous-matrix"></a>

### `reg.affine_to_homogeneous_matrix`

```python
affine_to_homogeneous_matrix(transform: SimpleITK.SimpleITK.AffineTransform) -> numpy.ndarray
```

Convert a SimpleITK affine transform to a 4x4 homogeneous matrix.

The bridge from a registration result to something DICOM can store.

**Parameters**

- **`transform`** — *sitk.AffineTransform*
  A 3D linear transform exposing `GetMatrix()` and `GetTranslation()`.
  A non-zero rotation centre is folded into the matrix offset automatically.

**Returns**

- **`np.ndarray`**
  4x4, `float64`, row-major: rotation in `[:3, :3]`, translation in `[:3, 3]`,
  bottom row `[0, 0, 0, 1]`.

**See Also**

- `to_centre_free_affine` — The equivalent conversion when a SimpleITK affine is needed.
- [`SpatialRegistrationBuilder.add_registration`](#regspatialregistrationbuilder) — Consumes the flattened result.

**Notes**

Direction is preserved, not flipped. Registration functions return *fixed to moving*
and Spatial REG stores *moving to fixed*, so for that IOD the usual call is
`affine_to_homogeneous_matrix(transform.GetInverse())`.

Flatten the matrix into 16 row-major values for the builder, for example:
`matrix.astype("float32").ravel().tolist()`.

<a id="regsitk-displacement-field-to-deformable-registration-grid"></a>

### `reg.sitk_displacement_field_to_deformable_registration_grid`

```python
sitk_displacement_field_to_deformable_registration_grid(transform: SimpleITK.SimpleITK.DisplacementFieldTransform) -> pydicom.dataset.Dataset
```

Convert a SimpleITK displacement field into a Deformable Registration Grid item.

**Parameters**

- **`transform`** — *sitk.DisplacementFieldTransform*
  The displacement field. Its geometry -- origin, spacing, direction and size --
  becomes the grid's, so the field defines the region the deformation is described
  over.

**Returns**

- **`Dataset`**
  One Deformable Registration Grid Sequence item, carrying GridDimensions,
  GridResolution, ImagePositionPatient, ImageOrientationPatient and VectorGridData.

**See Also**

- [`DeformableSpatialRegistrationBuilder.add_registration`](#regdeformablespatialregistrationbuilder) — Uses this internally.
- [`get_deformable_registrations`](#regget-deformable-registrations) — Read the result back.

**Notes**

VectorGridData is written as `float32` in DICOM's x-fastest order, the reverse of
the numpy `(z, y, x, 3)` layout the field is held in.

---

<a id="regmethod"></a>

## reg.method

Registration algorithms.

<a id="regmethoddemons-registration"></a>

### `reg.method.demons_registration`

```python
demons_registration(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image, resolution_staging: Sequence[int] = (8, 4, 1), iteration_staging: Sequence[int] = (10, 10, 10), isotropic_resample: bool = False, initial_displacement_field: Optional[SimpleITK.SimpleITK.Image] = None, regularization_kernel_mm: Union[float, Sequence[float]] = 1.5, smoothing_sigma_factor: float = 1.0, smoothing_sigmas: Union[float, Sequence[float], bool, NoneType] = False, default_value: Optional[float] = None, ncores: int = 1, interp_order: int = 2, verbose: bool = False) -> Tuple[SimpleITK.SimpleITK.Image, SimpleITK.SimpleITK.Transform, SimpleITK.SimpleITK.Image]
```

Deformable registration via Fast Symmetric-Forces Demons (SimpleITK) with multi-resolution strategy.

The function registers `moving_image` to `fixed_image` by estimating a dense displacement
field using a coarse-to-fine (image pyramid) schedule. It then applies the resulting
displacement-field transform to resample the moving image into the fixed image space.
Inputs may have different grids; moving data are sampled onto the fixed pyramid grid
through the current transform at each level. This does not estimate a rigid alignment.

**Parameters**

- **`fixed_image`** — *sitk.Image*
  Reference (target) image. Output is defined on this image grid.

- **`moving_image`** — *sitk.Image*
  Image to be warped into the fixed image space.

- **`resolution_staging`** — *Sequence[int or float], default = (8, 4, 1)*
  Pyramid schedule. Length of this sequence defines the number of pyramid levels.
  - If `isotropic_resample=False`: values are downsampling shrink factors (typically integers, e.g., [8, 4, 1]).
  - If `isotropic_resample=True`: values are target isotropic voxel sizes in millimeters (floats, e.g., [2.5, 1.0, 0.5]).

- **`iteration_staging`** — *Sequence[int], default = (20, 20, 20)*
  Number of iterations per pyramid level. Must match the length of `resolution_staging`.

- **`isotropic_resample`** — *bool, default = False*
  If True, interpret `resolution_staging` as isotropic voxel sizes (mm); otherwise as shrink factors.

- **`initial_displacement_field`** — *Optional[sitk.Image], default = None*
  Initial DVF (vector image). If provided, it overrides any `initial_transform` usage
  inside `multiscale_demons` and is resampled onto the fixed grid when needed.

- **`regularization_kernel_mm`** — *float or Sequence[float], default = 1.0*
  Standard deviation(s) in millimeters used by demons for smoothing the update and DVF.
  Scalar broadcasts to all axes; sequence length must be 1 or image dimension.

- **`smoothing_sigma_factor`** — *float, default = 1.0*
  Fallback factor to derive `smoothing_sigmas` per level if `smoothing_sigmas` is not provided.
  Each level's sigma = `resolution_staging[level] * smoothing_sigma_factor`.

- **`smoothing_sigmas`** — *float or Sequence[float] or bool, default = False*
  Per-level Gaussian sigmas (in physical units, mm) used for image smoothing
  prior to resampling at each pyramid level. This is separate from the demons
  regularization kernel (`regularization_kernel_mm`).
  - False/None: derived from `resolution_staging` and `smoothing_sigma_factor`
  - float: same sigma for all levels
  - sequence: length must match number of levels

- **`default_value`** — *Optional[float], default = None*
  Out-of-bounds value during final resampling. If None, use 0. For CT-like images
  (min <= -1000), default to -1000.

- **`ncores`** — *int, default = 1*
  Number of threads for demons filter.

- **`interp_order`** — *int, default = sitk.sitkLinear*
  Interpolator for resampling (e.g., sitk.sitkNearestNeighbor, sitk.sitkLinear, sitk.sitkBSpline).

- **`verbose`** — *bool, default = False*
  If True and a global callback `deformable_registration_command_iteration` is defined,
  prints per-iteration metric values.

**Returns**

- **`registered_image`** — *sitk.Image*
  The moving image warped into the fixed image space.

- **`output_transform`** — *sitk.Transform*
  DisplacementFieldTransform built from the final DVF.

- **`deformation_field`** — *sitk.Image*
  The final dense displacement field (vector image) aligned to `fixed_image`.

**Raises**

- **`ValueError`**
  On invalid arguments (length mismatches, non-positive values, dimension mismatch, etc.).

- **`RuntimeError`**
  If helper `multiscale_demons(...)` is not available in the current scope.

**Notes**

- Internally casts images to Float32 for demons, then casts result back to the original
  moving pixel type after resampling.
- Regularization sigmas are specified in *voxels* to SimpleITK; this function converts
  millimeter inputs to voxel units using the fixed image spacing.

<a id="regmethodrigid-registration"></a>

### `reg.method.rigid_registration`

```python
rigid_registration(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image, histogram_bins: int = 100, learning_rate: float = 2.0, iterations: int = 300, shrink_factors=(4, 2, 1), smoothing_sigmas=(2.0, 1.0, 0.0), optimizer: str = 'regular_step', relaxation_factor: float = 0.7, min_step_mm: float = 0.0001, convergence_minimum_value: float = 1e-06, convergence_window_size: int = 10, max_step_size_mm: float = 2.0) -> SimpleITK.SimpleITK.Transform
```

Rigidly align a moving image to a fixed image.

Estimates translation and rotation only -- no scaling, no shear -- by maximising
Mattes mutual information over a multi-resolution pyramid. Mutual information rather
than intensity difference, so CT-to-CBCT and CT-to-MR work without matching intensity
scales.

**Parameters**

- **`fixed_image`** — *sitk.Image*
  The reference. The result is expressed on this image's grid.

- **`moving_image`** — *sitk.Image*
  The image to align.

- **`histogram_bins`** — *int, optional*
  Bins for Mattes mutual information. Default 100.

- **`learning_rate`** — *float, optional*
  Initial optimizer step. Default 2.0.

- **`iterations`** — *int, optional*
  Maximum optimizer iterations per resolution level. Default 300.

- **`shrink_factors`** — *sequence of int or None, optional*
  Downsampling factor per level, coarsest first. `None` gives a single-resolution
  registration. Default `(4, 2, 1)`.

- **`smoothing_sigmas`** — *sequence of float or None, optional*
  Gaussian sigma in mm per level, matching `shrink_factors` element for element.
  Default `(2.0, 1.0, 0.0)`.

- **`optimizer`** — *{"regular_step", "gradient_descent"}, optional*
  Default `"regular_step"`. See Notes -- the default is slower and much more
  reliable.

- **`relaxation_factor`** — *float, optional*
  regular_step only. The step is multiplied by this whenever the gradient direction
  reverses. Default 0.7.

- **`min_step_mm`** — *float, optional*
  regular_step only. Stop once the step falls below this. Default 1e-4.

- **`convergence_minimum_value`** — *float, optional*
  gradient_descent only. Default 1e-6.

- **`convergence_window_size`** — *int, optional*
  gradient_descent only. Default 10.

- **`max_step_size_mm`** — *float, optional*
  gradient_descent only. Upper bound on a single step in mm; 0.0 leaves it
  unbounded. Default 2.0.

**Returns**

- **`sitk.Transform`**
  An affine transform holding rotation and translation only, pointing **fixed to
  moving** -- a *resampling* transform, the direction `sitk.Resample` expects.
  Invert it before storing in a DICOM Spatial REG, which runs moving to fixed.
  Keep its direction when using it as a Deformable REG post-matrix.

**Raises**

- **`ValueError`**
  If `shrink_factors` and `smoothing_sigmas` differ in length, or the optimizer
  name is unknown.

**See Also**

- [`registration_pipeline`](#regpipelineregistration-pipeline) — Adds preprocessing and an optional deformable stage.
- [`demons_registration`](#regmethoddemons-registration) — Deformable refinement, to run after this.
- [`SpatialRegistrationBuilder`](#regspatialregistrationbuilder) — Export the result as DICOM.

**Notes**

`"gradient_descent"` can return a *worse* alignment than its own initialisation and
still report convergence: it takes an estimated first step, overshoots into a flat
region, and the convergence window then sees the metric stop changing. On a real
CT/CBCT pair it finished 37 mm from the clinical registration, almost all of it
superior-inferior.

`"regular_step"` shrinks its step by `relaxation_factor` whenever the gradient
reverses, so it cannot run away from a good position. On the same pair it finished
7.5 mm out, with the residual concentrated where the metric itself disagrees with the
clinical answer. It costs roughly 8x the wall time -- about 9 minutes rather than 1 on
full-size volumes -- which is why the fast path is still selectable.

Registration is only bit-reproducible single-threaded: ITK reduces the metric in
thread-completion order, and on that real pair the thread count moved the answer by
more than 15 mm. Call
`sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)` when you need determinism.

Starting position is estimated with `CenteredTransformInitializer` on the two
images' geometry. Resampling both onto a common grid beforehand destroys that estimate
-- there is then nothing to centre -- so register in the images' own grids.

**Examples**

```python
>>> transform = rigid_registration(fixed_image, moving_image)
>>> registered = sitk.Resample(
...     moving_image, fixed_image, transform, sitk.sitkLinear, -1000.0)
```

A quick preview, at the cost of reliability:

```python
>>> transform = rigid_registration(
...     fixed_image, moving_image, optimizer="gradient_descent")
```

<a id="regmethodbspline-registration"></a>

### `reg.method.bspline_registration`

```python
bspline_registration(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image, fixed_mask: Optional[SimpleITK.SimpleITK.Image] = None, moving_mask: Optional[SimpleITK.SimpleITK.Image] = None, resolution_staging: Sequence[int] = (4, 2, 1), smoothing_sigmas: Sequence[float] = (4, 2, 0), sampling_rate: Optional[float] = None, optimizer: str = 'LBFGS', metric: str = 'mean_squares', initial_grid_spacing: float = 64, grid_scale_factors: Optional[Sequence[int]] = (1, 2, 4), interpolator: int = 23, default_value: Optional[float] = None, number_of_iterations: int = 20, isotropic_resample: bool = False, initial_isotropic_size: float = 1, histogram_bins: int = 30, verbose: bool = False, ncores: int = 1) -> Tuple[SimpleITK.SimpleITK.Image, SimpleITK.SimpleITK.Transform, SimpleITK.SimpleITK.Image]
```

Deformable registration with a B-spline transform.

Returns the same triple as `demons_registration` so the two are interchangeable.
The returned transform is a *resampling* transform: it maps a point on the fixed grid
back into the moving image.

**Parameters**

- **`fixed_image, moving_image`** — *sitk.Image*
  Images to register. Bring them into rough alignment first -- this stage refines, it
  does not search.

- **`fixed_mask, moving_mask`** — *sitk.Image, optional*
  Binary masks restricting where the metric is evaluated, in the fixed and moving
  image respectively.

- **`resolution_staging`** — *sequence[int]*
  Downsampling factor per level, coarsest first. The last entry should be 1 so the
  finest level sees the image at full resolution; the previous default stopped at 2,
  which meant it never did.

- **`smoothing_sigmas`** — *sequence[float]*
  Gaussian sigma per level, in **millimetres** -- not multiples of the shrink factor.

- **`sampling_rate`** — *float, optional*
  Fraction of voxels the metric samples. `None` (the default) selects ITK's dense
  path and is the only setting that registers well here -- any sampled configuration
  loses most of the accuracy regardless of rate. Reduce `number_of_iterations` or
  coarsen `initial_grid_spacing` for speed instead; see the module docstring.

- **`optimizer`** — *str*
  One of `OPTIMIZERS`. `"LBFGSB"` requires `grid_scale_factors=None`.

- **`metric`** — *str*
  One of `METRICS`.

- **`initial_grid_spacing`** — *float*
  Control point spacing at the coarsest level, in mm.

- **`grid_scale_factors`** — *sequence[int], optional*
  Mesh multiplier per level, so the grid refines as the resolution does. `None`
  keeps one grid throughout, which is what `"LBFGSB"` needs.

- **`interpolator`** — *int*
  Interpolator for the final resampling.

- **`default_value`** — *float, optional*
  Padding value. Inferred from the moving image when omitted (CT-like images use
  -1000).

- **`number_of_iterations`** — *int*
  Maximum optimizer iterations per level. Past roughly twenty the fit begins chasing
  noise, so raising this costs time and accuracy together -- see the module docstring.

- **`isotropic_resample`** — *bool*
  Resample both images to isotropic voxels before registering.

- **`initial_isotropic_size`** — *float*
  Voxel size in mm when `isotropic_resample` is set.

- **`histogram_bins`** — *int*
  Bins for `metric="mutual_information"`.

- **`verbose`** — *bool*
  Log per-iteration metric values.

- **`ncores`** — *int*
  Registration threads. Speedup is close to linear, but ITK reduces the metric in
  thread completion order, so only `ncores=1` is bit-reproducible.

**Returns**

- **`registered_image`** — *sitk.Image*
  Moving image resampled onto the fixed grid.

- **`output_transform`** — *sitk.Transform*
  Displacement field transform, fixed -> moving.

- **`deformation_field`** — *sitk.Image*
  The displacement field itself, in mm.

**Raises**

- **`ValueError`**
  If the optimizer or metric name is unknown, the resolution schedules disagree, or
  `"LBFGSB"` is combined with `grid_scale_factors`.

<a id="regmethoddemons-with-soft-mask"></a>

### `reg.method.demons_with_soft_mask`

```python
demons_with_soft_mask(fixed: SimpleITK.SimpleITK.Image, moving: SimpleITK.SimpleITK.Image, demons_fn, *, fixed_mask: SimpleITK.SimpleITK.Image, moving_mask: Optional[SimpleITK.SimpleITK.Image] = None, fixed_falloff_mm: float = 8.0, moving_falloff_mm: Optional[float] = None, background_fixed: Optional[float] = None, background_moving: Optional[float] = None, propagate_fixed_weight_to_moving: bool = False, initial_transform: Optional[SimpleITK.SimpleITK.Transform] = None, interp_order: int = 2, **demons_kwargs)
```

Run Demons registration weighted towards a region of interest.

Concentrates the registration on one structure without discarding its surroundings. A
hard mask would put a step edge at the boundary, which Demons reads as enormous
displacement; instead the mask is turned into a weight that falls off smoothly over
`fixed_falloff_mm`, and both images are blended towards a background value outside
it. Content well outside the ROI stops driving the result while the boundary stays
smooth.

**Parameters**

- **`fixed`** — *sitk.Image*
  The reference image.

- **`moving`** — *sitk.Image*
  The image to align.

- **`demons_fn`** — *callable*
  The registration to run, called as `demons_fn(weighted_fixed, weighted_moving,
  **demons_kwargs)`. Normally `demons_registration`. Taking it as an argument
  keeps the weighting independent of which Demons variant is used.

- **`fixed_mask`** — *sitk.Image*
  Required, keyword-only. Non-zero inside the ROI, on the `fixed` grid.

- **`moving_mask`** — *sitk.Image, optional*
  Same, on the `moving` grid. Without it the moving image is weighted only if
  `propagate_fixed_weight_to_moving` is set.

- **`fixed_falloff_mm`** — *float, optional*
  Distance in mm over which the weight decays from 1 inside the ROI to 0 outside.
  Default 8.0. Too small reintroduces the edge the soft mask exists to avoid; too
  large dilutes the focus.

- **`moving_falloff_mm`** — *float, optional*
  Falloff for `moving_mask`. Defaults to `fixed_falloff_mm`.

- **`background_fixed`** — *float, optional*
  Value the fixed image is blended towards outside the ROI. Inferred when omitted:
  -1000 if the image's minimum reaches -1000 (so it looks like CT), otherwise 0.

- **`background_moving`** — *float, optional*
  Same for the moving image, inferred the same way. Also becomes the resampling
  pad value handed to `demons_fn` unless `default_value` is given explicitly.

- **`propagate_fixed_weight_to_moving`** — *bool, optional*
  Warp the fixed weight map into moving space through `initial_transform` and use
  it as the moving weight. Useful when only the fixed image is segmented. Combined
  with `moving_mask`, the elementwise minimum is used. Default False.

- **`initial_transform`** — *sitk.Transform, optional*
  Required for propagation; identity when omitted, which only makes sense if the
  two images are already roughly aligned.

- **`interp_order`** — *int, optional*
  SimpleITK interpolator, passed on to `demons_fn`. Default `sitk.sitkLinear`.

- **`demons_kwargs`** — *dict, optional*
  Extra keyword arguments passed straight through to `demons_fn` --
  `iteration_staging`, `resolution_staging` and so on.

**Returns**

- **`tuple`**
  Whatever `demons_fn` returns -- for `demons_registration`, that is
  `(registered_image, transform, deformation_field)`.

**Raises**

- **`ValueError`**
  If `fixed_mask` is None.

**See Also**

- [`demons_registration`](#regmethoddemons-registration) — The unweighted version, and the usual `demons_fn`.

**Notes**

The returned image is a resample of the *weighted* moving image, so it carries the
background blending. Apply the returned transform to your original moving image when
you need untouched intensities.

The deformation outside the ROI is extrapolation, not measurement -- the metric was
weighted away from there. Do not read displacements far from the structure as
meaningful.

**Examples**

```python
>>> registered, transform, field = demons_with_soft_mask(
...     fixed_image, moving_image, demons_registration,
...     fixed_mask=prostate_mask, fixed_falloff_mm=10.0)
```

---

<a id="regpipeline"></a>

## reg.pipeline

Staged registration and transform composition.

<a id="regpipelineregistration-pipeline"></a>

### `reg.pipeline.registration_pipeline`

```python
registration_pipeline(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image, perform_rigid: bool = True, perform_deformable: bool = True, preprocess_config: Optional[Dict] = None, rigid_kwargs: Optional[Dict] = None, deformable_kwargs: Optional[Dict] = None, resample_interpolator: int = 2, default_value: Optional[float] = None) -> Tuple[SimpleITK.SimpleITK.Image, Optional[SimpleITK.SimpleITK.Transform], Optional[SimpleITK.SimpleITK.Transform], Optional[SimpleITK.SimpleITK.Image]]
```

Registration pipeline integrating rigid and deformable registration.

Stages:
1. Stage 0: Optional preprocessing before rigid (e.g. window clipping).
2. Stage 1: Rigid registration (if enabled).
3. Stage 2: Preprocessing before deform in rigid-aligned space (overlap crop + clip).
4. Stage 3: Deformable registration (if enabled).
5. Stage 4: Apply all transforms to the original moving_image.

**Parameters**

- **`fixed_image`** — *sitk.Image*
  Reference (fixed) image.

- **`moving_image`** — *sitk.Image*
  Image to be registered (moving image).

- **`perform_rigid`** — *bool, default = True*
  Whether to run rigid registration.

- **`perform_deformable`** — *bool, default = True*
  Whether to run deformable registration.

- **`preprocess_config`** — *Optional[Dict], default = None*
  Preprocessing config. Two formats supported:

  **Nested (recommended)**:
  {
  "rigid": {"window_clip": [-10, 500]},
  "deform": {
  "window_clip": [-10, 500],
  "align_extents": True,
  }
  }

  **Flat (backward compatible)**:
  {"window_clip": [-10, 500], "align_extents": True}


- **`rigid_kwargs`** — *Optional[Dict], default = None*
  Extra arguments for rigid_registration.

- **`deformable_kwargs`** — *Optional[Dict], default = None*
  Extra arguments for demons_registration.

- **`resample_interpolator`** — *int, default = sitk.sitkLinear*
  Interpolation for resampling.

- **`default_value`** — *Optional[float], default = None*
  Padding value for resampling. If None, inferred per image from its intensity range
  (CT-like images use -1000, everything else 0). Pass it explicitly for MR/PET.

**Returns**

- **`registered_image`** — *sitk.Image*
  Final registered image, always on the `fixed_image` grid and carrying the
  `moving_image` pixel type. Produced by a single resampling of the untouched
  `moving_image` through the composed transform, so no preprocessing (window
  clipping, extent cropping) leaks into the output.

- **`rigid_transform`** — *Optional[sitk.Transform]*
  Rigid transform, or None if rigid was not run.

- **`deformable_transform`** — *Optional[sitk.Transform]*
  Deformable transform, or None if deformable was not run.

- **`deformation_field`** — *Optional[sitk.Image]*
  Deformation field, or None if deformable was not run.

**Notes**

`rigid_transform` and `deformable_transform` are resampling transforms: they map a
point on the fixed grid back into the moving image. Compose them with
`compose_transforms` (stage order) rather than resampling once per stage.

The rigid stage registers the images in their own grids. Resampling them onto a shared
grid first -- which this used to do -- leaves `CenteredTransformInitializer` with two
identical grids and therefore nothing to estimate, so any offset larger than the
optimizer's own capture range is lost. With no preprocessing configured, the result is
identical to calling `rigid_registration` directly.

**Examples**

Rigid then deformable:

```python
>>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
...     fixed_image=ct_image, moving_image=cbct_image)
```

Rigid only:

```python
>>> registered, rigid_tfm, _, _ = registration_pipeline(
...     fixed_image=ct_image, moving_image=cbct_image,
...     perform_deformable=False)
```

Preprocessing applied to every stage (flat form):

```python
>>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
...     fixed_image=ct_image, moving_image=cbct_image,
...     preprocess_config={'window_clip': [-10, 500]})
```

Preprocessing configured per stage (nested form):

```python
>>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
...     fixed_image=ct_image, moving_image=cbct_image,
...     preprocess_config={
...         'rigid': {'window_clip': [-10, 500]},
...         'deform': {'window_clip': [-10, 500], 'align_extents': True},
...     })
```

<a id="regpipelinecompose-transforms"></a>

### `reg.pipeline.compose_transforms`

```python
compose_transforms(transforms: List[Optional[SimpleITK.SimpleITK.Transform]]) -> Optional[SimpleITK.SimpleITK.Transform]
```

Compose a stage-ordered list of transforms into a single resampling transform.

`transforms` is given in the order the stages were computed (rigid first, then
deformable). A resampling transform maps points from the *fixed* grid back into the
moving image, so the stages must be applied in reverse: the deformable transform maps
a fixed-space point into the rigid-aligned space, and the rigid transform then maps
that into the original moving space.

`SimpleITK.CompositeTransform` evaluates its list back-to-front, which matches
the stage order exactly, so the list is passed through unchanged.

**Parameters**

- **`transforms`** — *List[Optional[sitk.Transform]]*
  Stage-ordered transforms. `None` entries are skipped.

**Returns**

- **`Optional[sitk.Transform]`**
  A single transform equivalent to applying every stage, or None if the list holds
  no transform. A single-element list is returned as-is rather than wrapped.

<a id="regpipelinepreprocess-image"></a>

### `reg.pipeline.preprocess_image`

```python
preprocess_image(image: SimpleITK.SimpleITK.Image, preprocess_config: Optional[Dict] = None) -> SimpleITK.SimpleITK.Image
```

Preprocess the image by applying multiple preprocessing steps according to the config.

**Parameters**

- **`image`** — *sitk.Image*
  Input image.

- **`preprocess_config`** — *Optional[Dict], default = None*
  Preprocessing config dict. Supported options:
  - 'window_clip': Union[List[float], Tuple[float, float]] - Window clip range [lower, upper].

**Returns**

- **`sitk.Image`**
  Preprocessed image.

**Examples**

```python
>>> processed = preprocess_image(
...     image=ct_image,
...     preprocess_config={'window_clip': [-10, 500]})
```

<a id="regpipelineinfer-default-pixel-value"></a>

### `reg.pipeline.infer_default_pixel_value`

```python
infer_default_pixel_value(image: SimpleITK.SimpleITK.Image, default_value: Optional[float] = None) -> float
```

Resolve the padding value to use when resampling outside an image's extent.

An explicit `default_value` always wins. Otherwise the value is inferred from the
image intensities: an image reaching down to air (`<= -1000`) is treated as CT and
padded with `CT_AIR_HU`, anything else is padded with `0`.

Note that the inference only works on raw intensities. Once an image has been window
clipped its minimum no longer reaches `-1000`, so pass `default_value` explicitly
for preprocessed or non-CT (MR/PET) images instead of relying on the inference.

**Parameters**

- **`image`** — *sitk.Image*
  Image whose intensity range is inspected.

- **`default_value`** — *Optional[float], default = None*
  Explicit padding value. Returned unchanged when not None.

**Returns**

- **`float`**
  Padding value to hand to SimpleITK resampling.

<a id="regpipelinewindow-clip"></a>

### `reg.pipeline.window_clip`

```python
window_clip(image: SimpleITK.SimpleITK.Image, clip_range: Union[List[float], Tuple[float, float]]) -> SimpleITK.SimpleITK.Image
```

Apply window clipping (clamping) to the image.

**Parameters**

- **`image`** — *sitk.Image*
  Input image.

- **`clip_range`** — *Union[List[float], Tuple[float, float]]*
  Clipping range [lower, upper].

**Returns**

- **`sitk.Image`**
  Clipped image.

<a id="regpipelinealign-image-extents"></a>

### `reg.pipeline.align_image_extents`

```python
align_image_extents(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image, default_value: Optional[float] = None, use_initial_rigid: bool = True) -> Tuple[SimpleITK.SimpleITK.Image, SimpleITK.SimpleITK.Image, Optional[SimpleITK.SimpleITK.Transform]]
```

Align two images' physical extents by cropping them to their intersection.

If use_initial_rigid=True, an initial rigid transform is applied first to roughly align
the images, then intersection is computed and cropping is applied for a more accurate extent.

**Parameters**

- **`fixed_image`** — *sitk.Image*
  Reference image.

- **`moving_image`** — *sitk.Image*
  Image to be registered.

- **`default_value`** — *Optional[float], default = None*
  Default pixel value. If None, auto-detected (CT images use -1000).

- **`use_initial_rigid`** — *bool, default = True*
  Whether to use initial rigid transform to align images before computing intersection.

**Returns**

- **`fixed_aligned`** — *sitk.Image*
  Aligned fixed_image (cropped to intersection).

- **`moving_aligned`** — *sitk.Image*
  Aligned moving_image (cropped to intersection).

- **`initial_transform`** — *Optional[sitk.Transform]*
  Initial rigid transform if used, else None.

<a id="regpipelineget-image-physical-extent"></a>

### `reg.pipeline.get_image_physical_extent`

```python
get_image_physical_extent(image: SimpleITK.SimpleITK.Image) -> Tuple[numpy.ndarray, numpy.ndarray]
```

Compute the physical extent (bounding box) of the image in physical space.

**Parameters**

- **`image`** — *sitk.Image*
  Input image.

**Returns**

- **`min_corner`** — *np.ndarray*
  Minimum corner coordinates (3D) in physical space.

- **`max_corner`** — *np.ndarray*
  Maximum corner coordinates (3D) in physical space.

<a id="regpipelineget-intersection-extent"></a>

### `reg.pipeline.get_intersection_extent`

```python
get_intersection_extent(image1: SimpleITK.SimpleITK.Image, image2: SimpleITK.SimpleITK.Image) -> Tuple[Optional[numpy.ndarray], Optional[numpy.ndarray], bool]
```

Compute the intersection of two images' physical extents.

**Parameters**

- **`image1`** — *sitk.Image*
  First image.

- **`image2`** — *sitk.Image*
  Second image.

**Returns**

- **`min_corner`** — *Optional[np.ndarray]*
  Minimum corner of the intersection in physical space, or None if no intersection.

- **`max_corner`** — *Optional[np.ndarray]*
  Maximum corner of the intersection in physical space, or None if no intersection.

- **`has_intersection`** — *bool*
  Whether an intersection exists.

<a id="regpipelinecrop-image-to-extent"></a>

### `reg.pipeline.crop_image_to_extent`

```python
crop_image_to_extent(image: SimpleITK.SimpleITK.Image, min_corner: numpy.ndarray, max_corner: numpy.ndarray, default_value: Optional[float] = None) -> SimpleITK.SimpleITK.Image
```

Crop the image to the given physical extent.

**Parameters**

- **`image`** — *sitk.Image*
  Input image.

- **`min_corner`** — *np.ndarray*
  Minimum corner coordinates (3D) in physical space.

- **`max_corner`** — *np.ndarray*
  Maximum corner coordinates (3D) in physical space.

- **`default_value`** — *Optional[float], default = None*
  Default pixel value. If None, auto-detected (CT images use -1000).

**Returns**

- **`sitk.Image`**
  Image on a patient-axis-aligned grid covering the requested voxel-centre box.
  Oblique inputs are resampled; spacing is retained and direction becomes identity.

<a id="regpipelineget-initial-rigid-transform"></a>

### `reg.pipeline.get_initial_rigid_transform`

```python
get_initial_rigid_transform(fixed_image: SimpleITK.SimpleITK.Image, moving_image: SimpleITK.SimpleITK.Image) -> SimpleITK.SimpleITK.Transform
```

Get initial rigid transform using SimpleITK's CenteredTransformInitializer.

Used to roughly align two images before computing their intersection.

**Parameters**

- **`fixed_image`** — *sitk.Image*
  Reference image.

- **`moving_image`** — *sitk.Image*
  Image to be registered.

**Returns**

- **`sitk.Transform`**
  Initial rigid transform.

<a id="regpipelineget-image-center"></a>

### `reg.pipeline.get_image_center`

```python
get_image_center(image: SimpleITK.SimpleITK.Image) -> numpy.ndarray
```

Compute the center point of the image in physical space.

**Parameters**

- **`image`** — *sitk.Image*
  Input image.

**Returns**

- **`np.ndarray`**
  Center point coordinates (3D) in physical space.

<a id="regpipelineget-images-distance"></a>

### `reg.pipeline.get_images_distance`

```python
get_images_distance(image1: SimpleITK.SimpleITK.Image, image2: SimpleITK.SimpleITK.Image) -> float
```

Compute the Euclidean distance between the centers of two images.

**Parameters**

- **`image1`** — *sitk.Image*
  First image.

- **`image2`** — *sitk.Image*
  Second image.

**Returns**

- **`float`**
  Euclidean distance (mm) between the two image centers.

<a id="regpipelinecreate-reference-image-from-extent"></a>

### `reg.pipeline.create_reference_image_from_extent`

```python
create_reference_image_from_extent(min_corner: numpy.ndarray, max_corner: numpy.ndarray, spacing: Tuple[float, float, float], direction: Tuple[float, ...], pixel_id: int = 8) -> SimpleITK.SimpleITK.Image
```

Create a reference image from the given physical extent.

**Parameters**

- **`min_corner`** — *np.ndarray*
  Minimum corner coordinates (3D) in physical space.

- **`max_corner`** — *np.ndarray*
  Maximum corner coordinates (3D) in physical space.

- **`spacing`** — *Tuple[float, float, float]*
  Pixel spacing.

- **`direction`** — *Tuple[float, ...]*
  Direction matrix (9 elements).

- **`pixel_id`** — *int, default = sitk.sitkFloat32*
  Pixel type.

**Returns**

- **`sitk.Image`**
  Created reference image.

---

<a id="dose"></a>

## dose

RT Dose: build and read dose grids.

<a id="dosertdosebuilder"></a>

### `dose.RTDoseBuilder`

```python
RTDoseBuilder(reference_ds: pydicom.dataset.Dataset)
```

Build an RT Dose object from a dose grid.

Patient, study and frame-of-reference context is taken from a reference dataset -- a
slice of the planning CT, or the plan -- so the dose files into the same study.

**Parameters**

- **`reference_ds`** — *Dataset*
  A single dataset, not a series: one planning CT slice, or the RT Plan. Supplies
  patient identity, study and frame of reference.

**Methods**

- **`set_dose_grid(dose_image, ...)`**
  The dose values and their geometry. Required.

- **`add_referenced_plan(plan_ds)`**
  Reference the RT Plan. Required for most Dose Summation Types.

- **`set_series_number(series_number)`**
  Set SeriesNumber.

- **`set_uid_prefix(uid_prefix)`**
  Set the UID root.

- **`build()`**
  Assemble the `FileDataset`.

**Raises**

- **`ValueError`**
  If `reference_ds` is None.

**See Also**

- [`get_dose_image`](#doseget-dose-image) — Read the result back.
- [`check_rtdose_iod`](#dosecheck-rtdose-iod) — Validate, including the conditional plan reference.

**Notes**

Dose is stored as scaled unsigned 32-bit integers. Negative values are rejected --
unsigned storage would turn them into very large positive ones -- so clip or offset
beforehand if your grid dips below zero numerically.

**Examples**

```python
>>> dose_ds = (
...     RTDoseBuilder(planning_ct_series[0])
...     .set_dose_grid(dose_image)
...     .add_referenced_plan(plan_ds)
...     .build()
... )
>>> dose_ds.save_as("rtdose.dcm", enforce_file_format=True)
```

#### `RTDoseBuilder.add_referenced_plan`

```python
add_referenced_plan(self, plan_ds: pydicom.dataset.Dataset) -> 'RTDoseBuilder'
```

Reference the RT Plan this dose belongs to.

**Parameters**

- **`plan_ds`** — *Dataset*
  An RT Plan (SOP class `1.2.840.10008.5.1.4.1.1.481.5`). Only its SOP Class
  and Instance UIDs are stored.

**Returns**

- **`RTDoseBuilder`**
  self, so calls chain.

**Notes**

Referenced RT Plan Sequence is Type 1C, required for every Dose Summation Type
except `"FRACTION"` -- including the `"PLAN"` default. Omitting it produces a
file `check_rtdose_iod` reports and most TPSs refuse.

Can be called more than once; each plan is appended.

**Raises**

- **`ValueError`**
  If `plan_ds` is not an RT Plan.

#### `RTDoseBuilder.build`

```python
build(self) -> pydicom.dataset.FileDataset
```

Assemble the RT Dose object.

**Returns**

- **`FileDataset`**
  Ready to `save_as(path, enforce_file_format=True)`.

**Raises**

- **`ValueError`**
  If no dose grid has been set, or a referenced plan is not an RT Plan.

#### `RTDoseBuilder.set_dose_grid`

```python
set_dose_grid(self, dose_image: SimpleITK.SimpleITK.Image, dose_grid_scaling: float = None, dose_type: str = 'PHYSICAL', dose_summation_type: str = 'PLAN', dose_units: str = 'GY') -> 'RTDoseBuilder'
```

Set the dose grid and how it is described.

**Parameters**

- **`dose_image`** — *sitk.Image*
  Dose values in `dose_units`. Geometry -- origin, spacing, direction --
  becomes the grid's; it need not match the CT grid, and typically is coarser.

- **`dose_grid_scaling`** — *float, optional*
  Factor relating stored integers to real dose. Default 1e-7, covering about
  429 Gy. Raised automatically if the grid exceeds what the current factor can
  represent, so the default rarely needs changing; lower it for finer
  quantisation over a smaller range.

- **`dose_type`** — *str, optional*
  DoseType: `"PHYSICAL"` (default), `"EFFECTIVE"` or `"ERROR"`.

- **`dose_summation_type`** — *str, optional*
  DoseSummationType: `"PLAN"` (default), `"FRACTION"`, `"BEAM"`,
  `"BRACHY"` and so on. Everything except `"FRACTION"` makes Referenced RT
  Plan Sequence mandatory -- add one with `add_referenced_plan`.

- **`dose_units`** — *str, optional*
  DoseUnits: `"GY"` (default) or `"RELATIVE"`.

**Returns**

- **`RTDoseBuilder`**
  self, so calls chain.

**See Also**

- [`add_dose_grid_to_ds`](#doseadd-dose-grid-to-ds) — The underlying function.
- `add_referenced_plan` — Required for most summation types.

**Notes**

Calling this twice replaces the grid rather than accumulating; the last call wins.

#### `RTDoseBuilder.set_series_number`

```python
set_series_number(self, series_number: int) -> 'RTDoseBuilder'
```

Set SeriesNumber (0020,0011).

**Parameters**

- **`series_number`** — *int*
  Default 1. Give each dose its own number when writing several into one study
  -- beam doses alongside a plan dose, say -- so a viewer can tell them apart.

**Returns**

- **`RTDoseBuilder`**
  self, so calls chain.

#### `RTDoseBuilder.set_uid_prefix`

```python
set_uid_prefix(self, uid_prefix: str) -> 'RTDoseBuilder'
```

Set the root under which generated UIDs are minted.

**Parameters**

- **`uid_prefix`** — *str*
  Organisation UID root, dot-terminated. Defaults to the `DICOM_UID_PREFIX`
  environment variable.

**Returns**

- **`RTDoseBuilder`**
  self, so calls chain.

<a id="dosegenerate-base-dataset"></a>

### `dose.generate_base_dataset`

```python
generate_base_dataset(uid_prefix: str = None) -> pydicom.dataset.FileDataset
```

Create an empty RT Dose dataset with its required elements.

File meta, SOP class, the mandatory Type 1 and Type 2 elements, and empty sequences.
Carries no patient, study or dose data. Used by `RTDoseBuilder`.

**Parameters**

- **`uid_prefix`** — *str, optional*
  UID root. Defaults to the `DICOM_UID_PREFIX` environment variable.

**Returns**

- **`FileDataset`**
  An RT Dose skeleton.

**See Also**

- [`RTDoseBuilder`](#dosertdosebuilder) — The supported entry point.
- [`cp_information_from_ds`](#dosecp-information-from-ds) — Add patient and study context to the result.

<a id="dosecp-information-from-ds"></a>

### `dose.cp_information_from_ds`

```python
cp_information_from_ds(ds: pydicom.dataset.FileDataset, reference_ds: pydicom.dataset.Dataset, series_number: int = 1, uid_prefix: str = None)
```

Copy patient, study and frame-of-reference context into a dose dataset.

What makes the dose file into the same study as the images it was computed on, rather
than appearing as an unrelated study in the archive.

**Parameters**

- **`ds`** — *FileDataset*
  The dose dataset to populate, from `generate_base_dataset`. Modified in
  place.

- **`reference_ds`** — *Dataset*
  A planning CT slice or the RT Plan. Patient identity, StudyInstanceUID and
  FrameOfReferenceUID are read from it.

- **`series_number`** — *int, optional*
  SeriesNumber for the dose. Default 1.

- **`uid_prefix`** — *str, optional*
  UID root for the newly minted Series and SOP Instance UIDs. Defaults to the
  `DICOM_UID_PREFIX` environment variable.

**Returns**

- **`FileDataset`**
  `ds`, modified in place.

**See Also**

- [`RTDoseBuilder`](#dosertdosebuilder) — The supported entry point.

**Notes**

Study identity is inherited; series identity is new. A dose is a new series within an
existing study, so reusing the reference's SeriesInstanceUID would file it as more
slices of the CT.

<a id="doseadd-dose-grid-to-ds"></a>

### `dose.add_dose_grid_to_ds`

```python
add_dose_grid_to_ds(ds: pydicom.dataset.FileDataset, dose_sitk_image: SimpleITK.SimpleITK.Image, dose_grid_scaling: float = None, dose_type: str = 'PHYSICAL', dose_summation_type: str = 'PLAN', dose_units: str = 'GY')
```

Attach a dose grid to an RT Dose dataset.

**Parameters**

- **`ds`** — *FileDataset*
  Dataset from `generate_base_dataset`, with patient context already copied.

- **`dose_sitk_image`** — *sitk.Image*
  The dose grid. Values are in `dose_units`.

- **`dose_grid_scaling`** — *float, optional*
  Stored value -> dose factor. Defaults to 1e-7, raised automatically if the maximum
  dose would not otherwise fit in unsigned 32-bit.

- **`dose_type`** — *str*
  `"PHYSICAL"`, `"EFFECTIVE"` or `"ERROR"`. The previous default was EFFECTIVE,
  which means biologically weighted -- rarely what a plain dose grid holds.

- **`dose_summation_type`** — *str*
  `"PLAN"`, `"BEAM"`, `"FRACTION"` and so on. Several values make Referenced RT
  Plan Sequence required; see `check_dose_plan_reference`.

- **`dose_units`** — *str*
  `"GY"` or `"RELATIVE"`.

**Raises**

- **`ValueError`**
  If the grid holds negative dose, or is not 3D.

<a id="dosecheck-rtdose-iod"></a>

### `dose.check_rtdose_iod`

```python
check_rtdose_iod(dose_ds: pydicom.dataset.Dataset) -> dict
```

Validate an RT Dose dataset against its CIOD.

Also covers the conditional requirement most Dose Summation Types place on Referenced
RT Plan Sequence, which a field-presence check alone would miss.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object.

**Returns**

- **`dict`**
  `{"result": bool, "content": list of str}`. `result` is True when conformant
  and `content` empty; otherwise one message per problem.

**See Also**

- [`RTDoseBuilder`](#dosertdosebuilder) — Produces conformant objects.
- `check_dose_plan_reference` — The Type 1C check on its own.

**Notes**

The plan-reference rule is the one most commonly tripped: Referenced RT Plan Sequence
is required for every Dose Summation Type except `"FRACTION"`, and the default
Summation Type is `"PLAN"`. Structural conformance only -- dose values and geometry
are not checked.

<a id="doseget-dose-image"></a>

### `dose.get_dose_image`

```python
get_dose_image(dose_ds)
```

Read an RT Dose object into a `sitk.Image`.

The read path matching `RTDoseBuilder`: scaling applied, geometry attached.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object.

**Returns**

- **`sitk.Image`**
  Dose in DoseUnits (Gy for a conventional plan dose), on the dose grid -- which is
  typically coarser than the CT and need not cover the same extent. Use
  `resample_to_reference_image` with `default_value=0.0` to put it on the CT
  grid; the -1000 default would invent negative dose.

**See Also**

- [`RTDoseBuilder`](#dosertdosebuilder) — The write path.
- [`get_dose_array`](#doseget-dose-array) — Values without geometry.

**Examples**

```python
>>> dose_image = get_dose_image(dcmread("rtdose.dcm"))
>>> on_ct = resample_to_reference_image(ct_image, dose_image, default_value=0.0)
```

<a id="doseget-dose-array"></a>

### `dose.get_dose_array`

```python
get_dose_array(dose_ds)
```

Dose values as a numpy array, scaling applied.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object.

**Returns**

- **`np.ndarray`**
  `(frame, row, column)` in DoseUnits -- Gy for a conventional plan dose. Stored
  dose is scaled unsigned integers; multiplying by DoseGridScaling is what makes the
  numbers physical, so never read `dose_ds.pixel_array` directly.

**See Also**

- [`get_dose_image`](#doseget-dose-image) — The same values with geometry attached.

<a id="doseget-dose-spacing"></a>

### `dose.get_dose_spacing`

```python
get_dose_spacing(dose_ds)
```

Voxel size of an RT Dose grid, in SimpleITK order.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object carrying PixelSpacing and GridFrameOffsetVector.

**Returns**

- **`list of float`**
  `(x, y, z)` in mm. The in-plane pair is swapped out of DICOM's
  `[row, column]` order; `z` is the mean step in GridFrameOffsetVector.

**Warns**

- **`Logs a warning on the `pydicomrt.dose.sitk_transform` logger when the frame offsets`**

- **`are not uniformly spaced (standard deviation above 0.001 mm). A dose grid is allowed`**

- **`to be non-uniform in z; SimpleITK images are not, so the mean is used and the`**

- **`departure is worth knowing about.`**

**See Also**

- [`get_dose_image`](#doseget-dose-image) — Assembles this with the array, origin and direction.

<a id="doseget-dose-origin"></a>

### `dose.get_dose_origin`

```python
get_dose_origin(dose_ds)
```

Patient coordinates of the dose grid's first voxel.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object.

**Returns**

- **`list of float`**
  ImagePositionPatient, `[x, y, z]` in mm. Returned as pydicom's `DSfloat`
  values; wrap in `[float(v) for v in ...]` if you need plain floats.

<a id="doseget-dose-direction"></a>

### `dose.get_dose_direction`

```python
get_dose_direction(dose_ds)
```

Axis directions of an RT Dose grid, in SimpleITK convention.

**Parameters**

- **`dose_ds`** — *Dataset*
  RT Dose object carrying ImageOrientationPatient.

**Returns**

- **`np.ndarray`**
  3x3 with the axis vectors in **columns**. The slice direction is the cross product
  of the two stored vectors.

**Notes**

No orthonormality check is made here, unlike `get_slice_directions` for images.

---

<a id="ct"></a>

## ct

CT Image: emit a series from a volume.

<a id="ctctbuilder"></a>

### `ct.CTBuilder`

```python
CTBuilder(image_series: List[pydicom.dataset.Dataset])
```

Emit a CT Image series from a volume, borrowing patient and study context from a
reference series.

The output is DERIVED/SECONDARY: it is a rendering of a volume, not an acquisition.
It gets its own SeriesInstanceUID, so it cannot be mistaken for -- or merged into --
the acquired series it borrowed context from.

**Parameters**

- **`image_series`** — *list of Dataset*
  Reference slices supplying patient identity, study and frame of reference. Only
  the first is read for context unless `set_copy_all_attributes` is on.

**Methods**

- **`set_volume(sitk_image)`**
  The volume to slice, in HU. Required (or `set_volume_from_array`).

- **`set_volume_from_array(volume, origin, spacing, direction)`**
  The same, from numpy plus geometry.

- **`set_plane(plan_view)`**
  AXIAL, CORONAL or SAGITTAL. Default AXIAL.

- **`set_copy_all_attributes(copy_all_attributes)`**
  Inherit every attribute from the reference series.

- **`set_uid_prefix(uid_prefix)`**
  Set the UID root.

- **`build()`**
  Slice the volume into one `Dataset` per slice.

**Raises**

- **`ValueError`**
  If `image_series` is empty.

**See Also**

- [`check_ct_iod`](#ctcheck-ct-iod) — Validate the emitted slices.
- [`image_series_to_sitk_image`](#utilsimage-series-to-sitk-image) — The reverse direction.

**Notes**

Values are stored as unsigned 16-bit with a -1024 intercept, so the representable
range is -1024 to 64511 HU. Values outside it are clipped, not wrapped.

**Examples**

```python
>>> slices = (
...     CTBuilder(reference_series)
...     .set_volume(sitk_image)
...     .set_plane("AXIAL")
...     .build()
... )
>>> for i, s in enumerate(slices):
...     s.save_as(f"ct_{i:04d}.dcm", enforce_file_format=True)
```

#### `CTBuilder.build`

```python
build(self) -> List[pydicom.dataset.Dataset]
```

Slice the volume into a CT Image series.

**Returns**

- **`list[Dataset]`**
  One dataset per slice, ordered along the slice normal.

**Raises**

- **`ValueError`**
  If no volume has been set.

#### `CTBuilder.set_copy_all_attributes`

```python
set_copy_all_attributes(self, copy_all_attributes: bool) -> 'CTBuilder'
```

Copy every attribute from the reference series instead of building fresh.

**Parameters**

- **`copy_all_attributes`** — *bool*
  Default False, which writes only the elements the CT Image IOD requires.

**Returns**

- **`CTBuilder`**
  self, so calls chain.

**Notes**

Useful to inherit acquisition parameters -- KVP, exposure, convolution kernel --
that a fresh dataset would not carry. It also inherits everything else the source
had, including private tags and any burned-in identifiers, and those inherited
values now describe a volume they were not measured from. Geometry and pixel
elements are still overwritten from the volume.

#### `CTBuilder.set_plane`

```python
set_plane(self, plan_view: str) -> 'CTBuilder'
```

Choose the plane each slice lies in.

**Parameters**

- **`plan_view`** — *{"AXIAL", "CORONAL", "SAGITTAL"}*
  Case-insensitive. Named for the plane each slice lies *in*, so the axis the
  slices stack along is the one **not** in the name: axial stacks along z,
  coronal along y, sagittal along x. Default AXIAL.

**Returns**

- **`CTBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If `plan_view` is not one of the three.

**Notes**

Reslicing does not resample. A coronal series from an anisotropic axial volume
inherits the slice spacing as its in-plane row spacing, so it will look coarse in
one direction. Resample to isotropic first if that matters.

#### `CTBuilder.set_uid_prefix`

```python
set_uid_prefix(self, uid_prefix: str) -> 'CTBuilder'
```

Set the root under which generated UIDs are minted.

**Parameters**

- **`uid_prefix`** — *str*
  Organisation UID root, dot-terminated. Unlike the RTSTRUCT and dose builders,
  this one does **not** read `DICOM_UID_PREFIX`; it defaults to pydicom's
  root, which you do not own.

**Returns**

- **`CTBuilder`**
  self, so calls chain.

#### `CTBuilder.set_volume`

```python
set_volume(self, sitk_image: SimpleITK.SimpleITK.Image) -> 'CTBuilder'
```

Set the volume to slice.

**Parameters**

- **`sitk_image`** — *sitk.Image*
  3D. Intensities are treated as Hounsfield units and written with a -1024
  intercept, so values outside -1024 to 64511 HU are clipped. Geometry --
  origin, spacing, direction -- is taken from the image and becomes the output
  series' geometry.

**Returns**

- **`CTBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If the image is not 3D.

#### `CTBuilder.set_volume_from_array`

```python
set_volume_from_array(self, volume: numpy.ndarray, origin: Sequence[float], spacing: Sequence[float], direction: Sequence[float]) -> 'CTBuilder'
```

Set the volume from a numpy array plus its geometry.

**Parameters**

- **`volume`** — *np.ndarray*
  3D, `(slice, row, column)`. Treated as HU.

- **`origin`** — *sequence of float*
  Patient coordinates of voxel `[0, 0, 0]`, in mm.

- **`spacing`** — *sequence of float*
  `(x, y, z)` in mm -- SimpleITK's order, not DICOM's `[row, column]`.

- **`direction`** — *sequence of float*
  3x3 or flat 9-sequence, axis vectors in **columns**.

**Returns**

- **`CTBuilder`**
  self, so calls chain.

**Raises**

- **`ValueError`**
  If the volume is not 3D or any geometry argument is the wrong shape.

**See Also**

- `set_volume` — When you already have a `sitk.Image`.

<a id="ctcheck-ct-iod"></a>

### `ct.check_ct_iod`

```python
check_ct_iod(ct_ds: pydicom.dataset.Dataset) -> dict
```

Validate a CT Image dataset against its CIOD.

**Parameters**

- **`ct_ds`** — *Dataset*
  A single slice. CT is a single-frame IOD, so each slice is checked separately.

**Returns**

- **`dict`**
  `{"result": bool, "content": list of str}`. `result` is True when conformant
  and `content` empty; otherwise one message per missing or empty required
  element.

**See Also**

- [`CTBuilder`](#ctctbuilder) — Produces conformant slices.

**Notes**

Structural conformance only. Check every slice, not just the first: pixel-data
elements are per-slice, so a series can pass on slice 0 and fail later.

**Examples**

```python
>>> slices = CTBuilder(reference_series).set_volume(image).build()
>>> bad = [i for i, s in enumerate(slices)
...        if not check_ct_iod(s)["result"]]
```
