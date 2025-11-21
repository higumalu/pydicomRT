# pydicomRT

**pydicomRT** is a Python library for handling Radiation Therapy DICOM files. It provides utilities to create, modify, parse, and validate RTSTRUCT datasets; convert between RTSTRUCT and volumetric masks; and handle spatial/deformable registration and dose. It integrates smoothly with `pydicom`, `numpy`, and `SimpleITK`.

---

## Project Goals

- **Lower the development barrier for RT applications**  
  Provide intuitive APIs and tools that allow researchers and engineers to work with radiation therapy–related DICOM files more easily, without requiring deep knowledge of the complex DICOM standard.  

- **Enable seamless integration between Python 3D libraries and pydicom**  
  Build a robust bridge so that common Python 3D image processing libraries (e.g., `numpy`, `SimpleITK`) can work seamlessly with `pydicom`, accelerating medical imaging and radiotherapy application development.  

---

## Features

- Build RTSTRUCT datasets and manage Regions of Interest (ROIs)  
- Convert between 3D numpy masks and DICOM contours  
- Validate RTSTRUCT, Spatial REG, Deformable REG, and RTDOSE datasets  
- Load and sort DICOM image series with orientation/spacing sanity checks  
- Generate RTDOSE datasets and convert RTDOSE DICOM files to SimpleITK images  
- Export Spatial (`REG`) and Deformable Spatial (`REG-DR`) registrations from SimpleITK transforms  
- Parse deformable registration grids into NumPy displacement fields  
- SimpleITK helpers for image building, resampling, and registration (`rigid`, `demons`, `bspline`, `soft_demons`)  
- Coordinate transformation utilities between pixel and patient spaces  
- CT modality IOD helpers  

---

## Quick Links

- Examples: `example/try_demon_reg.py`, `example/try_sort_dcms.py`  
- Deformable REG builder demo: `test/reg/df_reg_build_test.py`  
- RTSTRUCT API: `src/pydicomrt/rs`  
- Registration API: `src/pydicomrt/reg`  
- Utilities: `src/pydicomrt/utils`  
- Architecture Doc: `docs/architecture.md`

---

## Installation

### Dependencies

- Python >= 3.8  
- pydicom >= 2.0.0  
- numpy >= 1.26.4  
- opencv-python >= 4.10.0  
- scipy >= 1.10.3  
- simpleitk >= 2.5.0  

### Install via pip

```bash
pip install pydicomrt
```

### Install from source

```bash
git clone https://github.com/higumalu/pydicomRT.git
cd pydicomRT
pip install .
```

---

## Usage Examples

### Create an RTSTRUCT Dataset and Add ROI

```python
import numpy as np
from pydicomrt.rs.make_contour_sequence import add_contour_sequence_from_mask3d
from pydicomrt.rs.add_new_roi import create_roi_into_rs_ds
from pydicomrt.rs.builder import create_rtstruct_dataset
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# Load DICOM image series
ds_list = load_sorted_image_series("path/to/dicom/images")

# Create an empty RTSTRUCT dataset
rs_ds = create_rtstruct_dataset(ds_list)

# Create an ROI (Region of Interest)
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "CTV")

# Create a 3D mask
mask = np.zeros((len(ds_list), 512, 512))
mask[100:200, 100:400, 100:400] = 1
mask[120:180, 200:300, 200:300] = 0

# Add 3D mask to RTSTRUCT dataset
rs_ds = add_contour_sequence_from_mask3d(rs_ds, ds_list, 1, mask)

# Save the RTSTRUCT dataset
rs_ds.save_as("path/to/output.dcm", write_like_original=False)
```

---

### Spatial Registration (Rigid) and DICOM REG Export

This example estimates a rigid transform between two CT series using SimpleITK and stores it in a DICOM Spatial Registration (REG) object.

```python
import numpy as np
import SimpleITK as sitk
from pydicomrt.utils.image_series_loader import load_sorted_image_series
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder
from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.type_transform import affine_to_homogeneous_matrix
from pydicomrt.reg.builder import SpatialRegistrationBuilder

# Load CT series as pydicom datasets
fixed_ds = load_sorted_image_series("/path/to/CT_fixed")
moving_ds = load_sorted_image_series("/path/to/CT_moving")

# Convert to SimpleITK images
fixed_img = SimpleITKImageBuilder().from_ds_list(fixed_ds)
moving_img = SimpleITKImageBuilder().from_ds_list(moving_ds)

# Run rigid registration in physical space (returns sitk.Transform)
rigid_tfm = rigid_registration(fixed_img, moving_img)

# Convert to 4x4 row-major list for DICOM REG
tfm_4x4 = affine_to_homogeneous_matrix(rigid_tfm).astype(np.float32).ravel().tolist()

# Build a DICOM Spatial Registration dataset and save
builder = SpatialRegistrationBuilder(fixed_ds)
builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")  # Optional but helps keep UIDs consistent
builder.add_rigid_registration(moving_ds, tfm_4x4)
reg_ds = builder.build()
reg_ds.save_as("/path/to/output_reg.dcm", write_like_original=False)
```

Notes:
- DICOM stores transforms as a 4x4 row-major matrix in the fixed image frame. Ensure transform directions match your use-case.
- You can set a custom UID root globally via the `DICOM_UID_PREFIX` environment variable.

---

### Deformable Spatial Registration (Demons) Export

This workflow performs a rigid pre-align, runs demons registration, and exports the resulting displacement field into a DICOM Deformable Spatial Registration dataset.

```python
import numpy as np
import SimpleITK as sitk
from pydicomrt.utils.image_series_loader import load_sorted_image_series
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder, resample_to_reference_image
from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.method.demons import demons_registration
from pydicomrt.reg.type_transform import affine_to_homogeneous_matrix
from pydicomrt.reg.builder import DeformableSpatialRegistrationBuilder

fixed_ds = load_sorted_image_series("/path/to/CT_fixed")
moving_ds = load_sorted_image_series("/path/to/CT_moving")
fixed_img = SimpleITKImageBuilder().from_ds_list(fixed_ds)
moving_img = SimpleITKImageBuilder().from_ds_list(moving_ds)

# Ensure voxel grids match before registration
moving_img = resample_to_reference_image(fixed_img, moving_img)

# Rigid pre-alignment
rigid_tfm = rigid_registration(fixed_img, moving_img)
rigid_matrix = affine_to_homogeneous_matrix(rigid_tfm).astype(np.float32).ravel().tolist()
moving_rigid = sitk.Resample(
    moving_img,
    fixed_img,
    rigid_tfm,
    sitk.sitkLinear,
    -1000.0,
    moving_img.GetPixelIDValue(),
)

# Demons deformable registration (returns registered image, transform, displacement field)
reg_img, deform_tfm, dvf = demons_registration(fixed_img, moving_rigid, verbose=False)

# Export to DICOM Deformable Spatial Registration
identity = np.eye(4, dtype=np.float32).ravel().tolist()

builder = DeformableSpatialRegistrationBuilder(fixed_ds)
builder.add_deformable_registration(
    moving_ds_list=moving_ds,
    vectorial_field_transform=deform_tfm,
    pre_transform=rigid_matrix,
    post_transform=identity,
)
reg_dr = builder.build()
reg_dr.save_as("/path/to/output_dr.dcm", write_like_original=False)
```

`demons_registration` returns a SimpleITK displacement-field transform (`deform_tfm`) and a DVF image (`dvf`). The builder converts the transform into the DICOM `VectorGridData` representation automatically.

---

### Parse a Deformable Registration Dataset

```python
from pydicom import dcmread
from pydicomrt.reg.parser import get_deformable_reg_list

reg_ds = dcmread("/path/to/output_dr.dcm")
reg_entries = get_deformable_reg_list(reg_ds)
field = reg_entries[0]["DeformableRegistrationGrid"]["VectorGridData"]
print(field.shape)  # (z, y, x, 3) float32 displacement vectors in mm
```

---

### Create an RTDOSE Dataset from a SimpleITK Image

```python
import SimpleITK as sitk
from pydicom import dcmread
from pydicomrt.dose.builder import generate_base_dataset, cp_information_from_ds, add_dose_grid_to_ds

reference = dcmread("path/to/reference_ct_or_plan.dcm")
dose_img = sitk.ReadImage("path/to/dose_image.nii.gz")

dose_ds = generate_base_dataset()
dose_ds = cp_information_from_ds(dose_ds, reference)
dose_ds = add_dose_grid_to_ds(dose_ds, dose_img)

dose_ds.save_as("path/to/output_dose.dcm", write_like_original=False)
```

To convert an RTDOSE DICOM dataset back to SimpleITK, use `pydicomrt.dose.sitk_transform.get_dose_sitk_image`.

---

### Extract Contour Information from RTSTRUCT Dataset

```python
from pydicomrt.rs.parser import get_roi_number_to_name, get_contour_dict

# Get ROI mapping
roi_map = get_roi_number_to_name(rs_ds)
print(roi_map)  # Output: {1: 'CTV'}

# Get contour dictionary
ctr_dict = get_contour_dict(rs_ds)
```

---

### Validate RTSTRUCT Dataset

```python
from pydicomrt.rs.checker import check_rs_iod

# Check whether the RTSTRUCT dataset conforms to IOD specification
result = check_rs_iod(rs_ds)
print(result)  # Output: {'result': True, 'content': []}
```

---

### Convert RTSTRUCT to 3D Mask

```python
from pydicomrt.rs.rs_to_volume import rtstruct_to_mask_dict, calc_image_series_affine_mapping
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# Load DICOM image series
ds_list = load_sorted_image_series("path/to/dicom/images")

# Calculate affine mapping and mask volume shape
affine_mapping, mask_volume_shape = calc_image_series_affine_mapping(ds_list)

# Convert RTSTRUCT to 3D mask dictionary
mask_dict = rtstruct_to_mask_dict(rs_ds, affine_mapping, mask_volume_shape)
```

---

## Module Structure

- **rs**: RTSTRUCT-related functionalities  
  - `builder`: Create RTSTRUCT datasets  
  - `add_new_roi`: Add new ROIs  
  - `make_contour_sequence`: Create contour sequences  
  - `parser`: Parse RTSTRUCT datasets  
  - `checker`: Validate RTSTRUCT datasets  
  - `rs_to_volume`: Convert between RTSTRUCT and volume data  
  - `packer`: Pack contour data  
  - `contour_process_method`: Contour processing methods  
  - `rs_ds_iod`: RTSTRUCT IOD definitions  

- **reg**: Spatial/deformable registration  
  - `builder`: Build DICOM REG / Deformable REG datasets  
  - `parser`: Extract spatial/deformable transforms (e.g., `get_deformable_reg_list`)  
  - `check`: Validate registration datasets  
  - `method`: SimpleITK registration helpers (`rigid`, `bspline`, `demons`, `soft_demons`)  
  - `ds_reg_ds_iod`: Deformable spatial registration IOD definitions  
  - `s_reg_ds_iod`: Spatial registration IOD definitions  
  - `type_transform`: Transform conversions (affine → homogeneous, displacement fields → DICOM grids)  

- **dose**: Dose distribution functionalities  
  - `builder`: Create dose datasets  
  - `sitk_transform`: Convert dose datasets to SimpleITK images  
  - `dose_ds_iod`: Dose IOD definitions  

- **ct**: CT image data functionalities  
  - `ct_ds_iod`: CT IOD definitions  

- **utils**: Utility tools  
  - `image_series_loader`: Load and sort DICOM image series  
  - `coordinate_transform`: Coordinate transformation utilities  
  - `validate_dcm_info`: Validate DICOM metadata  
  - `sitk_transform`: SimpleITK conversions, `SimpleITKImageBuilder`, and resampling helpers  
  - `rs_from_altas`: Build RTSTRUCT datasets from atlas inputs  

---

## Contributing

Issues and pull requests are welcome!

---

## Reference

- [SimpleITK](https://simpleitk.org/)
- [pydicom](https://pydicom.github.io/)
- [RT-Utils](https://github.com/qurit/rt-utils)
- [PlatiPy](https://github.com/pyplati/platipy)

---

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

---

## Author

- Higumalu (higuma.lu@gmail.com)
