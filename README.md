# pydicomRT

**pydicomRT** is a Python library for handling Radiation Therapy DICOM files. It provides capabilities for creating, modifying, and validating RTSTRUCT datasets, as well as tools for converting between RTSTRUCT and volumetric masks. Additionally, it supports handling other radiation therapy related DICOM files such as dose distributions and registration data.

---

## Features

- Create RTSTRUCT datasets  
- Add and manage Regions of Interest (ROIs)  
- Convert 3D masks to DICOM contours  
- Convert DICOM contours to 3D masks  
- Validate RTSTRUCT dataset compliance  
- Handle and sort DICOM image series  
- Coordinate transformation utilities  
- Create and validate DICOM dose distributions  
- Handle spatial registration data  
- Support for CT image data  

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

- **reg**: Spatial registration functionalities  
  - `builder`: Create registration datasets  
  - `parser`: Parse registration datasets  
  - `check`: Validate registration datasets  
  - `sitk_reg_method`: Registration methods using SimpleITK  

- **dose**: Dose distribution functionalities  
  - `builder`: Create dose datasets  

- **ct**: CT image data functionalities  
  - `ct_ds_iod`: CT IOD definitions  

- **utils**: Utility tools  
  - `image_series_loader`: Load and sort DICOM image series  
  - `coordinate_transform`: Coordinate transformation utilities  
  - `validate_dcm_info`: Validate DICOM metadata  
  - `sitk_transform`: Transformations using SimpleITK  

---

## Contributing

Issues and pull requests are welcome!

---

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

---

## Author

- Higumalu (higuma.lu@gmail.com)
