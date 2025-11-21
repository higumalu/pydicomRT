# pydicomRT
pydicomRT 是一個用於處理放射治療（Radiation Therapy）DICOM 檔案的 Python 函式庫。它提供建立、修改、解析與驗證 RTSTRUCT 資料集的工具，並支援 RTSTRUCT 與體素遮罩之間的轉換、空間/可變形註冊與劑量檔案的處理。與 `pydicom`、`numpy`、`SimpleITK` 等常用函式庫可以無縫整合。

---

## 專案目標

- **降低放射治療應用的開發門檻**  
  提供直觀的 API 與工具，讓研究人員與工程師不需深入 DICOM 標準，也能輕鬆處理放射治療相關的 DICOM 檔案。  

- **實現 Python 3D 函式庫與 pydicom 的無縫整合**  
  建立穩健的橋樑，讓常見的 3D 影像處理函式庫（例如 `numpy`、`SimpleITK`）能與 `pydicom` 協同運作，加速醫學影像與放射治療應用開發。  

---

## 功能特點

- 建構 RTSTRUCT 資料集並管理感興趣區域（ROI）  
- 在 3D numpy 遮罩與 DICOM 輪廓之間互轉  
- 驗證 RTSTRUCT、空間 REG、可變形 REG 以及 RTDOSE 資料集  
- 讀取與排序 DICOM 影像序列，並檢查方向與體素間距  
- 建立 RTDOSE 資料集，或將 RTDOSE DICOM 轉換為 SimpleITK 影像  
- 從 SimpleITK 轉換匯出空間（`REG`）與可變形空間（`REG-DR`）註冊資料集  
- 將可變形註冊網格解析為 NumPy 形變向量場  
- 提供 SimpleITK 建影、重採樣與註冊（`rigid`、`demons`、`bspline`、`soft_demons`）等輔助函式  
- 提供像素座標與病患座標之間的轉換工具  
- 提供 CT 模態 IOD 相關輔助  

---

## 快速連結

- 範例：`example/try_demon_reg.py`、`example/try_sort_dcms.py`  
- 可變形 REG 建置示範：`test/reg/df_reg_build_test.py`  
- RTSTRUCT API：`src/pydicomrt/rs`  
- 註冊 API：`src/pydicomrt/reg`  
- 工具集：`src/pydicomrt/utils`  
- 架構說明：`docs/architecture.md`

---

## 安裝

### 依賴套件

- Python >= 3.8  
- pydicom >= 2.0.0  
- numpy >= 1.26.4  
- opencv-python >= 4.10.0  
- scipy >= 1.10.3  
- simpleitk >= 2.5.0  

### 使用 pip 安裝

```bash
pip install pydicomrt
```

### 從原始碼安裝

```bash
git clone https://github.com/higumalu/pydicomRT.git
cd pydicomRT
pip install .
```

---

## 使用範例

### 建立 RTSTRUCT 資料集並新增 ROI

```python
import numpy as np
from pydicomrt.rs.make_contour_sequence import add_contour_sequence_from_mask3d
from pydicomrt.rs.add_new_roi import create_roi_into_rs_ds
from pydicomrt.rs.builder import create_rtstruct_dataset
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# 讀取 DICOM 影像序列
ds_list = load_sorted_image_series("path/to/dicom/images")

# 建立空的 RTSTRUCT 資料集
rs_ds = create_rtstruct_dataset(ds_list)

# 建立 ROI（感興趣區域）
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "CTV")

# 建立 3D 遮罩
mask = np.zeros((len(ds_list), 512, 512))
mask[100:200, 100:400, 100:400] = 1
mask[120:180, 200:300, 200:300] = 0

# 將 3D 遮罩寫入 RTSTRUCT
rs_ds = add_contour_sequence_from_mask3d(rs_ds, ds_list, 1, mask)

# 儲存 RTSTRUCT
rs_ds.save_as("path/to/output.dcm", write_like_original=False)
```

---

### 空間註冊（Rigid）並匯出 DICOM REG

以下範例示範如何使用 SimpleITK 計算兩組 CT 之間的剛體轉換，並儲存為 DICOM Spatial Registration（REG）物件。

```python
import numpy as np
import SimpleITK as sitk
from pydicomrt.utils.image_series_loader import load_sorted_image_series
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder
from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.type_transform import affine_to_homogeneous_matrix
from pydicomrt.reg.builder import SpatialRegistrationBuilder

# 讀取固定與移動的 CT 影像序列
fixed_ds = load_sorted_image_series("/path/to/CT_fixed")
moving_ds = load_sorted_image_series("/path/to/CT_moving")

# 轉換成 SimpleITK 影像
fixed_img = SimpleITKImageBuilder().from_ds_list(fixed_ds)
moving_img = SimpleITKImageBuilder().from_ds_list(moving_ds)

# 執行剛體註冊（回傳 sitk.Transform）
rigid_tfm = rigid_registration(fixed_img, moving_img)

# 將轉換轉為 DICOM 需要的 4x4 列主序矩陣
rigid_matrix = affine_to_homogeneous_matrix(rigid_tfm).astype(np.float32).ravel().tolist()

# 建立 REG 資料集並儲存
builder = SpatialRegistrationBuilder(fixed_ds)
builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")  # 可選，用以保持 UID 一致
builder.add_rigid_registration(moving_ds, rigid_matrix)
reg_ds = builder.build()
reg_ds.save_as("/path/to/output_reg.dcm", write_like_original=False)
```

注意：
- DICOM 以固定影像座標系中的 4x4 列主序矩陣儲存轉換，請確認方向與使用情境一致。
- 亦可透過環境變數 `DICOM_UID_PREFIX` 全域設定自訂 UID 根。  

---

### 可變形空間註冊（Demons）匯出

此流程先進行剛體對齊，再執行 Demons 可變形註冊，最後將形變場匯出為 DICOM Deformable Spatial Registration 資料集。

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

# 先確保體素網格一致
moving_img = resample_to_reference_image(fixed_img, moving_img)

# 剛體預對齊
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

# 執行 Demons 可變形註冊（回傳註冊影像、transform、位移場）
reg_img, deform_tfm, dvf = demons_registration(fixed_img, moving_rigid, verbose=False)

# 匯出為 DICOM Deformable Spatial Registration
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

`demons_registration` 會回傳 SimpleITK 的位移欄位 transform（`deform_tfm`）以及 DVF 影像（`dvf`），建置器會自動轉成 DICOM 所需的 `VectorGridData` 格式。

---

### 解析可變形註冊資料集

```python
from pydicom import dcmread
from pydicomrt.reg.parser import get_deformable_reg_list

reg_ds = dcmread("/path/to/output_dr.dcm")
reg_entries = get_deformable_reg_list(reg_ds)
field = reg_entries[0]["DeformableRegistrationGrid"]["VectorGridData"]
print(field.shape)  # (z, y, x, 3)，以毫米為單位的位移向量
```

---

### 以 SimpleITK 影像建立 RTDOSE 資料集

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

若要將 RTDOSE DICOM 轉回 SimpleITK 影像，可使用 `pydicomrt.dose.sitk_transform.get_dose_sitk_image`。

---

### 從 RTSTRUCT 資料集中提取輪廓

```python
from pydicomrt.rs.parser import get_roi_number_to_name, get_contour_dict

# 取得 ROI 編號與名稱的對應
roi_map = get_roi_number_to_name(rs_ds)
print(roi_map)  # 輸出: {1: 'CTV'}

# 取得輪廓資訊字典
ctr_dict = get_contour_dict(rs_ds)
```

---

### 驗證 RTSTRUCT 資料集

```python
from pydicomrt.rs.checker import check_rs_iod

# 檢查資料集是否符合 IOD 規範
result = check_rs_iod(rs_ds)
print(result)  # 輸出: {'result': True, 'content': []}
```

---

### 將 RTSTRUCT 轉換為 3D 遮罩

```python
from pydicomrt.rs.rs_to_volume import rtstruct_to_mask_dict, calc_image_series_affine_mapping
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# 讀取 DICOM 影像序列
ds_list = load_sorted_image_series("path/to/dicom/images")

# 計算仿射映射與遮罩體積形狀
affine_mapping, mask_volume_shape = calc_image_series_affine_mapping(ds_list)

# 轉換為 3D 遮罩字典
mask_dict = rtstruct_to_mask_dict(rs_ds, affine_mapping, mask_volume_shape)
```

---

## 模組結構

- **rs**：RTSTRUCT 相關功能  
  - `builder`：建立 RTSTRUCT 資料集  
  - `add_new_roi`：新增 ROI  
  - `make_contour_sequence`：建立輪廓序列  
  - `parser`：解析 RTSTRUCT 資料集  
  - `checker`：驗證 RTSTRUCT 資料集  
  - `rs_to_volume`：RTSTRUCT 與體素資料互轉  
  - `packer`：封裝輪廓資料  
  - `contour_process_method`：輪廓處理方法  
  - `rs_ds_iod`：RTSTRUCT IOD 定義  

- **reg**：空間 / 可變形註冊  
  - `builder`：建立 DICOM REG / Deformable REG 資料集  
  - `parser`：解析空間 / 可變形形變（例如 `get_deformable_reg_list`）  
  - `check`：驗證註冊資料集  
  - `method`：SimpleITK 註冊輔助函式（`rigid`、`bspline`、`demons`、`soft_demons`）  
  - `ds_reg_ds_iod`：可變形空間註冊 IOD 定義  
  - `s_reg_ds_iod`：空間註冊 IOD 定義  
  - `type_transform`：轉換工具（仿射→齊次、位移場→DICOM 格式）  

- **dose**：劑量功能  
  - `builder`：建立劑量資料集  
  - `sitk_transform`：劑量 DICOM 與 SimpleITK 影像互轉  
  - `dose_ds_iod`：劑量 IOD 定義  

- **ct**：CT 影像功能  
  - `ct_ds_iod`：CT IOD 定義  

- **utils**：實用工具  
  - `image_series_loader`：讀取與排序 DICOM 影像序列  
  - `coordinate_transform`：像素 / 病患座標轉換  
  - `validate_dcm_info`：驗證 DICOM 中繼資料  
  - `sitk_transform`：SimpleITK 轉換、`SimpleITKImageBuilder` 與重採樣工具  
  - `rs_from_altas`：由圖譜建立 RTSTRUCT  

---

## 貢獻
歡迎提交 Issue 與 Pull Request！

---

## 參考資源

- [SimpleITK](https://simpleitk.org/)
- [pydicom](https://pydicom.github.io/)
- [RT-Utils](https://github.com/qurit/rt-utils)
- [PlatiPy](https://github.com/pyplati/platipy)

---

## 許可證
本專案採用 MIT License，詳情請參考 `LICENSE` 檔案。

---

## 作者
- Higumalu (higuma.lu@gmail.com)
