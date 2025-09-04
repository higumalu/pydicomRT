# pydicomRT
pydicomRT 是一個用於處理放射治療結構集（RTSTRUCT）DICOM文件的Python庫。它提供了創建、修改和驗證RTSTRUCT數據集的功能，以及在RTSTRUCT和體積掩碼之間進行轉換的工具。

---

## 專案目標

- **降低放射治療應用的開發門檻**  
  提供直觀的 API 和工具，使研究人員與工程師能更輕鬆地處理與放射治療相關的 DICOM 檔案，而不需要深入了解複雜的 DICOM 標準。  

- **實現 Python 3D 函式庫與 pydicom 的無縫整合**  
  建立穩健的橋樑，讓常見的 Python 3D 影像處理函式庫（例如 `numpy`、`SimpleITK`）能與 `pydicom` 無縫協作，加速醫學影像與放射治療應用的開發。  


---

## 功能特點
- 創建RTSTRUCT數據集
- 添加和管理感興趣區域(ROI)
- 將3D掩碼轉換為DICOM輪廓
- 將DICOM輪廓轉換為3D掩碼
- 驗證RTSTRUCT數據集的合規性
- 處理和排序DICOM圖像序列
- 坐標轉換工具
- 創建和驗證DICOM劑量分佈
- 處理空間註冊數據
- 支持CT影像數據

## 安裝
### 依賴項
- Python >= 3.8
- pydicom >= 2.0.0
- numpy >= 1.26.4
- opencv-python >= 4.10.0
- scipy >= 1.10.3
- simpleitk >= 2.5.0

### 使用pip安裝
```bash
pip install pydicomrt
```

### 從源碼安裝
```bash
git clone https://github.com/higumalu/pydicomRT.git
cd pydicomRT
pip install .
```

## 使用示例
### 創建RTSTRUCT數據集並添加ROI

```python
import numpy as np
from pydicomrt.rs.make_contour_sequence import add_contour_sequence_from_mask3d
from pydicomrt.rs.add_new_roi import create_roi_into_rs_ds
from pydicomrt.rs.builder import create_rtstruct_dataset
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# 加載圖像序列
ds_list = load_sorted_image_series("path/to/dicom/images")

# 創建空的RTSTRUCT數據集
rs_ds = create_rtstruct_dataset(ds_list)

# 創建ROI（感興趣區域）
rs_ds = create_roi_into_rs_ds(rs_ds, [0, 255, 0], 1, "CTV", "CTV")

# 創建3D掩碼
mask = np.zeros((len(ds_list), 512, 512))
mask[100:200, 100:400, 100:400] = 1
mask[120:180, 200:300, 200:300] = 0

# 將3D掩碼添加到RTSTRUCT數據集
rs_ds = add_contour_sequence_from_mask3d(rs_ds, ds_list, 1, mask)

# 保存RTSTRUCT數據集
rs_ds.save_as("path/to/output.dcm", write_like_original=False)
```

### 從RTSTRUCT數據集中提取輪廓信息
```python
from pydicomrt.rs.parser import get_roi_number_to_name, get_contour_dict

# 獲取ROI映射
roi_map = get_roi_number_to_name(rs_ds)
print(roi_map)  # 輸出: {1: 'CTV'}

# 獲取輪廓字典
ctr_dict = get_contour_dict(rs_ds)
```

### 驗證RTSTRUCT數據集
```python
from pydicomrt.rs.checker import check_rs_iod

# 檢查RTSTRUCT數據集是否符合IOD規範
result = check_rs_iod(rs_ds)
print(result)  # 輸出: {'result': True, 'content': []}
```

### 將RTSTRUCT轉換為3D掩碼
```python
from pydicomrt.rs.rs_to_volume import rtstruct_to_mask_dict
from pydicomrt.utils.image_series_loader import load_sorted_image_series

# 加載圖像序列
ds_list = load_sorted_image_series("path/to/dicom/images")

# 將RTSTRUCT轉換為3D掩碼字典
mask_dict = rtstruct_to_mask_dict(rs_ds, ds_list)
```

## 模塊結構
- **rs**: 放射治療結構集相關功能
  - `builder`: 創建RTSTRUCT數據集
  - `add_new_roi`: 添加新的ROI
  - `make_contour_sequence`: 創建輪廓序列
  - `parser`: 解析RTSTRUCT數據集
  - `checker`: 驗證RTSTRUCT數據集
  - `rs_to_volume`: RTSTRUCT與體積數據轉換
  - `packer`: 輪廓數據打包
  - `contour_process_method`: 輪廓處理方法
  - `rs_ds_iod`: RTSTRUCT IOD定義

- **reg**: 空間註冊功能
  - `builder`: 創建註冊數據集
  - `parser`: 解析註冊數據集
  - `check`: 驗證註冊數據集
  - `method`: 使用SimpleITK的註冊方法
  - `ds_reg_ds_iod`: 可變形空間註冊IOD定義
  - `s_reg_ds_iod`: 空間註冊IOD定義
  - `type_transform`: 類型轉換

- **dose**: 劑量分佈功能
  - `builder`: 創建劑量數據集
  - `dose_ds_iod`: 劑量IOD定義

- **ct**: CT影像數據功能
  - `ct_ds_iod`: CT IOD定義

- **utils**: 實用工具
  - `image_series_loader`: 加載和排序DICOM圖像序列
  - `coordinate_transform`: 坐標轉換工具
  - `validate_dcm_info`: 驗證DICOM信息
  - `sitk_transform`: 使用SimpleITK的轉換
  - `rs_from_altas`: 從圖集創建RTSTRUCT

## 貢獻
歡迎提交問題和拉取請求。

## 許可證
此項目採用MIT許可證 - 請參閱 `LICENSE` 文件了解詳細信息。

## 作者
- Higumalu (higuma.lu@gmail.com)