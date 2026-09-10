# pydicomRT

以 Python 慣用的資料型別（`numpy` 與 `SimpleITK`）讀寫放射治療 DICOM 物件 —— RTSTRUCT、
空間與可變形註冊（Registration）、RTDOSE、CT Image。

目標是讓常見工作（把遮罩轉成輪廓、把輪廓轉回遮罩、匯出註冊結果）不必先讀懂 DICOM 標準，
同時把那些一旦弄錯就會**無聲地毀掉結果**的座標慣例講清楚。

**English: [README.md](README.md) · 搭配 AI coding agent 使用：[AGENTS.md](AGENTS.md)**

這一頁是導覽。逐項工作的完整手冊 —— 每個參數、每則錯誤訊息，以及 0.8 → 0.9 的遷移對照表 ——
請看 **[docs/user-guide_zh.md](docs/user-guide_zh.md)**。每個函式的逐項參考（由 docstring 生成）
請看 **[docs/api-reference.md](docs/api-reference.md)**。

---

## 安裝

```bash
pip install pydicomrt
```

從原始碼安裝：

```bash
git clone https://github.com/higumalu/pydicomRT.git
cd pydicomRT
pip install -e ".[dev]"
```

需要 Python ≥ 3.10，以及 `pydicom`、`numpy`、`SimpleITK`、`opencv-python`、`scipy`。
本函式庫是針對 pydicom 3 撰寫的。

---

## 能做什麼

| 工作 | 進入點 |
|---|---|
| 讀取並排序 DICOM 影像序列 | `utils.load_sorted_image_series` |
| DICOM 序列 → `sitk.Image` | `utils.sitk_transform.SimpleITKImageBuilder` |
| 從遮罩建立 RTSTRUCT | `rs.RTStructBuilder` |
| RTSTRUCT → 3D 遮罩 | `rs.rtstruct_to_masks` |
| 驗證 RTSTRUCT / REG / 劑量 / CT | `rs.check_rtstruct_iod`、`reg.check_spatial_reg_iod`、`dose.check_rtdose_iod`、`ct.check_ct_iod` |
| 剛體 / 可變形註冊 | `reg.method.rigid_registration`、`reg.method.demons_registration` |
| 將註冊結果匯出成 DICOM REG | `reg.SpatialRegistrationBuilder`、`reg.DeformableSpatialRegistrationBuilder` |
| 讀取 DICOM REG | `reg.get_spatial_registrations`、`reg.get_deformable_registrations` |
| 建立 / 讀回 RTDOSE | `dose.RTDoseBuilder`、`dose.get_dose_image` |
| 從體資料產生 CT 切片 | `ct.CTBuilder` |
| 降噪、N4 偏壓場校正 | `utils.sitk_image_process` |

---

## 你必須知道的慣例

慣例 1–4 是出錯時**不會拋例外**的地方 —— 程式跑得過、影像看起來正常，只有物理座標是錯的。
請在看範例之前先讀這一節。

### 0. 所有 builder 的用法都一樣

五個物件 builder 共用同一套語彙，學會一個就等於學會其餘四個：

```python
Builder(參考來源)          # 病患／檢查資訊從影像繼承
    .set_something(...)   # 單一值的設定
    .add_something(...)   # 可重複的內容
    .build()              # -> FileDataset；CT 回傳 list[FileDataset]
```

| Builder | 建構參數 | 內容加入方式 |
|---|---|---|
| `rs.RTStructBuilder` | 影像序列 | `.add_roi(mask, name, ...)` |
| `dose.RTDoseBuilder` | 參考資料集 | `.set_dose_grid(...)`、`.add_referenced_plan(...)` |
| `ct.CTBuilder` | 影像序列 | `.set_volume(...)`、`.set_plane(...)` |
| `reg.SpatialRegistrationBuilder` | fixed 序列 | `.add_registration(moving, matrix)` |
| `reg.DeformableSpatialRegistrationBuilder` | fixed 序列 | `.add_registration(moving, dvf, pre, post)` |

`set_*` 與 `add_*` 都回傳 builder 本身，因此可以鏈式呼叫。`build()` 會驗證並拋錯，
而不是產生一個缺少 Type 1 欄位的物件。

驗證函式也是同一個形狀：`check_rtstruct_iod`、`check_spatial_reg_iod`、
`check_deformable_reg_iod`、`check_rtdose_iod`、`check_ct_iod` 全部回傳
`{"result": bool, "content": [...]}`。

### 1. 軸順序會翻轉三次

| 對象 | 順序 |
|---|---|
| 由序列堆疊出的 numpy 體資料 | `(slice, row, column)` = `(z, y, x)` |
| `sitk.Image` 的 size / spacing / origin | `(x, y, z)` |
| DICOM `PixelSpacing` | `[row spacing, column spacing]` = `(y, x)` |

`PixelSpacing[0]` 是**列與列之間**的距離，所以它屬於 **y** 軸而非 x 軸。
傳給 `RTStructBuilder.add_roi()` 的遮罩是 `(slice, row, column)`。

### 2. 註冊回傳的 transform 是 fixed → moving

本函式庫所有註冊函式回傳的都是**重採樣用的 transform**：它把 *fixed* 網格上的點映射回
*moving* 影像。這正是 `sitk.Resample` 需要的方向，所以可以直接用：

```python
transform = rigid_registration(fixed_image, moving_image)
registered = sitk.Resample(moving_image, fixed_image, transform, sitk.sitkLinear, -1000.0)
```

但 DICOM Spatial REG 的矩陣方向**相反** —— 是 moving → fixed。匯出前必須反轉：

```python
matrix = affine_to_homogeneous_matrix(transform.GetInverse())
```

兩個方向產生的檔案都能正常載入，所以這個錯誤只會表現為「位移方向相反」，
往往要到治療機上才被發現。參見 `example/02_rigid_registration_to_reg.py`。

### 3. 先合成 transform，只重採樣一次

`sitk.CompositeTransform([a, b])` 是**先套用 b**。因此依階段順序排列的清單不需要反轉 ——
請用 `compose_transforms([rigid, deformable])`，並且只呼叫**一次** `Resample`。
若每個階段各重採樣一次，每次都會裁切到 fixed 網格，前一階段推出視野的內容就永久變成 padding，
後續階段再也救不回來。

### 4. Padding 值不是隨便填的

重採樣到來源影像之外時會憑空產生體素。CT 的正確值是 **-1000**（空氣）。
若用 0 填充，等於在病人外面包了一層軟組織，相似度指標會努力去對齊這層假影。
`infer_default_pixel_value()` 會依強度最小值推測，但這只在**原始強度**上成立 ——
經過 window clipping 的 CT 就不再像 CT 了。MR、PET 或已前處理的影像請明確傳入 `default_value`。

---

## 範例

`example/` 底下每個腳本都能獨立執行（使用合成資料，不需任何病人檔案），
下面的程式片段都取自這些已實際執行驗證過的腳本。

```bash
python example/01_rtstruct_from_mask.py
python example/02_rigid_registration_to_reg.py
python example/03_deformable_registration_to_reg.py
python example/04_rtdose.py
```

### 從 3D 遮罩建立 RTSTRUCT，再讀回來

```python
import numpy as np
from pydicomrt.rs import RTStructBuilder, calc_image_series_affine_mapping, rtstruct_to_masks
from pydicomrt.utils import load_sorted_image_series

series = load_sorted_image_series("path/to/ct")     # 已依切片法線方向排序

mask = np.zeros((len(series), series[0].Rows, series[0].Columns), dtype=np.uint8)
mask[8:16, 30:60, 30:60] = 1                        # (slice, row, column)

rs_ds = (
    RTStructBuilder(series)
    .add_roi(mask=mask, name="CTV", color=[0, 255, 0])
    .build()
)

rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)

# 讀回來 —— affine 與形狀來自影像序列，不是來自 RTSTRUCT
affine_mapping, mask_shape = calc_image_series_affine_mapping(series)
masks = rtstruct_to_masks(rs_ds, affine_mapping, mask_shape)
recovered = np.asarray(masks["CTV"]["mask_volume"])
```

### 剛體註冊並匯出成 DICOM REG

注意 `.GetInverse()` —— 原因見上方慣例 2。

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
    SpatialRegistrationBuilder(fixed_series)           # 用 FIXED 序列建構
    .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    .add_registration(moving_series, reg_matrix.astype(np.float32).ravel().tolist())
    .build()
)
reg_ds.save_as("registration.dcm", enforce_file_format=True)
```

Builder 在寫入前會驗證矩陣 —— 16 個元素、最後一列為 `[0, 0, 0, 1]`、
以及型別為 `RIGID` 時上方 3×3 必須是真正的旋轉矩陣。
若確實需要縮放或剪切，請傳入 `matrix_type="RIGID_SCALE"` 或 `"AFFINE"`。
它同時會為 fixed Frame of Reference 附加一個 identity item —— 這是標示「其他矩陣以哪個座標系為基準」的方式；
若下游系統不接受，可用 `build(include_identity=False)` 關閉。

讀回來：

```python
from pydicom import dcmread
from pydicomrt.reg import get_spatial_registrations

reg = get_spatial_registrations(dcmread("registration.dcm"))
matrix = reg[fixed_frame_of_reference_uid][moving_frame_of_reference_uid]   # moving -> fixed
```

### 可變形註冊並匯出成 DICOM REG

Deformable REG 的方向是 fixed → moving。對 rigid 後接 deformable 的流程，
pre matrix 使用單位矩陣，grid 放 residual deformation，post matrix 放原方向的 rigid；
這裡不要對 rigid 取反矩陣。

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
    RTDoseBuilder(reference_ds)          # 病患／檢查資訊
    .set_dose_grid(dose_sitk_image)      # 幾何資訊 + 縮放後的像素資料
    .add_referenced_plan(plan_ds)        # DoseSummationType 為 PLAN 時這是 Type 1C
    .build()
)
dose_ds.save_as("rtdose.dcm", enforce_file_format=True)

dose_image = get_dose_image(dose_ds)     # 轉回 sitk，已套用 DoseGridScaling
```

---

## 註冊功能：哪些驗證過、哪些還沒

註冊品質最容易讓人失望，所以這裡直接寫實測狀態，而不是宣稱。

**已用真實的計劃 CT / CBCT 配對與 TPS 產生的註冊檔驗證**（一組骨盆案例、單執行緒）：

- 序列讀取與 SimpleITK 自家的 GDCM reader **完全一致** —— 幾何與體素值都相同。
- `get_spatial_registrations` 能正確讀取 TPS 的註冊檔，方向與文件所述相符。
- `rigid_registration` 與臨床註冊整體相差 7.5mm，而原本最差的上下方向（z）已在 1mm 以內。

與 TPS transform 的實測偏差，單位 mm，格式 `(x, y, z)`：

| `optimizer=` | 偏差 | 距離 | 耗時 |
|---|---|---|---|
| `"regular_step"`（預設） | `(+2.8, −7.0, −0.9)` | **7.5** | 約 9 分鐘 |
| `"gradient_descent"` | `(+1.4, −8.0, −36.7)` | 37.5 | 約 1 分鐘 |

`"gradient_descent"` 可能在估計出的第一步就衝過頭、落進平坦區，然後收斂視窗判定「成功」——
回傳一個**比自己起點還差**的結果，而且完全不會拋出任何錯誤。
`"regular_step"` 在梯度方向反轉時會縮小步長，因此不會跑掉。
預設採用慢但可靠的那個；快速路徑仍可選用，適合預覽或需要反覆重跑的流程。

**已知落差：**

- **剩下的約 7mm 幾乎全在前後方向（y），而且不一定是本函式庫的誤差。**
  在臨床答案附近掃描 y 可以看到，全視野互資訊的最佳點落在偏離臨床答案約 10mm 處 ——
  在這組資料上，指標與 TPS 本身就不一致。要再縮小就得把指標限制在 ROI 內，目前尚未實作。
  **一個案例不等於驗證**：請用你自己的容差標準檢查殘差。
- **結果只有在單執行緒下可重現。** ITK 依執行緒完成順序歸約指標值；
  在真實配對上，執行緒數會讓答案差異超過 15mm。需要確定性時請呼叫
  `sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)`。
- **只有軸向掃描做過端到端驗證。** 傾斜與非等向的幾何處理已實作並有單元測試，
  但尚未在真實的斜向檢查上跑過。

原本列在這裡的兩項落差都已解決，而且兩次都先被寫上了聽起來合理的解釋。
上下方向的誤差曾歸因於 CBCT 視野受限，真正的原因是 optimizer 跑離了正確位置；
pipeline 的問題曾歸因於視野不匹配，真正的原因是它在配準前把兩張影像重採樣到同一個網格，
讓置中初始化無事可做。下次再遇到「這組資料就是難」這種解釋時，值得先想起這件事。

---

## 模組結構

| 模組 | 內容 |
|---|---|
| `rs` | `builder`、`add_new_roi`、`make_contour_sequence`、`parser`、`check`、`iod`、`rs_to_volume`、`contour_process_method` |
| `reg` | `builder`、`parser`、`check`、`iod`、`type_transform` |
| `reg.method` | `rigid`、`demons`、`soft_demons`、`bspline`、`common` |
| `reg.pipeline` | `registration_pipeline`、`compose_transforms`、`preprocessing` |
| `dose` | `builder`、`check`、`iod`、`sitk_transform` |
| `ct` | `builder`、`check`、`iod` |
| `utils` | `image_series_loader`、`coordinate_transform`、`sitk_transform`、`validate_dcm_info`、`sitk_image_process` |

每個 package 都有匯出公開 API，所以 `from pydicomrt.rs import ...` 就夠了，
不需要深入到子模組。

完整手冊請見 [docs/user-guide_zh.md](docs/user-guide_zh.md)；
每個函式的參數請見 [docs/api-reference.md](docs/api-reference.md)；
資料流與擴充點請見 [docs/architecture.md](docs/architecture.md)。

UID 根節點：RTSTRUCT 與 RTDOSE 會讀取 `DICOM_UID_PREFIX` 環境變數；
註冊與 CT 的 builder 則使用 `set_uid_prefix()`。

---

## 開發

```bash
pip install -e ".[dev]"
pytest                  # 約 295 個測試，約 34 秒
pytest -m slow          # 真實資料的完整註冊，約 22 分鐘
```

測試會在記憶體中自行合成 DICOM 序列與體模，**不需要任何病人資料**。
建構函式在 `test/synthetic.py`，fixture 在 `test/conftest.py`。

### 用真實臨床資料測試

`test/real_data/` 是選用的測試層，會針對真實的計劃 CT、治療時的 CBCT，
以及 TPS 產生的註冊檔執行。資料不存在時會自動跳過。目錄結構如下：

```
<testdata>/CTSIM/*.dcm          計劃 CT 序列
<testdata>/CBCT/*.dcm           治療時的 CBCT 序列
<testdata>/registration.dcm     TPS 產生的 Spatial Registration
```

```bash
PYDICOMRT_TESTDATA=/path/to/testdata pytest test/real_data
PYDICOMRT_TESTDATA=/path/to/testdata pytest test/real_data -m slow
```

**請勿把資料放進版本庫** —— `testdata/` 已列入 .gitignore。

### 若你要新增測試

- 斷言請針對 **映射結果**（transform 把探測點送到哪裡），而非 transform 的原始參數。
  同一個映射有多種參數表示法，而且帶旋轉中心的 transform，其 `GetTranslation()` 不等於 offset。
- 修改幾何相關程式時，請加上**非等向間距**或**旋轉方位**的案例。
  方形像素的軸向 fixture 無論程式對錯都會通過。
- `conftest.py` 把 SimpleITK 固定為單執行緒，原因如上所述。

---

## 參考資源

- [SimpleITK](https://simpleitk.org/) · [pydicom](https://pydicom.github.io/)
- [RT-Utils](https://github.com/qurit/rt-utils) · [PlatiPy](https://github.com/pyplati/platipy)

## 授權

MIT —— 見 [LICENSE](LICENSE)。

## 作者

Higumalu (higuma.lu@gmail.com)
