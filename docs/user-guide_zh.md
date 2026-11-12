# pydicomRT 使用手冊

pydicomRT 0.9 的工作導向參考文件。[README_zh.md](../README_zh.md) 是導覽，這份是你之後會
反覆回來翻的手冊。想了解函式庫內部怎麼組起來的，請看
[architecture.md](architecture.md)（英文）。

**English: [user-guide.md](user-guide.md)** · **逐函式參考：
[api-reference.md](api-reference.md)**

**0.9 重新命名了大部分公開 API。** 從 0.8 升上來的話，請先看
[從 0.8 遷移](#從-08-遷移)。

- [安裝](#安裝)
- [四個慣例](#四個慣例)
- [載入影像](#載入影像)
- [RT Structure Set](#rt-structure-set)
- [註冊（Registration）](#註冊registration)
- [RT Dose](#rt-dose)
- [CT Image](#ct-image)
- [驗證](#驗證)
- [選擇與調校註冊方法](#選擇與調校註冊方法)
- [疑難排解](#疑難排解)
- [從 0.8 遷移](#從-08-遷移)

---

## 安裝

```bash
pip install pydicomrt
```

需要 Python ≥ 3.10 與 **pydicom ≥ 3.0**。pydicom 3 是硬性需求而非偏好：builder 從
`file_meta` 的 transfer syntax 取得編碼方式，而 pydicom 2 需要的
`is_little_endian` / `is_implicit_VR` 屬性已被 pydicom 4 移除。

GPU 加速的 demons 註冊是選用的額外套件。PyTorch 安裝動輒數 GB，所以預設不安裝，而且除非
指定 `backend="torch"`，否則任何地方都不會 import 它：

```bash
pip install "pydicomrt[gpu]"
```

若 pip 預設的 PyTorch wheel 在你的平台上只有 CPU 版，請先自行安裝 CUDA 版的 PyTorch（見
pytorch.org）；這個 extra 只是把 `torch` 加進需求。

wheel 內含一個小型的 C++ 編譯函式庫，是 `rigid_registration` 預設演算法的內層迴圈。它本身
沒有任何依賴，所以不是額外套件，每次安裝都會有。從原始碼安裝（`pip install -e .`，或沒有
現成 wheel 的平台）時會編譯它，需要 C++17 編譯器；CMake 會自動下載。沒有編譯器也能裝好，
`rigid_registration` 改用 numpy 跑同一個演算法，約慢十倍。

存檔一律使用 pydicom 3 的寫法：

```python
ds.save_as("out.dcm", enforce_file_format=True)   # write_like_original 已移除
```

---

## 四個慣例

這四件事弄錯**不會拋例外** —— 檔案讀得進去、影像看起來正常，只有物理座標是錯的。
README 有完整說明，這裡是精簡版。

**1. 軸順序會翻轉三次。** 由序列堆疊出來的 numpy 體資料是 `(slice, row, column)`；
`sitk.Image` 回報的是 `(x, y, z)`；DICOM 的 `PixelSpacing` 是 `[row 間距, column 間距]`，
也就是 `(y, x)`。傳給 `RTStructBuilder.add_roi()` 的遮罩是 `(slice, row, column)`。

**2. 註冊轉換的方向是 fixed → moving。** 函式庫所有註冊回傳的都是*重取樣*用的轉換，
也就是 `sitk.Resample` 需要的方向。DICOM Spatial REG 矩陣方向相反，匯出前要先取反；Deformable REG 則維持 fixed → moving。

**3. 先合成轉換，最後只重取樣一次。** `sitk.CompositeTransform([a, b])` 會先套用 `b`，
所以依階段順序排列的 list 不需要反轉。請用 `compose_transforms()`。

**4. 補值不是免費的。** 重取樣到來源影像之外會憑空造出體素，CT 的正確值是 `-1000`。
MR、PET 或做過 window clip 的影像請明確傳入 `default_value` —— 自動推論是讀強度最小值，
而被 clip 過的 CT 已經不像 CT 了。

---

## 載入影像

```python
from pydicomrt.utils import load_sorted_image_series, image_series_to_sitk_image

series = load_sorted_image_series("/path/to/CT")        # 資料夾，或路徑 list
image = image_series_to_sitk_image(series)              # -> sitk.Image，已套用 HU
```

`load_sorted_image_series` 沿著切片法向量排序。不要依賴檔名順序或 `InstanceNumber`：
未排序的序列一樣會產生看起來正常的體資料，只是 z 軸翻轉，而建立在它之上的每一個輪廓與
註冊都會繼承這個翻轉。

轉換時會套用 rescale slope 與 intercept，所以 `sitk.Image` 是真實單位。

需要分開拿各部分，或想把自己的體資料放到某序列的幾何上：

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

# 或直接沿用既有影像的幾何
image = SimpleITKImageBuilder().from_reference_image(my_array, reference_image)
```

---

## RT Structure Set

### 遮罩 → RTSTRUCT

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

沒有傳 `number=` 時，ROI 編號會依序配發。遮罩形狀會與序列核對，因為形狀不符的遮罩會被切到
錯誤的影像上，產出位置錯誤的輪廓，而不是報錯。

`add_roi()` 另外接受 `description` 與 `contour_config`（輪廓萃取的調校參數：雜訊尺寸與
低通比例，見 `rs.builder.DEFAULT_CONTOUR_CONFIG`）。

### RTSTRUCT → 遮罩

```python
from pydicomrt.rs import calc_image_series_affine_mapping, rtstruct_to_masks

affine, shape = calc_image_series_affine_mapping(series)   # 來自影像，不是來自 RTSTRUCT
masks = rtstruct_to_masks(rs_ds, affine, shape)

ctv = np.asarray(masks["CTV"]["mask_volume"])              # (slice, row, column)
```

affine 與體資料形狀都來自影像序列。RTSTRUCT 存的是病患座標系下的輪廓點，它並不知道你想把
這些點放到哪個網格上。

### 讀取結構的中繼資料

```python
from pydicomrt.rs import get_roi_names, get_contours, is_rtstruct_matching_series

get_roi_names(rs_ds)                    # {roi_number: roi_name}
get_contours(rs_ds)                     # 依 ROI 與參考影像分組的輪廓
is_rtstruct_matching_series(rs_ds, series)
```

`get_roi_names` 的 key 是 `pydicom.valuerep.IS`（int 的子類別）。`names[1]` 可以，
`names["1"]` 會拋 `KeyError`，儘管 repr 顯示成 `{'1': ...}`。

### 已經有輪廓點的情況

如果你手上已經是病患座標系的輪廓點而不是遮罩：

```python
RTStructBuilder(series).add_roi_from_contours(
    contours={sop_instance_uid: {"sop_class_uid": ..., "contours": [[x, y, z, ...], ...]}},
    name="CTV",
).build()
```

---

## 註冊（Registration）

### 剛體註冊並匯出成 DICOM REG

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

# DICOM 要的是 moving -> fixed，所以取反
matrix = affine_to_homogeneous_matrix(transform.GetInverse())
reg_ds = (
    SpatialRegistrationBuilder(fixed_series)                   # 用 FIXED 序列建立
    .set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    .add_registration(moving_series, matrix.ravel().tolist())
    .build()
)
reg_ds.save_as("registration.dcm", enforce_file_format=True)
```

builder 在寫入前會驗證矩陣：16 個元素、最後一列是 `[0, 0, 0, 1]`、數值有限，以及當型別是
`RIGID` 時上 3×3 必須是真正的旋轉。`RIGID_SCALE` 也會檢查，並且會擋掉剪切。如果本來就想要
剪切，請明確宣告 `matrix_type="AFFINE"`，不要讓它掛著較窄的型別名義混過去。

`build()` 會為 fixed 的 Frame of Reference 補上一個單位矩陣項目 —— 這正是標示其他矩陣以哪
個 frame 為基準的東西。若下游系統不接受，可用 `build(include_identity=False)`。接著它會跑該
物件自己的 IOD checker，不符合就丟出 `Invalid Spatial Registration: ...`（或
`Invalid Deformable Registration: ...`），而不是交回一個過不了 conformance 的 dataset。
建構子這邊則要求所有 reference 影像共用同一個 Frame of Reference，所以由兩個 study 拼起來的
序列 list 在還沒加入任何 registration 前就會被擋下。

兩種 REG 物件都帶有四個 Type 1 的設備欄位 —— Manufacturer、Manufacturer's Model Name、
Device Serial Number 與 Software Versions —— 只要有一個缺漏或為空，TPS 就會拒絕匯入。它們的
預設值是 `pydicomRT`（Software Versions 另外記錄安裝的版本）。若接收端系統要求填院內自己的
值，可用 `set_equipment()` 設定；沒傳入的參數維持預設值：

```python
builder.set_equipment(manufacturer="My Hospital", device_serial_number="SN-001")
```

### 可變形註冊並匯出成 DICOM REG

Deformable REG 的方向是 fixed → moving。對 rigid 後接 deformable 的流程，grid 放 residual
deformation，post matrix 放原方向的 rigid，不要取反矩陣。`pre_transform` 與 `post_transform`
的預設值都是 `None`，會直接省略該 sequence —— 不需要傳單位矩陣進去。

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

`demons_registration` 與 `bspline_registration` 都回傳
`(registered_image, transform, deformation_field)`，可以互換。

### 讀取既有的 REG

```python
from pydicom import dcmread
from pydicomrt.reg import get_spatial_registrations, get_deformable_registrations

reg = get_spatial_registrations(dcmread("registration.dcm"))
matrix = reg[fixed_frame_uid][moving_frame_uid]        # moving -> fixed，16 個浮點數

deformable = get_deformable_registrations(dcmread("deformable.dcm"))
field = deformable[0]["DeformableRegistrationGrid"]["VectorGridData"]   # (z, y, x, 3)
```

**以下只適用於 Spatial REG 的矩陣。** 要把它當作重取樣轉換來*套用*，必須先取反 ——
Spatial REG 存的是 moving → fixed，而 `sitk.Resample` 要的是 fixed → moving：

```python
matrix_4x4 = np.array(matrix, dtype=float).reshape(4, 4)
inverse = np.linalg.inv(matrix_4x4)
transform = sitk.AffineTransform(3)
transform.SetMatrix(inverse[:3, :3].ravel().tolist())
transform.SetTranslation(inverse[:3, 3].tolist())
```

Deformable 物件的 pre/post 矩陣本身就是 fixed → moving，直接使用即可 —— 對它們取反正是
這個方向區分要防的錯誤。缺少的 pre/post sequence 會以單位矩陣回傳，而且只會回傳帶 grid
的項目。

### Pipeline

```python
from pydicomrt.reg.pipeline import registration_pipeline

registered, rigid, deformable, field = registration_pipeline(
    fixed_image, moving_image,
    perform_rigid=True,
    perform_deformable=True,
    preprocess_config={"rigid": {"window_clip": [-200, 800]}},
)
```

回傳的影像一定在 `fixed_image` 的網格上、帶有 moving 影像的像素型別，而且是用合成後的轉換
對**未經修改的** moving 影像做單次重取樣 —— 所以前處理不會洩漏到輸出，影像也與它回傳的轉換
完全一致。

沒有設定任何前處理時，pipeline 的行為與直接呼叫 `rigid_registration` 完全相同。

---

## RT Dose

```python
from pydicomrt.dose import RTDoseBuilder, get_dose_image

dose_ds = (
    RTDoseBuilder(reference_ds)              # 一張定位 CT 切片，或計畫檔
    .set_dose_grid(dose_sitk_image)          # 單位為 Gy
    .add_referenced_plan(plan_ds)
    .build()
)
dose_ds.save_as("rtdose.dcm", enforce_file_format=True)

dose_image = get_dose_image(dose_ds)         # 轉回 sitk，已套用 DoseGridScaling
```

`DoseSummationType` 預設是 `"PLAN"`，這會讓 Referenced RT Plan Sequence 成為 Type 1C —— 也
就是必須提供計畫參考。缺少時 `check_rtdose_iod()` 會回報。`set_dose_grid()` 另外接受
`dose_type`（預設 `"PHYSICAL"`）、`dose_summation_type`、`dose_units` 與 `dose_grid_scaling`。

劑量是以縮放後的無號 32 位元整數儲存。預設縮放係數 1e-7 可涵蓋約 429 Gy，超過時 builder 會
自動提高縮放係數而不是溢位繞回。負值劑量會被拒絕，因為無號儲存會把它變成極大的正值。

---

## CT Image

```python
from pydicomrt.ct import CTBuilder

slices = (
    CTBuilder(reference_series)              # 病患／檢查資訊來源
    .set_volume(sitk_image)                  # 數值視為 HU
    .set_plane("AXIAL")                      # 或 CORONAL、SAGITTAL
    .build()
)
for index, slice_ds in enumerate(slices):
    slice_ds.save_as(f"ct_{index:04d}.dcm", enforce_file_format=True)
```

平面名稱指的是**每張切片所在的平面**，所以切片堆疊的方向是名稱裡*沒有*出現的那個軸：
axial 沿 z 堆疊、coronal 沿 y、sagittal 沿 x。

輸出會標記為 `DERIVED\SECONDARY` 並取得自己的 Series Instance UID，因此不會被誤認為原始
掃描序列。數值會被裁切到無號 16 位元搭配 -1024 intercept 所能表示的範圍，即 -1024 至
64511 HU。

`set_volume_from_array(volume, origin, spacing, direction)` 改為接受 numpy 體資料加幾何資訊。
`set_copy_all_attributes(True)` 會從參考序列繼承所有屬性，而不是建立全新的 dataset —— 對保留
掃描參數很方便，但來源帶的其他東西也會一起繼承。

---

## 驗證

每個模態都有形狀一致的檢查器：

```python
from pydicomrt.rs import check_rtstruct_iod
from pydicomrt.reg import check_spatial_reg_iod, check_deformable_reg_iod
from pydicomrt.dose import check_rtdose_iod
from pydicomrt.ct import check_ct_iod

result = check_rtstruct_iod(rs_ds)
# {"result": True, "content": []}          -- 符合規範
# {"result": False, "content": ["Missing InstanceNumber in root", ...]}
```

它們依物件的 CIOD 檢查必填欄位，`check_rtdose_iod` 另外檢查大多數 Dose Summation Type 對
Referenced RT Plan Sequence 施加的條件式要求。

兩個 REG checker 比單純的欄位存在檢查再多做一些：驗證 SOP Class、轉換矩陣是否符合它自己宣告
的型別（正交性），以及 deformation grid 內部是否自洽 —— 維度、間距、正交的方向矩陣，還有
`VectorGridData` 的位元組數要對得上 `GridDimensions`。它們都看不出來的是註冊的**方向**對不對：
一個方向反了的轉換在結構上是完美的。

CT 是 single-frame IOD，所以 `check_ct_iod` 一次檢查一張切片。

---

## 選擇與調校註冊方法

### 選哪一種

| | 適用時機 | 成本 |
|---|---|---|
| `rigid_registration` | 只需要位置與方向 | 腦部 MR/CT 8 執行緒約 1 秒 |
| `demons_registration` | 可變形，預設首選 —— 同模態（CT↔CT、CT↔CBCT） | 快，且對參數寬容 |
| `demons_registration(metric="mind")` | 可變形，不同模態（MR→CT） | GPU 上數秒；需要 `[gpu]` |
| `bspline_registration` | 可變形，需要參數化轉換時 | 同等精度下約為 demons 的 40 倍 |
| `demons_with_soft_mask` | 可變形，限制在特定結構內 | 同 demons |

**可變形註冊請優先用 demons。** 在還原 2.11 mm 形變的合成假體上，demons 為 0.80 mm／0.4 秒，
B-spline 為 0.82 mm／17 秒。當你確實需要參數化轉換 —— 可儲存、可操作的控制點 —— 而不是位移
場時，才選 B-spline。demons 也可以在 GPU 上跑，見下文。

**demons 比較的是亮度**，所以預設只在同一組織在兩張影像亮度相同時才成立。MR 對 CT 請用
`metric="mind"`，見 [MR 對 CT](#mr-對-ctmetricmind)。

### 剛體

```python
transform = rigid_registration(fixed_image, moving_image)
```

本套件自己實作的 elastix 策略：每次迭代在新的隨機點上計算 Mattes 互資訊，在四層平滑金字塔上
做自適應隨機梯度上升，最後在一組固定的 50 萬個點上做確定性的收尾。它用的是線性內插本身的精確
導數，排程是固定的。

與 elastix 預設的 rigid 比較（平均誤差，mm；超過 5 mm 算失敗）：

| | `rigid_registration` | elastix |
|---|---:|---:|
| RIRE 腦部 MR → CT，骨釘真值（6 個 MR 序列） | 1.38 | 1.32 |
| RIRE，同一次檢查 T1/T2/PD 的一致性（16 位病人） | 1.33 | 1.31 |
| RIRE，MR header 隨機移動最多 ±10°、±20 mm（10 個起點） | 1.24 | 1.35 |
| SynthRAD2023 腦部 MR → CT，header 移動（9 組） | 0.53 | 1.00 |
| 骨盆 CT/CBCT，對 TPS 註冊（1 組） | 5.1 | 16.7 |
| 時間，腦部 MR/CT，8 執行緒 | 約 1.2 秒 | 約 1.2 秒 |

兩者都在同一位 RIRE 病人上失敗（該次檢查的序列之間相差約 8 mm），比較可能是序列之間病人真的
有移動。骨盆資料上，單執行緒約 15 秒。

- `backend`：`"auto"`（預設）使用編譯核心；若是在沒有 C++ 編譯器的環境從原始碼安裝，則用
  numpy，約慢十倍。兩者抽到相同的隨機點；但捨入誤差偶爾會讓個別點落到相鄰的直方圖 bin，
  隨機步驟會把它延續下去，所以兩個後端的答案最多相差約 0.2 mm。
- `seed`：隨機點的種子。骨盆資料上種子 0–2 的結果是 5.09–5.16 mm。
- `threads`：編譯核心的執行緒數，`None` 時取 SimpleITK 的全域預設值。答案與它無關。

`algorithm="itk"` 保留先前的預設：SimpleITK 的 `ImageRegistrationMethod`，以及它的參數
（`histogram_bins`、`optimizer`、`shrink_factors` 等）。它只有在單執行緒下才逐位元可重現，
因為 ITK 依執行緒完成順序歸約 metric。

### 在 GPU 上跑 demons

`backend="torch"` 用 PyTorch 跑同一個多解析度 demons，有 GPU 時就用 GPU。需要安裝 `[gpu]`
extra（見[安裝](#安裝)），回傳的三個物件和原本相同：

```python
registered, transform, field = demons_registration(
    fixed_image, moving_rigid, backend="torch"       # device="cuda:1" 可指定顯卡
)

# 或透過 pipeline
registration_pipeline(fixed_image, moving_image,
                      deformable_kwargs={"backend": "torch"})
```

它是 ITK filter 的移植，而不是照概念重寫 —— 邊界處理、平滑 kernel 和提前停止都依照 ITK
原始碼 —— 測試會在非等向、旋轉、moving 在不同格點等情況下和 SimpleITK backend 比對，兩者
差距在千分之一毫米等級。在一組真實的 418 × 417 × 126 定位 CT（2200 萬 voxel；RTX 5090、
24 核 CPU）上，使用預設參數：

| backend | 時間 |
|---|---|
| `"simpleitk"`，單執行緒（`ncores` 預設值） | 約 63 秒 |
| `"simpleitk"`，`ncores=24` | 約 8 秒 |
| `"torch"`，CUDA | 約 0.7 秒 |

兩者的時間都隨迭代次數線性增加，所以迭代調高後倍數不變。

- **記憶體**取決於金字塔最細一層，每個 voxel 約 180 bytes —— 那組 2200 萬 voxel 的 CT 約
  4 GB。1024 × 1024 × 200 的影像在全解析度下需要約 39 GB；改用
  `isotropic_resample=True, resolution_staging=(8, 4, 1.5)` 後只用 5.2 GB、3.7 秒（24 執行緒
  要 36.7 秒），準確度相同，回傳的形變場仍在 fixed 格點上。記憶體不足時會丟出附帶預估值的
  `torch.OutOfMemoryError`。
- **`backend="simpleitk"` 的 RAM** 一樣取決於最細一層，每個 voxel 約 150 bytes —— 512 × 512 ×
  160 的 CT 約 6 GB，大部分花在 ITK 的 demons filter 裡。RAM 不夠的電腦用同樣的方式把最細一層
  變粗：改成 2 mm（`isotropic_resample=True, resolution_staging=(8, 4, 2)`）後同一組 CT 只要
  2.5 GB，其中 2 GB 是回傳的 transform 和形變場，兩者都是全解析度的 `float64`。
- **多出來的迭代要加在最細一層。** 在真實 CT 加上已知合成形變的測試裡，每一層都加會變差 ——
  預設 `(10, 10, 10)` 是 1.15 mm，`(100, 100, 100)` 是 4.45 mm —— 只加最後一層則變好：
  `(10, 10, 100)` 是 0.55 mm。兩個 backend 都是如此，GPU 只是讓長排程變得負擔得起。
- **影像邊緣在兩個 backend 上都不穩定。** 解剖被影像邊界截斷的地方，demons 是混沌的：輸入
  只要擾動千分之一 HU，SimpleITK 在那裡的結果就會移動數毫米。兩個 backend 在那裡的差異也是
  這個量級；影像內部則一致。
- 同一個裝置上重複執行，結果完全相同。torch backend 只支援 3D、只用線性內插，並忽略
  `ncores` 和 `verbose`。沒有 CUDA 裝置時會在 CPU 上執行並記錄警告。

### MR 對 CT：`metric="mind"`

demons 的力是兩張影像的亮度差，這在 MR 和 CT 之間沒有意義 —— 骨頭在一邊亮、在另一邊暗。在
沒有形變的腦部 MR/CT 配對上，預設 metric 會把腦平均移動 9 mm；`registration_pipeline` 的
可變形階段就是 demons，所以也一樣。`metric="mind"` 把亮度換成 MIND-SSC 描述子（Heinrich
et al., 2013）：每個體素周圍的小區塊彼此有多像。這取決於邊在哪裡，而不是哪一側比較亮，
所以只要 MR 和 CT 都看得到同一條邊界，描述子就相同。

```python
rigid = rigid_registration(ct_image, mr_image)            # 互資訊：MR 沒問題
mr_on_ct = sitk.Resample(mr_image, ct_image, rigid, sitk.sitkLinear, 0.0)  # MR 用 0 補邊

registered, transform, field = demons_registration(
    ct_image, mr_on_ct, backend="torch", metric="mind", default_value=0.0,
)

# 或透過 pipeline
registration_pipeline(ct_image, mr_image, default_value=0.0,
                      deformable_kwargs={"backend": "torch", "metric": "mind"})
```

`example/05_mr_to_ct_deformable.py` 用一組已知答案的合成配對完整跑一次，並寫出 Deformable
REG；匯出方式和 CT 對 CT 相同。

- **需要 `backend="torch"`**（`[gpu]` 這個 extra；CUDA 顯卡不是必要的，CPU 也能跑，只是慢）。
- **直接輸入原始數值。** MIND 不受增益與偏移影響，CT 值和 MR 訊號照原樣傳入即可。MR 配對
  不要在 pipeline 設 `window_clip`：它會把同一個窗位套到兩張影像上，CT 的窗位會把 MR 壓平。
- **預設值不同**：`resolution_staging=(4, 2, 1)`、`iteration_staging=(100, 100, 200)`、
  `regularization_kernel_mm=3.0`、`smoothing_sigma_factor=0.5`。intensity 的預設值
  （(8, 4, 1)、1.5 mm）在這個 metric 下各自會讓憑空產生的位移升到 2–4 mm。明確傳入的值則照原樣使用。
- **記憶體**：最細一層每個體素約 360 bytes，是 intensity metric 的兩倍。原始解析度的定位
  CT 請在 1–1.5 mm 下註冊：`isotropic_resample=True, resolution_staging=(6, 3, 1.5)`。回傳的
  場仍在 fixed 的格點上。
- 它不會重現 SimpleITK 的結果：它讓同一個場貫穿所有金字塔層級，而不是每層組合一個修正量，
  這樣較細的層級才能修正較粗層級的錯。

在九組公開的腦部 MR/CT 配對上量測（SynthRAD2023：放療定位 CT 與 3D T1，主辦方只做了剛體
對位），已知形變加在 CT 上：

| | 未形變配對上腦的移動量 | 已知 2.7 mm 形變的誤差 | 時間 |
|---|---|---|---|
| 不做可變形 | 0 | 2.74 mm | — |
| `metric="intensity"` | 9.1 mm | 9.3 mm | 約 1 秒 |
| `bspline_registration(metric="mutual_information")`（單例） | 6.6 mm | 6.6 mm | 13 分鐘 |
| `metric="mind"` | 0.62 mm | 0.86 mm | 約 9 秒 |

互資訊適用於剛體註冊 —— 只有六個參數 —— 但密集或 B-spline 形變可以靠扭曲腦外的結構把它
拉高，而它真的會這麼做。未形變配對上的 0.62 mm 有一部分是真的（MR 失真、姿勢；這些配對只做
過剛體對位）：在直接由 CT 做出、完全沒有殘差的合成 MR 上，`"mind"` 讓腦移動約 0.1 mm。

**腦部標靶通常只要剛體就夠了。** 在這些量測中，可變形階段讓腦移動的量遠小於 1 mm；真正的
形變在頸部和顱底，8° 的頸部屈曲只還原了約三分之一。而且可變形結果會把真實的解剖變化 ——
腫瘤長大、水腫、術後腔室 —— 一起「對齊」掉。在上面畫輪廓之前請先檢查。MIND 看的是邊，所以
層厚 3–5 mm 的 2D MR 序列在層面方向上提供的資訊很少；請用 3D T1 對 CT 註冊，同一次 MR 檢查的
其他序列再透過它們與該 T1 的剛體對位帶過來。

### B-spline

決定一次 B-spline 擬合有沒有意義的參數，是 **metric 取樣點數與轉換係數量的比值**。取樣太少
時，optimizer 會靠扭曲網格把 metric 壓到近乎 0、宣告收斂，然後回傳一個在 metric 沒看過的地
方完全不受約束的形變。過程中不會有任何例外。

`describe_sampling()` 會逐層回報這個比值，而 `bspline_registration` 在某層欠定時會發出警告。

**取樣率不是速度旋鈕。** 它看起來很像：

| 取樣 | 誤差 | 時間 |
|---|---|---|
| `None`（稠密，預設） | 0.89 mm | 37 秒 |
| `REGULAR` 100% | 1.50 mm | 47 秒 |
| `REGULAR` 10% | 2.04 mm | 4 秒 |
| `RANDOM` 1% | 1.64 mm | 1 秒 |

`REGULAR` 100% 取樣的體素與 `None` 完全相同，結果仍是 1.50 mm，所以差距來自「一旦進入取樣
點集合模式，ITK 走的導數計算路徑不同」，而不是取樣數量。無論比率多少，所有取樣設定都落在
1.5 至 2.0 mm 之間。保持 `None` 就好。

要拿精度換時間，請用那些**平滑變化**的旋鈕：

| 旋鈕 | 效果 |
|---|---|
| `ncores` | 近乎線性；單執行緒 111 秒 → 十六執行緒 19 秒 |
| `number_of_iterations` | 各上限的平均誤差：5 → 1.09、10 → 1.01、**20 → 0.82**、30 → 0.88、50 → 1.07 mm |
| `initial_grid_spacing` | 越粗越快；96 mm → 1.00 mm／11 秒，64 mm → 0.89 mm／37 秒 |
| `resolution_staging` | 層數越少越快；`(2, 1)` → 0.95 mm／11 秒 |

注意迭代上限有**最佳值**而非天花板：超過約 20 次之後，擬合開始追逐雜訊，調高它同時損失時間
*與*精度。

`optimizer="LBFGSB"` 需要搭配 `grid_scale_factors=None` —— ITK 只會設定一次 optimizer scale，
而縮放的網格會讓各層之間的係數數量改變。

### 哪些驗證過、哪些沒有

註冊精度是在**一組**真實定位 CT／CBCT 骨盆資料上與 TPS 註冊比對，以及在合成假體上與已知
形變比對得出的。

- 序列載入與 SimpleITK 自己的 GDCM reader 完全一致，包含幾何與體素值。
- `rigid_registration` 整體重現臨床註冊到 5.1 mm，`(+1.6, −4.7, −1.1)` mm。殘差主要在前後軸，
  該處全視野互資訊的峰值偏離臨床答案 —— 也就是說在這組資料上，metric 與 TPS 是真的不一致。
  另外在公開的腦部 MR/CT 上與骨釘真值比對過（RIRE：一位病人，1990 年代 4 mm 切片）。
- **可變形註冊從未與臨床標準答案比對過。** 上面那些 B-spline 數字是用 B-spline 還原
  B-spline 形變，是最有利的測試條件。真實的解剖形變並不是 B-spline 形狀的。
- `metric="mind"` 是在九組公開的腦部 MR/CT 配對上，以加在真實影像上的*已知合成*形變量測
  的 —— 只有 3D T1，沒有 T2/FLAIR，也沒有腫瘤變化。

一組資料不構成驗證。請用你自己的容許值去檢查殘差。

---

## 疑難排解

### 載入

| 訊息 | 原因 |
|---|---|
| `No DICOM files found in ...` | 該目錄沒有 `*.dcm`，或你傳了檔案路徑而非目錄 |
| `Error reading DICOM file ...` | 不是 DICOM，或檔案截斷 |
| `ImageOrientationPatient [...] is not two orthogonal unit vectors` | orientation 損毀或不符規範。以 `InvalidImageOrientationError`（`ValueError` 子類別）拋出 |

### Builder

| 訊息 | 處理方式 |
|---|---|
| `no ROI added; call add_roi() before build()` | builder 沒有任何內容可寫 |
| `no registration added; call add_registration() before build()` | 同上，REG |
| `no dose grid set; call set_dose_grid() before build()` | 同上，劑量 |
| `no volume set; call set_volume() before build()` | 同上，CT |
| `mask shape (...) does not match the image series (...)` | 遮罩必須是 `(slice, row, column)` 且與序列相符 |
| `roi_number N was already added` | 兩個 ROI 用了同一個編號 |
| `the bottom row of the matrix must be [0, 0, 0, 1]` | 矩陣是 row-major，平移放在最後一**行** |
| `matrix_type is RIGID but the upper 3x3 is not orthonormal` | 改宣告 `RIGID_SCALE` 或 `AFFINE` |
| `RIGID_SCALE requires orthogonal, non-zero axes; use AFFINE for shear` | 矩陣含剪切，`RIGID_SCALE` 不涵蓋 |
| `Invalid Spatial Registration: ...` / `Invalid Deformable Registration: ...` | `build()` 自己的 conformance 檢查沒過，訊息會列出缺什麼 |
| `reference images must share one non-empty FrameOfReferenceUID` | 傳給 REG builder 的序列混到了不同的 Frame of Reference |
| `DICOM deformation grids require a 3D field with three components` | 位移場不是 3D、每點三分量 |
| `DICOM deformation grid direction must be right-handed and orthonormal` | 位移場的方向矩陣是鏡射或非正交 |
| `reference_ds is not an RT Plan` | `add_referenced_plan()` 需要 RT Plan（SOP class `...481.5`） |
| `dose grid holds negative values` | 儲存的劑量是無號的，請先裁切或位移 |

### 註冊

| 症狀 | 可能原因 |
|---|---|
| 警告 "the fit is underdetermined" | B-spline 的係數多於 metric 取樣點。調高 `sampling_rate`、加粗 `initial_grid_spacing`，或去掉最粗那層 |
| `optimizer must be one of (...)` | 檢查拼字，清單在 `bspline.OPTIMIZERS` |
| `optimizer='LBFGSB' cannot be combined with grid_scale_factors` | 傳入 `grid_scale_factors=None` |
| 每次跑結果都不同 | 用 `algorithm="itk"` 時，ITK 依執行緒完成順序歸約 metric：固定為單執行緒。預設演算法與執行緒數無關；檢查 `seed` |
| 註冊結果比什麼都不做還糟 | 若用了 `algorithm="itk", optimizer="gradient_descent"`，改回預設值 |
| `ImportError: backend='torch' needs PyTorch` | `pip install "pydicomrt[gpu]"` |
| `ImportError: backend='cpp' needs the compiled core` | pydicomrt 是在沒有 C++ 編譯器的環境從原始碼安裝的。改裝 wheel，或在有 C++17 編譯器的環境重裝；在那之前 `backend="auto"` 會用 numpy |
| `backend="torch"` 丟出 `torch.OutOfMemoryError` | 最細一層對這張顯卡太大。用 `isotropic_resample=True, resolution_staging=(8, 4, 1.5)` 把它變粗（`metric="mind"` 用 `(6, 3, 1.5)`） |
| `demons_registration` 時出現 `MemoryError` 或程序被系統砍掉 | 最細一層放不進 RAM（每個 voxel 約 150 bytes）。用 `isotropic_resample=True, resolution_staging=(8, 4, 2)` 把它變粗；回傳的形變場仍在 fixed 格點上 |
| `metric='mind' runs on backend='torch' only` | 加上 `backend="torch"`；MIND 沒有 SimpleITK 的實作 |
| MR 對到 CT 後變形扭曲，不該動的地方動了數毫米 | 在多模態配對上用了預設的 intensity metric。改用 `metric="mind"` |

### 輸出看起來位移或變形

幾乎一定是四個慣例之一：

- **位移方向相反** —— 轉換用錯方向。註冊回傳 fixed → moving；DICOM Spatial REG 存的是 moving → fixed；Deformable REG 則是 fixed → moving。
- **平面內被拉長或壓扁** —— 某處把 row/column 間距交換了。`PixelSpacing` 是 `[row, column]`，
  SimpleITK spacing 是 `(x, y, z)`。
- **z 軸鏡像** —— 序列沒有沿切片法向量排序。
- **病人外圍多一層軟組織殼** —— 補值用了 0 而不是 -1000。

---

## 從 0.8 遷移

0.9 為了一致性重新命名了公開 API，且**不保留任何別名**。行為也有變動，就算你的 import 還能
解析，也請讀[行為上的變動](#行為上的變動)。

### 物件建構

五個 builder 現在共用同一套形狀：`Builder(參考來源).set_*(...).add_*(...).build()`。

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

0.8 那些函式仍然存在於各自的子模組中、也仍然可用；它們正是 builder 組裝所用的積木，只是不再
是文件建議的進入點。

### 更名

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

### 搬移與移除的模組

| 0.8 | 0.9 |
|---|---|
| `rs.checker` | `rs.check` |
| `rs.rs_ds_iod`（`RT_STRUCTURE_SET_IOD`） | `rs.iod`（`RTSTRUCT_IOD`） |
| `reg.s_reg_ds_iod`（`SPATIAL_REGSITRATION_IOD`） | `reg.iod`（`SPATIAL_REGISTRATION_IOD`） |
| `reg.ds_reg_ds_iod` | `reg.iod`（`DEFORMABLE_SPATIAL_REGISTRATION_IOD`） |
| `dose.dose_ds_iod`（`RT_DOSE_IOD`） | `dose.iod`（`RTDOSE_IOD`） |
| `ct.ct_ds_iod` | `ct.iod` |
| `rs.packer` | 已移除 —— 是 `make_contour_sequence` + `add_new_roi` 未使用的舊副本 |
| `utils.rs_from_altas` | 已移除 —— 該檔案是空的 |

序列參數一律命名為 `image_series`，或在角色重要時用 `fixed_series` / `moving_series`。
`optimiser` 統一拼成 `optimizer`。

### 行為上的變動

你的 import 可能照常解析，但輸出不一樣了。

**幾何修正。** 這些原本是錯的、現在是對的，也就是說受影響的資料輸出會改變。方形像素的軸向
CT 不受影響。

- `PixelSpacing` 曾被直接放進 SimpleITK 的 `(x, y)`，對非等向像素而言等於交換了平面內間距。
  由 1.5 × 3.0 mm 體資料寫出的序列，讀回來會變成 3.0 × 1.5。
- 方向矩陣曾把軸向量放在**列**而非**行**，導致斜切與旋轉過的序列落在鏡像位置。
- `ImageOrientationPatient` 曾取自 row-major 攤平後的前六個元素，而不是前兩**行**。

**像素編碼。** `CTBuilder` 與劑量 builder 寫入無號資料卻宣告 `PixelRepresentation = 1`
（有號），因此超過有號中點的數值讀回來會變成負數 —— CT 上的高密度植入物、以及約 215 Gy 以上
的劑量。兩者現在都宣告為無號。

**註冊。** `rigid_registration` 預設改用 regular-step optimizer；先前無界的梯度下降可能回傳
比自己初始化還差的結果。在真實 CT/CBCT 上這讓誤差從 37.5 mm 降到 7.5 mm。
`registration_pipeline` 不再於剛體階段前把兩張影像重取樣到共用的 union-extent 網格 —— 那讓
置中初始化形同無效：同一組資料從 194 mm 降到 7.5 mm。`bspline_registration` 現在回傳三個值
而非兩個，與 demons 一致。

**Deformable REG 方向 —— 已經寫出去的檔案請重新檢查。** Deformable Spatial Registration 的
方向是 fixed → moving（PS3.3 C.20.3.1.1），與註冊函式回傳的方向*相同*，和 Spatial REG 矩陣相反。
先前的版本把它當成 moving → fixed 來記錄與匯出，還把 rigid 階段取反塞進 pre matrix。因此
0.9.0 之前寫出的每一個 deformable 物件方向都是反的，需要重新匯出；正確做法是 grid 放 residual
field、`post_transform` 放原方向的 rigid 矩陣。Spatial REG 不受影響 —— 它確實是 moving → fixed。

同批還有兩個相關修正。`affine_to_homogeneous_matrix` 現在會把非零的旋轉中心折進矩陣位移項，
而不是靜默丟棄，所以帶中心的 rigid transform 不必再先過 `to_centre_free_affine`。
`get_spatial_registrations` 則會依 DICOM 順序合成 Matrix Sequence 裡的每一項，先前只讀最後一項。

**DICOM 規範符合度。** REG builder 先前缺少 `InstanceNumber`、`ContentDate` 與 `ContentTime`
（皆為 Type 1）、把 Registration Type Code Sequence 寫在錯誤的巢狀層級，也遺漏了 fixed Frame
of Reference 的單位矩陣項目。0.8 寫出的檔案讀得進去，但在這幾點上不符規範。兩個 REG builder
現在都會在 `build()` 裡跑自己的 checker，畸形的物件會在建構當下就失敗，而不是流到下游。

**pydicom 下限提高到 3.0。** 見[安裝](#安裝)。
