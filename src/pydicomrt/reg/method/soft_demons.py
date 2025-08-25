import SimpleITK as sitk
import numpy as np
from typing import Optional


def _infer_background_value(img: sitk.Image, default_fallback: float = 0.0) -> float:
    """If the image looks like a CT (min <= -1000), return -1000. Otherwise, return fallback."""
    try:
        arr_min = float(np.min(sitk.GetArrayViewFromImage(img)))
        return -1000.0 if arr_min <= -1000.0 else float(default_fallback)
    except Exception:
        return float(default_fallback)

def _ensure_binary_mask(mask: sitk.Image) -> sitk.Image:
    """Ensure binary mask (0/1)."""
    return sitk.Cast(mask > 0, sitk.sitkUInt8)

def _soft_weight_from_mask(mask: sitk.Image, falloff_mm: float = 8.0) -> sitk.Image:
    """
    Generate soft weights (0..1):
    - Inside ROI, far from boundary → ~1
    - Near boundary → 0
    - Outside ROI → 0

    Uses SignedMaurerDistanceMap with linear falloff in mm.
    """
    mask = _ensure_binary_mask(mask)
    dt = sitk.SignedMaurerDistanceMap(
        mask, squaredDistance=False, insideIsPositive=True, useImageSpacing=True
    )
    w = sitk.Clamp(dt / float(falloff_mm), lowerBound=0.0, upperBound=1.0)
    return sitk.Cast(w, sitk.sitkFloat32)

def _apply_soft_mask(img: sitk.Image, w: sitk.Image, background_value: float) -> sitk.Image:
    """Apply soft mask: weighted_img = w * img + (1-w) * C."""
    img_f = sitk.Cast(img, sitk.sitkFloat32)
    return w * img_f + (1.0 - w) * float(background_value)

def _resample_to_reference(
    img: sitk.Image,
    reference: sitk.Image,
    transform: sitk.Transform = sitk.Transform(),
    interp: int = sitk.sitkLinear,
    default_value: float = 0.0,
    pixel_id: Optional[int] = None,
) -> sitk.Image:
    """Resample `img` onto the reference grid (direction, spacing, size, origin)."""
    if pixel_id is None:
        pixel_id = img.GetPixelID()
    return sitk.Resample(img, reference, transform, interp, default_value, pixel_id)


def demons_with_soft_mask(
    fixed: sitk.Image,
    moving: sitk.Image,
    demons_fn,  # demons_registration(fixed, moving, **kwargs)
    *,
    fixed_mask: sitk.Image,
    moving_mask: Optional[sitk.Image] = None,
    fixed_falloff_mm: float = 8.0,
    moving_falloff_mm: Optional[float] = None,
    background_fixed: Optional[float] = None,
    background_moving: Optional[float] = None,
    propagate_fixed_weight_to_moving: bool = False,
    initial_transform: Optional[sitk.Transform] = None,
    interp_order: int = sitk.sitkLinear,
    **demons_kwargs,
):
    """
    Run Demons registration with soft mask weighting, focusing on ROI.

    - `fixed_mask` is required; `moving_mask` is optional.
    - If `propagate_fixed_weight_to_moving=True`, the fixed weight map (w_f) will be
      warped to moving space using `initial_transform` and used as the moving weight.
    - If both `moving_mask` and propagation are provided, their minimum is used as the moving weight.
    """

    # 1) Weight field (fixed image)
    if fixed_mask is None:
        raise ValueError("fixed_mask is required.")
    w_f = _soft_weight_from_mask(fixed_mask, falloff_mm=float(fixed_falloff_mm))  # in fixed space

    # 2) Weight field (moving image)
    w_m_list = []

    # 2a) From moving_mask (in moving space)
    if moving_mask is not None:
        mf = float(moving_falloff_mm) if moving_falloff_mm is not None else float(fixed_falloff_mm)
        w_m_from_own = _soft_weight_from_mask(moving_mask, falloff_mm=mf)  # in moving space
        w_m_list.append(w_m_from_own)

    # 2b) Propagate fixed weight to moving space (if requested)
    if propagate_fixed_weight_to_moving:
        # Need to map fixed weights (w_f) into moving space:
        # moving(x) corresponds to fixed(T(x)) → weight = w_f(T(x))
        # If initial_transform maps moving→fixed, we can directly resample w_f.
        T = initial_transform if initial_transform is not None else sitk.Transform()
        w_f_on_m = _resample_to_reference(
            w_f, reference=moving, transform=T, interp=sitk.sitkLinear, default_value=0.0, pixel_id=sitk.sitkFloat32
        )
        w_m_list.append(w_f_on_m)

    # 2c) Combine moving weights: use min if multiple sources exist
    if len(w_m_list) == 1:
        w_m = w_m_list[0]
    elif len(w_m_list) > 1:
        w_m = sitk.Minimum(w_m_list[0], w_m_list[1]) if len(w_m_list) == 2 else w_m_list[0]
        for i in range(2, len(w_m_list)):
            w_m = sitk.Minimum(w_m, w_m_list[i])
    else:
        w_m = None  # No moving weight → no soft mask applied to moving

    # 3) Infer background values
    if background_fixed is None:
        background_fixed = _infer_background_value(fixed, default_fallback=0.0)
    if background_moving is None:
        background_moving = _infer_background_value(moving, default_fallback=0.0)

    # 4) Apply soft mask to images
    fixed_soft = _apply_soft_mask(fixed, w_f, background_fixed)
    moving_soft = _apply_soft_mask(moving, w_m, background_moving) if w_m is not None else moving

    # 5) Call user-provided Demons function
    #    This may internally perform multi-resolution pyramids and resampling.
    #    Pass initial_transform, interp_order, and other parameters downstream if supported.
    result_img, out_tfm, dvf = demons_fn(
        fixed_image=fixed_soft,
        moving_image=moving_soft,
        initial_displacement_field=demons_kwargs.pop("initial_displacement_field", None),
        isotropic_resample=demons_kwargs.pop("isotropic_resample", False),
        resolution_staging=demons_kwargs.pop("resolution_staging", (8, 4, 1)),
        iteration_staging=demons_kwargs.pop("iteration_staging", (20, 20, 20)),
        regularisation_kernel_mm=demons_kwargs.pop("regularisation_kernel_mm", 1.0),
        smoothing_sigma_factor=demons_kwargs.pop("smoothing_sigma_factor", 1.0),
        smoothing_sigmas=demons_kwargs.pop("smoothing_sigmas", False),
        default_value=demons_kwargs.pop("default_value", background_moving),
        ncores=demons_kwargs.pop("ncores", 1),
        interp_order=interp_order,
        verbose=demons_kwargs.pop("verbose", False),
        # Other kwargs (if demons_fn supports them):
        **demons_kwargs,
    )

    return result_img, out_tfm, dvf


# example: only use fixed_mask
# reg_img, tfm, dvf = demons_with_soft_mask(
#     fixed=fixed_img,
#     moving=moving_img,
#     demons_fn=demons_registration,
#     fixed_mask=fixed_mask,
#     fixed_falloff_mm=10.0,
# )

# example: use fixed_mask and moving_mask
# reg_img, tfm, dvf = demons_with_soft_mask(
#     fixed=fixed_img,
#     moving=moving_img,
#     demons_fn=demons_registration,
#     fixed_mask=fixed_mask,
#     moving_mask=moving_mask,
#     propagate_fixed_weight_to_moving=True,
#     initial_transform=initial_affine,
#     fixed_falloff_mm=8.0,
#     moving_falloff_mm=6.0,
#     regularisation_kernel_mm=1.5,
#     resolution_staging=(8, 4, 2, 1),
#     iteration_staging=(50, 40, 30, 20),
#     isotropic_resample=True,
#     smoothing_sigmas=False,
#     ncores=8,
#     verbose=True,
# )
