"""
Registration pipeline module.

Provides a packaged registration pipeline integrating rigid and deformable registration workflows.
"""

from typing import Optional, Dict, Tuple
import numpy as np
import SimpleITK as sitk

from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.method.demons import demons_registration
from pydicomrt.utils.sitk_transform import resample_to_reference_image
from pydicomrt.reg.pipeline.preprocessing import (
    preprocess_image,
    align_image_extents,
    get_image_physical_extent,
    get_images_distance,
    get_initial_rigid_transform,
    create_reference_image_from_extent,
)


def _apply_transforms(
    image: sitk.Image,
    reference_image: sitk.Image,
    transforms: list,
    interpolator: int = sitk.sitkLinear,
    default_value: Optional[float] = None
) -> sitk.Image:
    """
    Apply a sequence of transforms to the image.

    Parameters
    ----------
    image : sitk.Image
        Image to be transformed.
    reference_image : sitk.Image
        Reference image (defines output space).
    transforms : list
        List of transforms to apply in order.
    interpolator : int, default = sitk.sitkLinear
        Interpolation method.
    default_value : Optional[float], default = None
        Default pixel value. If None, auto-detected (CT images use -1000).

    Returns
    -------
    sitk.Image
        Transformed image.
    """
    if default_value is None:
        # Auto-detect CT image
        default_value = 0.0
        try:
            arr_min = float(np.min(sitk.GetArrayViewFromImage(image)))
            if arr_min <= -1000.0:
                default_value = -1000.0
        except Exception:
            pass

    result_image = image

    # Apply each transform in order
    for transform in transforms:
        if transform is not None:
            result_image = sitk.Resample(
                result_image,
                reference_image,
                transform,
                interpolator,
                default_value,
                reference_image.GetPixelID()
            )

    return result_image


def registration_pipeline(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    perform_rigid: bool = True,
    perform_deformable: bool = True,
    preprocess_config: Optional[Dict] = None,
    rigid_kwargs: Optional[Dict] = None,
    deformable_kwargs: Optional[Dict] = None,
    resample_interpolator: int = sitk.sitkLinear,
    default_value: Optional[float] = None,
) -> Tuple[sitk.Image, Optional[sitk.Transform], Optional[sitk.Transform], Optional[sitk.Image]]:
    """
    Registration pipeline integrating rigid and deformable registration.

    Stages:
    1. Stage 0: Optional preprocessing before rigid (e.g. window clipping).
    2. Stage 1: Rigid registration (if enabled).
    3. Stage 2: Preprocessing before deform in rigid-aligned space (overlap crop + clip).
    4. Stage 3: Deformable registration (if enabled).
    5. Stage 4: Apply all transforms to the original moving_image.

    Parameters
    ----------
    fixed_image : sitk.Image
        Reference (fixed) image.
    moving_image : sitk.Image
        Image to be registered (moving image).
    perform_rigid : bool, default = True
        Whether to run rigid registration.
    perform_deformable : bool, default = True
        Whether to run deformable registration.
    preprocess_config : Optional[Dict], default = None
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

    rigid_kwargs : Optional[Dict], default = None
        Extra arguments for rigid_registration.
    deformable_kwargs : Optional[Dict], default = None
        Extra arguments for demons_registration.
    resample_interpolator : int, default = sitk.sitkLinear
        Interpolation for resampling.
    default_value : Optional[float], default = None
        Default value for resampling. If None, auto-detected (CT uses -1000).

    Returns
    -------
    registered_image : sitk.Image
        Final registered image in fixed_image space.
    rigid_transform : Optional[sitk.Transform]
        Rigid transform, or None if rigid was not run.
    deformable_transform : Optional[sitk.Transform]
        Deformable transform, or None if deformable was not run.
    deformation_field : Optional[sitk.Image]
        Deformation field, or None if deformable was not run.

    Examples
    --------
    >>> from pydicomrt.reg.pipeline import registration_pipeline
    >>>
    >>> # Basic: rigid + deformable
    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
    ...     fixed_image=ct_b_image,
    ...     moving_image=ct_a_image
    ... )
    >>>
    >>> # Rigid only
    >>> registered, rigid_tfm, _, _ = registration_pipeline(
    ...     fixed_image=ct_b_image,
    ...     moving_image=ct_a_image,
    ...     perform_deformable=False
    ... )
    >>>
    >>> # With preprocessing (flat format)
    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
    ...     fixed_image=ct_b_image,
    ...     moving_image=ct_a_image,
    ...     preprocess_config={'window_clip': [-10, 500]}
    ... )
    >>>
    >>> # With preprocessing (nested format)
    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(
    ...     fixed_image=ct_b_image,
    ...     moving_image=ct_a_image,
    ...     preprocess_config={
    ...         'rigid': {'window_clip': [-10, 500]},
    ...         'deform': {'window_clip': [-10, 500], 'align_extents': True}
    ...     }
    ... )
    """
    if not perform_rigid and not perform_deformable:
        raise ValueError("At least one of rigid or deformable registration must be performed")

    # Keep original images for final output
    original_moving_image = moving_image
    original_fixed_image = fixed_image

    # ===================================================================
    # Stage 0: Parse preprocess_config and optional preprocessing before rigid
    # ===================================================================
    cfg = preprocess_config or {}

    if isinstance(cfg, dict) and ("rigid" in cfg or "deform" in cfg):
        rigid_cfg = cfg.get("rigid")
        deform_cfg = cfg.get("deform")
    else:
        rigid_cfg = None
        deform_cfg = cfg if cfg else None

    # Preprocessing before rigid (intensity only, no spatial crop)
    if rigid_cfg is not None:
        fixed_for_rigid = preprocess_image(fixed_image, rigid_cfg)
        moving_for_rigid = preprocess_image(moving_image, rigid_cfg)
    else:
        fixed_for_rigid = fixed_image
        moving_for_rigid = moving_image

    # ===================================================================
    # Stage 1: Rigid alignment
    # ===================================================================
    rigid_transform = None
    initial_alignment_transform = None
    used_union_extent = False
    reference_for_rigid = fixed_image

    if perform_rigid:
        rigid_params = rigid_kwargs if rigid_kwargs is not None else {}

        fixed_min, fixed_max = get_image_physical_extent(fixed_for_rigid)
        moving_min, moving_max = get_image_physical_extent(moving_for_rigid)

        fixed_size = fixed_max - fixed_min
        moving_size = moving_max - moving_min
        max_image_diagonal = max(
            np.linalg.norm(fixed_size),
            np.linalg.norm(moving_size)
        )

        center_distance = get_images_distance(fixed_for_rigid, moving_for_rigid)
        distance_threshold = max_image_diagonal * 1.5
        need_pre_alignment = center_distance > distance_threshold

        if need_pre_alignment:
            initial_alignment_transform = get_initial_rigid_transform(
                fixed_for_rigid,
                moving_for_rigid
            )

            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(fixed_for_rigid)
            resampler.SetInterpolator(sitk.sitkLinear)
            
            if default_value is None:
                default_val = 0.0
                try:
                    arr_min = float(np.min(sitk.GetArrayViewFromImage(moving_for_rigid)))
                    if arr_min <= -1000.0:
                        default_val = -1000.0
                except Exception:
                    pass
            else:
                default_val = default_value
            
            resampler.SetDefaultPixelValue(default_val)
            resampler.SetTransform(initial_alignment_transform)
            
            moving_pre_aligned = resampler.Execute(moving_for_rigid)
            moving_min, moving_max = get_image_physical_extent(moving_pre_aligned)
        else:
            moving_pre_aligned = moving_for_rigid

        union_min = np.minimum(fixed_min, moving_min)
        union_max = np.maximum(fixed_max, moving_max)
        union_size = union_max - union_min
        max_reasonable_size = max(
            np.linalg.norm(fixed_size),
            np.linalg.norm(moving_size)
        ) * 2.0
        use_union = np.linalg.norm(union_size) <= max_reasonable_size

        if use_union:
            reference_for_rigid = create_reference_image_from_extent(
                union_min,
                union_max,
                fixed_image.GetSpacing(),
                fixed_image.GetDirection(),
                fixed_image.GetPixelID()
            )
            used_union_extent = True

            fixed_in_union_space = resample_to_reference_image(
                reference_for_rigid,
                fixed_for_rigid
            )
            moving_in_union_space = resample_to_reference_image(
                reference_for_rigid,
                moving_pre_aligned
            )
            rigid_transform = rigid_registration(
                fixed_in_union_space,
                moving_in_union_space,
                **rigid_params
            )
        else:
            reference_for_rigid = fixed_image
            rigid_transform = rigid_registration(
                fixed_for_rigid,
                moving_pre_aligned,
                **rigid_params
            )
        
        if initial_alignment_transform is not None:
            composite_transform = sitk.CompositeTransform([
                initial_alignment_transform,
                rigid_transform
            ])
            rigid_transform = composite_transform

        if default_value is None:
            try:
                arr_min = float(np.min(sitk.GetArrayViewFromImage(moving_image)))
                default_val = -1000.0 if arr_min <= -1000.0 else 0.0
            except Exception:
                default_val = 0.0
        else:
            default_val = default_value

        moving_rigid = sitk.Resample(
            moving_image,
            reference_for_rigid,
            rigid_transform,
            resample_interpolator,
            default_val,
            fixed_image.GetPixelID()
        )
    else:
        moving_rigid = moving_image

    # ===================================================================
    # Stage 2: Preprocessing before deform (overlap crop + clip)
    # ===================================================================
    if used_union_extent:
        fixed_in_union_space = resample_to_reference_image(
            reference_for_rigid,
            fixed_image
        )
        fixed_for_deform = fixed_in_union_space
    else:
        fixed_for_deform = fixed_image
    
    moving_for_deform = moving_rigid
    used_align_extents = False

    if deform_cfg is not None:
        if deform_cfg.get("align_extents", False):
            fixed_for_deform, moving_for_deform, _ = align_image_extents(
                fixed_for_deform,
                moving_for_deform,
                default_value=default_value,
                use_initial_rigid=False,
            )
            used_align_extents = True

        deform_cfg_for_clip = {k: v for k, v in deform_cfg.items() if k != "align_extents"}
        if deform_cfg_for_clip:
            fixed_for_deform = preprocess_image(fixed_for_deform, deform_cfg_for_clip)
            moving_for_deform = preprocess_image(moving_for_deform, deform_cfg_for_clip)

    # ===================================================================
    # Stage 3: Deformable registration
    # ===================================================================
    deformable_transform = None
    deformation_field = None

    if perform_deformable:
        deformable_params = deformable_kwargs if deformable_kwargs is not None else {}
        _, deformable_transform, deformation_field = demons_registration(
            fixed_for_deform,
            moving_for_deform,
            **deformable_params
        )

    # ===================================================================
    # Stage 4: Final output and transform composition
    # ===================================================================
    if used_align_extents:
        reference_for_final = fixed_for_deform
    elif used_union_extent:
        reference_for_final = reference_for_rigid
    else:
        reference_for_final = fixed_image

    original_moving_resampled = resample_to_reference_image(
        reference_for_final,
        original_moving_image
    )

    transforms_to_apply = []
    if rigid_transform is not None:
        transforms_to_apply.append(rigid_transform)
    if deformable_transform is not None:
        transforms_to_apply.append(deformable_transform)

    if transforms_to_apply:
        registered_image = _apply_transforms(
            original_moving_resampled,
            reference_for_final,
            transforms_to_apply,
            resample_interpolator,
            default_value
        )
    else:
        registered_image = original_moving_resampled

    if used_align_extents or used_union_extent:
        if fixed_image.GetSize() != reference_for_final.GetSize():
            registered_image = resample_to_reference_image(
                fixed_image,
                registered_image
            )

    return registered_image, rigid_transform, deformable_transform, deformation_field
