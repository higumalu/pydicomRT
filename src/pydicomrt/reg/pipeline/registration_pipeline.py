"""
Registration pipeline module.

Provides a packaged registration pipeline integrating rigid and deformable registration workflows.
"""

import logging
from typing import Optional, Dict, List, Tuple

import SimpleITK as sitk

from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.method.demons import demons_registration
from pydicomrt.reg.pipeline.preprocessing import (
    preprocess_image,
    align_image_extents,
    infer_default_pixel_value,
)

logger = logging.getLogger(__name__)


def compose_transforms(transforms: List[Optional[sitk.Transform]]) -> Optional[sitk.Transform]:
    """
    Compose a stage-ordered list of transforms into a single resampling transform.

    ``transforms`` is given in the order the stages were computed (rigid first, then
    deformable). A resampling transform maps points from the *fixed* grid back into the
    moving image, so the stages must be applied in reverse: the deformable transform maps
    a fixed-space point into the rigid-aligned space, and the rigid transform then maps
    that into the original moving space.

    :class:`SimpleITK.CompositeTransform` evaluates its list back-to-front, which matches
    the stage order exactly, so the list is passed through unchanged.

    Parameters
    ----------
    transforms : List[Optional[sitk.Transform]]
        Stage-ordered transforms. ``None`` entries are skipped.

    Returns
    -------
    Optional[sitk.Transform]
        A single transform equivalent to applying every stage, or None if the list holds
        no transform. A single-element list is returned as-is rather than wrapped.
    """
    applicable = [transform for transform in transforms if transform is not None]

    if not applicable:
        return None
    if len(applicable) == 1:
        return applicable[0]
    return sitk.CompositeTransform(applicable)


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
        Padding value for resampling. If None, inferred per image from its intensity range
        (CT-like images use -1000, everything else 0). Pass it explicitly for MR/PET.

    Returns
    -------
    registered_image : sitk.Image
        Final registered image, always on the ``fixed_image`` grid and carrying the
        ``moving_image`` pixel type. Produced by a single resampling of the untouched
        ``moving_image`` through the composed transform, so no preprocessing (window
        clipping, extent cropping) leaks into the output.
    rigid_transform : Optional[sitk.Transform]
        Rigid transform, or None if rigid was not run.
    deformable_transform : Optional[sitk.Transform]
        Deformable transform, or None if deformable was not run.
    deformation_field : Optional[sitk.Image]
        Deformation field, or None if deformable was not run.

    Notes
    -----
    ``rigid_transform`` and ``deformable_transform`` are resampling transforms: they map a
    point on the fixed grid back into the moving image. Compose them with
    :func:`compose_transforms` (stage order) rather than resampling once per stage.

    The rigid stage registers the images in their own grids. Resampling them onto a shared
    grid first -- which this used to do -- leaves ``CenteredTransformInitializer`` with two
    identical grids and therefore nothing to estimate, so any offset larger than the
    optimizer's own capture range is lost. With no preprocessing configured, the result is
    identical to calling :func:`rigid_registration` directly.

    Examples
    --------
    Rigid then deformable:

    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(  # doctest: +SKIP
    ...     fixed_image=ct_image, moving_image=cbct_image)

    Rigid only:

    >>> registered, rigid_tfm, _, _ = registration_pipeline(             # doctest: +SKIP
    ...     fixed_image=ct_image, moving_image=cbct_image,
    ...     perform_deformable=False)

    Preprocessing applied to every stage (flat form):

    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(  # doctest: +SKIP
    ...     fixed_image=ct_image, moving_image=cbct_image,
    ...     preprocess_config={'window_clip': [-10, 500]})

    Preprocessing configured per stage (nested form):

    >>> registered, rigid_tfm, deform_tfm, dvf = registration_pipeline(  # doctest: +SKIP
    ...     fixed_image=ct_image, moving_image=cbct_image,
    ...     preprocess_config={
    ...         'rigid': {'window_clip': [-10, 500]},
    ...         'deform': {'window_clip': [-10, 500], 'align_extents': True},
    ...     })
    """
    if not perform_rigid and not perform_deformable:
        raise ValueError("At least one of rigid or deformable registration must be performed")

    # Keep original images for final output
    original_moving_image = moving_image

    # Resolve padding values once, on the raw images. Inference reads the intensity
    # minimum, so it must happen before window clipping pushes it above the air value.
    moving_default_value = infer_default_pixel_value(moving_image, default_value)

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

    if perform_rigid:
        rigid_params = rigid_kwargs if rigid_kwargs is not None else {}

        # Register the images in their own grids.
        #
        # rigid_registration starts from a CenteredTransformInitializer, whose entire
        # contribution is the offset between the two images' geometric centres. Resampling
        # either image onto a grid shared with the other -- a union extent, or simply the
        # fixed grid -- makes that offset identically zero, and the optimizer is left to
        # discover the whole misalignment unaided. This stage used to do exactly that, and
        # on a real CT/CBCT pair with a ~180 mm couch offset it finished further from the
        # answer than doing nothing at all.
        #
        # ITK does not need a shared grid: the metric samples in fixed space and maps
        # through the transform, so images in different physical spaces register directly.
        rigid_transform = rigid_registration(
            fixed_for_rigid,
            moving_for_rigid,
            **rigid_params
        )

        moving_rigid = sitk.Resample(
            moving_image,
            fixed_image,
            rigid_transform,
            resample_interpolator,
            moving_default_value,
            moving_image.GetPixelID()
        )
    else:
        moving_rigid = moving_image

    # ===================================================================
    # Stage 2: Preprocessing before deform (overlap crop + clip)
    # ===================================================================
    # The rigid stage has already brought the moving image onto the fixed grid, which is
    # also the domain the displacement field needs to cover and the grid the result is
    # returned on.
    fixed_for_deform = fixed_image
    moving_for_deform = moving_rigid

    if deform_cfg is not None:
        if deform_cfg.get("align_extents", False):
            fixed_for_deform, moving_for_deform, _ = align_image_extents(
                fixed_for_deform,
                moving_for_deform,
                default_value=default_value,
                use_initial_rigid=False,
            )

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
        deformable_params = dict(deformable_kwargs or {})
        deformable_params.setdefault("default_value", moving_default_value)
        _, deformable_transform, deformation_field = demons_registration(
            fixed_for_deform,
            moving_for_deform,
            **deformable_params
        )

    # ===================================================================
    # Stage 4: Final output and transform composition
    # ===================================================================
    # Every stage transform lives in physical space, so they compose into one transform that
    # maps the fixed grid straight back into the original moving image. Applying it in a
    # single resampling step matters twice over: resampling the moving image onto the fixed
    # grid *first* would throw away the very voxels the transform is about to bring into
    # view, and each additional resampling adds another round of interpolation blur.
    final_transform = compose_transforms([rigid_transform, deformable_transform])

    if final_transform is None:
        # Unreachable while the guard above requires at least one registration stage.
        final_transform = sitk.Transform()

    registered_image = sitk.Resample(
        original_moving_image,
        fixed_image,
        final_transform,
        resample_interpolator,
        moving_default_value,
        original_moving_image.GetPixelID()
    )

    return registered_image, rigid_transform, deformable_transform, deformation_field
