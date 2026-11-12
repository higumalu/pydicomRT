import logging
from typing import Optional, Sequence, Tuple, Union, List
import numpy as np
import SimpleITK as sitk
from .common import deformable_registration_command_iteration, smooth_and_resample

logger = logging.getLogger(__name__)




def multiscale_demons(
    registration_algorithm,
    fixed_image,
    moving_image,
    initial_transform=None,
    initial_displacement_field=None,
    isotropic_resample=None,
    resolution_staging=None,
    smoothing_sigmas=None,
    iteration_staging=None,
    interp_order=sitk.sitkLinear,
    default_value=None,
    regularization_kernel_mm=None,
):
    """
    Perform deformable image registration using a multi-scale (pyramid) demons strategy.

    This function applies a given registration algorithm (e.g., Demons) in a coarse-to-fine
    hierarchical fashion. At each pyramid level, the fixed and moving images are smoothed and
    resampled, and the registration is run iteratively. The displacement vector field (DVF)
    is refined across scales until the final resolution is reached.

    Parameters
    ----------
    registration_algorithm : SimpleITK registration filter
        A demons-based registration algorithm object with
        `Execute(fixed_image, moving_image) -> displacement_field` interface.
    fixed_image : sitk.Image
        The reference image. The resulting transform maps points from this image
        to the moving image's spatial domain.
    moving_image : sitk.Image
        The moving image. The algorithm estimates a transformation from
        fixed_image → moving_image space.
    initial_transform : Optional[sitk.Transform], default = None
        Initial transform to initialize the displacement field. Ignored if
        `initial_displacement_field` is provided.
    initial_displacement_field : Optional[sitk.Image], default = None
        Optional displacement field to initialize registration. If provided,
        it overrides `initial_transform`.
    isotropic_resample : bool, default = None
        If True, interpret `resolution_staging` as isotropic voxel sizes (mm).
        If False, interpret as shrink factors (relative downsampling).
    resolution_staging : Sequence[int or float], default = None
        Per-level downsampling factors or target isotropic resolutions (mm),
        depending on `isotropic_resample`.
        - If `isotropic_resample=False`: values are shrink factors (typically integers, e.g., [8, 4, 1]).
        - If `isotropic_resample=True`: values are target isotropic voxel sizes in mm (floats, e.g., [2.5, 1.0, 0.5]).
    smoothing_sigmas : Sequence[float], default = None
        Per-level Gaussian smoothing sigmas (in physical units, mm) applied to images
        before resampling at each pyramid level. This is separate from the demons
        regularization kernel. Length must match `resolution_staging`.
    iteration_staging : Sequence[int], default = None
        Number of iterations to run at each pyramid level.
    interp_order : int, default = sitk.sitkLinear
        Interpolator used for image resampling.
    default_value : float, optional
        Padding used when warping the moving image at each level. Inferred from raw
        moving intensities when omitted.
    regularization_kernel_mm : sequence of float, optional
        Displacement smoothing widths in mm, converted to voxels at each pyramid level.
        When omitted, use the filter's existing voxel widths.

    Returns
    -------
    dvf_total : sitk.Image
        Final displacement vector field (DVF) in the fixed image space.

    Notes
    -----
    - Uses Gaussian smoothing at each scale to stabilize registration.
    - The final DVF is smoothed at each iteration and accumulated across levels.
    - The DVF is accumulated as sitkVectorFloat32 and returned as sitkVectorFloat64, the type
      SimpleITK's DisplacementFieldTransform requires.
    """
    if default_value is None:
        default_value = -1000.0 if np.min(sitk.GetArrayViewFromImage(moving_image)) <= -1000 else 0.0
    # The accumulated DVF is held as float32. Only the demons filter and DisplacementFieldTransform
    # need float64, and a full-resolution vector field is the largest thing held here: 12 rather
    # than 24 bytes per voxel, at a rounding error of ~1e-5 mm. Every step warps with _warp,
    # which takes the field directly -- Resample would need a DisplacementFieldTransform, i.e. a
    # float64 copy of it.
    field_type = sitk.sitkVectorFloat32
    if initial_displacement_field is not None:
        dvf_total = sitk.Cast(initial_displacement_field, field_type)
    elif initial_transform:
        dvf_total = sitk.TransformToDisplacementField(
            initial_transform,
            field_type,
            fixed_image.GetSize(),
            fixed_image.GetOrigin(),
            fixed_image.GetSpacing(),
            fixed_image.GetDirection(),
        )
    else:
        dvf_total = None  # zero, created on the coarsest grid below

    # Multi-scale loop, coarsest → finest. Each level's images are built when the level runs
    # rather than as a pyramid up front, so the full-size ones are not held throughout.
    for resolution, smoothing_sigma, iters in zip(
        resolution_staging, smoothing_sigmas, iteration_staging, strict=True
    ):
        # `resolution` is an isotropic voxel size (mm) or a shrink factor
        level = dict(
            isotropic_voxel_size_mm=resolution if isotropic_resample else None,
            shrink_factor=None if isotropic_resample else resolution,
            smoothing_sigma=smoothing_sigma,
            interpolator=interp_order,
        )
        f_image = smooth_and_resample(fixed_image, **level)
        if regularization_kernel_mm is not None:
            registration_algorithm.SetStandardDeviations(
                (np.asarray(regularization_kernel_mm) / np.asarray(f_image.GetSpacing())).tolist()
            )

        # Resample DVF to current resolution
        if dvf_total is None:
            dvf_total = sitk.Image(f_image.GetSize(), field_type, f_image.GetDimension())
            dvf_total.CopyInformation(f_image)
        else:
            dvf_total = sitk.Resample(dvf_total, f_image)

        # Warp the moving image through the DVF and run demons at this scale
        m_image = _warp(
            smooth_and_resample(moving_image, **level), dvf_total, f_image, interp_order, float(default_value)
        )
        registration_algorithm.SetNumberOfIterations(iters)
        dvf_iter = registration_algorithm.Execute(f_image, m_image)

        # Moving was warped through T; the refinement D therefore gives T(D(x)).
        # u_total(x) = v(x) + u(x + v(x)), not u(x) + v(x + u(x)).
        dvf_total = _warp(dvf_total, dvf_iter, f_image, sitk.sitkLinear, 0.0)
        dvf_total += sitk.Cast(dvf_iter, field_type)

        # Regularize/smooth the DVF
        sigma_mm = np.asarray(registration_algorithm.GetStandardDeviations()) * np.asarray(f_image.GetSpacing())
        # The filter takes voxel sigmas; the standalone Gaussian takes physical mm.
        dvf_total = sitk.SmoothingRecursiveGaussian(dvf_total, sigma_mm.tolist())

    # Resample final DVF to fixed image grid
    dvf_total = sitk.Resample(dvf_total, fixed_image)

    return sitk.Cast(dvf_total, sitk.sitkVectorFloat64)


def _warp(image, displacement_field, reference, interpolator, default_value):
    """``image`` sampled at ``x + displacement_field(x)`` for every point ``x`` of ``reference``.

    The same as ``sitk.Resample`` through a ``DisplacementFieldTransform`` of the field, without
    the float64 copy that transform requires.
    """
    warp = sitk.WarpImageFilter()
    warp.SetOutputParameteresFromImage(reference)  # sic: SimpleITK's spelling
    warp.SetInterpolator(int(interpolator))
    warp.SetEdgePaddingValue(default_value)
    return warp.Execute(image, displacement_field)




#: Defaults that depend on ``metric``: the intensity values are the long-standing ones; the
#: MIND values were chosen on brain MR -> planning CT pairs, where the intensity ones invent
#: millimetres of motion (an 8x coarsest level and a 1.5 mm kernel let the descriptor
#: difference pull the field around).
METRIC_DEFAULTS = {
    "intensity": {"resolution_staging": (8, 4, 1), "iteration_staging": (10, 10, 10),
                  "regularization_kernel_mm": 1.5, "smoothing_sigma_factor": 1.0},
    "mind": {"resolution_staging": (4, 2, 1), "iteration_staging": (100, 100, 200),
             "regularization_kernel_mm": 3.0, "smoothing_sigma_factor": 0.5},
}


def demons_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    resolution_staging: Optional[Sequence[int]] = None,
    iteration_staging: Optional[Sequence[int]] = None,
    isotropic_resample: bool = False,
    initial_displacement_field: Optional[sitk.Image] = None,
    regularization_kernel_mm: Optional[Union[float, Sequence[float]]] = None,
    smoothing_sigma_factor: Optional[float] = None,
    smoothing_sigmas: Optional[Union[float, Sequence[float], bool]] = False,
    default_value: Optional[float] = None,
    ncores: int = 1,
    interp_order: int = sitk.sitkLinear,
    verbose: bool = False,
    backend: str = "simpleitk",
    device: Optional[str] = None,
    metric: str = "intensity",
) -> Tuple[sitk.Image, sitk.Transform, sitk.Image]:
    """
    Deformable registration via Fast Symmetric-Forces Demons (SimpleITK) with multi-resolution strategy.

    The function registers `moving_image` to `fixed_image` by estimating a dense displacement
    field using a coarse-to-fine (image pyramid) schedule. It then applies the resulting
    displacement-field transform to resample the moving image into the fixed image space.
    Inputs may have different grids; moving data are sampled onto the fixed pyramid grid
    through the current transform at each level. This does not estimate a rigid alignment.

    Demons compares intensities by default, which is only meaningful when the same tissue has
    the same intensity in both images: CT to CT, CT to CBCT. For MR to CT use
    ``metric="mind"`` (with ``backend="torch"``); the intensity metric produces several
    millimetres of motion that is not there.

    Parameters
    ----------
    fixed_image : sitk.Image
        Reference (target) image. Output is defined on this image grid.
    moving_image : sitk.Image
        Image to be warped into the fixed image space.
    resolution_staging : Sequence[int or float], optional
        Pyramid schedule. Length of this sequence defines the number of pyramid levels.
        ``None`` means (8, 4, 1) for ``metric="intensity"`` and (4, 2, 1) for ``"mind"``.
        - If `isotropic_resample=False`: values are downsampling shrink factors (typically integers, e.g., [8, 4, 1]).
        - If `isotropic_resample=True`: values are target isotropic voxel sizes in millimeters (floats, e.g., [2.5, 1.0, 0.5]).
    iteration_staging : Sequence[int], optional
        Number of iterations per pyramid level. Must match the length of `resolution_staging`.
        ``None`` means (10, 10, 10) for ``metric="intensity"`` and (100, 100, 200) for
        ``"mind"``. Either way a level stops early once its updates become negligible.
    isotropic_resample : bool, default = False
        If True, interpret `resolution_staging` as isotropic voxel sizes (mm); otherwise as shrink factors.
    initial_displacement_field : Optional[sitk.Image], default = None
        Initial DVF (vector image). If provided, it overrides any `initial_transform` usage
        inside `multiscale_demons` and is resampled onto the fixed grid when needed.
    regularization_kernel_mm : float or Sequence[float], optional
        Standard deviation(s) in millimeters used by demons for smoothing the update and DVF.
        Scalar broadcasts to all axes; sequence length must be 1 or image dimension.
        ``None`` means 1.5 for ``metric="intensity"`` and 3.0 for ``"mind"``.
    smoothing_sigma_factor : float, optional
        Fallback factor to derive `smoothing_sigmas` per level if `smoothing_sigmas` is not provided.
        Each level's sigma = `resolution_staging[level] * smoothing_sigma_factor`.
        ``None`` means 1.0 for ``metric="intensity"`` and 0.5 for ``"mind"``.
    smoothing_sigmas : float or Sequence[float] or bool, default = False
        Per-level Gaussian sigmas (in physical units, mm) used for image smoothing
        prior to resampling at each pyramid level. This is separate from the demons
        regularization kernel (`regularization_kernel_mm`).
        - False/None: derived from `resolution_staging` and `smoothing_sigma_factor`
        - float: same sigma for all levels
        - sequence: length must match number of levels
    default_value : Optional[float], default = None
        Out-of-bounds value during final resampling. If None, use 0. For CT-like images
        (min <= -1000), default to -1000.
    ncores : int, default = 1
        Number of threads for demons filter. Ignored by the torch backend.
    interp_order : int, default = sitk.sitkLinear
        Interpolator for resampling (e.g., sitk.sitkNearestNeighbor, sitk.sitkLinear, sitk.sitkBSpline).
        The torch backend interpolates linearly inside the registration and accepts only
        `sitk.sitkLinear`.
    verbose : bool, default = False
        If True and a global callback `deformable_registration_command_iteration` is defined,
        prints per-iteration metric values. Ignored by the torch backend.
    backend : {"simpleitk", "torch"}, default = "simpleitk"
        Which implementation runs the demons iterations. ``"torch"`` runs the same algorithm
        on PyTorch -- on a GPU when one is available -- and needs the optional dependency:
        ``pip install "pydicomrt[gpu]"``. It reproduces the SimpleITK field to within a few
        hundredths of a millimetre and supports 3D images only. Everything else, including
        the returned objects, is the same for both.
    device : str, optional
        PyTorch device for ``backend="torch"``, e.g. ``"cuda"``, ``"cuda:1"`` or ``"cpu"``.
        ``None`` uses CUDA when available and otherwise falls back to the CPU with a warning.
        Must be ``None`` for the SimpleITK backend.
    metric : {"intensity", "mind"}, default = "intensity"
        What the demons force compares. ``"intensity"`` is ITK's demons: the intensity
        difference, so single-modality only. ``"mind"`` compares MIND-SSC descriptors
        (Heinrich et al., 2013) -- at every voxel, how alike the neighbouring patches are --
        which depend on where the edges are and not on which side is brighter, so the images
        may be of different modalities (MR to CT) and need no intensity normalisation.
        ``"mind"`` needs ``backend="torch"`` and changes the defaults noted above; see Notes.

    Returns
    -------
    registered_image : sitk.Image
        The moving image warped into the fixed image space.
    output_transform : sitk.Transform
        DisplacementFieldTransform built from the final DVF.
    deformation_field : sitk.Image
        The final dense displacement field (vector image) aligned to `fixed_image`.

    Raises
    ------
    ValueError
        On invalid arguments (length mismatches, non-positive values, dimension mismatch, etc.).
    ImportError
        If ``backend="torch"`` is requested and PyTorch is not installed.
    RuntimeError
        If helper `multiscale_demons(...)` is not available in the current scope.

    Notes
    -----
    - Internally casts images to Float32 for demons, then casts result back to the original
      moving pixel type after resampling.
    - Regularization sigmas are specified in *voxels* to SimpleITK; this function converts
      millimeter inputs to voxel units using the fixed image spacing.
    - GPU memory for ``backend="torch"`` grows linearly with the voxel count of the finest
      pyramid level, at roughly 180 bytes per voxel: about 4.5 GB for 512 x 512 x 100. For a
      volume that does not fit, make the finest level coarser -- e.g. ``isotropic_resample=True``
      with ``resolution_staging=(8, 4, 1.5)`` -- the returned field is still on the fixed grid.
      Running out raises ``torch.OutOfMemoryError`` with an estimate.
    - RAM for ``backend="simpleitk"`` also follows the finest level, at roughly 150 bytes per
      voxel -- about 6 GB for 512 x 512 x 160, most of it inside ITK's demons filter. The same
      coarser finest level is the remedy: at 2 mm that volume needs 2.5 GB.
    - Iterations are best spent at the finest level: on synthetic tests, raising the coarse
      levels' counts made the result worse, raising the finest level's made it better.
    - ``metric="mind"`` carries one displacement field through all pyramid levels instead of
      composing a refinement per level as ITK's multiscale scheme does, so it does not
      reproduce any SimpleITK result. It needs about 360 bytes of device memory per voxel of
      the finest level -- twice the intensity metric -- so register a planning CT at
      1-1.5 mm (``isotropic_resample=True``, e.g. ``resolution_staging=(6, 3, 1.5)``) rather
      than at its native in-plane resolution. MIND sees edges, so the result is only as good
      as the edges the two images share; on 2D MR series with 3-5 mm slices, motion along the
      slice direction is mostly filled in by the regularisation.
    - Validated on public brain MR/CT pairs (SynthRAD2023) against known deformations only. A
      deformable result between MR and CT also moves real anatomical change -- tumour growth,
      oedema, resection cavities -- out of sight; review it before contouring on it.
    """

    if metric not in METRIC_DEFAULTS:
        raise ValueError(f"`metric` must be 'intensity' or 'mind', got {metric!r}")
    if metric == "mind" and backend != "torch":
        raise ValueError("metric='mind' runs on backend='torch' only; pass backend='torch' "
                         "(pip install \"pydicomrt[gpu]\")")
    defaults = METRIC_DEFAULTS[metric]
    if resolution_staging is None:
        resolution_staging = defaults["resolution_staging"]
    if iteration_staging is None:
        iteration_staging = defaults["iteration_staging"]
    if regularization_kernel_mm is None:
        regularization_kernel_mm = defaults["regularization_kernel_mm"]
    if smoothing_sigma_factor is None:
        smoothing_sigma_factor = defaults["smoothing_sigma_factor"]

    # ---- Basic validation on lists/levels/cores ----
    if not isinstance(resolution_staging, Sequence) or len(resolution_staging) == 0:
        raise ValueError("`resolution_staging` must be a non-empty sequence of positive numbers.")
    if not isinstance(iteration_staging, Sequence) or len(iteration_staging) == 0:
        raise ValueError("`iteration_staging` must be a non-empty sequence of positive integers.")
    if len(iteration_staging) != len(resolution_staging):
        raise ValueError("`iteration_staging` length must match `resolution_staging` length.")
    if any(float(x) <= 0 for x in resolution_staging):
        raise ValueError("All values in `resolution_staging` must be positive numbers.")
    if any(int(x) <= 0 for x in iteration_staging):
        raise ValueError("All values in `iteration_staging` must be positive integers.")
    if ncores <= 0:
        raise ValueError("`ncores` must be a positive integer.")
    if backend not in ("simpleitk", "torch"):
        raise ValueError(f"`backend` must be 'simpleitk' or 'torch', got {backend!r}")
    if backend == "simpleitk" and device is not None:
        raise ValueError("`device` applies to backend='torch' only")
    if backend == "torch":
        if fixed_image.GetDimension() != 3:
            raise ValueError("backend='torch' supports 3D images only")
        if int(interp_order) != sitk.sitkLinear:
            raise ValueError("backend='torch' interpolates linearly; pass interp_order=sitk.sitkLinear")

    # ---- Dimensionality checks ----
    dim = fixed_image.GetDimension()
    if moving_image.GetDimension() != dim:
        raise ValueError("`fixed_image` and `moving_image` must have the same dimension.")
    if initial_displacement_field is not None and initial_displacement_field.GetDimension() != dim:
        raise ValueError("`initial_displacement_field` must have the same dimension as the images.")
    if initial_displacement_field is not None and initial_displacement_field.GetNumberOfComponentsPerPixel() != dim:
        raise ValueError("initial_displacement_field must have one component per spatial dimension")

    # ---- Prepare per-level smoothing sigmas (in physical units) ----
    num_levels = len(resolution_staging)
    if smoothing_sigmas is False or smoothing_sigmas is None:
        # Derive from schedule × factor, e.g., (8,4,1) × 1.0 → (8.0,4.0,1.0)
        smoothing_sigmas_list: List[float] = [
            float(level) * float(smoothing_sigma_factor) for level in resolution_staging
        ]
    elif isinstance(smoothing_sigmas, (int, float)):
        smoothing_sigmas_list = [float(smoothing_sigmas)] * num_levels
    elif isinstance(smoothing_sigmas, Sequence):
        if len(smoothing_sigmas) != num_levels:
            raise ValueError("When providing a sequence, `smoothing_sigmas` length must match the number of levels.")
        smoothing_sigmas_list = [float(s) for s in smoothing_sigmas]
    else:
        raise ValueError("`smoothing_sigmas` must be False/None, a float, or a sequence of floats.")

    # ---- Convert regularisation sigma from mm → voxels for demons filter ----
    if isinstance(regularization_kernel_mm, (int, float)):
        kernel_mm = np.full(dim, float(regularization_kernel_mm), dtype=np.float64)
    else:
        ker = np.asarray(regularization_kernel_mm, dtype=np.float64).ravel()
        if ker.size not in (1, dim):
            raise ValueError(
                f"`regularization_kernel_mm` must be a float or a sequence of length 1 or {dim}."
            )
        kernel_mm = np.full(dim, float(ker[0]), dtype=np.float64) if ker.size == 1 else ker
    if not np.isfinite(kernel_mm).all() or np.any(kernel_mm <= 0):
        raise ValueError("regularization_kernel_mm must contain positive finite values")

    # Keep the original moving image pixel type to cast result back later
    moving_image_type = moving_image.GetPixelID()

    # Demons prefers float images; cast to Float32 if necessary
    if fixed_image.GetPixelID() != sitk.sitkFloat32:
        fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    if moving_image.GetPixelID() != sitk.sitkFloat32:
        moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # ---- Configure demons filter ----
    registration_method = sitk.FastSymmetricForcesDemonsRegistrationFilter()
    registration_method.SetNumberOfThreads(int(ncores))
    registration_method.SetSmoothUpdateField(True)          # smooth incremental updates
    registration_method.SetSmoothDisplacementField(True)    # smooth accumulated DVF

    # Convert mm to voxels using fixed image spacing (what the filter expects)
    spacing = np.asarray(fixed_image.GetSpacing(), dtype=np.float64)
    if spacing.size != dim:
        raise ValueError("Unexpected spacing dimensionality mismatch.")
    kernel_vox = (kernel_mm / spacing).tolist()
    registration_method.SetStandardDeviations(kernel_vox)

    # Optional iteration callback for progress printing
    if verbose and "deformable_registration_command_iteration" in globals():
        registration_method.AddCommand(
            sitk.sitkIterationEvent,
            lambda: globals()["deformable_registration_command_iteration"](registration_method),
        )

    # ---- Run multi-scale demons to obtain the final DVF ----
    if backend == "torch":
        try:
            from ._demons_torch import multiscale_demons_torch
        except ImportError as error:
            raise ImportError(
                "backend='torch' needs PyTorch, which pydicomrt does not install by default: "
                "pip install \"pydicomrt[gpu]\""
            ) from error
        deformation_field = multiscale_demons_torch(
            fixed_image,
            moving_image,
            resolution_staging=tuple(float(x) for x in resolution_staging),
            smoothing_sigmas=smoothing_sigmas_list,
            iteration_staging=tuple(int(x) for x in iteration_staging),
            isotropic_resample=bool(isotropic_resample),
            regularization_kernel_mm=kernel_mm,
            initial_displacement_field=initial_displacement_field,
            default_value=default_value,
            device=device,
            metric=metric,
        )
    else:
        try:
            # Preserve float values for isotropic_resample mode, but convert to appropriate types
            if isotropic_resample:
                # For isotropic mode, keep as floats (mm values)
                resolution_staging_processed = tuple(float(x) for x in resolution_staging)
            else:
                # For shrink factor mode, can use integers
                resolution_staging_processed = tuple(float(x) for x in resolution_staging)
        
            deformation_field: sitk.Image = multiscale_demons(
                registration_algorithm=registration_method,
                fixed_image=fixed_image,
                moving_image=moving_image,
                resolution_staging=resolution_staging_processed,
                smoothing_sigmas=smoothing_sigmas_list,
                iteration_staging=tuple(int(x) for x in iteration_staging),
                isotropic_resample=bool(isotropic_resample),
                initial_displacement_field=initial_displacement_field,
                interp_order=int(interp_order),
                default_value=default_value,
                regularization_kernel_mm=kernel_mm,
            )
        except NameError as e:
            # Provide a clearer error if helper function is missing
            raise RuntimeError(
                "Required helper `multiscale_demons(...)` is not defined in the current scope."
            ) from e
    # The filter keeps its last run's buffers, about two full-resolution fields, until it is
    # destroyed; it is not needed past this point.
    del registration_method

    # ---- Build a resampler to warp the moving image using the DVF ----
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed_image)      # output grid = fixed image grid
    resampler.SetInterpolator(int(interp_order))  # interpolation scheme

    # Heuristic default value: use -1000 for CT-like images, else 0
    if default_value is None:
        default_value = 0.0
        try:
            arr_min = float(np.min(sitk.GetArrayViewFromImage(moving_image)))
            if arr_min <= -1000.0:
                default_value = -1000.0
        except Exception:
            # If array view is not available, keep default 0.0
            pass
    resampler.SetDefaultPixelValue(float(default_value))

    # Create a displacement-field transform from the DVF (cast to Float64 for SITK transform)
    df64 = sitk.Cast(deformation_field, sitk.sitkVectorFloat64)
    output_transform = sitk.DisplacementFieldTransform(df64)
    resampler.SetTransform(output_transform)

    # Execute resampling to get the registered (warped) image
    registered_image = resampler.Execute(moving_image)

    # Copy spatial metadata from the fixed image; cast back to the original moving pixel type
    registered_image.CopyInformation(fixed_image)
    registered_image = sitk.Cast(registered_image, moving_image_type)

    return registered_image, output_transform, deformation_field
