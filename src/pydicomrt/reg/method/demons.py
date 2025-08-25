from typing import Optional, Sequence, Tuple, Union, List
import numpy as np
import SimpleITK as sitk


def smooth_and_resample(
    image,
    isotropic_voxel_size_mm=None,
    shrink_factor=None,
    smoothing_sigma=None,
    interpolator=sitk.sitkLinear,
):
    """
    Smooth (optional Gaussian) and resample an image to a lower resolution or
    isotropic voxel size.

    This function is typically used to build an image pyramid for multi-resolution
    registration. It first applies Gaussian smoothing (if requested), then resamples
    the image using either:
      - a shrink factor (downsampling), or
      - a target isotropic voxel size.

    Parameters
    ----------
    image : sitk.Image
        Input image to smooth and resample.
    isotropic_voxel_size_mm : float, optional
        Desired isotropic voxel size (mm). If provided, overrides `shrink_factor`.
    shrink_factor : float or list of floats, optional
        Downsampling factor(s). If scalar, applied equally to all dimensions.
        If list/tuple, must match the number of image dimensions.
        Mutually exclusive with `isotropic_voxel_size_mm`.
    smoothing_sigma : float or list of floats, optional
        Gaussian smoothing sigma(s), in physical units (mm). If scalar, same
        sigma is applied in all dimensions; if sequence, must match the number
        of dimensions.
    interpolator : int, default = sitk.sitkLinear
        Interpolator enum used by SimpleITK's Resample function.

    Returns
    -------
    sitk.Image
        Smoothed and resampled image.

    Raises
    ------
    AttributeError
        If both `isotropic_voxel_size_mm` and `shrink_factor` are specified.
    """

    # ---- Step 1. Optional Gaussian smoothing ----
    if smoothing_sigma:
        if hasattr(smoothing_sigma, "__iter__"):
            # Variance = sigma^2 per dimension
            smoothing_variance = [i * i for i in smoothing_sigma]
        else:
            smoothing_variance = (smoothing_sigma ** 2,) * image.GetDimension()

        # Kernel width = ~8*sigma/spacing (rounded)
        maximum_kernel_width = int(
            max([8 * j * i for i, j in zip(image.GetSpacing(), smoothing_variance)])
        )

        # Apply smoothing in physical space
        image = sitk.DiscreteGaussian(image, smoothing_variance, maximum_kernel_width)

    # ---- Step 2. Retrieve current metadata ----
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()

    # ---- Step 3. Check for conflicting arguments ----
    if shrink_factor and isotropic_voxel_size_mm:
        raise AttributeError(
            "Function must be called with either isotropic_voxel_size_mm or "
            "shrink_factor, not both."
        )

    # ---- Step 4. Compute new size ----
    if isotropic_voxel_size_mm:
        # Target isotropic resolution → scale factor = target/original spacing
        scale_factor = (
            isotropic_voxel_size_mm * np.ones_like(image.GetSize()) / np.array(image.GetSpacing())
        )
        # Compute new size = old_size / scale_factor (rounded)
        new_size = [int(sz / float(sf) + 0.5) for sz, sf in zip(original_size, scale_factor)]

    elif shrink_factor:
        if isinstance(shrink_factor, (list, tuple)):
            # Per-dimension shrink factor
            new_size = [int(sz / float(sf) + 0.5) for sz, sf in zip(original_size, shrink_factor)]
        else:
            # Same shrink factor in all dimensions
            new_size = [int(sz / float(shrink_factor) + 0.5) for sz in original_size]

    else:
        # Neither shrink nor isotropic resampling → return unchanged
        return image

    # ---- Step 5. Compute new spacing from new size ----
    # Keep same physical extent → spacing = (extent / (new_size-1))
    new_spacing = [
        ((size_o_i - 1) * spacing_o_i) / (size_n_i - 1)
        for size_o_i, spacing_o_i, size_n_i in zip(original_size, original_spacing, new_size)
    ]

    # ---- Step 6. Resample with new size/spacing ----
    return sitk.Resample(
        image,
        new_size,
        sitk.Transform(),          # identity transform
        interpolator,              # chosen interpolator
        image.GetOrigin(),         # preserve origin
        new_spacing,               # computed spacing
        image.GetDirection(),      # preserve direction cosines
        0.0,                       # default background value
        image.GetPixelID(),        # preserve pixel type
    )


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
    smoothing_sigmas : Sequence[float], default = None
        Per-level Gaussian smoothing sigmas (in physical units). Length must
        match `resolution_staging`.
    iteration_staging : Sequence[int], default = None
        Number of iterations to run at each pyramid level.
    interp_order : int, default = sitk.sitkLinear
        Interpolator used for image resampling.

    Returns
    -------
    dvf_total : sitk.Image
        Final displacement vector field (DVF) in the fixed image space.

    Notes
    -----
    - Uses Gaussian smoothing at each scale to stabilize registration.
    - The final DVF is smoothed at each iteration and accumulated across levels.
    - Pixel type of DVF is constrained to sitkVectorFloat64 for compatibility
      with Demons filters in SimpleITK.
    """
    # Lists to hold multi-resolution fixed and moving images
    fixed_images = []
    moving_images = []

    # Construct the image pyramid: apply smoothing + resampling per scale
    for resolution, smoothing_sigma in zip(resolution_staging, smoothing_sigmas):
        isotropic_voxel_size_mm = None
        shrink_factor = None

        # Choose between isotropic voxel resampling or shrink factor
        if isotropic_resample:
            isotropic_voxel_size_mm = resolution
        else:
            shrink_factor = resolution

        # Resample fixed image
        fixed_images.append(
            smooth_and_resample(
                fixed_image,
                isotropic_voxel_size_mm=isotropic_voxel_size_mm,
                shrink_factor=shrink_factor,
                smoothing_sigma=smoothing_sigma,
                interpolator=interp_order,
            )
        )

        # Resample moving image
        moving_images.append(
            smooth_and_resample(
                moving_image,
                isotropic_voxel_size_mm=isotropic_voxel_size_mm,
                shrink_factor=shrink_factor,
                smoothing_sigma=smoothing_sigma,
                interpolator=interp_order,
            )
        )

    # Initialize displacement field
    if not initial_displacement_field:
        if initial_transform:
            # Convert transform into a displacement field image
            initial_displacement_field = sitk.TransformToDisplacementField(
                initial_transform,
                sitk.sitkVectorFloat64,
                fixed_image.GetSize(),
                fixed_image.GetOrigin(),
                fixed_image.GetSpacing(),
                fixed_image.GetDirection(),
            )
        else:
            # Create an empty displacement field image matching fixed_image size
            if len(moving_image.GetSize()) == 2:
                initial_displacement_field = sitk.Image(
                    fixed_image.GetWidth(),
                    fixed_image.GetHeight(),
                    sitk.sitkVectorFloat64,
                )
            elif len(moving_image.GetSize()) == 3:
                initial_displacement_field = sitk.Image(
                    fixed_image.GetWidth(),
                    fixed_image.GetHeight(),
                    fixed_image.GetDepth(),
                    sitk.sitkVectorFloat64,
                )
            # Copy metadata (origin, spacing, direction) from fixed_image
            initial_displacement_field.CopyInformation(fixed_image)
    else:
        # Resample given initial DVF to fixed image grid
        initial_displacement_field = sitk.Resample(initial_displacement_field, fixed_image)

    # Initialize total deformation vector field (DVF) aligned to fixed image
    dvf_total = sitk.Resample(initial_displacement_field, fixed_image)

    # Multi-scale loop: process from coarsest → finest resolution
    for i in range(len(fixed_images)):
        f_image = fixed_images[i]
        m_image = moving_images[i]

        # Resample DVF to current resolution
        dvf_total = sitk.Resample(dvf_total, f_image)

        # Convert DVF into a transform and warp the moving image
        tfm_total = sitk.DisplacementFieldTransform(sitk.Cast(dvf_total, sitk.sitkVectorFloat64))
        m_image = sitk.Resample(m_image, tfm_total, interp_order)

        # Configure number of iterations at this scale
        iters = iteration_staging[i]
        registration_algorithm.SetNumberOfIterations(iters)

        # Run demons registration at current scale
        dvf_iter = registration_algorithm.Execute(f_image, m_image)

        # Accumulate update into running DVF (compose transforms)
        dvf_total = dvf_total + sitk.Resample(dvf_iter, tfm_total)

        # Regularize/smooth the DVF
        sigma = registration_algorithm.GetStandardDeviations()
        dvf_total = sitk.SmoothingRecursiveGaussian(dvf_total, sigma)
        dvf_total = sitk.Cast(dvf_total, sitk.sitkVectorFloat64)

    # Resample final DVF to fixed image grid
    dvf_total = sitk.Resample(dvf_total, fixed_image)

    return dvf_total


def deformable_registration_command_iteration(method):
    """
    Utility function to print information during demons registration
    """
    print("{0:3} = {1:10.5f}".format(method.GetElapsedIterations(), method.GetMetric()))


def demons_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    resolution_staging: Sequence[int] = (8, 4, 1),
    iteration_staging: Sequence[int] = (10, 10, 10),
    isotropic_resample: bool = False,
    initial_displacement_field: Optional[sitk.Image] = None,
    regularisation_kernel_mm: Union[float, Sequence[float]] = 1.5,
    smoothing_sigma_factor: float = 1.0,
    smoothing_sigmas: Optional[Union[float, Sequence[float], bool]] = False,
    default_value: Optional[float] = None,
    ncores: int = 1,
    interp_order: int = sitk.sitkLinear,
    verbose: bool = False,
) -> Tuple[sitk.Image, sitk.Transform, sitk.Image]:
    """
    Deformable registration via Fast Symmetric-Forces Demons (SimpleITK) with multi-resolution strategy.

    The function registers `moving_image` to `fixed_image` by estimating a dense displacement
    field using a coarse-to-fine (image pyramid) schedule. It then applies the resulting
    displacement-field transform to resample the moving image into the fixed image space.

    Parameters
    ----------
    fixed_image : sitk.Image
        Reference (target) image. Output is defined on this image grid.
    moving_image : sitk.Image
        Image to be warped into the fixed image space.
    resolution_staging : Sequence[int], default = (8, 4, 1)
        Pyramid schedule. If `isotropic_resample=False`, these are downsampling factors.
        If `isotropic_resample=True`, these are target isotropic voxel sizes in millimeters.
        Length of this sequence defines the number of pyramid levels.
    iteration_staging : Sequence[int], default = (20, 20, 20)
        Number of iterations per pyramid level. Must match the length of `resolution_staging`.
    isotropic_resample : bool, default = False
        If True, interpret `resolution_staging` as isotropic voxel sizes (mm); otherwise as shrink factors.
    initial_displacement_field : Optional[sitk.Image], default = None
        Initial DVF (vector image). If provided, it overrides any `initial_transform` usage
        inside `multiscale_demons` and is resampled onto the fixed grid when needed.
    regularisation_kernel_mm : float or Sequence[float], default = 1.0
        Standard deviation(s) in millimeters used by demons for smoothing the update and DVF.
        Scalar broadcasts to all axes; sequence length must be 1 or image dimension.
    smoothing_sigma_factor : float, default = 1.0
        Fallback factor to derive `smoothing_sigmas` per level if `smoothing_sigmas` is not provided.
        Each level's sigma = `resolution_staging[level] * smoothing_sigma_factor`.
    smoothing_sigmas : float or Sequence[float] or bool, default = False
        Per-level Gaussian sigmas (in physical units) used prior to resampling at each level.
        - False/None: derived from `resolution_staging` and `smoothing_sigma_factor`
        - float: same sigma for all levels
        - sequence: length must match number of levels
    default_value : Optional[float], default = None
        Out-of-bounds value during final resampling. If None, use 0. For CT-like images
        (min <= -1000), default to -1000.
    ncores : int, default = 1
        Number of threads for demons filter.
    interp_order : int, default = sitk.sitkLinear
        Interpolator for resampling (e.g., sitk.sitkNearestNeighbor, sitk.sitkLinear, sitk.sitkBSpline).
    verbose : bool, default = False
        If True and a global callback `deformable_registration_command_iteration` is defined,
        prints per-iteration metric values.

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
    RuntimeError
        If helper `multiscale_demons(...)` is not available in the current scope.

    Notes
    -----
    - Internally casts images to Float32 for demons, then casts result back to the original
      moving pixel type after resampling.
    - Regularization sigmas are specified in *voxels* to SimpleITK; this function converts
      millimeter inputs to voxel units using the fixed image spacing.
    """

    # ---- Basic validation on lists/levels/cores ----
    if not isinstance(resolution_staging, Sequence) or len(resolution_staging) == 0:
        raise ValueError("`resolution_staging` must be a non-empty sequence of positive integers.")
    if not isinstance(iteration_staging, Sequence) or len(iteration_staging) == 0:
        raise ValueError("`iteration_staging` must be a non-empty sequence of positive integers.")
    if len(iteration_staging) != len(resolution_staging):
        raise ValueError("`iteration_staging` length must match `resolution_staging` length.")
    if any(int(x) <= 0 for x in resolution_staging):
        raise ValueError("All values in `resolution_staging` must be positive integers.")
    if any(int(x) <= 0 for x in iteration_staging):
        raise ValueError("All values in `iteration_staging` must be positive integers.")
    if ncores <= 0:
        raise ValueError("`ncores` must be a positive integer.")

    # ---- Dimensionality checks ----
    dim = fixed_image.GetDimension()
    if moving_image.GetDimension() != dim:
        raise ValueError("`fixed_image` and `moving_image` must have the same dimension.")
    if initial_displacement_field is not None and initial_displacement_field.GetDimension() != dim:
        raise ValueError("`initial_displacement_field` must have the same dimension as the images.")

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
    if isinstance(regularisation_kernel_mm, (int, float)):
        kernel_mm = np.full(dim, float(regularisation_kernel_mm), dtype=np.float64)
    else:
        ker = np.asarray(regularisation_kernel_mm, dtype=np.float64).ravel()
        if ker.size not in (1, dim):
            raise ValueError(
                f"`regularisation_kernel_mm` must be a float or a sequence of length 1 or {dim}."
            )
        kernel_mm = np.full(dim, float(ker[0]), dtype=np.float64) if ker.size == 1 else ker

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
    try:
        deformation_field: sitk.Image = multiscale_demons(
            registration_algorithm=registration_method,
            fixed_image=fixed_image,
            moving_image=moving_image,
            resolution_staging=tuple(int(x) for x in resolution_staging),
            smoothing_sigmas=smoothing_sigmas_list,
            iteration_staging=tuple(int(x) for x in iteration_staging),
            isotropic_resample=bool(isotropic_resample),
            initial_displacement_field=initial_displacement_field,
            interp_order=int(interp_order),
        )
    except NameError as e:
        # Provide a clearer error if helper function is missing
        raise RuntimeError(
            "Required helper `multiscale_demons(...)` is not defined in the current scope."
        ) from e

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
