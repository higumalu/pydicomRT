import logging
import numpy as np
import SimpleITK as sitk
from .common import apply_transform, registration_command_iteration

logger = logging.getLogger(__name__)


def to_centre_free_affine(transform) -> sitk.AffineTransform:
    """
    Rewrite a centred rigid transform as an equivalent affine with its centre at the origin.

    A centred transform maps ``y = R (x - C) + C + T``. Copying only ``R`` and ``T`` into a
    centre-free affine silently drops the ``(I - R) C`` term, which is zero only while the
    rotation is identity -- so the error appears exactly once the optimizer starts rotating.
    Folding the centre into the translation keeps the mapping identical, which matters
    because the DICOM REG matrix has nowhere to record a rotation centre.

    Args:
        transform: A SimpleITK transform exposing GetMatrix/GetTranslation/GetCenter.

    Returns:
        sitk.AffineTransform: Transform with centre (0, 0, 0) and the same point mapping.
    """
    rotation = np.array(transform.GetMatrix()).reshape(3, 3)
    centre = np.array(transform.GetCenter())
    translation = np.array(transform.GetTranslation())

    offset = centre + translation - rotation @ centre

    affine_transform = sitk.AffineTransform(3)
    affine_transform.SetMatrix(transform.GetMatrix())
    affine_transform.SetTranslation(offset.tolist())
    return affine_transform


def rigid_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    histogram_bins: int = 100,
    learning_rate: float = 2.0,
    iterations: int = 300,
    shrink_factors=(4, 2, 1),
    smoothing_sigmas=(2.0, 1.0, 0.0),
    optimizer: str = "regular_step",
    relaxation_factor: float = 0.7,
    min_step_mm: float = 1e-4,
    convergence_minimum_value: float = 1e-6,
    convergence_window_size: int = 10,
    max_step_size_mm: float = 2.0,
    algorithm: str = "stochastic",
    backend: str = "auto",
    seed: int = 0,
    threads=None,
    ) -> sitk.Transform:
    """
    Rigidly align a moving image to a fixed image.

    Estimates translation and rotation only -- no scaling, no shear -- by maximising
    Mattes mutual information over a multi-resolution pyramid. Mutual information rather
    than intensity difference, so CT-to-CBCT and CT-to-MR work without matching intensity
    scales.

    Parameters
    ----------
    fixed_image : sitk.Image
        The reference. The result is expressed on this image's grid.
    moving_image : sitk.Image
        The image to align.
    histogram_bins : int, optional
        itk only. Bins for Mattes mutual information. Default 100.
    learning_rate : float, optional
        itk only. Initial optimizer step. Default 2.0.
    iterations : int, optional
        itk only. Maximum optimizer iterations per resolution level. Default 300.
    shrink_factors : sequence of int or None, optional
        itk only. Downsampling factor per level, coarsest first. ``None`` gives a
        single-resolution registration. Default ``(4, 2, 1)``.
    smoothing_sigmas : sequence of float or None, optional
        itk only. Gaussian sigma in mm per level, matching ``shrink_factors`` element for
        element. Default ``(2.0, 1.0, 0.0)``.
    optimizer : {"regular_step", "gradient_descent"}, optional
        itk only. Default ``"regular_step"``. See Notes.
    relaxation_factor : float, optional
        itk regular_step only. The step is multiplied by this whenever the gradient
        direction reverses. Default 0.7.
    min_step_mm : float, optional
        itk regular_step only. Stop once the step falls below this. Default 1e-4.
    convergence_minimum_value : float, optional
        itk gradient_descent only. Default 1e-6.
    convergence_window_size : int, optional
        itk gradient_descent only. Default 10.
    max_step_size_mm : float, optional
        itk gradient_descent only. Upper bound on a single step in mm; 0.0 leaves it
        unbounded. Default 2.0.
    algorithm : {"stochastic", "itk"}, optional
        ``"stochastic"`` (default) is the library's own implementation, with a fixed
        schedule -- see Notes. ``"itk"`` is SimpleITK's ``ImageRegistrationMethod``,
        configured by the parameters marked "itk only", which ``"stochastic"`` ignores.
    backend : {"auto", "numpy", "cpp"}, optional
        stochastic only. ``"cpp"`` runs the inner loops in the compiled core that ships
        with the package, about ten times faster; ``"numpy"`` runs them in numpy. Both use
        the same random points; rounding can still tip single points into another histogram
        bin, which a thousand stochastic steps carry forward, so their answers differ by up
        to about 0.2 mm. ``"auto"`` (default) uses the compiled core, or numpy in an
        installation built from source without a C++ compiler.
    seed : int, optional
        stochastic only. Seed of the random sample points. Default 0.
    threads : int or None, optional
        stochastic with the compiled core only. ``None`` (default) takes SimpleITK's global
        default thread count. The answer does not depend on it.

    Returns
    -------
    sitk.Transform
        An affine transform holding rotation and translation only, pointing **fixed to
        moving** -- a *resampling* transform, the direction ``sitk.Resample`` expects.
        Invert it before storing in a DICOM Spatial REG, which runs moving to fixed.
        Keep its direction when using it as a Deformable REG post-matrix.

    Raises
    ------
    ValueError
        If the algorithm, backend or optimizer name is unknown, an image is not 3D
        (stochastic), or ``shrink_factors`` and ``smoothing_sigmas`` differ in length (itk).
    ImportError
        If ``backend="cpp"`` and this installation has no compiled core.

    See Also
    --------
    registration_pipeline : Adds preprocessing and an optional deformable stage.
    demons_registration : Deformable refinement, to run after this.
    SpatialRegistrationBuilder : Export the result as DICOM.

    Notes
    -----
    ``"stochastic"`` follows elastix's strategy: Mattes mutual information on fresh random
    points every iteration, adaptive stochastic gradient ascent over a four-level smoothing
    pyramid, then a regular-step ascent on one fixed set of 500 000 points, so the answer is
    the optimum of a single estimate rather than wherever the stochastic walk stopped. It
    uses the exact derivative of the linear interpolant. Mean error in mm, elastix's default
    rigid in brackets: 1.38 (1.32) against the RIRE fiducial gold standard, 1.24 (1.35) from
    ten starts moved by up to 10 degrees and 20 mm, 0.53 (1.00) on nine SynthRAD2023 brain
    MR/CT pairs, and 5.1 (16.7) from the TPS registration of one pelvis CT/CBCT pair. Its
    answer does not depend on the thread count. It starts from the geometric centres
    aligned, as ``CenteredTransformInitializer(GEOMETRY)`` does.

    ``"itk"`` is the previous default. Its precomputed image gradient disagrees with the
    linear interpolant on thick slices, which walked it 30 mm away from a correct start on
    a 4 mm MR; on the pelvis pair it took about nine minutes single-threaded for 7.5 mm. It
    is only bit-reproducible single-threaded -- ITK reduces the metric in thread-completion
    order. Its ``optimizer="gradient_descent"`` can return an alignment worse than its own
    initialisation and still report convergence.

    Both register in the images' own grids: resampling both onto a common grid beforehand
    leaves no geometry to initialise from.

    Examples
    --------
    >>> transform = rigid_registration(fixed_image, moving_image)      # doctest: +SKIP
    >>> registered = sitk.Resample(                                    # doctest: +SKIP
    ...     moving_image, fixed_image, transform, sitk.sitkLinear, -1000.0)

    The previous SimpleITK algorithm:

    >>> transform = rigid_registration(                                # doctest: +SKIP
    ...     fixed_image, moving_image, algorithm="itk")
    """
    if algorithm not in ("itk", "stochastic"):
        raise ValueError(f"algorithm must be 'itk' or 'stochastic', got {algorithm!r}")
    if algorithm == "stochastic":
        if backend not in ("auto", "numpy", "cpp"):
            raise ValueError(f"backend must be 'auto', 'numpy' or 'cpp', got {backend!r}")
        from . import _stochastic_rigid

        if threads is None:
            threads = sitk.ProcessObject.GetGlobalDefaultNumberOfThreads()
        rotation, translation, centre = _stochastic_rigid.register(
            fixed_image, moving_image, seed=seed, backend=backend, threads=threads
        )
        centred = sitk.AffineTransform(3)
        centred.SetCenter(centre.tolist())
        centred.SetMatrix(rotation.ravel().tolist())
        centred.SetTranslation(translation.tolist())
        return to_centre_free_affine(centred)

    if optimizer not in ("regular_step", "gradient_descent"):
        raise ValueError(
            f"optimizer must be 'regular_step' or 'gradient_descent', got {optimizer!r}"
        )

    if shrink_factors is not None and smoothing_sigmas is not None:
        if len(shrink_factors) != len(smoothing_sigmas):
            raise ValueError(
                "shrink_factors and smoothing_sigmas must describe the same number of "
                f"resolution levels (got {len(shrink_factors)} and {len(smoothing_sigmas)})"
            )

    # Ensure both images are float32 for numerical stability in registration
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # Create the registration method object
    registration_method = sitk.ImageRegistrationMethod()

    # Use Mattes Mutual Information (robust for multimodal image registration)
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=histogram_bins)

    # Use linear interpolation for resampling the moving image
    registration_method.SetInterpolator(sitk.sitkLinear)

    # A step that shrinks on gradient reversal cannot walk away from a good alignment,
    # which plain gradient descent demonstrably does -- see the docstring note.
    if optimizer == "regular_step":
        registration_method.SetOptimizerAsRegularStepGradientDescent(
            learningRate=learning_rate,
            minStep=min_step_mm,
            numberOfIterations=iterations,
            relaxationFactor=relaxation_factor,
            gradientMagnitudeTolerance=1e-8,
            maximumStepSizeInPhysicalUnits=0.0,
        )
    else:
        registration_method.SetOptimizerAsGradientDescent(
            learningRate=learning_rate,
            numberOfIterations=iterations,
            convergenceMinimumValue=convergence_minimum_value,
            convergenceWindowSize=convergence_window_size,
            estimateLearningRate=registration_method.Once,
            maximumStepSizeInPhysicalUnits=max_step_size_mm,
        )

    # Scale optimizer step sizes according to physical units of the image
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Coarse-to-fine: the coarse levels carry the transform most of the way before the
    # fine level refines it, which is what gives the optimizer a usable capture range.
    if shrink_factors is not None:
        registration_method.SetShrinkFactorsPerLevel(list(shrink_factors))
        registration_method.SetSmoothingSigmasPerLevel(list(smoothing_sigmas))
        registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Initialize with a centered rigid transform (rotation + translation). The centre stays
    # at the image centre: rotating about a far-away world origin makes the versor
    # parameters badly conditioned, since a small angle becomes a large displacement.
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image,
        moving_image,
        sitk.VersorRigid3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )

    # Assign the initial transform to the registration method
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    # Run the registration (this optimizes the transform parameters)
    final_transform = registration_method.Execute(fixed_image, moving_image)

    # If the result is a composite transform, extract the first (rigid) transform
    if isinstance(final_transform, sitk.CompositeTransform):
        transform = final_transform.GetNthTransform(0)
    else:
        transform = final_transform

    # Convert rigid transform into an affine transform (matrix + translation only)
    return to_centre_free_affine(transform)


def alignment_registration(
    fixed_image, 
    moving_image, 
    moments=True, 
    interpolator=sitk.sitkLinear, 
    default_value=0,
    iterations=0,
    learning_rate=1.0,
    histogram_bins=50,
    convergence_minimum_value=1e-6,
    convergence_window_size=10
):
    """
    A simple registration procedure that can align images in a single step.
    Uses the image centres-of-mass (and optionally second moments) to
    estimate the shift (and rotation) needed for alignment.
    
    If iterations > 0, performs iterative optimization using gradient descent
    with Mattes Mutual Information metric.

    Args:
        fixed_image ([SimpleITK.Image]): The fixed (target/primary) image.
        moving_image ([SimpleITK.Image]): The moving (secondary) image.
        moments (bool, optional): Option to align images using the second moment. Defaults to True.
        interpolator (int, optional): The interpolation order.
                                Available options:
                                    - SimpleITK.sitkNearestNeighbor
                                    - SimpleITK.sitkLinear
                                    - SimpleITK.sitkBSpline
                                Defaults to SimpleITK.sitkLinear.
        default_value (int, optional): Default (background) value. Defaults to 0.
        iterations (int, optional): Number of optimization iterations. 
                                    If 0, only performs initialization without optimization (original behavior).
                                    Defaults to 0.
        learning_rate (float, optional): Step size for gradient descent optimizer. Defaults to 1.0.
        histogram_bins (int, optional): Number of histogram bins for Mattes Mutual Information metric. Defaults to 50.
        convergence_minimum_value (float, optional): Minimum convergence value for optimizer stopping criterion. Defaults to 1e-6.
        convergence_window_size (int, optional): Window size for convergence checking. Defaults to 10.

    Returns:
        [SimpleITK.Image]: The registered moving (secondary) image.
        [SimpleITK.Transform]: The linear transformation.
    """

    moving_image_type = moving_image.GetPixelIDValue()
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)
    
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image, moving_image, sitk.VersorRigid3DTransform(), moments
    )
    
    # If iterations > 0, perform iterative optimization
    if iterations > 0:
        # Create the registration method object
        registration_method = sitk.ImageRegistrationMethod()
        
        # Use Mattes Mutual Information (robust for multimodal image registration)
        registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=histogram_bins)
        
        # Use the specified interpolator
        registration_method.SetInterpolator(interpolator)
        
        # Configure the optimizer as gradient descent with given parameters
        registration_method.SetOptimizerAsGradientDescent(
            learningRate=learning_rate,
            numberOfIterations=iterations,
            convergenceMinimumValue=convergence_minimum_value,
            convergenceWindowSize=convergence_window_size,
        )
        
        # Scale optimizer step sizes according to physical units of the image
        registration_method.SetOptimizerScalesFromPhysicalShift()
        
        # Assign the initial transform to the registration method
        registration_method.SetInitialTransform(initial_transform, inPlace=False)
        
        # Run the registration (this optimizes the transform parameters)
        final_transform = registration_method.Execute(fixed_image, moving_image)
        
        # Use the optimized transform
        transform = final_transform
    else:
        # Original behavior: use initial transform without optimization
        transform = initial_transform
    
    # Resample the moving image using the transform
    aligned_image = sitk.Resample(moving_image, fixed_image, transform, interpolator, default_value)
    aligned_image = sitk.Cast(aligned_image, moving_image_type)
    
    return aligned_image, transform






def linear_registration(
    fixed_image,
    moving_image,
    fixed_structure=None,
    moving_structure=None,
    reg_method="similarity",
    metric="mean_squares",
    optimizer="gradient_descent",
    shrink_factors=[8, 2, 1],
    smooth_sigmas=[4, 2, 0],
    sampling_rate=0.25,
    final_interp=2,
    number_of_iterations=100,
    default_value=None,
    verbose=False,
):
    """
    Initial linear registration between two images.
    The images are not required to be in the same space.
    There are several transforms available, with options for the metric and optimizer to be used.
    Note the default_value, which should be set to match the image modality.

    Args:
        fixed_image ([SimpleITK.Image]): The fixed (target/primary) image.
        moving_image ([SimpleITK.Image]): The moving (secondary) image.
        fixed_structure (bool, optional): If defined, a binary SimpleITK.Image used to mask metric
                                          evaluation for the moving image. Defaults to False.
        moving_structure (bool, optional): If defined, a binary SimpleITK.Image used to mask metric
                                           evaluation for the fixed image. Defaults to False.
        reg_method (str, optional): The linear transformtation model to be used for image
                                    registration.
                                    Available options:
                                     - translation
                                     - rigid
                                     - similarity
                                     - affine
                                     - scale
                                     - scaleversor
                                     - scaleskewversor
                                    Defaults to "Similarity".
        metric (str, optional): The metric to be optimised during image registration.
                                Available options:
                                 - correlation
                                 - mean_squares
                                 - mattes_mi
                                 - joint_hist_mi
                                Defaults to "mean_squares".
        optimizer (str, optional): The optimizer algorithm used for image registration.
                                   Available options:
                                    - lbfgsb
                                      (limited-memory Broyden–Fletcher–Goldfarb–Shanno (bounded).)
                                    - gradient_descent
                                    - gradient_descent_line_search
                                   Defaults to "gradient_descent".
        shrink_factors (list, optional): The multi-resolution downsampling factors.
                                         Defaults to [8, 2, 1].
        smooth_sigmas (list, optional): The multi-resolution smoothing kernel scale (Gaussian).
                                        Defaults to [4, 2, 0].
        sampling_rate (float, optional): The fraction of voxels sampled during each iteration.
                                         Defaults to 0.25.
        ants_radius (int, optional): Used is the metric is set as "ants_radius". Defaults to 3.
        final_interp (int, optional): The final interpolation order. Defaults to 2 (linear).
        number_of_iterations (int, optional): Number of iterations in each multi-resolution step.
                                              Defaults to 50.
        default_value (int, optional): Default voxel value. Defaults to 0 unless image is CT-like.
        verbose (bool, optional): Print image registration process information. Defaults to False.

    Returns:
        [SimpleITK.Image]: The registered moving (secondary) image.
        [SimleITK.Transform]: The linear transformation.
    """

    # Re-cast
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)

    moving_image_type = moving_image.GetPixelIDValue()
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # Initialise using a VersorRigid3DTransform
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image, moving_image, sitk.Euler3DTransform(), False
    )
    # Set up image registration method
    registration = sitk.ImageRegistrationMethod()

    registration.SetShrinkFactorsPerLevel(shrink_factors)
    registration.SetSmoothingSigmasPerLevel(smooth_sigmas)
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    registration.SetMovingInitialTransform(initial_transform)

    if metric.lower() == "correlation":
        registration.SetMetricAsCorrelation()
    elif metric.lower() == "mean_squares":
        registration.SetMetricAsMeanSquares()
    elif metric.lower() == "mattes_mi":
        registration.SetMetricAsMattesMutualInformation()
    elif metric.lower() == "joint_hist_mi":
        registration.SetMetricAsJointHistogramMutualInformation()
    # to do: add the rest

    registration.SetInterpolator(sitk.sitkLinear)  # Perhaps a small gain in improvement
    registration.SetMetricSamplingPercentage(sampling_rate, seed=42)
    registration.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.REGULAR)

    # This is only necessary if using a transform comprising changes with different units
    # e.g. rigid (rotation: radians, translation: mm)
    # It can safely be left on
    registration.SetOptimizerScalesFromPhysicalShift()

    if moving_structure:
        registration.SetMetricMovingMask(moving_structure)

    if fixed_structure:
        registration.SetMetricFixedMask(fixed_structure)

    if isinstance(reg_method, str):
        if reg_method.lower() == "translation":
            registration.SetInitialTransform(sitk.TranslationTransform(3))
        elif reg_method.lower() == "similarity":
            registration.SetInitialTransform(sitk.Similarity3DTransform())
        elif reg_method.lower() == "affine":
            registration.SetInitialTransform(sitk.AffineTransform(3))
        elif reg_method.lower() == "rigid":
            registration.SetInitialTransform(sitk.VersorRigid3DTransform())
        elif reg_method.lower() == "scale":
            registration.SetInitialTransform(sitk.ScaleTransform(3))
        elif reg_method.lower() == "scaleversor":
            registration.SetInitialTransform(sitk.ScaleVersor3DTransform())
        elif reg_method.lower() == "scaleskewversor":
            registration.SetInitialTransform(sitk.ScaleSkewVersor3DTransform())
        else:
            raise ValueError(
                "You have selected a registration method that does not exist.\n Please select from"
                " Translation, Similarity, Affine, Rigid, ScaleVersor, ScaleSkewVersor"
            )
    elif isinstance(
        reg_method,
        (
            sitk.CompositeTransform,
            sitk.Transform,
            sitk.TranslationTransform,
            sitk.Similarity3DTransform,
            sitk.AffineTransform,
            sitk.VersorRigid3DTransform,
            sitk.ScaleVersor3DTransform,
            sitk.ScaleSkewVersor3DTransform,
        ),
    ):
        registration.SetInitialTransform(reg_method)
    else:
        raise ValueError(
            "'reg_method' must be either a string (see docs for acceptable registration names), "
            "or a custom sitk.CompositeTransform."
        )

    if optimizer.lower() == "lbfgsb":
        registration.SetOptimizerAsLBFGSB(
            gradientConvergenceTolerance=1e-5,
            numberOfIterations=number_of_iterations,
            maximumNumberOfCorrections=50,
            maximumNumberOfFunctionEvaluations=1024,
            costFunctionConvergenceFactor=1e7,
            trace=verbose,
        )
    elif optimizer.lower() == "exhaustive":
        """
        This isn't well implemented
        Needs some work to give options for sampling rates
        Use is not currently recommended
        """
        samples = [10, 10, 10, 10, 10, 10]
        registration.SetOptimizerAsExhaustive(samples)
    elif optimizer.lower() == "gradient_descent_line_search":
        registration.SetOptimizerAsGradientDescentLineSearch(
            learningRate=1.0, numberOfIterations=number_of_iterations
        )
    elif optimizer.lower() == "gradient_descent":
        registration.SetOptimizerAsGradientDescent(
            learningRate=2.0, numberOfIterations=number_of_iterations
        )

    if verbose:
        registration.AddCommand(
            sitk.sitkIterationEvent,
            lambda: registration_command_iteration(registration),
        )

    output_transform = registration.Execute(fixed=fixed_image, moving=moving_image)
    # Combine initial and optimised transform
    combined_transform = sitk.CompositeTransform([initial_transform, output_transform])

    # Try to find default value
    if default_value is None:
        default_value = 0

        # Test if image is CT-like
        if sitk.GetArrayViewFromImage(moving_image).min() <= -1000:
            default_value = -1000

    registered_image = apply_transform(
        input_image=moving_image,
        reference_image=fixed_image,
        transform=combined_transform,
        default_value=default_value,
        interpolator=final_interp,
    )

    registered_image = sitk.Cast(registered_image, moving_image_type)

    return registered_image, combined_transform
