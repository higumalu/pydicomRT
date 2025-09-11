import SimpleITK as sitk


def rigid_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    histogram_bins: int = 100,
    learning_rate: float = 2.0,
    iterations: int = 100,
    convergence_minimum_value: float = 1e-6,
    convergence_window_size: int = 10,
    ) -> sitk.Transform:
    """
    Perform rigid registration between two 3D images using SimpleITK.

    This function aligns a moving image to a fixed image by estimating
    a rigid transformation (translation + rotation, no scaling/shearing).
    The method uses Mattes Mutual Information as a similarity metric,
    a gradient descent optimizer, and returns an affine transform
    containing only the rigid components.

    Args:
        fixed_image (sitk.Image): The reference image to which the moving image is aligned.
        moving_image (sitk.Image): The image to be registered (transformed).
        histogram_bins (int, optional): Number of histogram bins for Mattes Mutual Information metric. Default = 50.
        learning_rate (float, optional): Step size for gradient descent optimizer. Default = 1.0.
        iterations (int, optional): Maximum number of optimizer iterations. Default = 100.
        convergence_minimum_value (float, optional): Minimum convergence value for optimizer stopping criterion. Default = 1e-6.
        convergence_window_size (int, optional): Window size for convergence checking. Default = 10.

    Returns:
        sitk.Transform: An affine transform containing only rotation and translation
                        that aligns the moving image to the fixed image.
    """

    # Ensure both images are float32 for numerical stability in registration
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # Create the registration method object
    registration_method = sitk.ImageRegistrationMethod()

    # Use Mattes Mutual Information (robust for multimodal image registration)
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=histogram_bins)

    # Use linear interpolation for resampling the moving image
    registration_method.SetInterpolator(sitk.sitkLinear)

    # Configure the optimizer as gradient descent with given parameters
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=learning_rate,
        numberOfIterations=iterations,
        convergenceMinimumValue=convergence_minimum_value,
        convergenceWindowSize=convergence_window_size,
    )

    # Scale optimizer step sizes according to physical units of the image
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Initialize with a centered rigid transform (rotation + translation)
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image,
        moving_image,
        sitk.VersorRigid3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )

    # Explicitly set the center to (0,0,0) – prevents unwanted shifting
    initial_transform.SetCenter((0.0, 0.0, 0.0))

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
    affine_transform = sitk.AffineTransform(3)
    affine_transform.SetMatrix(transform.GetMatrix())       # Rotation part
    affine_transform.SetTranslation(transform.GetTranslation())  # Translation part

    return affine_transform


def alignment_registration(fixed_image, moving_image, moments=True, interpolator=sitk.sitkLinear, default_value=0):
    """
    A simple registration procedure that can align images in a single step.
    Uses the image centres-of-mass (and optionally second moments) to
    estimate the shift (and rotation) needed for alignment.

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

    Returns:
        [SimpleITK.Image]: The registered moving (secondary) image.
        [SimleITK.Transform]: The linear transformation.
    """

    moving_image_type = moving_image.GetPixelIDValue()
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image, moving_image, sitk.VersorRigid3DTransform(), moments
    )
    aligned_image = sitk.Resample(moving_image, fixed_image, initial_transform, interpolator, default_value)
    aligned_image = sitk.Cast(aligned_image, moving_image_type)
    return aligned_image, initial_transform


def registration_command_iteration(method):
    """
    Utility function to print information during (rigid, similarity, translation, B-splines)
    registration
    """
    print("{0:3} = {1:10.5f}".format(method.GetOptimizerIteration(), method.GetMetricValue()))


def apply_transform(
    input_image,
    reference_image=None,
    transform=None,
    default_value=0,
    interpolator=sitk.sitkNearestNeighbor,
):
    """
    Transform a volume of structure with the given deformation field.

    Args
        input_image (SimpleITK.Image): The image, to which the transform is applied
        reference_image (SimpleITK.Image): The image will be resampled into this reference space.
        transform (SimpleITK.Transform): The transformation
        default_value: Default (background) value. Defaults to 0.
        interpolator (int, optional): The interpolation order.
                                Available options:
                                    - SimpleITK.sitkNearestNeighbor
                                    - SimpleITK.sitkLinear
                                    - SimpleITK.sitkBSpline
                                Defaults to SimpleITK.sitkNearestNeighbor

    Returns
        (SimpleITK.Image): the transformed image

    """
    original_image_type = input_image.GetPixelID()

    resampler = sitk.ResampleImageFilter()

    if reference_image:
        resampler.SetReferenceImage(reference_image)
    else:
        resampler.SetReferenceImage(input_image)

    if transform:
        resampler.SetTransform(transform)

    resampler.SetDefaultPixelValue(default_value)
    resampler.SetInterpolator(interpolator)

    output_image = resampler.Execute(input_image)
    output_image = sitk.Cast(output_image, original_image_type)

    return output_image


def linear_registration(
    fixed_image,
    moving_image,
    fixed_structure=None,
    moving_structure=None,
    reg_method="similarity",
    metric="mean_squares",
    optimiser="gradient_descent",
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
    There are several transforms available, with options for the metric and optimiser to be used.
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
        optimiser (str, optional): The optimiser algorithm used for image registration.
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

    if optimiser.lower() == "lbfgsb":
        registration.SetOptimizerAsLBFGSB(
            gradientConvergenceTolerance=1e-5,
            numberOfIterations=number_of_iterations,
            maximumNumberOfCorrections=50,
            maximumNumberOfFunctionEvaluations=1024,
            costFunctionConvergenceFactor=1e7,
            trace=verbose,
        )
    elif optimiser.lower() == "exhaustive":
        """
        This isn't well implemented
        Needs some work to give options for sampling rates
        Use is not currently recommended
        """
        samples = [10, 10, 10, 10, 10, 10]
        registration.SetOptimizerAsExhaustive(samples)
    elif optimiser.lower() == "gradient_descent_line_search":
        registration.SetOptimizerAsGradientDescentLineSearch(
            learningRate=1.0, numberOfIterations=number_of_iterations
        )
    elif optimiser.lower() == "gradient_descent":
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