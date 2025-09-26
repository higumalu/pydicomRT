import numpy as np
import SimpleITK as sitk


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

def registration_command_iteration(method):
    """
    Utility function to print information during (rigid, similarity, translation, B-splines)
    registration
    """
    print("{0:3} = {1:10.5f}".format(method.GetOptimizerIteration(), method.GetMetricValue()))


def stage_iteration(method):
    """
    Utility function to print information during stage change in registration
    """
    print(f"Number of parameters = {method.GetInitialTransform().GetNumberOfParameters()}")


def deformable_registration_command_iteration(method):
    """
    Utility function to print information during demons registration
    """
    print("{0:3} = {1:10.5f}".format(method.GetElapsedIterations(), method.GetMetric()))


def control_point_spacing_distance_to_number(image, grid_spacing):
    """
    Convert grid spacing specified in distance to number of control points
    """
    image_spacing = np.array(image.GetSpacing())
    image_size = np.array(image.GetSize())
    number_points = image_size * image_spacing / np.array(grid_spacing)
    return (number_points + 0.5).astype(int)


def smooth_and_resample(
    image,
    isotropic_voxel_size_mm=None,
    shrink_factor=None,
    smoothing_sigma=None,
    interpolator=sitk.sitkLinear,
):
    """
    Args:
        image (SimpleITK.Image): The image we want to resample.
        isotropic_voxel_size_mm (float | None): New voxel size in millimetres
        shrink_factor (list | float): A number greater than one, such that the new image's size is
            original_size/shrink_factor. Can also be specified independently for each
            dimension (sagittal, coronal, axial).
        smoothing_sigma (list | float): Scale for Gaussian smoothing, this is in physical
            (image spacing) units, not pixels. Can also be specified independently for
            each dimension (sagittal, coronal, axial).
    Return:
        SimpleITK.Image: Image which is a result of smoothing the input and then resampling
        it using the specified Gaussian kernel and shrink factor.
    """
    if smoothing_sigma:
        if hasattr(smoothing_sigma, "__iter__"):
            smoothing_variance = [i * i for i in smoothing_sigma]
        else:
            smoothing_variance = (smoothing_sigma ** 2,) * 3

        maximum_kernel_width = int(
            max([8 * j * i for i, j in zip(image.GetSpacing(), smoothing_variance)])
        )

        image = sitk.DiscreteGaussian(image, smoothing_variance, maximum_kernel_width)

    original_spacing = image.GetSpacing()
    original_size = image.GetSize()

    if shrink_factor and isotropic_voxel_size_mm:
        raise AttributeError(
            "Function must be called with either isotropic_voxel_size_mm or "
            "shrink_factor, not both."
        )

    elif isotropic_voxel_size_mm:
        scale_factor = (
            isotropic_voxel_size_mm * np.ones_like(image.GetSize()) / np.array(image.GetSpacing())
        )
        new_size = [int(sz / float(sf) + 0.5) for sz, sf in zip(original_size, scale_factor)]

    elif shrink_factor:
        if isinstance(shrink_factor, list):
            new_size = [int(sz / float(sf) + 0.5) for sz, sf in zip(original_size, shrink_factor)]
        else:
            new_size = [int(sz / float(shrink_factor) + 0.5) for sz in original_size]

    else:
        return image

    new_spacing = [
        ((size_o_i - 1) * spacing_o_i) / (size_n_i - 1)
        for size_o_i, spacing_o_i, size_n_i in zip(original_size, original_spacing, new_size)
    ]

    return sitk.Resample(
        image,
        new_size,
        sitk.Transform(),
        interpolator,
        image.GetOrigin(),
        new_spacing,
        image.GetDirection(),
        0.0,
        image.GetPixelID(),
    )


def bspline_registration(
    fixed_image,
    moving_image,
    fixed_structure=False,
    moving_structure=False,
    resolution_staging=[8, 4, 2],
    smooth_sigmas=[4, 2, 1],
    sampling_rate=0.1,
    optimiser="LBFGS",
    metric="mean_squares",
    initial_grid_spacing=64,
    grid_scale_factors=[1, 2, 4],
    interp_order=sitk.sitkBSpline,
    default_value=-1000,
    number_of_iterations=20,
    isotropic_resample=False,
    initial_isotropic_size=1,
    number_of_histogram_bins_mi=30,
    verbose=False,
    ncores=8,
):
    """
    B-Spline image registration using ITK

    IMPORTANT - THIS IS UNDER ACTIVE DEVELOPMENT

    Args:
        fixed_image ([SimpleITK.Image]): The fixed (target/primary) image.
        moving_image ([SimpleITK.Image]): The moving (secondary) image.
        fixed_structure (bool, optional): If defined, a binary SimpleITK.Image used to mask metric
                                          evaluation for the moving image. Defaults to False.
        moving_structure (bool, optional): If defined, a binary SimpleITK.Image used to mask metric
                                           evaluation for the fixed image. Defaults to False.
        resolution_staging (list, optional): The multi-resolution downsampling factors.
                                             Defaults to [8, 4, 2].
        smooth_sigmas (list, optional): The multi-resolution smoothing kernel scale (Gaussian).
                                        Defaults to [4, 2, 1].
        sampling_rate (float, optional): The fraction of voxels sampled during each iteration.
                                         Defaults to 0.1.
        optimiser (str, optional): The optimiser algorithm used for image registration.
                                   Available options:
                                    - LBFSGS
                                      (limited-memory Broyden–Fletcher–Goldfarb–Shanno (bounded).)
                                    - LBFSG
                                      (limited-memory Broyden–Fletcher–Goldfarb–Shanno
                                      (unbounded).)
                                    - CGLS (conjugate gradient line search)
                                    - gradient_descent
                                    - gradient_descent_line_search
                                   Defaults to "LBFGS".
        metric (str, optional): The metric to be optimised during image registration.
                                Available options:
                                 - correlation
                                 - mean_squares
                                 - demons
                                 - mutual_information
                                   (used with parameter number_of_histogram_bins_mi)
                                Defaults to "mean_squares".
        initial_grid_spacing (int, optional): Grid spacing of lower resolution stage (in mm).
                                              Defaults to 64.
        grid_scale_factors (list, optional): Factors to determine grid spacing at each
                                             multiresolution stage.
                                             Defaults to [1, 2, 4].
        interp_order (int, optional): Interpolation order of final resampling.
                                      Defaults to sitk.sitkBSpline (cubic).
        default_value (int, optional): Default image value. Defaults to -1000.
        number_of_iterations (int, optional): Number of iterations at each resolution stage.
                                              Defaults to 20.
        isotropic_resample (bool, optional): Flag whether to resample to isotropic resampling
                                             prior to registration.
                                             Defaults to False.
        initial_isotropic_size (int, optional): Voxel size (in mm) of resampled isotropic image
                                                (if used). Defaults to 1.
        number_of_histogram_bins_mi (int, optional): Number of histogram bins used when calculating
                                                     mutual information. Defaults to 30.
        verbose (bool, optional): Print image registration process information. Defaults to False.
        ncores (int, optional): Number of CPU cores used. Defaults to 8.

    Returns:
        [SimpleITK.Image]: The registered moving (secondary) image.
        [SimleITK.Transform]: The linear transformation.

    Notes:
     - smooth_sigmas are relative to resolution staging
        e.g. for image spacing of 1x1x1 mm^3, with smooth sigma=2 and resolution_staging=4, the
        scale of the Gaussian filter would be 2x4 = 8mm (i.e. 8x8x8 mm^3)

    """

    # Re-cast input images
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)

    moving_image_type = moving_image.GetPixelID()
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # (Optional) isotropic resample
    # This changes the behaviour, so care should be taken
    # For highly anisotropic images may be preferable

    if isotropic_resample:
        # First, copy the fixed image so we can resample back into this space at the end
        fixed_image_original = fixed_image
        fixed_image_original.MakeUnique()

        fixed_image = smooth_and_resample(
            fixed_image,
            isotropic_voxel_size_mm=initial_isotropic_size,
        )

        moving_image = smooth_and_resample(
            moving_image,
            isotropic_voxel_size_mm=initial_isotropic_size,
        )

    else:
        fixed_image_original = fixed_image

    # Set up image registration method
    registration = sitk.ImageRegistrationMethod()
    registration.SetNumberOfThreads(ncores)

    registration.SetShrinkFactorsPerLevel(resolution_staging)
    registration.SetSmoothingSigmasPerLevel(smooth_sigmas)
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Choose optimiser
    if optimiser.lower() == "lbfgsb":
        registration.SetOptimizerAsLBFGSB(
            gradientConvergenceTolerance=1e-5,
            numberOfIterations=number_of_iterations,
            maximumNumberOfCorrections=5,
            maximumNumberOfFunctionEvaluations=1024,
            costFunctionConvergenceFactor=1e7,
            trace=verbose,
        )
    elif optimiser.lower() == "lbfgs":
        registration.SetOptimizerAsLBFGS2(
            numberOfIterations=number_of_iterations,
            solutionAccuracy=1e-2,
            hessianApproximateAccuracy=6,
            deltaConvergenceDistance=0,
            deltaConvergenceTolerance=0.01,
            lineSearchMaximumEvaluations=40,
            lineSearchMinimumStep=1e-20,
            lineSearchMaximumStep=1e20,
            lineSearchAccuracy=0.01,
        )
    elif optimiser.lower() == "cgls":
        registration.SetOptimizerAsConjugateGradientLineSearch(
            learningRate=0.05, numberOfIterations=number_of_iterations
        )
        registration.SetOptimizerScalesFromPhysicalShift()
    elif optimiser.lower() == "gradient_descent":
        registration.SetOptimizerAsGradientDescent(
            learningRate=5.0,
            numberOfIterations=number_of_iterations,
            convergenceMinimumValue=1e-6,
            convergenceWindowSize=10,
        )
        registration.SetOptimizerScalesFromPhysicalShift()
    elif optimiser.lower() == "gradient_descent_line_search":
        registration.SetOptimizerAsGradientDescentLineSearch(
            learningRate=1.0, numberOfIterations=number_of_iterations
        )
        registration.SetOptimizerScalesFromPhysicalShift()

    # Set metric
    if metric == "correlation":
        registration.SetMetricAsCorrelation()
    elif metric == "mean_squares":
        registration.SetMetricAsMeanSquares()
    elif metric == "demons":
        registration.SetMetricAsDemons()
    elif metric == "mutual_information":
        registration.SetMetricAsMattesMutualInformation(
            numberOfHistogramBins=number_of_histogram_bins_mi
        )

    registration.SetInterpolator(sitk.sitkLinear)

    # Set sampling
    if isinstance(sampling_rate, float):
        registration.SetMetricSamplingPercentage(sampling_rate)
    elif type(sampling_rate) in [np.ndarray, list]:
        registration.SetMetricSamplingPercentagePerLevel(sampling_rate)

    registration.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.REGULAR)

    # Set masks
    if moving_structure is not False:
        registration.SetMetricMovingMask(moving_structure)

    if fixed_structure is not False:
        registration.SetMetricFixedMask(fixed_structure)

    # Set control point spacing
    transform_domain_mesh_size = control_point_spacing_distance_to_number(
        fixed_image, initial_grid_spacing
    )

    if verbose:
        print(f"Initial grid size: {transform_domain_mesh_size}")

    # Initialise transform
    initial_transform = sitk.BSplineTransformInitializer(
        fixed_image,
        transformDomainMeshSize=[int(i) for i in transform_domain_mesh_size],
    )
    registration.SetInitialTransformAsBSpline(
        initial_transform, inPlace=True, scaleFactors=grid_scale_factors
    )

    # (Optionally) add iteration commands
    if verbose:
        registration.AddCommand(
            sitk.sitkIterationEvent,
            lambda: registration_command_iteration(registration),
        )
        registration.AddCommand(
            sitk.sitkMultiResolutionIterationEvent,
            lambda: stage_iteration(registration),
        )

    # Run the registration
    output_transform = registration.Execute(fixed=fixed_image, moving=moving_image)

    # Resample moving image
    registered_image = apply_transform(
        input_image=moving_image,
        reference_image=fixed_image_original,
        transform=output_transform,
        default_value=default_value,
        interpolator=interp_order,
    )

    registered_image = sitk.Cast(registered_image, moving_image_type)

    dvf_image = sitk.TransformToDisplacementField(
        output_transform,
        sitk.sitkVectorFloat32,
        fixed_image_original.GetSize(),
        fixed_image_original.GetOrigin(),
        fixed_image_original.GetSpacing(),
        fixed_image_original.GetDirection(),
    )
    output_transform = sitk.DisplacementFieldTransform(dvf_image)

    # Return outputs
    return registered_image, output_transform
