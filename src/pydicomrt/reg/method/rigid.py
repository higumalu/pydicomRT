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
