import numpy as np
import SimpleITK as sitk


def rigid_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    histogram_bins: int = 50,
    learning_rate: float = 1.0,
    iterations: int = 100,
    convergence_minimum_value: float = 1e-6,
    convergence_window_size: int = 10,
    ) -> sitk.Transform:
    """
    Perform rigid registration (translation and rotation only) between two 3D SimpleITK images and return an affine transform.
    Args:
        fixed_image (sitk.Image)
        moving_image (sitk.Image)
        histogram_bins (int)
        learning_rate (float)
        iterations (int)
        convergence_minimum_value (float)
        convergence_window_size (int)
    Returns:
        sitk.Transform: The resulting affine transform (containing only rotation + translation)
    """
    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=histogram_bins)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=learning_rate,
        numberOfIterations=iterations,
        convergenceMinimumValue=convergence_minimum_value,
        convergenceWindowSize=convergence_window_size,
    )
    registration_method.SetOptimizerScalesFromPhysicalShift()
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image,
        moving_image,
        sitk.VersorRigid3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    initial_transform.SetCenter((0.0, 0.0, 0.0))
    registration_method.SetInitialTransform(initial_transform, inPlace=False)
    final_transform = registration_method.Execute(fixed_image, moving_image)
    if isinstance(final_transform, sitk.CompositeTransform):
        transform = final_transform.GetNthTransform(0)
    else:
        transform = final_transform
    affine_transform = sitk.AffineTransform(3)
    affine_transform.SetMatrix(transform.GetMatrix())
    affine_transform.SetTranslation(transform.GetTranslation())
    return affine_transform


def deformable_registration(fixed_image, moving_image, mesh_size=8):
    """
    Perform deformable registration and output displacement field (numpy array with spatial information) and SimpleITK transform object.

    Args:
        fixed_image (sitk.Image): Fixed image (reference)
        moving_image (sitk.Image): Moving image (to be registered)
        mesh_size (int or list): B-spline mesh size (default 8, can be adjusted based on image size)

    Returns:
        dfield_transform (sitk.DisplacementFieldTransform): SimpleITK displacement transform
    """
    if isinstance(mesh_size, int):
        mesh_size = [mesh_size] * fixed_image.GetDimension()
    tx = sitk.BSplineTransformInitializer(fixed_image, mesh_size)

    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=100)
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=50,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10
    )
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetInitialTransform(tx, inPlace=False)
    # registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])       # 4 2 1
    # registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])   # 2 1 0
    # registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    outTx = registration_method.Execute(fixed_image, moving_image)

    displacement_field = sitk.TransformToDisplacementField(
        outTx,
        sitk.sitkVectorFloat64,
        fixed_image.GetSize(),
        fixed_image.GetOrigin(),
        fixed_image.GetSpacing(),
        fixed_image.GetDirection()
    )

    dfield_transform = sitk.DisplacementFieldTransform(displacement_field)
    return dfield_transform


def affine_to_homogeneous_matrix(transform: sitk.AffineTransform) -> np.ndarray:
    """
    Convert a SimpleITK AffineTransform to a 4x4 homogeneous transformation matrix.

    Args:
        transform (sitk.AffineTransform): The affine transform to convert.

    Returns:
        np.ndarray: 4x4 homogeneous matrix (dtype float64).
    """
    matrix3x3 = np.array(transform.GetMatrix()).reshape((3, 3))
    translation = np.array(transform.GetTranslation())
    hom_mat = np.eye(4)
    hom_mat[:3, :3] = matrix3x3
    hom_mat[:3, 3] = translation
    return hom_mat


def displacement_field_to_dict(transform: sitk.DisplacementFieldTransform) -> dict:
    """
    Convert a SimpleITK DisplacementFieldTransform to a dictionary.
    Args:
        transform (sitk.DisplacementFieldTransform): The displacement field transform to convert.
    Returns:
        dict: Dictionary containing displacement field information.
    """
    displacement_field = sitk.GetArrayFromImage(transform.GetDisplacementField())
    origin = transform.GetOrigin()
    spacing = transform.GetSpacing()
    direction = transform.GetDirection()
    size = transform.GetSize()
    return {
        "vector_grid": displacement_field,
        "origin": origin,
        "spacing": spacing,
        "direction": direction,
        "size": size,
    }
