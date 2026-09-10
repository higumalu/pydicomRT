"""
Helpers shared by the registration algorithms.

These existed as byte-identical copies in two or three modules each, which is how the
pyramid helper came to be correct in ``demons`` and wrong in ``bspline``: the latter
computed its Gaussian kernel width as ``8 * variance * spacing`` rather than
``8 * sigma / spacing``, which is not even dimensionally right. The version kept here is
the corrected one.
"""

import logging

import numpy as np
import SimpleITK as sitk

logger = logging.getLogger(__name__)

__all__ = [
    "smooth_and_resample",
    "apply_transform",
    "registration_command_iteration",
    "stage_iteration",
    "deformable_registration_command_iteration",
]


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
            sigmas = list(smoothing_sigma)
        else:
            smoothing_variance = (smoothing_sigma ** 2,) * image.GetDimension()
            sigmas = [smoothing_sigma] * image.GetDimension()

        # Kernel width = ~8*sigma/spacing (in voxels, rounded)
        # Convert sigma from mm to voxels: sigma_vox = sigma_mm / spacing
        maximum_kernel_width = int(
            max([8.0 * sigma / spacing for sigma, spacing in zip(sigmas, image.GetSpacing())])
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

    # ---- Step 4.5. Validate new_size to avoid division by zero ----
    # Ensure no dimension becomes less than 2 (which would cause division by zero in spacing calculation)
    for i, size_n_i in enumerate(new_size):
        if size_n_i < 1:
            raise ValueError(
                f"Computed new size for dimension {i} is {size_n_i}, which is less than 1. "
                f"This may occur if shrink_factor is too large or isotropic_voxel_size_mm is too small."
            )
        if size_n_i == 1:
            raise ValueError(
                f"Computed new size for dimension {i} is 1, which would cause division by zero "
                f"in spacing calculation. Please use a smaller shrink_factor or larger "
                f"isotropic_voxel_size_mm to avoid this issue."
            )

    # ---- Step 5. Compute new spacing from new size ----
    # Keep same physical extent → spacing = (extent / (new_size-1))
    # Note: new_size is guaranteed to be >= 2 at this point, so no division by zero
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
    logger.info("%3d = %10.5f", method.GetOptimizerIteration(), method.GetMetricValue())


def stage_iteration(method):
    """
    Utility function to print information during stage change in registration
    """
    logger.info("Number of parameters = %s", method.GetInitialTransform().GetNumberOfParameters())


def deformable_registration_command_iteration(method):
    """
    Utility function to print information during demons registration
    """
    logger.info("%3d = %10.5f", method.GetElapsedIterations(), method.GetMetric())
