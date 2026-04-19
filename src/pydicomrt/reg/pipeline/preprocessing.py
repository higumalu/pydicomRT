"""
Registration preprocessing module.

Provides image preprocessing utilities such as window clipping and normalization.
"""

from typing import Optional, Dict, List, Tuple, Union
import numpy as np
import SimpleITK as sitk

from pydicomrt.utils.sitk_transform import resample_to_reference_image


def window_clip(
    image: sitk.Image,
    clip_range: Union[List[float], Tuple[float, float]]
) -> sitk.Image:
    """
    Apply window clipping (clamping) to the image.

    Parameters
    ----------
    image : sitk.Image
        Input image.
    clip_range : Union[List[float], Tuple[float, float]]
        Clipping range [lower, upper].

    Returns
    -------
    sitk.Image
        Clipped image.
    """
    if not isinstance(clip_range, (list, tuple)) or len(clip_range) != 2:
        raise ValueError("clip_range must be a list or tuple of two elements [lower, upper]")
    
    lower = min(clip_range)
    upper = max(clip_range)
    return sitk.Clamp(image, lowerBound=lower, upperBound=upper)


def preprocess_image(
    image: sitk.Image,
    preprocess_config: Optional[Dict] = None
) -> sitk.Image:
    """
    Preprocess the image by applying multiple preprocessing steps according to the config.

    Parameters
    ----------
    image : sitk.Image
        Input image.
    preprocess_config : Optional[Dict], default = None
        Preprocessing config dict. Supported options:
        - 'window_clip': Union[List[float], Tuple[float, float]] - Window clip range [lower, upper].

    Returns
    -------
    sitk.Image
        Preprocessed image.

    Examples
    --------
    >>> import SimpleITK as sitk
    >>> from pydicomrt.reg.pipeline.preprocessing import preprocess_image
    >>>
    >>> # Window clipping
    >>> processed = preprocess_image(
    ...     image=ct_image,
    ...     preprocess_config={'window_clip': [-10, 500]}
    ... )
    """
    if preprocess_config is None:
        return image

    processed_image = image

    # Window clipping
    if 'window_clip' in preprocess_config:
        clip_range = preprocess_config['window_clip']
        processed_image = window_clip(processed_image, clip_range)

    return processed_image


def get_image_physical_extent(image: sitk.Image) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the physical extent (bounding box) of the image in physical space.

    Parameters
    ----------
    image : sitk.Image
        Input image.

    Returns
    -------
    min_corner : np.ndarray
        Minimum corner coordinates (3D) in physical space.
    max_corner : np.ndarray
        Maximum corner coordinates (3D) in physical space.
    """
    # Get image attributes
    origin = np.array(image.GetOrigin())
    spacing = np.array(image.GetSpacing())
    size = np.array(image.GetSize())
    direction = np.array(image.GetDirection()).reshape(image.GetDimension(), image.GetDimension())

    # Compute 8 corner positions in image coordinates
    # Corner indices: (0,0,0), (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1), (1,1,0), (1,1,1)
    corners_image_coords = np.array([
        [0, 0, 0],
        [0, 0, size[2] - 1],
        [0, size[1] - 1, 0],
        [0, size[1] - 1, size[2] - 1],
        [size[0] - 1, 0, 0],
        [size[0] - 1, 0, size[2] - 1],
        [size[0] - 1, size[1] - 1, 0],
        [size[0] - 1, size[1] - 1, size[2] - 1],
    ])

    # Transform corners from image coordinates to physical coordinates
    # physical_coord = origin + direction @ (spacing * image_coord)
    corners_physical_coords = []
    for corner in corners_image_coords:
        physical_coord = origin + direction @ (spacing * corner)
        corners_physical_coords.append(physical_coord)

    corners_physical_coords = np.array(corners_physical_coords)

    # Find min and max coordinates
    min_corner = np.min(corners_physical_coords, axis=0)
    max_corner = np.max(corners_physical_coords, axis=0)

    return min_corner, max_corner


def get_image_center(image: sitk.Image) -> np.ndarray:
    """
    Compute the center point of the image in physical space.

    Parameters
    ----------
    image : sitk.Image
        Input image.

    Returns
    -------
    np.ndarray
        Center point coordinates (3D) in physical space.
    """
    min_corner, max_corner = get_image_physical_extent(image)
    center = (min_corner + max_corner) / 2.0
    return center


def get_images_distance(image1: sitk.Image, image2: sitk.Image) -> float:
    """
    Compute the Euclidean distance between the centers of two images.

    Parameters
    ----------
    image1 : sitk.Image
        First image.
    image2 : sitk.Image
        Second image.

    Returns
    -------
    float
        Euclidean distance (mm) between the two image centers.
    """
    center1 = get_image_center(image1)
    center2 = get_image_center(image2)
    distance = np.linalg.norm(center1 - center2)
    return distance


def create_reference_image_from_extent(
    min_corner: np.ndarray,
    max_corner: np.ndarray,
    spacing: Tuple[float, float, float],
    direction: Tuple[float, ...],
    pixel_id: int = sitk.sitkFloat32
) -> sitk.Image:
    """
    Create a reference image from the given physical extent.

    Parameters
    ----------
    min_corner : np.ndarray
        Minimum corner coordinates (3D) in physical space.
    max_corner : np.ndarray
        Maximum corner coordinates (3D) in physical space.
    spacing : Tuple[float, float, float]
        Pixel spacing.
    direction : Tuple[float, ...]
        Direction matrix (9 elements).
    pixel_id : int, default = sitk.sitkFloat32
        Pixel type.

    Returns
    -------
    sitk.Image
        Created reference image.
    """
    # Compute physical size
    size_physical = max_corner - min_corner
    
    # Compute pixel size
    spacing_array = np.array(spacing)
    size_pixels = np.ceil(size_physical / spacing_array).astype(int)
    size_pixels = np.maximum(size_pixels, [1, 1, 1])

    # Create image
    ref_image = sitk.Image(size_pixels.tolist(), pixel_id)
    ref_image.SetOrigin(min_corner)
    ref_image.SetSpacing(spacing)
    ref_image.SetDirection(direction)
    
    return ref_image


def get_initial_rigid_transform(
    fixed_image: sitk.Image,
    moving_image: sitk.Image
) -> sitk.Transform:
    """
    Get initial rigid transform using SimpleITK's CenteredTransformInitializer.

    Used to roughly align two images before computing their intersection.

    Parameters
    ----------
    fixed_image : sitk.Image
        Reference image.
    moving_image : sitk.Image
        Image to be registered.

    Returns
    -------
    sitk.Transform
        Initial rigid transform.
    """
    # Get initial transform via CenteredTransformInitializer
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image,
        moving_image,
        sitk.VersorRigid3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )

    # Set center to (0, 0, 0) to avoid unnecessary offset
    initial_transform.SetCenter((0.0, 0.0, 0.0))

    return initial_transform


def get_intersection_extent(
    image1: sitk.Image,
    image2: sitk.Image
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], bool]:
    """
    Compute the intersection of two images' physical extents.

    Parameters
    ----------
    image1 : sitk.Image
        First image.
    image2 : sitk.Image
        Second image.

    Returns
    -------
    min_corner : Optional[np.ndarray]
        Minimum corner of the intersection in physical space, or None if no intersection.
    max_corner : Optional[np.ndarray]
        Maximum corner of the intersection in physical space, or None if no intersection.
    has_intersection : bool
        Whether an intersection exists.
    """
    # Compute physical extents of both images
    min1, max1 = get_image_physical_extent(image1)
    min2, max2 = get_image_physical_extent(image2)

    # Compute intersection
    intersection_min = np.maximum(min1, min2)
    intersection_max = np.minimum(max1, max2)

    # Check if intersection exists (min < max in each dimension)
    has_intersection = np.all(intersection_min < intersection_max)

    if has_intersection:
        return intersection_min, intersection_max, True
    else:
        return None, None, False


def crop_image_to_extent(
    image: sitk.Image,
    min_corner: np.ndarray,
    max_corner: np.ndarray,
    default_value: Optional[float] = None
) -> sitk.Image:
    """
    Crop the image to the given physical extent.

    Parameters
    ----------
    image : sitk.Image
        Input image.
    min_corner : np.ndarray
        Minimum corner coordinates (3D) in physical space.
    max_corner : np.ndarray
        Maximum corner coordinates (3D) in physical space.
    default_value : Optional[float], default = None
        Default pixel value. If None, auto-detected (CT images use -1000).

    Returns
    -------
    sitk.Image
        Cropped image.
    """
    # Auto-detect default value
    if default_value is None:
        default_value = 0.0
        try:
            arr_min = float(np.min(sitk.GetArrayViewFromImage(image)))
            if arr_min <= -1000.0:
                default_value = -1000.0
        except Exception:
            pass

    # Compute size of crop extent
    size_physical = max_corner - min_corner

    # Get image attributes
    spacing = np.array(image.GetSpacing())
    direction = np.array(image.GetDirection()).reshape(image.GetDimension(), image.GetDimension())

    # Compute cropped image size in pixels (physical size -> pixel size)
    size_pixels = np.ceil(size_physical / spacing).astype(int)
    size_pixels = np.maximum(size_pixels, [1, 1, 1])

    # Create new image with crop extent spatial info (origin = min_corner)
    cropped_image = sitk.Image(size_pixels.tolist(), image.GetPixelID())
    cropped_image.SetOrigin(min_corner)
    cropped_image.SetSpacing(spacing)
    cropped_image.SetDirection(image.GetDirection())

    # Resample original image onto crop extent
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(cropped_image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(default_value)
    resampler.SetTransform(sitk.Transform())

    result = resampler.Execute(image)

    return result


def align_image_extents(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    default_value: Optional[float] = None,
    use_initial_rigid: bool = True
) -> Tuple[sitk.Image, sitk.Image, Optional[sitk.Transform]]:
    """
    Align two images' physical extents by cropping them to their intersection.

    If use_initial_rigid=True, an initial rigid transform is applied first to roughly align
    the images, then intersection is computed and cropping is applied for a more accurate extent.

    Parameters
    ----------
    fixed_image : sitk.Image
        Reference image.
    moving_image : sitk.Image
        Image to be registered.
    default_value : Optional[float], default = None
        Default pixel value. If None, auto-detected (CT images use -1000).
    use_initial_rigid : bool, default = True
        Whether to use initial rigid transform to align images before computing intersection.

    Returns
    -------
    fixed_aligned : sitk.Image
        Aligned fixed_image (cropped to intersection).
    moving_aligned : sitk.Image
        Aligned moving_image (cropped to intersection).
    initial_transform : Optional[sitk.Transform]
        Initial rigid transform if used, else None.
    """
    initial_transform = None
    current_moving_image = moving_image

    # If allowed, apply initial rigid transform before computing intersection
    if use_initial_rigid:
        initial_transform = get_initial_rigid_transform(fixed_image, current_moving_image)

        # Apply transform to moving_image (resample moving_image into fixed_image space)
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(fixed_image)
        resampler.SetInterpolator(sitk.sitkLinear)

        # Auto-detect default value
        if default_value is None:
            default_val = 0.0
            try:
                arr_min = float(np.min(sitk.GetArrayViewFromImage(current_moving_image)))
                if arr_min <= -1000.0:
                    default_val = -1000.0
            except Exception:
                pass
        else:
            default_val = default_value
        
        resampler.SetDefaultPixelValue(default_val)
        resampler.SetTransform(initial_transform)

        current_moving_image = resampler.Execute(current_moving_image)

    # Compute intersection (based on aligned images)
    min_corner, max_corner, has_intersection = get_intersection_extent(
        fixed_image,
        current_moving_image
    )

    if not has_intersection:
        raise ValueError(
            "The two images have no intersection in physical space. "
            "No intersection found even with initial rigid transform. "
            "Check that the images' spatial extents overlap."
        )

    # Crop both images to intersection; use fixed_image spacing for consistency
    fixed_spacing = fixed_image.GetSpacing()
    fixed_aligned = crop_image_to_extent(fixed_image, min_corner, max_corner, default_value)
    moving_aligned = crop_image_to_extent(current_moving_image, min_corner, max_corner, default_value)

    # Ensure moving_aligned has same spacing and spatial attributes as fixed_aligned
    if (moving_aligned.GetSpacing() != fixed_spacing or 
        moving_aligned.GetSize() != fixed_aligned.GetSize() or
        moving_aligned.GetOrigin() != fixed_aligned.GetOrigin()):
        moving_aligned = resample_to_reference_image(fixed_aligned, moving_aligned)

    return fixed_aligned, moving_aligned, initial_transform
