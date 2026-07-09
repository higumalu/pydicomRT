import SimpleITK as sitk

def median_denoise(input_image: sitk.Image, radius: int = 1) -> sitk.Image:
    """
    Denoise using a median filter.
    
    Args:
        input_image (sitk.Image): Input image.
        radius (int): Filter radius.
            radius = 1 -> 3×3×3 neighborhood (median of 27 voxels)
            radius = 2 -> 5×5×5 neighborhood
            
    Returns:
        sitk.Image: Denoised image.
    """
    # Median filtering often works on integer types (e.g. Int16); float avoids some overflow edge cases
    original_type = input_image.GetPixelID()
    img_float = sitk.Cast(input_image, sitk.sitkFloat32)

    # Build median filter
    median_filter = sitk.MedianImageFilter()
    median_filter.SetRadius(radius)
    
    # Execute
    output = median_filter.Execute(img_float)
    
    # Cast back to original type
    return sitk.Cast(output, original_type)
