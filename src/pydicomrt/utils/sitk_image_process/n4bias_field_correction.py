import SimpleITK as sitk

def n4bfc(input_image: sitk.Image, shrink_factor: int = 1) -> tuple[sitk.Image, sitk.Image]:
    """
    Apply N4 bias field (shading) correction to a CT image, with optional downsampling for speed.

    Args:
        input_image (sitk.Image): The input 3D CT image.
        shrink_factor (int): Downsampling factor. Default is 1 (no scaling).
            Use 2, 3, or 4 to significantly speed up the correction.

    Returns:
        tuple[sitk.Image, sitk.Image]: (corrected full-resolution image, 3D low-frequency bias map)
    """
    
    # 1. Ensure float precision (N4 requires Float32)
    input_image_float = sitk.Cast(input_image, sitk.sitkFloat32)
    
    # 2. Downsampling
    if shrink_factor > 1:
        shrink_filter = sitk.ShrinkImageFilter()
        # Shrink all three dimensions (X, Y, Z) of the 3D image
        shrink_filter.SetShrinkFactor(shrink_factor)
        work_image = shrink_filter.Execute(input_image_float)
    else:
        work_image = input_image_float
        
    # 3. Build foreground mask (on downsampled image for speed)
    # Otsu thresholding separates body from air (air background often interferes with correction)
    mask_image = sitk.OtsuThreshold(work_image, 0, 1, 200)
    
    # 4. Run N4 on the downsampled image
    corrector = sitk.N4BiasFieldCorrectionImageFilter()
    # We only need the internal B-spline model, not the downsampled corrected image
    corrector.Execute(work_image, mask_image)
    
    # 5. Reconstruct full-resolution correction map
    # Pass the original image as reference; SimpleITK interpolates the log bias field
    # to full size using the computed B-spline control points.
    log_bias_field = corrector.GetLogBiasFieldAsImage(input_image_float)
    
    # Convert log field to multiplicative map (Map = exp(log_map))
    bias_map = sitk.Exp(log_bias_field)
    
    # 6. Apply correction: divide original image by bias map (I_corrected = I_original / Map)
    corrected_image = sitk.Divide(input_image_float, bias_map)
    
    # Cast back to original pixel type (e.g. Int16) to save memory
    corrected_image = sitk.Cast(corrected_image, input_image.GetPixelID())
    
    return corrected_image, bias_map