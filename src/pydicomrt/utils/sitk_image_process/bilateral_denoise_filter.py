import SimpleITK as sitk

def bilateral_denoise(
    input_image: sitk.Image, 
    domain_sigma: float = 1.0, 
    range_sigma: float = 50.0,
    shrink_factor: int = 1
) -> sitk.Image:
    """
    Denoise the volume with a bilateral filter while preserving edge detail.

    Args:
        input_image (sitk.Image): Input 3D image.
        domain_sigma (float): Spatial standard deviation; larger values smooth over a wider neighborhood.
        range_sigma (float): Intensity standard deviation; larger values behave more like Gaussian blur;
            smaller values are more edge-preserving.
        shrink_factor (int): If greater than 1, process a shrunk image then resample back to full size
            (for testing or speed). Note: denoising is usually run at full resolution unless compute is very limited.

    Returns:
        sitk.Image: Denoised image with the original pixel type.
    """

    # 1. Use floating-point input (bilateral filter expects float)
    original_pixel_type = input_image.GetPixelID()
    input_image_float = sitk.Cast(input_image, sitk.sitkFloat32)

    # 2. Optional downsampling
    if shrink_factor > 1:
        work_image = sitk.Shrink(input_image_float, [shrink_factor] * input_image.GetDimension())
    else:
        work_image = input_image_float

    # 3. Configure and run bilateral filter
    # BilateralImageFilter is invoked via Execute in SimpleITK
    # Domain sigma is often set to ~0.5–2.0× the image spacing
    # Range sigma is often set to ~1–2× the noise standard deviation
    bf_filter = sitk.BilateralImageFilter()
    bf_filter.SetDomainSigma(domain_sigma)
    bf_filter.SetRangeSigma(range_sigma)
    
    denoised_image = bf_filter.Execute(work_image)

    # 4. Resample back to original size if shrink was used
    if shrink_factor > 1:
        # Linear interpolation to upsample to original grid
        denoised_image = sitk.Resample(
            denoised_image, 
            input_image_float, 
            sitk.Transform(), 
            sitk.sitkLinear, 
            0.0, 
            input_image_float.GetPixelID()
        )

    # 5. Cast back to original pixel type (e.g. Int16)
    final_image = sitk.Cast(denoised_image, original_pixel_type)

    return final_image
