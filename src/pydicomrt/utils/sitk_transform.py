"""
A SimpleITK image Builder.
Date: 2025-07-18
Author: higumalu
"""

import numpy as np
import SimpleITK as sitk

from pydicomrt.utils.image_series_loader import sort_ds_list, get_slice_position, load_sorted_image_series, get_slice_directions

class SimpleITKImageBuilder:
    """
    Build a SimpleITK image
    """
    def __init__(self):
        self.volume = None
        self.origin = None
        self.spacing = None
        self.direction = None
        self.image = None

    def set_volume(self, volume: np.ndarray):
        assert len(volume.shape) == 3, "Volume must be a 3D array"
        self.volume = volume

    def set_origin(self, origin: np.ndarray):
        assert origin.shape == (3,), "Origin must be a 3-element vector"
        self.origin = origin

    def set_spacing(self, spacing: np.ndarray):
        assert spacing.shape == (3,), "Spacing must be a 3-element vector"
        self.spacing = spacing

    def set_direction(self, direction: np.ndarray):
        assert direction.shape == (3, 3), "Direction must be a 3x3 matrix"
        self.direction = direction

    def from_ds_list(self, ds_list):
        self.volume, self.origin, self.spacing, self.direction = parse_ds_list(ds_list)
        return self.build()

    def from_ref_sitk_image(self, volume: np.ndarray, ref_image: sitk.Image):
        self.volume = volume
        self.origin = np.array(ref_image.GetOrigin())
        self.spacing = np.array(ref_image.GetSpacing())
        self.direction = np.array(ref_image.GetDirection()).reshape(3, 3)
        return self.build()

    def from_dcms_dir(self, dcms_dir: str):
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(dcms_dir)
        reader.SetFileNames(dicom_names)
        image = reader.Execute()
        self.volume = sitk.GetArrayFromImage(image)
        self.origin = np.array(image.GetOrigin())
        self.spacing = np.array(image.GetSpacing())
        self.direction = np.array(image.GetDirection()).reshape(3, 3)
        return self.build()

    def build(self):
        assert self.volume is not None, "Volume must be set"
        assert self.origin is not None, "Origin must be set"
        assert self.spacing is not None, "Spacing must be set"
        assert self.direction is not None, "Direction must be set"
        self.image = sitk.GetImageFromArray(self.volume)
        self.image.SetOrigin(self.origin)
        self.image.SetSpacing(self.spacing)
        self.image.SetDirection(self.direction.ravel())
        return self.image

def get_spacing_between_slices(series_data):
    if len(series_data) > 1:
        first = get_slice_position(series_data[0])
        last = get_slice_position(series_data[-1])
        return (last - first) / (len(series_data) - 1)
    # Return nonzero value for one slice just to make the transformation matrix invertible
    return 1.0

def parse_ds_list(ds_list):
    """
    Parse a list of DICOM datasets
    Limitations:
        Only works for axial slices.
        Only works for linear transforms.
        Only works for DICOM datasets with ImagePositionPatient and PixelSpacing attributes.
        Only works for CT/MR DICOM datasets.
    Args:
        ds_list (list): List of DICOM datasets.
    Returns:
        tuple: Tuple of image volume, origin, spacing, and direction.
    """
    ds_list = sort_ds_list(ds_list)
    image_volume = np.stack([ds.pixel_array for ds in ds_list], axis=0)
    slope = getattr(ds_list[0], "RescaleSlope", 1)
    intercept = getattr(ds_list[0], "RescaleIntercept", 0)
    image_volume = image_volume * slope + intercept
    slice_spacing = get_spacing_between_slices(ds_list)
    origin = np.array([float(ds_list[0].ImagePositionPatient[0]), float(ds_list[0].ImagePositionPatient[1]), float(ds_list[0].ImagePositionPatient[2])])
    spacing = np.array([float(ds_list[0].PixelSpacing[0]), float(ds_list[0].PixelSpacing[1]), slice_spacing])
    direction = np.array(get_slice_directions(ds_list[0])).reshape(3, 3)
    return image_volume, origin, spacing, direction

def ds_list_to_sitk_image(ds_list):
    """
    Convert a list of DICOM datasets to a SimpleITK image.
    Args:
        ds_list (list): List of DICOM datasets.
    Returns:
        sitk.Image: SimpleITK image.
    Limitations:
        - Only works for axial slices.
    """
    ds_list = sort_ds_list(ds_list)
    image_volume = np.stack([ds.pixel_array for ds in ds_list], axis=0)
    slope = getattr(ds_list[0], "RescaleSlope", 1)
    intercept = getattr(ds_list[0], "RescaleIntercept", 0)
    image_volume = image_volume * slope + intercept
    slice_spacing = get_spacing_between_slices(ds_list)
    spacing = [float(ds_list[0].PixelSpacing[0]), float(ds_list[0].PixelSpacing[1]), slice_spacing]
    origin = [float(ds_list[0].ImagePositionPatient[0]), float(ds_list[0].ImagePositionPatient[1]), float(ds_list[0].ImagePositionPatient[2])]
    direction = np.array(get_slice_directions(ds_list[0])).reshape(3, 3).ravel()
    sitk_image = sitk.GetImageFromArray(image_volume)
    sitk_image.SetOrigin(origin)
    sitk_image.SetSpacing(spacing)
    sitk_image.SetDirection(direction)
    return sitk_image

def resample_to_reference_image(reference_image, source_image, min_val=-1000.0):
    """
    Resample the input image to the reference image space.
    Args:
        input_image (sitk.Image): Input image.
        target_image (sitk.Image): Target image.
    Returns:
        sitk.Image: Input image with synchronized spacing, origin, and direction.
    """
    resample_image = sitk.Resample(
        source_image,
        reference_image,
        sitk.Transform(),
        sitk.sitkLinear,
        min_val,
        reference_image.GetPixelID())
    return resample_image


if __name__ == "__main__":
    ct_path = "example/data/Mirror/CT"
    ct_ds_list = load_sorted_image_series(ct_path)

    ct_builder = SimpleITKImageBuilder()
    ct_builder.from_dcms_dir(ct_path)
    ct_image = ct_builder.build()
    print(ct_image)

    ct_builder.from_ds_list(ct_ds_list)
    ct_image = ct_builder.build()
    print(ct_image)

    ct_volume = np.random.rand(10, 10, 10)
    ct_builder.from_ref_sitk_image(ct_volume, ct_image)
    ct_image = ct_builder.build()
    print(ct_image)

    # ct_ds_list = load_sorted_image_series(ct_path)
    # image = ds_list_to_sitk_image(ct_ds_list)
    # print(image)
