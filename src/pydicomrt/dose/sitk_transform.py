
import numpy as np
import SimpleITK as sitk

def get_dose_array_spacing(dose_ds):
    x_spacing, y_spacing = dose_ds.PixelSpacing
    frame_offset_vector_list = dose_ds.GridFrameOffsetVector
    frame_distance_list = [np.linalg.norm(frame_offset_vector_list[i + 1] - frame_offset_vector_list[i])
                           for i in range(len(frame_offset_vector_list) - 1)]
    z_spacing = np.mean(frame_distance_list)
    std_distance = np.std(frame_distance_list)
    if std_distance > 0.001:
        print("Warning: dose array spacing is not uniform")
    return [x_spacing, y_spacing, z_spacing]

def get_dose_direction(dose_ds):
    ori = dose_ds.ImageOrientationPatient
    row = ori[:3]
    col = ori[3:]
    cross = np.cross(row, col)
    dose_direction = np.asarray([row, col, cross]).T
    return dose_direction

def get_dose_origin(dose_ds):
    return dose_ds.ImagePositionPatient

def get_dose_array(dose_ds):
    return dose_ds.pixel_array * dose_ds.DoseGridScaling

def get_dose_sitk_image(dose_ds):
    dose_array = get_dose_array(dose_ds)
    spacing = get_dose_array_spacing(dose_ds)
    origin = get_dose_origin(dose_ds)
    direction = get_dose_direction(dose_ds)
    sitk_image = sitk.GetImageFromArray(dose_array)
    sitk_image.SetSpacing(spacing)
    sitk_image.SetOrigin(origin)
    sitk_image.SetDirection(direction.flatten())
    return sitk_image
