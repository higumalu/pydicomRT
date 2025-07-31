import numpy as np
import glob
import os

from pathlib import Path
from pydicom import dcmread
from pydicom.dataset import Dataset


def sort_ds_list(ds_list: list[Dataset]) -> list[Dataset]:
    return sorted(ds_list, key=get_slice_position, reverse=False)

def get_slice_position(series_slice: Dataset):
    _, _, slice_direction = get_slice_directions(series_slice)
    return np.dot(slice_direction, series_slice.ImagePositionPatient)

def get_slice_directions(series_slice: Dataset):
    orientation = series_slice.ImageOrientationPatient
    row_direction = np.array(orientation[:3])
    column_direction = np.array(orientation[3:])
    slice_direction = np.cross(row_direction, column_direction)
    if not np.allclose(
        np.dot(row_direction, column_direction), 0.0, atol=1e-3
    ) or not np.allclose(np.linalg.norm(slice_direction), 1.0, atol=1e-3):
        raise Exception("Invalid Image Orientation (Patient) attribute")
    return row_direction, column_direction, slice_direction

def load_sorted_image_series(
    image_series_path: str | list[str],
    ) -> list[Dataset]:
    """
    Load a sorted list of DICOM image series from a directory or a list of paths.
    Args:
        image_series_path (str | list[str]): The path to the directory containing DICOM images or a list of paths to DICOM files.
    Returns:
        list[Dataset]: A sorted list of DICOM datasets.
    """
    dcm_path_list = []
    if isinstance(image_series_path, str) and Path(image_series_path).is_dir():
        # Check if directory exists
        if not Path(image_series_path).exists():
            raise FileNotFoundError(f"Directory {image_series_path} does not exist")
        dcm_path_list = glob.glob(os.path.join(image_series_path, "*.dcm"))
    # Check if image_series_path is a list of strings
    if isinstance(image_series_path, list) and all(isinstance(path, str) for path in image_series_path):
        # Validate each path in the list exists
        for path in image_series_path:
            if not Path(path).exists():
                raise FileNotFoundError(f"File {path} does not exist")
        dcm_path_list = image_series_path

    if len(dcm_path_list) == 0:
        raise FileNotFoundError(f"No DICOM files found in {image_series_path}")

    image_ds_list = []
    try:
        for dcm_path in dcm_path_list:
            image_ds_list.append(dcmread(dcm_path))
    except Exception as e:
        raise ValueError(f"Error reading DICOM file {dcm_path}: {e}")

    image_ds_list = sort_ds_list(image_ds_list)
    return image_ds_list


if __name__ == "__main__":
    orientation = [-1, 0, 0, 0, -1, 0]
    row_direction = np.array(orientation[:3])
    column_direction = np.array(orientation[3:])
    slice_direction = np.cross(row_direction, column_direction)
    print(slice_direction)
    image_position = [2, 2, 2]
    slice_idx = np.dot(slice_direction, image_position)
    print(slice_idx)
    print(sorted([-1, -2, -3, -4], reverse=False))
