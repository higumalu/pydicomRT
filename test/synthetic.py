"""
Synthetic DICOM and image builders shared by the test suite.

Kept out of conftest.py so test modules can import these by name: there is a second
conftest.py under test/real_data/, and two files with that basename cannot both be imported
as the top-level ``conftest`` module.
"""

import numpy as np
import SimpleITK as sitk
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid


# Axial, head-first-supine: rows run +x, columns run +y, slices stack along +z.
AXIAL_ORIENTATION = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]


def make_ct_slice(
    pixel_array: np.ndarray,
    image_position: tuple,
    pixel_spacing: tuple,
    slice_thickness: float,
    study_uid: str,
    series_uid: str,
    frame_of_reference_uid: str,
    instance_number: int,
    orientation=AXIAL_ORIENTATION,
    rescale_slope: float = 1.0,
    rescale_intercept: float = -1024.0,
) -> FileDataset:
    """
    Build one CT slice as a decodable pydicom ``FileDataset``.

    ``pixel_array`` holds stored values (uint16); the HU an application sees is
    ``stored * rescale_slope + rescale_intercept``.
    """
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = CTImageStorage
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset("slice.dcm", {}, file_meta=file_meta, preamble=b"\0" * 128)

    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.FrameOfReferenceUID = frame_of_reference_uid

    ds.Modality = "CT"
    ds.PatientID = "TEST-0001"
    ds.PatientName = "Synthetic^Phantom"
    ds.PatientBirthDate = "19700101"
    ds.PatientSex = "O"
    ds.PatientPosition = "HFS"
    ds.StudyID = "1"
    ds.StudyDate = "20260101"
    ds.StudyTime = "120000"
    ds.SeriesNumber = 1
    ds.InstanceNumber = instance_number
    ds.AccessionNumber = ""
    ds.ReferringPhysicianName = ""
    ds.PositionReferenceIndicator = ""

    ds.ImagePositionPatient = [float(v) for v in image_position]
    ds.ImageOrientationPatient = [float(v) for v in orientation]
    # DICOM order is [row spacing (along the column direction), column spacing].
    ds.PixelSpacing = [float(pixel_spacing[0]), float(pixel_spacing[1])]
    ds.SliceThickness = float(slice_thickness)

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.Rows, ds.Columns = pixel_array.shape
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.RescaleSlope = rescale_slope
    ds.RescaleIntercept = rescale_intercept
    ds.PixelData = pixel_array.astype(np.uint16).tobytes()

    return ds


def make_ct_series(
    volume: np.ndarray = None,
    shape: tuple = (8, 16, 16),
    origin: tuple = (-50.0, -60.0, -20.0),
    pixel_spacing: tuple = (2.0, 2.0),
    slice_spacing: float = 3.0,
    orientation=AXIAL_ORIENTATION,
    rescale_intercept: float = -1024.0,
) -> list:
    """
    Build a synthetic axial CT series, ordered from the first slice upward.

    Parameters
    ----------
    volume : np.ndarray, optional
        Stored values shaped ``(z, rows, cols)``. A smooth gradient phantom is generated
        when omitted.
    shape : tuple
        ``(slices, rows, cols)`` used when ``volume`` is None.
    pixel_spacing : tuple
        DICOM ``PixelSpacing``: ``(row_spacing, column_spacing)``. Deliberately allowed to
        be anisotropic so tests can pin down which axis each value lands on.
    """
    if volume is None:
        slices, rows, cols = shape
        z_idx, y_idx, x_idx = np.meshgrid(
            np.arange(slices), np.arange(rows), np.arange(cols), indexing="ij"
        )
        volume = (1024 + 10 * x_idx + 20 * y_idx + 40 * z_idx).astype(np.uint16)

    study_uid = generate_uid()
    series_uid = generate_uid()
    frame_of_reference_uid = generate_uid()

    ds_list = []
    for index in range(volume.shape[0]):
        position = (origin[0], origin[1], origin[2] + index * slice_spacing)
        ds_list.append(
            make_ct_slice(
                pixel_array=volume[index],
                image_position=position,
                pixel_spacing=pixel_spacing,
                slice_thickness=slice_spacing,
                study_uid=study_uid,
                series_uid=series_uid,
                frame_of_reference_uid=frame_of_reference_uid,
                instance_number=index + 1,
                orientation=orientation,
                rescale_intercept=rescale_intercept,
            )
        )
    return ds_list


def make_blob_image(
    size=(40, 40, 24),
    spacing=(2.0, 2.0, 3.0),
    origin=(0.0, 0.0, 0.0),
    centre_offset=(0.0, 0.0, 0.0),
    radius_mm: float = 22.0,
    background: float = -1000.0,
    foreground: float = 300.0,
) -> sitk.Image:
    """
    A CT-like sphere on an air background, positioned in physical space.

    Registration needs actual gradients to lock onto, so the sphere edge is smoothed rather
    than left as a hard binary step.
    """
    grid = np.meshgrid(
        *[np.arange(n) * s for n, s in zip(size[::-1], spacing[::-1])], indexing="ij"
    )
    extent_mm = [(n - 1) * s for n, s in zip(size, spacing)]
    centre = [e / 2.0 + o for e, o in zip(extent_mm, centre_offset)]

    # grid is (z, y, x); centre/extent are (x, y, z)
    distance = np.sqrt(
        (grid[0] - centre[2]) ** 2 + (grid[1] - centre[1]) ** 2 + (grid[2] - centre[0]) ** 2
    )

    # Smooth falloff over one voxel-ish so the metric has a usable gradient
    edge = np.clip((radius_mm - distance) / max(spacing), 0.0, 1.0)
    array = background + (foreground - background) * edge

    image = sitk.GetImageFromArray(array.astype(np.float32))
    image.SetSpacing([float(s) for s in spacing])
    image.SetOrigin([float(o) for o in origin])
    return image


def make_textured_phantom(
    size=(48, 48, 32),
    spacing=(2.0, 2.0, 3.0),
    origin=(0.0, 0.0, 0.0),
    noise_sigma: float = 15.0,
) -> sitk.Image:
    """
    A CT-like phantom with enough structure to register against.

    A single uniform blob is a degenerate target for mutual information -- only its edge
    carries any signal. This adds an oval "body" plus off-centre inserts at distinct HU
    values and light noise, so the metric has a well-conditioned optimum in every axis.
    """
    rng = np.random.default_rng(0)

    axes = [np.arange(n) * s for n, s in zip(size[::-1], spacing[::-1])]
    z_mm, y_mm, x_mm = np.meshgrid(*axes, indexing="ij")

    extent = [(n - 1) * s for n, s in zip(size, spacing)]
    centre = [e / 2.0 for e in extent]

    array = np.full(z_mm.shape, -1000.0)
    body = ((x_mm - centre[0]) / (centre[0] * 0.8)) ** 2 + \
           ((y_mm - centre[1]) / (centre[1] * 0.8)) ** 2 < 1
    array[body] = 0.0

    inserts = [(-18, -10, -8, 12, 300), (20, 14, 10, 10, -700),
               (4, -22, 14, 9, 900), (-12, 18, -16, 8, 150)]
    for dx, dy, dz, radius, value in inserts:
        distance = np.sqrt(
            (x_mm - centre[0] - dx) ** 2
            + (y_mm - centre[1] - dy) ** 2
            + (z_mm - centre[2] - dz) ** 2
        )
        array[distance < radius] = value

    array += rng.normal(0.0, noise_sigma, array.shape)

    image = sitk.GetImageFromArray(array.astype(np.float32))
    image.SetSpacing([float(s) for s in spacing])
    image.SetOrigin([float(o) for o in origin])
    return image


def phantom_centre(size=(48, 48, 32), spacing=(2.0, 2.0, 3.0)):
    return [(n - 1) * s / 2.0 for n, s in zip(size, spacing)]
