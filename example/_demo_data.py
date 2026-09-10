"""
Synthetic CT series for the examples.

The examples are runnable without any patient data, so they can be executed as-is to check
that a change to the library did not break the documented workflow. Swap
``make_demo_series()`` for ``load_sorted_image_series("/path/to/dicom")`` to use real data.
"""

import numpy as np
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

AXIAL_ORIENTATION = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]


def make_demo_series(
    shape=(24, 96, 96),
    origin=(-60.0, -60.0, -30.0),
    pixel_spacing=(1.25, 1.25),
    slice_spacing=2.5,
    shift_mm=(0.0, 0.0, 0.0),
    noise_sigma=15.0,
):
    """
    Build a CT-like axial series containing a phantom with a few distinct structures.

    ``shift_mm`` moves the series origin, which is how the examples produce a "moving"
    series that is misaligned from the "fixed" one by a known amount.

    ``noise_sigma`` is not decoration. A noiseless phantom is piecewise constant, so its
    joint histogram collapses into a handful of bins and mutual information becomes a poor,
    spiky objective -- intensity-based registration on it behaves far worse than on real CT.
    """
    slices, rows, cols = shape
    z_mm, y_mm, x_mm = np.meshgrid(
        np.arange(slices) * slice_spacing,
        np.arange(rows) * pixel_spacing[0],
        np.arange(cols) * pixel_spacing[1],
        indexing="ij",
    )
    extent = [(cols - 1) * pixel_spacing[1], (rows - 1) * pixel_spacing[0], (slices - 1) * slice_spacing]
    centre = [e / 2.0 for e in extent]

    volume = np.full(z_mm.shape, -1000.0)
    # Deliberately elliptical, wider than it is tall. A circular body is rotationally
    # symmetric enough that a 90-degree rotation becomes a competing optimum and the
    # registration examples land on it instead of the intended alignment.
    body = ((x_mm - centre[0]) / (centre[0] * 0.85)) ** 2 + \
           ((y_mm - centre[1]) / (centre[1] * 0.55)) ** 2 < 1
    volume[body] = 0.0
    for dx, dy, dz, radius, value in [(-20, -8, -6, 14, 300), (18, 9, 8, 11, -700), (2, -14, 10, 9, 900)]:
        distance = np.sqrt(
            (x_mm - centre[0] - dx) ** 2 + (y_mm - centre[1] - dy) ** 2 + (z_mm - centre[2] - dz) ** 2
        )
        volume[distance < radius] = value

    if noise_sigma:
        volume = volume + np.random.default_rng(0).normal(0.0, noise_sigma, volume.shape)

    stored = np.clip(volume + 1024, 0, 4095).astype(np.uint16)

    study_uid, series_uid, frame_uid = generate_uid(), generate_uid(), generate_uid()
    ds_list = []
    for index in range(slices):
        file_meta = FileMetaDataset()
        file_meta.MediaStorageSOPClassUID = CTImageStorage
        file_meta.MediaStorageSOPInstanceUID = generate_uid()
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        ds = FileDataset("slice.dcm", {}, file_meta=file_meta, preamble=b"\0" * 128)
        ds.SOPClassUID = CTImageStorage
        ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
        ds.StudyInstanceUID, ds.SeriesInstanceUID = study_uid, series_uid
        ds.FrameOfReferenceUID = frame_uid
        ds.Modality, ds.PatientID, ds.PatientName = "CT", "DEMO-001", "Demo^Phantom"
        ds.PatientBirthDate, ds.PatientSex, ds.PatientPosition = "19700101", "O", "HFS"
        ds.StudyID, ds.StudyDate, ds.StudyTime = "1", "20260101", "120000"
        ds.SeriesNumber, ds.InstanceNumber = 1, index + 1
        ds.AccessionNumber = ds.ReferringPhysicianName = ds.PositionReferenceIndicator = ""

        ds.ImagePositionPatient = [
            origin[0] + shift_mm[0],
            origin[1] + shift_mm[1],
            origin[2] + shift_mm[2] + index * slice_spacing,
        ]
        ds.ImageOrientationPatient = list(AXIAL_ORIENTATION)
        ds.PixelSpacing = [float(pixel_spacing[0]), float(pixel_spacing[1])]
        ds.SliceThickness = float(slice_spacing)

        ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
        ds.Rows, ds.Columns = rows, cols
        ds.BitsAllocated = ds.BitsStored = 16
        ds.HighBit, ds.PixelRepresentation = 15, 0
        ds.RescaleSlope, ds.RescaleIntercept = 1, -1024
        ds.PixelData = stored[index].tobytes()

        ds_list.append(ds)
    return ds_list
