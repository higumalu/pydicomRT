"""
Synthetic CT and MR series for the examples.

The examples are runnable without any patient data, so they can be executed as-is to check
that a change to the library did not break the documented workflow. Swap
``make_demo_series()`` for ``load_sorted_image_series("/path/to/dicom")`` to use real data.
"""

import numpy as np
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, MRImageStorage, generate_uid

AXIAL_ORIENTATION = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]


# Phantom inserts: (dx, dy, dz) from the centre in mm, radius in mm, HU.
INSERTS = [(-20, -8, -6, 14, 300), (18, 9, 8, 11, -700), (2, -14, 10, 9, 900)]

# CT number -> MR signal for the demo MR: not monotonic, the way CT -> T1 is not. Air stays
# dark, the -700 (lung-like) insert turns bright, the 900 (bone-like) one nearly black.
MR_SIGNAL = ([-1000, -700, 0, 300, 900], [0, 900, 500, 650, 40])


def _phantom_hu(shape, pixel_spacing, slice_spacing, bulge_mm=0.0):
    """
    The phantom in HU on a (slice, row, column) grid, in mm from the first voxel.

    ``bulge_mm`` deforms it: the anatomy at each point is taken from a point displaced along
    x by a smooth bump -- ``bulge_mm`` at the centre, fading over ~25 mm -- so a deformable
    registration has something to find.
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
    if bulge_mm:
        x_mm = x_mm + bulge_displacement(x_mm - centre[0], y_mm - centre[1], z_mm - centre[2], bulge_mm)

    volume = np.full(z_mm.shape, -1000.0)
    # Deliberately elliptical, wider than it is tall. A circular body is rotationally
    # symmetric enough that a 90-degree rotation becomes a competing optimum and the
    # registration examples land on it instead of the intended alignment.
    body = ((x_mm - centre[0]) / (centre[0] * 0.85)) ** 2 + \
           ((y_mm - centre[1]) / (centre[1] * 0.55)) ** 2 < 1
    volume[body] = 0.0
    for dx, dy, dz, radius, value in INSERTS:
        distance = np.sqrt(
            (x_mm - centre[0] - dx) ** 2 + (y_mm - centre[1] - dy) ** 2 + (z_mm - centre[2] - dz) ** 2
        )
        volume[distance < radius] = value
    return volume


def bulge_displacement(dx, dy, dz, bulge_mm, width_mm=25.0):
    """x-displacement of the demo deformation at an offset (dx, dy, dz) from the centre."""
    return bulge_mm * np.exp(-(dx ** 2 + dy ** 2 + dz ** 2) / (2.0 * width_mm ** 2))


def _series(stored, origin, pixel_spacing, slice_spacing, shift_mm, modality, sop_class, rescale):
    slices, rows, cols = stored.shape
    study_uid, series_uid, frame_uid = generate_uid(), generate_uid(), generate_uid()
    ds_list = []
    for index in range(slices):
        file_meta = FileMetaDataset()
        file_meta.MediaStorageSOPClassUID = sop_class
        file_meta.MediaStorageSOPInstanceUID = generate_uid()
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        ds = FileDataset("slice.dcm", {}, file_meta=file_meta, preamble=b"\0" * 128)
        ds.SOPClassUID = sop_class
        ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
        ds.StudyInstanceUID, ds.SeriesInstanceUID = study_uid, series_uid
        ds.FrameOfReferenceUID = frame_uid
        ds.Modality, ds.PatientID, ds.PatientName = modality, "DEMO-001", "Demo^Phantom"
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
        if rescale:
            ds.RescaleSlope, ds.RescaleIntercept = 1, -1024
        ds.PixelData = stored[index].tobytes()

        ds_list.append(ds)
    return ds_list


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
    volume = _phantom_hu(shape, pixel_spacing, slice_spacing)
    if noise_sigma:
        volume = volume + np.random.default_rng(0).normal(0.0, noise_sigma, volume.shape)
    stored = np.clip(volume + 1024, 0, 4095).astype(np.uint16)
    return _series(stored, origin, pixel_spacing, slice_spacing, shift_mm, "CT", CTImageStorage, True)


def make_demo_mr_series(
    shape=(24, 96, 96),
    origin=(-60.0, -60.0, -30.0),
    pixel_spacing=(1.25, 1.25),
    slice_spacing=2.5,
    shift_mm=(0.0, 0.0, 0.0),
    bulge_mm=0.0,
    noise_sigma=20.0,
):
    """
    The same phantom as an MR series: a different, non-monotonic contrast (``MR_SIGNAL``),
    optionally deformed by ``bulge_mm`` (see ``_phantom_hu``), and shifted by ``shift_mm``.
    MR has no physical unit, so there is no rescale -- stored values are the signal.
    """
    hu = _phantom_hu(shape, pixel_spacing, slice_spacing, bulge_mm)
    signal = np.interp(hu, *MR_SIGNAL)
    if noise_sigma:
        signal = signal + np.random.default_rng(1).normal(0.0, noise_sigma, signal.shape)
    stored = np.clip(signal, 0, 4095).astype(np.uint16)
    return _series(stored, origin, pixel_spacing, slice_spacing, shift_mm, "MR", MRImageStorage, False)


def true_mr_point(ct_point, shape, origin, pixel_spacing, slice_spacing, shift_mm, bulge_mm):
    """
    Where a CT point's anatomy sits in the demo MR -- the answer a CT -> MR registration should
    give. Inverts the MR's shift and bulge by fixed-point iteration (the bulge is small and
    smooth, so a handful of steps converge to well under a micrometre).
    """
    slices, rows, cols = shape
    centre = np.array([(cols - 1) * pixel_spacing[1], (rows - 1) * pixel_spacing[0],
                       (slices - 1) * slice_spacing]) / 2.0
    target = np.asarray(ct_point, dtype=float) - np.asarray(origin)  # mm from the first voxel
    local = target.copy()
    for _ in range(50):
        offset = local - centre
        local = target - np.array([bulge_displacement(*offset, bulge_mm), 0.0, 0.0])
    return local + np.asarray(origin) + np.asarray(shift_mm)
