"""
CT Image Builder
Date: 2026-02-14
Author: higumalu
"""

import copy
import datetime
import logging
from typing import List, Sequence

import numpy as np
import SimpleITK as sitk
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.uid import ImplicitVRLittleEndian, generate_uid

from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

logger = logging.getLogger(__name__)

#: Stored values are unsigned 16-bit with this intercept, so HU = stored + RESCALE_INTERCEPT.
RESCALE_INTERCEPT = -1024
_STORED_MAX = 2 ** 16 - 1

#: Which SimpleITK axis plays which DICOM role, per plane.
#:
#: DICOM counts a *row* along the column index and a *column* along the row index, and
#: SimpleITK reports (x, y, z) while numpy arrays come back (z, y, x). Writing the index
#: juggling out once, per plane, keeps the three views from drifting apart -- SAGITTAL and
#: CORONAL used to be swapped because each was hand-indexed separately.
#:
#: (column axis, row axis, slice-normal axis) as SimpleITK axis numbers.
PLANE_AXES = {
    "AXIAL": (0, 1, 2),      # normal along z
    "CORONAL": (0, 2, 1),    # normal along y
    "SAGITTAL": (1, 2, 0),   # normal along x
}


class CTBuilder:
    """
    Emit a CT Image series from a volume, borrowing patient and study context from a
    reference series.

    The output is DERIVED/SECONDARY: it is a rendering of a volume, not an acquisition.
    It gets its own SeriesInstanceUID, so it cannot be mistaken for -- or merged into --
    the acquired series it borrowed context from.

    Parameters
    ----------
    image_series : list of Dataset
        Reference slices supplying patient identity, study and frame of reference. Only
        the first is read for context unless :meth:`set_copy_all_attributes` is on.

    Methods
    -------
    set_volume(sitk_image)
        The volume to slice, in HU. Required (or ``set_volume_from_array``).
    set_volume_from_array(volume, origin, spacing, direction)
        The same, from numpy plus geometry.
    set_plane(plan_view)
        AXIAL, CORONAL or SAGITTAL. Default AXIAL.
    set_copy_all_attributes(copy_all_attributes)
        Inherit every attribute from the reference series.
    set_uid_prefix(uid_prefix)
        Set the UID root.
    build()
        Slice the volume into one ``Dataset`` per slice.

    Raises
    ------
    ValueError
        If ``image_series`` is empty.

    See Also
    --------
    check_ct_iod : Validate the emitted slices.
    image_series_to_sitk_image : The reverse direction.

    Notes
    -----
    Values are stored as unsigned 16-bit with a -1024 intercept, so the representable
    range is -1024 to 64511 HU. Values outside it are clipped, not wrapped.

    Examples
    --------
    >>> slices = (                                                  # doctest: +SKIP
    ...     CTBuilder(reference_series)
    ...     .set_volume(sitk_image)
    ...     .set_plane("AXIAL")
    ...     .build()
    ... )
    >>> for i, s in enumerate(slices):                              # doctest: +SKIP
    ...     s.save_as(f"ct_{i:04d}.dcm", enforce_file_format=True)
    """

    def __init__(self, image_series: List[Dataset]):
        if not image_series:
            raise ValueError("image_series must be non-empty")
        self.image_series = list(image_series)
        self.uid_prefix = "1.2.826.0.1.3680043.8.498."  # PYDICOM_ROOT_UID
        self.ref_ds = self.image_series[0]
        self._base_file_name = "ct_image"
        self.template_ds = copy.deepcopy(self.image_series[0])
        self.volume_image = None
        self.plan_view = "AXIAL"
        self.copy_all_attributes = False

    def set_uid_prefix(self, uid_prefix: str) -> "CTBuilder":
        """
        Set the root under which generated UIDs are minted.

        Parameters
        ----------
        uid_prefix : str
            Organisation UID root, dot-terminated. Unlike the RTSTRUCT and dose builders,
            this one does **not** read ``DICOM_UID_PREFIX``; it defaults to pydicom's
            root, which you do not own.

        Returns
        -------
        CTBuilder
            self, so calls chain.
        """
        self.uid_prefix = uid_prefix
        return self

    def set_volume(self, sitk_image: sitk.Image) -> "CTBuilder":
        """
        Set the volume to slice.

        Parameters
        ----------
        sitk_image : sitk.Image
            3D. Intensities are treated as Hounsfield units and written with a -1024
            intercept, so values outside -1024 to 64511 HU are clipped. Geometry --
            origin, spacing, direction -- is taken from the image and becomes the output
            series' geometry.

        Returns
        -------
        CTBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If the image is not 3D.
        """
        if sitk_image.GetDimension() != 3:
            raise ValueError(f"volume must be 3D, got {sitk_image.GetDimension()}D")
        self.volume_image = sitk_image
        return self

    def set_volume_from_array(
        self,
        volume: np.ndarray,
        origin: Sequence[float],
        spacing: Sequence[float],
        direction: Sequence[float],
        ) -> "CTBuilder":
        """
        Set the volume from a numpy array plus its geometry.

        Parameters
        ----------
        volume : np.ndarray
            3D, ``(slice, row, column)``. Treated as HU.
        origin : sequence of float
            Patient coordinates of voxel ``[0, 0, 0]``, in mm.
        spacing : sequence of float
            ``(x, y, z)`` in mm -- SimpleITK's order, not DICOM's ``[row, column]``.
        direction : sequence of float
            3x3 or flat 9-sequence, axis vectors in **columns**.

        Returns
        -------
        CTBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If the volume is not 3D or any geometry argument is the wrong shape.

        See Also
        --------
        set_volume : When you already have a ``sitk.Image``.
        """
        return self.set_volume(
            SimpleITKImageBuilder()
            .set_volume(volume)
            .set_origin(origin)
            .set_spacing(spacing)
            .set_direction(direction)
            .build()
        )

    def set_plane(self, plan_view: str) -> "CTBuilder":
        """
        Choose the plane each slice lies in.

        Parameters
        ----------
        plan_view : {"AXIAL", "CORONAL", "SAGITTAL"}
            Case-insensitive. Named for the plane each slice lies *in*, so the axis the
            slices stack along is the one **not** in the name: axial stacks along z,
            coronal along y, sagittal along x. Default AXIAL.

        Returns
        -------
        CTBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``plan_view`` is not one of the three.

        Notes
        -----
        Reslicing does not resample. A coronal series from an anisotropic axial volume
        inherits the slice spacing as its in-plane row spacing, so it will look coarse in
        one direction. Resample to isotropic first if that matters.
        """
        plan_view = plan_view.upper()
        if plan_view not in PLANE_AXES:
            raise ValueError(
                f"Invalid plan view: {plan_view!r}; expected one of {sorted(PLANE_AXES)}"
            )
        self.plan_view = plan_view
        return self

    def set_copy_all_attributes(self, copy_all_attributes: bool) -> "CTBuilder":
        """
        Copy every attribute from the reference series instead of building fresh.

        Parameters
        ----------
        copy_all_attributes : bool
            Default False, which writes only the elements the CT Image IOD requires.

        Returns
        -------
        CTBuilder
            self, so calls chain.

        Notes
        -----
        Useful to inherit acquisition parameters -- KVP, exposure, convolution kernel --
        that a fresh dataset would not carry. It also inherits everything else the source
        had, including private tags and any burned-in identifiers, and those inherited
        values now describe a volume they were not measured from. Geometry and pixel
        elements are still overwritten from the volume.
        """
        self.copy_all_attributes = copy_all_attributes
        return self

    def build(self) -> List[Dataset]:
        """
        Slice the volume into a CT Image series.

        Returns
        -------
        list[Dataset]
            One dataset per slice, ordered along the slice normal.

        Raises
        ------
        ValueError
            If no volume has been set.
        """
        if self.volume_image is None:
            raise ValueError("no volume set; call set_volume() before build()")
        return self._build_series(
            self.volume_image, self.plan_view, self.copy_all_attributes
        )

    def _generate_file_meta(self):
        # https://dicom.nema.org/dicom/2013/output/chtml/part04/sect_i.4.html
        # Group length is recomputed on write, so it is not set here.
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        file_meta.MediaStorageSOPInstanceUID = generate_uid(prefix=self.uid_prefix)
        file_meta.ImplementationClassUID = self.uid_prefix + "1"
        file_meta.ImplementationVersionName = "pydicomRT"
        return file_meta

    def _generate_base_dataset(self):
        # Following the DICOM standard for CT Image
        # https://dicom.innolitics.com/ciods/ct-image
        file_meta = self._generate_file_meta()
        ds = FileDataset(
                self._base_file_name,
                dataset={},
                file_meta=file_meta,
                preamble=b"\0" * 128,
            )
        ds = self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        ds = self._add_general_study_information(ds, self.ref_ds)
        ds = self._add_general_series_information(ds, self.ref_ds)
        ds = self._add_frame_of_reference_information(ds, self.ref_ds)
        ds = self._add_general_equipment_information(ds, self.ref_ds)
        ds = self._add_contrast_bolus_information(ds, self.ref_ds)
        return ds

    def _add_patient_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.PatientID = getattr(ref_ds, "PatientID", "")
        ds.PatientName = getattr(ref_ds, "PatientName", "")
        ds.PatientBirthDate = getattr(ref_ds, "PatientBirthDate", "")
        ds.PatientBirthTime = getattr(ref_ds, "PatientBirthTime", "")
        ds.PatientSex = getattr(ref_ds, "PatientSex", "")
        ds.PatientAge = getattr(ref_ds, "PatientAge", "")
        ds.PatientWeight = getattr(ref_ds, "PatientWeight", "")
        ds.PatientPosition = getattr(ref_ds, "PatientPosition", "")
        return ds

    def _add_general_study_information(self, ds: Dataset, ref_ds: Dataset):
        ds.StudyInstanceUID = getattr(ref_ds, "StudyInstanceUID", "")
        ds.StudyID = getattr(ref_ds, "StudyID", "")
        ds.StudyDescription = getattr(ref_ds, "StudyDescription", "")
        ds.StudyDate = getattr(ref_ds, "StudyDate", "")
        ds.StudyTime = getattr(ref_ds, "StudyTime", "")
        ds.AccessionNumber = getattr(ref_ds, "AccessionNumber", "")
        ds.ReferringPhysicianName = getattr(ref_ds, "ReferringPhysicianName", "")
        return ds

    def _add_general_series_information(
        self,
        ds: Dataset,
        ref_ds: Dataset,
        series_number: str = "1",
        series_desc_prefix: str = "CT"
        ):
        dt = datetime.datetime.now()
        ds.Modality = "CT"
        ds.SpecificCharacterSet = "ISO_IR 192"
        ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        ds.SeriesNumber = series_number
        ds.SeriesDescription = f"{series_desc_prefix} {dt.strftime('%Y%m%d%H%M')}"
        ds.SeriesDate = dt.strftime("%Y%m%d")
        ds.SeriesTime = dt.strftime("%H%M%S")
        ds.ContentCreatorName = "pydicomRT"
        ds.ProtocolName = getattr(ref_ds, "ProtocolName", "")
        return ds

    def _add_frame_of_reference_information(self, ds: Dataset, ref_ds: Dataset):
        ds.FrameOfReferenceUID = getattr(ref_ds, "FrameOfReferenceUID", "")
        ds.PositionReferenceIndicator = getattr(ref_ds, "PositionReferenceIndicator", "")
        return ds

    def _add_general_equipment_information(self, ds: Dataset, ref_ds: Dataset):
        ds.Manufacturer = getattr(ref_ds, "Manufacturer", "pydicomRT")
        ds.InstitutionName = getattr(ref_ds, "InstitutionName", "pydicomRT")
        ds.InstitutionAddress = getattr(ref_ds, "InstitutionAddress", "")
        ds.StationName = getattr(ref_ds, "StationName", "")
        ds.InstitutionalDepartmentName = getattr(ref_ds, "InstitutionalDepartmentName", "")
        ds.ManufacturerModelName = getattr(ref_ds, "ManufacturerModelName", "modelv1")
        ds.DeviceSerialNumber = getattr(ref_ds, "DeviceSerialNumber", "")
        ds.SoftwareVersions = getattr(ref_ds, "SoftwareVersions", "v1.0")
        spatial_resolution = getattr(ref_ds, "SpatialResolution", None)
        if spatial_resolution: ds.SpatialResolution = spatial_resolution
        return ds

    def _add_contrast_bolus_information(self, ds: Dataset, ref_ds: Dataset):
        ds.ContrastBolusAgent = getattr(ref_ds, "ContrastBolusAgent", "")
        return ds

    def _build_series(
        self,
        sitk_image: sitk.Image,
        plan_view: str,
        copy_all_attributes: bool,
        ) -> List[Dataset]:
        ds_list = []

        now = datetime.datetime.now()
        if copy_all_attributes:
            base_ds = copy.deepcopy(self.template_ds)
            base_ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        else:
            base_ds = self._generate_base_dataset()
            base_ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
            base_ds.InstanceCreationDate = now.strftime("%Y%m%d")

        base_ds.BitsAllocated = 16
        base_ds.BitsStored = 16
        base_ds.HighBit = 15
        # The stored values below are unsigned; declaring them signed makes every value
        # above 32767 read back negative.
        base_ds.PixelRepresentation = 0
        base_ds.PhotometricInterpretation = "MONOCHROME2"
        base_ds.SamplesPerPixel = 1
        # Air, expressed in the stored (pre-rescale) representation. VR is US here, so it is
        # an integer -- it used to be written as raw bytes.
        base_ds.PixelPaddingValue = 0
        base_ds.ImageType = ["DERIVED", "SECONDARY", plan_view]
        base_ds.KVP = getattr(self.ref_ds, "KVP", 120)
        base_ds.AcquisitionNumber = getattr(self.ref_ds, "AcquisitionNumber", 1)
        base_ds.ContentDate = now.strftime("%Y%m%d")
        base_ds.ContentTime = now.strftime("%H%M%S")
        base_ds.RescaleType = "HU"

        column_axis, row_axis, normal_axis = PLANE_AXES[plan_view]

        origin = np.array(sitk_image.GetOrigin(), dtype=float)
        voxel_spacing = np.array(sitk_image.GetSpacing(), dtype=float)
        # Column j of a SimpleITK direction matrix is the unit vector travelled when index
        # axis j increases.
        direction = np.array(sitk_image.GetDirection(), dtype=float).reshape(3, 3)
        image_size = sitk_image.GetSize()
        image_array = sitk.GetArrayFromImage(sitk_image)

        slice_number = image_size[normal_axis]
        rows = image_size[row_axis]
        columns = image_size[column_axis]

        # ImageOrientationPatient is the row direction (column index increasing) followed by
        # the column direction (row index increasing).
        image_orientation = (
            direction[:, column_axis].tolist() + direction[:, row_axis].tolist()
        )
        # PixelSpacing is ordered [between rows, between columns] -- the reverse of the
        # index order, which is why it is not simply (spacing[column], spacing[row]).
        pixel_spacing = [float(voxel_spacing[row_axis]), float(voxel_spacing[column_axis])]
        thickness = float(voxel_spacing[normal_axis])
        slice_step = direction[:, normal_axis] * thickness

        # image_array is (z, y, x); the SimpleITK axis n is numpy axis 2 - n.
        numpy_slice_axis = 2 - normal_axis

        for slice_idx in range(slice_number):
            slice_ds = copy.deepcopy(base_ds)
            slice_ds.InstanceCreationTime = datetime.datetime.now().strftime("%H%M%S.%f")

            # Generate a new SOP Instance UID for each slice
            instance_uid = generate_uid(prefix=self.uid_prefix)
            slice_ds.file_meta.MediaStorageSOPInstanceUID = instance_uid
            slice_ds.SOPInstanceUID = instance_uid

            slice_2d = np.take(image_array, slice_idx, axis=numpy_slice_axis)

            # Stored value = HU - RESCALE_INTERCEPT, clipped to what uint16 can hold. The
            # upper bound used to be 65536 *before* adding the offset, which overflowed.
            stored = np.clip(
                slice_2d, RESCALE_INTERCEPT, _STORED_MAX + RESCALE_INTERCEPT
            ) - RESCALE_INTERCEPT
            slice_ds.PixelData = np.rint(stored).astype(np.uint16).tobytes()
            slice_ds.RescaleIntercept = RESCALE_INTERCEPT
            slice_ds.RescaleSlope = 1

            slice_ds.ImagePositionPatient = (origin + slice_idx * slice_step).tolist()
            slice_ds.ImageOrientationPatient = image_orientation
            slice_ds.PixelSpacing = pixel_spacing
            slice_ds.SliceThickness = thickness
            slice_ds.SpacingBetweenSlices = thickness
            slice_ds.Rows = rows
            slice_ds.Columns = columns
            slice_ds.InstanceNumber = slice_idx + 1
            ds_list.append(slice_ds)

        return ds_list
