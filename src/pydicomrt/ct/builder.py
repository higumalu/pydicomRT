"""
CT Image Builder
Date: 2026-02-14
Author: higumalu
"""

from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

from typing import Tuple, List
from abc import ABC, abstractmethod
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, ImplicitVRLittleEndian

import copy
import datetime
import numpy as np
import SimpleITK as sitk

class CTBuilder(ABC):
    def __init__(self, referenced_ds_list: List[Dataset]):
        self.referenced_ds_list = referenced_ds_list
        self.uid_prefix = "1.2.826.0.1.3680043.8.498."  # PYDICOM_ROOT_UID
        self.ref_ds = referenced_ds_list[0]
        self._base_file_name = "ct_image"
        self.template_ds = copy.deepcopy(referenced_ds_list[0])

    def set_uid_prefix(self, uid_prefix: str):
        self.uid_prefix = uid_prefix

    def _generate_file_meta(self):
        # DICOM File Meta Information Group Length: 202
        # https://dicom.nema.org/dicom/2013/output/chtml/part04/sect_i.4.html
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationGroupLength = 202
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        file_meta.MediaStorageSOPInstanceUID = generate_uid()
        file_meta.ImplementationClassUID = "1.2.826.0.1.3680043.8.498.1"
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
                is_implicit_VR=True,
                is_little_endian=True
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
        ds.PixelPaddingValue = int(getattr(ref_ds, "PixelPaddingValue", -1024)).to_bytes(2, byteorder='little', signed=True)
        return ds

    def _add_contrast_bolus_information(self, ds: Dataset, ref_ds: Dataset):
        ds.ContrastBolusAgent = getattr(ref_ds, "ContrastBolusAgent", "")
        return ds

    def build_from_sitk_image(
        self, 
        sitk_image: sitk.Image,
        plan_view: str = "AXIAL",
        cp_all_attr: bool = False,
        ) -> List[Dataset]:
        self.plan_view = plan_view
        ds_list = []

        if cp_all_attr:
            base_ds = copy.deepcopy(self.template_ds)
            base_ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        else:
            base_ds = self._generate_base_dataset()
            base_ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
            base_ds.InstanceCreationDate = datetime.datetime.now().strftime("%Y%m%d")

        base_ds.BitsAllocated = 16
        base_ds.BitsStored = 16
        base_ds.HighBit = 15
        base_ds.PixelRepresentation = 1
        base_ds.PhotometricInterpretation = "MONOCHROME2"
        base_ds.SamplesPerPixel = 1
        base_ds.PixelPaddingValue = int(-1024).to_bytes(2, byteorder='little', signed=True)
        base_ds.ImageType = [f"DERIVED", f"SECONDARY", f"{plan_view.upper()}"]
        base_ds.KVP = getattr(self.ref_ds, "KVP", 120)
        base_ds.AcquisitionNumber = 51


        origin = list(sitk_image.GetOrigin())
        voxel_spacing = list(sitk_image.GetSpacing())
        image_direction = list(sitk_image.GetDirection())
        image_size = list(sitk_image.GetSize())
        image_array = sitk.GetArrayFromImage(sitk_image)

        match plan_view:
            case "AXIAL":
                slice_number = image_size[2]
                row = image_size[0]
                col = image_size[1]
                pixel_spacing = [voxel_spacing[0], voxel_spacing[1]]
                thickness = voxel_spacing[2]
                image_orientation = [image_direction[0], image_direction[1], image_direction[2], image_direction[3], image_direction[4], image_direction[5]]

            case "SAGITTAL":
                slice_number = image_size[1]
                row = image_size[0]
                col = image_size[2]
                pixel_spacing = [voxel_spacing[0], voxel_spacing[2]]
                thickness = voxel_spacing[1]
                # row=axis0, col=axis2 → (d0,d1,d2, d6,d7,d8)
                image_orientation = [image_direction[0], image_direction[1], image_direction[2], image_direction[6], image_direction[7], image_direction[8]]

            case "CORONAL":
                slice_number = image_size[0]
                row = image_size[1]
                col = image_size[2]
                pixel_spacing = [voxel_spacing[1], voxel_spacing[2]]
                thickness = voxel_spacing[0]
                # row=axis1, col=axis2 → (d3,d4,d5, d6,d7,d8)
                image_orientation = [image_direction[3], image_direction[4], image_direction[5], image_direction[6], image_direction[7], image_direction[8]]
            case _:
                raise ValueError(f"Invalid plan view: {plan_view}")

        # image_array from SimpleITK is (z, y, x)
        for slice_idx in range(slice_number):
            slice_ds = copy.deepcopy(base_ds)
            
            slice_ds.InstanceCreationTime = datetime.datetime.now().strftime("%H%M%S.%f")

            # Generate a new SOP Instance UID for each slice
            instance_uid = generate_uid(prefix=self.uid_prefix)
            slice_ds.file_meta.MediaStorageSOPInstanceUID = instance_uid
            slice_ds.SOPInstanceUID = instance_uid

            match plan_view:
                case "AXIAL":
                    slice_2d = image_array[slice_idx, :, :]  # (row, col) = (y, x)
                case "SAGITTAL":
                    slice_2d = image_array[:, slice_idx, :]  # (z, x)
                case "CORONAL":
                    slice_2d = image_array[:, :, slice_idx]  # (z, y)
            
            # Clip the values in slice_2d to the range [-1024, 1024]
            row = slice_2d.shape[0]
            col = slice_2d.shape[1]
            slice_2d = np.clip(slice_2d, -1024, 65536) + 1024
            slice_2d = slice_2d.astype(np.uint16)
            slice_ds.PixelData = slice_2d.tobytes()
            slice_ds.RescaleIntercept = -1024
            slice_ds.RescaleSlope = 1

            match plan_view:
                case "AXIAL":
                    slice_ds.ImagePositionPatient = [origin[0], origin[1], origin[2] + slice_idx * thickness]
                case "SAGITTAL":
                    slice_ds.ImagePositionPatient = [origin[0], origin[1] + slice_idx * thickness, origin[2]]
                case "CORONAL":
                    slice_ds.ImagePositionPatient = [origin[0] + slice_idx * thickness, origin[1], origin[2]]

            slice_ds.ImageOrientationPatient = image_orientation
            slice_ds.PixelSpacing = pixel_spacing
            slice_ds.SliceThickness = thickness
            slice_ds.Rows = row
            slice_ds.Columns = col
            slice_ds.InstanceNumber = slice_idx + 1
            ds_list.append(slice_ds)

        return ds_list

    def build_from_np_array(
        self,
        volume: np.ndarray,
        origin: np.ndarray,
        spacing: np.ndarray,
        direction: np.ndarray,
        plan_view: str = "AXIAL",
        ) -> List[Dataset]:

        image_builder = SimpleITKImageBuilder()
        image_builder.set_volume(volume)
        image_builder.set_origin(origin)
        image_builder.set_spacing(spacing)
        image_builder.set_direction(direction)
        sitk_image = image_builder.build()
        print(sitk_image)
        ds_list = self.build_from_sitk_image(sitk_image, plan_view)

        return ds_list

if __name__ == "__main__":
    from pydicomrt.utils.image_series_loader import load_sorted_image_series
    ct_path = "example/data/Mirror/CT"
    ct_ds_list = load_sorted_image_series(ct_path)
    ct_builder = CTBuilder(ct_ds_list)
    volume = np.random.rand(10, 10, 10) * 1000 - 500
    origin = np.array([0, 0, 0])
    spacing = np.array([1, 1, 1])
    direction = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1]).reshape(3, 3)
    ds_list = ct_builder.build_from_np_array(volume, origin, spacing, direction, plan_view="AXIAL")
    for ds in ds_list:
        ds.save_as(f"example/data/Mirror/CT_new/ct_image_{ds.InstanceNumber}.dcm")