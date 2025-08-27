import os
import datetime

from pydicom.uid import generate_uid
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ImplicitVRLittleEndian


class SpatialRegistrationBuilder:
    def __init__(
        self,
        fixed_ds_list,
        moving_ds_list,
        transform_matrix):
        assert len(fixed_ds_list) > 0, "fixed_ds_list must be non-empty"
        assert len(moving_ds_list) > 1, "moving_ds_list must be non-empty"

        self.fixed_ds_list = fixed_ds_list
        self.moving_ds_list = moving_ds_list
        self.transform_matrix = transform_matrix
        self.uid_prefix = None
        self.ref_ds = moving_ds_list[0]

    def set_uid_prefix(self, uid_prefix: str):
        self.uid_prefix = uid_prefix

    def build(self):
        pass

    def _generate_file_meta(self):
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationGroupLength = 202
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.66.1"
        file_meta.MediaStorageSOPInstanceUID = (
            generate_uid(prefix=self.uid_prefix)
        )
        file_meta.ImplementationClassUID = self.uid_prefix + "1"
        file_meta.ImplementationVersionName = "pydicomRT"
        return file_meta

    def _generate_base_dataset(self):
        file_name = "spatial_registration"
        file_meta = self._generate_file_meta()
        ds = FileDataset(file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
        return ds

    def _add_required_elements(self, ds: Dataset):
        dt = datetime.datetime.now()
        ds.SpecificCharacterSet = "ISO_IR 192"
        ds.InstanceCreationDate = dt.strftime("%Y%m%d")
        ds.InstanceCreationTime = dt.strftime("%H%M%S")
        ds.Modality = "REG"
        ds.Manufacturer = ""
        ds.ManufacturerModelName = ""
        ds.InstitutionName = ""
        # Set the transfer syntax
        ds.is_little_endian = True
        ds.is_implicit_VR = True
        # Set values already defined in the file meta
        ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID

        ds.RegistrationSequence = Sequence()
        return ds

    def _add_patient_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.PatientID = ref_ds.getattr("PatientID", "")
        ds.PatientName = ref_ds.getattr("PatientName", "")
        ds.PatientBirthDate = ref_ds.getattr("PatientBirthDate", "")
        ds.PatientBirthTime = ref_ds.getattr("PatientBirthTime", "")
        ds.PatientSex = ref_ds.getattr("PatientSex", "")
        ds.PatientAge = ref_ds.getattr("PatientAge", "")
        ds.PatientWeight = ref_ds.getattr("PatientWeight", "")
        ds.PatientPosition = ref_ds.getattr("PatientPosition", "")
        return ds

    def _add_study_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.StudyInstanceUID = ref_ds.getattr("StudyInstanceUID", "")
        ds.StudyID = ref_ds.getattr("StudyID", "")
        ds.StudyDescription = ref_ds.getattr("StudyDescription", "")
        ds.StudyDate = ref_ds.getattr("StudyDate", "")
        ds.StudyTime = ref_ds.getattr("StudyTime", "")
        ds.AccessionNumber = ref_ds.getattr("AccessionNumber", "")
        ds.ReferringPhysicianName = ref_ds.getattr("ReferringPhysicianName", "")
        return ds

    def _add_series_information(self, ds: Dataset):
        dt = datetime.datetime.now()
        ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        ds.SeriesNumber = "50"
        ds.SeriesDescription = "Spatial Registration " + dt.strftime("%Y%m%d%H%M")
        ds.SeriesDate = dt.strftime("%Y%m%d")
        ds.SeriesTime = dt.strftime("%H%M%S")
        ds.ContentLabel = "SPATIAL REGISTRATION"
        ds.ContentDescription = "Spatial Registration"
        ds.ContentCreatorName = "pydicomRT"
        return ds

###################################################################################

class DeformableSpatialRegistrationBuilder:
    def __init__(
        self,
        fixed_ds_list,
        moving_ds_list,
        pre_transform,
        vectorial_field,
        post_transform
        ):
        self.fixed_ds_list = fixed_ds_list
        self.moving_ds_list = moving_ds_list
        self.pre_transform = pre_transform
        self.vectorial_field = vectorial_field
        self.post_transform = post_transform
        self.uid_prefix = None

    def set_uid_prefix(self, uid_prefix: str):
        self.uid_prefix = uid_prefix

    def build(self):
        pass

    def _generate_file_meta(self):
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationGroupLength = 202
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.66.3"
        file_meta.MediaStorageSOPInstanceUID = (
            generate_uid(prefix=self.uid_prefix)
        )
        file_meta.ImplementationClassUID = self.uid_prefix + "1"
        file_meta.ImplementationVersionName = "pydicomRT"
        return file_meta

    def _generate_base_dataset(self):
        file_name = "deformable_spatial_registration"
        file_meta = self._generate_file_meta()
        ds = FileDataset(file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
        return ds
