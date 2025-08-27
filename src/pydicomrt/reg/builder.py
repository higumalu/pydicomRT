import datetime
from abc import ABC, abstractmethod
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, ImplicitVRLittleEndian

# --------------------- Basic Registration Builder ---------------------
class BaseRegistrationBuilder(ABC):
    def __init__(self, fixed_ds_list):
        assert len(fixed_ds_list) > 0, "fixed_ds_list must be non-empty"
        self.fixed_ds_list = fixed_ds_list
        self.uid_prefix = None
        self.ref_ds = fixed_ds_list[0]

    def set_uid_prefix(self, uid_prefix: str):
        self.uid_prefix = uid_prefix

    @property
    @abstractmethod
    def _sop_class_uid(self) -> str:
        raise NotImplementedError

    @property
    @abstractmethod
    def _base_file_name(self) -> str:
        raise NotImplementedError

    def _generate_file_meta(self):
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationGroupLength = 202
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = self._sop_class_uid
        file_meta.MediaStorageSOPInstanceUID = generate_uid(prefix=self.uid_prefix)
        file_meta.ImplementationClassUID = (self.uid_prefix or "") + "1"
        file_meta.ImplementationVersionName = "pydicomRT"
        return file_meta

    def _generate_base_dataset(self):
        file_meta = self._generate_file_meta()
        ds = FileDataset(self._base_file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
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
        ds.is_little_endian = True
        ds.is_implicit_VR = True
        ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID

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

    def _add_study_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.StudyInstanceUID = getattr(ref_ds, "StudyInstanceUID", "")
        ds.StudyID = getattr(ref_ds, "StudyID", "")
        ds.StudyDescription = getattr(ref_ds, "StudyDescription", "")
        ds.StudyDate = getattr(ref_ds, "StudyDate", "")
        ds.StudyTime = getattr(ref_ds, "StudyTime", "")
        ds.AccessionNumber = getattr(ref_ds, "AccessionNumber", "")
        ds.ReferringPhysicianName = getattr(ref_ds, "ReferringPhysicianName", "")
        return ds

    def _add_series_information(
        self,
        ds: Dataset,
        series_number: str = "50",
        series_desc_prefix: str = "Spatial Registration"
        ):
        dt = datetime.datetime.now()
        ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        ds.SeriesNumber = series_number
        ds.SeriesDescription = f"{series_desc_prefix} {dt.strftime('%Y%m%d%H%M')}"
        ds.SeriesDate = dt.strftime("%Y%m%d")
        ds.SeriesTime = dt.strftime("%H%M%S")
        ds.ContentLabel = "SPATIAL REGISTRATION"
        ds.ContentDescription = "Spatial Registration"
        ds.ContentCreatorName = "pydicomRT"
        return ds


# --------------------- Spatial Registration Builder ---------------------

class SpatialRegistrationBuilder(BaseRegistrationBuilder):
    def __init__(self, fixed_ds_list):
        super().__init__(fixed_ds_list)

    @property
    def _sop_class_uid(self) -> str:
        # Spatial Registration (rigid/affine) SOP Class UID
        return "1.2.840.10008.5.1.4.1.1.66.1"

    @property
    def _base_file_name(self) -> str:
        return "spatial_registration"

    def build(self):
        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(ds)
        # TODO: RegistrationSequence
        ds.RegistrationSequence = Sequence()

        return ds

    def add_rigid_registration(self, moving_ds_list, rigid_transform_matrix):
        pass


# --------------------- Deformable Spatial Registration Builder ---------------------

class DeformableSpatialRegistrationBuilder(BaseRegistrationBuilder):
    def __init__(self, fixed_ds_list, moving_ds_list, pre_transform, vectorial_field, post_transform):
        super().__init__(fixed_ds_list)
        self.moving_ds_list = moving_ds_list
        self.pre_transform = pre_transform
        self.vectorial_field = vectorial_field
        self.post_transform = post_transform

    @property
    def _sop_class_uid(self) -> str:
        # Deformable Spatial Registration SOP Class UID
        return "1.2.840.10008.5.1.4.1.1.66.3"

    @property
    def _base_file_name(self) -> str:
        return "deformable_spatial_registration"

    def build(self):
        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(ds, series_desc_prefix="Deformable Spatial Registration")
        # TODO: DeformableRegistrationSequence
        return ds
