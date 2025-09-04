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
        self.registration_dataset_list = []

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
        ref_ds: Dataset,
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
        ds.FrameOfReferenceUID = getattr(ref_ds, "FrameOfReferenceUID", "")
        ds.PositionReferenceIndicator = getattr(ref_ds, "PositionReferenceIndicator", "")
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

    def add_rigid_registration(self, moving_ds_list, rigid_transform_matrix):
        rigid_registration_block = Dataset()
        # TODO
        # ReferencedImageSequence   (moving)
        # FrameOfReferenceUID   (moving)
        # MatrixRegistrationSequence
        #     MatrixSequence
        #         FrameOfReferenceTransformationMatrixType
        #         FrameOfReferenceTransformationMatrix    (moving to fixed)
        #     RegistrationTypeCodeSequence
        # https://dicom.nema.org/medical/dicom/2016a/output/chtml/part16/chapter_D.html#DCM_125024
        #         CodeValue     125024
        #         CodingSchemeDesignator     DCM
        #         CodingSchemeVersion       01
        #         CodeMeaning     Image Content-based Alignment

        frame_of_reference_uid = getattr(moving_ds_list[0], "FrameOfReferenceUID", "")
        rigid_registration_block.FrameOfReferenceUID = frame_of_reference_uid
        rigid_registration_block.ReferencedImageSequence = Sequence()
        for moving_ds in moving_ds_list:
            referenced_image = Dataset()
            referenced_image.ReferencedSOPClassUID = moving_ds.SOPClassUID
            referenced_image.ReferencedSOPInstanceUID = moving_ds.SOPInstanceUID
            rigid_registration_block.ReferencedImageSequence.append(referenced_image)

        rigid_registration_block.MatrixRegistrationSequence = Sequence()
        matrix_registration_sequence_block = Dataset()
        matrix_registration_sequence_block.MatrixSequence = Sequence()
        matrix_sequence = Sequence()
        matrix_sequence_block = Dataset()
        matrix_sequence_block.FrameOfReferenceTransformationMatrixType = "RIGID"
        matrix_sequence_block.FrameOfReferenceTransformationMatrix = rigid_transform_matrix
        matrix_sequence.append(matrix_sequence_block)
        matrix_registration_sequence_block.MatrixSequence = matrix_sequence
        rigid_registration_block.MatrixRegistrationSequence.append(matrix_registration_sequence_block)

        rigid_registration_block.RegistrationTypeCodeSequence = Sequence()
        registration_type_code = Dataset()
        registration_type_code.CodeValue = "125024"
        registration_type_code.CodingSchemeDesignator = "DCM"
        registration_type_code.CodingSchemeVersion = "01"
        registration_type_code.CodeMeaning = "Image Content-based Alignment"
        rigid_registration_block.RegistrationTypeCodeSequence.append(registration_type_code)

        self.registration_dataset_list.append(rigid_registration_block)
        return rigid_registration_block

    def build(self):
        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(ds, self.ref_ds)
        ds.RegistrationSequence = Sequence()
        for reg_ds in self.registration_dataset_list:
            ds.RegistrationSequence.append(reg_ds)
        return ds


# --------------------- Deformable Spatial Registration Builder ---------------------

class DeformableSpatialRegistrationBuilder(BaseRegistrationBuilder):
    def __init__(self, fixed_ds_list):
        super().__init__(fixed_ds_list)

    @property
    def _sop_class_uid(self) -> str:
        # Deformable Spatial Registration SOP Class UID
        return "1.2.840.10008.5.1.4.1.1.66.3"

    @property
    def _base_file_name(self) -> str:
        return "deformable_spatial_registration"

    def add_deformable_registration(
        self,
        moving_ds_list,
        vectorial_field,
        pre_transform,
        post_transform):
        deformable_registration_block = Dataset()

        # TODO
        # ReferencedImageSequence   (moving)
        # SourceFrameOfReferenceUID   (moving)
        # DeformableRegistrationGridSequence
        #     ImagePositionPatient
        #     ImageOrientationPatient
        #     GridDimensions
        #     GridResolution
        #     VectorGridData
        # PreDeformationMatrixRegistrationSequence
        #     FrameOfReferenceTransformationMatrixType
        #     FrameOfReferenceTransformationMatrix
        # PostDeformationMatrixRegistrationSequence
        #     FrameOfReferenceTransformationMatrixType
        #     FrameOfReferenceTransformationMatrix
        # RegistrationTypeCodeSequence
        # https://dicom.nema.org/medical/dicom/2016a/output/chtml/part16/chapter_D.html#DCM_125024
        #     CodeValue     125024
        #     CodingSchemeDesignator     DCM
        #     CodingSchemeVersion       01
        #     CodeMeaning     Image Content-based Alignment
        deformable_registration_block.ReferencedImageSequence = Sequence()
        deformable_registration_block.DeformableRegistrationGridSequence = Sequence()
        deformable_registration_block.PreDeformationMatrixRegistrationSequence = Sequence()
        deformable_registration_block.PostDeformationMatrixRegistrationSequence = Sequence()
        deformable_registration_block.RegistrationTypeCodeSequence = Sequence()
        deformable_registration_block.SourceFrameOfReferenceUID = getattr(moving_ds_list[0], "FrameOfReferenceUID", "")

        for moving_ds in moving_ds_list:
            referenced_image = Dataset()
            referenced_image.ReferencedSOPClassUID = moving_ds.SOPClassUID
            referenced_image.ReferencedSOPInstanceUID = moving_ds.SOPInstanceUID
            deformable_registration_block.ReferencedImageSequence.append(referenced_image)

        deformable_registration_grid = Dataset()
        deformable_registration_grid.ImagePositionPatient = vectorial_field.image_position_patient
        deformable_registration_grid.ImageOrientationPatient = vectorial_field.image_orientation_patient
        deformable_registration_grid.GridDimensions = vectorial_field.grid_dimensions
        deformable_registration_grid.GridResolution = vectorial_field.grid_resolution
        deformable_registration_grid.VectorGridData = vectorial_field.vector_grid_data
        deformable_registration_block.DeformableRegistrationGridSequence.append(deformable_registration_grid)

        pre_deformation_matrix_registration = Dataset()
        pre_deformation_matrix_registration.FrameOfReferenceTransformationMatrixType = "RIGID"
        pre_deformation_matrix_registration.FrameOfReferenceTransformationMatrix = pre_transform
        deformable_registration_block.PreDeformationMatrixRegistrationSequence.append(pre_deformation_matrix_registration)

        post_deformation_matrix_registration = Dataset()
        post_deformation_matrix_registration.FrameOfReferenceTransformationMatrixType = "RIGID"
        post_deformation_matrix_registration.FrameOfReferenceTransformationMatrix = post_transform
        deformable_registration_block.PostDeformationMatrixRegistrationSequence.append(post_deformation_matrix_registration)

        self.registration_dataset_list.append(deformable_registration_block)
        pass

    def build(self):
        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(ds, series_desc_prefix="Deformable Spatial Registration")
        ds.DeformableRegistrationSequence = Sequence()
        for reg_ds in self.registration_dataset_list:
            ds.DeformableRegistrationSequence.append(reg_ds)
        return ds


if __name__ == "__main__":
    from pydicomrt.utils import load_sorted_image_series
    fixed_dcm_path = "example/data/CT_001"
    moving_dcm_path = "example/data/CT_001"
    rigid_transform_matrix = [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
    ]

    fixed_ds_list = load_sorted_image_series(fixed_dcm_path)
    moving_ds_list = load_sorted_image_series(moving_dcm_path)

    spatial_registration_builder = SpatialRegistrationBuilder(fixed_ds_list)
    spatial_registration_builder.add_rigid_registration(moving_ds_list, rigid_transform_matrix)
    ds = spatial_registration_builder.build()
    ds.save_as("./CT_001_rigid_registration.dcm")
