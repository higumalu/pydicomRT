import os
import numpy as np
import SimpleITK as sitk

from datetime import datetime
from pydicom.tag import Tag
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid
from pydicom.uid import ImplicitVRLittleEndian


DICOM_UID_PREFIX = os.getenv('DICOM_UID_PREFIX', "1.2.826.0.1.3680043.8.498.")  # PYDICOM_ROOT_UID

def generate_base_dataset() -> FileDataset:
    file_name = "bear_dose"
    file_meta = get_file_meta()
    ds = FileDataset(file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
    add_required_elements_to_ds(ds)
    add_sequence_lists_to_ds(ds)
    return ds

def get_file_meta() -> FileMetaDataset:
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationGroupLength = 190
    file_meta.FileMetaInformationVersion = b"\x00\x01"
    file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2"
    file_meta.MediaStorageSOPInstanceUID = (
        generate_uid(prefix=DICOM_UID_PREFIX)
    )
    file_meta.ImplementationClassUID = DICOM_UID_PREFIX + "1"
    # TODO file_meta.ImplementationVersionName = ""
    return file_meta

def add_required_elements_to_ds(ds: FileDataset):
    dt = datetime.now()
    # Append data elements required by the DICOM standarad
    ds.SpecificCharacterSet = "ISO_IR 192"
    ds.InstanceCreationDate = dt.strftime("%Y%m%d")
    ds.InstanceCreationTime = dt.strftime("%H%M%S")
    ds.Modality = "RTDOSE"
    ds.Manufacturer = "higumalu"
    ds.ManufacturerModelName = "modelv1"
    ds.InstitutionName = "higumalu"
    # Set the transfer syntax
    ds.is_little_endian = True
    ds.is_implicit_VR = True
    # Set values already defined in the file meta
    ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID

    ds.ApprovalStatus = "UNAPPROVED"

def add_sequence_lists_to_ds(ds: FileDataset):
    ds.ReferencedRTPlanSequence = Sequence()

# ---------------------------------------------------------------------------------------- #

def add_patient_information(ds: FileDataset, reference_ds: Dataset):
    ds.PatientName = getattr(reference_ds, "PatientName", "Unknown")
    ds.PatientID = getattr(reference_ds, "PatientID", "Unknown")
    ds.PatientBirthDate = getattr(reference_ds, "PatientBirthDate", "")
    ds.PatientSex = getattr(reference_ds, "PatientSex", "")
    ds.PatientAge = getattr(reference_ds, "PatientAge", "")
    ds.PatientSize = getattr(reference_ds, "PatientSize", "")
    ds.PatientWeight = getattr(reference_ds, "PatientWeight", "")
    return ds

def add_study_information(ds: FileDataset, reference_ds: Dataset):
    dt = datetime.now()
    ds.StudyInstanceUID = getattr(reference_ds, "StudyInstanceUID", generate_uid(prefix=DICOM_UID_PREFIX))
    ds.StudyDate = getattr(reference_ds, "StudyDate", dt.strftime("%Y%m%d"))
    ds.StudyTime = getattr(reference_ds, "StudyTime", dt.strftime("%H%M%S"))
    ds.StudyID = getattr(reference_ds, "StudyID", "")
    ds.AccessionNumber = getattr(reference_ds, "AccessionNumber", "")
    ds.ReferringPhysicianName = getattr(reference_ds, "ReferringPhysicianName", "")
    return ds

def add_series_information(ds: FileDataset, reference_ds: Dataset, series_number: int = 1):
    dt = datetime.now()
    ds.Modality = "RTDOSE"
    ds.SeriesInstanceUID = generate_uid(prefix=DICOM_UID_PREFIX)
    ds.SeriesDate = dt.strftime("%Y%m%d")
    ds.SeriesTime = dt.strftime("%H%M%S")
    ds.SeriesDescription = getattr(reference_ds, "SeriesDescription", "") + "_dose" + dt.strftime("%Y%m%d%H%M%S")
    ds.SeriesNumber = getattr(reference_ds, "SeriesNumber", series_number)
    ds.OperatorsName = getattr(reference_ds, "OperatorsName", "")
    return ds

def add_frame_of_reference_information(ds: FileDataset, reference_ds: Dataset):
    ds.FrameOfReferenceUID = getattr(reference_ds, "FrameOfReferenceUID", generate_uid(prefix=DICOM_UID_PREFIX))
    ds.PositionReferenceIndicator = getattr(reference_ds, "PositionReferenceIndicator", "")
    return ds

# ---------------------------------------------------------------------------------------- #

def cp_information_from_ds(
    ds: FileDataset,
    reference_ds: Dataset,
    ):
    ds = add_patient_information(ds, reference_ds)
    ds = add_study_information(ds, reference_ds)
    ds = add_series_information(ds, reference_ds)
    ds = add_frame_of_reference_information(ds, reference_ds)
    return ds

# ---------------------------------------------------------------------------------------- #

def add_dose_grid_to_ds(
    ds: FileDataset,
    dose_sitk_image: sitk.Image,
    ):
    scaling = 1e-7
    array = sitk.GetArrayFromImage(dose_sitk_image)
    size = list(dose_sitk_image.GetSize())
    spacing = list(dose_sitk_image.GetSpacing())
    direction = list(dose_sitk_image.GetDirection())
    origin = list(dose_sitk_image.GetOrigin())

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = True
    ds.DoseUnits = "GY"
    ds.DoseType = "EFFECTIVE"
    ds.DoseSummationType = "PLAN"
    ds.NumberOfFrames = size[2]
    ds.FrameIncrementPointer = Tag(0x3004, 0x000C)

    ds.ImagePositionPatient = origin
    ds.ImageOrientationPatient = direction
    ds.PixelSpacing = spacing[0:2]
    ds.SliceThickness = spacing[2]
    ds.Columns = size[0]
    ds.Rows = size[1]
    ds.GridFrameOffsetVector = [ov * spacing[2] for ov in range(size[2])]
    ds.DoseGridScaling = scaling
    ds.PixelData = (array // scaling).astype(np.uint32).tobytes()
    return ds


# -------------------------------- example --------------------------------------------------- #

if __name__ == "__main__":
    import pydicom
    import SimpleITK as sitk
    dose_sitk_image = sitk.ReadImage("example/data/DF_001/RD_003/CT_001_dose.dcm")
    ds = pydicom.dcmread("example/data/DF_001/RD_003/CT_001.dcm")

    dose_ds = generate_base_dataset()
    dose_ds = cp_information_from_ds(dose_ds, ds)
    dose_ds = add_dose_grid_to_ds(dose_ds, dose_sitk_image)
    dose_ds.save_as("example/data/DF_001/RD_003/dose.dcm")
