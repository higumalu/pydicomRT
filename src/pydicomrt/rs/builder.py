import os
import datetime

from pydicom.uid import generate_uid
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ImplicitVRLittleEndian


DICOM_UID_PREFIX = os.getenv('DICOM_UID_PREFIX', "1.2.826.0.1.3680043.8.498.")  # PYDICOM_ROOT_UID

def create_rtstruct_dataset(series_data) -> FileDataset:
    ds = generate_base_dataset()
    add_study_and_series_information(ds, series_data)
    add_patient_information(ds, series_data)
    add_refd_frame_of_ref_sequence(ds, series_data)
    add_rs_series_information(ds, series_data)
    return ds

def generate_base_dataset() -> FileDataset:
    file_name = "rs_test"
    file_meta = get_file_meta()
    ds = FileDataset(file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
    add_required_elements_to_ds(ds)
    add_sequence_lists_to_ds(ds)
    return ds

def get_file_meta() -> FileMetaDataset:
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationGroupLength = 202
    file_meta.FileMetaInformationVersion = b"\x00\x01"
    file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
    file_meta.MediaStorageSOPInstanceUID = (
        generate_uid(prefix=DICOM_UID_PREFIX)
    )
    file_meta.ImplementationClassUID = DICOM_UID_PREFIX + "1"
    # TODO file_meta.ImplementationVersionName = ""
    return file_meta

def add_required_elements_to_ds(ds: FileDataset):
    dt = datetime.datetime.now()
    # Append data elements required by the DICOM standarad
    ds.SpecificCharacterSet = "ISO_IR 192"
    ds.InstanceCreationDate = dt.strftime("%Y%m%d")
    ds.InstanceCreationTime = dt.strftime("%H%M%S")
    ds.StructureSetLabel = "RS_" + dt.strftime("%Y%m%d")
    ds.StructureSetDate = dt.strftime("%Y%m%d")
    ds.StructureSetTime = dt.strftime("%H%M%S.%f")
    ds.Modality = "RTSTRUCT"
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
    ds.StructureSetROISequence = Sequence()
    ds.ROIContourSequence = Sequence()
    ds.RTROIObservationsSequence = Sequence()

# ---------------------------------------------------------------------------------------- #

def add_study_and_series_information(ds: FileDataset, series_data):
    reference_ds = series_data[0]  # All elements in series should have the same data
    ds.StudyDate = reference_ds.StudyDate
    ds.SeriesDate = getattr(reference_ds, "SeriesDate", "")
    ds.StudyTime = reference_ds.StudyTime
    ds.SeriesTime = getattr(reference_ds, "SeriesTime", "")
    ds.StudyDescription = getattr(reference_ds, "StudyDescription", "")
    ds.SeriesDescription = getattr(reference_ds, "SeriesDescription", "")
    ds.StudyInstanceUID = reference_ds.StudyInstanceUID
    ds.SeriesInstanceUID = generate_uid(prefix=DICOM_UID_PREFIX)  # TODO: find out if random generation is ok
    ds.StudyID = reference_ds.StudyID
    ds.AccessionNumber = getattr(reference_ds, "AccessionNumber", "")
    ds.SeriesNumber = "1"  # TODO: find out if we can just use 1 (Should be fine since its a new series)
    ds.OperatorsName = getattr(reference_ds, "OperatorsName", "")

def add_patient_information(ds: FileDataset, series_data):
    reference_ds = series_data[0]  # All elements in series should have the same data
    ds.PatientName = getattr(reference_ds, "PatientName", "")
    ds.PatientID = getattr(reference_ds, "PatientID", "")
    ds.PatientBirthDate = getattr(reference_ds, "PatientBirthDate", "")
    ds.PatientSex = getattr(reference_ds, "PatientSex", "")
    ds.PatientAge = getattr(reference_ds, "PatientAge", "")
    ds.PatientSize = getattr(reference_ds, "PatientSize", "")
    ds.PatientWeight = getattr(reference_ds, "PatientWeight", "")

def add_refd_frame_of_ref_sequence(ds: FileDataset, series_data):
    refd_frame_of_ref = Dataset()
    refd_frame_of_ref.FrameOfReferenceUID = getattr(series_data[0], 'FrameOfReferenceUID', generate_uid(prefix=DICOM_UID_PREFIX))
    refd_frame_of_ref.RTReferencedStudySequence = create_frame_of_ref_study_sequence(series_data)

    ds.ReferencedFrameOfReferenceSequence = Sequence()
    ds.ReferencedFrameOfReferenceSequence.append(refd_frame_of_ref)
    ds.FrameOfReferenceUID = getattr(series_data[0], 'FrameOfReferenceUID', generate_uid(prefix=DICOM_UID_PREFIX))

def create_frame_of_ref_study_sequence(series_data) -> Sequence:
    reference_ds = series_data[0]
    rt_refd_series = Dataset()
    rt_refd_series.SeriesInstanceUID = reference_ds.SeriesInstanceUID
    rt_refd_series.ContourImageSequence = create_contour_image_sequence(series_data)
    rt_refd_series_sequence = Sequence()
    rt_refd_series_sequence.append(rt_refd_series)
    rt_refd_study = Dataset()
    rt_refd_study.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"
    rt_refd_study.ReferencedSOPInstanceUID = reference_ds.StudyInstanceUID
    rt_refd_study.RTReferencedSeriesSequence = rt_refd_series_sequence
    rt_refd_study_sequence = Sequence()
    rt_refd_study_sequence.append(rt_refd_study)
    return rt_refd_study_sequence

def create_contour_image_sequence(series_data) -> Sequence:
    contour_image_sequence = Sequence()
    for series in series_data:
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = series.SOPClassUID
        contour_image.ReferencedSOPInstanceUID = series.SOPInstanceUID
        contour_image_sequence.append(contour_image)
    return contour_image_sequence

def add_rs_series_information(ds: FileDataset, series_data):
    dt = datetime.datetime.now()
    ds.StructureSetDescription = dt.strftime("%Y%m%d") + "MBB"
