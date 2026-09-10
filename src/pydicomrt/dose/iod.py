from pydicom.uid import UID

# RT Dose CIOD
RTDOSE_IOD = {
    "PatientID": {},
    "PatientName": {},
    "PatientBirthDate": {},
    "PatientSex": {},
    "StudyInstanceUID": {},
    "StudyID": {},
    "StudyDate": {},
    "StudyTime": {},
    "AccessionNumber": {},
    "SeriesInstanceUID": {"type": UID},
    "SeriesNumber": {},
    "Modality": {"value": "RTDOSE"},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID},
    "SOPClassUID": {"type": UID},
    "SOPInstanceUID": {"type": UID},
    "InstanceNumber": {},

    "NumberOfFrames": {},
    "Rows": {},
    "Columns": {},
    "PixelSpacing": {},
    "ImageOrientationPatient": {},
    "ImagePositionPatient": {},     # if value == (3004,000C) GridFrameOffsetVector is required
    "SliceThickness": {"optional": True},

    "BitsAllocated": {},
    "BitsStored": {},
    "HighBit": {},
    "PixelRepresentation": {},

    "DoseUnits": {"value": "GY"},
    "DoseType": {},
    "DoseSummationType": {},
    "GridFrameOffsetVector": {},

    "PixelData": {},
}
