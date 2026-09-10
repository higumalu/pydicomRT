from pydicom.uid import UID


CT_IMAGE_IOD = {
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
    "Modality": {"value": "CT"},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID},
    "SOPClassUID": {"type": UID},
    "SOPInstanceUID": {"type": UID},
    "InstanceNumber": {},

    "SliceThickness": {},
    "ImagePositionPatient": {},
    "ImageOrientationPatient": {},
    "PixelSpacing": {},

    "Rows": {},
    "Columns": {},

    "ImageType": {},
    "BitsAllocated": {},
    "BitsStored": {},
    "HighBit": {},
    "RescaleIntercept": {},
    "RescaleSlope": {},

    "PixelData": {},

}
