from pydicom.uid import UID

# Spatial Registration CIOD
SPATIAL_REGSITRATION_IOD = {
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
    "Modality": {"value": "REG"},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID},
    "SOPClassUID": {"type": UID},
    "SOPInstanceUID": {"type": UID},
    "InstanceNumber": {},
    "ContentLabel": {},
    "ContentDescription": {},
    "RegistrationSequence": {
        "FrameOfReferenceUID": {"type": UID},
        "submap": {
            "MatrixRegistrationSequence": {
                "submap": {
                    "MatrixSequence": {
                        "FrameOfReferenceTransformationMatrixType": {},
                        "FrameOfReferenceTransformationMatrix": {},
                    }
                }
            },
        }
    },

}
