from pydicom.uid import UID

# Deformable Spatial Registration CIOD
DEFORMABLE_SPATIAL_REGSITRATION_IOD = {
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
    "DeformableRegistrationSequence": {
        "SourceFrameOfReferenceUID": {"type": UID},
        "submap": {
            "DeformableRegistrationGridSequence": {
                "ImagePositionPatient": {},
                "ImageOrientationPatient": {},
                "GridDimensions": {},
                "GridResolution": {},
                "VectorGridData": {}
            },
            "PreDeformationMatrixRegistrationSequence": {},
            "PostDeformationMatrixRegistrationSequence": {},
            "RegistrationTypeCodeSequence": {},
        }
    },

}
