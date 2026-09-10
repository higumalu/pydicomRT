"""
Spatial Registration CIOD field map.

Follows PS3.3 A.39.1 (Spatial Registration IOD) and C.20.2 (Spatial Registration Module).

Nesting matters here. Registration Type Code Sequence (0070,030D) belongs *inside* Matrix
Registration Sequence, alongside Matrix Sequence -- a reader looking for it one level up in
the Registration Sequence item will not find it. Sequence contents must sit under a
``submap`` key or ``check_iod`` never descends into them.
"""

from pydicom.uid import UID

# Spatial Registration CIOD
SPATIAL_REGISTRATION_IOD = {
    "PatientID": {},
    "PatientName": {},
    "PatientBirthDate": {},
    "PatientSex": {},
    "StudyInstanceUID": {"nonempty": True},
    "StudyID": {},
    "StudyDate": {},
    "StudyTime": {},
    "AccessionNumber": {},
    "SeriesInstanceUID": {"type": UID, "nonempty": True},
    "SeriesNumber": {},
    "Modality": {"value": "REG"},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID, "nonempty": True},
    "PositionReferenceIndicator": {},
    "SOPClassUID": {"type": UID, "nonempty": True},
    "SOPInstanceUID": {"type": UID, "nonempty": True},
    "InstanceNumber": {"nonempty": True},
    "ContentDate": {"nonempty": True},
    "ContentTime": {"nonempty": True},
    "ContentLabel": {"nonempty": True},
    "ContentDescription": {},
    "ContentCreatorName": {},
    "RegistrationSequence": {
        "min_items": 1,
        "submap": {
            "FrameOfReferenceUID": {"type": UID, "nonempty": True, "optional": True},
            "ReferencedImageSequence": {"optional": True, "min_items": 1, "submap": {
                "ReferencedSOPClassUID": {"nonempty": True},
                "ReferencedSOPInstanceUID": {"nonempty": True},
            }},
            "MatrixRegistrationSequence": {
                "min_items": 1, "max_items": 1,
                "submap": {
                    "MatrixSequence": {
                        "min_items": 1,
                        "submap": {
                            "FrameOfReferenceTransformationMatrixType": {"nonempty": True},
                            "FrameOfReferenceTransformationMatrix": {"nonempty": True},
                        }
                    },
                    # Type 2: required, but zero items is allowed.
                    "RegistrationTypeCodeSequence": {
                        "max_items": 1,
                        "submap": {
                            "CodeValue": {},
                            "CodingSchemeDesignator": {},
                            "CodeMeaning": {},
                        },
                    },
                }
            },
        }
    },
    # Common Instance Reference module: User Optional in the IOD, but it is how a viewer
    # finds the images a registration belongs to.
    "ReferencedSeriesSequence": {
        "optional": True,
        "submap": {
            "SeriesInstanceUID": {"type": UID, "nonempty": True},
            "ReferencedInstanceSequence": {
                "submap": {
                    "ReferencedSOPClassUID": {"type": UID},
                    "ReferencedSOPInstanceUID": {"type": UID},
                }
            },
        },
    },
}

_FRAME_OF_REFERENCE_TRANSFORMATION_MATRIX = {
    "optional": True, "min_items": 1, "max_items": 1,
    "submap": {
        "FrameOfReferenceTransformationMatrixType": {"nonempty": True},
        "FrameOfReferenceTransformationMatrix": {"nonempty": True},
    }
}

# Deformable Spatial Registration CIOD
DEFORMABLE_SPATIAL_REGISTRATION_IOD = {
    "PatientID": {},
    "PatientName": {},
    "PatientBirthDate": {},
    "PatientSex": {},
    "StudyInstanceUID": {"nonempty": True},
    "StudyID": {},
    "StudyDate": {},
    "StudyTime": {},
    "AccessionNumber": {},
    "SeriesInstanceUID": {"type": UID, "nonempty": True},
    "SeriesNumber": {},
    "Modality": {"value": "REG"},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID, "nonempty": True},
    "PositionReferenceIndicator": {},
    "SOPClassUID": {"type": UID, "nonempty": True},
    "SOPInstanceUID": {"type": UID, "nonempty": True},
    "InstanceNumber": {"nonempty": True},
    "ContentDate": {"nonempty": True},
    "ContentTime": {"nonempty": True},
    "ContentLabel": {"nonempty": True},
    "ContentDescription": {},
    "ContentCreatorName": {},
    "DeformableRegistrationSequence": {
        "min_items": 1,
        "submap": {
            "SourceFrameOfReferenceUID": {"type": UID, "nonempty": True},
            "ReferencedImageSequence": {"optional": True, "min_items": 1, "submap": {
                "ReferencedSOPClassUID": {"nonempty": True},
                "ReferencedSOPInstanceUID": {"nonempty": True},
            }},
            "DeformableRegistrationGridSequence": {
                "optional": True, "max_items": 1,
                "submap": {
                    "ImagePositionPatient": {"nonempty": True},
                    "ImageOrientationPatient": {"nonempty": True},
                    "GridDimensions": {"nonempty": True},
                    "GridResolution": {"nonempty": True},
                    "VectorGridData": {"nonempty": True},
                }
            },
            "PreDeformationMatrixRegistrationSequence":
                _FRAME_OF_REFERENCE_TRANSFORMATION_MATRIX,
            "PostDeformationMatrixRegistrationSequence":
                _FRAME_OF_REFERENCE_TRANSFORMATION_MATRIX,
            "RegistrationTypeCodeSequence": {
                "max_items": 1,
                "submap": {
                    "CodeValue": {},
                    "CodingSchemeDesignator": {},
                    "CodeMeaning": {},
                },
            },
        }
    },
    "ReferencedSeriesSequence": {
        "optional": True,
        "submap": {
            "SeriesInstanceUID": {"type": UID, "nonempty": True},
            "ReferencedInstanceSequence": {
                "submap": {
                    "ReferencedSOPClassUID": {"type": UID},
                    "ReferencedSOPInstanceUID": {"type": UID},
                }
            },
        },
    },
}
