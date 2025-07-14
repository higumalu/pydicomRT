from pydicom.multival import MultiValue
from pydicom.uid import UID

RT_STRUCTURE_SET_IOD = {
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
    "OperatorsName": {},
    "Modality": {},
    "Manufacturer": {},
    "FrameOfReferenceUID": {"type": UID},
    "SOPClassUID": {"type": UID},
    "SOPInstanceUID": {"type": UID},
    "StructureSetLabel": {"type": str},
    "StructureSetDate": {},
    "StructureSetTime": {},
    "StructureSetROISequence": {
        "submap": {
            "ROINumber": {"type": int},
            "ROIName": {"type": str},
            "ROIGenerationAlgorithm": {},
            "ReferencedFrameOfReferenceUID": {},
        }
    },
    "ROIContourSequence": {
        "submap": {
            "ReferencedROINumber": {"type": int},
            "ROIDisplayColor": {},
            "ContourSequence": {
                "optional": True,
                "submap": {
                    "ContourImageSequence": {
                        "submap": {
                            "ReferencedSOPClassUID": {"type": UID},
                            "ReferencedSOPInstanceUID": {"type": UID},
                        }
                    },
                    "ContourGeometricType": {},
                    "NumberOfContourPoints": {"type": int},
                    "ContourData": {
                        "type": MultiValue,
                        "validator": ["ContourDataValidator"]
                        }
                }
            }
        }
    },
    "RTROIObservationsSequence": {
        "submap": {
            "ReferencedROINumber": {"type": int},
            "ObservationNumber": {"type": int},
            "RTROIInterpretedType": {},
        }
    }
}
