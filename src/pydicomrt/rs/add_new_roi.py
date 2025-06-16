from pydicom.dataset import Dataset
from pydicom.sequence import Sequence


def create_roi_into_rs_ds(
    rs_ds: Dataset,
    roi_color: list,
    roi_number: int,
    roi_name: str,
    roi_description: str,
    roi_interpreted_type: str = "ORGAN",
    ) -> Dataset:
    frame_of_ref_uid = rs_ds.FrameOfReferenceUID
    rs_ds.StructureSetROISequence.append(create_empty_structure_set_roi(roi_number, roi_name, roi_description, frame_of_ref_uid))
    rs_ds.ROIContourSequence.append(create_empty_roi_contour_sequence(roi_number, roi_color))
    rs_ds.RTROIObservationsSequence.append(create_empty_rtroi_observation(roi_number, roi_interpreted_type))
    return rs_ds

def create_empty_structure_set_roi(roi_number,
                                   roi_name,
                                   roi_description,
                                   frame_of_ref_uid,
                                   roi_generation_algorithm="AUTOMATIC",
                                   roi_generation_description="KumaGenerated",
                                   ) -> Dataset:
    structure_set_roi = Dataset()
    structure_set_roi.ROINumber = roi_number
    structure_set_roi.ROIName = roi_name
    structure_set_roi.ROIDescription = roi_description
    structure_set_roi.ROIGenerationAlgorithm = roi_generation_algorithm
    structure_set_roi.ROIGenerationDescription = roi_generation_description
    structure_set_roi.ReferencedFrameOfReferenceUID = frame_of_ref_uid
    return structure_set_roi

def create_empty_roi_contour_sequence(roi_number, roi_color) -> Dataset:
    roi_contour_sequence = Dataset()
    roi_contour_sequence.ReferencedROINumber = roi_number    # also in create_structure_set_roi
    roi_contour_sequence.ROIDisplayColor = roi_color
    roi_contour_sequence.ContourSequence = Sequence()
    return roi_contour_sequence

def create_empty_rtroi_observation(roi_number, roi_interpreted_type='ORGAN') -> Dataset:
    rtroi_observation = Dataset()
    rtroi_observation.ObservationNumber = roi_number
    rtroi_observation.ReferencedROINumber = roi_number
    rtroi_observation.ROIObservationDescription = "Type:Soft,Range:*/*,Fill:0,Opacity:0.0,Thickness:1,LineThickness:2,read-only:false"
    rtroi_observation.private_creators = "higumalu"
    rtroi_observation.RTROIInterpretedType = roi_interpreted_type
    rtroi_observation.ROIInterpreter = ""
    return rtroi_observation
