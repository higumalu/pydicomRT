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
    """
    Declare an ROI in all three sequences an RTSTRUCT needs it in.

    An ROI is not one element but three parallel entries, cross-referenced by ROINumber:
    StructureSetROISequence (name and frame of reference), ROIContourSequence (colour and,
    later, the contours), and RTROIObservationsSequence (interpreted type). Adding it to
    only some of them produces a file most TPSs silently ignore.

    This creates all three, with an empty ContourSequence for the contour extraction to
    fill. Called by :meth:`RTStructBuilder.add_roi`, which is the supported entry point.

    Parameters
    ----------
    rs_ds : Dataset
        RT Structure Set with the three sequences present, from
        :func:`create_rtstruct_dataset`. Modified in place.
    roi_color : list of int
        ROIDisplayColor as RGB, 0-255.
    roi_number : int
        ROINumber. Must be unique within the dataset; nothing here enforces that.
    roi_name : str
        ROIName, as the TPS displays it.
    roi_description : str
        ROIDescription. May be empty.
    roi_interpreted_type : str, optional
        RTROIInterpretedType, e.g. ``"ORGAN"``, ``"PTV"``, ``"CTV"``, ``"EXTERNAL"``.
        Default ``"ORGAN"``.

    Returns
    -------
    Dataset
        ``rs_ds``, modified in place.

    See Also
    --------
    RTStructBuilder.add_roi : Adds the ROI *and* its contours, and rejects duplicate
        numbers.
    add_contour_sequence_from_mask3d : Fills in the contours afterwards.
    """
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
