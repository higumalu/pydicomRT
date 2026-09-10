"""
Check RTSTRUCT dataset Class Information Object Definition (IOD)

Author: Higumalu
Date: 2025-06-13
"""
from pydicom.dataset import Dataset

from pydicomrt.utils.validate_dcm_info import build_check_result, check_iod, ValidationError
from .iod import RTSTRUCT_IOD


def contour_data_validator(value):
    if hasattr(value, '__len__'):
        if len(value) == 0:
            raise ValidationError("Contour data is empty")
        elif len(value) < 9:
            raise ValidationError("Contour data less then 3 points")
        elif len(value) % 3 != 0:
            raise ValidationError("Contour data length is not multiple of 3")
    else:
        raise ValidationError("Contour data must be a sequence type that supports len()")


RS_VALIDATORS_MAP = {
    "ContourDataValidator": contour_data_validator,
}


def check_rtstruct_iod(rs_ds: Dataset) -> dict:
    """
    Validate an RT Structure Set against its CIOD.

    Parameters
    ----------
    rs_ds : Dataset
        RT Structure Set, as built or as read from disk.

    Returns
    -------
    dict
        ``{"result": bool, "content": list of str}``. ``result`` is True when the object
        is conformant and ``content`` is empty; otherwise ``content`` holds one message
        per missing or empty required element, each naming its path in the dataset.

    See Also
    --------
    is_rtstruct_matching_series : Whether the structure set belongs to a given series.

    Notes
    -----
    Structural conformance only: this reports elements that are absent or empty, not
    contours in the wrong place. An RTSTRUCT can pass here and still be geometrically
    wrong.

    Examples
    --------
    >>> result = check_rtstruct_iod(rs_ds)                          # doctest: +SKIP
    >>> if not result["result"]:                                    # doctest: +SKIP
    ...     print("\\n".join(result["content"]))
    """
    return build_check_result(
        check_iod(ds=rs_ds, config_map=RTSTRUCT_IOD, validators=RS_VALIDATORS_MAP)
    )

# --------------------------------------------------------------------------------------------------------------------- #

def is_rtstruct_matching_series(
    rs_ds: Dataset,
    image_series: list) -> bool:
    """
    Check whether an RT Structure Set belongs to a given image series.

    Worth calling before :func:`rtstruct_to_masks`, since rasterising a structure set
    against the wrong series produces masks in the wrong place rather than an error.

    Parameters
    ----------
    rs_ds : Dataset
        RT Structure Set.
    image_series : list of Dataset
        Candidate image series.

    Returns
    -------
    bool
        True if the two are linked and no contour references an image outside the series.

    Notes
    -----
    The test is deliberately lenient in one direction and strict in the other. Evidence
    *for* a match is either a shared FrameOfReferenceUID or a matching SeriesInstanceUID
    in ReferencedFrameOfReferenceSequence -- one is enough, since real exports frequently
    carry only one. Evidence *against* is decisive: a contour referencing a SOP Instance
    UID not present in the series returns False regardless.

    A structure set covering a *subset* of the series still matches; only references
    pointing outside it fail. An empty ``image_series`` is False.
    """
    if not image_series:
        return False

    status = 0
    first_slice = image_series[0]

    rs_for_uid = getattr(rs_ds, "FrameOfReferenceUID", None)
    series_for_uid = getattr(first_slice, "FrameOfReferenceUID", None)
    if rs_for_uid is not None and series_for_uid is not None and rs_for_uid == series_for_uid:
        status += 1

    referenced_for_seq = getattr(rs_ds, "ReferencedFrameOfReferenceSequence", None)
    if referenced_for_seq:
        try:
            rt_referenced_series_seq = getattr(referenced_for_seq[0], "RTReferencedSeriesSequence", None)
            if rt_referenced_series_seq:
                rs_series_uid = getattr(rt_referenced_series_seq[0], "SeriesInstanceUID", None)
                series_instance_uid = getattr(first_slice, "SeriesInstanceUID", None)
                if rs_series_uid is not None and series_instance_uid is not None and rs_series_uid == series_instance_uid:
                    status += 1
        except (AttributeError, IndexError, TypeError):
            pass

    sop_instance_uid_list = [
        uid for uid in (getattr(ds, "SOPInstanceUID", None) for ds in image_series)
        if uid is not None
    ]
    for roi_contour_sequence in getattr(rs_ds, "ROIContourSequence", []):
        for contour in getattr(roi_contour_sequence, "ContourSequence", []):
            for contour_image in getattr(contour, "ContourImageSequence", []):
                referenced_sop_uid = getattr(contour_image, "ReferencedSOPInstanceUID", None)
                if referenced_sop_uid is not None and referenced_sop_uid not in sop_instance_uid_list:
                    status = 0
                    break

    return status > 0
