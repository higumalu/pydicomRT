from pydicom.dataset import Dataset


def get_contours(
    rs_ds: Dataset,
    ) -> dict:
    """
    Read every contour out of an RT Structure Set, grouped by ROI and by image.

    Returns the points as stored -- patient coordinates in mm, no rasterisation and no
    image series required. Use this to inspect, edit or copy contours; use
    :func:`rtstruct_to_masks` when you want voxels.

    Parameters
    ----------
    rs_ds : Dataset
        RT Structure Set. Missing sequences are tolerated and yield an empty result rather
        than raising.

    Returns
    -------
    dict
        Keyed by ROI *number*::

            {roi_number: {
                "color": [R, G, B],
                "name": str,
                "dcm_contour": {
                    sop_instance_uid: {
                        "sop_class_uid": str,
                        "contours": [[x1, y1, z1, x2, y2, z2, ...], ...],
                    },
                },
            }}

        Each contour is one flat list of interleaved xyz triples, the DICOM ContourData
        layout. Reshape with ``np.reshape(contour, (-1, 3))``.

        Keys are ``pydicom.valuerep.IS``, an ``int`` subclass: ``contours[1]`` works,
        ``contours["1"]`` raises ``KeyError`` despite the repr showing ``{'1': ...}``.

    See Also
    --------
    rtstruct_to_masks : Rasterise onto an image grid instead.
    get_roi_names : Just the number-to-name mapping.
    RTStructBuilder.add_roi_from_contours : Feed this structure back into a new RTSTRUCT.

    Notes
    -----
    Contours are grouped under the SOP Instance UID of the image they reference. A contour
    whose ContourImageSequence is absent or empty is skipped, since it cannot be attributed
    to an image.

    Examples
    --------
    >>> contours = get_contours(rs_ds)                              # doctest: +SKIP
    >>> {n: d["name"] for n, d in contours.items()}                 # doctest: +SKIP
    {1: 'CTV', 2: 'SpinalCord'}
    """
    contour_dict = {}
    number_name_map = get_roi_names(rs_ds)
    for roi_contour_sequence in getattr(rs_ds, "ROIContourSequence", []):
        roi_number = getattr(roi_contour_sequence, "ReferencedROINumber", None)
        if roi_number is None:
            continue

        roi_color = getattr(roi_contour_sequence, "ROIDisplayColor", [255, 255, 255])
        roi_name = number_name_map.get(roi_number, f"ROI_{roi_number}")
        contour_dict[roi_number] = {
            'color': roi_color,
            'name': roi_name,
            'dcm_contour': {},
        }

        for contour_sequence in getattr(roi_contour_sequence, "ContourSequence", []):
            try:
                contour_image_sequence = getattr(contour_sequence, "ContourImageSequence", None)
                if not contour_image_sequence:
                    continue
                contour_image = contour_image_sequence[0]
                sop_instance_uid = getattr(contour_image, "ReferencedSOPInstanceUID", None)
                sop_class_uid = getattr(contour_image, "ReferencedSOPClassUID", None)
                contour_data = getattr(contour_sequence, "ContourData", None)
                if sop_instance_uid is None or contour_data is None:
                    continue

                if sop_instance_uid not in contour_dict[roi_number]['dcm_contour']:
                    contour_dict[roi_number]['dcm_contour'][sop_instance_uid] = {
                        'sop_class_uid': sop_class_uid,
                        'contours': [],
                    }

                contour_dict[roi_number]['dcm_contour'][sop_instance_uid]['contours'].append(contour_data)
            except (AttributeError, IndexError, TypeError):
                continue

    for roi_observations_sequence in getattr(rs_ds, "RTROIObservationsSequence", []):
        roi_number = getattr(roi_observations_sequence, "ReferencedROINumber", None)
        if roi_number is None or roi_number not in contour_dict:
            continue
        contour_dict[roi_number]['interpreted_type'] = getattr(
            roi_observations_sequence, "RTROIInterpretedType", "NOTAG"
        )

    return contour_dict


def get_roi_names(
    rs_ds: Dataset,
    ) -> dict:
    """
    Map ROI numbers to ROI names.

    The cheapest way to see what is in a structure set: it reads StructureSetROISequence
    only and never touches contour data.

    Parameters
    ----------
    rs_ds : Dataset
        RT Structure Set.

    Returns
    -------
    dict
        ``{roi_number: roi_name}``. An ROI with no ROIName is reported as
        ``"ROI_<number>"`` rather than ``None``, so the value is always usable as a label.

        Keys are ``pydicom.valuerep.IS``, an ``int`` subclass. ``names[1]`` works;
        ``names["1"]`` raises ``KeyError`` even though the repr prints ``{'1': 'CTV'}``.

    See Also
    --------
    get_contours : The same ROIs with their contour points.

    Notes
    -----
    DICOM does not require ROI names to be unique. Compare ``len(names)`` against
    ``len(set(names.values()))`` before keying anything by name --
    :func:`rtstruct_to_masks` does key by name, and collapses duplicates.

    Examples
    --------
    >>> get_roi_names(rs_ds)                                        # doctest: +SKIP
    {1: 'CTV', 2: 'SpinalCord', 3: 'BODY'}
    """
    roi_number_to_name = {}
    for ssroi in getattr(rs_ds, "StructureSetROISequence", []):
        roi_number = getattr(ssroi, "ROINumber", None)
        if roi_number is None:
            continue
        roi_name = getattr(ssroi, "ROIName", None)
        roi_number_to_name[roi_number] = roi_name if roi_name is not None else f"ROI_{roi_number}"
    return roi_number_to_name
