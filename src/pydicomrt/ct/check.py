"""IOD validation for CT Image."""

from pydicom.dataset import Dataset

from pydicomrt.utils.validate_dcm_info import build_check_result, check_iod

from .iod import CT_IMAGE_IOD


def check_ct_iod(ct_ds: Dataset) -> dict:
    """
    Validate a CT Image dataset against its CIOD.

    Parameters
    ----------
    ct_ds : Dataset
        A single slice. CT is a single-frame IOD, so each slice is checked separately.

    Returns
    -------
    dict
        ``{"result": bool, "content": list of str}``. ``result`` is True when conformant
        and ``content`` empty; otherwise one message per missing or empty required
        element.

    See Also
    --------
    CTBuilder : Produces conformant slices.

    Notes
    -----
    Structural conformance only. Check every slice, not just the first: pixel-data
    elements are per-slice, so a series can pass on slice 0 and fail later.

    Examples
    --------
    >>> slices = CTBuilder(reference_series).set_volume(image).build()   # doctest: +SKIP
    >>> bad = [i for i, s in enumerate(slices)                           # doctest: +SKIP
    ...        if not check_ct_iod(s)["result"]]
    """
    return build_check_result(check_iod(ds=ct_ds, config_map=CT_IMAGE_IOD))
