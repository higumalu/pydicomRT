"""IOD validation for RT Dose."""

from pydicom.dataset import Dataset

from pydicomrt.utils.validate_dcm_info import build_check_result, check_iod

from .builder import check_dose_plan_reference
from .iod import RTDOSE_IOD


def check_rtdose_iod(dose_ds: Dataset) -> dict:
    """
    Validate an RT Dose dataset against its CIOD.

    Also covers the conditional requirement most Dose Summation Types place on Referenced
    RT Plan Sequence, which a field-presence check alone would miss.

    Parameters
    ----------
    dose_ds : Dataset
        RT Dose object.

    Returns
    -------
    dict
        ``{"result": bool, "content": list of str}``. ``result`` is True when conformant
        and ``content`` empty; otherwise one message per problem.

    See Also
    --------
    RTDoseBuilder : Produces conformant objects.
    check_dose_plan_reference : The Type 1C check on its own.

    Notes
    -----
    The plan-reference rule is the one most commonly tripped: Referenced RT Plan Sequence
    is required for every Dose Summation Type except ``"FRACTION"``, and the default
    Summation Type is ``"PLAN"``. Structural conformance only -- dose values and geometry
    are not checked.
    """
    errors = check_iod(ds=dose_ds, config_map=RTDOSE_IOD)
    errors += check_dose_plan_reference(dose_ds)
    return build_check_result(errors)
