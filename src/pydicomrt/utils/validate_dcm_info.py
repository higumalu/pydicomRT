"""
Check Dicom Class Information Object Definition (IOD)

Author: Higumalu
Date: 2025-06-13
"""
from typing import Callable, Dict, List
from pydicom.tag import Tag
from pydicom.sequence import Sequence
from pydicom.datadict import tag_for_keyword

class ValidationError(Exception):
    pass


def build_check_result(errors: List[str]) -> Dict[str, object]:
    """
    Wrap a list of IOD problems in the result shape every ``check_*_iod`` returns.

    Centralised so the four modality checkers cannot drift into reporting differently.

    Returns
    -------
    dict
        ``{"result": bool, "content": list[str]}`` -- ``result`` is True when clean.
    """
    errors = list(errors)
    return {"result": not errors, "content": errors}

def resolve_tag(key):
    if isinstance(key, str) and not key.lower().startswith("0x"):
        tag = tag_for_keyword(key)
        if tag is None:
            raise ValueError(f"Unknown DICOM keyword: {key}")
        return tag
    return Tag(key)

def check_iod(
    ds,
    config_map: Dict[str, dict],
    validators: Dict[str, Callable] = None,
    path: str = ""
    ) -> List[str]:
    """
    Check a DICOM dataset against an IOD description.

    The engine behind every ``check_*_iod`` function. Call those instead unless you are
    describing a new object type; this one takes the IOD as data so a new modality needs
    an ``iod.py`` rather than new validation code.

    Parameters
    ----------
    ds : Dataset
        The dataset to check. Nested sequences are recursed into via ``submap``.
    config_map : dict
        Keyword to requirement, nestable::

            {
                "PatientID":     {},
                "ROIContourSequence": {
                    "min_items": 1,
                    "submap": {
                        "ReferencedROINumber": {"nonempty": True},
                    },
                },
                "DoseSummationType": {"nonempty": True, "validator": ["dose_summation"]},
            }

        Elements are required unless ``optional=True``. ``nonempty=True`` rejects empty
        values; ``min_items`` and ``max_items`` constrain sequence cardinality. ``type``
        is a Python value class, not the DICOM requirement number. ``validator`` names
        entries in ``validators``. ``submap`` describes sequence items.
    validators : dict of str to callable, optional
        Named checks beyond presence, as ``{name: func}``. Each is called with the element
        value and raises ``ValidationError`` on failure. Used
        for the conditional (Type 1C) requirements presence alone cannot express.
    path : str, optional
        Prefix used when reporting nested elements; set by the recursion. Leave it alone
        at the top level, where it produces messages ending "in root".

    Returns
    -------
    list of str
        One message per problem, empty when the dataset conforms. Pass it through
        :func:`build_check_result` to get the ``{"result", "content"}`` shape the public
        checkers return.

    See Also
    --------
    build_check_result : Wrap the returned list in the public result shape.

    Notes
    -----
    This checks that required elements are *present*, not that their values are correct.
    A conformant object can still describe the wrong geometry.
    """
    # print(validators)
    errors = []
    for key, cfg in config_map.items():
        tag = resolve_tag(key)
        elem = ds.get(tag, None)
        loc = path or "root"

        if elem is None:
            if cfg.get("optional", False):
                continue
            else:
                errors.append(f"Missing {key} in {loc}")
                continue

        val = elem.value
        if cfg.get("nonempty") and (val is None or (hasattr(val, "__len__") and len(val) == 0)):
            errors.append(f"Empty {key} in {loc}")
        if elem.VR == "SQ":
            if len(val) < cfg.get("min_items", 0):
                errors.append(f"{key} in {loc} requires at least {cfg['min_items']} item(s)")
            if "max_items" in cfg and len(val) > cfg["max_items"]:
                errors.append(f"{key} in {loc} permits at most {cfg['max_items']} item(s)")

        # Check type
        expected_type = cfg.get("type")
        if expected_type and not isinstance(val, expected_type):
            errors.append(f"{key} in {loc} should be {expected_type.__name__}, got {type(val).__name__} {val}")

        # Check value
        expected_value = cfg.get("value")
        if expected_value and val != expected_value:
            errors.append(f"{key} in {loc} should be {expected_value}, got {val}")

        # ------------------------------ Custom validator --------------------------------------- #
        for vname in cfg.get("validator", []):
            func = validators.get(vname) if validators else None
            if func is None:
                # raise ValueError(f"Unknown validator '{vname}', provide a dict of validators")
                errors.append(f"No validator named '{vname}' provided for {key} in {loc}, provide a dict of validators")
                continue
            try:
                func(val)
            except ValidationError as e:
                errors.append(f"{key} in {loc} failed '{vname}': {e}")
        # ---------------------------------------------------------------------------------------- #

        # If there is a submap, recursively check the Sequence
        submap = cfg.get("submap")
        if submap:
            if elem.VR != "SQ" or not isinstance(val, Sequence):
                errors.append(f"{key} in {loc} should be SQ, but VR={elem.VR}")
            else:
                for idx, item in enumerate(val):
                    errors += check_iod(item, submap, path=f"{loc}.{key}[{idx}]", validators=validators)

    return errors


if __name__ == "__main__":
    # Example validator function
    def example_string_len_16_validator(value):
        if len(value) > 16:
            raise ValidationError(f"String length should be <= 16, but got {len(value)}")

    EXAMPLE_VAILDATORS = {
        "string_len_16": example_string_len_16_validator,
    }

    EXAMPLE_CONFIG = {
        "PatientName": {
            "type": str,
            "validator": ["string_len_16"],
        },
    }
