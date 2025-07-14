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
    Check if a DICOM dataset conforms to a given IOD configuration.
    :param ds: The DICOM dataset to check.
    :param config_map: The IOD configuration to check against.
        {
        "tag_name": {
            "type": type,
            "validator": ["validator_name_1", "validator_name_2"],
            "submap": {
                "tag_name": {
                    "type": type,
                    "validator": ["validator_name_1", "validator_name_2"],
                    "submap": {
                        ...
                    }
                }
            }
        }
    :param validators: A dictionary of custom validators.
        {
        "validator_name_1": validator_function_1,
        "validator_name_2": validator_function_2,
        ...
        }
    :param path: The path to the dataset in the IOD configuration.
    :return: A list of errors.
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
