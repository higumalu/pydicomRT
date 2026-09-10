"""
Cross-module conventions.

The package's usability rests on every module reading the same way, and that is not
something the per-module tests can check -- each one only sees its own module. These
assertions are what stop the five builders drifting back apart.
"""

import importlib
import inspect

import pytest

BUILDERS = [
    ("pydicomrt.rs", "RTStructBuilder"),
    ("pydicomrt.dose", "RTDoseBuilder"),
    ("pydicomrt.ct", "CTBuilder"),
    ("pydicomrt.reg", "SpatialRegistrationBuilder"),
    ("pydicomrt.reg", "DeformableSpatialRegistrationBuilder"),
]

CHECKERS = [
    ("pydicomrt.rs", "check_rtstruct_iod"),
    ("pydicomrt.reg", "check_spatial_reg_iod"),
    ("pydicomrt.reg", "check_deformable_reg_iod"),
    ("pydicomrt.dose", "check_rtdose_iod"),
    ("pydicomrt.ct", "check_ct_iod"),
]

MODALITY_PACKAGES = ["pydicomrt.rs", "pydicomrt.reg", "pydicomrt.dose", "pydicomrt.ct"]


def resolve(module_name, attribute):
    return getattr(importlib.import_module(module_name), attribute)


class TestBuilderProtocol:
    @pytest.mark.parametrize("module_name, name", BUILDERS)
    def test_has_build_and_uid_prefix(self, module_name, name):
        builder = resolve(module_name, name)

        assert callable(getattr(builder, "build", None))
        assert callable(getattr(builder, "set_uid_prefix", None))

    @pytest.mark.parametrize("module_name, name", BUILDERS)
    def test_public_methods_are_set_add_or_build(self, module_name, name):
        """
        The whole vocabulary is set_* for single values, add_* for repeatable content, and
        build() to finish. Anything else is a method a reader has to learn separately.
        """
        builder = resolve(module_name, name)
        public = [
            method for method, _ in inspect.getmembers(builder, inspect.isfunction)
            if not method.startswith("_")
        ]

        unexpected = [
            method for method in public
            if not (method.startswith(("set_", "add_")) or method == "build")
        ]
        assert unexpected == [], f"{name} exposes {unexpected}"

    @pytest.mark.parametrize("module_name, name", BUILDERS)
    def test_set_and_add_methods_are_annotated_to_chain(self, module_name, name):
        """
        Accepts either the class itself or a self-type TypeVar: an inherited setter is
        correctly annotated with the base class, and a TypeVar preserves the subclass.
        """
        builder = resolve(module_name, name)
        allowed = {cls.__name__ for cls in builder.__mro__} | {"BuilderT"}

        for method, function in inspect.getmembers(builder, inspect.isfunction):
            if not method.startswith(("set_", "add_")):
                continue
            annotation = inspect.signature(function).return_annotation
            assert annotation is not inspect.Signature.empty, f"{name}.{method}"
            rendered = str(annotation).strip("~'\"")
            assert rendered in allowed, (
                f"{name}.{method} should return the builder so calls chain, "
                f"annotated {annotation!r}"
            )


class TestCheckerProtocol:
    @pytest.mark.parametrize("module_name, name", CHECKERS)
    def test_named_check_something_iod(self, module_name, name):
        assert name.startswith("check_") and name.endswith("_iod")
        assert callable(resolve(module_name, name))

    @pytest.mark.parametrize("module_name, name", CHECKERS)
    def test_takes_one_dataset(self, module_name, name):
        parameters = list(inspect.signature(resolve(module_name, name)).parameters)

        assert len(parameters) == 1, f"{name}{tuple(parameters)}"

    def test_all_return_the_same_shape(self):
        """
        Verified against a real dataset rather than by reading the source: an empty Dataset
        fails every IOD, which is exactly the branch that has to report consistently.
        """
        from pydicom.dataset import Dataset

        for module_name, name in CHECKERS:
            result = resolve(module_name, name)(Dataset())

            assert set(result) == {"result", "content"}, name
            assert result["result"] is False, name
            assert isinstance(result["content"], list), name
            assert all(isinstance(item, str) for item in result["content"]), name


class TestPackageLayout:
    @pytest.mark.parametrize("package", MODALITY_PACKAGES)
    def test_each_modality_has_builder_check_and_iod(self, package):
        """One layout per package, so finding anything in one tells you where it is in the rest."""
        for module in ("builder", "check", "iod"):
            importlib.import_module(f"{package}.{module}")

    @pytest.mark.parametrize("package", MODALITY_PACKAGES + ["pydicomrt.utils"])
    def test_public_api_is_declared(self, package):
        module = importlib.import_module(package)

        assert getattr(module, "__all__", None), f"{package} declares no __all__"
        for name in module.__all__:
            assert hasattr(module, name), f"{package}.__all__ lists missing {name}"
