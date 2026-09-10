"""
Docstring conventions for the public API.

The reference manual is generated from these docstrings, so a docstring that drifts from
the signature does not merely read badly -- it publishes something false. These assertions
are what keep ``docs/api-reference.md`` honest without anyone re-reading it.

The style is numpydoc, as used by numpy and scipy: underlined section headings, one
``name : type`` entry per parameter.
"""

import importlib
import inspect
import re

import pytest

PUBLIC_MODULES = [
    "pydicomrt.utils",
    "pydicomrt.rs",
    "pydicomrt.reg",
    "pydicomrt.reg.method",
    "pydicomrt.reg.pipeline",
    "pydicomrt.dose",
    "pydicomrt.ct",
]

SECTION = re.compile(r"^([A-Z][a-zA-Z ]+)\n-{3,}$", re.MULTILINE)

# numpydoc's own set. A heading outside it is a typo or an invented section, and either way
# a generator will not render it as a section.
KNOWN_SECTIONS = {
    "Parameters", "Returns", "Yields", "Raises", "Warns", "Warnings", "Other Parameters",
    "See Also", "Notes", "References", "Examples", "Attributes", "Methods",
}


def public_callables():
    """Every callable a user reaches through a package's ``__all__``, methods included."""
    for module_name in PUBLIC_MODULES:
        module = importlib.import_module(module_name)
        for name in getattr(module, "__all__", []):
            obj = getattr(module, name, None)
            if not (inspect.isfunction(obj) or inspect.isclass(obj)):
                continue
            yield module_name, name, obj
            if inspect.isclass(obj):
                for method_name, method in inspect.getmembers(obj, inspect.isfunction):
                    if method_name.startswith("_"):
                        continue
                    yield module_name, f"{name}.{method_name}", method


CALLABLES = list(public_callables())
IDS = [f"{module}:{name}" for module, name, _ in CALLABLES]

# Exceptions carry no parameters or return value of their own; prose describing when they
# are raised is the correct shape for them.
EXCEPTIONS = {name for _, name, obj in CALLABLES
              if inspect.isclass(obj) and issubclass(obj, BaseException)}


def sections(obj):
    return set(SECTION.findall(inspect.getdoc(obj) or ""))


class TestEveryPublicCallableIsDocumented:
    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_has_a_docstring(self, module_name, name, obj):
        assert inspect.getdoc(obj), f"{module_name}.{name} has no docstring"

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_opens_with_a_one_line_summary(self, module_name, name, obj):
        """
        A generated index shows the first line and nothing else, so a docstring opening
        mid-sentence -- or with a section heading -- produces an index entry that says
        nothing.
        """
        first = inspect.getdoc(obj).strip().split("\n")[0]

        assert first, f"{module_name}.{name} opens with a blank line"
        assert not first.startswith("-"), f"{module_name}.{name} opens with an underline"
        # "4x4 affine mapping ..." is a fine opening; a lowercase word is a continuation.
        assert not first[0].islower(), (
            f"{module_name}.{name} summary is not capitalised: {first!r}"
        )

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_uses_numpydoc_sections(self, module_name, name, obj):
        if name in EXCEPTIONS:
            pytest.skip("an exception class documents when it is raised, not parameters")

        found = sections(obj)

        assert found & {"Parameters", "Returns", "Methods"}, (
            f"{module_name}.{name} has no numpydoc section; found {sorted(found) or 'none'}"
        )

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_section_headings_are_real(self, module_name, name, obj):
        unknown = sections(obj) - KNOWN_SECTIONS

        assert not unknown, f"{module_name}.{name} has unknown section(s) {sorted(unknown)}"


class TestParametersMatchSignatures:
    """
    A documented parameter that no longer exists is worse than an undocumented one: it
    reads as usable and silently becomes a TypeError.
    """

    @staticmethod
    def documented_parameters(obj):
        doc = inspect.getdoc(obj) or ""
        match = re.search(r"^Parameters\n-{3,}\n(.*?)(?=\n[A-Z][a-zA-Z ]+\n-{3,}|\Z)",
                          doc, re.MULTILINE | re.DOTALL)
        if not match:
            return set()
        names = set()
        for line in match.group(1).split("\n"):
            # A parameter entry is unindented and of the form "name : type".
            entry = re.match(r"^(\*{0,2}\w[\w, ]*?)\s*:", line)
            if entry:
                names.update(part.strip().lstrip("*") for part in entry.group(1).split(","))
        return names

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_no_documented_parameter_is_missing_from_the_signature(self, module_name, name, obj):
        try:
            signature = inspect.signature(obj)
        except (ValueError, TypeError):
            pytest.skip("no introspectable signature")

        actual = {p for p in signature.parameters if p != "self"}
        documented = self.documented_parameters(obj)
        # Classes document __init__'s parameters; methods, their own.
        phantom = documented - actual

        assert not phantom, (
            f"{module_name}.{name} documents {sorted(phantom)}, "
            f"which the signature does not accept: {sorted(actual)}"
        )

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_every_parameter_is_documented(self, module_name, name, obj):
        if name in EXCEPTIONS:
            pytest.skip("an exception class documents when it is raised, not parameters")
        try:
            signature = inspect.signature(obj)
        except (ValueError, TypeError):
            pytest.skip("no introspectable signature")

        actual = {p for p in signature.parameters if p != "self"}
        if not actual:
            pytest.skip("takes no parameters")
        documented = self.documented_parameters(obj)
        if not documented and "Methods" in sections(obj):
            pytest.skip("a builder documenting its methods rather than its constructor")

        assert not actual - documented, (
            f"{module_name}.{name} does not document {sorted(actual - documented)}"
        )


class TestExamplesAreRunnable:
    """
    Every ``>>>`` line either runs or is explicitly marked as illustrative. The middle
    case -- an example that looks executable and raises NameError -- is the one that
    reaches users, because nothing checks it.
    """

    @pytest.mark.parametrize("module_name, name, obj", CALLABLES, ids=IDS)
    def test_examples_run_or_are_marked_skip(self, module_name, name, obj):
        import doctest

        examples = doctest.DocTestParser().get_examples(inspect.getdoc(obj) or "")
        # One namespace for the whole docstring, and seeded with the defining module's
        # globals -- the scope doctest itself gives these examples. So a later line may
        # use a name an earlier line bound, and only a genuinely undefined name fails.
        defining_module = inspect.getmodule(obj)
        namespace = dict(vars(defining_module)) if defining_module else {}
        for example in examples:
            if doctest.SKIP in example.options:
                continue
            try:
                exec(compile(example.source, "<docstring>", "single"), namespace)
            except Exception as error:  # noqa: BLE001 - the failure is the finding
                pytest.fail(
                    f"{module_name}.{name} example is neither runnable nor marked "
                    f"# doctest: +SKIP\n    {example.source.strip()}\n    "
                    f"{type(error).__name__}: {error}"
                )


class TestGeneratedReferenceIsCurrent:
    """
    ``docs/api-reference.md`` is generated, so it can silently fall behind the docstrings
    it was generated from. Regenerating and comparing is the only check that catches that;
    reading the file cannot.
    """

    def test_regenerating_produces_no_diff(self, tmp_path):
        import importlib.util
        import shutil
        from pathlib import Path

        docs = Path(__file__).resolve().parent.parent / "docs"
        generator = docs / "generate_api_reference.py"
        reference = docs / "api-reference.md"
        if not generator.exists():
            pytest.skip("generator not present")

        assert reference.exists(), "docs/api-reference.md is missing; run the generator"

        spec = importlib.util.spec_from_file_location("_api_reference_generator", generator)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        committed = reference.read_text()
        backup = tmp_path / "api-reference.md"
        shutil.copy(reference, backup)
        try:
            module.main()          # exec_module alone runs no main(), it is __main__-guarded
            regenerated = reference.read_text()
        finally:
            shutil.copy(backup, reference)

        assert regenerated == committed, (
            "docs/api-reference.md is stale -- a docstring changed without the reference "
            "being regenerated. Run: python docs/generate_api_reference.py"
        )
