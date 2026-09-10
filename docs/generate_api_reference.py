#!/usr/bin/env python
"""
Generate docs/api-reference.md from the package's docstrings.

The reference is generated rather than written so it cannot drift from the code. Edit the
docstrings, then re-run this; editing api-reference.md directly loses the change.

    python docs/generate_api_reference.py

``test_docstring_conventions.py`` enforces the shape this relies on, so a docstring that
would render badly fails the suite before it reaches here.
"""

import importlib
import inspect
import re
import sys
from pathlib import Path

PACKAGES = [
    ("pydicomrt.utils", "utils", "Series loading, geometry, and SimpleITK conversion."),
    ("pydicomrt.rs", "rs", "RT Structure Sets: masks to contours and back."),
    ("pydicomrt.reg", "reg", "DICOM registration objects: build, read, validate."),
    ("pydicomrt.reg.method", "reg.method", "Registration algorithms."),
    ("pydicomrt.reg.pipeline", "reg.pipeline", "Staged registration and transform composition."),
    ("pydicomrt.dose", "dose", "RT Dose: build and read dose grids."),
    ("pydicomrt.ct", "ct", "CT Image: emit a series from a volume."),
]

SECTION = re.compile(r"^([A-Z][a-zA-Z ]+)\n(-{3,})$", re.MULTILINE)


def anchor(label):
    return label.lower().replace(".", "").replace("_", "-").replace(" ", "-")


# Sections whose body is a list of "name : type" entries rather than prose.
PARAMETER_SECTIONS = {"Parameters", "Returns", "Yields", "Raises", "Warns",
                      "Other Parameters", "Attributes", "Methods"}

# RST roles pydicom-style docstrings use, which markdown would show verbatim.
RST_ROLE = re.compile(r":(?:func|meth|class|mod|data|attr|obj):`~?([^`]+)`")


def demote_rst(text):
    """Turn the RST inline markup numpydoc allows into its markdown equivalent."""
    text = RST_ROLE.sub(r"`\1`", text)
    text = text.replace("``", "`")
    return text


def summary(obj):
    """First line of the docstring -- what the index entry shows."""
    doc = inspect.getdoc(obj) or ""
    return demote_rst(doc.strip().split("\n")[0]) if doc else ""


def signature_of(obj):
    try:
        return f"{obj.__name__}{inspect.signature(obj)}"
    except (ValueError, TypeError):
        return obj.__name__


def split_sections(doc):
    """Split a numpydoc docstring into (title, body) pairs; the summary has title None."""
    lines = doc.split("\n")
    sections, title, body = [], None, []
    index = 0
    while index < len(lines):
        following = lines[index + 1] if index + 1 < len(lines) else ""
        if (re.fullmatch(r"[A-Z][a-zA-Z ]+", lines[index].strip())
                and re.fullmatch(r"-{3,}", following.strip())):
            sections.append((title, body))
            title, body = lines[index].strip(), []
            index += 2
            continue
        body.append(lines[index])
        index += 1
    sections.append((title, body))
    return sections


def render_prose(body):
    """
    Prose, with RST literal blocks turned into fenced code.

    ``foo::`` followed by an indented block is RST's literal block; markdown needs a
    fence, and without one the block renders as a run-on paragraph.
    """
    out, index = [], 0
    while index < len(body):
        line = body[index]
        if line.rstrip().endswith("::"):
            out.append(demote_rst(line.rstrip()[:-2].rstrip() + ":"))
            index += 1
            while index < len(body) and not body[index].strip():
                index += 1
            block = []
            while index < len(body) and (not body[index].strip() or body[index].startswith("    ")):
                block.append(body[index][4:] if body[index].startswith("    ") else body[index])
                index += 1
            while block and not block[-1].strip():
                block.pop()
            out += ["", "```python", *block, "```", ""]
            continue
        out.append(demote_rst(line))
        index += 1
    return "\n".join(out).strip()


def render_parameters(body):
    """
    Render "name : type" entries as a markdown definition list.

    Markdown collapses the indentation numpydoc relies on, so an unconverted block runs
    the parameter names and their descriptions together into one paragraph.
    """
    out, pending = [], []
    index = 0
    while index < len(body):
        line = body[index]
        is_entry = bool(line.strip()) and not line.startswith(" ")

        if is_entry:
            if pending:
                out += pending + [""]
            name, _, kind = line.partition(" : ")
            head = f"- **`{demote_rst(name.strip())}`**"
            if kind.strip():
                head += f" — *{demote_rst(kind.strip())}*"
            pending = [head]
            index += 1
            continue

        text = line.strip()
        if text.endswith("::"):
            # An RST literal block inside a description. Fence it, indented to stay
            # inside the bullet -- unfenced, markdown reflows it into the prose.
            pending += [f"  {demote_rst(text[:-2].rstrip())}:", ""]
            index += 1
            while index < len(body) and not body[index].strip():
                index += 1
            base = len(body[index]) - len(body[index].lstrip()) if index < len(body) else 0
            block = []
            while index < len(body) and (not body[index].strip()
                                         or len(body[index]) - len(body[index].lstrip()) >= base):
                block.append(body[index][base:] if body[index].strip() else "")
                index += 1
            while block and not block[-1].strip():
                block.pop()
            pending += ["  ```python", *[f"  {b}" for b in block], "  ```", ""]
            continue

        if text:
            pending.append(f"  {demote_rst(text)}")
        elif pending:
            pending.append("")
        index += 1

    if pending:
        out += pending
    return "\n".join(out).strip()


SKIP_DIRECTIVE = re.compile(r"\s*#\s*doctest:\s*\+SKIP\s*$")


def render_examples(body):
    """
    Wrap each run of >>> lines in a python fence, leaving prose between them.

    Expected output follows its `>>>` line with no marker of its own, so it belongs
    inside the same fence -- a blank line is what ends an example. The `+SKIP` directives
    are stripped: they tell the test suite the snippet needs data it has not got, and
    mean nothing to a reader.
    """
    out, buffer = [], []

    def flush():
        nonlocal buffer
        if buffer:
            out.extend(["```python", *buffer, "```", ""])
            buffer = []

    for line in body:
        if line.startswith(">>>") or line.startswith("..."):
            buffer.append(SKIP_DIRECTIVE.sub("", line))
        elif buffer and line.strip():
            buffer.append(line)          # expected output belongs to the example above
        else:
            flush()
            if line.strip():
                out.append(demote_rst(line))
            elif out and out[-1]:
                out.append("")
    flush()
    return "\n".join(out).strip()


def render_see_also(body, index):
    """
    Render See Also as a linked list.

    numpydoc's format is ``name, name : description``. Names present in this reference
    become links, which is most of the value of the section -- an unlinked name is just
    a string the reader has to search for.
    """
    out = []
    for line in body:
        text = line.strip()
        if not text:
            continue
        if line.startswith("    ") and out:
            out[-1] += " " + demote_rst(text)      # wrapped description
            continue
        names, _, description = text.partition(" : ")
        linked = []
        for name in names.split(","):
            name = name.strip()
            if not name:
                continue
            target = index.get(name) or index.get(name.split(".")[-1])
            linked.append(f"[`{name}`](#{target})" if target else f"`{name}`")
        entry = "- " + ", ".join(linked)
        if description:
            entry += f" — {demote_rst(description.strip())}"
        out.append(entry)
    return "\n".join(out)


def render_docstring(doc, index=None):
    """Render one numpydoc docstring as markdown."""
    if not doc:
        return "_Undocumented._"

    parts = []
    for title, body in split_sections(doc):
        if title is None:
            rendered = render_prose(body)
            if rendered:
                parts.append(rendered)
            continue
        if title == "Examples":
            rendered = render_examples(body)
        elif title == "See Also":
            rendered = render_see_also(body, index or {})
        elif title in PARAMETER_SECTIONS:
            rendered = render_parameters(body)
        else:
            rendered = render_prose(body)
        if rendered:
            parts += [f"**{title}**", rendered]
    return "\n\n".join(parts)


def public_members(module):
    for name in getattr(module, "__all__", []):
        obj = getattr(module, name, None)
        if inspect.isfunction(obj) or inspect.isclass(obj):
            yield name, obj


def main():
    parts = [
        "# pydicomRT API Reference",
        "",
        "Every public function and class, with parameters, return values and the behaviour",
        "that is not visible from the signature.",
        "",
        "This page is **generated from the docstrings** by `docs/generate_api_reference.py`.",
        "Edit the docstrings and re-run it; edits made here are lost on the next run. The",
        "same text is what `help(...)` prints and what an IDE shows on hover.",
        "",
        "For task-oriented instructions see [user-guide.md](user-guide.md); for how the",
        "library is put together see [architecture.md](architecture.md).",
        "",
        "---",
        "",
        "## Contents",
        "",
    ]

    modules = [(importlib.import_module(path), path, label, blurb)
               for path, label, blurb in PACKAGES]

    # Symbol -> anchor, so See Also entries can link. Methods are indexed under both
    # "Builder.method" and "method", since docstrings use whichever reads better.
    index = {}
    for module, _, label, _ in modules:
        for name, obj in public_members(module):
            index.setdefault(name, anchor(label + name))
            if inspect.isclass(obj):
                for method_name, _ in inspect.getmembers(obj, inspect.isfunction):
                    if not method_name.startswith("_"):
                        index.setdefault(f"{name}.{method_name}", anchor(label + name))

    for module, _, label, blurb in modules:
        parts.append(f"**[`{label}`](#{anchor(label)})** — {blurb}")
        parts.append("")
        for name, obj in public_members(module):
            kind = "class" if inspect.isclass(obj) else "func"
            parts.append(f"- [`{name}`](#{anchor(label + name)}) <sub>{kind}</sub> — {summary(obj)}")
        parts.append("")

    for module, _, label, blurb in modules:
        parts += ["---", "", f'<a id="{anchor(label)}"></a>', "",
                  f"## {label}", "", blurb, ""]

        for name, obj in public_members(module):
            parts.append(f'<a id="{anchor(label + name)}"></a>')
            parts.append("")
            parts.append(f"### `{label}.{name}`")
            parts.append("")
            if inspect.isclass(obj):
                parts += ["```python", signature_of(obj), "```", ""]
            else:
                parts += ["```python", signature_of(obj), "```", ""]
            parts.append(render_docstring(inspect.getdoc(obj), index))
            parts.append("")

            if inspect.isclass(obj):
                methods = [(n, m) for n, m in inspect.getmembers(obj, inspect.isfunction)
                           if not n.startswith("_")]
                for method_name, method in methods:
                    parts.append(f"#### `{name}.{method_name}`")
                    parts.append("")
                    parts += ["```python", signature_of(method), "```", ""]
                    parts.append(render_docstring(inspect.getdoc(method), index))
                    parts.append("")

    output = Path(__file__).parent / "api-reference.md"
    output.write_text("\n".join(parts).rstrip() + "\n")

    entries = sum(1 for module, _, _, _ in modules for _ in public_members(module))
    print(f"wrote {output} — {entries} top-level entries, {len(output.read_text())} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
