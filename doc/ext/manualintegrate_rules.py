"""
Table of the rules of :mod:`sympy.integrals.manualintegrate`.

The ``manualintegrate-rules`` directive collects the concrete subclasses of
``Rule`` and builds a table from their docstrings. The ``Explanation`` section
of a rule docstring holds a ``.. math::`` block with the integral before and
after the rule, separated by ``\\longrightarrow``, followed by the conditions
under which the rule applies; the ``.. math::`` blocks of the ``Examples``
section are the examples. The full docstrings are rendered after the table.
"""
from __future__ import annotations

import importlib
import inspect
from typing import Any, TYPE_CHECKING

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import StringList

if TYPE_CHECKING:
    from sphinx.application import Sphinx

MODULE = "sympy.integrals.manualintegrate"
ARROW = r"\longrightarrow"


def _rule_classes():
    module = importlib.import_module(MODULE)
    classes = [obj for obj in vars(module).values()
               if inspect.isclass(obj) and issubclass(obj, module.Rule)
               and obj.__module__ == MODULE and not inspect.isabstract(obj)]
    classes.sort(key=lambda cls: inspect.getsourcelines(cls)[1])
    return classes


def _sections(docstring):
    """Split a docstring into its sections: {name: lines}, summary under ''."""
    lines = inspect.cleandoc(docstring or "").splitlines()
    sections: dict[str, list[str]] = {"": []}
    name = ""
    i = 0
    while i < len(lines):
        line = lines[i]
        underline = lines[i + 1].strip() if i + 1 < len(lines) else ""
        if (line.strip() and underline and set(underline) <= {"=", "-"}
                and len(underline) >= len(line.strip())):
            name = line.strip()
            sections[name] = []
            i += 2
            continue
        sections[name].append(line)
        i += 1
    return sections


def _math_blocks(lines):
    """Return the ``.. math::`` blocks (each joined into one line) and the
    remaining lines."""
    blocks, other = [], []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith(".. math::"):
            indent = len(line) - len(line.lstrip())
            content = [line.split(".. math::", 1)[1].strip()]
            i += 1
            while i < len(lines) and (not lines[i].strip() or
                                      len(lines[i]) - len(lines[i].lstrip()) > indent):
                content.append(lines[i].strip())
                i += 1
            blocks.append(" ".join(c for c in content if c))
        else:
            other.append(line)
            i += 1
    return blocks, other


def _row(cls):
    sections = _sections(cls.__doc__)
    blocks, prose = _math_blocks(sections.get("Explanation", []))
    before = after = ""
    if blocks and ARROW in blocks[0]:
        before, after = (s.strip() for s in blocks[0].split(ARROW, 1))
    conditions = " ".join(l.strip() for l in prose if l.strip())
    examples, _ = _math_blocks(sections.get("Examples", []))
    summary = " ".join(l.strip() for l in sections[""] if l.strip())
    return summary, before, after, conditions, examples


class ManualintegrateRules(Directive):

    has_content = False

    def run(self):
        lines = [
            ".. list-table::",
            "   :header-rows: 1",
            "   :widths: 14 28 30 28",
            "",
            "   * - Rule",
            "     - Integral",
            "     - Result",
            "     - Examples",
        ]
        classes = _rule_classes()
        for cls in classes:
            summary, before, after, conditions, examples = _row(cls)
            lines.append("   * - :class:`~%s.%s`" % (MODULE, cls.__name__))
            if before:
                lines.append("     - .. math:: %s" % before)
                lines.append("     - .. math:: %s" % after)
                if conditions:
                    lines.append("")
                    lines.append("       %s" % conditions)
            else:
                lines.append("     - %s" % summary)
                lines.append("     -")
            if examples:
                for k, example in enumerate(examples):
                    lines.append("     - .. math:: %s" % example if k == 0
                                 else "       .. math:: %s" % example)
                    lines.append("")
            else:
                lines.append("     -")
        lines.append("")
        lines.append(".. autoclass:: %s.Rule" % MODULE)
        for cls in classes:
            lines.append(".. autoclass:: %s.%s" % (MODULE, cls.__name__))
        node = nodes.section()
        node.document = self.state.document
        self.state.nested_parse(StringList(lines, source=__file__),
                                self.content_offset, node)
        return node.children


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_directive("manualintegrate-rules", ManualintegrateRules)
    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
