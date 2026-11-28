from __future__ import annotations
from sympy.printing.mathml import mathml
from sympy.utilities.decorator import deprecated
from sympy.utilities.mathml import _c2p
import tempfile
import subprocess


@deprecated("""\
print_gtk() is deprecated. Use preview() to display an expression in an
external viewer, or use

    print_mathml(expr, printer='presentation')

to get presentation MathML that can be displayed in a web browser.""",
    deprecated_since_version="1.15",
    active_deprecations_target="deprecated-mathml-xsl")
def print_gtk(x, start_viewer=True):
    """Print to Gtkmathview, a gtk widget capable of rendering MathML.

    .. deprecated:: 1.15

       ``print_gtk`` is deprecated. Use :func:`~sympy.printing.preview.preview`
       to display an expression in an external viewer, or use
       ``print_mathml(expr, printer='presentation')`` to get presentation
       MathML that can be displayed in a web browser. See
       :ref:`deprecated-mathml-xsl` for details.

    Needs libgtkmathview-bin"""
    with tempfile.NamedTemporaryFile('w') as file:
        file.write(_c2p(mathml(x), simple=True))
        file.flush()

        if start_viewer:
            subprocess.check_call(('mathmlviewer', file.name))
