from sympy.external import import_module

from .errors import LaTeXSyntaxError  # noqa


def parse_latex(latex_str):
    _latex = import_module(
        'sympy.parsing.latex._parse_latex_antlr',
        __import__kwargs={'fromlist': ['X']})

    if _latex is not None:
        return _latex.parse_latex(latex_str)
