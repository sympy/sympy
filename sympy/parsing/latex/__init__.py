from sympy.external import import_module

from .errors import LaTeXSyntaxError  # noqa


def parse_latex(s):
    """Converts the string ``s`` to a SymPy ``Expr``

    Parameters
    ==========

    s : str
    The LaTeX string to parse.

    >>> from sympy.parsing.latex import parse_latex
    >>> parse_latex("1 + 1")  # doctest: +SKIP
    2
    """

    _latex = import_module(
        'sympy.parsing.latex._parse_latex_antlr',
        __import__kwargs={'fromlist': ['X']})

    if _latex is not None:
        return _latex.parse_latex(s)
