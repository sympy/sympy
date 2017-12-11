__all__ = ['parse_latex', 'LaTeXSyntaxError']


class LaTeXSyntaxError(Exception):
    pass


def parse_latex(latex_str):
    """ Generate a SymPy Expression from a LaTeX math string.

    Examples
    ========

    >>> from sympy.abc import a, b
    >>> from sympy.parsing.latex import parse_latex
    >>> parse_latex("a + 1")
    a + 1
    >>> parse_latex("a + \\sin b")
    a + sin(b)
    """
    from ._latex import (parse_latex as _parse_latex)
    return _parse_latex(latex_str)
