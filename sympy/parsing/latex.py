__all__ = ['parse_latex', 'LaTeXSyntaxError']


class LaTeXSyntaxError(Exception):
    pass


def parse_latex(latex_str, parser="antlr", debug=False):
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
    if parser == "ply":
        from ._latex._ply.parser import parse_latex as _parse_latex
    elif parser == "antlr":
        from ._latex._parse_latex_antlr import parse_latex as _parse_latex

    return _parse_latex(latex_str)
