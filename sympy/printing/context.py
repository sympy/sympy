"""
A Printer which converts an expression into its ConTeXt equivalent.

ConTeXt is a document preparation system based on TeX, similar to LaTeX.
This printer extends the LaTeX printer since most mathematical expressions
have the same syntax in both systems.
"""
from __future__ import annotations
from typing import Any

from sympy.printing.latex import LatexPrinter, latex_escape


class ContextPrinter(LatexPrinter):
    """
    A printer for converting SymPy expressions to ConTeXt code.

    ConTeXt is a document preparation system based on TeX. Most mathematical
    expressions are identical to LaTeX, but ConTeXt uses different environment
    delimiters.

    Examples
    ========

    >>> from sympy import symbols, sqrt
    >>> from sympy.printing.context import context
    >>> x, y = symbols('x y')
    >>> context(x**2 + sqrt(2))
    'x^{2} + \\\\sqrt{2}'
    >>> context(x**2 + sqrt(2), mode='inline')
    '$x^{2} + \\\\sqrt{2}$'
    >>> context(x**2 + sqrt(2), mode='equation')
    '\\\\startformula x^{2} + \\\\sqrt{2} \\\\stopformula'
    """

    # printmethod inherited from LatexPrinter ("_latex")
    # We use all the LaTeX printer's methods for printing expressions

    _default_settings: dict[str, Any] = {
        **LatexPrinter._default_settings,
        "mode": "plain",
    }

    def __init__(self, settings=None):
        # Call parent constructor
        super().__init__(settings)

        # Validate mode settings for ConTeXt
        if 'mode' in self._settings:
            valid_modes = ['inline', 'plain', 'equation', 'equation*']
            if self._settings['mode'] not in valid_modes:
                raise ValueError("'mode' must be one of 'inline', 'plain', "
                                 "'equation' or 'equation*'")

    def doprint(self, expr) -> str:
        """
        Print an expression with ConTeXt delimiters.

        ConTeXt uses different environment delimiters than LaTeX:
        - LaTeX: \\begin{equation}...\\end{equation}
        - ConTeXt: \\startformula...\\stopformula

        Parameters
        ==========
        expr : Expr
            The expression to print

        Returns
        =======
        str
            The ConTeXt representation of the expression
        """
        # Get the core tex representation from the parent class
        from sympy.printing.printer import Printer
        tex = Printer.doprint(self, expr)

        mode = self._settings['mode']

        if mode == 'plain':
            return tex
        elif mode == 'inline':
            # Inline math is the same in both LaTeX and ConTeXt
            return r"$%s$" % tex
        elif mode == 'equation':
            # ConTeXt uses \startformula...\stopformula
            return r"\startformula %s \stopformula" % tex
        elif mode == 'equation*':
            # ConTeXt uses \startformula...\stopformula (unnumbered by default)
            # For numbered equations, use \placeformula\startformula
            return r"\startformula %s \stopformula" % tex
        else:
            # Fallback to plain mode
            return tex


def context(expr, **settings):
    r"""Convert the given expression to ConTeXt string representation.

    ConTeXt is a document preparation system based on TeX, similar to LaTeX.
    Most mathematical expressions have identical syntax in LaTeX and ConTeXt,
    but ConTeXt uses different environment delimiters.

    Parameters
    ==========

    expr : Expr
        A SymPy expression to be converted to ConTeXt.

    mode : string, optional
        Specifies how the generated code will be delimited. ``mode`` can be one
        of ``'plain'``, ``'inline'``, ``'equation'`` or ``'equation*'``. If
        ``mode`` is set to ``'plain'``, then the resulting code will not be
        delimited at all (this is the default). If ``mode`` is set to
        ``'inline'`` then inline ConTeXt ``$...$`` will be used. If ``mode`` is
        set to ``'equation'`` or ``'equation*'``, the resulting code will be
        enclosed in the ConTeXt formula environment using
        ``\startformula...\stopformula``.

    **settings : dict
        Additional settings are passed to the LatexPrinter. See the
        documentation for the ``latex()`` function for more details on
        available settings.

    Examples
    ========

    >>> from sympy import context, symbols, sqrt, Integral, pi
    >>> from sympy.abc import x, y, tau, mu
    >>> from sympy import Rational

    Basic usage:

    >>> print(context((2*tau)**Rational(7,2)))
    8 \sqrt{2} \tau^{\frac{7}{2}}

    ``mode`` option:

    >>> print(context((2*tau)**Rational(7,2), mode='plain'))
    8 \sqrt{2} \tau^{\frac{7}{2}}
    >>> print(context((2*tau)**Rational(7,2), mode='inline'))
    $8 \sqrt{2} \tau^{7 / 2}$
    >>> print(context((2*mu)**Rational(7,2), mode='equation'))
    \startformula 8 \sqrt{2} \mu^{\frac{7}{2}} \stopformula
    >>> print(context((2*mu)**Rational(7,2), mode='equation*'))
    \startformula 8 \sqrt{2} \mu^{\frac{7}{2}} \stopformula

    Other settings from ``latex()`` are also supported:

    >>> print(context((2*tau)**Rational(7,2), fold_frac_powers=True))
    8 \sqrt{2} \tau^{7/2}
    >>> print(context(3*x**2/y))
    \frac{3 x^{2}}{y}
    >>> print(context(3*x**2/y, fold_short_frac=True))
    3 x^{2} / y
    >>> print(context(Integral(x, x)))
    \int x\, dx

    Notes
    =====

    ConTeXt is a document preparation system that is an alternative to LaTeX.
    While the mathematical syntax is largely the same, ConTeXt uses different
    commands for document structure and environments.

    Key differences from LaTeX:
    - ConTeXt uses ``\startformula...\stopformula`` instead of
      ``\begin{equation}...\end{equation}``
    - ConTeXt uses ``\startformulas...\stopformulas`` instead of
      ``\begin{align}...\end{align}``
    - For numbered equations, ConTeXt uses ``\placeformula\startformula...``

    Most mathematical expressions will have identical output to the LaTeX
    printer since the mathematical syntax is the same.

    See Also
    ========

    latex : Convert expression to LaTeX
    """
    return ContextPrinter(settings).doprint(expr)


def print_context(expr, **settings):
    """Prints ConTeXt representation of the given expression. Takes the same
    settings as ``context()``."""
    print(context(expr, **settings))


# For backward compatibility and convenience
ctext = context
print_ctext = print_context
