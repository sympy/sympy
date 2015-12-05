# -*- coding: utf-8 -*-

from sympy.core.function import UndefinedFunction
from sympy.core.symbol import Symbol
from sympy.printing.conventions import split_super_sub
from sympy.printing.latex import LatexPrinter, translate
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.str import StrPrinter

__all__ = ['vprint', 'vsstrrepr', 'vsprint', 'vpprint', 'vlatex']


class VectorStrPrinter(StrPrinter):
    """String Printer for vector expressions. """

    def _print_Function(self, e):
        return e.func.__name__ + "(%s)" % self.stringify(e.args, ", ")


class VectorStrReprPrinter(VectorStrPrinter):
    """String repr printer for vector expressions."""

    def _print_str(self, s):
        return repr(s)


class VectorLatexPrinter(LatexPrinter):
    """Latex Printer for vector expressions. """

    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__
        
        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [str(self._print(arg)) for arg in expr.args]
            # How inverse trig functions should be displayed, formats are:
            # abbreviated: asin, full: arcsin, power: sin^-1
            inv_trig_style = self._settings['inv_trig_style']
            # If we are dealing with a power-style inverse trig function
            inv_trig_power_case = False
            # If it is applicable to fold the argument brackets
            can_fold_brackets = self._settings['fold_func_brackets'] and \
                len(args) == 1 and \
                not self._needs_function_brackets(expr.args[0])

            inv_trig_table = ["asin", "acos", "atan", "acot"]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    func = func
                elif inv_trig_style == "full":
                    func = "arc" + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                name = r"\operatorname{%s}^{%s}" % (func, exp)
            else:
                name = r"\operatorname{%s}" % func

            if can_fold_brackets:
                name += r"%s"
            else:
                name += r"\left(%s\right)"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            return name % ",".join(args)


class VectorPrettyPrinter(PrettyPrinter):
    """Pretty Printer for vectorialexpressions. """

    def _print_Function(self, e):
        # XXX works only for applied functions
        func_name = func.__name__
        pform = self._print_Symbol(Symbol(func_name))
        return pform


def vprint(expr, **settings):
    r"""Function for printing of expressions generated in the
    sympy.vector package.

    Extends SymPy's StrPrinter, takes the same setting accepted by SymPy's
    `sstr()`, and is equivalent to `print(sstr(foo))`.

    Parameters
    ==========

    expr : valid SymPy object
        SymPy expression to print.
    settings : args
        Same as the settings accepted by SymPy's sstr().

    Examples
    ========

    >>> from sympy.vector import vprint, CoordSysCartesian, Vector
    >>> N = CoordSysCartesian('N')
    >>> u1 = N.i + N.j + Vector.zero
    >>> vprint(u1)
    N.i + N.j

    """

    outstr = vsprint(expr, **settings)

    from sympy.core.compatibility import builtins
    if (outstr != 'None'):
        builtins._ = outstr
        print(outstr)


def vsstrrepr(expr, **settings):
    """Function for displaying expression representation's with vector
    printing enabled.

    Parameters
    ==========

    expr : valid SymPy object
        SymPy expression to print.
    settings : args
        Same as the settings accepted by SymPy's sstrrepr().

    """
    p = VectorStrReprPrinter(settings)
    return p.doprint(expr)


def vsprint(expr, **settings):
    r"""Function for displaying expressions generated in the
    sympy.vector package.

    Returns the output of vprint() as a string.

    Parameters
    ==========

    expr : valid SymPy object
        SymPy expression to print
    settings : args
        Same as the settings accepted by SymPy's sstr().

    Examples
    ========

    >>> from sympy.vector import vsprint, CoordSysCartesian
    >>> from sympy.core.symbol import Symbol
    >>> N = CoordSysCartesian('N')
    >>> u1, u2, u3 = Symbol('x')*N.i, Symbol('y')*N.j, Symbol('z')*N.i
    >>> print("%s = %s" % (vsprint(u1), vsprint(u2 + u3)))
    x*N.i = z*N.i + y*N.j

    """

    string_printer = VectorStrPrinter(settings)
    return string_printer.doprint(expr)


def vpprint(expr, **settings):
    r"""Function for pretty printing of expressions generated in the
    sympy.vector package.

    Mainly used for expressions not inside a vector; the output of running
    scripts and generating equations of motion. Takes the same options as
    SymPy's pretty_print(); see that function for more information.

    Parameters
    ==========

    expr : valid SymPy object
        SymPy expression to pretty print
    settings : args
        Same as those accepted by SymPy's pretty_print.
    """

    pp = VectorPrettyPrinter(settings)

    # Note that this is copied from sympy.printing.pretty.pretty_print:

    # XXX: this is an ugly hack, but at least it works
    use_unicode = pp._settings['use_unicode']
    from sympy.printing.pretty.pretty_symbology import pretty_use_unicode
    uflag = pretty_use_unicode(use_unicode)

    try:
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)


def vlatex(expr, **settings):
    r"""Function for printing latex representation of sympy.vector
    objects.

    For latex representation of Vectors and Dyadics. Takes the
    same options as SymPy's latex(); see that function for more information;

    Parameters
    ==========

    expr : valid SymPy object
        SymPy expression to represent in LaTeX form
    settings : args
        Same as latex()

    Examples
    ========

    >>> from sympy.vector import vlatex, CoordSysCartesian
    >>> N = CoordSysCartesian('N')
    >>> vlatex(N.i + N.j)
    '\\mathbf{\\hat{i}_{N}} + \\mathbf{\\hat{j}_{N}}'
    """
    latex_printer = VectorLatexPrinter(settings)

    return latex_printer.doprint(expr)
