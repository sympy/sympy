"""
A Printer which converts an expression into its LaTeX equivalent.
"""

from sympy.core import S, C, Basic, Symbol
from printer import Printer
from sympy.simplify import fraction
import re

class LatexPrinter(Printer):
    printmethod = "_latex_"

    def __init__(self, inline=True):
        Printer.__init__(self)
        self._inline = inline

    def doprint(self, expr):
        tex = Printer.doprint(self, expr)

        if self._inline:
            return r"$%s$" % tex
        else:
            return r"\begin{equation*}%s\end{equation*}" % tex

    def _needs_brackets(self, expr):
        return not ((expr.is_Integer and expr.is_nonnegative) or expr.is_Atom)

    def _do_exponent(self, expr, exp):
        if exp is not None:
            return r"\left(%s\right)^{%s}" % (expr, exp)
        else:
            return expr

    def _print_Add(self, expr):
        args = list(expr.args)
        args.sort(Basic._compare_pretty)

        tex = str(self._print(args[0]))

        for term in args[1:]:
            coeff = term.as_coeff_terms()[0]

            if coeff.is_negative:
                tex += r" %s" % self._print(term)
            else:
                tex += r" + %s" % self._print(term)

        return tex

    def _print_Mul(self, expr):
        coeff, terms = expr.as_coeff_terms()

        if not coeff.is_negative:
            tex = ""
        else:
            coeff = -coeff
            tex = "- "

        numer, denom = fraction(C.Mul(*terms))

        def convert(terms):
            product = []

            if not terms.is_Mul:
                return str(self._print(terms))
            else:
                for term in terms.args:
                    pretty = self._print(term)

                    if term.is_Add:
                        product.append(r"\left(%s\right)" % pretty)
                    else:
                        product.append(str(pretty))

                return r" ".join(product)

        if denom is S.One:
            if coeff is not S.One:
                tex += str(self._print(coeff)) + " "

            if numer.is_Add:
                tex += r"\left(%s\right)" % convert(numer)
            else:
                tex += r"%s" % convert(numer)
        else:
            if numer is S.One:
                if coeff.is_Integer:
                    numer *= coeff.p
                elif coeff.is_Rational:
                    if coeff.p != 1:
                        numer *= coeff.p

                    denom *= coeff.q
                elif coeff is not S.One:
                    tex += str(self._print(coeff)) + " "
            else:
                if coeff.is_Rational and coeff.p == 1:
                    denom *= coeff.q
                elif coeff is not S.One:
                    tex += str(self._print(coeff)) + " "

            tex += r"\frac{%s}{%s}" % \
                (convert(numer), convert(denom))

        return tex

    def _print_Pow(self, expr):
        if expr.exp.is_Rational and expr.exp.q == 2:
            base, exp = self._print(expr.base), abs(expr.exp.p)

            if exp == 1:
                tex = r"\sqrt{%s}" % base
            else:
                tex = r"\sqrt[%s]{%s}" % (exp, base)

            if expr.exp.is_negative:
                return r"\frac{1}{%s}" % tex
            else:
                return tex
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                if expr.exp == S.NegativeOne:
                    #solves issue 1030
                    #As Mul always simplify 1/x to x**-1
                    #The objective is achieved with this hack
                    #first we get the latex for -1 * expr,
                    #which is a Mul expression
                    tex = self._print(S.NegativeOne * expr).strip()
                    #the result comes with a minus and a space, so we remove
                    if tex[:1] == "-":
                        return tex[1:].strip()
                if self._needs_brackets(expr.base):
                    tex = r"\left(%s\right)^{%s}"
                else:
                    tex = r"{%s}^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Derivative(self, expr):
        dim = len(expr.symbols)

        if dim == 1:
            tex = r"\frac{\partial}{\partial %s}" % \
                self._print(expr.symbols[0])
        else:
            multiplicity, i, tex = [], 1, ""
            current = expr.symbols[0]

            for symbol in expr.symbols[1:]:
                if symbol == current:
                    i = i + 1
                else:
                    multiplicity.append((current, i))
                    current, i = symbol, 1
            else:
                multiplicity.append((current, i))

            for x, i in multiplicity:
                if i == 1:
                    tex += r"\partial %s" % self._print(x)
                else:
                    tex += r"\partial^{%s} %s" % (i, self._print(x))

            tex = r"\frac{\partial^{%s}}{%s} " % (dim, tex)

        if isinstance(expr.expr, C.AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            return r"%s %s" % (tex, self._print(expr.expr))

    def _print_Integral(self, expr):
        tex, symbols = "", []

        for symbol, limits in reversed(expr.limits):
            tex += r"\int"

            if limits is not None:
                if not self._inline:
                    tex += r"\limits"

                tex += "_{%s}^{%s}" % (self._print(limits[0]),
                                       self._print(limits[1]))

            symbols.insert(0, "d%s" % self._print(symbol))

        return r"%s %s\,%s" % (tex,
            str(self._print(expr.function)), " ".join(symbols))

    def _print_Limit(self, expr):
        tex = r"\lim_{%s \to %s}" % (self._print(expr.var),
                                     self._print(expr.varlim))

        if isinstance(expr.expr, C.AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            return r"%s %s" % (tex, self._print(expr.expr))

    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__

        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [ str(self._print(arg)) for arg in expr.args ]

            if exp is not None:
                name = r"\operatorname{%s}^{%s}" % (func, exp)
            else:
                name = r"\operatorname{%s}" % func

            return name + r"\left(%s\right)" % ",".join(args)

    def _print_floor(self, expr, exp=None):
        tex = r"\lfloor{%s}\rfloor" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_ceiling(self, expr, exp=None):
        tex = r"\lceil{%s}\rceil" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_abs(self, expr, exp=None):
        tex = r"\lvert{%s}\rvert" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_re(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Re\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\Re{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_im(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Im\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\Im{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_conjugate(self, expr, exp=None):
        tex = r"\overline{%s}" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_exp(self, expr, exp=None):
        tex = r"{e}^{%s}" % self._print(expr.args[0])
        return self._do_exponent(tex, exp)

    def _print_gamma(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])

        if exp is not None:
            return r"\operatorname{\Gamma}^{%s}%s" % (exp, tex)
        else:
            return r"\operatorname{\Gamma}%s" % tex

    def _print_Factorial(self, expr, exp=None):
        x = expr.args[0]
        if self._needs_brackets(x):
            tex = r"\left(%s\right)!" % self._print(x)
        else:
            tex = self._print(x) + "!"

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_Binomial(self, expr, exp=None):
        tex = r"{{%s}\choose{%s}}" % (self._print(expr[0]),
                                      self._print(expr[1]))

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_RisingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}^{\left(%s\right)}" % \
            (self._print(expr[0]), self._print(expr[1]))

        return self._do_exponent(tex, exp)

    def _print_FallingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}_{\left(%s\right)}" % \
            (self._print(expr[0]), self._print(expr[1]))

        return self._do_exponent(tex, exp)

    def _print_Rational(self, expr):
        if expr.q != 1:
            sign = ""
            p = expr.p
            if expr.p < 0:
                sign = "- "
                p = -p
            return r"%s\frac{%d}{%d}" % (sign, p, expr.q)
        else:
            return self._print(expr.p)

    def _print_Infinity(self, expr):
        return r"\infty"

    def _print_NegativeInfinity(self, expr):
        return r"-\infty"

    def _print_ComplexInfinity(self, expr):
        return r"\tilde{\infty}"

    def _print_ImaginaryUnit(self, expr):
        return r"\mathbf{\imath}"

    def _print_NaN(self, expr):
        return r"\bot"

    def _print_Pi(self, expr):
        return r"\pi"

    def _print_Exp1(self, expr):
        return r"e"

    def _print_EulerGamma(self, expr):
        return r"\gamma"

    def _print_Order(self, expr):
        return r"\operatorname{\mathcal{O}}\left(%s\right)" % \
            self._print(expr.args[0])

    def _print_Symbol(self, expr):
        if len(expr.name) == 1:
            return expr.name

        #convert trailing digits to subscript
        m = re.match('(^[a-zA-Z]+)([0-9]+)$', expr.name)
        if m is not None:
            name, sub=m.groups()
            tex=self._print_Symbol(Symbol(name))
            tex="%s_{%s}" %(tex, sub)
            return tex

        # insert braces to expresions containing '_' or '^'
        m = re.match('(^[a-zA-Z0-9]+)([_\^]{1})([a-zA-Z0-9]+)$', expr.name)
        if m is not None:
            name, sep, rest=m.groups()
            tex=self._print_Symbol(Symbol(name))
            tex="%s%s{%s}" %(tex, sep, rest)
            return tex

        greek = set([ 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                      'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                      'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                      'phi', 'chi', 'psi', 'omega' ])

        other = set( ['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                      'hbar', 'hslash', 'mho' ])

        if expr.name.lower() in greek:
            return "\\" + expr.name
        elif expr.name in other:
            return "\\" + expr.name
        else:
            return expr.name

    def _print_Relational(self, expr):
        charmap = {
            "==" : "=",
            "<"  : "<",
            "<=" : r"\leq",
            "!=" : r"\neq",
        }

        return "%s %s %s" % (self._print(expr.lhs),
            charmap[expr.rel_op], self._print(expr.rhs))

    def _print_Piecewise(self, expr):
        ecpairs = [r"%s & for %s" % (self._print(e), self._print(c)) \
                       for e, c in expr.args[:-1]]
        if expr.args[-1].cond is S.One:
            ecpairs.append(r"%s & \textrm{otherwise}" % \
                               self._print(expr.args[-1].expr))
        else:
            ecpairs.append(r"%s & for %s" % \
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        tex = r"\left\{\begin{array}{cl} %s \end{array}\right."
        return tex % r" \\".join(ecpairs)

    def _print_Matrix(self, expr):
        lines = []

        for line in range(expr.lines): # horrible, should be 'rows'
            lines.append(" & ".join([ self._print(i) for i in expr[line,:] ]))

        if self._inline:
            tex = r"\left(\begin{smallmatrix}%s\end{smallmatrix}\right)"
        else:
            tex = r"\begin{pmatrix}%s\end{pmatrix}"

        return tex % r"\\".join(lines)

    def _print_tuple(self, expr):
        return r"\begin{pmatrix}%s\end{pmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_list(self, expr):
        return r"\begin{bmatrix}%s\end{bmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_dict(self, expr):
        items = []

        keys = expr.keys()
        keys.sort(Basic.compare_pretty)
        for key in keys:
            val = expr[key]
            items.append("%s : %s" % (self._print(key), self._print(val)))

        return r"\begin{Bmatrix}%s\end{Bmatrix}" % r", & ".join(items)

    def _print_DiracDelta(self, expr):
        if len(expr.args) == 1 or expr.args[1] == 0:
            tex = r"\delta\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\delta^{\left( %s \right)}\left( %s \right)" % (\
            self._print(expr.args[1]), self._print(expr.args[0]))
        return tex

def latex(expr, inline=True):
    r"""Convert the given expression to LaTeX representation.

        You can specify how the generated code will be delimited.
        If the 'inline' keyword is set then inline LaTeX $ $ will
        be used. Otherwise the resulting code will be enclosed in
        'equation*' environment (remember to import 'amsmath').

        >>> from sympy import *
        >>> from sympy.abc import *

        >>> latex((2*tau)**Rational(7,2))
        '$8 \\sqrt{2} \\sqrt[7]{\\tau}$'

        >>> latex((2*mu)**Rational(7,2), inline=False)
        '\\begin{equation*}8 \\sqrt{2} \\sqrt[7]{\\mu}\\end{equation*}'

        Besides all Basic based expressions, you can recursively
        convert Pyhon containers (lists, tuples and dicts) and
        also SymPy matrices:

        >>> latex([2/x, y])
        '$\\begin{bmatrix}\\frac{2}{x}, & y\\end{bmatrix}$'

    """

    return LatexPrinter(inline).doprint(expr)

def print_latex(expr):
    """Prints LaTeX representation of the given expression."""
    print latex(expr)
