from sympy.core import Basic
from printer import Printer


class LatexPrinter(Printer):
    """A printer which converts an expression into its LaTeX equivalent."""
    def doprint(self, e):
        return "$%s$" % Printer.doprint(self, e)

    def _print_Mul(self, e):
        # TODO Improve this to collect numerator/denominator and possibly
        #      use \frac{num}{den}?
        f = []
        coeff, a = e.as_coeff_terms()
        if isinstance(coeff, Basic.NegativeOne):
            f = ["-"]

        multsymb = r"\cdot"
        for x in a:
            xs = self._print(x)
            if isinstance(x, Basic.AssocOp):
                f.append("(%s)" % xs)

            # Insert explicit multiplication sign in some cases
            elif (isinstance(x, Basic.Symbol) and "mathrm" in xs) or \
                 (isinstance(x, Basic.Pow) and "mathrm" in self._print(x.base)):
                f.extend([multsymb, xs, multsymb])
            else:
                f.append(xs)

        # Remove extra multiplication signs
        for i in range(len(f)-1):
            if f[i] == f[i+1] == multsymb:
                f[i] = ""
        f = [x for x in f if x]
        if f[0] == multsymb: f = f[1:]
        if f[-1] == multsymb: f = f[:-1]
        return str.join(" ", f)

    def _print_Add(self, e):
        f = "%s" % self._print(e[0])
        for term in e[1:]:
            num_part,other = term.as_coeff_terms()
            if num_part < 0:
              f += "%s" % self._print(term)
            else:
              f += "+%s" % self._print(term)
        return f

    def _print_Integral(self, e):
        if e.a is None:
            # if this is an indefinite integral
            return "\int %s\,d%s" % (self._print(e.f), self._print(e.x))
        else:
            return "\int^%s_%s %s\,d%s" % (self._print(e.a), self._print(e.b),
                                           self._print(e.f), self._print(e.x))

    def _print_Apply(self, e):
        # Check to see if there is something here for this func first
        func = e.func.__class__.__name__
        if hasattr(self, '_print_'+func):
            return getattr(self, '_print_'+func)(e)

        # No handler function, do a generic apply
        args = []
        for arg in e.args:
            args.append(self._print(arg))

        s = self._print(Basic.Symbol(e.func.name))
        s += r"\left(%s\right)" % str.join(',', args)
        return s

    def _print_ApplyExp(self, e):
        return "{e}^{%s}" % self._print(e.args[0])

    def _print_Derivative(self, e):
        # TODO Upgrade for multiple symbols used in differentiation
        x = e.symbols[0]
        s = r"\frac{\partial}{\partial %s} " % self._print(x)
        if isinstance(e.expr, Basic.Add):
            s += r"\left(" + self._print(e.expr) + r"\right)"
        else:
            s += self._print(e.expr)
        return s

    def _print_Factorial(self, e):
        x = e.args[0]
        if (isinstance(x, Basic.Integer) and x.is_nonnegative) or \
            isinstance(x, Basic.Symbol):
            s = self._print(x)
        else:
            s = "(" + self._print(x) + ")"
        return s + "!"

    def _print_RisingFactorial(self, e):
        x, n = e.args
        return "{(%s)}^{(%s)}" % (self._print(x), self._print(n))

    def _print_FallingFactorial(self, e):
        x, n = e.args
        return "{(%s)}_{(%s)}" % (self._print(x), self._print(n))

    def _print_Binomial2(self, e):
        n, k = e.args
        return r"{{%s}\choose{%s}}" % (self._print(x), self._print(k))

    def _print_Infinity(self, e):
        return r"\infty"

    def _print_NegativeInfinity(self, e):
        return r"-\infty"

    def _print_ImaginaryUnit(self, e):
        return r"\mathrm{i}"

    def _print_Pi(self, e):
        return r"\pi"

    def _print_Pow(self, e):
        f = ""
        if isinstance(e.base, Basic.AssocOp) or isinstance(e.base, Basic.Pow):
            f += "{(%s)}"
        else:
            f += "{%s}"
        f += "^"
        if isinstance(e.exp, Basic.AssocOp) or isinstance(e.exp, Basic.Pow) \
            or (isinstance(e.exp, Basic.Rational) and \
            (not e.exp.is_integer or (e.exp.is_integer and \
            int(e.exp) < 0)) ):
            f += "{(%s)}"
        else:
            f += "{%s}"
        return f % (self._print(e.base), self._print(e.exp))

    def _print_Symbol(self, e):
        if len(e.name) == 1:
            return e.name
        greek = set(['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
          'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi',
          'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi',
          'psi', 'omega'])
        if e.name.lower() in greek:
            return "\\" + e.name
        return r"\mathrm{%s}" % e.name

    def _print_Relational(self, e):
        charmap = {
            '==': '=',
            '<':  '<',
            '<=': '\leq',
            '!=': '\neq'
        }

        rsym = charmap[e.rel_op]
        return "%s %s %s" % (self._print(e.lhs), rsym, self._print(e.rhs))

def latex(expr):
    """
    Usage
    =====
        Returns the latex code representing the current object.

    Notes
    =====
        @param x: a sympy object. It can be any expression as long 
            as it inherits from basic
        @return: a string with TeX code (something like $\int \left( 1 +yx\right)dx$)
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> from sympy.printing.latex import latex
        >>> print latex( integrate(x*y-2, x, evaluate=False))
        $\int -2+x y\,dx$
    """
    lp = LatexPrinter()
    return lp.doprint(expr)

def print_latex(expr):
    """
    Prints expr in pretty form.

    pprint is just a shortcut for this function
    """
    print latex(expr)