
from sympy.core import Basic, S, Symbol
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.integrals.risch import risch_norman
from sympy.series import limit

class Integral(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated integral."""

    precedence = Basic.Apply_precedence

    def __new__(cls, function, *symbols, **assumptions):
        function = Basic.sympify(function)

        if isinstance(function, Basic.Number):
            if isinstance(function, Basic.NaN):
                return S.NaN
            elif isinstance(function, Basic.Infinity):
                return S.Infinity
            elif isinstance(function, Basic.NegativeInfinity):
                return S.NegativeInfinity

        if symbols:
            limits = []

            for V in symbols:
                if isinstance(V, Symbol):
                    limits.append((V,None))
                    continue
                elif isinstance(V, (tuple, list)):
                    if len(V) == 3:
                        limits.append( (V[0],tuple(V[1:])) )
                        continue
                    elif len(V) == 1:
                        if isinstance(V[0], Symbol):
                            limits.append((V[0],None))
                            continue

                raise ValueError("Invalid integration variable or limits")
        else:
            limits = func.atoms(Symbol)

            if not limits:
                return function

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (function, tuple(limits))

        return obj

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1]

    @property
    def variables(self):
        variables = []

        for x,ab in self.limits:
            variables.append(x)

        return variables

    @staticmethod
    def _xab_tostr(xab):
        """str representation of integration variable with optional limits"""
        x,ab = xab
        if ab is None:
            return str(x)
        else:
            return str(xab)


    def tostr(self, level=0):
        L = ', '.join([ self._xab_tostr(l) for l in self.limits ])
        return 'Integral(%s, %s)' % (self.function.tostr(), L)

    def doit(self, **hints):
        if not hints.get('integrals', True):
            return self

        function = self.function

        for x,ab in self.limits:
            antideriv = self._eval_integral(function, x)

            if antideriv is None:
                return self
            else:
                if ab is None:
                    function = antideriv
                else:
                    a,b = ab
                    A = antideriv.subs(x, a)

                    if isinstance(A, Basic.NaN):
                        A = limit(antideriv, x, a)
                    if isinstance(A, Basic.NaN):
                        return self

                    B = antideriv.subs(x, b)

                    if isinstance(B, Basic.NaN):
                        B = limit(antideriv, x, b)
                    if isinstance(B, Basic.NaN):
                        return self

                    function = B - A

        return function

    def _eval_integral(self, f, x):
        # TODO : add table lookup for logarithmic and sine/cosine integrals
        # and for some elementary special cases for speed improvement.
        return risch_norman(f, x)

def integrate(*args, **kwargs):
    """Compute definite or indefinite integral of one or more variables
       using Risch-Norman algorithm and table lookup. This procedure is
       able to handle elementary algebraic and transcendental functions
       and also a huge class of special functions, including Airy,
       Bessel, Whittaker and Lambert.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> integrate(log(x), x)
       -x + x*log(x)

    """
    integral = Integral(*args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit()
    else:
        return integral
