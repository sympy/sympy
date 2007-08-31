
from sympy.core import Basic, S, Symbol
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.integrals.risch import risch_norman
from sympy.series import limit

class Integral(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated integral."""

    precedence = Basic.Apply_precedence

    def __new__(cls, func, *symbols, **assumptions):
        func = Basic.sympify(func)

        if isinstance(func, Basic.Number):
            if isinstance(func, Basic.NaN):
                return S.NaN
            elif isinstance(func, Basic.Infinity):
                return S.Infinity
            elif isinstance(func, Basic.NegativeInfinity):
                return S.NegativeInfinity

        variables, limits = [], {}

        if symbols:
            for var in symbols:
                if isinstance(var, Symbol):
                    variables.append(var)
                elif isinstance(var, (tuple, list)) and len(var) == 3:
                    variables.append(var[0])
                    limits[var[0]] = var[1:]
                else:
                    raise ValueError("Invalid arguments")
        else:
            variables = list(func.atoms(Symbol))

        if not variables:
            return func

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (func, variables, limits)

        return obj

    @property
    def function(self):
        return self._args[0]

    @property
    def variables(self):
        return self._args[1]

    @property
    def limits(self):
        return self._args[2]

    @property
    def is_definite(self):
        return not self.is_indefinite

    @property
    def is_indefinite(self):
        return not self.limits

    def tostr(self, level=0):
        elems = []

        for var in self.variables:
            try:
                a, b = self.limits[var]
                elems.append(repr((var, a, b)))
            except KeyError:
                elems.append(var.tostr())

        r = 'Integral(%s)' % ', '.join([self.function.tostr()] + elems)

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

    def doit(self):
        func = self.function

        for i, var in enumerate(self.variables):
            antideriv = self._eval_integral(func, var)

            if antideriv is None:
                return self
            else:
                try:
                    a, b = self.limits[var]

                    # TODO: find singularities

                    Fb = limit(antideriv, var, b)
                    Fa = limit(antideriv, var, a)

                    if isinstance(Fb, Basic.NaN):
                        return self
                    if isinstance(Fa, Basic.NaN):
                        return self

                    func = Fb - Fa
                except KeyError:
                    func = antideriv

        if self.is_indefinite:
            return func# + Symbol('C')
        else:
            return func

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

       #>>> integrate(log(x), x)
       #C - x + x*log(x)

       #>>> integrate(tan(x), x)
       #-log(cos(x)) + C

    """
    integral = Integral(*args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit()
    else:
        return integral
