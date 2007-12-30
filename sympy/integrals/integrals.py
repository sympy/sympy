
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
            # no symbols provided -- let's compute full antiderivative
            limits = [(symb,None) for symb in function.atoms(Symbol)]

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
        """Calculate the antiderivative to the function f(x).

        This is a powerful function that should in theory be able to integrate
        everything that can be integrated. If you find something, that it
        doesn't, it is easy to implement it.

        (1) Simple heuristics (based on pattern matching and integral table):

         - most frequently used functions (eg. polynomials)
         - functions non-integrable by any of the following algorithms (eg.
           exp(-x**2))

        (2) Integration of rational functions:

         (a) using apart() - apart() is full partial fraction decomposition
         procedure based on Bronstein-Salvy algorithm. It gives formal
         decomposition with no polynomial factorization at all (so it's fast
         and gives the most general results). However it needs much better
         implementation of RootsOf class (if fact any implementation).
         (b) using Trager's algorithm - possibly faster than (a) but needs
         implementation :)

        (3) Whichever implementation of pmInt (Mateusz, Kirill's or a
        combination of both).

          - this way we can handle efficiently huge class of elementary and
            special functions

        (4) Recursive Risch algorithm as described in Bronstein's integration
        tutorial.

          - this way we can handle those integrable functions for which (3)
            fails

        (5) Powerful heuristics based mostly on user defined rules.

         - handle complicated, rarely used cases
        """

        # Let's first try some simple functions, that we know fast how to
        # integrate.

        # simple powers:
        from sympy import Pow, log
        if isinstance(f, Pow) and isinstance(f[0], Symbol) and f[1].is_number:
            if f[1] == -1:
                return log(f[0])
            else:
                return f[0]**(f[1]+1)/(f[1]+1) 

        # polynomials:
        from sympy import Polynomial, PolynomialException
        p = None
        try:
            p = Polynomial(f)
        except PolynomialException:
            pass
        if p != None:
            return p.integrate(x)

        # f is not a simple function, let's try the risch norman (that can
        # btw. integrate all the functions above, but slower):
        r = risch_norman(f, x)
        return r

def integrate(*args, **kwargs):
    """integrate(f, var, ...)

       Compute definite or indefinite integral of one or more variables
       using Risch-Norman algorithm and table lookup. This procedure is
       able to handle elementary algebraic and transcendental functions
       and also a huge class of special functions, including Airy,
       Bessel, Whittaker and Lambert.

       var can be:

       - a symbol                   -- indefinite integration
       - a tuple (symbol, a, b)     -- definite integration

       Several variables can be specified, in which case the result is multiple
       integration.

       Also, if no var is specified at all, then full-antiderivative of f is
       returned. This is equivalent of integrating f over all it's variables.

       Examples
       --------

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> integrate(x*y, x)
       (1/2)*y*x**2

       >>> integrate(log(x), x)
       -x + x*log(x)

       >>> integrate(x)
       (1/2)*x**2

       >>> integrate(x*y)
       (1/4)*x**2*y**2

       See also the doctest of Integral._eval_integral(), which explains
       thoroughly the strategy that SymPy uses for integration.

    """
    integral = Integral(*args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit()
    else:
        return integral
