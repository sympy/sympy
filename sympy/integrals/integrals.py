
from sympy.core import Basic, S, C, Symbol, Wild, Pow, sympify
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.integrals.risch import heurisch
from sympy.integrals.trigonometry import trigintegrate
from sympy.simplify import apart
from sympy.series import limit
from sympy.polys import Poly

class Integral(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated integral."""

    precedence = Basic.Apply_precedence

    def __new__(cls, function, *symbols, **assumptions):
        function = sympify(function)

        if function.is_Number:
            if function is S.NaN:
                return S.NaN
            elif function is S.Infinity:
                return S.Infinity
            elif function is S.NegativeInfinity:
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

                    if A is S.NaN:
                        A = limit(antideriv, x, a)
                    if A is S.NaN:
                        return self

                    B = antideriv.subs(x, b)

                    if B is S.NaN:
                        B = limit(antideriv, x, b)
                    if B is S.NaN:
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

        # if it is a poly(x) then let the polynomial integrate itself (fast)
        #
        # It is important to make this check first, otherwise the other code
        # will return a sympy expression instead of a Polynomial.
        #
        # see Polynomial for details.
        if isinstance(f, Poly):
            return f.integrate(x)

        # let's cut it short if `f` does not depend on `x`
        if not f.has(x):
            return f*x

        # try to convert to poly(x) and then integrate if successful (fast)
        poly = f.as_poly(x)

        if poly is not None:
            return poly.integrate(x).as_basic()

        # since Integral(f=g1+g2+...) == Integral(g1) + Integral(g2) + ...
        # we are going to handle Add terms separately,
        # if `f` is not Add -- we only have one term
        if not f.is_Add:
            f = [f]

        parts = []

        if isinstance(f, Basic):
            f = f.args
        for g in f:
            coeff, g = g.as_independent(x)

            # g(x) = const
            if g is S.One:
                parts.append(coeff * x)
                continue

            #               c
            # g(x) = (a*x+b)
            if g.is_Pow and not g.exp.has(x):
                a = Wild('a', exclude=[x])
                b = Wild('b', exclude=[x])

                M = g.base.match(a*x + b)

                if M is not None:
                    if g.exp == -1:
                        h = C.log(g.base)
                    else:
                        h = g.base**(g.exp+1) / (g.exp+1)

                    parts.append(coeff * h / M[a])
                    continue

            #        poly(x)
            # g(x) = -------
            #        poly(x)
            if g.is_fraction(x):
                h = self._eval_integral(apart(g, x), x)
                parts.append(coeff * h)
                continue

            # g(x) = Mul(trig)
            h = trigintegrate(g, x)
            if h is not None:
                parts.append(coeff * h)
                continue

            # fall back to the more general algorithm
            h = heurisch(g, x, hints=[])

            if h is not None:
                parts.append(coeff * h)
            else:
                return None

        return C.Add(*parts)

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
