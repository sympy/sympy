from __future__ import print_function, division

from sympy.core import S, Symbol, Add, sympify, Expr, PoleError, Mul
from sympy.core.compatibility import string_types
from sympy.core.symbol import Dummy
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.numbers import GoldenRatio, Float, Rational
from sympy.functions.combinatorial.numbers import fibonacci
from sympy.functions.special.gamma_functions import gamma
from sympy.series.order import Order
from .gruntz import gruntz, limitinf
from sympy.core.exprtools import factor_terms
from sympy.simplify.ratsimp import ratsimp
from sympy.polys import PolynomialError

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-).  For infinite
    z0 (oo or -oo), the dir argument is determined from the direction
    of the infinity (i.e., dir="-" for oo).

    Examples
    ========

    >>> from sympy import limit, sin, Symbol, oo
    >>> from sympy.abc import x
    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, 0, dir="+")
    oo
    >>> limit(1/x, x, 0, dir="-")
    -oo
    >>> limit(1/x, x, oo)
    0

    Notes
    =====

    First we try some heuristics for easy and frequent cases like "x", "1/x",
    "x**2" and similar, so that it's fast. For all other cases, we use the
    Gruntz algorithm (see the gruntz() function).
    """

    return Limit(e, z, z0, dir).doit(deep=False)


def heuristics(e, z, z0, dir):
    rv = None
    if abs(z0) is S.Infinity:
        rv = limit(e.subs(z, 1/z), z, S.Zero, "+" if z0 is S.Infinity else "-")
        if isinstance(rv, Limit):
            return
    elif e.is_Mul or e.is_Add or e.is_Pow or e.is_Function:
        r = []
        for a in e.args:
            l = limit(a, z, z0, dir)
            if l.has(S.Infinity) and l.is_finite is None:
                return
            elif isinstance(l, Limit):
                return
            elif l is S.NaN:
                return
            else:
                r.append(l)
        if r:
            rv = e.func(*r)
            if rv is S.NaN:
                try:
                    rat_e = ratsimp(e)
                except PolynomialError:
                    return
                if rat_e is S.NaN or rat_e == e:
                    return
                return limit(rat_e, z, z0, dir)
    return rv


class Limit(Expr):
    """Represents an unevaluated limit.

    Examples
    ========

    >>> from sympy import Limit, sin
    >>> from sympy.abc import x
    >>> Limit(sin(x)/x, x, 0)
    Limit(sin(x)/x, x, 0)
    >>> Limit(1/x, x, 0, dir="-")
    Limit(1/x, x, 0, dir='-')

    """

    def __new__(cls, e, z, z0, dir="+"):
        e = sympify(e)
        z = sympify(z)
        z0 = sympify(z0)

        if z0 is S.Infinity:
            dir = "-"
        elif z0 is S.NegativeInfinity:
            dir = "+"

        if isinstance(dir, string_types):
            dir = Symbol(dir)
        elif not isinstance(dir, Symbol):
            raise TypeError("direction must be of type basestring or Symbol, not %s" % type(dir))
        if str(dir) not in ('+', '-'):
            raise ValueError(
                "direction must be either '+' or '-', not %s" % dir)

        obj = Expr.__new__(cls)
        obj._args = (e, z, z0, dir)
        return obj

    @property
    def free_symbols(self):
        e = self.args[0]
        isyms = e.free_symbols
        isyms.difference_update(self.args[1].free_symbols)
        isyms.update(self.args[2].free_symbols)
        return isyms

    def doit(self, **hints):
        """Evaluates limit.

        Notes
        =====

        First we handle some trivial cases (i.e. constant), then try
        Gruntz algorithm (see the :py:mod:`~sympy.series.gruntz` module).
        """
        from sympy.series.limitseq import limit_seq
        from sympy.functions import RisingFactorial

        e, z, z0, dir = self.args

        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)

        has_Floats = e.has(Float)
        if has_Floats:
            e = e.subs({k: Rational(k) for k in e.atoms(Float)},
                       simultaneous=True)

        if z0.has(z):
            newz = z.as_dummy()
            r = limit(e.subs(z, newz), newz, z0, dir)
            if r.func is Limit:
                r = r.subs(newz, z)
            return r

        if e == z:
            return z0

        if not e.has(z):
            return e

        if z0 is S.NaN:
            return S.NaN

        # gruntz fails on factorials but works with the gamma function
        # If no factorial term is present, e should remain unchanged.
        # factorial is defined to be zero for negative inputs (which
        # differs from gamma) so only rewrite for positive z0.
        if z0.is_positive:
            e = e.rewrite([factorial, RisingFactorial], gamma)

        if e.has(Order):
            e = e.expand()
            order = e.getO()
            if order:
                if (z, z0) in zip(order.variables, order.point):
                    order = limit(order.expr, z, z0, dir)
                    e = e.removeO() + order

        try:
            # Convert to the limit z->oo and use Gruntz algorithm.
            newe, newz = e, z
            if z0 == S.NegativeInfinity:
                newe = e.subs(z, -z)
            elif z0 != S.Infinity:
                if str(dir) == "+":
                    newe = e.subs(z, z0 + 1/z)
                else:
                    newe = e.subs(z, z0 - 1/z)

            newe = factor_terms(newe)
            newe = newe.rewrite(fibonacci, GoldenRatio)

            if not z.is_positive or not z.is_finite:
                # We need a fresh variable here to simplify expression further.
                newz = Dummy(z.name, positive=True, finite=True)
                newe = newe.subs(z, newz)

            r = limitinf(newe, newz)
        except (PoleError, ValueError):
            r = heuristics(e, z, z0, dir)
            if r is None:
                return self
        except NotImplementedError:
            # Trying finding limits of sequences
            if hints.get('sequence', True) and z0 is S.Infinity:
                trials = hints.get('trials', 5)
                r = limit_seq(e, z, trials)
                if r is None:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()

        if has_Floats:
            r = r.evalf()

        return r
