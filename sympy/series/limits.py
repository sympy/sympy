from __future__ import print_function, division

from sympy.core import S, Symbol, Add, sympify, Expr, PoleError, Mul
from sympy.core.compatibility import string_types
from sympy.core.symbol import Dummy
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.numbers import GoldenRatio
from sympy.functions.combinatorial.numbers import fibonacci
from sympy.functions.special.gamma_functions import gamma
from sympy.series.order import Order
from .gruntz import gruntz
from sympy.core.exprtools import factor_terms
from sympy.simplify.ratsimp import ratsimp
from sympy.polys import PolynomialError, factor
from sympy.simplify.simplify import together

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of ``e(z)`` at the point ``z0``.

    ``z0`` can be any expression, including ``oo`` and ``-oo``.

    For ``dir="+-"`` it calculates the bi-directional limit; for
    ``dir="+"`` (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-).
    For infinite ``z0`` (``oo`` or ``-oo``), the ``dir`` argument is
    determined from the direction of the infinity (i.e.,
    ``dir="-"`` for ``oo``).

    Examples
    ========

    >>> from sympy import limit, sin, Symbol, oo
    >>> from sympy.abc import x
    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, 0) # default dir='+'
    oo
    >>> limit(1/x, x, 0, dir="-")
    -oo
    >>> limit(1/x, x, 0, dir='+-')
    Traceback (most recent call last):
        ...
    ValueError: The limit does not exist since left hand limit = -oo and right hand limit = oo

    >>> limit(1/x, x, oo)
    0

    Notes
    =====

    First we try some heuristics for easy and frequent cases like "x", "1/x",
    "x**2" and similar, so that it's fast. For all other cases, we use the
    Gruntz algorithm (see the gruntz() function).
    """

    if dir == "+-":
        llim = Limit(e, z, z0, dir="-").doit(deep=False)
        rlim = Limit(e, z, z0, dir="+").doit(deep=False)
        if llim == rlim:
            return rlim
        else:
            # TODO: choose a better error?
            raise ValueError("The limit does not exist since "
                    "left hand limit = %s and right hand limit = %s"
                    % (llim, rlim))
    else:
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
                if isinstance(e, Add):
                    m = factor_terms(e)
                    if not isinstance(m, Mul): # try together
                        m = together(m)
                    if not isinstance(m, Mul): # try factor if the previous methods failed
                        m = factor(e)
                    if isinstance(m, Mul):
                        return heuristics(m, z, z0, dir)
                    return
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

    >>> from sympy import Limit, sin, Symbol
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
            raise TypeError("direction must be of type basestring or "
                    "Symbol, not %s" % type(dir))
        if str(dir) not in ('+', '-', '+-'):
            raise ValueError("direction must be one of '+', '-' "
                    "or '+-', not %s" % dir)

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
        """Evaluates limit"""
        from sympy.series.limitseq import limit_seq
        from sympy.functions import RisingFactorial

        e, z, z0, dir = self.args

        if z0 is S.ComplexInfinity:
            raise NotImplementedError("Limits at complex "
                                    "infinity are not implemented")

        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)

        if e == z:
            return z0

        if not e.has(z):
            return e

        # gruntz fails on factorials but works with the gamma function
        # If no factorial term is present, e should remain unchanged.
        # factorial is defined to be zero for negative inputs (which
        # differs from gamma) so only rewrite for positive z0.
        if z0.is_positive:
            e = e.rewrite([factorial, RisingFactorial], gamma)

        if e.is_Mul:
            if abs(z0) is S.Infinity:
                e = factor_terms(e)
                e = e.rewrite(fibonacci, GoldenRatio)
                ok = lambda w: (z in w.free_symbols and
                                any(a.is_polynomial(z) or
                                    any(z in m.free_symbols and m.is_polynomial(z)
                                        for m in Mul.make_args(a))
                                    for a in Add.make_args(w)))
                if all(ok(w) for w in e.as_numer_denom()):
                    u = Dummy(positive=True)
                    if z0 is S.NegativeInfinity:
                        inve = e.subs(z, -1/u)
                    else:
                        inve = e.subs(z, 1/u)
                    r = limit(inve.as_leading_term(u), u, S.Zero, "+")
                    if isinstance(r, Limit):
                        return self
                    else:
                        return r

        if e.is_Order:
            return Order(limit(e.expr, z, z0), *e.args[1:])

        try:
            r = gruntz(e, z, z0, dir)
            if r is S.NaN:
                raise PoleError()
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

        return r
