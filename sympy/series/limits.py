from __future__ import print_function, division

from sympy.core import S, Symbol, Add, sympify, Expr, PoleError, Mul, oo, C
from sympy.core.compatibility import string_types
from sympy.functions import tan, cot, factorial, gamma
from .gruntz import gruntz


def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

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
    e = sympify(e)
    z = sympify(z)
    z0 = sympify(z0)

    if e == z:
        return z0

    if not e.has(z):
        return e

    # gruntz fails on factorials but works with the gamma function
    # If no factorial term is present, e should remain unchanged.
    # factorial is defined to be zero for negative inputs (which
    # differs from gamma) so only rewrite for positive z0.
    if z0.is_positive:
        e = e.rewrite(factorial, gamma)

    if e.is_Mul:
        if abs(z0) is S.Infinity:
            # XXX todo: this should probably be stated in the
            # negative -- i.e. to exclude expressions that should
            # not be handled this way but I'm not sure what that
            # condition is; when ok is True it means that the leading
            # term approach is going to succeed (hopefully)
            ok = lambda w: (z in w.free_symbols and
                 any(a.is_polynomial(z) or
                 any(z in m.free_symbols and m.is_polynomial(z)
                 for m in Mul.make_args(a))
                 for a in Add.make_args(w)))
            if all(ok(w) for w in e.as_numer_denom()):
                u = C.Dummy(positive=(z0 is S.Infinity))
                inve = e.subs(z, 1/u)
                return limit(inve.as_leading_term(u), u,
                    S.Zero, "+" if z0 is S.Infinity else "-")

    if e.is_Order:
        return C.Order(limit(e.expr, z, z0), *e.args[1:])

    try:
        r = gruntz(e, z, z0, dir)
        if r is S.NaN:
            raise PoleError()
    except (PoleError, ValueError):
        r = heuristics(e, z, z0, dir)
    return r


def heuristics(e, z, z0, dir):
    if abs(z0) is S.Infinity:
        return limit(e.subs(z, 1/z), z, S.Zero, "+" if z0 is S.Infinity else "-")

    rv = None
    bad = (S.NaN, None)

    if e.is_Mul or e.is_Add or e.is_Pow or e.is_Function:
        r = []
        for a in e.args:
            try:
                r.append(limit(a, z, z0, dir))
            except PoleError:
                break
            if r[-1] in bad:
                break
        else:
            if r:
                rv = e.func(*r)

    if rv in bad:
        msg = "Don't know how to calculate the limit(%s, %s, %s, dir=%s), sorry."
        raise PoleError(msg % (e, z, z0, dir))

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

    def doit(self, **hints):
        """Evaluates limit"""
        e, z, z0, dir = self.args
        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)
        return limit(e, z, z0, str(dir))
