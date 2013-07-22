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
    from sympy import Wild, log

    e = sympify(e)
    z = sympify(z)
    z0 = sympify(z0)

    if e == z:
        return z0

    if e.is_Rational:
        return e

    if not e.has(z):
        return e

    # gruntz fails on factorials but works with the gamma function
    # If no factorial term is present, e should remain unchanged.
    # factorial is defined to be zero for negative inputs (which
    # differs from gamma) so only rewrite for positive z0.
    if z0.is_positive:
        e = e.rewrite(factorial, gamma)

    if e.func is tan:
        # discontinuity at odd multiples of pi/2; 0 at even
        disc = S.Pi/2
        sign = 1
        if dir == '-':
            sign *= -1
        i = limit(sign*e.args[0], z, z0)/disc
        if i.is_integer:
            if i.is_even:
                return S.Zero
            elif i.is_odd:
                if dir == '+':
                    return S.NegativeInfinity
                else:
                    return S.Infinity

    if e.func is cot:
        # discontinuity at multiples of pi; 0 at odd pi/2 multiples
        disc = S.Pi
        sign = 1
        if dir == '-':
            sign *= -1
        i = limit(sign*e.args[0], z, z0)/disc
        if i.is_integer:
            if dir == '-':
                return S.NegativeInfinity
            else:
                return S.Infinity
        elif (2*i).is_integer:
            return S.Zero

    if e.is_Pow:
        b, ex = e.args
        c = None  # records sign of b if b is +/-z or has a bounded value
        if b.is_Mul:
            c, b = b.as_two_terms()
            if c is S.NegativeOne and b == z:
                c = '-'
        elif b == z:
            c = '+'

        if ex.is_number:
            if c is None:
                base = b.subs(z, z0)
                if base != 0 and (ex.is_bounded or base is not S.One):
                    return base**ex
            else:
                if z0 == 0 and ex < 0:
                    if dir != c:
                        # integer
                        if ex.is_even:
                            return S.Infinity
                        elif ex.is_odd:
                            return S.NegativeInfinity
                        # rational
                        elif ex.is_Rational:
                            return (S.NegativeOne**ex)*S.Infinity
                        else:
                            return S.ComplexInfinity
                    return S.Infinity
                return z0**ex

    if e.is_Mul or not z0 and e.is_Pow and b.func is log:
        if e.is_Mul:
            if abs(z0) is S.Infinity:
                n, d = e.as_numer_denom()
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
                if all(ok(w) for w in (n, d)):
                    u = C.Dummy(positive=(z0 is S.Infinity))
                    inve = (n/d).subs(z, 1/u)
                    return limit(inve.as_leading_term(u), u,
                        S.Zero, "+" if z0 is S.Infinity else "-")

            # weed out the z-independent terms
            i, d = e.as_independent(z)
            if i is not S.One and i.is_bounded:
                return i*limit(d, z, z0, dir)
        else:
            i, d = S.One, e
        if not z0:
            # look for log(z)**q or z**p*log(z)**q
            p, q = Wild("p"), Wild("q")
            r = d.match(z**p * log(z)**q)
            if r:
                p, q = [r.get(w, w) for w in [p, q]]
                if q and q.is_number and p.is_number:
                    if q > 0:
                        if p > 0:
                            return S.Zero
                        else:
                            return -oo*i
                    else:
                        if p >= 0:
                            return S.Zero
                        else:
                            return -oo*i

    if e.is_Add:
        if e.is_polynomial():
            if not z0.is_unbounded:
                return Add(*[limit(term, z, z0, dir) for term in e.args])
        elif e.is_rational_function(z):
            rval = Add(*[limit(term, z, z0, dir) for term in e.args])
            if rval != S.NaN:
                return rval
        if not any([a.is_unbounded for a in e.args]):
            e = e.normal() # workaround for issue 3744

    if e.is_Order:
        args = e.args
        return C.Order(limit(args[0], z, z0), *args[1:])

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
    bad = (S.Infinity, S.NegativeInfinity, S.NaN, None)
    if e.is_Mul:
        r = []
        for a in e.args:
            if not a.is_bounded:
                r.append(a.limit(z, z0, dir))
                if r[-1] in bad:
                    break
        else:
            if r:
                rv = Mul(*r)
    if rv is None and (e.is_Add or e.is_Pow or e.is_Function):
        rv = e.func(*[limit(a, z, z0, dir) for a in e.args])
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
