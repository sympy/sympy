from sympy.core import S, Add, sympify, Basic, PoleError, Mul, oo, C
from gruntz import gruntz

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

    Examples:

    >>> from sympy import limit, sin, Symbol
    >>> x = Symbol('x')
    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, 0, dir="+")
    oo
    >>> limit(1/x, x, 0, dir="-")
    -oo
    >>> limit(1/x, x, oo)
    0

    Strategy:

    First we try some heuristics for easy and frequent cases like "x", "1/x",
    "x**2" and similar, so that it's fast. For all other cases, we use the
    Gruntz algorithm (see the gruntz() function).
    """
    e = sympify(e)
    z = sympify(z)
    z0 = sympify(z0)

    if e == z:
        return z0

    if e.is_Rational:
        return e

    if e.is_Pow:
        if e.args[0] == z:
            if e.args[1].is_Rational:
                if e.args[1] > 0:
                    return z0**e.args[1]
                else:
                    if z0 == 0:
                        if dir == "+":
                            return S.Infinity
                        else:
                            return -S.Infinity
                    else:
                        return z0**e.args[1]
            if e.args[1].is_number:
                if e.args[1].evalf() > 0:
                    return S.Zero
                else:
                    if dir == "+":
                        return S.Infinity
                    else:
                        return -S.Infinity

    if e.is_Add:
        if e.is_polynomial() and z0.is_finite:
            return Add(*[limit(term, z, z0, dir) for term in e.args])
        else:
            # this is a case like limit(x*y+x*z, z, 2) == x*y+2*x
            # but we need to make sure, that the general gruntz() algorithm is
            # executed for a case like "limit(sqrt(x+1)-sqrt(x),x,oo)==0"
            unbounded = []; unbounded_result=[]
            finite = []
            for term in e.args:
                result = term.subs(z, z0)
                if result.is_unbounded or result is S.NaN:
                    unbounded.append(term)
                    unbounded_result.append(result)
                else:
                    finite.append(result)
            if unbounded:
                inf_limit = Add(*unbounded_result)
                if inf_limit is not S.NaN:
                    return inf_limit
                if finite:
                    return Add(*finite) + limit(Add(*unbounded), z, z0, dir)
            else:
                return Add(*finite)

    try:
        r = gruntz(e, z, z0, dir)
    except PoleError:
        r = heuristics(e, z, z0, dir)
    return r

def heuristics(e, z, z0, dir):
    if z0 == oo:
        return heuristics(e.subs(z, 1/z), z, sympify(0), "+")
    elif e.is_Mul:
        r = []
        for a in e.args:
            if not a.is_bounded:
                r.append(a.limit(z, z0, dir))
        if not (r is []):
            return Mul(*r)
    elif e.is_Add:
        r = []
        for a in e.args:
            r.append(a.limit(z, z0, dir))
        return Add(*r)
    elif isinstance(e, C.Function):
        return e.subs(e.args[0], heuristics(e.args[0], z, z0, dir))
    msg = "Don't know how to calculate the limit(%s, %s, %s, dir=%s), sorry."
    raise PoleError(msg % (e, z, z0, dir))


class Limit(Basic):
    """Represents unevaluated limit.

    Examples:

    >>> from sympy import limit, sin, Symbol
    >>> x = Symbol('x')
    >>> Limit(sin(x)/x, x, 0)
    Limit(sin(x)/x, x, 0)
    >>> Limit(1/x, x, 0, dir="-")
    Limit(1/x, x, 0, dir='-')

    """

    def __new__(cls, e, z, z0, dir="+"):
        e = sympify(e)
        z = sympify(z)
        z0 = sympify(z0)
        obj = Basic.__new__(cls)
        obj._args = (e, z, z0, dir)
        return obj

    def doit(self):
        e, z, z0, dir = self.args
        return limit(e, z, z0, dir)
