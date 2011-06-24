from sympy.core import S, Add, sympify, Expr, PoleError, Mul, oo, C
from gruntz import gruntz
from sympy.functions import sign, tan, cot

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

    Examples:

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

    Strategy:

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
        c = None # records sign of b if b is +/-z or has a bounded value
        if b.is_Mul:
            c, b = b.as_two_terms()
            if c is S.NegativeOne and b == z:
                c = '-'
        elif b == z:
            c = '+'

        if ex.is_number:
            if c is None:
                base = b.subs(z, z0)
                if base.is_bounded and (ex.is_bounded or base is not S.One):
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
        if e.is_polynomial() and not z0.is_unbounded:
            return Add(*[limit(term, z, z0, dir) for term in e.args])

        # this is a case like limit(x*y+x*z, z, 2) == x*y+2*x
        # but we need to make sure, that the general gruntz() algorithm is
        # executed for a case like "limit(sqrt(x+1)-sqrt(x),x,oo)==0"
        unbounded = []; unbounded_result=[]
        finite = []; unknown = []
        ok = True
        for term in e.args:
            if not term.has(z) and not term.is_unbounded:
                finite.append(term)
                continue
            result = term.subs(z, z0)
            bounded = result.is_bounded
            if bounded is False or result is S.NaN:
                if unknown:
                    ok = False
                    break
                unbounded.append(term)
                if result != S.NaN:
                    # take result from direction given
                    result = limit(term, z, z0, dir)
                unbounded_result.append(result)
            elif bounded:
                finite.append(result)
            else:
                if unbounded:
                    ok = False
                    break
                unknown.append(result)
        if not ok:
            # we won't be able to resolve this with unbounded
            # terms, e.g. Sum(1/k, (k, 1, n)) - log(n) as n -> oo:
            # since the Sum is unevaluated it's boundedness is
            # unknown and the log(n) is oo so you get Sum - oo
            # which is unsatisfactory.
            raise NotImplementedError('unknown boundedness for %s' %
                                      (unknown or result))
        u = Add(*unknown)
        if unbounded:
            inf_limit = Add(*unbounded_result)
            if inf_limit is not S.NaN:
                return inf_limit + u
            if finite:
                return Add(*finite) + limit(Add(*unbounded), z, z0, dir) + u
        else:
            return Add(*finite) + u

    if e.is_Order:
        args = e.args
        return C.Order(limit(args[0], z, z0), *args[1:])

    try:
        r = gruntz(e, z, z0, dir)
        if r is S.NaN:
            raise PoleError()
    except PoleError:
        r = heuristics(e, z, z0, dir)
    return r

def heuristics(e, z, z0, dir):
    if z0 == oo:
        return limit(e.subs(z, 1/z), z, sympify(0), "+")
    elif e.is_Mul:
        r = []
        for a in e.args:
            if not a.is_bounded:
                r.append(a.limit(z, z0, dir))
        if r:
            return Mul(*r)
    elif e.is_Add:
        r = []
        for a in e.args:
            r.append(a.limit(z, z0, dir))
        return Add(*r)
    elif e.is_Function:
        return e.subs(e.args[0], limit(e.args[0], z, z0, dir))
    msg = "Don't know how to calculate the limit(%s, %s, %s, dir=%s), sorry."
    raise PoleError(msg % (e, z, z0, dir))


class Limit(Expr):
    """Represents an unevaluated limit.

    Examples:

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
        obj = Expr.__new__(cls)
        obj._args = (e, z, z0, dir)
        return obj

    def doit(self, **hints):
        e, z, z0, dir = self.args
        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)
        return limit(e, z, z0, dir)
