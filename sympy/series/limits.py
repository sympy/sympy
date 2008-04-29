from sympy.core import S, Add, sympify, Basic
from sympy.core.methods import NoRelMeths, ArithMeths
from gruntz import gruntz

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

    Examples:

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
        if e.is_polynomial():
            return Add(*[limit(term, z, z0, dir) for term in e.args])
        else:
            # this is a case like limit(x*y+x*z, z, 2) == x*y+2*x
            # but we need to make sure, that the general gruntz() algorithm is
            # executed for a case like "limit(sqrt(x+1)-sqrt(x),x,oo)==0"
            r = e.subs(z, z0)
            if r is not S.NaN:
                return r

    return gruntz(e, z, z0, dir)

class Limit(Basic, NoRelMeths, ArithMeths):
    """Represents unevaluated limit.

    Examples:

    >>> Limit(sin(x)/x, x, 0)
    Limit(1/x*sin(x), x, 0, dir='+')
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

    def tostr(self, level=0):
        e, z, z0, dir = self.args
        return "Limit(%s, %s, %s, dir='%s')" % (e, z, z0, dir)

