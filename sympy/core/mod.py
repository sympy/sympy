from function import Function


class Mod(Function):
    """Represents a modulo operation on symbolic expressions.

    Receives two arguments, dividend p and divisor q.

    The convention used is the same as Python's: the remainder always has the
    same sign as the divisor.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> x**2 % y
    Mod(x**2, y)
    >>> _.subs({x: 5, y: 6})
    1

    """
    nargs = 2

    @classmethod
    def eval(cls, p, q):
        from sympy.core.singleton import S
        from sympy.core.add import Add
        from sympy.core.exprtools import gcd_terms
        from sympy.polys.polytools import gcd
        from sympy.utilities.iterables import sift

        if p == q:
            return S.Zero

        def doit(p, q):
            """Try to return p % q if both are numbers or +/-p is known
            to be less than q.
            """
            if p.is_Number and q.is_Number:
                return (p % q)
            diff = p - q
            if diff.is_negative:
                if p.is_negative:
                    if (-p < q) is True:
                        return (p + q)
                elif p.is_positive:
                    if (p < q) is True:
                        return p

        rv = doit(p, q)
        if rv is not None:
            return rv

        # remove known terms
        # (x + y + 2) % x -> Mod(y + 2, x)
        if p.is_Add:
            args = [Mod(i, q) for i in p.args]
            if any(a.func is not cls for a in args):
                d = sift(args, lambda x: x.args[1] if x.func is cls else None)
                rv = Add(*d.pop(None, []))
                for k, v in d.iteritems():
                    rv += Mod(Add(*[i.args[0] for i in v]), k, evaluate=False)
                print p, q, rv
                return rv

        # extract gcd; don't do anything else heroic since that should
        # be done as simplification by the user
        G = gcd(p, q)
        if G is not S.One:
            p, q = [
                gcd_terms(i/G, clear=False, fraction=False) for i in (p, q)]

        # handle coefficients if they are not Rational
        # since those are not handled by factor_terms
        # e.g. Mod(.6*x, .3*y) -> 0.3*Mod(2*x, y)
        cp, p = p.as_coeff_Mul()
        cq, q = q.as_coeff_Mul()
        ok = False
        if not cp.is_Rational or not cq.is_Rational:
            r = cp % cq
            if r == 0:
                G *= cq
                p *= int(cp/cq)
                ok = True
        if not ok:
            p = cp*p
            q = cq*q

        # check again to see if p and q can now be handled as
        # numbers
        rv = doit(p, q)
        if rv is not None:
            return G*rv

        # simple -1 extraction
        if p.could_extract_minus_sign() and q.could_extract_minus_sign():
            G, p, q = [-i for i in (G, p, q)]

        return G*Mod(p, q, evaluate=False)
