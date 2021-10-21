from .function import Function
from .kind import NumberKind
from .singleton import S


class Rem(Function):
    """
    Parameters
    ==========

    p : Expr
        Dividend.

    q : Expr
        Divisor.

    Notes
    =====

    the remainder always has the same sign as the dividend.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> Rem(x**3,y)
    Rem(x, y)
    >>> _.subs({x: -5, y: 6})
    -5

    """

    kind = NumberKind

    @classmethod
    def eval(cls,p,q):

        def doit(p,q):
            """ the function remainder if both p,q are numbers
                and q is not zero
            """

            if q.is_zero:
                raise ZeroDivisionError("Division by zero")
            if p is S.NaN or q is S.NaN or p.is_finite is False or q.is_finite is False:
                return S.NaN
            if p is S.Zero or p in (q, -q) or (p.is_integer and q == 1):
                return S.Zero

            if q.is_Number:
                if p.is_Number:
                    return p - int(p/q)*q
        rv=doit(p,q)
        if rv is not None:
            return rv

      