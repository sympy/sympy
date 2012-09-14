from function import Function
from sympy.core.numbers import Float
from sympy.core.function import expand_mul

class Mod(Function):
    """Represents a modulo operation on symbolic expressions.

    Receives two arguments, dividend p and divisor q.

    The convention used is the same as python's: the remainder always has the
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
        from sympy.simplify.simplify import nsimplify
        if q.is_Number:
            float = not q.is_Rational
            pnew = expand_mul(p)
            if pnew.is_Number:
                float = float or not pnew.is_Rational
                if not float:
                    return pnew % q
                return Float(nsimplify(pnew) % nsimplify(q))
            elif pnew.is_Add and pnew.args[0].is_Number:
                r, p = pnew.as_two_terms()
                p += Mod(r, q)
        return Mod(p, q, evaluate=False)
