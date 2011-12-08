from function import Function

class Mod(Function):
    '''Represents a modulo operation on symbolic expressions.

    Receives two arguments, dividend p and divisor q.

    The convention used is the same as python's: the remainder always has the
    same sign as the divisor.

    Examples
    --------

    >>> from sympy.abc import x, y
    >>> x**2 % y
    Mod(x**2, y)
    >>> _.subs({x: 5, y: 6})
    1

    '''
    nargs = 2

    @classmethod
    def eval(cls, p, q):
        if p.is_Number and q.is_Number:
            return p % q
        return Mod(p, q, evaluate=False)
