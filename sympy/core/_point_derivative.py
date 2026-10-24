"""Experimental exact point derivatives from truncated Taylor coefficients.

This is a small starting point for arithmetic expressions, not a replacement
for ordinary symbolic differentiation or the general series framework.
"""

from math import factorial

from .add import Add
from .singleton import S


class _TaylorNotSupported(Exception):
    """The limited Taylor evaluator cannot handle the expression or point."""


class _TaylorJet:
    """Coefficients of f(x0 + t), truncated after t**n."""

    __slots__ = ('coeffs', 'order')

    def __init__(self, coeffs):
        self.coeffs = coeffs
        self.order = len(coeffs) - 1

    @classmethod
    def constant(cls, value, order):
        return cls([value] + [S.Zero]*order)

    @classmethod
    def variable(cls, value, order):
        if not order:
            return cls([value])
        return cls([value, S.One] + [S.Zero]*(order - 1))

    def __add__(self, other):
        return _TaylorJet([a + b for a, b in zip(self.coeffs, other.coeffs)])

    def __mul__(self, other):
        a, b = self.coeffs, other.coeffs
        return _TaylorJet([
            Add(*(a[k]*b[j - k] for k in range(j + 1)))
            for j in range(self.order + 1)
        ])

    def __pow__(self, exponent):
        if exponent.is_Integer and exponent >= 0:
            n = int(exponent)
            base = self
            result = self.constant(S.One, self.order)
            while n:
                if n & 1:
                    result = result*base
                n >>= 1
                if n:
                    base = base*base
            return result

        a = self.coeffs
        if a[0].is_zero is not False:
            raise _TaylorNotSupported('nonregular power at evaluation point')
        result = [a[0]**exponent] + [S.Zero]*self.order
        # a*h' = exponent*a'*h, where h = a**exponent.
        for j in range(1, self.order + 1):
            result[j] = Add(*(
                ((exponent + 1)*k - j)*a[k]*result[j - k]
                for k in range(1, j + 1)
            ))/(j*a[0])
        return _TaylorJet(result)


def _derivative_at_point(expr, x, x0, n):
    """Return d**n(expr)/dx**n at x0 without constructing the derivative.

    Supports Add/Mul/Pow with numeric exponents and no free parameters apart
    from x. All Taylor coefficients are computed with SymPy arithmetic.
    """
    if not (n.is_Integer and n.is_nonnegative and n.is_number):
        raise _TaylorNotSupported('requires explicit nonnegative order')
    n = int(n)
    if not (x.is_Symbol and x0.is_number and x0.is_finite is True):
        raise _TaylorNotSupported('requires a finite numerical point')
    if expr.is_commutative is not True:
        raise _TaylorNotSupported('noncommutative expression')
    if expr.free_symbols - {x}:
        raise _TaylorNotSupported('unsubstituted parameters')

    cache = {}

    def eval_node(e):
        if e in cache:
            return cache[e]
        if e == x:
            result = _TaylorJet.variable(x0, n)
        elif e.is_number and e.is_finite is True:
            result = _TaylorJet.constant(e, n)
        elif e.is_Add:
            result = _TaylorJet.constant(S.Zero, n)
            for arg in e.args:
                result += eval_node(arg)
        elif e.is_Mul:
            result = _TaylorJet.constant(S.One, n)
            for arg in e.args:
                result *= eval_node(arg)
        elif e.is_Pow and e.exp.is_number:
            result = eval_node(e.base)**e.exp
        else:
            raise _TaylorNotSupported('unsupported expression node')
        cache[e] = result
        return result

    return factorial(n)*eval_node(expr).coeffs[n]


def _subs_derivative_at_point(expr, variables, point):
    """Try a Subs(Derivative(...)) point evaluation; return None to defer."""
    # Delayed import because Derivative is defined in core.function.
    from .function import Derivative

    if not isinstance(expr, Derivative) or len(expr.variable_count) != 1:
        return None
    x, n = expr.variable_count[0]
    if not x.is_Symbol:
        return None
    replacements = dict(zip(variables, point))
    if x not in replacements:
        return None
    x0 = replacements.pop(x)
    if x0.has(x) or any(value.has(x) for value in replacements.values()):
        return None
    try:
        e = expr.expr.subs(replacements)
        n = n.subs(replacements)
        return _derivative_at_point(e, x, x0, n)
    except (_TaylorNotSupported, ValueError, TypeError, ZeroDivisionError):
        return None
