from __future__ import print_function, division

from sympy.functions import SingularityFunction, DiracDelta, Heaviside
from sympy.core import sympify
from sympy.integrals import integrate


def singularityintegrate(f, x):

    if not f.has(SingularityFunction):
        return None

    if f.func == SingularityFunction:
        x = sympify(f.args[0])
        a = sympify(f.args[1])
        n = sympify(f.args[2])
        if n.is_positive or n.is_zero:
            return SingularityFunction(x, a, n + 1)/(n + 1)
        elif n == -1 or n == -2:
            return SingularityFunction(x, a, n + 1)

    if f.is_Mul or f.is_Pow:

        expr = f.rewrite(DiracDelta)
        expr = integrate(expr, x)
        return expr.rewrite(SingularityFunction)
    return None
