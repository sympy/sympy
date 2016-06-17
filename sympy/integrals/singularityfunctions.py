from __future__ import print_function, division

from sympy.functions import SingularityFunction, DiracDelta, Heaviside
from sympy.core import sympify


def singularityintegrate(f, x):
    if not f.has(SingularityFunction):
        return None

    if f.func == SingularityFunction:
        # I dont know how to handle the case for x*SingularityFunction(x, a, n)
        # or ever more complex one where f.is_Mul or f.is_Pow is True. This way
        # we can handle only the most simple form of the Singularity Functions.
        x = sympify(f.args[0])
        a = sympify(f.args[1])
        n = sympify(f.args[2])
        if n.is_positive or n.is_zero:
            return SingularityFunction(x, a, n + 1)/(n + 1)
        elif n == -1 or n == -2:
            # For this case I am using the inbuilt DiracDelta integration of Sympy.
            from sympy.integrals import integrate

            expr = f.rewrite(DiracDelta)
            expr = integrate(expr, x)
            return expr.rewrite(SingularityFunction)

    return None
