from sympy.integrals.singularityfunctions import singularityintegrate
from sympy import SingularityFunction, symbols, DiracDelta, Heaviside

x, a, n = symbols('x a n')


def test_singularityintegrate():
    assert singularityintegrate(SingularityFunction(x, a, 3), x) == SingularityFunction(x, a, 4)/4
    assert singularityintegrate(SingularityFunction(x, a, -1), x) == Heaviside(x - a)
