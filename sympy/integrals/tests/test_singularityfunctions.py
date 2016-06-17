from sympy.integrals.singularityfunctions import singularityintegrate
from sympy import SingularityFunction, symbols, DiracDelta, Heaviside

x, a, n = symbols('x a n')


def test_singularityintegrate():
    assert singularityintegrate(x, x) is None

    assert 4*singularityintegrate(SingularityFunction(x, a, 3), x) == 4*SingularityFunction(x, a, 4)/4
    assert 5*singularityintegrate(SingularityFunction(x, 5, -2), x) == 5*SingularityFunction(x, 5, -1)
    assert 6*singularityintegrate(SingularityFunction(x, 5, -1), x) == 6*SingularityFunction(x, 5, 0)
