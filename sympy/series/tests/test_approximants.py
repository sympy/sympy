from sympy.core.compatibility import range

from sympy.series import approximants
from sympy import fibonacci, lucas


def test_approximants():
    g = [lucas(k) for k in range(16)]
    assert len([x for x in approximants(g)]) == 4
    g = [lucas(k)+fibonacci(k+2) for k in range(16)]
    assert len([x for x in approximants(g)]) == 4
    g = [lucas(k)**2 for k in range(16)]
    assert len([x for x in approximants(g)]) == 6
