from sympy import Symbol, Q
from sympy.assumptions.wrapper import (AssumptionsWrapper, is_infinite,
    is_extended_real)


def test_AssumptionsWrapper():
    x = Symbol('x', positive=True)
    y = Symbol('y')
    assert AssumptionsWrapper(x).is_positive
    assert AssumptionsWrapper(y).is_positive is None
    assert AssumptionsWrapper(y, Q.positive(y)).is_positive


def test_is_infinite():
    x = Symbol('x', infinite=True)
    y = Symbol('y', infinite=False)
    z = Symbol('z')
    assert is_infinite(x)
    assert not is_infinite(y)
    assert is_infinite(z) is None
    assert is_infinite(z, Q.infinite(z))


def test_is_extended_real():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=False)
    z = Symbol('z')
    assert is_extended_real(x)
    assert not is_extended_real(y)
    assert is_extended_real(z) is None
    assert is_extended_real(z, Q.extended_real(z))
