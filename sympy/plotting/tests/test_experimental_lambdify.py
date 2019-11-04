from sympy.core.symbol import symbols
from sympy.plotting.experimental_lambdify import experimental_lambdify
from sympy.plotting.intervalmath.interval_arithmetic import \
    interval, intervalMembership


def test_composite_boolean_region():
    x, y = symbols('x y')

    r1 = (x - 1)**2 + y**2 < 2
    r2 = (x + 1)**2 + y**2 < 2

    f = experimental_lambdify((x, y), r1 & r2)
    a = (interval(-0.1, 0.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(-1.1, -0.9), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(0.9, 1.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-0.1, 0.1), interval(1.9, 2.1))
    assert f(*a) == intervalMembership(False, True)

    f = experimental_lambdify((x, y), r1 | r2)
    a = (interval(-0.1, 0.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(-1.1, -0.9), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(0.9, 1.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(-0.1, 0.1), interval(1.9, 2.1))
    assert f(*a) == intervalMembership(False, True)

    f = experimental_lambdify((x, y), r1 & ~r2)
    a = (interval(-0.1, 0.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-1.1, -0.9), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(0.9, 1.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(-0.1, 0.1), interval(1.9, 2.1))
    assert f(*a) == intervalMembership(False, True)

    f = experimental_lambdify((x, y), ~r1 & r2)
    a = (interval(-0.1, 0.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-1.1, -0.9), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(True, True)
    a = (interval(0.9, 1.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-0.1, 0.1), interval(1.9, 2.1))
    assert f(*a) == intervalMembership(False, True)

    f = experimental_lambdify((x, y), ~r1 & ~r2)
    a = (interval(-0.1, 0.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-1.1, -0.9), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(0.9, 1.1), interval(-0.1, 0.1))
    assert f(*a) == intervalMembership(False, True)
    a = (interval(-0.1, 0.1), interval(1.9, 2.1))
    assert f(*a) == intervalMembership(True, True)
