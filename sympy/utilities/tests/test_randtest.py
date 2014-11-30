from sympy.utilities.randtest import numerically_zero
from sympy import cos, sin, sqrt, solve, pi, Add, S
from sympy.abc import a, b, c, x, y, z

from sympy.utilities.pytest import slow


def is_zero(x, **args):
    args.setdefault('number', 1)
    return numerically_zero(x, **args)


def test_numerically_zero():
    assert is_zero(0) is True
    assert is_zero(1) is False
    assert is_zero(pi) is False
    assert is_zero(S.NaN) is False
    assert is_zero(sin(x)**2 + cos(x)**2 - 1, maxprec=10) is None
    eps = lambda e: Add(pi, 10**S(-e), -pi, evaluate=False)
    assert is_zero(eps(-9)) is False
    # the expression is <= tol
    assert is_zero(eps(-10)) is False
    assert is_zero(S(10)**(-10)) is True  # a Number can return True
    # one must not trust this blindly
    eq = 100 - x - sqrt((100 - x)**2)  # nonzero for x > 100
    assert is_zero(eq, maxprec=10) is None
    assert is_zero(eq, a=101, c=200, b=0, d=0, maxprec=10) is False
    eq = -x - 1 + (-x**2 + 1)/(-x + 1)  # undefined at x = 1
    assert is_zero(eq, maxprec=10) is None
    assert is_zero(eq, prec=None, maxprec=10) is None


@slow
def test_numerically_zero_slow():
    # 8516
    eqs, syms = [x + y + z - a, x*y + y*z + x*z - b, x*y*z - c], [x, y, z]
    sol = solve(eqs, syms, dict=True, check=False, simplify=False, manual=True)
    assert [is_zero(eqs[0].xreplace(i), number=1, maxprec=10, b=0, d=0) for i in sol] == \
        [False, False] + [True]*6
    assert is_zero(eqs[1].xreplace(sol[2]), number=1, maxprec=10, b=0, d=0) is None
    assert is_zero(eqs[2].xreplace(sol[2]), number=1, maxprec=10, b=0, d=0) is None
