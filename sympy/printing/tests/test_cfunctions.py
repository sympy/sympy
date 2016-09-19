from sympy import Symbol, exp, log, pi, S
from sympy.printing.cfunctions import expm1, log1p, exp2, log2

def test_expm1():
    # Eval
    assert expm1(0) == 0

    x = Symbol('x', real=True, finite=True)

    # Expand and rewrite
    assert expm1(x).expand(func=True) - exp(x) == -1
    assert expm1(x).rewrite('tractable') - exp(x) == -1
    assert expm1(x).rewrite('exp') - exp(x) == -1

    # Precision
    assert not ((exp(1e-10).evalf() - 1) - 1e-10 - 5e-21) < 1e-22  # for comparison
    assert abs(expm1(1e-10).evalf() - 1e-10 - 5e-21) < 1e-22

    # Properties
    assert expm1(x).is_real
    assert expm1(x).is_finite

    # Diff
    assert expm1(42*x).diff(x) - 42*exp(42*x) == 0


def test_log1p():
    # Eval
    assert log1p(0) == 0
    d = S(10)
    assert log1p(d**-1000) - log(d**1000 + 1) + log(d**1000) == 0

    x = Symbol('x', real=True, finite=True)

    # Expand and rewrite
    assert log1p(x).expand(func=True) - log(x + 1) == 0
    assert log1p(x).rewrite('tractable') - log(x + 1) == 0
    assert log1p(x).rewrite('log') - log(x + 1) == 0

    # Precision
    assert not abs(log(1e-99 + 1).evalf() - 1e-99) < 1e-100  # for comparison
    assert abs(log1p(1e-99).evalf() - 1e-99) < 1e-100

    # Properties
    assert log1p(-2**(-S(1)/2)).is_real

    assert not log1p(-1).is_finite
    assert log1p(pi).is_finite

    assert not log1p(x).is_positive
    assert log1p(Symbol('y', positive=True)).is_positive

    assert not log1p(x).is_zero
    assert log1p(Symbol('z', zero=True)).is_zero

    assert not log1p(x).is_nonnegative
    assert log1p(Symbol('o', nonnegative=True)).is_nonnegative

    # Diff
    assert log1p(42*x).diff(x) - 42/(42*x + 1) == 0


def test_exp2():
    # Eval
    assert exp2(2) == 4

    x = Symbol('x', real=True, finite=True)

    # Expand
    assert exp2(x).expand(func=True) - 2**x == 0

    # Diff
    assert exp2(42*x).diff(x) - 42*exp2(42*x)*log(2) == 0


def test_log2():
    # Eval
    assert log2(8) == 3
    assert log2(pi) != log(pi)/log(2)  # log2 should *save* (CPU) operations

    x = Symbol('x', real=True, finite=True)
    assert log2(x) != log(x)/log(2)
    assert log2(2**x) == x

    # Expand
    assert log2(x).expand(func=True) - log(x)/log(2) == 0

    # Diff
    assert log2(42*x).diff() - 1/(log(2)*x) == 0
