"""Tests for the implementation of RootOf class and related tools. """

from sympy.polys.rootoftools import RootOf

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
    PolynomialError,
    DomainError,
)

from sympy import sqrt, I, Real

from sympy.utilities.pytest import raises
from sympy.abc import x, y

def test_RootOf___new__():
    assert RootOf(x, 0) == 0
    assert RootOf(x,-1) == 0

    assert RootOf(x - 1, 0) == 1
    assert RootOf(x - 1,-1) == 1

    assert RootOf(x + 1, 0) ==-1
    assert RootOf(x + 1,-1) ==-1

    assert RootOf((x - 1)*(x + 1), 0) ==-1
    assert RootOf((x - 1)*(x + 1), 1) == 1
    assert RootOf((x - 1)*(x + 1),-1) == 1
    assert RootOf((x - 1)*(x + 1),-2) ==-1

    assert RootOf((x - 1)*(x**2 - 2), 0) == RootOf(x**2 - 2, 0)
    assert RootOf((x - 1)*(x**2 - 2), 1) == 1
    assert RootOf((x - 1)*(x**2 - 2), 2) == RootOf(x**2 - 2, 1)
    assert RootOf((x - 1)*(x**2 - 2),-1) == RootOf(x**2 - 2, 1)
    assert RootOf((x - 1)*(x**2 - 2),-2) == 1
    assert RootOf((x - 1)*(x**2 - 2),-3) == RootOf(x**2 - 2, 0)

    raises(GeneratorsNeeded, "RootOf(0, 0)")
    raises(GeneratorsNeeded, "RootOf(1, 0)")

    raises(PolynomialError, "RootOf(x - y, 0)")

    raises(DomainError, "RootOf(x - sqrt(2), 0)")
    raises(DomainError, "RootOf(x - I, 0)")

    raises(IndexError, "RootOf(x**2 - 1,-4)")
    raises(IndexError, "RootOf(x**2 - 1,-3)")
    raises(IndexError, "RootOf(x**2 - 1, 2)")
    raises(IndexError, "RootOf(x**2 - 1, 3)")

def test_RootOf___eq__():
    assert (RootOf(x**3 + 2, 0) == RootOf(x**3 + 2, 0)) == True
    assert (RootOf(x**3 + 2, 0) == RootOf(x**3 + 2, 1)) == False
    assert (RootOf(x**3 + 2, 1) == RootOf(x**3 + 2, 1)) == True
    assert (RootOf(x**3 + 2, 1) == RootOf(x**3 + 2, 2)) == False
    assert (RootOf(x**3 + 2, 2) == RootOf(x**3 + 2, 2)) == True

def test_RootOf_is_real():
    assert RootOf(x**3 + 2, 0).is_real == True
    assert RootOf(x**3 + 2, 1).is_real == False
    assert RootOf(x**3 + 2, 2).is_real == False

def test_RootOf_is_complex():
    assert RootOf(x**3 + 2, 0).is_complex == False
    assert RootOf(x**3 + 2, 1).is_complex == True
    assert RootOf(x**3 + 2, 2).is_complex == True

def test_RootOf_evalf():
    real = RootOf(x**3 + 2, 0).evalf(n=20)

    assert real.epsilon_eq(Real("-1.2599210498948731648"))

    re, im = RootOf(x**3 + 2, 1).evalf(n=20).as_real_imag()

    assert re.epsilon_eq(Real("0.62996052494743658238"))
    assert im.epsilon_eq(Real("1.09112363597172140360"))

    re, im = RootOf(x**3 + 2, 2).evalf(n=20).as_real_imag()

    assert re.epsilon_eq( Real("0.62996052494743658238"))
    assert im.epsilon_eq(-Real("1.09112363597172140360"))
