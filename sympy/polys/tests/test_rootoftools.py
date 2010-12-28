"""Tests for the implementation of RootOf class and related tools. """

from sympy.polys.polytools import Poly

from sympy.polys.rootoftools import (
    RootOf, RootSum, _rootof_preprocess,
)

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
    PolynomialError,
    DomainError,
)

from sympy import (
    symbols, sqrt, I, Rational, Real, Lambda,
)

from sympy.utilities.pytest import raises
from sympy.abc import x, y

def test_RootOf___new__():
    assert RootOf(x, 0) == 0
    assert RootOf(x,-1) == 0

    assert RootOf(x - 1, 0) == 1
    assert RootOf(x - 1,-1) == 1

    assert RootOf(x + 1, 0) ==-1
    assert RootOf(x + 1,-1) ==-1

    assert RootOf(x**2 + 2*x + 3, 0) == -1 + I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3, 1) == -1 - I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3,-1) == -1 - I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3,-2) == -1 + I*sqrt(2)

    r = RootOf(x**2 + 2*x + 3, 0, radicals=False)
    assert isinstance(r, RootOf) == True

    r = RootOf(x**2 + 2*x + 3, 1, radicals=False)
    assert isinstance(r, RootOf) == True

    r = RootOf(x**2 + 2*x + 3,-1, radicals=False)
    assert isinstance(r, RootOf) == True

    r = RootOf(x**2 + 2*x + 3,-2, radicals=False)
    assert isinstance(r, RootOf) == True

    assert RootOf((x - 1)*(x + 1), 0, radicals=False) ==-1
    assert RootOf((x - 1)*(x + 1), 1, radicals=False) == 1
    assert RootOf((x - 1)*(x + 1),-1, radicals=False) == 1
    assert RootOf((x - 1)*(x + 1),-2, radicals=False) ==-1

    assert RootOf((x - 1)*(x + 1), 0, radicals=True) ==-1
    assert RootOf((x - 1)*(x + 1), 1, radicals=True) == 1
    assert RootOf((x - 1)*(x + 1),-1, radicals=True) == 1
    assert RootOf((x - 1)*(x + 1),-2, radicals=True) ==-1

    assert RootOf((x - 1)*(x**3 + x + 3), 0) == RootOf(x**3 + x + 3, 0)
    assert RootOf((x - 1)*(x**3 + x + 3), 1) == 1
    assert RootOf((x - 1)*(x**3 + x + 3), 2) == RootOf(x**3 + x + 3, 1)
    assert RootOf((x - 1)*(x**3 + x + 3), 3) == RootOf(x**3 + x + 3, 2)
    assert RootOf((x - 1)*(x**3 + x + 3),-1) == RootOf(x**3 + x + 3, 2)
    assert RootOf((x - 1)*(x**3 + x + 3),-2) == RootOf(x**3 + x + 3, 1)
    assert RootOf((x - 1)*(x**3 + x + 3),-3) == 1
    assert RootOf((x - 1)*(x**3 + x + 3),-4) == RootOf(x**3 + x + 3, 0)

    assert RootOf(x**4 + 3*x**3, 0) ==-3
    assert RootOf(x**4 + 3*x**3, 1) == 0
    assert RootOf(x**4 + 3*x**3, 2) == 0
    assert RootOf(x**4 + 3*x**3, 3) == 0

    raises(GeneratorsNeeded, "RootOf(0, 0)")
    raises(GeneratorsNeeded, "RootOf(1, 0)")

    raises(PolynomialError, "RootOf(Poly(0, x), 0)")
    raises(PolynomialError, "RootOf(Poly(1, x), 0)")

    raises(PolynomialError, "RootOf(x - y, 0)")

    raises(DomainError, "RootOf(x**3 - x + sqrt(2), 0)")
    raises(DomainError, "RootOf(x**3 - x + I, 0)")

    raises(IndexError, "RootOf(x**2 - 1,-4)")
    raises(IndexError, "RootOf(x**2 - 1,-3)")
    raises(IndexError, "RootOf(x**2 - 1, 2)")
    raises(IndexError, "RootOf(x**2 - 1, 3)")

    assert RootOf(Poly(x - y, x), 0) == y

    assert RootOf(Poly(x**2 - y, x), 0) == +sqrt(y)
    assert RootOf(Poly(x**2 - y, x), 1) == -sqrt(y)

    assert RootOf(Poly(x**3 - y, x), 0) == y**Rational(1,3)

    raises(DomainError, "RootOf(Poly(x**3 + x - y, x), 0)")

def test_RootOf___new___indices():
    r0 = RootOf(x**3 + x + 3, 0)
    r1 = RootOf(x**3 + x + 3, 1)
    r2 = RootOf(x**3 + x + 3, 2)

    assert RootOf(x**3 + x + 3) == [r0]

    assert RootOf(x**3 + x + 3, (0,)) == [r0]
    assert RootOf(x**3 + x + 3, (0,1)) == [r0, r1]
    assert RootOf(x**3 + x + 3, (0,1,2)) == [r0, r1, r2]

    assert RootOf(x**3 + x + 3, (-3,)) == [r0]
    assert RootOf(x**3 + x + 3, (-3,-2)) == [r0, r1]
    assert RootOf(x**3 + x + 3, (-3,-2,-1)) == [r0, r1, r2]

def test_RootOf___eq__():
    assert (RootOf(x**3 + x + 3, 0) == RootOf(x**3 + x + 3, 0)) == True
    assert (RootOf(x**3 + x + 3, 0) == RootOf(x**3 + x + 3, 1)) == False
    assert (RootOf(x**3 + x + 3, 1) == RootOf(x**3 + x + 3, 1)) == True
    assert (RootOf(x**3 + x + 3, 1) == RootOf(x**3 + x + 3, 2)) == False
    assert (RootOf(x**3 + x + 3, 2) == RootOf(x**3 + x + 3, 2)) == True

def test_RootOf_is_real():
    assert RootOf(x**3 + x + 3, 0).is_real == True
    assert RootOf(x**3 + x + 3, 1).is_real == False
    assert RootOf(x**3 + x + 3, 2).is_real == False

def test_RootOf_is_complex():
    assert RootOf(x**3 + x + 3, 0).is_complex == False
    assert RootOf(x**3 + x + 3, 1).is_complex == True
    assert RootOf(x**3 + x + 3, 2).is_complex == True

def test_RootOf_evalf():
    real = RootOf(x**3 + x + 3, 0).evalf(n=20)

    assert real.epsilon_eq(Real("-1.2134116627622296341"))

    re, im = RootOf(x**3 + x + 3, 1).evalf(n=20).as_real_imag()

    assert re.epsilon_eq(Real("0.60670583138111481707"))
    assert im.epsilon_eq(Real("1.45061224918844152650"))

    re, im = RootOf(x**3 + x + 3, 2).evalf(n=20).as_real_imag()

    assert re.epsilon_eq( Real("0.60670583138111481707"))
    assert im.epsilon_eq(-Real("1.45061224918844152650"))

def test_RootOf_preprocessing():
    E, F, J, L = symbols("E,F,J,L")

    f = -21601054687500000000*E**8*J**8/L**16 + \
         508232812500000000*F*x*E**7*J**7/L**14 - \
         4269543750000000*E**6*F**2*J**6*x**2/L**12 + \
         16194716250000*E**5*F**3*J**5*x**3/L**10 - \
         27633173750*E**4*F**4*J**4*x**4/L**8 + \
         14840215*E**3*F**5*J**3*x**5/L**6 + \
         54794*E**2*F**6*J**2*x**6/(5*L**4) - \
         1153*E*J*F**7*x**7/(80*L**2) + \
         633*F**8*x**8/160000

    coeff, poly = _rootof_preprocess(Poly(f, x))

    assert coeff == 20*E*J/(F*L**2)
    assert poly == 633*x**8 - 115300*x**7 + 4383520*x**6 + 296804300*x**5 - 27633173750*x**4 + \
        809735812500*x**3 - 10673859375000*x**2 + 63529101562500*x - 135006591796875

def test_RootSum___new__():
    f = x**3 + x + 3

    assert isinstance(RootSum(f), RootSum) == True
    assert RootSum(f).doit() == RootOf(f, 0) + RootOf(f, 1) + RootOf(f, 2)

    assert RootSum(f**2) == 2*RootSum(f)
    assert RootSum(f**2).doit() == 2*(RootOf(f, 0) + RootOf(f, 1) + RootOf(f, 2))

    assert RootSum((x - 7)*f**3) == 7 + 3*RootSum(f)
    assert RootSum((x - 7)*f**3).doit() == 7 + 3*(RootOf(f, 0) + RootOf(f, 1) + RootOf(f, 2))

    g = Lambda(x, x**2)

    assert isinstance(RootSum(f, g), RootSum) == True
    assert RootSum(f, g).doit() == RootOf(f, 0)**2 + RootOf(f, 1)**2 + RootOf(f, 2)**2

    assert RootSum(f**2, g) == 2*RootSum(f, Lambda(x, x**2))
    assert RootSum(f**2, g).doit() == 2*(RootOf(f, 0)**2 + RootOf(f, 1)**2 + RootOf(f, 2)**2)

    assert RootSum((x - 7)*f**3, g) == 49 + 3*RootSum(f, Lambda(x, x**2))
    assert RootSum((x - 7)*f**3, g).doit() == 49 + 3*(RootOf(f, 0)**2 + RootOf(f, 1)**2 + RootOf(f, 2)**2)

    assert RootSum(f, formal=False) == RootOf(f, 0) + RootOf(f, 1) + RootOf(f, 2)

    raises(PolynomialError, "RootSum(x**3 + x + y)")
    raises(TypeError, "RootSum(x**2 + 3, 10)")
