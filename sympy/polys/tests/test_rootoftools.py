"""Tests for the implementation of RootOf class and related tools. """

from sympy.polys.polytools import Poly

from sympy.polys.rootoftools import (
    RootOf, RootSum, _rootof_preprocess,
)

from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    GeneratorsNeeded,
    PolynomialError,
    DomainError,
)

from sympy import (
    S, symbols, sqrt, I, Rational, Real, Lambda, log, exp,
)

from sympy.utilities.pytest import raises

from sympy.abc import x, y, z, r

def test_RootOf___new__():
    assert RootOf(x, 0) == 0
    assert RootOf(x,-1) == 0

    assert RootOf(x, S.Zero) == 0

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

    assert RootOf(y*x**3 + y*x + 2*y, x, 0) == -1
    raises(DomainError, "RootOf(x**3 + x + 2*y, x, 0)")

    assert RootOf(x**3 + x + 1, 0).is_commutative == True

def test_RootOf___new___indices():
    f = x**3 + x + 3
    r0 = RootOf(f, 0)
    r1 = RootOf(f, 1)
    r2 = RootOf(f, 2)

    assert RootOf(f) == [r0]

    assert RootOf(f, (0,)) == [r0]
    assert RootOf(f, (0,1)) == [r0, r1]
    assert RootOf(f, (0,1,2)) == [r0, r1, r2]

    assert RootOf(f, (-3,)) == [r0]
    assert RootOf(f, (-3,-2)) == [r0, r1]
    assert RootOf(f, (-3,-2,-1)) == [r0, r1, r2]

    assert RootOf(f, True) == [r0, r1, r2]

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

def test_RootOf_subs():
    assert RootOf(x**3 + x + 1, 0).subs(x, y) == RootOf(y**3 + y + 1, 0)

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

    g = Lambda(r, log(r*x))
    s = RootSum(f, g)

    rootofs = sum(log(RootOf(f, i)*x) for i in (0, 1, 2))

    assert isinstance(s, RootSum) == True
    assert s.doit() == rootofs

    assert RootSum(f**2, g) == 2*RootSum(f, g)
    assert RootSum(f**2, g).doit() == 2*rootofs

    assert RootSum((x - 7)*f**3, g) == log(7*x) + 3*RootSum(f, g)
    assert RootSum((x - 7)*f**3, g).doit() == log(7*x) + 3*rootofs

    raises(MultivariatePolynomialError, "RootSum(x**3 + x + y)")
    raises(ValueError, "RootSum(x**2 + 3, lambda x: x)")

    assert RootSum(f, exp) == RootSum(f, Lambda(x, exp(x)))
    assert RootSum(f, log) == RootSum(f, Lambda(x, log(x)))

    assert isinstance(RootSum(f, auto=False), RootSum) == True

    assert RootSum(f) == 0
    assert RootSum(f, Lambda(x, x)) == 0
    assert RootSum(f, Lambda(x, x**2)) == -2

    assert RootSum(f, Lambda(x, 1)) == 3
    assert RootSum(f, Lambda(x, 2)) == 6

    assert RootSum(f, auto=False).is_commutative == True

    assert RootSum(f, Lambda(x, 1/(x + x**2))) == S(11)/3
    assert RootSum(f, Lambda(x, y/(x + x**2))) == S(11)/3*y

    assert RootSum(x**2 - 1, Lambda(x, 3*x**2), x) == 6
    assert RootSum(x**2 - y, Lambda(x, 3*x**2), x) == 6*y

    assert RootSum(x**2 - 1, Lambda(x, z*x**2), x) == 2*z
    assert RootSum(x**2 - y, Lambda(x, z*x**2), x) == 2*z*y

def test_RootSum_diff():
    f = x**3 + x + 3

    g = Lambda(r,   exp(r*x))
    h = Lambda(r, r*exp(r*x))

    assert RootSum(f, g).diff(x) == RootSum(f, h)

def test_RootSum_subs():
    f = x**3 + x + 3
    g = Lambda(r, exp(r*x))

    F = y**3 + y + 3
    G = Lambda(r, exp(r*y))

    assert RootSum(f, g).subs(y, 1) == RootSum(f, g)
    assert RootSum(f, g).subs(x, y) == RootSum(F, G)

def test_RootSum_rational():
    assert RootSum(z**5 - z + 1, Lambda(z, z/(x - z))) == (4*x - 5)/(x**5 - x + 1)

    f = 161*z**3 + 115*z**2 + 19*z + 1
    g = Lambda(z, z*log(-3381*z**4/4 - 3381*z**3/4 - 625*z**2/2 - 125*z/2 - 5 + exp(x)))

    assert RootSum(f, g).diff(x) == -((5*exp(2*x) - 6*exp(x) + 4)*exp(x)/(exp(3*x) - exp(2*x) + 1))/7
