from sympy import (
    Abs, Dummy, Eq, Gt,
    LambertW, Piecewise, Poly, Rational, S, Symbol,
    acos, atan, atanh, cos, erf, erfinv, erfc, erfcinv,
    exp, log, pi, sin, sinh, sqrt, symbols,
    tan, tanh, atan2, arg,
    Lambda, imageset, cot, acot, I, EmptySet, Union, E, Interval, oo)

from sympy.core.function import nfloat
from sympy.functions.elementary.complexes import im, re
from sympy.functions.elementary.hyperbolic import HyperbolicFunction
from sympy.functions.elementary.trigonometric import TrigonometricFunction

from sympy.polys.rootoftools import CRootOf

from sympy.sets import FiniteSet

from sympy.utilities.pytest import XFAIL, raises, skip
from sympy.utilities.randtest import verify_numerically as tn
from sympy.physics.units import cm


from sympy.solvers.solveset import (
    solveset_real, domain_check, solveset_complex,
    _is_function_class_equation, invert_real, invert_complex, solveset)

a = Symbol('a', real=True)
b = Symbol('b', real=True)
c = Symbol('c', real=True)
x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
q = Symbol('q', real=True)
m = Symbol('m', real=True)
n = Symbol('n', real=True)


def test_invert_real():
    x = Symbol('x', real=True)
    x = Dummy(real=True)
    n = Symbol('n')
    d = Dummy()
    assert solveset(abs(x) - n, x) == solveset(abs(x) - d, x) == EmptySet()

    n = Symbol('n', real=True)
    assert invert_real(x + 3, y, x) == (x, FiniteSet(y - 3))
    assert invert_real(x*3, y, x) == (x, FiniteSet(y / 3))

    assert invert_real(exp(x), y, x) == (x, FiniteSet(log(y)))
    assert invert_real(exp(3*x), y, x) == (x, FiniteSet(log(y) / 3))
    assert invert_real(exp(x + 3), y, x) == (x, FiniteSet(log(y) - 3))

    assert invert_real(exp(x) + 3, y, x) == (x, FiniteSet(log(y - 3)))
    assert invert_real(exp(x)*3, y, x) == (x, FiniteSet(log(y / 3)))

    assert invert_real(log(x), y, x) == (x, FiniteSet(exp(y)))
    assert invert_real(log(3*x), y, x) == (x, FiniteSet(exp(y) / 3))
    assert invert_real(log(x + 3), y, x) == (x, FiniteSet(exp(y) - 3))

    assert invert_real(Abs(x), y, x) == (x, FiniteSet(-y, y))

    assert invert_real(2**x, y, x) == (x, FiniteSet(log(y)/log(2)))
    assert invert_real(2**exp(x), y, x) == (x, FiniteSet(log(log(y)/log(2))))

    assert invert_real(x**2, y, x) == (x, FiniteSet(sqrt(y), -sqrt(y)))
    assert invert_real(x**Rational(1, 2), y, x) == (x, FiniteSet(y**2))

    raises(ValueError, lambda: invert_real(x, x, x))
    raises(ValueError, lambda: invert_real(x**pi, y, x))
    raises(ValueError, lambda: invert_real(S.One, y, x))

    assert invert_real(x**31 + x, y, x) == (x**31 + x, FiniteSet(y))

    assert invert_real(Abs(x**31 + x + 1), y, x) == (x**31 + x,
                                                     FiniteSet(-y - 1, y - 1))

    assert invert_real(tan(x), y, x) == \
        (x, imageset(Lambda(n, n*pi + atan(y)), S.Integers))

    assert invert_real(tan(exp(x)), y, x) == \
        (x, imageset(Lambda(n, log(n*pi + atan(y))), S.Integers))

    assert invert_real(cot(x), y, x) == \
        (x, imageset(Lambda(n, n*pi + acot(y)), S.Integers))
    assert invert_real(cot(exp(x)), y, x) == \
        (x, imageset(Lambda(n, log(n*pi + acot(y))), S.Integers))

    assert invert_real(tan(tan(x)), y, x) == \
        (tan(x), imageset(Lambda(n, n*pi + atan(y)), S.Integers))

    x = Symbol('x', positive=True)
    assert invert_real(x**pi, y, x) == (x, FiniteSet(y**(1/pi)))


def test_invert_complex():
    assert invert_complex(x + 3, y, x) == (x, FiniteSet(y - 3))
    assert invert_complex(x*3, y, x) == (x, FiniteSet(y / 3))

    assert invert_complex(exp(x), y, x) == \
        (x, imageset(Lambda(n, I*(2*pi*n + arg(y)) + log(Abs(y))), S.Integers))

    assert invert_complex(log(x), y, x) == (x, FiniteSet(exp(y)))

    raises(ValueError, lambda: invert_real(S.One, y, x))
    raises(ValueError, lambda: invert_complex(x, x, x))


def test_domain_check():
    assert domain_check(1/(1 + (1/(x+1))**2), x, -1) is False
    assert domain_check(x**2, x, 0) is True
    assert domain_check(x, x, oo) is False
    assert domain_check(0, x, oo) is False


def test_is_function_class_equation():
    from sympy.abc import x, a
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x), x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x) - 1, x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x) + sin(x), x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x) + sin(x) - a, x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       sin(x)*tan(x) + sin(x), x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       sin(x)*tan(x + a) + sin(x), x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       sin(x)*tan(x*a) + sin(x), x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       a*tan(x) - 1, x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x)**2 + sin(x) - 1, x) is True
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x**2), x) is False
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x**2) + sin(x), x) is False
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(x)**sin(x), x) is False
    assert _is_function_class_equation(TrigonometricFunction,
                                       tan(sin(x)) + sin(x), x) is False
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x), x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x) - 1, x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x) + sinh(x), x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x) + sinh(x) - a, x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       sinh(x)*tanh(x) + sinh(x), x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       sinh(x)*tanh(x + a) + sinh(x), x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       sinh(x)*tanh(x*a) + sinh(x), x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       a*tanh(x) - 1, x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x)**2 + sinh(x) - 1, x) is True
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x**2), x) is False
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x**2) + sinh(x), x) is False
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(x)**sinh(x), x) is False
    assert _is_function_class_equation(HyperbolicFunction,
                                       tanh(sinh(x)) + sinh(x), x) is False


def test_garbage_input():
    raises(ValueError, lambda: solveset_real([x], x))
    raises(ValueError, lambda: solveset_real(x, pi))
    raises(ValueError, lambda: solveset_real(x, x**2))

    raises(ValueError, lambda: solveset_complex([x], x))
    raises(ValueError, lambda: solveset_complex(x, pi))


def test_solve_mul():
    assert solveset_real((a*x + b)*(exp(x) - 3), x) == \
        FiniteSet(-b/a, log(3))
    assert solveset_real((2*x + 8)*(8 + exp(x)), x) == FiniteSet(S(-4))
    assert solveset_real(x/log(x), x) == EmptySet()


def test_solve_invert():
    assert solveset_real(exp(x) - 3, x) == FiniteSet(log(3))
    assert solveset_real(log(x) - 3, x) == FiniteSet(exp(3))

    assert solveset_real(3**(x + 2), x) == FiniteSet()
    assert solveset_real(3**(2 - x), x) == FiniteSet()

    b = Symbol('b', positive=True)
    y = Symbol('y', positive=True)
    assert solveset_real(y - b*exp(a/x), x) == FiniteSet(a/log(y/b))
    # issue 4504
    assert solveset_real(2**x - 10, x) == FiniteSet(log(10)/log(2))


def test_errorinverses():
    assert solveset_real(erf(x) - S.One/2, x) == \
        FiniteSet(erfinv(S.One/2))
    assert solveset_real(erfinv(x) - 2, x) == \
        FiniteSet(erf(2))
    assert solveset_real(erfc(x) - S.One, x) == \
        FiniteSet(erfcinv(S.One))
    assert solveset_real(erfcinv(x) - 2, x) == FiniteSet(erfc(2))


def test_solve_polynomial():
    assert solveset_real(3*x - 2, x) == FiniteSet(Rational(2, 3))

    assert solveset_real(x**2 - 1, x) == FiniteSet(-S(1), S(1))
    assert solveset_real(x - y**3, x) == FiniteSet(y ** 3)

    a11, a12, a21, a22, b1, b2 = symbols('a11, a12, a21, a22, b1, b2')

    assert solveset_real(x**3 - 15*x - 4, x) == FiniteSet(
        -2 + 3 ** Rational(1, 2),
        S(4),
        -2 - 3 ** Rational(1, 2))

    assert solveset_real(sqrt(x) - 1, x) == FiniteSet(1)
    assert solveset_real(sqrt(x) - 2, x) == FiniteSet(4)
    assert solveset_real(x**Rational(1, 4) - 2, x) == FiniteSet(16)
    assert solveset_real(x**Rational(1, 3) - 3, x) == FiniteSet(27)
    assert len(solveset_real(x**5 + x**3 + 1, x)) == 1
    assert len(solveset_real(-2*x**3 + 4*x**2 - 2*x + 6, x)) > 0


def test_return_root_of():
    f = x**5 - 15*x**3 - 5*x**2 + 10*x + 20
    s = list(solveset_complex(f, x))
    for root in s:
        assert root.func == CRootOf

    # if one uses solve to get the roots of a polynomial that has a CRootOf
    # solution, make sure that the use of nfloat during the solve process
    # doesn't fail. Note: if you want numerical solutions to a polynomial
    # it is *much* faster to use nroots to get them than to solve the
    # equation only to get CRootOf solutions which are then numerically
    # evaluated. So for eq = x**5 + 3*x + 7 do Poly(eq).nroots() rather
    # than [i.n() for i in solve(eq)] to get the numerical roots of eq.
    assert nfloat(list(solveset_complex(x**5 + 3*x**3 + 7, x))[0],
                  exponent=False) == CRootOf(x**5 + 3*x**3 + 7, 0).n()

    sol = list(solveset_complex(x**6 - 2*x + 2, x))
    assert all(isinstance(i, CRootOf) for i in sol) and len(sol) == 6

    f = x**5 - 15*x**3 - 5*x**2 + 10*x + 20
    s = list(solveset_complex(f, x))
    for root in s:
        assert root.func == CRootOf

    s = x**5 + 4*x**3 + 3*x**2 + S(7)/4
    assert solveset_complex(s, x) == \
        FiniteSet(*Poly(s*4, domain='ZZ').all_roots())

    # XXX: this comparison should work without converting the FiniteSet to list
    # See #7876
    eq = x*(x - 1)**2*(x + 1)*(x**6 - x + 1)
    assert list(solveset_complex(eq, x)) == \
        list(FiniteSet(-1, 0, 1, CRootOf(x**6 - x + 1, 0),
                       CRootOf(x**6 - x + 1, 1),
                       CRootOf(x**6 - x + 1, 2),
                       CRootOf(x**6 - x + 1, 3),
                       CRootOf(x**6 - x + 1, 4),
                       CRootOf(x**6 - x + 1, 5)))


def test__has_rational_power():
    from sympy.solvers.solveset import _has_rational_power
    assert _has_rational_power(sqrt(2), x)[0] is False
    assert _has_rational_power(x*sqrt(2), x)[0] is False

    assert _has_rational_power(x**2*sqrt(x), x) == (True, 2)
    assert _has_rational_power(sqrt(2)*x**(S(1)/3), x) == (True, 3)
    assert _has_rational_power(sqrt(x)*x**(S(1)/3), x) == (True, 6)


def test_solveset_sqrt_1():
    assert solveset_real(sqrt(5*x + 6) - 2 - x, x) == \
        FiniteSet(-S(1), S(2))
    assert solveset_real(sqrt(x - 1) - x + 7, x) == FiniteSet(10)
    assert solveset_real(sqrt(x - 2) - 5, x) == FiniteSet(27)
    assert solveset_real(sqrt(x) - 2 - 5, x) == FiniteSet(49)
    assert solveset_real(sqrt(x**3), x) == FiniteSet(0)
    assert solveset_real(sqrt(x - 1), x) == FiniteSet(1)


def test_solveset_sqrt_2():
    # http://tutorial.math.lamar.edu/Classes/Alg/SolveRadicalEqns.aspx#Solve_Rad_Ex2_a
    assert solveset_real(sqrt(2*x - 1) - sqrt(x - 4) - 2, x) == \
        FiniteSet(S(5), S(13))
    assert solveset_real(sqrt(x + 7) + 2 - sqrt(3 - x), x) == \
        FiniteSet(-6)

    # http://www.purplemath.com/modules/solverad.htm
    assert solveset_real(sqrt(17*x - sqrt(x**2 - 5)) - 7, x) == \
        FiniteSet(3)

    eq = x + 1 - (x**4 + 4*x**3 - x)**Rational(1, 4)
    assert solveset_real(eq, x) == FiniteSet(-S(1)/2, -S(1)/3)

    eq = sqrt(2*x + 9) - sqrt(x + 1) - sqrt(x + 4)
    assert solveset_real(eq, x) == FiniteSet(0)

    eq = sqrt(x + 4) + sqrt(2*x - 1) - 3*sqrt(x - 1)
    assert solveset_real(eq, x) == FiniteSet(5)

    eq = sqrt(x)*sqrt(x - 7) - 12
    assert solveset_real(eq, x) == FiniteSet(16)

    eq = sqrt(x - 3) + sqrt(x) - 3
    assert solveset_real(eq, x) == FiniteSet(4)

    eq = sqrt(2*x**2 - 7) - (3 - x)
    assert solveset_real(eq, x) == FiniteSet(-S(8), S(2))

    # others
    eq = sqrt(9*x**2 + 4) - (3*x + 2)
    assert solveset_real(eq, x) == FiniteSet(0)

    assert solveset_real(sqrt(x - 3) - sqrt(x) - 3, x) == FiniteSet()

    eq = (2*x - 5)**Rational(1, 3) - 3
    assert solveset_real(eq, x) == FiniteSet(16)

    assert solveset_real(sqrt(x) + sqrt(sqrt(x)) - 4, x) == \
        FiniteSet((-S.Half + sqrt(17)/2)**4)

    eq = sqrt(x) - sqrt(x - 1) + sqrt(sqrt(x))
    assert solveset_real(eq, x) == FiniteSet()

    eq = (sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5)
    ans = solveset_real(eq, x)
    ra = S('''-1484/375 - 4*(-1/2 + sqrt(3)*I/2)*(-12459439/52734375 +
    114*sqrt(12657)/78125)**(1/3) - 172564/(140625*(-1/2 +
    sqrt(3)*I/2)*(-12459439/52734375 + 114*sqrt(12657)/78125)**(1/3))''')
    rb = S(4)/5
    assert all(abs(eq.subs(x, i).n()) < 1e-10 for i in (ra, rb)) and \
        len(ans) == 2 and \
        set([i.n(chop=True) for i in ans]) == \
        set([i.n(chop=True) for i in (ra, rb)])

    assert solveset_real(sqrt(x) + x**Rational(1, 3) +
                                 x**Rational(1, 4), x) == FiniteSet(0)

    assert solveset_real(x/sqrt(x**2 + 1), x) == FiniteSet(0)

    eq = (x - y**3)/((y**2)*sqrt(1 - y**2))
    assert solveset_real(eq, x) == FiniteSet(y**3)

    # issue 4497
    assert solveset_real(1/(5 + x)**(S(1)/5) - 9, x) == \
        FiniteSet(-295244/S(59049))


@XFAIL
def test_solve_sqrt_fail():
    # this only works if we check real_root(eq.subs(x, S(1)/3))
    # but checksol doesn't work like that
    eq = (x**3 - 3*x**2)**Rational(1, 3) + 1 - x
    assert solveset_real(eq, x) == FiniteSet(S(1)/3)


def test_solve_sqrt_3():
    R = Symbol('R')
    eq = sqrt(2)*R*sqrt(1/(R + 1)) + (R + 1)*(sqrt(2)*sqrt(1/(R + 1)) - 1)
    sol = solveset_complex(eq, R)

    assert sol == FiniteSet(*[S(5)/3 + 4*sqrt(10)*cos(atan(3*sqrt(111)/251)/3)/3,
        -sqrt(10)*cos(atan(3*sqrt(111)/251)/3)/3 + 40*re(1/((-S(1)/2 -
        sqrt(3)*I/2)*(S(251)/27 + sqrt(111)*I/9)**(S(1)/3)))/9 +
        sqrt(30)*sin(atan(3*sqrt(111)/251)/3)/3 + S(5)/3 +
        I*(-sqrt(30)*cos(atan(3*sqrt(111)/251)/3)/3 -
        sqrt(10)*sin(atan(3*sqrt(111)/251)/3)/3 + 40*im(1/((-S(1)/2 -
        sqrt(3)*I/2)*(S(251)/27 + sqrt(111)*I/9)**(S(1)/3)))/9)])

    # the number of real roots will depend on the value of m: for m=1 there are 4
    # and for m=-1 there are none.
    eq = -sqrt((m - q)**2 + (-m/(2*q) + S(1)/2)**2) + sqrt((-m**2/2 - sqrt(
        4*m**4 - 4*m**2 + 8*m + 1)/4 - S(1)/4)**2 + (m**2/2 - m - sqrt(
            4*m**4 - 4*m**2 + 8*m + 1)/4 - S(1)/4)**2)
    raises(NotImplementedError, lambda: solveset_real(eq, q))


def test_solve_polynomial_symbolic_param():
    assert solveset_complex((x**2 - 1)**2 - a, x) == \
        FiniteSet(sqrt(1 + sqrt(a)), -sqrt(1 + sqrt(a)),
                  sqrt(1 - sqrt(a)), -sqrt(1 - sqrt(a)))

    # By attempt to make Set.contains behave symbolically SetDifference on
    # FiniteSet isn't working very well.
    # Simple operations like `FiniteSet(a) - FiniteSet(-b)` raises `TypeError`
    # The likely course of action will making such operations return
    # SetDifference object. That will also change the expected output of
    # the given tests. Till the SetDifference becomes well behaving again the
    # following tests are kept as comments.

    # # issue 4508
    # assert solveset_complex(y - b*x/(a + x), x) == \
    #     FiniteSet(-a*y/(y - b))
    #
    # # issue 4507
    # assert solveset_complex(y - b/(1 + a*x), x) == \
    #     FiniteSet((b - y)/(a*y))


def test_solve_rational():
    assert solveset_real(1/x + 1, x) == FiniteSet(-S.One)
    assert solveset_real(1/exp(x) - 1, x) == FiniteSet(0)
    assert solveset_real(x*(1 - 5/x), x) == FiniteSet(5)
    assert solveset_real(2*x/(x + 2) - 1, x) == FiniteSet(2)
    assert solveset_real((x**2/(7 - x)).diff(x), x) == \
        FiniteSet(S(0), S(14))


def test_solveset_real_gen_is_pow():
    assert solveset_real(sqrt(1) + 1, x) == EmptySet()


def test_no_sol():
    assert solveset_real(4, x) == EmptySet()
    assert solveset_real(exp(x), x) == EmptySet()
    assert solveset_real(x**2 + 1, x) == EmptySet()
    assert solveset_real(-3*a/sqrt(x), x) == EmptySet()
    assert solveset_real(1/x, x) == EmptySet()
    assert solveset_real(-(1 + x)/(2 + x)**2 + 1/(2 + x), x) == \
        EmptySet()


def test_sol_zero_real():
    assert solveset_real(0, x) == S.Reals
    assert solveset_real(-x**2 - 2*x + (x + 1)**2 - 1, x) == S.Reals


def test_no_sol_rational_extragenous():
    assert solveset_real((x/(x + 1) + 3)**(-2), x) == EmptySet()
    assert solveset_real((x - 1)/(1 + 1/(x - 1)), x) == EmptySet()


def test_solve_polynomial_cv_1a():
    """
    Test for solving on equations that can be converted to
    a polynomial equation using the change of variable y -> x**Rational(p, q)
    """
    assert solveset_real(sqrt(x) - 1, x) == FiniteSet(1)
    assert solveset_real(sqrt(x) - 2, x) == FiniteSet(4)
    assert solveset_real(x**Rational(1, 4) - 2, x) == FiniteSet(16)
    assert solveset_real(x**Rational(1, 3) - 3, x) == FiniteSet(27)
    assert solveset_real(x*(x**(S(1) / 3) - 3), x) == \
        FiniteSet(S(0), S(27))


def test_solveset_real_rational():
    """Test solveset_real for rational functions"""
    assert solveset_real((x - y**3) / ((y**2)*sqrt(1 - y**2)), x) \
        == FiniteSet(y**3)
    # issue 4486
    assert solveset_real(2*x/(x + 2) - 1, x) == FiniteSet(2)


def test_solveset_real_log():
    assert solveset_real(log((x-1)*(x+1)), x) == \
        FiniteSet(sqrt(2), -sqrt(2))


def test_poly_gens():
    assert solveset_real(4**(2*(x**2) + 2*x) - 8, x) == \
        FiniteSet(-Rational(3, 2), S.Half)


@XFAIL
def test_uselogcombine_1():
    assert solveset_real(log(x - 3) + log(x + 3), x) == \
        FiniteSet(sqrt(10))
    assert solveset_real(log(x + 1) - log(2*x - 1), x) == FiniteSet(2)
    assert solveset_real(log(x + 3) + log(1 + 3/x) - 3) == FiniteSet(
        -3 + sqrt(-12 + exp(3))*exp(S(3)/2)/2 + exp(3)/2,
        -sqrt(-12 + exp(3))*exp(S(3)/2)/2 - 3 + exp(3)/2)


@XFAIL
def test_uselogcombine_2():
    eq = z - log(x) + log(y/(x*(-1 + y**2/x**2)))
    assert solveset_real(eq, x) == \
        FiniteSet(-sqrt(y*(y - exp(z))), sqrt(y*(y - exp(z))))


def test_solve_abs():
    assert solveset_real(Abs(x) - 2, x) == FiniteSet(-2, 2)
    assert solveset_real(Abs(x + 3) - 2*Abs(x - 3), x) == \
        FiniteSet(1, 9)
    assert solveset_real(2*Abs(x) - Abs(x - 1), x) == \
        FiniteSet(-1, Rational(1, 3))

    assert solveset_real(Abs(x - 7) - 8, x) == FiniteSet(-S(1), S(15))


@XFAIL
def test_rewrite_trigh():
    # if this import passes then the test below should also pass
    from sympy import sech
    assert solveset_real(sinh(x) + sech(x), x) == FiniteSet(
        2*atanh(-S.Half + sqrt(5)/2 - sqrt(-2*sqrt(5) + 2)/2),
        2*atanh(-S.Half + sqrt(5)/2 + sqrt(-2*sqrt(5) + 2)/2),
        2*atanh(-sqrt(5)/2 - S.Half + sqrt(2 + 2*sqrt(5))/2),
        2*atanh(-sqrt(2 + 2*sqrt(5))/2 - sqrt(5)/2 - S.Half))


@XFAIL
def test_real_imag_splitting1():
    a, b = symbols('a b', real=True, finite=True)
    s = solveset_real(sqrt(a**2 + b**2) - 3, a)
    assert s != S.EmptySet
    # FiniteSet(-sqrt(-b**2 + 9), sqrt(-b**2 + 9))
    # fails now because whether it is real or not depends
    # on the value of b, e.g. b = 4 gives an imaginary value


def test_real_imag_splitting():
    a, b = symbols('a b', real=True, finite=True)
    assert solveset_real(sqrt(a**2 - b**2) - 3, a) == \
        FiniteSet(-sqrt(b**2 + 9), sqrt(b**2 + 9))


def test_units():
    assert solveset_real(1/x - 1/(2*cm), x) == FiniteSet(2*cm)


def test_solve_only_exp_1():
    y = Symbol('y', positive=True, finite=True)
    assert solveset_real(exp(x) - y, x) == FiniteSet(log(y))


@XFAIL
def test_only_exp_2():
    assert solveset_real(exp(x) + exp(-x) - 4, x) == \
        FiniteSet(log(-sqrt(3) + 2), log(sqrt(3) + 2))
    assert solveset_real(exp(x) + exp(-x) - y, x) != S.EmptySet
    # FiniteSet(log(y/2 - sqrt((y - 2)*(y + 2))/2),
    #           log(y/2 + sqrt((y - 2)*(y + 2))/2))
    # fails now because whether it is real or not depends
    # on whether y >= 2


@XFAIL
def test_only_exp_3():
    assert solveset_real(exp(x/y)*exp(-z/y) - 2, y) == \
        FiniteSet((x - z)/log(2))
    assert solveset_real(sqrt(exp(x)) + sqrt(exp(-x)) - 4, x) == \
        FiniteSet(2*log(-sqrt(3) + 2), 2*log(sqrt(3) + 2))


def test_atan2():
    # The .inverse() method on atan2 works only if x.is_real is True and the
    # second argument is a real constant
    assert solveset_real(atan2(x, 2) - pi/3, x) == FiniteSet(2*sqrt(3))


def test_piecewise():
    eq = Piecewise((x - 2, Gt(x, 2)), (2 - x, True)) - 3
    assert set(solveset_real(eq, x)) == set(FiniteSet(-1, 5))
    absxm3 = Piecewise(
        (x - 3, S(0) <= x - 3),
        (3 - x, S(0) > x - 3)
    )
    y = Symbol('y', positive=True)
    assert solveset_real(absxm3 - y, x) == FiniteSet(-y + 3, y + 3)


def test_solveset_complex_polynomial():
    from sympy.abc import x, a, b, c
    assert solveset_complex(a*x**2 + b*x + c, x) == \
        FiniteSet(-b/(2*a) - sqrt(-4*a*c + b**2)/(2*a),
                  -b/(2*a) + sqrt(-4*a*c + b**2)/(2*a))

    assert solveset_complex(x - y**3, y) == FiniteSet(
        (-x**Rational(1, 3))/2 + I*sqrt(3)*x**Rational(1, 3)/2,
        x**Rational(1, 3),
        (-x**Rational(1, 3))/2 - I*sqrt(3)*x**Rational(1, 3)/2)

    assert solveset_complex(x + 1/x - 1, x) == \
        FiniteSet(Rational(1, 2) + I*sqrt(3)/2, Rational(1, 2) - I*sqrt(3)/2)


def test_sol_zero_complex():
    # This should return the complex set after it is implemented
    raises(NotImplementedError, lambda: solveset_complex(0, x))


def test_solveset_complex_rational():
    assert solveset_complex((x - 1)*(x - I)/(x - 3), x) == \
        FiniteSet(1, I)

    assert solveset_complex((x - y**3)/((y**2)*sqrt(1 - y**2)), x) == \
        FiniteSet(y**3)
    assert solveset_complex(-x**2 - I, x) == \
        FiniteSet(-sqrt(2)/2 + sqrt(2)*I/2, sqrt(2)/2 - sqrt(2)*I/2)


def test_solve_quintics():
    skip("This test is too slow")
    f = x**5 - 110*x**3 - 55*x**2 + 2310*x + 979
    s = solveset_complex(f, x)
    for root in s:
        res = f.subs(x, root.n()).n()
        assert tn(res, 0)

    f = x**5 + 15*x + 12
    s = solveset_complex(f, x)
    for root in s:
        res = f.subs(x, root.n()).n()
        assert tn(res, 0)


def test_solveset_complex_exp():
    from sympy.abc import x, n
    assert solveset_complex(exp(x) - 1, x) == \
        imageset(Lambda(n, I*2*n*pi), S.Integers)
    assert solveset_complex(exp(x) - I, x) == \
        imageset(Lambda(n, I*(2*n*pi + pi/2)), S.Integers)


def test_solve_complex_log():
    assert solveset_complex(log(x), x) == FiniteSet(1)
    assert solveset_complex(1 - log(a + 4*x**2), x) == \
        FiniteSet(-sqrt(-a/4 + E/4), sqrt(-a/4 + E/4))


def test_solve_complex_sqrt():
    assert solveset_complex(sqrt(5*x + 6) - 2 - x, x) == \
        FiniteSet(-S(1), S(2))
    assert solveset_complex(sqrt(5*x + 6) - (2 + 2*I) - x, x) == \
        FiniteSet(-S(2), 3 - 4*I)
    assert solveset_complex(4*x*(1 - a * sqrt(x)), x) == \
        FiniteSet(S(0), 1 / a ** 2)


def test_solveset_complex_tan():
    s = solveset_complex(tan(x).rewrite(exp), x)
    assert s == imageset(Lambda(n, pi*n), S.Integers) - \
        imageset(Lambda(n, pi*n + pi/2), S.Integers)


def test_solve_trig():
    from sympy.abc import n
    assert solveset_real(sin(x), x) == \
        Union(imageset(Lambda(n, 2*pi*n), S.Integers),
              imageset(Lambda(n, 2*pi*n + pi), S.Integers))

    assert solveset_real(sin(x) - 1, x) == \
        imageset(Lambda(n, 2*pi*n + pi/2), S.Integers)

    assert solveset_real(cos(x), x) == \
        Union(imageset(Lambda(n, 2*pi*n - pi/2), S.Integers),
              imageset(Lambda(n, 2*pi*n + pi/2), S.Integers))

    assert solveset_real(sin(x) + cos(x), x) == \
        Union(imageset(Lambda(n, 2*n*pi - pi/4), S.Integers),
              imageset(Lambda(n, 2*n*pi + 3*pi/4), S.Integers))

    assert solveset_real(sin(x)**2 + cos(x)**2, x) == S.EmptySet


def test_solve_invalid_sol():
    assert 0 not in solveset_real(sin(x)/x, x)
    assert 0 not in solveset_complex((exp(x) - 1)/x, x)


def test_solve_complex_unsolvable():
    raises(NotImplementedError, lambda: solveset_complex(cos(x) - S.Half, x))


@XFAIL
def test_solve_trig_simplified():
    from sympy.abc import n
    assert solveset_real(sin(x), x) == \
        imageset(Lambda(n, n*pi), S.Integers)

    assert solveset_real(cos(x), x) == \
        imageset(Lambda(n, n*pi + pi/2), S.Integers)

    assert solveset_real(cos(x) + sin(x), x) == \
        imageset(Lambda(n, n*pi - pi/4), S.Integers)


@XFAIL
def test_solve_lambert():
    assert solveset_real(x*exp(x) - 1, x) == FiniteSet(LambertW(1))
    assert solveset_real(x + 2**x, x) == \
        FiniteSet(-LambertW(log(2))/log(2))

    # issue 4739
    assert solveset_real(exp(log(5)*x) - 2**x, x) == FiniteSet(0)
    ans = solveset_real(3*x + 5 + 2**(-5*x + 3), x)
    assert ans == FiniteSet(-Rational(5, 3) +
                            LambertW(-10240*2**(S(1)/3)*log(2)/3)/(5*log(2)))

    eq = 2*(3*x + 4)**5 - 6*7**(3*x + 9)
    result = solveset_real(eq, x)
    ans = FiniteSet((log(2401) +
                     5*LambertW(-log(7**(7*3**Rational(1, 5)/5))))/(3*log(7))/-1)
    assert result == ans
    assert solveset_real(eq.expand(), x) == result

    assert solveset_real(5*x - 1 + 3*exp(2 - 7*x), x) == \
        FiniteSet(Rational(1, 5) + LambertW(-21*exp(Rational(3, 5))/5)/7)

    assert solveset_real(2*x + 5 + log(3*x - 2), x) == \
        FiniteSet(Rational(2, 3) + LambertW(2*exp(-Rational(19, 3))/3)/2)

    assert solveset_real(3*x + log(4*x), x) == \
        FiniteSet(LambertW(Rational(3, 4))/3)

    assert solveset_complex(x**z*y**z - 2, z) == \
        FiniteSet(log(2)/(log(x) + log(y)))

    assert solveset_real(x**x - 2) == FiniteSet(exp(LambertW(log(2))))

    a = Symbol('a')
    assert solveset_real(-a*x + 2*x*log(x), x) == FiniteSet(exp(a/2))
    a = Symbol('a', real=True)
    assert solveset_real(a/x + exp(x/2), x) == \
        FiniteSet(2*LambertW(-a/2))
    assert solveset_real((a/x + exp(x/2)).diff(x), x) == \
        FiniteSet(4*LambertW(sqrt(2)*sqrt(a)/4))

    assert solveset_real(1/(1/x - y + exp(y)), x) == EmptySet()
    # coverage test
    p = Symbol('p', positive=True)
    w = Symbol('w')
    assert solveset_real((1/p + 1)**(p + 1), p) == EmptySet()
    assert solveset_real(tanh(x + 3)*tanh(x - 3) - 1, x) == EmptySet()
    assert solveset_real(2*x**w - 4*y**w, w) == \
        solveset_real((x/y)**w - 2, w)

    assert solveset_real((x**2 - 2*x + 1).subs(x, log(x) + 3*x), x) == \
        FiniteSet(LambertW(3*S.Exp1)/3)
    assert solveset_real((x**2 - 2*x + 1).subs(x, (log(x) + 3*x)**2 - 1), x) == \
        FiniteSet(LambertW(3*exp(-sqrt(2)))/3, LambertW(3*exp(sqrt(2)))/3)
    assert solveset_real((x**2 - 2*x - 2).subs(x, log(x) + 3*x), x) == \
        FiniteSet(LambertW(3*exp(1 + sqrt(3)))/3, LambertW(3*exp(-sqrt(3) + 1))/3)
    assert solveset_real(x*log(x) + 3*x + 1, x) == \
        FiniteSet(exp(-3 + LambertW(-exp(3))))
    eq = (x*exp(x) - 3).subs(x, x*exp(x))
    assert solveset_real(eq, x) == \
        FiniteSet(LambertW(3*exp(-LambertW(3))))

    assert solveset_real(3*log(a**(3*x + 5)) + a**(3*x + 5), x) == \
        FiniteSet(-((log(a**5) + LambertW(S(1)/3))/(3*log(a))))
    p = symbols('p', positive=True)
    assert solveset_real(3*log(p**(3*x + 5)) + p**(3*x + 5), x) == \
        FiniteSet(
        log((-3**(S(1)/3) - 3**(S(5)/6)*I)*LambertW(S(1)/3)**(S(1)/3)/(2*p**(S(5)/3)))/log(p),
        log((-3**(S(1)/3) + 3**(S(5)/6)*I)*LambertW(S(1)/3)**(S(1)/3)/(2*p**(S(5)/3)))/log(p),
        log((3*LambertW(S(1)/3)/p**5)**(1/(3*log(p)))),)  # checked numerically
    # check collection
    b = Symbol('b')
    eq = 3*log(a**(3*x + 5)) + b*log(a**(3*x + 5)) + a**(3*x + 5)
    assert solveset_real(eq, x) == FiniteSet(
        -((log(a**5) + LambertW(1/(b + 3)))/(3*log(a))))

    # issue 4271
    assert solveset_real((a/x + exp(x/2)).diff(x, 2), x) == FiniteSet(
        6*LambertW((-1)**(S(1)/3)*a**(S(1)/3)/3))

    assert solveset_real(x**3 - 3**x, x) == \
        FiniteSet(-3/log(3)*LambertW(-log(3)/3))
    assert solveset_real(x**2 - 2**x, x) == FiniteSet(2)
    assert solveset_real(-x**2 + 2**x, x) == FiniteSet(2)
    assert solveset_real(3**cos(x) - cos(x)**3) == FiniteSet(
        acos(-3*LambertW(-log(3)/3)/log(3)))

    assert solveset_real(4**(x/2) - 2**(x/3), x) == FiniteSet(0)
    assert solveset_real(5**(x/2) - 2**(x/3), x) == FiniteSet(0)
    b = sqrt(6)*sqrt(log(2))/sqrt(log(5))
    assert solveset_real(5**(x/2) - 2**(3/x), x) == FiniteSet(-b, b)


def test_solveset():
    x = Symbol('x', real=True)
    raises(ValueError, lambda: solveset(x + y))

    assert solveset(exp(x) - 1) == FiniteSet(0)
    assert solveset(exp(x) - 1, x) == FiniteSet(0)
    assert solveset(Eq(exp(x), 1), x) == FiniteSet(0)

    assert solveset(x - 1 >= 0, x) == Interval(1, oo)
    assert solveset(exp(x) - 1 >= 0, x) == Interval(0, oo)

    x = Symbol('x')
    assert solveset(exp(x) - 1, x) == imageset(Lambda(n, 2*I*pi*n), S.Integers)
    assert solveset(Eq(exp(x), 1), x) == imageset(Lambda(n, 2*I*pi*n),
                                                  S.Integers)


def test_improve_coverage():
    from sympy.solvers.solveset import _has_rational_power
    x = Symbol('x', real=True)
    y = exp(x+1/x**2)
    raises(NotImplementedError, lambda: solveset(y**2+y, x))

    assert _has_rational_power(sin(x)*exp(x) + 1, x) == (False, S.One)
    assert _has_rational_power((sin(x)**2)*(exp(x) + 1)**3, x) == (False, S.One)


def test_issue_9522():
    x = Symbol('x', real=True)
    expr1 = Eq(1/(x**2 - 4) + x, 1/(x**2 - 4) + 2)
    expr2 = Eq(1/x + x, 1/x)

    assert solveset(expr1, x) == EmptySet()
    assert solveset(expr2, x) == EmptySet()
