from sympy import (
    Abs, acos, Add, atan, Basic, besselsimp, binomial, collect,
    collect_const, combsimp, cos, cosh, cot, coth, count_ops,
    Derivative, diff, Dummy, E, Eq, erf, exp, exp_polar, expand,
    expand_multinomial, exptrigsimp, factor, factorial, FallingFactorial,
    Float, fraction, Function, gamma, GoldenRatio, hyper,
    hypersimp, I, Integral, integrate, log, logcombine, Matrix, MatrixSymbol,
    Mul, nsimplify, O, oo, pi, Piecewise,
    posify, powdenest, powsimp, rad, radsimp, Rational, ratsimp,
    ratsimpmodprime, rcollect, RisingFactorial, root, rootofsimp, S,
    separatevars, signsimp, simplify, sin, sinh, solve, sqrt, Subs,
    Symbol, symbols, sympify, tan, tanh, trigsimp, Wild, zoo, Sum, Lt)
from sympy.core.mul import _keep_coeff, _unevaluated_Mul as umul
from sympy.polys.polyerrors import NotInvertible
from sympy.polys.rootoftools import RootOf
from sympy.simplify.simplify import (
    collect_sqrt, fraction_expand, _unevaluated_Add, nthroot)
from sympy.utilities.pytest import XFAIL, slow, raises
from sympy.core.compatibility import range

from sympy.abc import x, y, z, t, a, b, c, d, e, f, g, h, i, k


def test_ratsimp():
    f, g = 1/x + 1/y, (x + y)/(x*y)

    assert f != g and ratsimp(f) == g

    f, g = 1/(1 + 1/x), 1 - 1/(x + 1)

    assert f != g and ratsimp(f) == g

    f, g = x/(x + y) + y/(x + y), 1

    assert f != g and ratsimp(f) == g

    f, g = -x - y - y**2/(x + y) + x**2/(x + y), -2*y

    assert f != g and ratsimp(f) == g

    f = (a*c*x*y + a*c*z - b*d*x*y - b*d*z - b*t*x*y - b*t*x - b*t*z +
         e*x)/(x*y + z)
    G = [a*c - b*d - b*t + (-b*t*x + e*x)/(x*y + z),
         a*c - b*d - b*t - ( b*t*x - e*x)/(x*y + z)]

    assert f != g and ratsimp(f) in G

    A = sqrt(pi)

    B = log(erf(x) - 1)
    C = log(erf(x) + 1)

    D = 8 - 8*erf(x)

    f = A*B/D - A*C/D + A*C*erf(x)/D - A*B*erf(x)/D + 2*A/D

    assert ratsimp(f) == A*B/8 - A*C/8 - A/(4*erf(x) - 4)


def test_ratsimpmodprime():
    a = y**5 + x + y
    b = x - y
    F = [x*y**5 - x - y]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (x**2 + x*y + x + y) / (x**2 - x*y)

    a = x + y**2 - 2
    b = x + y**2 - y - 1
    F = [x*y - 1]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (1 + y - x)/(y - x)

    a = 5*x**3 + 21*x**2 + 4*x*y + 23*x + 12*y + 15
    b = 7*x**3 - y*x**2 + 31*x**2 + 2*x*y + 15*y + 37*x + 21
    F = [x**2 + y**2 - 1]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (1 + 5*y - 5*x)/(8*y - 6*x)

    a = x*y - x - 2*y + 4
    b = x + y**2 - 2*y
    F = [x - 2, y - 3]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        Rational(2, 5)

    # Test a bug where denominators would be dropped
    assert ratsimpmodprime(x, [y - 2*x], order='lex') == \
        y/2


def test_rootofsimp():
    r = RootOf(x**5 - 5*x + 12, 0)
    s = RootOf(x**4 - x**2 + 1, 0)
    t = RootOf(x**4 - x**2 + 1, 1)

    f, g = r**5, 5*r - 12
    assert f != g and rootofsimp(f) == g

    f, g = r**10, 25*r**2 - 120*r + 144
    assert f != g and rootofsimp(f) == g

    f, g = r**4*(1 + r), r**4 + 5*r - 12
    assert f != g and rootofsimp(f) == g

    # test modular inversion
    f, g = 1/r, -r**4/12 + 5/12
    assert f != g and rootofsimp(f) == g

    f, g = -r**4/12 + 5/12, r
    assert f*g != 1 and rootofsimp(f*g) == 1

    # test multiple roots in expression
    f, g = r**5 * s**4, expand((5*r - 12)*(s**2 - 1))
    assert f != g and rootofsimp(f).expand() == g

    # test multiple roots of the same polynomial
    r0 = RootOf(x**5 - 5*x + 12, 3)
    r1 = RootOf(x**5 - 5*x + 12, 4)
    f, g = r0**5 - 5*r1 + 12, 5*r0 - 5*r1
    assert f != g and rootofsimp(f) == g

    f, g = (s**4 + 1)*(t**4 + 1), s**2*t**2
    assert f != g and rootofsimp(f) == g

    # test RootOfs with forced radicals=False
    u = RootOf(x**2 + 1, 0, radicals=False)
    v = RootOf(x**2 - 2, 0, radicals=False)

    f, g = u**2, -1
    assert f != g and rootofsimp(f) == g

    f, g = u**3, -u
    assert f != g and rootofsimp(f) == g

    f, g = v**2, 2
    assert f != g and rootofsimp(f) == g

    f, g = u**2*v**2, -2
    assert f != g and rootofsimp(f) == g

    # test no-op when expression is not polynomial the RootOf
    f, g = sin(u**2 + 1), 0
    assert f != g and rootofsimp(f) != g


def test_rootofsimp_notinvertible():
    r = RootOf(x**2 + 1, 0, radicals=False)
    expr = 1/(r**2 + 1)
    raises(NotInvertible, lambda: rootofsimp(expr))


def test_trigsimp1():
    x, y = symbols('x,y')

    assert trigsimp(1 - sin(x)**2) == cos(x)**2
    assert trigsimp(1 - cos(x)**2) == sin(x)**2
    assert trigsimp(sin(x)**2 + cos(x)**2) == 1
    assert trigsimp(1 + tan(x)**2) == 1/cos(x)**2
    assert trigsimp(1/cos(x)**2 - 1) == tan(x)**2
    assert trigsimp(1/cos(x)**2 - tan(x)**2) == 1
    assert trigsimp(1 + cot(x)**2) == 1/sin(x)**2
    assert trigsimp(1/sin(x)**2 - 1) == 1/tan(x)**2
    assert trigsimp(1/sin(x)**2 - cot(x)**2) == 1

    assert trigsimp(5*cos(x)**2 + 5*sin(x)**2) == 5
    assert trigsimp(5*cos(x/2)**2 + 2*sin(x/2)**2) == 3*cos(x)/2 + S(7)/2

    assert trigsimp(sin(x)/cos(x)) == tan(x)
    assert trigsimp(2*tan(x)*cos(x)) == 2*sin(x)
    assert trigsimp(cot(x)**3*sin(x)**3) == cos(x)**3
    assert trigsimp(y*tan(x)**2/sin(x)**2) == y/cos(x)**2
    assert trigsimp(cot(x)/cos(x)) == 1/sin(x)

    assert trigsimp(sin(x + y) + sin(x - y)) == 2*sin(x)*cos(y)
    assert trigsimp(sin(x + y) - sin(x - y)) == 2*sin(y)*cos(x)
    assert trigsimp(cos(x + y) + cos(x - y)) == 2*cos(x)*cos(y)
    assert trigsimp(cos(x + y) - cos(x - y)) == -2*sin(x)*sin(y)
    assert trigsimp(tan(x + y) - tan(x)/(1 - tan(x)*tan(y))) == \
        sin(y)/(-sin(y)*tan(x) + cos(y))  # -tan(y)/(tan(x)*tan(y) - 1)

    assert trigsimp(sinh(x + y) + sinh(x - y)) == 2*sinh(x)*cosh(y)
    assert trigsimp(sinh(x + y) - sinh(x - y)) == 2*sinh(y)*cosh(x)
    assert trigsimp(cosh(x + y) + cosh(x - y)) == 2*cosh(x)*cosh(y)
    assert trigsimp(cosh(x + y) - cosh(x - y)) == 2*sinh(x)*sinh(y)
    assert trigsimp(tanh(x + y) - tanh(x)/(1 + tanh(x)*tanh(y))) == \
        sinh(y)/(sinh(y)*tanh(x) + cosh(y))

    assert trigsimp(cos(0.12345)**2 + sin(0.12345)**2) == 1
    e = 2*sin(x)**2 + 2*cos(x)**2
    assert trigsimp(log(e)) == log(2)


def test_trigsimp1a():
    assert trigsimp(sin(2)**2*cos(3)*exp(2)/cos(2)**2) == tan(2)**2*cos(3)*exp(2)
    assert trigsimp(tan(2)**2*cos(3)*exp(2)*cos(2)**2) == sin(2)**2*cos(3)*exp(2)
    assert trigsimp(cot(2)*cos(3)*exp(2)*sin(2)) == cos(3)*exp(2)*cos(2)
    assert trigsimp(tan(2)*cos(3)*exp(2)/sin(2)) == cos(3)*exp(2)/cos(2)
    assert trigsimp(cot(2)*cos(3)*exp(2)/cos(2)) == cos(3)*exp(2)/sin(2)
    assert trigsimp(cot(2)*cos(3)*exp(2)*tan(2)) == cos(3)*exp(2)
    assert trigsimp(sinh(2)*cos(3)*exp(2)/cosh(2)) == tanh(2)*cos(3)*exp(2)
    assert trigsimp(tanh(2)*cos(3)*exp(2)*cosh(2)) == sinh(2)*cos(3)*exp(2)
    assert trigsimp(coth(2)*cos(3)*exp(2)*sinh(2)) == cosh(2)*cos(3)*exp(2)
    assert trigsimp(tanh(2)*cos(3)*exp(2)/sinh(2)) == cos(3)*exp(2)/cosh(2)
    assert trigsimp(coth(2)*cos(3)*exp(2)/cosh(2)) == cos(3)*exp(2)/sinh(2)
    assert trigsimp(coth(2)*cos(3)*exp(2)*tanh(2)) == cos(3)*exp(2)


def test_trigsimp2():
    x, y = symbols('x,y')
    assert trigsimp(cos(x)**2*sin(y)**2 + cos(x)**2*cos(y)**2 + sin(x)**2,
            recursive=True) == 1
    assert trigsimp(sin(x)**2*sin(y)**2 + sin(x)**2*cos(y)**2 + cos(x)**2,
            recursive=True) == 1
    assert trigsimp(
        Subs(x, x, sin(y)**2 + cos(y)**2)) == Subs(x, x, 1)


def test_issue_4373():
    x = Symbol("x")
    assert abs(trigsimp(2.0*sin(x)**2 + 2.0*cos(x)**2) - 2.0) < 1e-10


def test_trigsimp3():
    x, y = symbols('x,y')
    assert trigsimp(sin(x)/cos(x)) == tan(x)
    assert trigsimp(sin(x)**2/cos(x)**2) == tan(x)**2
    assert trigsimp(sin(x)**3/cos(x)**3) == tan(x)**3
    assert trigsimp(sin(x)**10/cos(x)**10) == tan(x)**10

    assert trigsimp(cos(x)/sin(x)) == 1/tan(x)
    assert trigsimp(cos(x)**2/sin(x)**2) == 1/tan(x)**2
    assert trigsimp(cos(x)**10/sin(x)**10) == 1/tan(x)**10

    assert trigsimp(tan(x)) == trigsimp(sin(x)/cos(x))


def test_issue_4661():
    a, x, y = symbols('a x y')
    eq = -4*sin(x)**4 + 4*cos(x)**4 - 8*cos(x)**2
    assert trigsimp(eq) == -4
    n = sin(x)**6 + 4*sin(x)**4*cos(x)**2 + 5*sin(x)**2*cos(x)**4 + 2*cos(x)**6
    d = -sin(x)**2 - 2*cos(x)**2
    assert simplify(n/d) == -1
    assert trigsimp(-2*cos(x)**2 + cos(x)**4 - sin(x)**4) == -1
    eq = (- sin(x)**3/4)*cos(x) + (cos(x)**3/4)*sin(x) - sin(2*x)*cos(2*x)/8
    assert trigsimp(eq) == 0


def test_issue_4494():
    a, b = symbols('a b')
    eq = sin(a)**2*sin(b)**2 + cos(a)**2*cos(b)**2*tan(a)**2 + cos(a)**2
    assert trigsimp(eq) == 1


def test_issue_5948():
    a, x, y = symbols('a x y')
    assert trigsimp(diff(integrate(cos(x)/sin(x)**7, x), x)) == \
           cos(x)/sin(x)**7


def test_issue_4775():
    a, x, y = symbols('a x y')
    assert trigsimp(sin(x)*cos(y)+cos(x)*sin(y)) == sin(x + y)
    assert trigsimp(sin(x)*cos(y)+cos(x)*sin(y)+3) == sin(x + y) + 3


def test_issue_4280():
    a, x, y = symbols('a x y')
    assert trigsimp(cos(x)**2 + cos(y)**2*sin(x)**2 + sin(y)**2*sin(x)**2) == 1
    assert trigsimp(a**2*sin(x)**2 + a**2*cos(y)**2*cos(x)**2 + a**2*cos(x)**2*sin(y)**2) == a**2
    assert trigsimp(a**2*cos(y)**2*sin(x)**2 + a**2*sin(y)**2*sin(x)**2) == a**2*sin(x)**2


def test_issue_3210():
    eqs = (sin(2)*cos(3) + sin(3)*cos(2),
        -sin(2)*sin(3) + cos(2)*cos(3),
        sin(2)*cos(3) - sin(3)*cos(2),
        sin(2)*sin(3) + cos(2)*cos(3),
        sin(2)*sin(3) + cos(2)*cos(3) + cos(2),
        sinh(2)*cosh(3) + sinh(3)*cosh(2),
        sinh(2)*sinh(3) + cosh(2)*cosh(3),
        )
    assert [trigsimp(e) for e in eqs] == [
        sin(5),
        cos(5),
        -sin(1),
        cos(1),
        cos(1) + cos(2),
        sinh(5),
        cosh(5),
        ]


def test_issue_7263():
    assert abs((simplify(30.8**2 - 82.5**2 * sin(rad(11.6))**2)).evalf() - \
            673.447451402970) < 1e-12


def test_trigsimp_issues():
    a, x, y = symbols('a x y')

    # issue 4625 - factor_terms works, too
    assert trigsimp(sin(x)**3 + cos(x)**2*sin(x)) == sin(x)

    # issue 5948
    assert trigsimp(diff(integrate(cos(x)/sin(x)**3, x), x)) == \
        cos(x)/sin(x)**3
    assert trigsimp(diff(integrate(sin(x)/cos(x)**3, x), x)) == \
        sin(x)/cos(x)**3

    # check integer exponents
    e = sin(x)**y/cos(x)**y
    assert trigsimp(e) == e
    assert trigsimp(e.subs(y, 2)) == tan(x)**2
    assert trigsimp(e.subs(x, 1)) == tan(1)**y

    # check for multiple patterns
    assert (cos(x)**2/sin(x)**2*cos(y)**2/sin(y)**2).trigsimp() == \
        1/tan(x)**2/tan(y)**2
    assert trigsimp(cos(x)/sin(x)*cos(x+y)/sin(x+y)) == \
        1/(tan(x)*tan(x + y))

    eq = cos(2)*(cos(3) + 1)**2/(cos(3) - 1)**2
    assert trigsimp(eq) == eq.factor()  # factor makes denom (-1 + cos(3))**2
    assert trigsimp(cos(2)*(cos(3) + 1)**2*(cos(3) - 1)**2) == \
        cos(2)*sin(3)**4

    # issue 6789; this generates an expression that formerly caused
    # trigsimp to hang
    assert cot(x).equals(tan(x)) is False

    # nan or the unchanged expression is ok, but not sin(1)
    z = cos(x)**2 + sin(x)**2 - 1
    z1 = tan(x)**2 - 1/cot(x)**2
    n = (1 + z1/z)
    assert trigsimp(sin(n)) != sin(1)
    eq = x*(n - 1) - x*n
    assert trigsimp(eq) is S.NaN
    assert trigsimp(eq, recursive=True) is S.NaN
    assert trigsimp(1).is_Integer

    assert trigsimp(-sin(x)**4 - 2*sin(x)**2*cos(x)**2 - cos(x)**4) == -1


def test_trigsimp_issue_2515():
    x = Symbol('x')
    assert trigsimp(x*cos(x)*tan(x)) == x*sin(x)
    assert trigsimp(-sin(x) + cos(x)*tan(x)) == 0


def test_trigsimp_issue_3826():
    assert trigsimp(tan(2*x).expand(trig=True)) == tan(2*x)


def test_trigsimp_issue_4032():
    n = Symbol('n', integer=True, positive=True)
    assert trigsimp(2**(n/2)*cos(pi*n/4)/2 + 2**(n - 1)/2) == \
        2**(n/2)*cos(pi*n/4)/2 + 2**n/4


def test_trigsimp_issue_7761():
    assert trigsimp(cosh(pi/4)) == cosh(pi/4)


def test_trigsimp_noncommutative():
    x, y = symbols('x,y')
    A, B = symbols('A,B', commutative=False)

    assert trigsimp(A - A*sin(x)**2) == A*cos(x)**2
    assert trigsimp(A - A*cos(x)**2) == A*sin(x)**2
    assert trigsimp(A*sin(x)**2 + A*cos(x)**2) == A
    assert trigsimp(A + A*tan(x)**2) == A/cos(x)**2
    assert trigsimp(A/cos(x)**2 - A) == A*tan(x)**2
    assert trigsimp(A/cos(x)**2 - A*tan(x)**2) == A
    assert trigsimp(A + A*cot(x)**2) == A/sin(x)**2
    assert trigsimp(A/sin(x)**2 - A) == A/tan(x)**2
    assert trigsimp(A/sin(x)**2 - A*cot(x)**2) == A

    assert trigsimp(y*A*cos(x)**2 + y*A*sin(x)**2) == y*A

    assert trigsimp(A*sin(x)/cos(x)) == A*tan(x)
    assert trigsimp(A*tan(x)*cos(x)) == A*sin(x)
    assert trigsimp(A*cot(x)**3*sin(x)**3) == A*cos(x)**3
    assert trigsimp(y*A*tan(x)**2/sin(x)**2) == y*A/cos(x)**2
    assert trigsimp(A*cot(x)/cos(x)) == A/sin(x)

    assert trigsimp(A*sin(x + y) + A*sin(x - y)) == 2*A*sin(x)*cos(y)
    assert trigsimp(A*sin(x + y) - A*sin(x - y)) == 2*A*sin(y)*cos(x)
    assert trigsimp(A*cos(x + y) + A*cos(x - y)) == 2*A*cos(x)*cos(y)
    assert trigsimp(A*cos(x + y) - A*cos(x - y)) == -2*A*sin(x)*sin(y)

    assert trigsimp(A*sinh(x + y) + A*sinh(x - y)) == 2*A*sinh(x)*cosh(y)
    assert trigsimp(A*sinh(x + y) - A*sinh(x - y)) == 2*A*sinh(y)*cosh(x)
    assert trigsimp(A*cosh(x + y) + A*cosh(x - y)) == 2*A*cosh(x)*cosh(y)
    assert trigsimp(A*cosh(x + y) - A*cosh(x - y)) == 2*A*sinh(x)*sinh(y)

    assert trigsimp(A*cos(0.12345)**2 + A*sin(0.12345)**2) == 1.0*A


def test_hyperbolic_simp():
    x, y = symbols('x,y')

    assert trigsimp(sinh(x)**2 + 1) == cosh(x)**2
    assert trigsimp(cosh(x)**2 - 1) == sinh(x)**2
    assert trigsimp(cosh(x)**2 - sinh(x)**2) == 1
    assert trigsimp(1 - tanh(x)**2) == 1/cosh(x)**2
    assert trigsimp(1 - 1/cosh(x)**2) == tanh(x)**2
    assert trigsimp(tanh(x)**2 + 1/cosh(x)**2) == 1
    assert trigsimp(coth(x)**2 - 1) == 1/sinh(x)**2
    assert trigsimp(1/sinh(x)**2 + 1) == 1/tanh(x)**2
    assert trigsimp(coth(x)**2 - 1/sinh(x)**2) == 1

    assert trigsimp(5*cosh(x)**2 - 5*sinh(x)**2) == 5
    assert trigsimp(5*cosh(x/2)**2 - 2*sinh(x/2)**2) == 3*cosh(x)/2 + S(7)/2

    assert trigsimp(sinh(x)/cosh(x)) == tanh(x)
    assert trigsimp(tanh(x)) == trigsimp(sinh(x)/cosh(x))
    assert trigsimp(cosh(x)/sinh(x)) == 1/tanh(x)
    assert trigsimp(2*tanh(x)*cosh(x)) == 2*sinh(x)
    assert trigsimp(coth(x)**3*sinh(x)**3) == cosh(x)**3
    assert trigsimp(y*tanh(x)**2/sinh(x)**2) == y/cosh(x)**2
    assert trigsimp(coth(x)/cosh(x)) == 1/sinh(x)

    e = 2*cosh(x)**2 - 2*sinh(x)**2
    assert trigsimp(log(e)) == log(2)

    assert trigsimp(cosh(x)**2*cosh(y)**2 - cosh(x)**2*sinh(y)**2 - sinh(x)**2,
            recursive=True) == 1
    assert trigsimp(sinh(x)**2*sinh(y)**2 - sinh(x)**2*cosh(y)**2 + cosh(x)**2,
            recursive=True) == 1

    assert abs(trigsimp(2.0*cosh(x)**2 - 2.0*sinh(x)**2) - 2.0) < 1e-10

    assert trigsimp(sinh(x)**2/cosh(x)**2) == tanh(x)**2
    assert trigsimp(sinh(x)**3/cosh(x)**3) == tanh(x)**3
    assert trigsimp(sinh(x)**10/cosh(x)**10) == tanh(x)**10
    assert trigsimp(cosh(x)**3/sinh(x)**3) == 1/tanh(x)**3

    assert trigsimp(cosh(x)/sinh(x)) == 1/tanh(x)
    assert trigsimp(cosh(x)**2/sinh(x)**2) == 1/tanh(x)**2
    assert trigsimp(cosh(x)**10/sinh(x)**10) == 1/tanh(x)**10

    assert trigsimp(x*cosh(x)*tanh(x)) == x*sinh(x)
    assert trigsimp(-sinh(x) + cosh(x)*tanh(x)) == 0

    assert tan(x) != 1/cot(x)  # cot doesn't auto-simplify

    assert trigsimp(tan(x) - 1/cot(x)) == 0
    assert trigsimp(3*tanh(x)**7 - 2/coth(x)**7) == tanh(x)**7


def test_trigsimp_groebner():
    from sympy.simplify.simplify import trigsimp_groebner

    c = cos(x)
    s = sin(x)
    ex = (4*s*c + 12*s + 5*c**3 + 21*c**2 + 23*c + 15)/(
        -s*c**2 + 2*s*c + 15*s + 7*c**3 + 31*c**2 + 37*c + 21)
    resnum = (5*s - 5*c + 1)
    resdenom = (8*s - 6*c)
    results = [resnum/resdenom, (-resnum)/(-resdenom)]
    assert trigsimp_groebner(ex) in results
    assert trigsimp_groebner(s/c, hints=[tan]) == tan(x)

    assert trigsimp((-s + 1)/c + c/(-s + 1),
                    method='groebner') == 2/c
    assert trigsimp((-s + 1)/c + c/(-s + 1),
                    method='groebner', polynomial=True) == 2/c

    # Test quick=False works
    assert trigsimp_groebner(ex, hints=[2]) in results

    # test "I"
    assert trigsimp_groebner(sin(I*x)/cos(I*x), hints=[tanh]) == I*tanh(x)

    # test hyperbolic / sums
    assert trigsimp_groebner((tanh(x)+tanh(y))/(1+tanh(x)*tanh(y)),
                             hints=[(tanh, x, y)]) == tanh(x + y)


@XFAIL
def test_factorial_simplify():
    # There are more tests in test_factorials.py. These are just to
    # ensure that simplify() calls factorial_simplify correctly
    from sympy.specfun.factorials import factorial
    x = Symbol('x')
    assert simplify(factorial(x)/x) == factorial(x - 1)
    assert simplify(factorial(factorial(x))) == factorial(factorial(x))


def test_simplify_expr():
    x, y, z, k, n, m, w, f, s, A = symbols('x,y,z,k,n,m,w,f,s,A')

    assert all(simplify(tmp) == tmp for tmp in [I, E, oo, x, -x, -oo, -E, -I])

    e = 1/x + 1/y
    assert e != (x + y)/(x*y)
    assert simplify(e) == (x + y)/(x*y)

    e = A**2*s**4/(4*pi*k*m**3)
    assert simplify(e) == e

    e = (4 + 4*x - 2*(2 + 2*x))/(2 + 2*x)
    assert simplify(e) == 0

    e = (-4*x*y**2 - 2*y**3 - 2*x**2*y)/(x + y)**2
    assert simplify(e) == -2*y

    e = -x - y - (x + y)**(-1)*y**2 + (x + y)**(-1)*x**2
    assert simplify(e) == -2*y

    e = (x + x*y)/x
    assert simplify(e) == 1 + y

    e = (f(x) + y*f(x))/f(x)
    assert simplify(e) == 1 + y

    e = (2 * (1/n - cos(n * pi)/n))/pi
    assert simplify(e) == (-cos(pi*n) + 1)/(pi*n)*2

    e = integrate(1/(x**3 + 1), x).diff(x)
    assert simplify(e) == 1/(x**3 + 1)

    e = integrate(x/(x**2 + 3*x + 1), x).diff(x)
    assert simplify(e) == x/(x**2 + 3*x + 1)

    A = Matrix([[2*k - m*w**2, -k], [-k, k - m*w**2]]).inv()
    assert simplify((A*Matrix([0, f]))[1]) == \
        -f*(2*k - m*w**2)/(k**2 - (k - m*w**2)*(2*k - m*w**2))

    f = -x + y/(z + t) + z*x/(z + t) + z*a/(z + t) + t*x/(z + t)
    assert simplify(f) == (y + a*z)/(z + t)

    A, B = symbols('A,B', commutative=False)

    assert simplify(A*B - B*A) == A*B - B*A
    assert simplify(A/(1 + y/x)) == x*A/(x + y)
    assert simplify(A*(1/x + 1/y)) == A/x + A/y  #(x + y)*A/(x*y)

    assert simplify(log(2) + log(3)) == log(6)
    assert simplify(log(2*x) - log(2)) == log(x)

    assert simplify(hyper([], [], x)) == exp(x)


def test_issue_3557():
    f_1 = x*a + y*b + z*c - 1
    f_2 = x*d + y*e + z*f - 1
    f_3 = x*g + y*h + z*i - 1

    solutions = solve([f_1, f_2, f_3], x, y, z, simplify=False)

    assert simplify(solutions[y]) == \
        (a*i + c*d + f*g - a*f - c*g - d*i)/ \
        (a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g)


def test_simplify_other():
    assert simplify(sin(x)**2 + cos(x)**2) == 1
    assert simplify(gamma(x + 1)/gamma(x)) == x
    assert simplify(sin(x)**2 + cos(x)**2 + factorial(x)/gamma(x)) == 1 + x
    assert simplify(
        Eq(sin(x)**2 + cos(x)**2, factorial(x)/gamma(x))) == Eq(x, 1)
    nc = symbols('nc', commutative=False)
    assert simplify(x + x*nc) == x*(1 + nc)
    # issue 6123
    # f = exp(-I*(k*sqrt(t) + x/(2*sqrt(t)))**2)
    # ans = integrate(f, (k, -oo, oo), conds='none')
    ans = I*(-pi*x*exp(-3*I*pi/4 + I*x**2/(4*t))*erf(x*exp(-3*I*pi/4)/
        (2*sqrt(t)))/(2*sqrt(t)) + pi*x*exp(-3*I*pi/4 + I*x**2/(4*t))/
        (2*sqrt(t)))*exp(-I*x**2/(4*t))/(sqrt(pi)*x) - I*sqrt(pi) * \
        (-erf(x*exp(I*pi/4)/(2*sqrt(t))) + 1)*exp(I*pi/4)/(2*sqrt(t))
    assert simplify(ans) == -(-1)**(S(3)/4)*sqrt(pi)/sqrt(t)
    # issue 6370
    assert simplify(2**(2 + x)/4) == 2**x


def test_simplify_complex():
    cosAsExp = cos(x)._eval_rewrite_as_exp(x)
    tanAsExp = tan(x)._eval_rewrite_as_exp(x)
    assert simplify(cosAsExp*tanAsExp).expand() == (
        sin(x))._eval_rewrite_as_exp(x).expand()  # issue 4341


def test_simplify_ratio():
    # roots of x**3-3*x+5
    roots = ['(1/2 - sqrt(3)*I/2)*(sqrt(21)/2 + 5/2)**(1/3) + 1/((1/2 - '
             'sqrt(3)*I/2)*(sqrt(21)/2 + 5/2)**(1/3))',
             '1/((1/2 + sqrt(3)*I/2)*(sqrt(21)/2 + 5/2)**(1/3)) + '
             '(1/2 + sqrt(3)*I/2)*(sqrt(21)/2 + 5/2)**(1/3)',
             '-(sqrt(21)/2 + 5/2)**(1/3) - 1/(sqrt(21)/2 + 5/2)**(1/3)']

    for r in roots:
        r = S(r)
        assert count_ops(simplify(r, ratio=1)) <= count_ops(r)
        # If ratio=oo, simplify() is always applied:
        assert simplify(r, ratio=oo) is not r


def test_simplify_measure():
    measure1 = lambda expr: len(str(expr))
    measure2 = lambda expr: -count_ops(expr)
                                       # Return the most complicated result
    expr = (x + 1)/(x + sin(x)**2 + cos(x)**2)
    assert measure1(simplify(expr, measure=measure1)) <= measure1(expr)
    assert measure2(simplify(expr, measure=measure2)) <= measure2(expr)

    expr2 = Eq(sin(x)**2 + cos(x)**2, 1)
    assert measure1(simplify(expr2, measure=measure1)) <= measure1(expr2)
    assert measure2(simplify(expr2, measure=measure2)) <= measure2(expr2)


def test_simplify_issue_1308():
    assert simplify(exp(-Rational(1, 2)) + exp(-Rational(3, 2))) == \
        (1 + E)*exp(-Rational(3, 2))


def test_issue_5652():
    assert simplify(E + exp(-E)) == exp(-E) + E
    n = symbols('n', commutative=False)
    assert simplify(n + n**(-n)) == n + n**(-n)


def test_simplify_fail1():
    x = Symbol('x')
    y = Symbol('y')
    e = (x + y)**2/(-4*x*y**2 - 2*y**3 - 2*x**2*y)
    assert simplify(e) == 1 / (-2*y)


def test_fraction():
    x, y, z = map(Symbol, 'xyz')
    A = Symbol('A', commutative=False)

    assert fraction(Rational(1, 2)) == (1, 2)

    assert fraction(x) == (x, 1)
    assert fraction(1/x) == (1, x)
    assert fraction(x/y) == (x, y)
    assert fraction(x/2) == (x, 2)

    assert fraction(x*y/z) == (x*y, z)
    assert fraction(x/(y*z)) == (x, y*z)

    assert fraction(1/y**2) == (1, y**2)
    assert fraction(x/y**2) == (x, y**2)

    assert fraction((x**2 + 1)/y) == (x**2 + 1, y)
    assert fraction(x*(y + 1)/y**7) == (x*(y + 1), y**7)

    assert fraction(exp(-x), exact=True) == (exp(-x), 1)

    assert fraction(x*A/y) == (x*A, y)
    assert fraction(x*A**-1/y) == (x*A**-1, y)

    n = symbols('n', negative=True)
    assert fraction(exp(n)) == (1, exp(-n))
    assert fraction(exp(-n)) == (exp(-n), 1)


def test_powsimp():
    x, y, z, n = symbols('x,y,z,n')
    f = Function('f')
    assert powsimp( 4**x * 2**(-x) * 2**(-x) ) == 1
    assert powsimp( (-4)**x * (-2)**(-x) * 2**(-x) ) == 1

    assert powsimp(
        f(4**x * 2**(-x) * 2**(-x)) ) == f(4**x * 2**(-x) * 2**(-x))
    assert powsimp( f(4**x * 2**(-x) * 2**(-x)), deep=True ) == f(1)
    assert exp(x)*exp(y) == exp(x)*exp(y)
    assert powsimp(exp(x)*exp(y)) == exp(x + y)
    assert powsimp(exp(x)*exp(y)*2**x*2**y) == (2*E)**(x + y)
    assert powsimp(exp(x)*exp(y)*2**x*2**y, combine='exp') == \
        exp(x + y)*2**(x + y)
    assert powsimp(exp(x)*exp(y)*exp(2)*sin(x) + sin(y) + 2**x*2**y) == \
        exp(2 + x + y)*sin(x) + sin(y) + 2**(x + y)
    assert powsimp(sin(exp(x)*exp(y))) == sin(exp(x)*exp(y))
    assert powsimp(sin(exp(x)*exp(y)), deep=True) == sin(exp(x + y))
    assert powsimp(x**2*x**y) == x**(2 + y)
    # This should remain factored, because 'exp' with deep=True is supposed
    # to act like old automatic exponent combining.
    assert powsimp((1 + E*exp(E))*exp(-E), combine='exp', deep=True) == \
        (1 + exp(1 + E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), deep=True) == \
        (1 + exp(1 + E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E)) == (1 + exp(1 + E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), combine='exp') == \
        (1 + exp(1 + E))*exp(-E)
    assert powsimp((1 + E*exp(E))*exp(-E), combine='base') == \
        (1 + E*exp(E))*exp(-E)
    x, y = symbols('x,y', nonnegative=True)
    n = Symbol('n', real=True)
    assert powsimp(y**n * (y/x)**(-n)) == x**n
    assert powsimp(x**(x**(x*y)*y**(x*y))*y**(x**(x*y)*y**(x*y)), deep=True) \
        == (x*y)**(x*y)**(x*y)
    assert powsimp(2**(2**(2*x)*x), deep=False) == 2**(2**(2*x)*x)
    assert powsimp(2**(2**(2*x)*x), deep=True) == 2**(x*4**x)
    assert powsimp(
        exp(-x + exp(-x)*exp(-x*log(x))), deep=False, combine='exp') == \
        exp(-x + exp(-x)*exp(-x*log(x)))
    assert powsimp(
        exp(-x + exp(-x)*exp(-x*log(x))), deep=False, combine='exp') == \
        exp(-x + exp(-x)*exp(-x*log(x)))
    assert powsimp((x + y)/(3*z), deep=False, combine='exp') == (x + y)/(3*z)
    assert powsimp((x/3 + y/3)/z, deep=True, combine='exp') == (x/3 + y/3)/z
    assert powsimp(exp(x)/(1 + exp(x)*exp(y)), deep=True) == \
        exp(x)/(1 + exp(x + y))
    assert powsimp(x*y**(z**x*z**y), deep=True) == x*y**(z**(x + y))
    assert powsimp((z**x*z**y)**x, deep=True) == (z**(x + y))**x
    assert powsimp(x*(z**x*z**y)**x, deep=True) == x*(z**(x + y))**x
    p = symbols('p', positive=True)
    assert powsimp((1/x)**log(2)/x) == (1/x)**(1 + log(2))
    assert powsimp((1/p)**log(2)/p) == p**(-1 - log(2))

    # coefficient of exponent can only be simplified for positive bases
    assert powsimp(2**(2*x)) == 4**x
    assert powsimp((-1)**(2*x)) == (-1)**(2*x)
    i = symbols('i', integer=True)
    assert powsimp((-1)**(2*i)) == 1
    assert powsimp((-1)**(-x)) != (-1)**x  # could be 1/((-1)**x), but is not
    # force=True overrides assumptions
    assert powsimp((-1)**(2*x), force=True) == 1

    # rational exponents allow combining of negative terms
    w, n, m = symbols('w n m', negative=True)
    e = i/a  # not a rational exponent if `a` is unknown
    ex = w**e*n**e*m**e
    assert powsimp(ex) == m**(i/a)*n**(i/a)*w**(i/a)
    e = i/3
    ex = w**e*n**e*m**e
    assert powsimp(ex) == (-1)**i*(-m*n*w)**(i/3)
    e = (3 + i)/i
    ex = w**e*n**e*m**e
    assert powsimp(ex) == (-1)**(3*e)*(-m*n*w)**e

    eq = x**(2*a/3)
    # eq != (x**a)**(2/3) (try x = -1 and a = 3 to see)
    assert powsimp(eq).exp == eq.exp == 2*a/3
    # powdenest goes the other direction
    assert powsimp(2**(2*x)) == 4**x

    assert powsimp(exp(p/2)) == exp(p/2)

    # issue 6368
    eq = Mul(*[sqrt(Dummy(imaginary=True)) for i in range(3)])
    assert powsimp(eq) == eq and eq.is_Mul

    assert all(powsimp(e) == e for e in (sqrt(x**a), sqrt(x**2)))

    # issue 8836
    assert str( powsimp(exp(I*pi/3)*root(-1,3)) ) == '(-1)**(2/3)'

def test_issue_6367():
    z = -5*sqrt(2)/(2*sqrt(2*sqrt(29) + 29)) + sqrt(-sqrt(29)/29 + S(1)/2)
    assert Mul(*[powsimp(a) for a in Mul.make_args(z.normal())]) == 0
    assert powsimp(z.normal()) == 0
    assert simplify(z) == 0
    assert powsimp(sqrt(2 + sqrt(3))*sqrt(2 - sqrt(3)) + 1) == 2
    assert powsimp(z) != 0


def test_powsimp_negated_base():
    assert powsimp((-x + y)/sqrt(x - y)) == -sqrt(x - y)
    assert powsimp((-x + y)*(-z + y)/sqrt(x - y)/sqrt(z - y)) == sqrt(x - y)*sqrt(z - y)
    p = symbols('p', positive=True)
    assert powsimp((-p)**a/p**a) == (-1)**a
    n = symbols('n', negative=True)
    assert powsimp((-n)**a/n**a) == (-1)**a
    # if x is 0 then the lhs is 0**a*oo**a which is not (-1)**a
    assert powsimp((-x)**a/x**a) != (-1)**a


def test_issue_6440():
    assert powsimp(16*2**a*8**b) == 2**(a + 3*b + 4)


def test_powsimp_polar():
    from sympy import polar_lift, exp_polar
    x, y, z = symbols('x y z')
    p, q, r = symbols('p q r', polar=True)

    assert (polar_lift(-1))**(2*x) == exp_polar(2*pi*I*x)
    assert powsimp(p**x * q**x) == (p*q)**x
    assert p**x * (1/p)**x == 1
    assert (1/p)**x == p**(-x)

    assert exp_polar(x)*exp_polar(y) == exp_polar(x)*exp_polar(y)
    assert powsimp(exp_polar(x)*exp_polar(y)) == exp_polar(x + y)
    assert powsimp(exp_polar(x)*exp_polar(y)*p**x*p**y) == \
        (p*exp_polar(1))**(x + y)
    assert powsimp(exp_polar(x)*exp_polar(y)*p**x*p**y, combine='exp') == \
        exp_polar(x + y)*p**(x + y)
    assert powsimp(
        exp_polar(x)*exp_polar(y)*exp_polar(2)*sin(x) + sin(y) + p**x*p**y) \
        == p**(x + y) + sin(x)*exp_polar(2 + x + y) + sin(y)
    assert powsimp(sin(exp_polar(x)*exp_polar(y))) == \
        sin(exp_polar(x)*exp_polar(y))
    assert powsimp(sin(exp_polar(x)*exp_polar(y)), deep=True) == \
        sin(exp_polar(x + y))


def test_powsimp_nc():
    x, y, z = symbols('x,y,z')
    A, B, C = symbols('A B C', commutative=False)

    assert powsimp(A**x*A**y, combine='all') == A**(x + y)
    assert powsimp(A**x*A**y, combine='base') == A**x*A**y
    assert powsimp(A**x*A**y, combine='exp') == A**(x + y)

    assert powsimp(A**x*B**x, combine='all') == (A*B)**x
    assert powsimp(A**x*B**x, combine='base') == (A*B)**x
    assert powsimp(A**x*B**x, combine='exp') == A**x*B**x

    assert powsimp(B**x*A**x, combine='all') == (B*A)**x
    assert powsimp(B**x*A**x, combine='base') == (B*A)**x
    assert powsimp(B**x*A**x, combine='exp') == B**x*A**x

    assert powsimp(A**x*A**y*A**z, combine='all') == A**(x + y + z)
    assert powsimp(A**x*A**y*A**z, combine='base') == A**x*A**y*A**z
    assert powsimp(A**x*A**y*A**z, combine='exp') == A**(x + y + z)

    assert powsimp(A**x*B**x*C**x, combine='all') == (A*B*C)**x
    assert powsimp(A**x*B**x*C**x, combine='base') == (A*B*C)**x
    assert powsimp(A**x*B**x*C**x, combine='exp') == A**x*B**x*C**x

    assert powsimp(B**x*A**x*C**x, combine='all') == (B*A*C)**x
    assert powsimp(B**x*A**x*C**x, combine='base') == (B*A*C)**x
    assert powsimp(B**x*A**x*C**x, combine='exp') == B**x*A**x*C**x


def test_nthroot():
    assert nthroot(90 + 34*sqrt(7), 3) == sqrt(7) + 3
    q = 1 + sqrt(2) - 2*sqrt(3) + sqrt(6) + sqrt(7)
    assert nthroot(expand_multinomial(q**3), 3) == q
    assert nthroot(41 + 29*sqrt(2), 5) == 1 + sqrt(2)
    assert nthroot(-41 - 29*sqrt(2), 5) == -1 - sqrt(2)
    expr = 1320*sqrt(10) + 4216 + 2576*sqrt(6) + 1640*sqrt(15)
    assert nthroot(expr, 5) == 1 + sqrt(6) + sqrt(15)
    q = 1 + sqrt(2) + sqrt(3) + sqrt(5)
    assert expand_multinomial(nthroot(expand_multinomial(q**5), 5)) == q
    q = 1 + sqrt(2) + 7*sqrt(6) + 2*sqrt(10)
    assert nthroot(expand_multinomial(q**5), 5, 8) == q
    q = 1 + sqrt(2) - 2*sqrt(3) + 1171*sqrt(6)
    assert nthroot(expand_multinomial(q**3), 3) == q
    assert nthroot(expand_multinomial(q**6), 6) == q


@slow
def test_nthroot1():
    q = 1 + sqrt(2) + sqrt(3) + S(1)/10**20
    p = expand_multinomial(q**5)
    assert nthroot(p, 5) == q
    q = 1 + sqrt(2) + sqrt(3) + S(1)/10**30
    p = expand_multinomial(q**5)
    assert nthroot(p, 5) == q


def test_collect_1():
    """Collect with respect to a Symbol"""
    x, y, z, n = symbols('x,y,z,n')
    assert collect( x + y*x, x ) == x * (1 + y)
    assert collect( x + x**2, x ) == x + x**2
    assert collect( x**2 + y*x**2, x ) == (x**2)*(1 + y)
    assert collect( x**2 + y*x, x ) == x*y + x**2
    assert collect( 2*x**2 + y*x**2 + 3*x*y, [x] ) == x**2*(2 + y) + 3*x*y
    assert collect( 2*x**2 + y*x**2 + 3*x*y, [y] ) == 2*x**2 + y*(x**2 + 3*x)

    assert collect( ((1 + y + x)**4).expand(), x) == ((1 + y)**4).expand() + \
        x*(4*(1 + y)**3).expand() + x**2*(6*(1 + y)**2).expand() + \
        x**3*(4*(1 + y)).expand() + x**4
    # symbols can be given as any iterable
    expr = x + y
    assert collect(expr, expr.free_symbols) == expr


def test_collect_2():
    """Collect with respect to a sum"""
    a, b, x = symbols('a,b,x')
    assert collect(a*(cos(x) + sin(x)) + b*(cos(x) + sin(x)),
        sin(x) + cos(x)) == (a + b)*(cos(x) + sin(x))


def test_collect_3():
    """Collect with respect to a product"""
    a, b, c = symbols('a,b,c')
    f = Function('f')
    x, y, z, n = symbols('x,y,z,n')

    assert collect(-x/8 + x*y, -x) == x*(y - S(1)/8)

    assert collect( 1 + x*(y**2), x*y ) == 1 + x*(y**2)
    assert collect( x*y + a*x*y, x*y) == x*y*(1 + a)
    assert collect( 1 + x*y + a*x*y, x*y) == 1 + x*y*(1 + a)
    assert collect(a*x*f(x) + b*(x*f(x)), x*f(x)) == x*(a + b)*f(x)

    assert collect(a*x*log(x) + b*(x*log(x)), x*log(x)) == x*(a + b)*log(x)
    assert collect(a*x**2*log(x)**2 + b*(x*log(x))**2, x*log(x)) == \
        x**2*log(x)**2*(a + b)

    # with respect to a product of three symbols
    assert collect(y*x*z + a*x*y*z, x*y*z) == (1 + a)*x*y*z


def test_collect_4():
    """Collect with respect to a power"""
    a, b, c, x = symbols('a,b,c,x')

    assert collect(a*x**c + b*x**c, x**c) == x**c*(a + b)
    # issue 6096: 2 stays with c (unless c is integer or x is positive0
    assert collect(a*x**(2*c) + b*x**(2*c), x**c) == x**(2*c)*(a + b)


def test_collect_5():
    """Collect with respect to a tuple"""
    a, x, y, z, n = symbols('a,x,y,z,n')
    assert collect(x**2*y**4 + z*(x*y**2)**2 + z + a*z, [x*y**2, z]) in [
        z*(1 + a + x**2*y**4) + x**2*y**4,
        z*(1 + a) + x**2*y**4*(1 + z) ]
    assert collect((1 + (x + y) + (x + y)**2).expand(),
                   [x, y]) == 1 + y + x*(1 + 2*y) + x**2 + y**2


def test_collect_D():
    D = Derivative
    f = Function('f')
    x, a, b = symbols('x,a,b')
    fx = D(f(x), x)
    fxx = D(f(x), x, x)

    assert collect(a*fx + b*fx, fx) == (a + b)*fx
    assert collect(a*D(fx, x) + b*D(fx, x), fx) == (a + b)*D(fx, x)
    assert collect(a*fxx + b*fxx, fx) == (a + b)*D(fx, x)
    # issue 4784
    assert collect(5*f(x) + 3*fx, fx) == 5*f(x) + 3*fx
    assert collect(f(x) + f(x)*diff(f(x), x) + x*diff(f(x), x)*f(x), f(x).diff(x)) == \
        (x*f(x) + f(x))*D(f(x), x) + f(x)
    assert collect(f(x) + f(x)*diff(f(x), x) + x*diff(f(x), x)*f(x), f(x).diff(x), exact=True) == \
        (x*f(x) + f(x))*D(f(x), x) + f(x)
    assert collect(1/f(x) + 1/f(x)*diff(f(x), x) + x*diff(f(x), x)/f(x), f(x).diff(x), exact=True) == \
        (1/f(x) + x/f(x))*D(f(x), x) + 1/f(x)


@XFAIL
def test_collect_issues():
    D = Derivative
    f = Function('f')
    e = (1 + x*D(f(x), x) + D(f(x), x))/f(x)
    assert collect(e.expand(), f(x).diff(x)) != e


def test_collect_D_0():
    D = Derivative
    f = Function('f')
    x, a, b = symbols('x,a,b')
    fxx = D(f(x), x, x)

    # collect does not distinguish nested derivatives, so it returns
    #                                           -- (a + b)*D(D(f, x), x)
    assert collect(a*fxx + b*fxx, fxx) == (a + b)*fxx


def test_collect_Wild():
    """Collect with respect to functions with Wild argument"""
    a, b, x, y = symbols('a b x y')
    f = Function('f')
    w1 = Wild('.1')
    w2 = Wild('.2')
    assert collect(f(x) + a*f(x), f(w1)) == (1 + a)*f(x)
    assert collect(f(x, y) + a*f(x, y), f(w1)) == f(x, y) + a*f(x, y)
    assert collect(f(x, y) + a*f(x, y), f(w1, w2)) == (1 + a)*f(x, y)
    assert collect(f(x, y) + a*f(x, y), f(w1, w1)) == f(x, y) + a*f(x, y)
    assert collect(f(x, x) + a*f(x, x), f(w1, w1)) == (1 + a)*f(x, x)
    assert collect(a*(x + 1)**y + (x + 1)**y, w1**y) == (1 + a)*(x + 1)**y
    assert collect(a*(x + 1)**y + (x + 1)**y, w1**b) == \
        a*(x + 1)**y + (x + 1)**y
    assert collect(a*(x + 1)**y + (x + 1)**y, (x + 1)**w2) == \
        (1 + a)*(x + 1)**y
    assert collect(a*(x + 1)**y + (x + 1)**y, w1**w2) == (1 + a)*(x + 1)**y


def test_collect_func():
    f = ((x + a + 1)**3).expand()

    assert collect(f, x) == a**3 + 3*a**2 + 3*a + x**3 + x**2*(3*a + 3) + \
        x*(3*a**2 + 6*a + 3) + 1
    assert collect(f, x, factor) == x**3 + 3*x**2*(a + 1) + 3*x*(a + 1)**2 + \
        (a + 1)**3

    assert collect(f, x, evaluate=False) == {
        S.One: a**3 + 3*a**2 + 3*a + 1,
        x: 3*a**2 + 6*a + 3, x**2: 3*a + 3,
        x**3: 1
    }


@XFAIL
def test_collect_func_xfail():
    # XXX: this test will pass when automatic constant distribution is removed (issue 4596)
    assert collect(f, x, factor, evaluate=False) == {S.One: (a + 1)**3,
                   x: 3*(a + 1)**2, x**2: 3*(a + 1), x**3: 1}


def test_collect_order():
    a, b, x, t = symbols('a,b,x,t')

    assert collect(t + t*x + t*x**2 + O(x**3), t) == t*(1 + x + x**2 + O(x**3))
    assert collect(t + t*x + x**2 + O(x**3), t) == \
        t*(1 + x + O(x**3)) + x**2 + O(x**3)

    f = a*x + b*x + c*x**2 + d*x**2 + O(x**3)
    g = x*(a + b) + x**2*(c + d) + O(x**3)

    assert collect(f, x) == g
    assert collect(f, x, distribute_order_term=False) == g

    f = sin(a + b).series(b, 0, 10)

    assert collect(f, [sin(a), cos(a)]) == \
        sin(a)*cos(b).series(b, 0, 10) + cos(a)*sin(b).series(b, 0, 10)
    assert collect(f, [sin(a), cos(a)], distribute_order_term=False) == \
        sin(a)*cos(b).series(b, 0, 10).removeO() + \
        cos(a)*sin(b).series(b, 0, 10).removeO() + O(b**10)


def test_rcollect():
    assert rcollect((x**2*y + x*y + x + y)/(x + y), y) == \
        (x + y*(1 + x + x**2))/(x + y)
    assert rcollect(sqrt(-((x + 1)*(y + 1))), z) == sqrt(-((x + 1)*(y + 1)))


def test_separatevars():
    x, y, z, n = symbols('x,y,z,n')
    assert separatevars(2*n*x*z + 2*x*y*z) == 2*x*z*(n + y)
    assert separatevars(x*z + x*y*z) == x*z*(1 + y)
    assert separatevars(pi*x*z + pi*x*y*z) == pi*x*z*(1 + y)
    assert separatevars(x*y**2*sin(x) + x*sin(x)*sin(y)) == \
        x*(sin(y) + y**2)*sin(x)
    assert separatevars(x*exp(x + y) + x*exp(x)) == x*(1 + exp(y))*exp(x)
    assert separatevars((x*(y + 1))**z).is_Pow  # != x**z*(1 + y)**z
    assert separatevars(1 + x + y + x*y) == (x + 1)*(y + 1)
    assert separatevars(y/pi*exp(-(z - x)/cos(n))) == \
        y*exp(x/cos(n))*exp(-z/cos(n))/pi
    assert separatevars((x + y)*(x - y) + y**2 + 2*x + 1) == (x + 1)**2
    # issue 4858
    p = Symbol('p', positive=True)
    assert separatevars(sqrt(p**2 + x*p**2)) == p*sqrt(1 + x)
    assert separatevars(sqrt(y*(p**2 + x*p**2))) == p*sqrt(y*(1 + x))
    assert separatevars(sqrt(y*(p**2 + x*p**2)), force=True) == \
        p*sqrt(y)*sqrt(1 + x)
    # issue 4865
    assert separatevars(sqrt(x*y)).is_Pow
    assert separatevars(sqrt(x*y), force=True) == sqrt(x)*sqrt(y)
    # issue 4957
    # any type sequence for symbols is fine
    assert separatevars(((2*x + 2)*y), dict=True, symbols=()) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    # separable
    assert separatevars(((2*x + 2)*y), dict=True, symbols=[x]) == \
        {'coeff': y, x: 2*x + 2}
    assert separatevars(((2*x + 2)*y), dict=True, symbols=[]) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    assert separatevars(((2*x + 2)*y), dict=True) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    assert separatevars(((2*x + 2)*y), dict=True, symbols=None) == \
        {'coeff': y*(2*x + 2)}
    # not separable
    assert separatevars(3, dict=True) is None
    assert separatevars(2*x + y, dict=True, symbols=()) is None
    assert separatevars(2*x + y, dict=True) is None
    assert separatevars(2*x + y, dict=True, symbols=None) == {'coeff': 2*x + y}
    # issue 4808
    n, m = symbols('n,m', commutative=False)
    assert separatevars(m + n*m) == (1 + n)*m
    assert separatevars(x + x*n) == x*(1 + n)
    # issue 4910
    f = Function('f')
    assert separatevars(f(x) + x*f(x)) == f(x) + x*f(x)
    # a noncommutable object present
    eq = x*(1 + hyper((), (), y*z))
    assert separatevars(eq) == eq


def test_separatevars_advanced_factor():
    x, y, z = symbols('x,y,z')
    assert separatevars(1 + log(x)*log(y) + log(x) + log(y)) == \
        (log(x) + 1)*(log(y) + 1)
    assert separatevars(1 + x - log(z) - x*log(z) - exp(y)*log(z) -
        x*exp(y)*log(z) + x*exp(y) + exp(y)) == \
        -((x + 1)*(log(z) - 1)*(exp(y) + 1))
    x, y = symbols('x,y', positive=True)
    assert separatevars(1 + log(x**log(y)) + log(x*y)) == \
        (log(x) + 1)*(log(y) + 1)


def test_hypersimp():
    n, k = symbols('n,k', integer=True)

    assert hypersimp(factorial(k), k) == k + 1
    assert hypersimp(factorial(k**2), k) is None

    assert hypersimp(1/factorial(k), k) == 1/(k + 1)

    assert hypersimp(2**k/factorial(k)**2, k) == 2/(k + 1)**2

    assert hypersimp(binomial(n, k), k) == (n - k)/(k + 1)
    assert hypersimp(binomial(n + 1, k), k) == (n - k + 1)/(k + 1)

    term = (4*k + 1)*factorial(k)/factorial(2*k + 1)
    assert hypersimp(term, k) == (S(1)/2)*((4*k + 5)/(3 + 14*k + 8*k**2))

    term = 1/((2*k - 1)*factorial(2*k + 1))
    assert hypersimp(term, k) == (k - S(1)/2)/((k + 1)*(2*k + 1)*(2*k + 3))

    term = binomial(n, k)*(-1)**k/factorial(k)
    assert hypersimp(term, k) == (k - n)/(k + 1)**2


def test_nsimplify():
    x = Symbol("x")
    assert nsimplify(0) == 0
    assert nsimplify(-1) == -1
    assert nsimplify(1) == 1
    assert nsimplify(1 + x) == 1 + x
    assert nsimplify(2.7) == Rational(27, 10)
    assert nsimplify(1 - GoldenRatio) == (1 - sqrt(5))/2
    assert nsimplify((1 + sqrt(5))/4, [GoldenRatio]) == GoldenRatio/2
    assert nsimplify(2/GoldenRatio, [GoldenRatio]) == 2*GoldenRatio - 2
    assert nsimplify(exp(5*pi*I/3, evaluate=False)) == \
        sympify('1/2 - sqrt(3)*I/2')
    assert nsimplify(sin(3*pi/5, evaluate=False)) == \
        sympify('sqrt(sqrt(5)/8 + 5/8)')
    assert nsimplify(sqrt(atan('1', evaluate=False))*(2 + I), [pi]) == \
        sqrt(pi) + sqrt(pi)/2*I
    assert nsimplify(2 + exp(2*atan('1/4')*I)) == sympify('49/17 + 8*I/17')
    assert nsimplify(pi, tolerance=0.01) == Rational(22, 7)
    assert nsimplify(pi, tolerance=0.001) == Rational(355, 113)
    assert nsimplify(0.33333, tolerance=1e-4) == Rational(1, 3)
    assert nsimplify(2.0**(1/3.), tolerance=0.001) == Rational(635, 504)
    assert nsimplify(2.0**(1/3.), tolerance=0.001, full=True) == \
        2**Rational(1, 3)
    assert nsimplify(x + .5, rational=True) == Rational(1, 2) + x
    assert nsimplify(1/.3 + x, rational=True) == Rational(10, 3) + x
    assert nsimplify(log(3).n(), rational=True) == \
        sympify('109861228866811/100000000000000')
    assert nsimplify(Float(0.272198261287950), [pi, log(2)]) == pi*log(2)/8
    assert nsimplify(Float(0.272198261287950).n(3), [pi, log(2)]) == \
        -pi/4 - log(2) + S(7)/4
    assert nsimplify(x/7.0) == x/7
    assert nsimplify(pi/1e2) == pi/100
    assert nsimplify(pi/1e2, rational=False) == pi/100.0
    assert nsimplify(pi/1e-7) == 10000000*pi
    assert not nsimplify(
        factor(-3.0*z**2*(z**2)**(-2.5) + 3*(z**2)**(-1.5))).atoms(Float)
    e = x**0.0
    assert e.is_Pow and nsimplify(x**0.0) == 1
    assert nsimplify(3.333333, tolerance=0.1, rational=True) == Rational(10, 3)
    assert nsimplify(3.333333, tolerance=0.01, rational=True) == Rational(10, 3)
    assert nsimplify(3.666666, tolerance=0.1, rational=True) == Rational(11, 3)
    assert nsimplify(3.666666, tolerance=0.01, rational=True) == Rational(11, 3)
    assert nsimplify(33, tolerance=10, rational=True) == Rational(33)
    assert nsimplify(33.33, tolerance=10, rational=True) == Rational(30)
    assert nsimplify(37.76, tolerance=10, rational=True) == Rational(40)
    assert nsimplify(-203.1) == -S(2031)/10
    assert nsimplify(.2, tolerance=0) == S.One/5
    assert nsimplify(-.2, tolerance=0) == -S.One/5
    assert nsimplify(.2222, tolerance=0) == S(1111)/5000
    assert nsimplify(-.2222, tolerance=0) == -S(1111)/5000
    # issue 7211, PR 4112
    assert nsimplify(S(2e-8)) == S(1)/50000000
    # issue 7322 direct test
    assert nsimplify(1e-42, rational=True) != 0

def test_issue_9448():
    tmp = sympify("1/(1 - (-1)**(2/3) - (-1)**(1/3)) + 1/(1 + (-1)**(2/3) + (-1)**(1/3))")
    assert nsimplify(tmp) == S(1)/2

def test_extract_minus_sign():
    x = Symbol("x")
    y = Symbol("y")
    a = Symbol("a")
    b = Symbol("b")
    assert simplify(-x/-y) == x/y
    assert simplify(-x/y) == -x/y
    assert simplify(x/y) == x/y
    assert simplify(x/-y) == -x/y
    assert simplify(-x/0) == zoo*x
    assert simplify(S(-5)/0) == zoo
    assert simplify(-a*x/(-y - b)) == a*x/(b + y)


def test_diff():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    g = Function("g")
    assert simplify(g(x).diff(x)*f(x).diff(x) - f(x).diff(x)*g(x).diff(x)) == 0
    assert simplify(2*f(x)*f(x).diff(x) - diff(f(x)**2, x)) == 0
    assert simplify(diff(1/f(x), x) + f(x).diff(x)/f(x)**2) == 0
    assert simplify(f(x).diff(x, y) - f(x).diff(y, x)) == 0


def test_logcombine_1():
    x, y = symbols("x,y")
    a = Symbol("a")
    z, w = symbols("z,w", positive=True)
    b = Symbol("b", real=True)
    assert logcombine(log(x) + 2*log(y)) == log(x) + 2*log(y)
    assert logcombine(log(x) + 2*log(y), force=True) == log(x*y**2)
    assert logcombine(a*log(w) + log(z)) == a*log(w) + log(z)
    assert logcombine(b*log(z) + b*log(x)) == log(z**b) + b*log(x)
    assert logcombine(b*log(z) - log(w)) == log(z**b/w)
    assert logcombine(log(x)*log(z)) == log(x)*log(z)
    assert logcombine(log(w)*log(x)) == log(w)*log(x)
    assert logcombine(cos(-2*log(z) + b*log(w))) in [cos(log(w**b/z**2)),
                                                   cos(log(z**2/w**b))]
    assert logcombine(log(log(x) - log(y)) - log(z), force=True) == \
        log(log(x/y)/z)
    assert logcombine((2 + I)*log(x), force=True) == (2 + I)*log(x)
    assert logcombine((x**2 + log(x) - log(y))/(x*y), force=True) == \
        (x**2 + log(x/y))/(x*y)
    # the following could also give log(z*x**log(y**2)), what we
    # are testing is that a canonical result is obtained
    assert logcombine(log(x)*2*log(y) + log(z), force=True) == \
        log(z*y**log(x**2))
    assert logcombine((x*y + sqrt(x**4 + y**4) + log(x) - log(y))/(pi*x**Rational(2, 3)*
            sqrt(y)**3), force=True) == (
            x*y + sqrt(x**4 + y**4) + log(x/y))/(pi*x**(S(2)/3)*y**(S(3)/2))
    assert logcombine(gamma(-log(x/y))*acos(-log(x/y)), force=True) == \
        acos(-log(x/y))*gamma(-log(x/y))

    assert logcombine(2*log(z)*log(w)*log(x) + log(z) + log(w)) == \
        log(z**log(w**2))*log(x) + log(w*z)
    assert logcombine(3*log(w) + 3*log(z)) == log(w**3*z**3)
    assert logcombine(x*(y + 1) + log(2) + log(3)) == x*(y + 1) + log(6)
    assert logcombine((x + y)*log(w) + (-x - y)*log(3)) == (x + y)*log(w/3)


def test_logcombine_complex_coeff():
    i = Integral((sin(x**2) + cos(x**3))/x, x)
    assert logcombine(i, force=True) == i
    assert logcombine(i + 2*log(x), force=True) == \
        i + log(x**2)


def test_posify():
    from sympy.abc import x

    assert str(posify(
        x +
        Symbol('p', positive=True) +
        Symbol('n', negative=True))) == '(_x + n + p, {_x: x})'

    # log(1/x).expand() should be log(1/x) but it comes back as -log(x)
    # when it is corrected, posify will allow the change to be made. The
    # force=True option can do so as well when it is implemented.
    eq, rep = posify(1/x)
    assert log(eq).expand().subs(rep) == -log(x)
    assert str(posify([x, 1 + x])) == '([_x, _x + 1], {_x: x})'

    x = symbols('x')
    p = symbols('p', positive=True)
    n = symbols('n', negative=True)
    orig = [x, n, p]
    modified, reps = posify(orig)
    assert str(modified) == '[_x, n, p]'
    assert [w.subs(reps) for w in modified] == orig

    assert str(Integral(posify(1/x + y)[0], (y, 1, 3)).expand()) == \
        'Integral(1/_x, (y, 1, 3)) + Integral(_y, (y, 1, 3))'
    assert str(Sum(posify(1/x**n)[0], (n,1,3)).expand()) == \
        'Sum(_x**(-n), (n, 1, 3))'


def test_powdenest():
    from sympy import powdenest
    from sympy.abc import x, y, z, a, b
    p, q = symbols('p q', positive=True)
    i, j = symbols('i,j', integer=True)

    assert powdenest(x) == x
    assert powdenest(x + 2*(x**(2*a/3))**(3*x)) == (x + 2*(x**(2*a/3))**(3*x))
    assert powdenest((exp(2*a/3))**(3*x))  # -X-> (exp(a/3))**(6*x)
    assert powdenest((x**(2*a/3))**(3*x)) == ((x**(2*a/3))**(3*x))
    assert powdenest(exp(3*x*log(2))) == 2**(3*x)
    assert powdenest(sqrt(p**2)) == p
    i, j = symbols('i,j', integer=True)
    eq = p**(2*i)*q**(4*i)
    assert powdenest(eq) == (p*q**2)**(2*i)
    # -X-> (x**x)**i*(x**x)**j == x**(x*(i + j))
    assert powdenest((x**x)**(i + j))
    assert powdenest(exp(3*y*log(x))) == x**(3*y)
    assert powdenest(exp(y*(log(a) + log(b)))) == (a*b)**y
    assert powdenest(exp(3*(log(a) + log(b)))) == a**3*b**3
    assert powdenest(((x**(2*i))**(3*y))**x) == ((x**(2*i))**(3*y))**x
    assert powdenest(((x**(2*i))**(3*y))**x, force=True) == x**(6*i*x*y)
    assert powdenest(((x**(2*a/3))**(3*y/i))**x) == \
        (((x**(2*a/3))**(3*y/i))**x)
    assert powdenest((x**(2*i)*y**(4*i))**z, force=True) == (x*y**2)**(2*i*z)
    assert powdenest((p**(2*i)*q**(4*i))**j) == (p*q**2)**(2*i*j)
    e = ((p**(2*a))**(3*y))**x
    assert powdenest(e) == e
    e = ((x**2*y**4)**a)**(x*y)
    assert powdenest(e) == e
    e = (((x**2*y**4)**a)**(x*y))**3
    assert powdenest(e) == ((x**2*y**4)**a)**(3*x*y)
    assert powdenest((((x**2*y**4)**a)**(x*y)), force=True) == \
        (x*y**2)**(2*a*x*y)
    assert powdenest((((x**2*y**4)**a)**(x*y))**3, force=True) == \
        (x*y**2)**(6*a*x*y)
    assert powdenest((x**2*y**6)**i) != (x*y**3)**(2*i)
    x, y = symbols('x,y', positive=True)
    assert powdenest((x**2*y**6)**i) == (x*y**3)**(2*i)

    assert powdenest((x**(2*i/3)*y**(i/2))**(2*i)) == (x**(S(4)/3)*y)**(i**2)
    assert powdenest(sqrt(x**(2*i)*y**(6*i))) == (x*y**3)**i

    assert powdenest(4**x) == 2**(2*x)
    assert powdenest((4**x)**y) == 2**(2*x*y)
    assert powdenest(4**x*y) == 2**(2*x)*y


def test_powdenest_polar():
    x, y, z = symbols('x y z', polar=True)
    a, b, c = symbols('a b c')
    assert powdenest((x*y*z)**a) == x**a*y**a*z**a
    assert powdenest((x**a*y**b)**c) == x**(a*c)*y**(b*c)
    assert powdenest(((x**a)**b*y**c)**c) == x**(a*b*c)*y**(c**2)


def test_issue_5805():
    arg = ((gamma(x)*hyper((), (), x))*pi)**2
    assert powdenest(arg) == (pi*gamma(x)*hyper((), (), x))**2
    assert arg.is_positive is None


def test_issue_4194():
    # simplify should call cancel
    from sympy.abc import x, y
    f = Function('f')
    assert simplify((4*x + 6*f(y))/(2*x + 3*f(y))) == 2


@XFAIL
def test_simplify_float_vs_integer():
    # Test for issue 4473:
    # https://github.com/sympy/sympy/issues/4473
    assert simplify(x**2.0 - x**2) == 0
    assert simplify(x**2 - x**2.0) == 0


def test_combsimp():
    from sympy.abc import n, k

    assert combsimp(factorial(n)) == factorial(n)
    assert combsimp(binomial(n, k)) == binomial(n, k)

    assert combsimp(factorial(n)/factorial(n - 3)) == n*(-1 + n)*(-2 + n)
    assert combsimp(binomial(n + 1, k + 1)/binomial(n, k)) == (1 + n)/(1 + k)

    assert combsimp(binomial(3*n + 4, n + 1)/binomial(3*n + 1, n)) == \
        S(3)/2*((3*n + 2)*(3*n + 4)/((n + 1)*(2*n + 3)))

    assert combsimp(factorial(n)**2/factorial(n - 3)) == \
        factorial(n)*n*(-1 + n)*(-2 + n)
    assert combsimp(factorial(n)*binomial(n + 1, k + 1)/binomial(n, k)) == \
        factorial(n)*(1 + n)/(1 + k)

    assert combsimp(binomial(n - 1, k)) == -((-n + k)*binomial(n, k))/n

    assert combsimp(binomial(n + 2, k + S(1)/2)) == 4*((n + 1)*(n + 2) *
        binomial(n, k + S(1)/2))/((2*k - 2*n - 1)*(2*k - 2*n - 3))
    assert combsimp(binomial(n + 2, k + 2.0)) == \
        -((1.0*n + 2.0)*binomial(n + 1.0, k + 2.0))/(k - n)

    # coverage tests
    assert combsimp(factorial(n*(1 + n) - n**2 - n)) == 1
    assert combsimp(binomial(n + k - 2, n)) == \
        k*(k - 1)*binomial(n + k, n)/((n + k)*(n + k - 1))
    i = Symbol('i', integer=True)
    e = gamma(i + 3)
    assert combsimp(e) == e
    e = gamma(exp(i))
    assert combsimp(e) == e
    e = gamma(n + S(1)/3)*gamma(n + S(2)/3)
    assert combsimp(e) == e
    assert combsimp(gamma(4*n + S(1)/2)/gamma(2*n - S(3)/4)) == \
        2**(4*n - S(5)/2)*(8*n - 3)*gamma(2*n + S(3)/4)/sqrt(pi)

    assert combsimp(6*FallingFactorial(-4, n)/factorial(n)) == \
        (-1)**n*(n + 1)*(n + 2)*(n + 3)
    assert combsimp(6*FallingFactorial(-4, n - 1)/factorial(n - 1)) == \
        (-1)**(n - 1)*n*(n + 1)*(n + 2)
    assert combsimp(6*FallingFactorial(-4, n - 3)/factorial(n - 3)) == \
        (-1)**(n - 3)*n*(n - 1)*(n - 2)
    assert combsimp(6*FallingFactorial(-4, -n - 1)/factorial(-n - 1)) == \
        -(-1)**(-n - 1)*n*(n - 1)*(n - 2)

    assert combsimp(6*RisingFactorial(4, n)/factorial(n)) == \
        (n + 1)*(n + 2)*(n + 3)
    assert combsimp(6*RisingFactorial(4, n - 1)/factorial(n - 1)) == \
        n*(n + 1)*(n + 2)
    assert combsimp(6*RisingFactorial(4, n - 3)/factorial(n - 3)) == \
        n*(n - 1)*(n - 2)
    assert combsimp(6*RisingFactorial(4, -n - 1)/factorial(-n - 1)) == \
        -n*(n - 1)*(n - 2)


def test_issue_5615():
    aA, Re, a, b, D = symbols('aA Re a b D')
    e = ((D**3*a + b*aA**3)/Re).expand()
    assert collect(e, [aA**3/Re, a]) == e


def test_issue_5728():
    b = x*sqrt(y)
    a = sqrt(b)
    c = sqrt(sqrt(x)*y)
    assert powsimp(a*b) == sqrt(b)**3
    assert powsimp(a*b**2*sqrt(y)) == sqrt(y)*a**5
    assert powsimp(a*x**2*c**3*y) == c**3*a**5
    assert powsimp(a*x*c**3*y**2) == c**7*a
    assert powsimp(x*c**3*y**2) == c**7
    assert powsimp(x*c**3*y) == x*y*c**3
    assert powsimp(sqrt(x)*c**3*y) == c**5
    assert powsimp(sqrt(x)*a**3*sqrt(y)) == sqrt(x)*sqrt(y)*a**3
    assert powsimp(Mul(sqrt(x)*c**3*sqrt(y), y, evaluate=False)) == \
        sqrt(x)*sqrt(y)**3*c**3
    assert powsimp(a**2*a*x**2*y) == a**7

    # symbolic powers work, too
    b = x**y*y
    a = b*sqrt(b)
    assert a.is_Mul is True
    assert powsimp(a) == sqrt(b)**3

    # as does exp
    a = x*exp(2*y/3)
    assert powsimp(a*sqrt(a)) == sqrt(a)**3
    assert powsimp(a**2*sqrt(a)) == sqrt(a)**5
    assert powsimp(a**2*sqrt(sqrt(a))) == sqrt(sqrt(a))**9


def test_as_content_primitive():
    # although the _as_content_primitive methods do not alter the underlying structure,
    # the as_content_primitive function will touch up the expression and join
    # bases that would otherwise have not been joined.
    assert ((x*(2 + 2*x)*(3*x + 3)**2)).as_content_primitive() == \
        (18, x*(x + 1)**3)
    assert (2 + 2*x + 2*y*(3 + 3*y)).as_content_primitive() == \
        (2, x + 3*y*(y + 1) + 1)
    assert ((2 + 6*x)**2).as_content_primitive() == \
        (4, (3*x + 1)**2)
    assert ((2 + 6*x)**(2*y)).as_content_primitive() == \
        (1, (_keep_coeff(S(2), (3*x + 1)))**(2*y))
    assert (5 + 10*x + 2*y*(3 + 3*y)).as_content_primitive() == \
        (1, 10*x + 6*y*(y + 1) + 5)
    assert ((5*(x*(1 + y)) + 2*x*(3 + 3*y))).as_content_primitive() == \
        (11, x*(y + 1))
    assert ((5*(x*(1 + y)) + 2*x*(3 + 3*y))**2).as_content_primitive() == \
        (121, x**2*(y + 1)**2)
    assert (y**2).as_content_primitive() == \
        (1, y**2)
    assert (S.Infinity).as_content_primitive() == (1, oo)
    eq = x**(2 + y)
    assert (eq).as_content_primitive() == (1, eq)
    assert (S.Half**(2 + x)).as_content_primitive() == (S(1)/4, 2**-x)
    assert ((-S.Half)**(2 + x)).as_content_primitive() == \
           (S(1)/4, (-S.Half)**x)
    assert ((-S.Half)**(2 + x)).as_content_primitive() == \
           (S(1)/4, (-S.Half)**x)
    assert (4**((1 + y)/2)).as_content_primitive() == (2, 4**(y/2))
    assert (3**((1 + y)/2)).as_content_primitive() == \
           (1, 3**(Mul(S(1)/2, 1 + y, evaluate=False)))
    assert (5**(S(3)/4)).as_content_primitive() == (1, 5**(S(3)/4))
    assert (5**(S(7)/4)).as_content_primitive() == (5, 5**(S(3)/4))
    assert Add(5*z/7, 0.5*x, 3*y/2, evaluate=False).as_content_primitive() == \
              (S(1)/14, 7.0*x + 21*y + 10*z)
    assert (2**(S(3)/4) + 2**(S(1)/4)*sqrt(3)).as_content_primitive(radical=True) == \
           (1, 2**(S(1)/4)*(sqrt(2) + sqrt(3)))


def test_radsimp():
    r2 = sqrt(2)
    r3 = sqrt(3)
    r5 = sqrt(5)
    r7 = sqrt(7)
    assert fraction(radsimp(1/r2)) == (sqrt(2), 2)
    assert radsimp(1/(1 + r2)) == \
        -1 + sqrt(2)
    assert radsimp(1/(r2 + r3)) == \
        -sqrt(2) + sqrt(3)
    assert fraction(radsimp(1/(1 + r2 + r3))) == \
        (-sqrt(6) + sqrt(2) + 2, 4)
    assert fraction(radsimp(1/(r2 + r3 + r5))) == \
        (-sqrt(30) + 2*sqrt(3) + 3*sqrt(2), 12)
    assert fraction(radsimp(1/(1 + r2 + r3 + r5))) == (
        (-34*sqrt(10) - 26*sqrt(15) - 55*sqrt(3) - 61*sqrt(2) + 14*sqrt(30) +
        93 + 46*sqrt(6) + 53*sqrt(5), 71))
    assert fraction(radsimp(1/(r2 + r3 + r5 + r7))) == (
        (-50*sqrt(42) - 133*sqrt(5) - 34*sqrt(70) - 145*sqrt(3) + 22*sqrt(105)
        + 185*sqrt(2) + 62*sqrt(30) + 135*sqrt(7), 215))
    z = radsimp(1/(1 + r2/3 + r3/5 + r5 + r7))
    assert len((3616791619821680643598*z).args) == 16
    assert radsimp(1/z) == 1/z
    assert radsimp(1/z, max_terms=20).expand() == 1 + r2/3 + r3/5 + r5 + r7
    assert radsimp(1/(r2*3)) == \
        sqrt(2)/6
    assert radsimp(1/(r2*a + r3 + r5 + r7)) == (
        (8*sqrt(2)*a**7 - 8*sqrt(7)*a**6 - 8*sqrt(5)*a**6 - 8*sqrt(3)*a**6 -
        180*sqrt(2)*a**5 + 8*sqrt(30)*a**5 + 8*sqrt(42)*a**5 + 8*sqrt(70)*a**5
        - 24*sqrt(105)*a**4 + 84*sqrt(3)*a**4 + 100*sqrt(5)*a**4 +
        116*sqrt(7)*a**4 - 72*sqrt(70)*a**3 - 40*sqrt(42)*a**3 -
        8*sqrt(30)*a**3 + 782*sqrt(2)*a**3 - 462*sqrt(3)*a**2 -
        302*sqrt(7)*a**2 - 254*sqrt(5)*a**2 + 120*sqrt(105)*a**2 -
        795*sqrt(2)*a - 62*sqrt(30)*a + 82*sqrt(42)*a + 98*sqrt(70)*a -
        118*sqrt(105) + 59*sqrt(7) + 295*sqrt(5) + 531*sqrt(3))/(16*a**8 -
        480*a**6 + 3128*a**4 - 6360*a**2 + 3481))
    assert radsimp(1/(r2*a + r2*b + r3 + r7)) == (
        (sqrt(2)*a*(a + b)**2 - 5*sqrt(2)*a + sqrt(42)*a + sqrt(2)*b*(a +
        b)**2 - 5*sqrt(2)*b + sqrt(42)*b - sqrt(7)*(a + b)**2 - sqrt(3)*(a +
        b)**2 - 2*sqrt(3) + 2*sqrt(7))/(2*a**4 + 8*a**3*b + 12*a**2*b**2 -
        20*a**2 + 8*a*b**3 - 40*a*b + 2*b**4 - 20*b**2 + 8))
    assert radsimp(1/(r2*a + r2*b + r2*c + r2*d)) == \
        sqrt(2)/(2*a + 2*b + 2*c + 2*d)
    assert radsimp(1/(1 + r2*a + r2*b + r2*c + r2*d)) == (
        (sqrt(2)*a + sqrt(2)*b + sqrt(2)*c + sqrt(2)*d - 1)/(2*a**2 + 4*a*b +
        4*a*c + 4*a*d + 2*b**2 + 4*b*c + 4*b*d + 2*c**2 + 4*c*d + 2*d**2 - 1))
    assert radsimp((y**2 - x)/(y - sqrt(x))) == \
        sqrt(x) + y
    assert radsimp(-(y**2 - x)/(y - sqrt(x))) == \
        -(sqrt(x) + y)
    assert radsimp(1/(1 - I + a*I)) == \
        (-I*a + 1 + I)/(a**2 - 2*a + 2)
    assert radsimp(1/((-x + y)*(x - sqrt(y)))) == \
        (-x - sqrt(y))/((x - y)*(x**2 - y))
    e = (3 + 3*sqrt(2))*x*(3*x - 3*sqrt(y))
    assert radsimp(e) == x*(3 + 3*sqrt(2))*(3*x - 3*sqrt(y))
    assert radsimp(1/e) == (
        (-9*x + 9*sqrt(2)*x - 9*sqrt(y) + 9*sqrt(2)*sqrt(y))/(9*x*(9*x**2 -
        9*y)))
    assert radsimp(1 + 1/(1 + sqrt(3))) == \
        Mul(S.Half, -1 + sqrt(3), evaluate=False) + 1
    A = symbols("A", commutative=False)
    assert radsimp(x**2 + sqrt(2)*x**2 - sqrt(2)*x*A) == \
        x**2 + sqrt(2)*x**2 - sqrt(2)*x*A
    assert radsimp(1/sqrt(5 + 2 * sqrt(6))) == -sqrt(2) + sqrt(3)
    assert radsimp(1/sqrt(5 + 2 * sqrt(6))**3) == -(-sqrt(3) + sqrt(2))**3

    # issue 6532
    assert fraction(radsimp(1/sqrt(x))) == (sqrt(x), x)
    assert fraction(radsimp(1/sqrt(2*x + 3))) == (sqrt(2*x + 3), 2*x + 3)
    assert fraction(radsimp(1/sqrt(2*(x + 3)))) == (sqrt(2*x + 6), 2*x + 6)

    # issue 5994
    e = S('-(2 + 2*sqrt(2) + 4*2**(1/4))/'
        '(1 + 2**(3/4) + 3*2**(1/4) + 3*sqrt(2))')
    assert radsimp(e).expand() == -2*2**(S(3)/4) - 2*2**(S(1)/4) + 2 + 2*sqrt(2)

    # issue 5986 (modifications to radimp didn't initially recognize this so
    # the test is included here)
    assert radsimp(1/(-sqrt(5)/2 - S(1)/2 + (-sqrt(5)/2 - S(1)/2)**2)) == 1

    # from issue 5934
    eq = (
        (-240*sqrt(2)*sqrt(sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) -
        360*sqrt(2)*sqrt(-8*sqrt(5) + 40)*sqrt(-sqrt(5) + 5) -
        120*sqrt(10)*sqrt(-8*sqrt(5) + 40)*sqrt(-sqrt(5) + 5) +
        120*sqrt(2)*sqrt(-sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) +
        120*sqrt(2)*sqrt(-8*sqrt(5) + 40)*sqrt(sqrt(5) + 5) +
        120*sqrt(10)*sqrt(-sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) +
        120*sqrt(10)*sqrt(-8*sqrt(5) + 40)*sqrt(sqrt(5) + 5))/(-36000 -
        7200*sqrt(5) + (12*sqrt(10)*sqrt(sqrt(5) + 5) +
        24*sqrt(10)*sqrt(-sqrt(5) + 5))**2))
    assert radsimp(eq) is S.NaN  # it's 0/0

    # work with normal form
    e = 1/sqrt(sqrt(7)/7 + 2*sqrt(2) + 3*sqrt(3) + 5*sqrt(5)) + 3
    assert radsimp(e) == (
        -sqrt(sqrt(7) + 14*sqrt(2) + 21*sqrt(3) +
        35*sqrt(5))*(-11654899*sqrt(35) - 1577436*sqrt(210) - 1278438*sqrt(15)
        - 1346996*sqrt(10) + 1635060*sqrt(6) + 5709765 + 7539830*sqrt(14) +
        8291415*sqrt(21))/1300423175 + 3)

    # obey power rules
    base = sqrt(3) - sqrt(2)
    assert radsimp(1/base**3) == (sqrt(3) + sqrt(2))**3
    assert radsimp(1/(-base)**3) == -(sqrt(2) + sqrt(3))**3
    assert radsimp(1/(-base)**x) == (-base)**(-x)
    assert radsimp(1/base**x) == (sqrt(2) + sqrt(3))**x
    assert radsimp(root(1/(-1 - sqrt(2)), -x)) == (-1)**(-1/x)*(1 + sqrt(2))**(1/x)

    # recurse
    e = cos(1/(1 + sqrt(2)))
    assert radsimp(e) == cos(-sqrt(2) + 1)
    assert radsimp(e/2) == cos(-sqrt(2) + 1)/2
    assert radsimp(1/e) == 1/cos(-sqrt(2) + 1)
    assert radsimp(2/e) == 2/cos(-sqrt(2) + 1)
    assert fraction(radsimp(e/sqrt(x))) == (sqrt(x)*cos(-sqrt(2)+1), x)

    # test that symbolic denominators are not processed
    r = 1 + sqrt(2)
    assert radsimp(x/r, symbolic=False) == -x*(-sqrt(2) + 1)
    assert radsimp(x/(y + r), symbolic=False) == x/(y + 1 + sqrt(2))
    assert radsimp(x/(y + r)/r, symbolic=False) == \
        -x*(-sqrt(2) + 1)/(y + 1 + sqrt(2))

    # issue 7408
    eq = sqrt(x)/sqrt(y)
    assert radsimp(eq) == umul(sqrt(x), sqrt(y), 1/y)
    assert radsimp(eq, symbolic=False) == eq

    # issue 7498
    assert radsimp(sqrt(x)/sqrt(y)**3) == umul(sqrt(x), sqrt(y**3), 1/y**3)

    # for coverage
    eq = sqrt(x)/y**2
    assert radsimp(eq) == eq


def test_radsimp_issue_3214():
    c, p = symbols('c p', positive=True)
    s = sqrt(c**2 - p**2)
    b = (c + I*p - s)/(c + I*p + s)
    assert radsimp(b) == -I*(c + I*p - sqrt(c**2 - p**2))**2/(2*c*p)


def test_collect_const():
    # coverage not provided by above tests
    assert collect_const(2*sqrt(3) + 4*a*sqrt(5)) == \
        2*(2*sqrt(5)*a + sqrt(3))  # let the primitive reabsorb
    assert collect_const(2*sqrt(3) + 4*a*sqrt(5), sqrt(3)) == \
        2*sqrt(3) + 4*a*sqrt(5)
    assert collect_const(sqrt(2)*(1 + sqrt(2)) + sqrt(3) + x*sqrt(2)) == \
        sqrt(2)*(x + 1 + sqrt(2)) + sqrt(3)

    # issue 5290
    assert collect_const(2*x + 2*y + 1, 2) == \
        collect_const(2*x + 2*y + 1) == \
        Add(S(1), Mul(2, x + y, evaluate=False), evaluate=False)
    assert collect_const(-y - z) == Mul(-1, y + z, evaluate=False)
    assert collect_const(2*x - 2*y - 2*z, 2) == \
        Mul(2, x - y - z, evaluate=False)
    assert collect_const(2*x - 2*y - 2*z, -2) == \
        _unevaluated_Add(2*x, Mul(-2, y + z, evaluate=False))

    # this is why the content_primitive is used
    eq = (sqrt(15 + 5*sqrt(2))*x + sqrt(3 + sqrt(2))*y)*2
    assert collect_sqrt(eq + 2) == \
        2*sqrt(sqrt(2) + 3)*(sqrt(5)*x + y) + 2


def test_issue_5933():
    from sympy import Polygon, RegularPolygon, denom
    x = Polygon(*RegularPolygon((0, 0), 1, 5).vertices).centroid.x
    assert abs(denom(x).n()) > 1e-12
    assert abs(denom(radsimp(x))) > 1e-12  # in case simplify didn't handle it


def test_fraction_expand():
    eq = (x + y)*y/x
    assert eq.expand(frac=True) == fraction_expand(eq) == (x*y + y**2)/x
    assert eq.expand() == y + y**2/x


def test_combsimp_gamma():
    from sympy.abc import x, y
    R = Rational

    assert combsimp(gamma(x)) == gamma(x)
    assert combsimp(gamma(x + 1)/x) == gamma(x)
    assert combsimp(gamma(x)/(x - 1)) == gamma(x - 1)
    assert combsimp(x*gamma(x)) == gamma(x + 1)
    assert combsimp((x + 1)*gamma(x + 1)) == gamma(x + 2)
    assert combsimp(gamma(x + y)*(x + y)) == gamma(x + y + 1)
    assert combsimp(x/gamma(x + 1)) == 1/gamma(x)
    assert combsimp((x + 1)**2/gamma(x + 2)) == (x + 1)/gamma(x + 1)
    assert combsimp(x*gamma(x) + gamma(x + 3)/(x + 2)) == \
        (x + 2)*gamma(x + 1)

    assert combsimp(gamma(2*x)*x) == gamma(2*x + 1)/2
    assert combsimp(gamma(2*x)/(x - S(1)/2)) == 2*gamma(2*x - 1)

    assert combsimp(gamma(x)*gamma(1 - x)) == pi/sin(pi*x)
    assert combsimp(gamma(x)*gamma(-x)) == -pi/(x*sin(pi*x))
    assert combsimp(1/gamma(x + 3)/gamma(1 - x)) == \
        sin(pi*x)/(pi*x*(x + 1)*(x + 2))

    assert powsimp(combsimp(
        gamma(x)*gamma(x + S(1)/2)*gamma(y)/gamma(x + y))) == \
        2**(-2*x + 1)*sqrt(pi)*gamma(2*x)*gamma(y)/gamma(x + y)
    assert combsimp(1/gamma(x)/gamma(x - S(1)/3)/gamma(x + S(1)/3)) == \
        3**(3*x - S(3)/2)/(2*pi*gamma(3*x - 1))
    assert simplify(
        gamma(S(1)/2 + x/2)*gamma(1 + x/2)/gamma(1 + x)/sqrt(pi)*2**x) == 1
    assert combsimp(gamma(S(-1)/4)*gamma(S(-3)/4)) == 16*sqrt(2)*pi/3

    assert powsimp(combsimp(gamma(2*x)/gamma(x))) == \
        2**(2*x - 1)*gamma(x + S(1)/2)/sqrt(pi)

    # issue 6792
    e = (-gamma(k)*gamma(k + 2) + gamma(k + 1)**2)/gamma(k)**2
    assert combsimp(e) == -k
    assert combsimp(1/e) == -1/k
    e = (gamma(x) + gamma(x + 1))/gamma(x)
    assert combsimp(e) == x + 1
    assert combsimp(1/e) == 1/(x + 1)
    e = (gamma(x) + gamma(x + 2))*(gamma(x - 1) + gamma(x))/gamma(x)
    assert combsimp(e) == (x**2 + x + 1)*gamma(x + 1)/(x - 1)
    e = (-gamma(k)*gamma(k + 2) + gamma(k + 1)**2)/gamma(k)**2
    assert combsimp(e**2) == k**2
    assert combsimp(e**2/gamma(k + 1)) == k/gamma(k)
    a = R(1, 2) + R(1, 3)
    b = a + R(1, 3)
    assert combsimp(gamma(2*k)/gamma(k)*gamma(k + a)*gamma(k + b))
    3*2**(2*k + 1)*3**(-3*k - 2)*sqrt(pi)*gamma(3*k + R(3, 2))/2

    A, B = symbols('A B', commutative=False)
    assert combsimp(e*B*A) == combsimp(e)*B*A

    # check iteration
    assert combsimp(gamma(2*k)/gamma(k)*gamma(-k - R(1, 2))) == (
        -2**(2*k + 1)*sqrt(pi)/(2*((2*k + 1)*cos(pi*k))))
    assert combsimp(
        gamma(k)*gamma(k + R(1, 3))*gamma(k + R(2, 3))/gamma(3*k/2)) == (
        3*2**(3*k + 1)*3**(-3*k - S.Half)*sqrt(pi)*gamma(3*k/2 + S.Half)/2)


def test_polarify():
    from sympy import polar_lift, polarify
    x = Symbol('x')
    z = Symbol('z', polar=True)
    f = Function('f')
    ES = {}

    assert polarify(-1) == (polar_lift(-1), ES)
    assert polarify(1 + I) == (polar_lift(1 + I), ES)

    assert polarify(exp(x), subs=False) == exp(x)
    assert polarify(1 + x, subs=False) == 1 + x
    assert polarify(f(I) + x, subs=False) == f(polar_lift(I)) + x

    assert polarify(x, lift=True) == polar_lift(x)
    assert polarify(z, lift=True) == z
    assert polarify(f(x), lift=True) == f(polar_lift(x))
    assert polarify(1 + x, lift=True) == polar_lift(1 + x)
    assert polarify(1 + f(x), lift=True) == polar_lift(1 + f(polar_lift(x)))

    newex, subs = polarify(f(x) + z)
    assert newex.subs(subs) == f(x) + z

    mu = Symbol("mu")
    sigma = Symbol("sigma", positive=True)

    # Make sure polarify(lift=True) doesn't try to lift the integration
    # variable
    assert polarify(
        Integral(sqrt(2)*x*exp(-(-mu + x)**2/(2*sigma**2))/(2*sqrt(pi)*sigma),
        (x, -oo, oo)), lift=True) == Integral(sqrt(2)*(sigma*exp_polar(0))**exp_polar(I*pi)*
        exp((sigma*exp_polar(0))**(2*exp_polar(I*pi))*exp_polar(I*pi)*polar_lift(-mu + x)**
        (2*exp_polar(0))/2)*exp_polar(0)*polar_lift(x)/(2*sqrt(pi)), (x, -oo, oo))


def test_unpolarify():
    from sympy import (exp_polar, polar_lift, exp, unpolarify,
                       principal_branch)
    from sympy import gamma, erf, sin, tanh, uppergamma, Eq, Ne
    from sympy.abc import x
    p = exp_polar(7*I) + 1
    u = exp(7*I) + 1

    assert unpolarify(1) == 1
    assert unpolarify(p) == u
    assert unpolarify(p**2) == u**2
    assert unpolarify(p**x) == p**x
    assert unpolarify(p*x) == u*x
    assert unpolarify(p + x) == u + x
    assert unpolarify(sqrt(sin(p))) == sqrt(sin(u))

    # Test reduction to principal branch 2*pi.
    t = principal_branch(x, 2*pi)
    assert unpolarify(t) == x
    assert unpolarify(sqrt(t)) == sqrt(t)

    # Test exponents_only.
    assert unpolarify(p**p, exponents_only=True) == p**u
    assert unpolarify(uppergamma(x, p**p)) == uppergamma(x, p**u)

    # Test functions.
    assert unpolarify(sin(p)) == sin(u)
    assert unpolarify(tanh(p)) == tanh(u)
    assert unpolarify(gamma(p)) == gamma(u)
    assert unpolarify(erf(p)) == erf(u)
    assert unpolarify(uppergamma(x, p)) == uppergamma(x, p)

    assert unpolarify(uppergamma(sin(p), sin(p + exp_polar(0)))) == \
        uppergamma(sin(u), sin(u + 1))
    assert unpolarify(uppergamma(polar_lift(0), 2*exp_polar(0))) == \
        uppergamma(0, 2)

    assert unpolarify(Eq(p, 0)) == Eq(u, 0)
    assert unpolarify(Ne(p, 0)) == Ne(u, 0)
    assert unpolarify(polar_lift(x) > 0) == (x > 0)

    # Test bools
    assert unpolarify(True) is True


def test_issue_6097():
    assert collect(a*y**(2.0*x) + b*y**(2.0*x), y**x) == y**(2.0*x)*(a + b)
    assert collect(a*2**(2.0*x) + b*2**(2.0*x), 2**x) == 2**(2.0*x)*(a + b)


def test_signsimp():
    e = x*(-x + 1) + x*(x - 1)
    assert signsimp(Eq(e, 0)) is S.true
    assert Abs(x - 1) == Abs(1 - x)


def test_besselsimp():
    from sympy import besselj, besseli, exp_polar, cosh, cosine_transform
    assert besselsimp(exp(-I*pi*y/2)*besseli(y, z*exp_polar(I*pi/2))) == \
        besselj(y, z)
    assert besselsimp(exp(-I*pi*a/2)*besseli(a, 2*sqrt(x)*exp_polar(I*pi/2))) == \
        besselj(a, 2*sqrt(x))
    assert besselsimp(sqrt(2)*sqrt(pi)*x**(S(1)/4)*exp(I*pi/4)*exp(-I*pi*a/2) *
                      besseli(-S(1)/2, sqrt(x)*exp_polar(I*pi/2)) *
                      besseli(a, sqrt(x)*exp_polar(I*pi/2))/2) == \
        besselj(a, sqrt(x)) * cos(sqrt(x))
    assert besselsimp(besseli(S(-1)/2, z)) == \
        sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(a, z*exp_polar(-I*pi/2))) == \
        exp(-I*pi*a/2)*besselj(a, z)
    assert cosine_transform(1/t*sin(a/t), t, y) == \
        sqrt(2)*sqrt(pi)*besselj(0, 2*sqrt(a)*sqrt(y))/2


def test_Piecewise():
    e1 = x*(x + y) - y*(x + y)
    e2 = sin(x)**2 + cos(x)**2
    e3 = expand((x + y)*y/x)
    s1 = simplify(e1)
    s2 = simplify(e2)
    s3 = simplify(e3)
    assert simplify(Piecewise((e1, x < e2), (e3, True))) == \
        Piecewise((s1, x < s2), (s3, True))

    # trigsimp tries not to touch non-trig containing args
    assert trigsimp(Piecewise((e1, e3 < e2), (e3, True))) == \
        Piecewise((e1, e3 < s2), (e3, True))


def test_polymorphism():
    class A(Basic):
        def _eval_simplify(x, **kwargs):
            return 1

    a = A(5, 2)
    assert simplify(a) == 1


def test_issue_from_PR1599():
    n1, n2, n3, n4 = symbols('n1 n2 n3 n4', negative=True)
    assert simplify(I*sqrt(n1)) == -sqrt(-n1)
    assert (powsimp(sqrt(n1)*sqrt(n2)*sqrt(n3)) ==
        -I*sqrt(-n1)*sqrt(-n2)*sqrt(-n3))
    assert (powsimp(root(n1, 3)*root(n2, 3)*root(n3, 3)*root(n4, 3)) ==
        -(-1)**(S(1)/3)*
        (-n1)**(S(1)/3)*(-n2)**(S(1)/3)*(-n3)**(S(1)/3)*(-n4)**(S(1)/3))


def test_issue_6811():
    eq = (x + 2*y)*(2*x + 2)
    assert simplify(eq) == (x + 1)*(x + 2*y)*2
    # reject the 2-arg Mul -- these are a headache for test writing
    assert simplify(eq.expand()) == \
        2*x**2 + 4*x*y + 2*x + 4*y


@XFAIL
def test_issue_6811_fail():
    # from doc/src/modules/physics/mechanics/examples.rst, the current `eq`
    # at Line 576 (in different variables) was formerly the equivalent and
    # shorter expression given below...it would be nice to get the short one
    # back again
    xp, y, x, z = symbols('xp, y, x, z')
    eq = 4*(-19*sin(x)*y + 5*sin(3*x)*y + 15*cos(2*x)*z - 21*z)*xp/(9*cos(x) - 5*cos(3*x))
    assert trigsimp(eq) == -2*(2*cos(x)*tan(x)*y + 3*z)*xp/cos(x)


def test_issue_6920():
    e = [cos(x) + I*sin(x), cos(x) - I*sin(x),
        cosh(x) - sinh(x), cosh(x) + sinh(x)]
    ok = [exp(I*x), exp(-I*x), exp(-x), exp(x)]
    # wrap in f to show that the change happens wherever ei occurs
    f = Function('f')
    assert [simplify(f(ei)).args[0] for ei in e] == ok


def test_issue_7001():
    from sympy.abc import r, R
    assert simplify(-(r*Piecewise((4*pi/3, r <= R),
        (-8*pi*R**3/(3*r**3), True)) + 2*Piecewise((4*pi*r/3, r <= R),
        (4*pi*R**3/(3*r**2), True)))/(4*pi*r)) == \
        Piecewise((-1, r <= R), (0, True))


def test_exptrigsimp():
    def valid(a, b):
        from sympy.utilities.randtest import verify_numerically as tn
        if not (tn(a, b) and a == b):
            return False
        return True

    assert exptrigsimp(exp(x) + exp(-x)) == 2*cosh(x)
    assert exptrigsimp(exp(x) - exp(-x)) == 2*sinh(x)
    e = [cos(x) + I*sin(x), cos(x) - I*sin(x),
         cosh(x) - sinh(x), cosh(x) + sinh(x)]
    ok = [exp(I*x), exp(-I*x), exp(-x), exp(x)]
    assert all(valid(i, j) for i, j in zip(
        [exptrigsimp(ei) for ei in e], ok))

    ue = [cos(x) + sin(x), cos(x) - sin(x),
          cosh(x) + I*sinh(x), cosh(x) - I*sinh(x)]
    assert [exptrigsimp(ei) == ei for ei in ue]

    res = []
    ok = [y*tanh(1), 1/(y*tanh(1)), I*y*tan(1), -I/(y*tan(1)),
        y*tanh(x), 1/(y*tanh(x)), I*y*tan(x), -I/(y*tan(x)),
        y*tanh(1 + I), 1/(y*tanh(1 + I))]
    for a in (1, I, x, I*x, 1 + I):
        w = exp(a)
        eq = y*(w - 1/w)/(w + 1/w)
        s = simplify(eq)
        assert s == exptrigsimp(eq)
        res.append(s)
        sinv = simplify(1/eq)
        assert sinv == exptrigsimp(1/eq)
        res.append(sinv)
    assert all(valid(i, j) for i, j in zip(res, ok))

    for a in range(1, 3):
        w = exp(a)
        e = w + 1/w
        s = simplify(e)
        assert s == exptrigsimp(e)
        assert valid(s, 2*cosh(a))
        e = w - 1/w
        s = simplify(e)
        assert s == exptrigsimp(e)
        assert valid(s, 2*sinh(a))


def test_issue_2827_trigsimp_methods():
    measure1 = lambda expr: len(str(expr))
    measure2 = lambda expr: -count_ops(expr)
                                       # Return the most complicated result
    expr = (x + 1)/(x + sin(x)**2 + cos(x)**2)
    ans = Matrix([1])
    M = Matrix([expr])
    assert trigsimp(M, method='fu', measure=measure1) == ans
    assert trigsimp(M, method='fu', measure=measure2) != ans
    # all methods should work with Basic expressions even if they
    # aren't Expr
    M = Matrix.eye(1)
    assert all(trigsimp(M, method=m) == M for m in
        'fu matching groebner old'.split())
    # watch for E in exptrigsimp, not only exp()
    eq = 1/sqrt(E) + E
    assert exptrigsimp(eq) == eq


def test_powsimp_on_numbers():
    assert 2**(S(1)/3 - 2) == 2**(S(1)/3)/4


def test_inequality_no_auto_simplify():
    # no simplify on creation but can be simplified
    lhs = cos(x)**2 + sin(x)**2
    rhs = 2;
    e = Lt(lhs, rhs)
    assert e == Lt(lhs, rhs, evaluate=False)
    assert simplify(e)

def test_issue_9398():
    from sympy import Number, cancel
    assert cancel(1e-14) != 0
    assert cancel(1e-14*I) != 0

    assert simplify(1e-14) != 0
    assert simplify(1e-14*I) != 0

    assert (I*Number(1.)*Number(10)**Number(-14)).simplify() != 0

    assert cancel(1e-20) != 0
    assert cancel(1e-20*I) != 0

    assert simplify(1e-20) != 0
    assert simplify(1e-20*I) != 0

    assert cancel(1e-100) != 0
    assert cancel(1e-100*I) != 0

    assert simplify(1e-100) != 0
    assert simplify(1e-100*I) != 0

    f = Float("1e-1000")
    assert cancel(f) != 0
    assert cancel(f*I) != 0

    assert simplify(f) != 0
    assert simplify(f*I) != 0

def test_issue_9324_simplify():
    M = MatrixSymbol('M', 10, 10)
    e = M[0, 0] + M[5, 4] + 1304
    assert simplify(e) == e


def test_issue_9324_powsimp_on_matrix_symbol():
    M = MatrixSymbol('M', 10, 10)
    expr = powsimp(M, deep=True)
    assert expr == M
    assert expr.args[0] == 'M'
