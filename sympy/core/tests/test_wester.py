""" Tests from Michael Wester's 1999 paper "Review of CAS mathematical
capabilities".

http://www.math.unm.edu/~wester/cas/book/Wester.pdf
See also http://math.unm.edu/~wester/cas_review.html for detailed output of each
tested system.
"""

from sympy import (Rational, symbols, factorial, sqrt, log, exp, oo, product,
    binomial, rf, pi, gamma, igcd, factorint, radsimp, combsimp,
    npartitions, totient, primerange, factor, simplify, gcd, resultant, expand,
    I, trigsimp, tan, sin, cos, cot, diff, nan, limit, EulerGamma, polygamma,
    bernoulli, hyper, hyperexpand, besselj, asin, assoc_legendre, Function, re,
    im, DiracDelta, chebyshevt, atan, sinh, cosh, floor, ceiling, solve, asinh,
    LambertW, N, apart, sqrtdenest, factorial2, powdenest, Mul, S, mpmath, ZZ,
    Poly, expand_func, E, Q, And, Ne, Or, Le, Lt, Ge, Gt, QQ, ask, refine, AlgebraicNumber,
    elliptic_e, elliptic_f, powsimp, hessian, wronskian, fibonacci)

from sympy.functions.combinatorial.numbers import stirling
from sympy.functions.special.zeta_functions import zeta
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.utilities.pytest import XFAIL, slow
from sympy.utilities.iterables import partitions
from sympy.mpmath import mpi, mpc
from sympy.matrices import Matrix, GramSchmidt, eye
from sympy.matrices.expressions.blockmatrix import (BlockMatrix, block_collapse)
from sympy.matrices.expressions import (MatrixSymbol, ZeroMatrix)
from sympy.galgebra.ga import MV
from sympy.physics.quantum import Commutator
from sympy.assumptions import assuming
from sympy.polys.rings import vring
from sympy.polys.fields import vfield
from sympy.polys.solvers import solve_lin_sys
from sympy.concrete import Sum
from sympy.concrete.products import Product

R = Rational
x, y, z = symbols('x y z')
i, j, k, l, m, n = symbols('i j k l m n', integer=True)
f = Function('f')
g = Function('g')

# A. Boolean Logic and Quantifier Elimination
#   Not implemented.

# B. Set Theory
#   Not implemented.

# C. Numbers


def test_C1():
    assert (factorial(50) ==
        30414093201713378043612608166064768844377641568960512000000000000)


def test_C2():
    assert (factorint(factorial(50)) == {2: 47, 3: 22, 5: 12, 7: 8,
        11: 4, 13: 3, 17: 2, 19: 2, 23: 2, 29: 1, 31: 1, 37: 1,
        41: 1, 43: 1, 47: 1})


def test_C3():
    assert (factorial2(10), factorial2(9)) == (3840, 945)


# Base conversions; not really implemented by sympy
# Whatever. Take credit!
def test_C4():
    assert 0xABC == 2748


def test_C5():
    assert 123 == int('234', 7)


def test_C6():
    assert int('677', 8) == int('1BF', 16) == 447


def test_C7():
    assert log(32768, 8) == 5


def test_C8():
    # Modular multiplicative inverse. Would be nice if divmod could do this.
    assert ZZ.invert(5, 7) == 3
    assert ZZ.invert(5, 6) == 5


def test_C9():
    assert igcd(igcd(1776, 1554), 5698) == 74


def test_C10():
    x = 0
    for n in range(2, 11):
        x += R(1, n)
    assert x == R(4861, 2520)


def test_C11():
    assert R(1, 7) == S('0.[142857]')


def test_C12():
    assert R(7, 11) * R(22, 7) == 2


def test_C13():
    test = R(10, 7) * (1 + R(29, 1000)) ** R(1, 3)
    good = 3 ** R(1, 3)
    assert test == good


def test_C14():
    assert sqrtdenest(sqrt(2*sqrt(3) + 4)) == 1 + sqrt(3)


def test_C15():
    test = sqrtdenest(sqrt(14 + 3*sqrt(3 + 2*sqrt(5 - 12*sqrt(3 - 2*sqrt(2))))))
    good = sqrt(2) + 3
    assert test == good


def test_C16():
    test = sqrtdenest(sqrt(10 + 2*sqrt(6) + 2*sqrt(10) + 2*sqrt(15)))
    good = sqrt(2) + sqrt(3) + sqrt(5)
    assert test == good


def test_C17():
    test = radsimp((sqrt(3) + sqrt(2)) / (sqrt(3) - sqrt(2)))
    good = 5 + 2*sqrt(6)
    assert test == good


def test_C18():
    assert simplify((sqrt(-2 + sqrt(-5)) * sqrt(-2 - sqrt(-5))).expand(complex=True)) == 3


@XFAIL
def test_C19():
    assert radsimp(simplify((90 + 35*sqrt(7)) ** R(1, 3))) == 3 + sqrt(7)


@XFAIL
def test_C20():
    inside = (135 + 78*sqrt(3))
    test = simplify((inside**R(2, 3) + 3) * sqrt(3) / inside**R(1, 3))
    assert test == 12


@XFAIL
def test_C21():
    assert simplify((41 + 29*sqrt(2)) ** R(1, 5)) == 1 + sqrt(2)


@XFAIL
def test_C22():
    test = simplify(((6 - 4*sqrt(2))*log(3 - 2*sqrt(2)) + (3 - 2*sqrt(2))*log(17
        - 12*sqrt(2)) + 32 - 24*sqrt(2)) / (48*sqrt(2) - 72))
    good = sqrt(2)/3 - log(sqrt(2) - 1)/3
    assert test == good


def test_C23():
    assert 2 * oo - 3 == oo


@XFAIL
def test_C24():
    raise NotImplementedError("2**aleph_null == aleph_1")

# D. Numerical Analysis


def test_D1():
    assert 0.0 / sqrt(2) == 0.0


def test_D2():
    assert str(exp(-1000000).evalf()) == '3.29683147808856e-434295'


def test_D3():
    assert exp(pi*sqrt(163)).evalf(50).num.ae(262537412640768744)


def test_D4():
    assert floor(R(-5, 3)) == -2
    assert ceiling(R(-5, 3)) == -1


@XFAIL
def test_D5():
    raise NotImplementedError("cubic_spline([1, 2, 4, 5], [1, 4, 2, 3], x)(3) == 27/8")


@XFAIL
def test_D6():
    raise NotImplementedError("translate sum(a[i]*x**i, (i,1,n)) to FORTRAN")


@XFAIL
def test_D7():
    raise NotImplementedError("translate sum(a[i]*x**i, (i,1,n)) to C")


@XFAIL
def test_D8():
    # One way is to cheat by converting the sum to a string,
    # and replacing the '[' and ']' with ''.
    # E.g., horner(S(str(_).replace('[','').replace(']','')))
    raise NotImplementedError("apply Horner's rule to sum(a[i]*x**i, (i,1,5))")


@XFAIL
def test_D9():
    raise NotImplementedError("translate D8 to FORTRAN")


@XFAIL
def test_D10():
    raise NotImplementedError("translate D8 to C")


@XFAIL
def test_D11():
    #Is there a way to use count_ops?
    raise NotImplementedError("flops(sum(product(f[i][k], (i,1,k)), (k,1,n)))")


@XFAIL
def test_D12():
    assert (mpi(-4, 2) * x + mpi(1, 3)) ** 2 == mpi(-8, 16)*x**2 + mpi(-24, 12)*x + mpi(1, 9)


@XFAIL
def test_D13():
    raise NotImplementedError("discretize a PDE: diff(f(x,t),t) == diff(diff(f(x,t),x),x)")

# E. Statistics
#   See scipy; all of this is numerical.

# F. Combinatorial Theory.


def test_F1():
    assert rf(x, 3) == x*(1 + x)*(2 + x)


def test_F2():
    assert expand_func(binomial(n, 3)) == n*(n - 1)*(n - 2)/6


@XFAIL
def test_F3():
    assert combsimp(2**n * factorial(n) * factorial2(2*n - 1)) == factorial(2*n)


@XFAIL
def test_F4():
    assert combsimp((2**n * factorial(n) * product(2*k - 1, (k, 1, n)))) == factorial(2*n)


@XFAIL
def test_F5():
    assert gamma(n + R(1, 2)) / sqrt(pi) / factorial(n) == factorial(2*n)/2**(2*n)/factorial(n)**2


def test_F6():
    partTest = [p.copy() for p in partitions(4)]
    partDesired = [{4: 1}, {1: 1, 3: 1}, {2: 2}, {1: 2, 2:1}, {1: 4}]
    assert partTest == partDesired


def test_F7():
    assert npartitions(4) == 5


def test_F8():
    assert stirling(5, 2, signed=True) == -50  # if signed, then kind=1


def test_F9():
    assert totient(1776) == 576

# G. Number Theory


def test_G1():
    assert list(primerange(999983, 1000004)) == [999983, 1000003]


@XFAIL
def test_G2():
    raise NotImplementedError("find the primitive root of 191 == 19")


@XFAIL
def test_G3():
    raise NotImplementedError("(a+b)**p mod p == a**p + b**p mod p; p prime")

# ... G20 Modular equations and continued fractions are not implemented.

# H. Algebra


def test_H1():
    assert simplify(2*2**n) == simplify(2**(n + 1))
    assert powdenest(2*2**n) == simplify(2**(n + 1))


def test_H2():
    assert powsimp(4 * 2**n) == 2**(n + 2)


def test_H3():
    assert (-1)**(n*(n + 1)) == 1


def test_H4():
    expr = factor(6*x - 10)
    assert type(expr) is Mul
    assert expr.args[0] == 2
    assert expr.args[1] == 3*x - 5

p1 = 64*x**34 - 21*x**47 - 126*x**8 - 46*x**5 - 16*x**60 - 81
p2 = 72*x**60 - 25*x**25 - 19*x**23 - 22*x**39 - 83*x**52 + 54*x**10 + 81
q = 34*x**19 - 25*x**16 + 70*x**7 + 20*x**3 - 91*x - 86


def test_H5():
    assert gcd(p1, p2, x) == 1


def test_H6():
    assert gcd(expand(p1 * q), expand(p2 * q)) == q


def test_H7():
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    assert gcd(p1, p2, x, y, z) == 1


def test_H8():
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    q = 11*x**12*y**7*z**13 - 23*x**2*y**8*z**10 + 47*x**17*y**5*z**8
    assert gcd(p1 * q, p2 * q, x, y, z) == q


def test_H9():
    p1 = 2*x**(n + 4) - x**(n + 2)
    p2 = 4*x**(n + 1) + 3*x**n
    assert gcd(p1, p2) == x**n


def test_H10():
    p1 = 3*x**4 + 3*x**3 + x**2 - x - 2
    p2 = x**3 - 3*x**2 + x + 5
    assert resultant(p1, p2, x) == 0


def test_H11():
    assert resultant(p1 * q, p2 * q, x) == 0


def test_H12():
    num = x**2 - 4
    den = x**2 + 4*x + 4
    assert simplify(num/den) == (x - 2)/(x + 2)


@XFAIL
def test_H13():
    assert simplify((exp(x) - 1) / (exp(x/2) + 1)) == exp(x/2) - 1


def test_H14():
    p = (x + 1) ** 20
    ep = expand(p)
    assert ep == (1 + 20*x + 190*x**2 + 1140*x**3 + 4845*x**4 + 15504*x**5
        + 38760*x**6 + 77520*x**7 + 125970*x**8 + 167960*x**9 + 184756*x**10
        + 167960*x**11 + 125970*x**12 + 77520*x**13 + 38760*x**14 + 15504*x**15
        + 4845*x**16 + 1140*x**17 + 190*x**18 + 20*x**19 + x**20)
    dep = diff(ep, x)
    assert dep == (20 + 380*x + 3420*x**2 + 19380*x**3 + 77520*x**4
        + 232560*x**5 + 542640*x**6 + 1007760*x**7 + 1511640*x**8 + 1847560*x**9
        + 1847560*x**10 + 1511640*x**11 + 1007760*x**12 + 542640*x**13
        + 232560*x**14 + 77520*x**15 + 19380*x**16 + 3420*x**17 + 380*x**18
        + 20*x**19)
    assert factor(dep) == 20*(1 + x)**19


def test_H15():
    assert simplify((Mul(*[x - r for r in solve(x**3 + x**2 - 7)]))) == x**3 + x**2 - 7


def test_H16():
    assert factor(x**100 - 1) == ((x - 1)*(x + 1)*(x**2 + 1)*(x**4 - x**3
        + x**2 - x + 1)*(x**4 + x**3 + x**2 + x + 1)*(x**8 - x**6 + x**4
        - x**2 + 1)*(x**20 - x**15 + x**10 - x**5 + 1)*(x**20 + x**15 + x**10
        + x**5 + 1)*(x**40 - x**30 + x**20 - x**10 + 1))


@slow
def test_H17():
    assert simplify(factor(expand(p1 * p2)) - p1*p2) == 0


@XFAIL
def test_H18():
    # Factor over complex rationals.
    test = factor(4*x**4 + 8*x**3 + 77*x**2 + 18*x + 53)
    good = (2*x + 3*I)*(2*x - 3*I)*(x + 1 - 4*I)(x + 1 + 4*I)
    assert test == good


def test_H19():
    a = symbols('a')
    # The idea is to let a**2 == 2, then solve 1/(a-1). Answer is a+1")
    assert Poly(a - 1).invert(Poly(a**2 - 2)) == a + 1


@XFAIL
def test_H20():
    raise NotImplementedError("let a**2==2; (x**3 + (a-2)*x**2 - "
        + "(2*a+3)*x - 3*a) / (x**2-2) = (x**2 - 2*x - 3) / (x-a)")


@XFAIL
def test_H21():
    raise NotImplementedError("evaluate (b+c)**4 assuming b**3==2, c**2==3. \
                              Answer is 2*b + 8*c + 18*b**2 + 12*b*c + 9")


def test_H22():
    assert factor(x**4 - 3*x**2 + 1, modulus=5) == (x - 2)**2 * (x + 2)**2


def test_H23():
    f = x**11 + x + 1
    g = (x**2 + x + 1) * (x**9 - x**8 + x**6 - x**5 + x**3 - x**2 + 1)
    assert factor(f, modulus=65537) == g


def test_H24():
    phi = AlgebraicNumber(S.GoldenRatio.expand(func=True), alias='phi')
    assert factor(x**4 - 3*x**2 + 1, extension=phi) == \
        (x - phi)*(x + 1 - phi)*(x - 1 + phi)*(x + phi)


@slow
def test_H25():
    e = (x - 2*y**2 + 3*z**3) ** 20
    assert factor(expand(e)) == e


@slow
def test_H26():
    g = expand((sin(x) - 2*cos(y)**2 + 3*tan(z)**3)**20)
    assert factor(g, expand=False) == (-sin(x) + 2*cos(y)**2 - 3*tan(z)**3)**20


@slow
def test_H27():
    f = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    g = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    h = -2*z*y**7 \
        *(6*x**9*y**9*z**3 + 10*x**7*z**6 + 17*y*x**5*z**12 + 40*y**7) \
        *(3*x**22 + 47*x**17*y**5*z**8 - 6*x**15*y**9*z**2 - 24*x*y**19*z**8 - 5)
    assert factor(expand(f*g)) == h


@XFAIL
def test_H28():
    raise NotImplementedError("expand ((1 - c**2)**5 * (1 - s**2)**5 * "
        + "(c**2 + s**2)**10) with c**2 + s**2 = 1. Answer is c**10*s**10.")


@XFAIL
def test_H29():
    assert factor(4*x**2 - 21*x*y + 20*y**2, modulus=3) == (x + y)*(x - y)


def test_H30():
    test = factor(x**3 + y**3, extension=sqrt(-3))
    answer = (x + y)*(x + y*(-R(1, 2) - sqrt(3)/2*I))*(x + y*(-R(1, 2) + sqrt(3)/2*I))
    assert answer == test


def test_H31():
    f = (x**2 + 2*x + 3)/(x**3 + 4*x**2 + 5*x + 2)
    g = 2 / (x + 1)**2 - 2 / (x + 1) + 3 / (x + 2)
    assert apart(f) == g


@XFAIL
def test_H32():  # issue 3459
    raise NotImplementedError("[A*B*C - (A*B*C)**(-1)]*A*C*B (product \
                              of a non-commuting product and its inverse)")


def test_H33():
    A, B, C = symbols('A, B, C', commutatative=False)
    assert (Commutator(A, Commutator(B, C))
        + Commutator(B, Commutator(C, A))
        + Commutator(C, Commutator(A, B))).doit().expand() == 0


# I. Trigonometry

@XFAIL
def test_I1():
    assert tan(7*pi/10) == -sqrt(1 + 2/sqrt(5))


@XFAIL
def test_I2():
    assert sqrt((1 + cos(6))/2) == -cos(3)


def test_I3():
    assert cos(n*pi) + sin((4*n - 1)*pi/2) == (-1)**n - 1


def test_I4():
    assert refine(cos(pi*cos(n*pi)) + sin(pi/2*cos(n*pi)), Q.integer(n)) == (-1)**n - 1


@XFAIL
def test_I5():
    assert sin((n**5/5 + n**4/2 + n**3/3 - n/30) * pi) == 0


@XFAIL
def test_I6():
    raise NotImplementedError("assuming -3*pi<x<-5*pi/2, abs(cos(x)) == -cos(x), abs(sin(x)) == -sin(x)")


@XFAIL
def test_I7():
    assert cos(3*x)/cos(x) == cos(x)**2 - 3*sin(x)**2


@XFAIL
def test_I8():
    assert cos(3*x)/cos(x) == 2*cos(2*x) - 1


@XFAIL
def test_I9():
    # Supposed to do this with rewrite rules.
    assert cos(3*x)/cos(x) == cos(x)**2 - 3*sin(x)**2


def test_I10():
    assert trigsimp((tan(x)**2 + 1 - cos(x)**-2) / (sin(x)**2 + cos(x)**2 - 1)) == nan

#@XFAIL
#def test_I11():
#    assert limit((tan(x)**2 + 1 - cos(x)**-2) / (sin(x)**2 + cos(x)**2 - 1), x, 0) != 0


@XFAIL
def test_I12():
    try:
        # This should fail or return nan or something.
        diff((tan(x)**2 + 1 - cos(x)**-2) / (sin(x)**2 + cos(x)**2 - 1), x)
    except:
        assert True
    else:
        assert False, "taking the derivative with a fraction equivalent to 0/0 should fail"

# J. Special functions.


def test_J1():
    assert bernoulli(16) == R(-3617, 510)


def test_J2():
    assert diff(elliptic_e(x, y**2), y) == (elliptic_e(x, y**2) - elliptic_f(x, y**2))/y


@XFAIL
def test_J3():
    raise NotImplementedError("Jacobi elliptic functions: diff(dn(u,k), u) == -k**2*sn(u,k)*cn(u,k)")


def test_J4():
    assert gamma(R(-1, 2)) == -2*sqrt(pi)


def test_J5():
    assert polygamma(0, R(1, 3)) == -EulerGamma - pi/2*sqrt(R(1, 3)) - R(3, 2)*log(3)


def test_J6():
    assert mpmath.besselj(2, 1 + 1j).ae(mpc('0.04157988694396212', '0.24739764151330632'))


def test_J7():
    assert simplify(besselj(R(-5,2), pi/2)) == 12/(pi**2)


def test_J8():
    p = besselj(R(3,2), z)
    q = (sin(z)/z - cos(z))/sqrt(pi*z/2)
    assert simplify(expand_func(p) -q) == 0


def test_J9():
    assert besselj(0, z).diff(z) == - besselj(1, z)


def test_J10():
    mu, nu = symbols('mu, nu', integer=True)
    assert assoc_legendre(nu, mu, 0) == 2**mu*sqrt(pi)/gamma((nu - mu)/2 + 1)/gamma((-nu - mu + 1)/2)


def test_J11():
    assert simplify(assoc_legendre(3, 1, x)) == simplify(-R(3, 2)*sqrt(1 - x**2)*(5*x**2 - 1))


@slow
def test_J12():
    assert simplify(chebyshevt(1008, x) - 2*x*chebyshevt(1007, x) + chebyshevt(1006, x)) == 0


def test_J13():
    a = symbols('a', integer=True, negative=False)
    assert chebyshevt(a, -1) == (-1)**a


def test_J14():
    p = hyper([S(1)/2, S(1)/2], [S(3)/2], z**2)
    assert hyperexpand(p) == asin(z)/z


@XFAIL
def test_J15():
    raise NotImplementedError("F((n+2)/2,-(n-2)/2,R(3,2),sin(z)**2) == sin(n*z)/(n*sin(z)*cos(z)); F(.) is hypergeometric function")


@XFAIL
def test_J16():
    raise NotImplementedError("diff(zeta(x), x) @ x=0 == -log(2*pi)/2")


@XFAIL
def test_J17():
    assert deltaintegrate(f((x + 2)/5)*DiracDelta((x - 2)/3) - g(x)*diff(DiracDelta(x - 1), x), (x, 0, 3))


@XFAIL
def test_J18():
    raise NotImplementedError("define an antisymmetric function")


# K. The Complex Domain

def test_K1():
    z1, z2 = symbols('z1, z2', complex=True)
    assert re(z1 + I*z2) == -im(z2) + re(z1)
    assert im(z1 + I*z2) == im(z1) + re(z2)


@XFAIL  # abs(...).n() does evaluate to 1.00000...
def test_K2():
    assert abs(3 - sqrt(7) + I*sqrt(6*sqrt(7) - 15)) == 1


@XFAIL
def test_K3():
    a, b = symbols('a, b', real=True)
    assert simplify(abs(1/(a + I/a + I*b))) == 1/sqrt(a**2 + (I/a + b)**2)


def test_K4():
    assert log(3 + 4*I).expand(complex=True) == log(5) + I*atan(R(4, 3))


def test_K5():
    x, y = symbols('x, y', real=True)
    assert tan(x + I*y).expand(complex=True) == sin(x)*cos(x) / (cos(x)**2 +
    sinh(y)**2) + I*sinh(y)*cosh(y) / (cos(x)**2 + sinh(y)**2)


def test_K6():
    assert sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z)) == sqrt(x*y)/sqrt(x)
    assert sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z)) != sqrt(y)


def test_K7():
    y = symbols('y', real=True, negative=False)
    expr = sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z))
    sexpr = simplify(expr)
    assert sexpr == sqrt(y)


@XFAIL
def test_K8():
    z = symbols('z', complex=True)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) != 0  # Passes
    z = symbols('z', complex=True, negative=False)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0  # Fails


def test_K9():
    z = symbols('z', real=True, positive=True)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0


def test_K10():
    z = symbols('z', real=True, negative=True)
    assert simplify(sqrt(1/z) + 1/sqrt(z)) == 0

# This goes up to K25

# L. Determining Zero Equivalence


def test_L1():
    assert sqrt(997) - (997**3)**R(1, 6) == 0


def test_L2():
    assert sqrt(999983) - (999983**3)**R(1, 6) == 0


def test_L3():
    assert simplify((2**R(1, 3) + 4**R(1, 3))**3 - 6*(2**R(1, 3) + 4**R(1, 3)) - 6) == 0


def test_L4():
    assert trigsimp(cos(x)**3 + cos(x)*sin(x)**2 - cos(x)) == 0


@XFAIL
def test_L5():
    assert log(tan(R(1, 2)*x + pi/4)) - asinh(tan(x)) == 0


def test_L6():
    assert (log(tan(x/2 + pi/4)) - asinh(tan(x))).diff(x).subs({x: 0}) == 0


@XFAIL
def test_L7():
    assert simplify(log((2*sqrt(x) + 1)/(sqrt(4*x + 4*sqrt(x) + 1)))) == 0


@XFAIL
def test_L8():
    assert simplify((4*x + 4*sqrt(x) + 1)**(sqrt(x)/(2*sqrt(x) + 1)) \
        *(2*sqrt(x) + 1)**(1/(2*sqrt(x) + 1)) - 2*sqrt(x) - 1) == 0


@XFAIL
def test_L9():
    z = symbols('z', complex=True)
    assert simplify(2**(1 - z)*gamma(z)*zeta(z)*cos(z*pi/2) - pi**2*zeta(1 - z)) == 0

# M. Equations


@XFAIL
def test_M1():
    assert Equality(x, 2)/2 + Equality(1, 1) == Equality(x/2 + 1, 2)


def test_M2():
    # The roots of this equation should all be real. Note that this doesn't test
    # that they are correct.
    sol = solve(3*x**3 - 18*x**2 + 33*x - 19, x)
    assert all(expand(x, complex=True).is_real for x in sol)


@XFAIL
def test_M5():
    assert solve(x**6 - 9*x**4 - 4*x**3 + 27*x**2 - 36*x - 23, x) == [2**(1/3) + sqrt(3), 2**(1/3) - sqrt(3), +sqrt(3) - 1/2**(2/3) + I*sqrt(3)/2**(2/3), +sqrt(3) - 1/2**(2/3) - I*sqrt(3)/2**(2/3), -sqrt(3) - 1/2**(2/3) + I*sqrt(3)/2**(2/3), -sqrt(3) - 1/2**(2/3) - I*sqrt(3)/2**(2/3)]


def test_M6():
    assert set(solve(x**7 - 1, x)) == set([cos(n*2*pi/7) + I*sin(n*2*pi/7) for n in range(0, 7)])
    # The paper asks for exp terms, but sin's and cos's may be acceptable


def test_M7():
    assert set(solve(x**8 - 8*x**7 + 34*x**6 - 92*x**5 + 175*x**4 - 236*x**3 +
        226*x**2 - 140*x + 46, x)) == set([
        1 + sqrt(2)*I*sqrt(sqrt(-3 + 4*sqrt(3)) + 3)/2,
        1 + sqrt(2)*sqrt(-3 + sqrt(-3 + 4*sqrt(3)))/2,
        1 - sqrt(2)*sqrt(-3 + I*sqrt(3 + 4*sqrt(3)))/2,
        1 - sqrt(2)*I*sqrt(sqrt(-3 + 4*sqrt(3)) + 3)/2,
        1 + sqrt(2)*sqrt(-3 - I*sqrt(3 + 4*sqrt(3)))/2,
        1 + sqrt(2)*sqrt(-3 + I*sqrt(3 + 4*sqrt(3)))/2,
        1 - sqrt(2)*sqrt(-3 - I*sqrt(3 + 4*sqrt(3)))/2,
        1 - sqrt(2)*sqrt(-3 + sqrt(-3 + 4*sqrt(3)))/2,
        ])


@XFAIL  # There are an infinite number of solutions.
def test_M8():
    z = symbols('z', complex=True)
    assert set(solve(exp(2*x) + 2*exp(x) + 1 - z, x)) == \
        set([log(1 + z - 2*sqrt(z))/2, log(1 + z + 2*sqrt(z))/2])
    # This one could be simplified better (the 1/2 could be pulled into the log
    # as a sqrt, and the function inside the log can be factored as a square,
    # giving [log(sqrt(z) - 1), log(sqrt(z) + 1)]). Also, there should be an
    # infinite number of solutions.
    # x = {log(sqrt(z) - 1), log(sqrt(z) + 1) + i pi} [+ n 2 pi i, + n 2 pi i]
    # where n is an arbitrary integer.  See url of detailed output above.


@XFAIL
def test_M9():
    x = symbols('x', complex=True)
    raise NotImplementedError("solve(exp(2-x**2)-exp(-x),x) has complex solutions.")


def test_M10():
    assert solve(exp(x) - x, x) == [-LambertW(-1)]


@XFAIL
def test_M11():
    assert solve(x**x - x, x) == [-1, 1]


def test_M12():
    # TODO: x = [-1, 2*(+/-asinh(1)*I + n*pi}, 3*(pi/6 + n*pi/3)]
    assert solve((x + 1)*(sin(x)**2 + 1)**2*cos(3*x)**3, x) == [
        -1, pi/6, pi/2,
           - I*log(1 + sqrt(2)),      I*log(1 + sqrt(2)),
        pi - I*log(1 + sqrt(2)), pi + I*log(1 + sqrt(2)),
    ]


def test_M13():
    assert solve(sin(x) - cos(x), x) == [-3*pi/4, pi/4]


def test_M14():
    assert solve(tan(x) - 1, x) == [pi/4]


def test_M15():
    assert solve(sin(x) - S.Half) == [pi/6, 5*pi/6]


def test_M16():
    assert solve(sin(x) - tan(x), x) == [0, 2*pi]


@XFAIL
def test_M17():
    assert solve(asin(x) - atan(x),x) == [0]


@XFAIL
def test_M18():
    assert solve(acos(x) - atan(x), x) == [sqrt((sqrt(5) - 1)/2)]


def test_M19():
    assert solve((x - 2)/x**R(1, 3), x) == [2]


def test_M20():
    assert solve(sqrt(x**2 + 1) - x + 2, x) == []


def test_M21():
    assert solve(x + sqrt(x) - 2) == [1]


def test_M22():
    assert solve(2*sqrt(x) + 3*x**R(1, 4) - 2) == [R(1, 16)]


def test_M23():
    x = symbols('x', complex=True)

    assert solve(x - 1/sqrt(1 + x**2)) == [
        simplify(-I*sqrt((sqrt(5) + 1)/2)),
        simplify(   sqrt((sqrt(5) - 1)/2)),
    ]


def test_M24():
    solution = solve(1 - binomial(m, 2)*2**k, k)
    answer = log(2/(m*(m - 1)), 2)
    assert solution[0].expand() == answer.expand()


def test_M25():
    a, b, c, d = symbols(':d', positive=True)
    x = symbols('x')
    assert solve(a*b**x - c*d**x, x)[0].expand() == (log(c/a)/log(b/d)).expand()


def test_M26():
    assert solve(sqrt(log(x)) - log(sqrt(x))) == [1, exp(4)]


@XFAIL
def test_M27():
    x = symbols('x', real=True)
    b = symbols('b', real=True)
    with assuming(Q.is_true(sin(cos(1/E**2) + 1) + b > 0)):
        solve(log(acos(asin(x**R(2,3) - b) - 1)) + 2, x) == [-b - sin(1 + cos(1/e**2))**R(3/2), b + sin(1 + cos(1/e**2))**R(3/2)]


@XFAIL
def test_M28():
    assert solve(5*x + exp((x - 5)/2) - 8*x**3, x, assume=Q.real(x)) == [-0.784966, -0.016291, 0.802557]


def test_M29():
    assert solve(abs(x - 1) - 2) == [-1, 3]


@XFAIL
def test_M30():
    assert solve(abs(2*x + 5) - abs(x - 2),x, assume=Q.real(x)) == [-1, -7]


@XFAIL
def test_M31():
    assert solve(1 - abs(x) - max(-x - 2, x - 2),x, assume=Q.real(x)) == [-3/2, 3/2]


@XFAIL
def test_M32():
    assert solve(max(2 - x**2, x)- max(-x, (x**3)/9), assume=Q.real(x)) == [-1, 3]


@XFAIL
def test_M33():
    # Second answer can be written in another form. The second answer is the root of x**3 + 9*x**2 - 18 = 0 in the interval (-2, -1).
    assert solve(max(2 - x**2, x) - x**3/9, assume=Q.real(x)) == [-3, -1.554894, 3]


@XFAIL
def test_M34():
    z = symbols('z', complex=True)
    assert solve((1 + I) * z + (2 - I) * conjugate(z) + 3*I, z) == [2 + 3*I]


def test_M35():
    x, y = symbols('x y', real=True)
    assert solve((3*x - 2*y - I*y + 3*I).as_real_imag()) == {y: 3, x: 2}


@XFAIL
def test_M36():
    assert solve(f**2 + f - 2, x) == [Eq(f(x), 1), Eq(f(x), -2)]


def test_M37():
    assert solve([x + y + z - 6, 2*x + y + 2*z - 10, x + 3*y + z - 10 ]) == {x: -z + 4, y: 2}


@slow
def test_M38():
    variabes = vring("k1:50", vfield("a,b,c", ZZ).to_domain())
    system = [
        -b*k8/a + c*k8/a, -b*k11/a + c*k11/a, -b*k10/a + c*k10/a + k2, -k3 - b*k9/a + c*k9/a,
        -b*k14/a + c*k14/a, -b*k15/a + c*k15/a, -b*k18/a + c*k18/a - k2, -b*k17/a + c*k17/a,
        -b*k16/a + c*k16/a + k4, -b*k13/a + c*k13/a - b*k21/a + c*k21/a + b*k5/a - c*k5/a,
        b*k44/a - c*k44/a, -b*k45/a + c*k45/a, -b*k20/a + c*k20/a, -b*k44/a + c*k44/a,
        b*k46/a - c*k46/a, b**2*k47/a**2 - 2*b*c*k47/a**2 + c**2*k47/a**2, k3, -k4,
        -b*k12/a + c*k12/a - a*k6/b + c*k6/b, -b*k19/a + c*k19/a + a*k7/c - b*k7/c,
        b*k45/a - c*k45/a, -b*k46/a + c*k46/a, -k48 + c*k48/a + c*k48/b - c**2*k48/(a*b),
        -k49 + b*k49/a + b*k49/c - b**2*k49/(a*c), a*k1/b - c*k1/b, a*k4/b - c*k4/b,
        a*k3/b - c*k3/b + k9, -k10 + a*k2/b - c*k2/b, a*k7/b - c*k7/b, -k9, k11,
        b*k12/a - c*k12/a + a*k6/b - c*k6/b, a*k15/b - c*k15/b, k10 + a*k18/b - c*k18/b,
        -k11 + a*k17/b - c*k17/b, a*k16/b - c*k16/b, -a*k13/b + c*k13/b + a*k21/b - c*k21/b + a*k5/b - c*k5/b,
        -a*k44/b + c*k44/b, a*k45/b - c*k45/b, a*k14/c - b*k14/c + a*k20/b - c*k20/b,
        a*k44/b - c*k44/b, -a*k46/b + c*k46/b, -k47 + c*k47/a + c*k47/b - c**2*k47/(a*b),
        a*k19/b - c*k19/b, -a*k45/b + c*k45/b, a*k46/b - c*k46/b, a**2*k48/b**2 - 2*a*c*k48/b**2 + c**2*k48/b**2,
        -k49 + a*k49/b + a*k49/c - a**2*k49/(b*c), k16, -k17, -a*k1/c + b*k1/c,
        -k16 - a*k4/c + b*k4/c, -a*k3/c + b*k3/c, k18 - a*k2/c + b*k2/c, b*k19/a - c*k19/a - a*k7/c + b*k7/c,
        -a*k6/c + b*k6/c, -a*k8/c + b*k8/c, -a*k11/c + b*k11/c + k17, -a*k10/c + b*k10/c - k18,
        -a*k9/c + b*k9/c, -a*k14/c + b*k14/c - a*k20/b + c*k20/b, -a*k13/c + b*k13/c + a*k21/c - b*k21/c - a*k5/c + b*k5/c,
        a*k44/c - b*k44/c, -a*k45/c + b*k45/c, -a*k44/c + b*k44/c, a*k46/c - b*k46/c,
        -k47 + b*k47/a + b*k47/c - b**2*k47/(a*c), -a*k12/c + b*k12/c, a*k45/c - b*k45/c,
        -a*k46/c + b*k46/c, -k48 + a*k48/b + a*k48/c - a**2*k48/(b*c),
        a**2*k49/c**2 - 2*a*b*k49/c**2 + b**2*k49/c**2, k8, k11, -k15, k10 - k18,
        -k17, k9, -k16, -k29, k14 - k32, -k21 + k23 - k31, -k24 - k30, -k35, k44,
        -k45, k36, k13 - k23 + k39, -k20 + k38, k25 + k37, b*k26/a - c*k26/a - k34 + k42,
        -2*k44, k45, k46, b*k47/a - c*k47/a, k41, k44, -k46, -b*k47/a + c*k47/a,
        k12 + k24, -k19 - k25, -a*k27/b + c*k27/b - k33, k45, -k46, -a*k48/b + c*k48/b,
        a*k28/c - b*k28/c + k40, -k45, k46, a*k48/b - c*k48/b, a*k49/c - b*k49/c,
        -a*k49/c + b*k49/c, -k1, -k4, -k3, k15, k18 - k2, k17, k16, k22, k25 - k7,
        k24 + k30, k21 + k23 - k31, k28, -k44, k45, -k30 - k6, k20 + k32, k27 + b*k33/a - c*k33/a,
        k44, -k46, -b*k47/a + c*k47/a, -k36, k31 - k39 - k5, -k32 - k38, k19 - k37,
        k26 - a*k34/b + c*k34/b - k42, k44, -2*k45, k46, a*k48/b - c*k48/b,
        a*k35/c - b*k35/c - k41, -k44, k46, b*k47/a - c*k47/a, -a*k49/c + b*k49/c,
        -k40, k45, -k46, -a*k48/b + c*k48/b, a*k49/c - b*k49/c, k1, k4, k3, -k8,
        -k11, -k10 + k2, -k9, k37 + k7, -k14 - k38, -k22, -k25 - k37, -k24 + k6,
        -k13 - k23 + k39, -k28 + b*k40/a - c*k40/a, k44, -k45, -k27, -k44, k46,
        b*k47/a - c*k47/a, k29, k32 + k38, k31 - k39 + k5, -k12 + k30, k35 - a*k41/b + c*k41/b,
        -k44, k45, -k26 + k34 + a*k42/c - b*k42/c, k44, k45, -2*k46, -b*k47/a + c*k47/a,
        -a*k48/b + c*k48/b, a*k49/c - b*k49/c, k33, -k45, k46, a*k48/b - c*k48/b,
        -a*k49/c + b*k49/c
        ]
    solution = {
        k49: 0, k48: 0, k47: 0, k46: 0, k45: 0, k44: 0, k41: 0, k40: 0,
        k38: 0, k37: 0, k36: 0, k35: 0, k33: 0, k32: 0, k30: 0, k29: 0,
        k28: 0, k27: 0, k25: 0, k24: 0, k22: 0, k21: 0, k20: 0, k19: 0,
        k18: 0, k17: 0, k16: 0, k15: 0, k14: 0, k13: 0, k12: 0, k11: 0,
        k10: 0, k9:  0, k8:  0, k7:  0, k6:  0, k5:  0, k4:  0, k3:  0,
        k2:  0, k1:  0,
        k34: b/c*k42, k31: k39, k26: a/c*k42, k23: k39
    }
    assert solve_lin_sys(system, variabes) == solution

def test_M39():
    x, y, z = symbols('x y z', complex=True)
    assert solve([x**2*y + 3*y*z - 4, -3*x**2*z + 2*y**2 + 1, 2*y*z**2 - z**2 - 1 ]) ==\
            [{y: 1, z: 1, x: -1}, {y: 1, z: 1, x: 1},\
             {y: sqrt(2)*I, z: R(1,3) - sqrt(2)*I/3, x: -sqrt(-1 - sqrt(2)*I)},\
             {y: sqrt(2)*I, z: R(1,3) - sqrt(2)*I/3, x: sqrt(-1 - sqrt(2)*I)},\
             {y: -sqrt(2)*I, z: R(1,3) + sqrt(2)*I/3, x: -sqrt(-1 + sqrt(2)*I)},\
             {y: -sqrt(2)*I, z: R(1,3) + sqrt(2)*I/3, x: sqrt(-1 + sqrt(2)*I)}]

# N. Inequalities


def test_N1():
    assert ask(Q.is_true(E**pi > pi**E))


@XFAIL
def test_N2():
    x = symbols('x', real=True)
    assert ask(Q.is_true(x**4 - x + 1 > 0))
    assert ask(Q.is_true(x**4 - x + 1 > 1)) == False


@XFAIL
def test_N3():
    x = symbols('x', real=True)
    assert ask(Q.is_true(And(Lt(-1, x), Lt(x, 1))), Q.is_true(abs(x) < 1 ))

@XFAIL
def test_N4():
    x, y = symbols('x y', real=True)
    assert ask(Q.is_true(2*x**2 > 2*y**2), Q.is_true((x > y) & (y > 0)))


@XFAIL
def test_N5():
    x, y, k = symbols('x y k', real=True)
    assert ask(Q.is_true(k*x**2 > k*y**2), Q.is_true((x > y) & (y > 0) & (k > 0)))


@XFAIL
def test_N6():
    x, y, k, n = symbols('x y k n', real=True)
    assert ask(Q.is_true(k*x**n > k*y**n), Q.is_true((x > y) & (y > 0) & (k > 0) & (n > 0)))


@XFAIL
def test_N7():
    x, y = symbols('x y', real=True)
    assert ask(Q.is_true(y > 0), Q.is_true((x > 1) & (y >= x - 1)))


@XFAIL
def test_N8():
    x, y, z = symbols('x y z', real=True)
    assert ask(Q.is_true((x == y) & (y == z)), Q.is_true((x >= y) & (y >= z) & (z >= x)))


def test_N9():
    with assuming(Q.real(x)):
        assert solve(abs(x-1) > 2) == Or(x < -1, x > 3)


def test_N10():
    p=(x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5)
    assert solve(expand(p) < 0, assume=Q.real(x)) == Or( And(Lt(2, x), Lt(x, 3)), And(Lt(4, x), Lt(x, 5)), Lt(x, 1))


def test_N11():
    assert solve(6/(x - 3) <= 3, assume=Q.real(x)) == Or(5 <= x, x < 3)


@XFAIL
def test_N12():
    assert solve(sqrt(x)<2, assume=Q.real(x)) == And(Le(0,x),Lt(x,4))

@XFAIL
def test_N13():
    assert solve(sin(x)<2, assume=Q.real(x)) == S.Reals # unsupported

@XFAIL
def test_N14():
    assert solve(sin(x)<1, assume=Q.real(x)) ==  Ne(x,pi/2) # unsupported should return

@XFAIL
def test_N15():
    r, t = symbols('r t', real=True)
    solve(abs(2*r*(cos(t)-1)+1)<=1,r) # unsupported

@XFAIL
def test_N16():
    r, t = symbols('r t', real=True)
    solve((r**2)*((cos(t) - 4)**2)*sin(t)**2 < 9, r)

@XFAIL
def test_N17():
    assert solve(x+y>0, x-y<0)


def test_O1():
    M = Matrix((1 + I, -2, 3*I))
    assert sqrt(expand(M.dot(M.H))) == sqrt(15)

def test_O2():
    assert Matrix((2,2,-3)).cross(Matrix((1,3,1))) == Matrix([[11, -5, 4]])

@slow
def test_O3():
    (va, vb, vc, vd)  = MV.setup('va vb vc vd')
    assert (va^vb)|(vc^vd) == -(va|vc)*(vb|vd) + (va|vd)*(vb|vc)

def test_O4():
    (ex,ey,ez,grad) = MV.setup('e*x|y|z',metric='[1,1,1]',coords=(x,y,z))
    F=ex*(x*y*z)+ey*((x*y*z)**2)+ez*((y**2)*(z**3))
    cu=grad^F
    assert cu|ex == (x*z*(-2*y**2*z + 1))*ey + x*y*ez
    assert cu|ey == (x*z*(2*y**2*z - 1))*ex + (2*y*z*(x**2*y - z**2))*ez
    assert cu|ez == -x*y*ex + (2*y*z*(-x**2*y + z**2))*ey
    #assert cu == (x*z*(2*y**2*z - 1))*ex^ey - x*y*ex^ez + (2*y*z*(-x**2*y + z**2))*ey^ez

@XFAIL
@slow
def test_O5():
    (ex,ey,ez,grad) = MV.setup('e*x|y|z',metric='[1,1,1]',coords=(x,y,z))
    f = MV('f','vector',fct=True)
    g = MV('g','vector',fct=True)
    assert grad|(f^g)-g|(grad^f)+f|(grad^g)  == 0

#O8-O9 MISSING!!
def test_O10():
    L = [Matrix([2,3,5]), Matrix([3,6,2]), Matrix([8,3,6])]
    assert GramSchmidt(L) == [Matrix([
                                [2],
                                [3],
                                [5]]), Matrix([
                                [ 23/19],
                                [ 63/19],
                                [-47/19]]), Matrix([
                                [ 1692/353],
                                [-1551/706],
                                [ -423/706]])]

@XFAIL
def test_P1():
    raise NotImplementedError("Matrix property/function to extract Nth diagnoal not implemented")


def test_P2():
    M = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    M.row_del(1)
    M.col_del(2)
    assert M == Matrix([
                    [1, 2],
                    [7, 8]])

@XFAIL
def test_P3():
    A = Matrix([
        [11, 12, 13, 14],
        [21, 22, 23, 24],
        [31, 32, 33, 34],
        [41, 42, 43, 44]])

    A11 = A[0:3,1:4]
    A12 = A[(0,1,3),(2,0,3)] # unsupported rasies exception
    A21 = A
    A221 = A[0:2,2:4]
    A222 = A[(3,0),(2,1)] # unsupported rasies exception
    A22 = BlockMatrix([A221,A222])
    B= BlockMatrix([[A11,A12],[A21,A22]])
    assert B  ==  Matrix([
        [12,13,14,13,11,14],
        [22,22,24,23,21,24],
        [32,33,34,43,41,44],
        [11,12,13,14,13,14],
        [21,22,23,24,23,24],
        [31,32,33,34,43,42],
        [41,42,43,44,13,12]])

@XFAIL
def test_P4():
    raise NotImplementedError("Block matrix diagonalization not supported")

@XFAIL
def test_P5():
    M = Matrix([[7,11],[3,8]])
    assert  M % 2 == Matrix([ # Raises exception % not supported for matrixes
                        [1, 1],
                        [1, 0]])

def test_P5_workaround():
    M = Matrix([[7,11],[3,8]])
    assert  M.applyfunc(lambda i:i%2) == Matrix([
                                            [1, 1],
                                            [1, 0]])
def test_P6():
    M = Matrix([[cos(x),sin(x)],[-sin(x),cos(x)]])
    assert  M.diff(x,2) == Matrix([
        [-cos(x), -sin(x)],
        [ sin(x), -cos(x)]])

def test_P7():
    M = Matrix([[x,y]])*(z*Matrix([
                    [1,3,5],
                    [2,4,6]])+Matrix([
                        [7,-9,11],
                        [-8,10,-12]]))
    assert M == Matrix([
        [x*(z + 7) + y*(2*z - 8), x*(3*z - 9) + y*(4*z + 10), x*(5*z + 11) + y*(6*z - 12)]])

@XFAIL
def test_P8():
    M=Matrix([[1,-2*I],[-3*I,4]])
    assert M.norm(ord=S.Infinity) == 7 # Matrix.norm(ord=inf) not implemented

def test_P9():
    a, b, c = symbols('a b c', real=True)
    M=Matrix([[a/(b*c), 1/c, 1/b], [1/c, b/(a*c), 1/a], [1/b, 1/a, c/(a*b)]])
    assert factor(M.norm('fro')) == (a**2 + b**2 + c**2)/(abs(a)*abs(b)*abs(c))

@XFAIL # conugate(f(4-5*i)) is not simplified to f(4+5*I)
def test_P10():
    M=Matrix([[1,2+3*I],[f(4-5*i),6]])
    assert M.H == Matrix([[1,f(4+5*I)],[2+3*I,6]])

@XFAIL
def test_P11():
    #raise NotImplementedError("Matrix([[x,y],[1,x*y]]).inv() not simplifying to extract common factor")
    assert Matrix([[x,y],[1,x*y]]).inv() ==  (1/(x**2-1))*Matrix([
                                                            [x, -1],
                                                            [-1/y, x/y]])

def test_P12():
    A11 = MatrixSymbol('A11',n,n)
    A12 = MatrixSymbol('A12',n,n)
    A22 = MatrixSymbol('A22',n,n)
    B = BlockMatrix([[A11,A12],[ZeroMatrix(n,n),A22]])
    assert block_collapse(B.I)  == BlockMatrix([
                                    [A11.I,          (-1)*A11.I*A12*A22.I],
                                    [ZeroMatrix(n,n), A22.I]])

def test_P13():
    M=Matrix([
        [ 1,     x-2,         x-3     ],
        [x-1, x**2-3*x+6,   x**2-3*x-2  ],
        [x-2,   x**2-8,   2*(x**2)-12*x+14]])
    L, U, _= M.LUdecomposition()
    assert simplify(L) == Matrix([
                            [    1,     0, 0],
                            [x - 1,     1, 0],
                            [x - 2, x - 3, 1]])
    assert simplify(U) == Matrix([
                            [1, x - 2, x - 3],
                            [0,     4, x - 5],
                            [0,     0, x - 7]])

def test_P14():
    M = Matrix([
            [1, 2, 3, 1, 3],
            [3, 2, 1, 1, 7],
            [0, 2, 4, 1, 1],
            [1, 1, 1, 1, 4]])
    R,_ = M.rref()
    assert R == Matrix([
                    [1, 0, -1, 0,  2],
                    [0, 1,  2, 0, -1],
                    [0, 0,  0, 1,  3],
                    [0, 0,  0, 0,  0]])

def test_P15():
    M = Matrix([
        [-1, 3, 7, -5],
        [4, -2, 1, 3],
        [2, 4, 15, -7]])
    assert M.rank() == 2

def test_P16():
    M = Matrix([
        [2*sqrt(2), 8],
        [6*sqrt(6), 24*sqrt(3)]])
    assert M.rank() == 1

@XFAIL
def test_P17():
    t = symbols('t', real=True)
    M=Matrix([
        [sin(2*t), cos(2*t)],
        [2*(1 - (cos(t)**2))*cos(t), (1 - 2*(sin(t)**2))*sin(t)]])
    assert M.rank() == 1

def test_P18():
    M = Matrix([
        [1, 0, -2, 0],
        [-2, 1, 0, 3],
        [-1, 2, -6, 6]])
    assert M.nullspace() == [   Matrix([
                                    [2],
                                    [4],
                                    [1],
                                    [0]]),
                                Matrix([
                                    [ 0],
                                    [-3],
                                    [ 0],
                                    [ 1]])]

def test_P19():
    w = symbols('w')
    M = Matrix([[1,   1,   1,   1  ],
                [w,   x,   y,   z  ],
                [w**2, x**2, y**2, z**2],
                [w**3, x**3, y**3, z**3]])
    assert M.det()  == w**3*x**2*y - w**3*x**2*z - w**3*x*y**2 + w**3*x*z**2 + w**3*y**2*z - w**3*y*z**2 - w**2*x**3*y + w**2*x**3*z + w**2*x*y**3 - w**2*x*z**3 - w**2*y**3*z + w**2*y*z**3 + w*x**3*y**2 - w*x**3*z**2 - w*x**2*y**3 + w*x**2*z**3 + w*y**3*z**2 - w*y**2*z**3 - x**3*y**2*z + x**3*y*z**2 + x**2*y**3*z - x**2*y*z**3 - x*y**3*z**2 + x*y**2*z**3

@XFAIL
def test_P20():
    raise NotImplementedError("Matrix minimal polynomial not supported")


def test_P21():
    M=Matrix([
        [ 5, -3, -7],
        [-2,  1,  2],
        [ 2, -3, -4]])
    assert M.charpoly(x).as_expr() == x**3 - 2*x**2 - 5*x + 6

@slow
def test_P22():
#   Wester test requires calculating eigenvalues for a matrix of dimension 100
#   This currently takes forever with sympy
#    M=(2-x)*eye(100);
#    assert M.eigenvals() == {-x + 2: 100}
#   So we will speed-up for the moment the test checking only for dimension 12
    M=(2-x)*eye(12)
    assert M.eigenvals() == {-x + 2: 12}

def test_P23():
    M = Matrix([
        [2, 1, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [0, 1, 2, 1, 0],
        [0, 0, 1, 2, 1],
        [0, 0, 0, 1, 2]])
    assert M.eigenvals() == {
        S('1'): 1,
        S('2'): 1,
        S('3'): 1,
        S('sqrt(3) + 2'): 1,
        S('-sqrt(3) + 2'): 1}

def test_P24():
    M = Matrix([
        [ 611,  196, -192,  407,   -8,  -52,  -49,   29],
        [ 196,  899,  113, -192,  -71,  -43,   -8,  -44],
        [-192,  113,  899,  196,   61,   49,    8,   52],
        [ 407, -192,  196,  611,    8,   44,   59,  -23],
        [  -8,  -71,   61,    8,  411, -599,  208,  208],
        [ -52,  -43,   49,   44, -599,  411,  208,  208],
        [ -49,   -8,    8,   59,  208,  208,   99, -911],
        [  29,  -44,   52,  -23,  208,  208, -911,   99]])
    assert M.eigenvals() == {
        S('0'): 1,
        S('10*sqrt(10405)'): 1,
        S('100*sqrt(26) + 510'): 1,
        S('1000'): 2,
        S('-100*sqrt(26) + 510'): 1,
        S('-10*sqrt(10405)'): 1,
        S('1020'): 1}

def test_P25():
    MF = N(Matrix([
        [ 611,  196, -192,  407,   -8,  -52,  -49,   29],
        [ 196,  899,  113, -192,  -71,  -43,   -8,  -44],
        [-192,  113,  899,  196,   61,   49,    8,   52],
        [ 407, -192,  196,  611,    8,   44,   59,  -23],
        [  -8,  -71,   61,    8,  411, -599,  208,  208],
        [ -52,  -43,   49,   44, -599,  411,  208,  208],
        [ -49,   -8,    8,   59,  208,  208,   99, -911],
        [  29,  -44,   52,  -23,  208,  208, -911,   99]]))
    assert [float(i) for i in sorted(MF.eigenvals())]  == [-1020.0490184299969, 0.0, 0.09804864072151699, 1000.0, 1019.9019513592784, 1020.0, 1020.0490184299969]

def test_P26():
    a0,a1,a2,a3,a4 = symbols('a0 a1 a2 a3 a4')
    M = Matrix([
    [-a4, -a3, -a2, -a1, -a0,  0,  0,  0,  0],
    [  1,   0,   0,   0,   0,  0,  0,  0,  0],
    [  0,   1,   0,   0,   0,  0,  0,  0,  0],
    [  0,   0,   1,   0,   0,  0,  0,  0,  0],
    [  0,   0,   0,   1,   0,  0,  0,  0,  0],
    [  0,   0,   0,   0,   0, -1, -1,  0,  0],
    [  0,   0,   0,   0,   0,  1,  0,  0,  0],
    [  0,   0,   0,   0,   0,  0,  1, -1, -1],
    [  0,   0,   0,   0,   0,  0,  0,  1,  0]])
    assert M.eigenvals() == {
        S('-1/2 - sqrt(3)*I/2'): 2,
        S('-1/2 + sqrt(3)*I/2'): 2}

def test_P27():
    a = symbols('a')
    M = Matrix([
    [a,  0, 0, 0, 0],
    [0,  0, 0, 0, 1],
    [0,  0, a, 0, 0],
    [0,  0, 0, a, 0],
    [0, -2, 0, 0, 2]])
    M.eigenvects() == [
                        (a, 3, [Matrix([
                        [1],
                        [0],
                        [0],
                        [0],
                        [0]]), Matrix([
                        [0],
                        [0],
                        [1],
                        [0],
                        [0]]), Matrix([
                        [0],
                        [0],
                        [0],
                        [1],
                        [0]])]), (1 - I, 1, [Matrix([
                        [          0],
                        [-1/(-1 + I)],
                        [          0],
                        [          0],
                        [          1]])]), (1 + I, 1, [Matrix([
                        [          0],
                        [-1/(-1 - I)],
                        [          0],
                        [          0],
                        [          1]])])]
@XFAIL
def test_P28():
    raise NotImplementedError("Generalized eigen vectors not supported https://code.google.com/p/sympy/issues/detail?id=2194")

@XFAIL
def test_P29():
    raise NotImplementedError("Generalized eigen vectors not supported https://code.google.com/p/sympy/issues/detail?id=2194")

def test_P30():
    M = Matrix([
        [1,  0,  0,  1, -1],
        [0,  1, -2,  3, -3],
        [0,  0, -1,  2, -2],
        [1, -1,  1,  0,  1],
        [1, -1,  1, -1,  2]])
    P,J = M.jordan_form()
    assert J == Matrix([
        [-1, 0, 0, 0, 0],
        [ 0, 1, 1, 0, 0],
        [ 0, 0, 1, 0, 0],
        [ 0, 0, 0, 1, 1],
        [ 0, 0, 0, 0, 1]])

@XFAIL
def test_P31():
    raise NotImplementedError("Smith normal form not implemented")

def test_P32():
    M=Matrix([
        [1, -2],
        [2, 1]])
    assert exp(M).rewrite(cos).simplify() == Matrix([
                                        [E*cos(2), -E*sin(2)],
                                        [E*sin(2),  E*cos(2)]])

def test_P33():
    w,t = symbols('w t')
    M = Matrix([
        [0, 1,    0,     0  ],
        [0, 0,    0,     2*w],
        [0, 0,    0,     1  ],
        [0, -2*w, 3*w**2, 0  ]])
    assert exp(M*t).rewrite(cos).expand() == Matrix([
        [1, -3*t + 4*sin(t*w)/w,  6*t*w - 6*sin(t*w), -2*cos(t*w)/w + 2/w],
        [0,      4*cos(t*w) - 3, -6*w*cos(t*w) + 6*w,          2*sin(t*w)],
        [0,  2*cos(t*w)/w - 2/w,     -3*cos(t*w) + 4,          sin(t*w)/w],
        [0,         -2*sin(t*w),        3*w*sin(t*w),            cos(t*w)]])


@XFAIL
def test_P34():
    a,b,c = symbols('a b c',real=True)
    M=Matrix([
    [a, 1, 0, 0, 0, 0],
    [0, a, 0, 0, 0, 0],
    [0, 0, b, 0, 0, 0],
    [0, 0, 0, c, 1, 0],
    [0, 0, 0, 0, c, 1],
    [0, 0, 0, 0, 0, c]])
    # raises exception, sin(M) not supported. exp(M*I) also not supported
    # https://code.google.com/p/sympy/issues/detail?id=3119
    assert sin(M) == Matrix([
                        [sin(a), cos(a), 0, 0, 0, 0],
                        [0, sin(a), 0, 0, 0, 0	],
                        [0, 0, sin(b), 0, 0, 0	],
                        [0, 0, 0, sin(c), cos(c), -sin(c)/2],
                        [0, 0, 0, 0, sin(c), cos(c)],
                        [0, 0, 0, 0, 0, sin(c)]])

@XFAIL
def test_P35():
    M = pi/2*Matrix([[2, 1, 1], [2, 3, 2], [1, 1, 2]])
    # raises exception, sin(M) not supported. exp(M*I) also not supported
    # https://code.google.com/p/sympy/issues/detail?id=3119
    assert sin(M) ==  eye(3)

@XFAIL
def test_P36():
    M=Matrix([
        [10, 7],
        [7, 17]])
    # sqrt(M) not performmed
    # sqrtdenest not simplifying sqrt(M)
    assert sqrtdenest(M**Rational(1,2)) == Matrix([[3, 1], [1, 4]])
@XFAIL
def test_P37():
    M=Matrix([
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1]])
    #raises NotImplementedError: Implemented only for diagonalizable matrices
    M**Rational(1,2)

@XFAIL
def test_P38():
    M=Matrix([
    [0, 1, 0],
    [0, 0, 0],
    [0, 0, 0]])
    #raises NotImplementedError: Implemented only for diagonalizable matrices
    M**Rational(1,2)

@XFAIL
def test_P39():
    '''
    M=Matrix([
        [1, 1],
        [2, 2],
        [3, 3]])
    M.SVD()
    '''
    raise NotImplementedError("Singular value decomposition not implemented normal form not implemented")

def test_P40():
    r,t = symbols('r t',real=True)
    M=Matrix([r*cos(t), r*sin(t)])
    assert M.jacobian(Matrix([r, t])) == Matrix([
                                [cos(t), -r*sin(t)],
                                [sin(t),  r*cos(t)]])

def test_P41():
    r,t = symbols('r t',real=True)
    assert hessian(r**2*sin(t),(r,t)) == Matrix([
                                            [  2*sin(t),   2*r*cos(t)],
                                            [2*r*cos(t), -r**2*sin(t)]])

def test_P42():
    assert wronskian([cos(x), sin(x)], x).simplify() == 1

def test_P43():
    def __my_jacobian(M,Y):
        return Matrix([M.diff(v).T for v in Y]).T
    r,t = symbols('r t',real=True)
    M=Matrix([r*cos(t), r*sin(t)])
    assert __my_jacobian(M,[r,t]) == Matrix([
                                [cos(t), -r*sin(t)],
                                [sin(t),  r*cos(t)]])

def test_P44():
    def __my_hessian(f,Y):
        V=Matrix([diff(f,v) for v in Y])
        return  Matrix([V.T.diff(v) for v in Y])
    r,t = symbols('r t',real=True)
    assert __my_hessian(r**2*sin(t),(r,t)) == Matrix([
                                            [  2*sin(t),   2*r*cos(t)],
                                            [2*r*cos(t), -r**2*sin(t)]])

def test_P45():
    def __my_wronskian(Y,v):
        return  Matrix([Matrix(Y).T.diff(x,n) for n in range(0,len(Y))]).det()
    assert __my_wronskian([cos(x), sin(x)], x).simplify() == 1

# Q1-Q6  Tensor tests missing

@XFAIL
def test_R1():
    i,n = symbols('i n', integer=True, positive=True)
    xn=MatrixSymbol('xn',n,1)
    Sm = Sum((xn[i,0]-Sum(xn[j,0],(j,0,n-1))/n)**2,(i,0,n-1))
    Sm.doit() # raises AttributeError: 'str' object has no attribute 'is_Piecewise'

@XFAIL
def test_R2():
    m,b = symbols('m b',real=True)
    i,n = symbols('i n', integer=True, positive=True)
    xn=MatrixSymbol('xn',n,1)
    yn=MatrixSymbol('yn',n,1)
    f=Sum((yn[i,0]-m*xn[i,0]-b)**2,(i,0,n-1))
    f1=diff(f,m)
    f2=diff(f,b)
    solve((f1,f2),m,b) # raises AttributeError: 'str' object has no attribute 'is_Piecewise'

@XFAIL
def test_R3():
    n,k = symbols('n k', integer=True, positive=True)
    sk = ((-1)**k) * (binomial(2*n, k))**2
    Sm = Sum(sk, (k,1,oo))
    T = Sm.doit()
    assert T.combsimp() == (-1)**n*binomial(2*n, n)# returns -((-1)**n*factorial(2*n) - (factorial(n))**2)*exp_polar(-I*pi)/(factorial(n))**2

@XFAIL
def test_R4():
    n,k = symbols('n k', integer=True, positive=True)
    sk = binomial(n, k)/(2**n) - binomial(n + 1, k)/(2**(n + 1))
    Sm = Sum(sk, (k,1,oo))
    T = Sm.doit()
    assert T.combsimp() == 2**(-n-1)*binomial(n,k-1) # returns -2**(-n)/2

@XFAIL
def test_R5():
    a,b,c,n,k = symbols('a b c n k', integer=True, positive=True)
    sk = ((-1)**k) * binomial(a+b, a+k) * binomial(b+c, b+k) * binomial(c+a, c+k)
    Sm = Sum(sk, (k,1,oo))
    T = Sm.doit() # hypergeometric series not calculated
    assert T == factorial(a+b+c)/(factorial(a)*factorial(b)*factorial(c))

@XFAIL
def test_R6():
    n,k = symbols('n k', integer=True, positive=True)
    gn = MatrixSymbol('gn',n+1,1)
    Sm = Sum(gn[k,0]-gn[k-1,0],(k,1,n+1))
    assert Sm.doit() == -gn[0, 0] + gn[n+1, 0] #raises AttributeError: 'str' object has no attribute 'is_Piecewise'

def test_R7():
    n,k = symbols('n k', integer=True, positive=True)
    T = Sum(k**3,(k,1,n)).doit()
    assert T.factor() == n**2*(n + 1)**2/4

@XFAIL
def test_R8():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(k**2*binomial(n,k),(k,1,n))
    T = Sm.doit() #returns Piecewise function
    #T.simplify() raisesAttributeError: 'Or' object has no attribute 'as_numer_denom'
    assert T.combsimp() == n*(n+1)*2**(n-2)


def test_R9():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(binomial(n,k-1)/k,(k,1,n+1))
    assert Sm.doit().simplify() ==  (2**(n + 1) - 1)/(n + 1)

@XFAIL
def test_R10():
    n,m,r,k = symbols('n m r k', integer=True, positive=True)
    Sm = Sum(binomial(n,k)*binomial(m,r-k),(k,0,r))
    T = Sm.doit()
    T2 = T.combsimp().rewrite(factorial)
    assert T2 == factorial(m + n)/(factorial(r)*factorial(m + n - r))
    assert T2 == binomial(m+n,r).rewrite(factorial)
    T3 = T2.rewrite(binomial) # rewrite(binomial) is not working. https://code.google.com/p/sympy/issues/detail?id=4036
    assert T3 ==  binomial(m+n,r)

@XFAIL
def test_R11():
    n,k = symbols('n k', integer=True, positive=True)
    sk = binomial(n,k)*fibonacci(k)
    Sm = Sum(sk,(k,0,n))
    T = Sm.doit()
    # Fibonnaci simplification not implemented https://code.google.com/p/sympy/issues/detail?id=4035
    assert T == fibonacci(2*n)

@XFAIL
def test_R12():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(fibonacci(k)**2,(k,0,n))
    T = Sm.doit()
    assert T == fibonacci(n)*fibonacci(n+1)

@XFAIL
def test_R13():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(sin(k*x),(k,1,n))
    T = Sm.doit() # Sum is not calculated
    assert T.simplify() ==  cot(x/2)/2 - cos(x*(2*n + 1)/2)/(2*sin(x/2))

@XFAIL
def test_R14():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(sin((2*k-1)*x),(k,1,n))
    T = Sm.doit() # Sum is not calculated
    assert T.simplify() == sin(n*x)**2/sin(x)

@XFAIL
def test_R15():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(binomial(n - k, k), (k, 0, floor(n/2)))
    T = Sm.doit() # Sum is not calculated
    assert T.simplify() == fibonacci(n+1)

def test_R16():
    k = symbols('k', integer=True, positive=True)
    Sm = Sum(1/k**2 + 1/k**3, (k,1,oo))
    assert Sm.doit() == zeta(3) + pi**2/6

def test_R17():
    k = symbols('k', integer=True, positive=True)
    assert float(Sum(1/k**2 + 1/k**3, (k,1,oo))) == 2.8469909700078206

@XFAIL
def test_R18():
    k = symbols('k', integer=True, positive=True)
    Sm = Sum(1/(2**k*k**2), (k,1,oo))
    T = Sm.doit() # returns polylog(2, 1/2), seems particular value for 1/2 is not known. https://code.google.com/p/sympy/issues/detail?id=4033
    assert T.simplify() == -log(2)**2/2 + pi**2/12

@XFAIL
def test_R19():
    k = symbols('k', integer=True, positive=True)
    Sm = Sum(1/((3*k+1)*(3*k+2)*(3*k+3)), (k,0,oo))
    T = Sm.doit()
    assert T.simplify() == -log(3)/4 + sqrt(3)*pi/12 # fails, no simplification

@XFAIL
def test_R20():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(binomial(n,4*k) , (k,0,oo))
    T = Sm.doit()
    assert T.simplify() == 2**(n/2)*cos(pi*n/4)/2 + 2**(n - 1)/2 # fails, no simplification

@XFAIL
def test_R21():
    k = symbols('k', integer=True, positive=True)
    Sm = Sum(1/(sqrt(k*(k + 1)) * (sqrt(k) + sqrt(k + 1))) , (k,1,oo))
    T = Sm.doit() # Sum not calculated
    assert T.simplify() == 1

@XFAIL
def test_R22():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(Sum(binomial(n, k)*binomial(n - k, n - 2*k)*x**n*y**(n - 2*k),(k,0,floor(n/2))),(n,0,oo))
    # How to express constraint abs(x*y)<1?
    T = Sm.doit()
    # Correct answer unknown, not possible to provide assert
    assert T == False

@XFAIL
def test_R23():
    n,k = symbols('n k', integer=True, positive=True)
    Sm = Sum(Sum(factorial(n)/(factorial(k)**2*factorial(n - 2*k))*(x/y)**k*(x*y)**(n - k), (n,2*k,oo)), (k,0,oo))
    # How to express constraint abs(x*y)<1?
    T = Sm.doit() # Sum not calculated
    assert T == -1/sqrt(x**2*y**2 - 4*x**2 - 2*x*y + 1)

def test_R24():
    m,k = symbols('m k', integer=True, positive=True)
    Sm = Sum(Product(k/(2*k - 1), (k,1,m)), (m,2,oo))
    assert Sm.doit() == pi/2

def test_S1():
    k = symbols('k', integer=True, positive=True)
    Pr = Product(gamma(k/3), (k,1,8))
    assert Pr.doit().simplify() == 640*sqrt(3)*pi**3/6561

def test_S2():
    n,k = symbols('n k', integer=True, positive=True)
    assert Product(k,(k,1,n)).doit(),simplify() == factorial(n)

def test_S3():
    n,k = symbols('n k', integer=True, positive=True)
    assert Product(x**k,(k,1,n)).doit().simplify() == x**(n*(n + 1)/2)

def test_S4():
    n,k = symbols('n k', integer=True, positive=True)
    assert Product(1+1/k,(k,1,n-1)).doit().simplify() == n

def test_S5():
    n,k = symbols('n k', integer=True, positive=True)
    assert Product((2*k-1)/(2*k),(k,1,n)).doit().combsimp() == factorial(n-Rational(1,2))/(sqrt(pi)*factorial(n))

@XFAIL
def test_S6():
    n,k = symbols('n k', integer=True, positive=True)
    # Product raises Infinite recursion error. https://code.google.com/p/sympy/issues/detail?id=4034
    assert Product(x**2-2*x*cos(k*pi/n)+1, (k,1,n-1)).doit().simplify() == (x**(2*n)-1)/(x**2-1)

@XFAIL
def test_S7():
    k = symbols('k', integer=True, positive=True)
    Pr = Product((k**3-1)/(k**3+1), (k,2,oo))
    T = Pr.doit()
    assert T.simplify() == Rational(2,3) #T simplifies incorrectly to 0

@XFAIL
def test_S8():
    k = symbols('k', integer=True, positive=True)
    Pr = Product(1 - 1/(2*k)**2, (k,1,oo))
    T = Pr.doit() # returns nan https://code.google.com/p/sympy/issues/detail?id=4037
    assert T.simplify() == 2/pi

@XFAIL
def test_S9():
    k = symbols('k', integer=True, positive=True)
    Pr = Product(1 + (-1)**(k + 1)/(2*k - 1), (k, 1, oo))
    T = Pr.doit() # Product raises Infinite recursion error. https://code.google.com/p/sympy/issues/detail?id=4034
    assert T.simplify() == sqrt(2)

@XFAIL
def test_S10():
    k = symbols('k', integer=True, positive=True)
    Pr = Product((k*(k +  1) + 1 + I)/(k*(k + 1) + 1 - I), (k,0,oo))
    T = Pr.doit()
    assert T.simplify() == -1 # raises OverflowError  https://code.google.com/p/sympy/issues/detail?id=4038

def test_T1():
    assert limit((1 + 1/n)**n, n, oo) == E
    assert limit((1 - cos(x))/x**2, x, 0) == Rational(1,2)

def test_T2():
    assert limit((3**x + 5**x)**(1/x), x, oo) == 5

@XFAIL
def test_T3():
    assert limit(log(x)/(log(x) + sin(x)), x, oo) == 1 #raises PoleError

def test_T4():
    assert limit((exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1)))) - exp(x))/x, x, oo) == -exp(2)

def test_T5():
    assert  limit(x*log(x)*log(x*exp(x) - x**2)**2/log(log(x**2 + 2*exp(exp(3*x**3*log(x))))),x,oo) == Rational(1,3)

def test_T6():
    assert limit(1/n * factorial(n)**(1/n), n, oo) == exp(-1)

def test_T7():
    limit(1/n * gamma(n + 1)**(1/n), n, oo)

def test_T8():
    a,z = symbols('a z', real=True,positive=True)
    assert limit(gamma(z + a)/gamma(z)*exp(-a*log(z)), z, oo) == 1

@XFAIL
def test_T9():
    z,k = symbols('z k', real=True,positive=True)
    assert limit(hyper((1, k), (1,), z/k), k,oo) == exp(z) # raises NotImplementedError: Don't know how to calculate the mrv of '(1, k)'

@XFAIL
def test_T10():
    limit(zeta(x) - 1/(x - 1), x, 1)# raises PoleError shouldreturn euler-mascheroni constant

@XFAIL
def test_T11():
    n,k = symbols('n k', integer=True, positive=True)
    limit(n**x/(x*product((1 + x/k), (k, 1, n))),n,oo) == gamma(x) #raises NotImplementedError
