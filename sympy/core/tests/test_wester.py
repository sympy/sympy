""" Tests from Michael Wester's 1999 paper "Review of CAS mathematical
capabilities".

http://www.math.unm.edu/~wester/cas/book/Wester.pdf
See also http://math.unm.edu/~wester/cas_review.html for detailed output of each
tested system.
"""

from sympy import (Rational, symbols, factorial, sqrt, log, exp, oo, product,
    binomial, rf, pi, gamma, igcd, factorint, nsimplify, radsimp, combsimp,
    npartitions, totient, primerange, factor, simplify, gcd, resultant, expand,
    I, trigsimp, tan, sin, cos, diff, nan, limit, EulerGamma, polygamma,
    bernoulli, assoc_legendre, Function, re, im, DiracDelta, chebyshevt, atan,
    sinh, cosh, floor, ceiling, solve, asinh, LambertW, N, apart, sqrtdenest,
    factorial2, powdenest, Mul, S, mpmath, ZZ, Poly, expand_func)

from sympy.functions.combinatorial.numbers import stirling
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.utilities.pytest import XFAIL, slow
from sympy.utilities.iterables import partitions
from sympy.mpmath import mpi, mpc
from sympy.physics.quantum import Commutator

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
    assert nsimplify(sqrt(2*sqrt(3) + 4)) == 1 + sqrt(3)


def test_C15():
    test = nsimplify(sqrt(14 + 3*sqrt(3 + 2*sqrt(5 - 12*sqrt(3 - 2*sqrt(2))))))
    good = sqrt(2) + 3
    assert test == good


def test_C16():
    test = sqrtdenest(sqrt(10 + 2*sqrt(6) + 2*sqrt(10) + 2*sqrt(15)))
    good = sqrt(2) + sqrt(3) + sqrt(5)
    assert test == good


def test_C17():
    test = nsimplify((sqrt(3) + sqrt(2)) / (sqrt(3) - sqrt(2)))
    good = 5 + 2*sqrt(6)
    assert test == good


def test_C18():
    assert nsimplify(sqrt(-2 + sqrt(-5)) * sqrt(-2 - sqrt(-5))) == 3


@XFAIL
def test_C19():
    assert radsimp(nsimplify((90 + 35*sqrt(7)) ** R(1, 3))) == 3 + sqrt(7)


def test_C20():
    inside = (135 + 78*sqrt(3))
    test = nsimplify((inside**R(2, 3) + 3) * sqrt(3) / inside**R(1, 3))
    assert test == 12


def test_C21():
    assert nsimplify((41 + 29*sqrt(2)) ** R(1, 5)) == 1 + sqrt(2)


@XFAIL
def test_C22():
    test = nsimplify(((6 - 4*sqrt(2))*log(3 - 2*sqrt(2)) + (3 - 2*sqrt(2))*log(17
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


@XFAIL
def test_H2():
    assert 4 * 2**n == 2 ** (n + 2)


@XFAIL
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


@XFAIL
def test_H24():
    raise NotImplementedError("factor x**4 - 3*x**2 + 1, GoldenRatio")


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


@XFAIL
def test_I3():
    assert cos(n*pi) + sin((4*n - 1)*pi/2) == (-1)**n - 1


@XFAIL
def test_I4():
    assert cos(pi*cos(n*pi)) + sin(pi/2*cos(n*pi)) == (-1)**n - 1


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


@XFAIL
def test_J2():
    raise NotImplementedError("diff(E(phi,k), k) == (E(phi,k) - F(phi,k)) / k; F() and E() are elliptic integrals of the 1st and 2nd kind, respectively")


@XFAIL
def test_J3():
    raise NotImplementedError("Jacobi elliptic functions: diff(dn(u,k), u) == -k**2*sn(u,k)*cn(u,k)")


def test_J4():
    assert gamma(R(-1, 2)) == -2*sqrt(pi)


def test_J5():
    assert polygamma(0, R(1, 3)) == -EulerGamma - pi/2*sqrt(R(1, 3)) - R(3, 2)*log(3)


def test_J6():
    assert mpmath.besselj(2, 1 + 1j).ae(mpc('0.04157988694396212', '0.24739764151330632'))


@XFAIL
def test_J7():
    raise NotImplementedError("jv(R(-5,2), pi/2) == 12/(pi**2)")


@XFAIL
def test_J8():
    raise NotImplementedError("jv(R(3,2), z) == sqrt(2/(pi*z))*(sin(z)/z - cos(z))")


@XFAIL
def test_J9():
    raise NotImplementedError("diff(j0(z), z) == -j1(z)")


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


@XFAIL
def test_J14():
    raise NotImplementedError("F(R(1,2),R(1,2),R(3,2),z**2) == asin(z)/z; F(.) is hypergeometric function")


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
    assert (4*x + 4*sqrt(x) + 1)**(sqrt(x)/(2*sqrt(x) + 1)) \
        *(2*sqrt(x) + 1)**(1/(2*sqrt(x) + 1)) - 2*sqrt(x) - 1 == 0


@XFAIL
def test_L9():
    z = symbols('z', complex=True)
    assert 2**(1 - z)*gamma(z)*zeta(z)*cos(z*pi/2) - pi**2*zeta(1 - z) == 0

# M. Equations


@XFAIL
def test_M1():
    assert Equality(x, 2)/2 + Equality(1, 1) == Equality(x/2 + 1, 2)


@XFAIL
def test_M2():
    # The roots of this equation should all be real. Note that this doesn't test
    # that they are correct.
    sol = solve(3*x**3 - 18*x**2 + 33*x - 19, x)
    for i in sol:
        assert im(i) == 0


@XFAIL
def test_M5():
    assert solve(x**6 - 9*x**4 - 4*x**3 + 27*x**2 - 36*x - 23, x) == [2**(1/3) + sqrt(3), 2**(1/3) - sqrt(3), +sqrt(3) - 1/2**(2/3) + I*sqrt(3)/2**(2/3), +sqrt(3) - 1/2**(2/3) - I*sqrt(3)/2**(2/3), -sqrt(3) - 1/2**(2/3) + I*sqrt(3)/2**(2/3), -sqrt(3) - 1/2**(2/3) - I*sqrt(3)/2**(2/3)]


def test_M6():
    assert set(solve(x**7 - 1, x)) == set([cos(n*2*pi/7) + I*sin(n*2*pi/7) for n in range(0, 7)])
    # The paper asks for exp terms, but sin's and cos's may be acceptable


def test_M7():
    assert set(solve(x**8 - 8*x**7 + 34*x**6 - 92*x**5 + 175*x**4 - 236*x**3 +
        226*x**2 - 140*x + 46, x)) == set(
            [1 + sqrt(-6 + 2*sqrt(-3 - 4*sqrt(3)))/2,
        1 + sqrt(-6 - 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 + sqrt(-6 + 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 - 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 + 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 - 2*sqrt(-3 - 4*sqrt(3)))/2,
        1 - sqrt(-6 + 2*sqrt(-3 - 4*sqrt(3)))/2,
        1 + sqrt(-6 - 2*sqrt(-3 - 4*sqrt(3)))/2])


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


@XFAIL
def test_M12():
    raise NotImplementedError("solve((x+1)*(sin(x)**2+1)**2*cos(3*x)**3,x)")


@XFAIL
def test_M13():
    raise NotImplementedError("solve(sin(x)-cos(x),x)")


def test_M14():
    assert solve(tan(x) - 1, x) == [pi/4]


@XFAIL
def test_M15():
    assert solve(sin(x) - 1/2) == [pi/6, 5*pi/6]


@XFAIL
def test_M16():
    raise NotImplementedError("solve(sin(x)-tan(x),x)")


@XFAIL
def test_M17():
    raise NotImplementedError("solve(asin(x)-atan(x),x)")


@XFAIL
def test_M18():
    raise NotImplementedError("solve(acos(x)-atan(x),x)")
    #assert solve(acos(x) - atan(x), x) == [sqrt((sqrt(5) - 1)/2)]


def test_M19():
    assert solve((x - 2)/x**R(1, 3), x) == [2]


def test_M20():
    assert solve(sqrt(x**2 + 1) - x + 2, x) == []
