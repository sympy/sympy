""" Tests from Michael Wester's 1999 paper "Review of CAS mathematical
capabilities".

http://www.math.unm.edu/~wester/cas/book/Wester.pdf
See also http://math.unm.edu/~wester/cas_review.html for detailed output of each
tested system.
"""

from sympy import (Rational, symbols, factorial, sqrt, log, exp, oo, product,
    binomial, rf, pi, gamma, igcd, factorint, nsimplify, radsimp, combsimp,
    npartitions, totient, primerange, factor, simplify, gcd, resultant, expand,
    normal, I, trigsimp, tan, sin, cos, diff, nan, limit, EulerGamma, polygamma,
    bernoulli, assoc_legendre, Function, re, im, DiracDelta, chebyshevt, atan,
    sinh, cosh, Symbol, floor, ceiling)
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.utilities.pytest import XFAIL, skip
from sympy.mpmath import mpi, mpc
from sympy import mpmath

R = Rational
x,y,z = symbols('x','y','z')
i,j,k,l,m,n = symbols('i', 'j', 'k', 'l', 'm', 'n', integer=True)
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
    assert (factorint(factorial(50)) == [(2, 47), (3, 22), (5, 12), (7, 8),
        (11, 4), (13, 3), (17, 2), (19, 2), (23, 2), (29, 1), (31, 1), (37, 1),
        (41, 1), (43, 1), (47, 1)])

@XFAIL
def test_C3():
    raise NotImplementedError("10!! == 3840\n"
        " 9!! == 945")

# Base conversions; not really implemented by sympy
# Whatever. Take credit!
def test_C4():
    assert 0xABC == 2748

def test_C5():
    assert 123 == int('234', 7)

def test_C6():
    assert int('677', 8) == int('1BF', 16) == 447

@XFAIL
def test_C7():
    # You can evalf, but that's not really the point.
    assert log(32768, 8) == 5

@XFAIL
def test_C8():
    raise NotImplementedError("modular arithmetic:\n"
        "  5 ** -1 mod 7 == 3\n"
        "  5 ** -1 mod 6 == 5"
    )

def test_C9():
    assert igcd(igcd(1776, 1554), 5698) == 74

def test_C10():
    x = 0
    for n in range(2, 11):
        x += R(1, n)
    assert x == R(4861, 2520)

@XFAIL
def test_C11():
    raise NotImplementedError("evalf(1/4) == 0.142857 *repeating*")

def test_C12():
    assert R(7, 11) * R(22, 7) == 2

def test_C13():
    test = R(10, 7) * (1 + R(29, 1000)) ** R(1,3)
    good = 3 ** R(1,3)
    assert test == good

def test_C14():
    assert nsimplify(sqrt(2*sqrt(3) + 4)) == 1 + sqrt(3)

def test_C15():
    test = nsimplify(sqrt(14 + 3*sqrt(3 + 2*sqrt(5 - 12*sqrt(3 - 2*sqrt(2))))))
    good = sqrt(2) + 3
    assert test == good

@XFAIL
def test_C16():
    test = radsimp(nsimplify(sqrt(10 + 2*sqrt(6) + 2*sqrt(10) + 2*sqrt(15))))
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
    assert radsimp(nsimplify((90 + 35*sqrt(7)) ** R(1,3))) == 3 + sqrt(7)

def test_C20():
    inside = (135 + 78*sqrt(3))
    test = nsimplify((inside**R(2,3) + 3) * sqrt(3) / inside**R(1,3))
    assert test == 12

def test_C21():
    assert nsimplify((41 + 29*sqrt(2)) ** R(1,5)) == 1 + sqrt(2)

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
    raise NotImplementedError("apply Horner's rule to sum(a[i]*x**i, (i,1,5))")

@XFAIL
def test_D9():
    raise NotImplementedError("translate D8 to FORTRAN")

@XFAIL
def test_D10():
    raise NotImplementedError("translate D8 to C")

@XFAIL
def test_D11():
    raise NotImplementedError("flops(sum(product(f[i][k], (i,1,k)), (k,1,n)))")

@XFAIL
def test_D12():
    assert (mpi(-4,2) * x + mpi(1,3)) ** 2 == mpi(-8, 16)*x**2 + mpi(-24, 12)*x + mpi(1, 9)

@XFAIL
def test_D13():
    raise NotImplementedError("discretize a PDE: diff(f(x,t),t) == diff(diff(f(x,t),x),x)")

# E. Statistics
#   See scipy; all of this is numerical.

# F. Combinatorial Theory.

def test_F1():
    assert rf(x, 3) == x*(1+x)*(2+x)

def test_F2():
    assert binomial(n, 3) == n*(1-n)*(2-n)/6

@XFAIL
def test_F3():
    raise NotImplementedError("2**n * n! * (2*n-1)!! == (2*n)! or gamma(2*n+1)")

@XFAIL
def test_F4():
    assert combsimp((2**n * factorial(n) * product(2*k-1, (k,1,n)))) == factorial(2*n)

@XFAIL
def test_F5():
    assert gamma(n+R(1,2)) / sqrt(pi) / factorial(n) == factorial(2*n)/2**(2*n)/factorial(n)**2

@XFAIL
def test_F6():
    raise NotImplementedError("find the partitions of 4: [4,2+2,1+3,1+1+2,1+1+1+1]")

def test_F7():
    assert npartitions(4) == 5

@XFAIL
def test_F8():
    raise NotImplementedError("S1(5,2) == -50; Stirling numbers")

def test_F9():
    assert totient(1776) == 576

# G. Number Theory

def test_G1():
    assert list(primerange(999983, 1000004)) ==  [999983, 1000003]

@XFAIL
def test_G2():
    raise NotImplementedError("find the primitive root of 191 == 19")

@XFAIL
def test_G3():
    raise NotImplementedError("(a+b)**p mod p == a**p + b**p mod p; p prime")

# ... G20 I don't think these are implemented.

# H. Algebra

@XFAIL
def test_H1():
    assert 2 * 2**n == 2 ** (n+1)

@XFAIL
def test_H2():
    assert 4 * 2**n == 2 ** (n+2)

@XFAIL
def test_H3():
    assert (-1) ** (n*(n+1)) == 1

@XFAIL
def test_H4():
    expr = factor(6*x - 10)
    assert type(expr) is Mul
    assert expr.args[0] == 2
    assert factor(6*x - 10) == 2 * (3*x-5)

p1 = 64*x**34 - 21*x**47 - 126*x**8 - 46*x**5 - 16*x**60 - 81
p2 = 72*x**60 - 25*x**25 - 19*x**23 - 22*x**39 - 83*x**52 + 54*x**10 + 81
q = 34*x**19 - 25*x**16 + 70*x**7 + 20*x**3 - 91*x - 86

def test_H5():
    assert gcd(p1, p2, x) == 1

@XFAIL
def test_H6():
    assert gcd(expand(p1 * q), expand(p2 * q), x) == q

@XFAIL
def test_H7():
    skip('takes too much time')
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    assert gcd(p1, p2, x, y, z) == 1

def test_H8():
    skip('takes too much time')
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    q = 11*x**12*y**7*z**13 - 23*x**2*y**8*z**10 + 47*x**17*y**5*z**8
    assert gcd(p1 * q, p2 * q, x, y, z) == q

@XFAIL
def test_H9():
    p1 = 2*x**(n+4) - x**(n+2)
    p2 = 4*x**(n+1) + 3*x**n
    assert gcd(p1, p2, x) == x**n

def test_H10():
    p1 = 3*x**4 + 3*x**3 + x**2 - x - 2
    p2 = x**3 - 3*x**2 + x + 5
    assert resultant(p1, p2, x) == 0

def test_H11():
    skip('takes too much time')
    assert resultant(p1 * q, p2 * q, x) == 0

def test_H12():
    num = x**2 - 4
    den = x**2 + 4*x + 4
    assert normal(num, den, x) == (x-2, x+2, 1)

@XFAIL
def test_H13():
    assert simplify((exp(x) - 1) / (exp(x/2) + 1)) == exp(x/2) - 1

def test_H14():
    p = (x+1) ** 20
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

@XFAIL
def test_H15():
    # You can use solve() to do the factorization, but putting it back is hard.
    # The coefficients are also much nastier than, say, Maple provides, even
    # after nsimplify.
    raise NotImplementedError("factor(x**3 + x**2 - 7) then expand it back to the original form")

def test_H16():
    assert factor(x**100 - 1) == (-(1 + x)*(1 + x**2)*(1 - x)*(1 + x + x**2
        + x**3 + x**4)*(1 - x + x**2 - x**3 + x**4)*(1 + x**5 + x**10 + x**15
            + x**20)*(1 - x**5 + x**10 - x**15 + x**20)*(1 - x**10 + x**20
                - x**30 + x**40)*(1 - x**2 + x**4 - x**6 + x**8))

# Takes too long.
@XFAIL
def test_H17():
    skip('takes too much time')
    assert factor(expand(p1 * p2)) == p1 * p2

@XFAIL
def test_H18():
    # Factor over complex rationals.
    test = factor(4*x**4 + 8*x**3 + 77*x**2 + 18*x + 53)
    good = (2*x + 3*I)*(2*x - 3*I)*(x + 1 - 4*I) (x + 1 + 4*I)
    assert test == good

@XFAIL
def test_H19():
    raise NotImplementedError("let a**2==2; 1/(a-1) == a+1")

@XFAIL
def test_H20():
    raise NotImplementedError("let a**2==2; (x**3 + (a-2)*x**2 - (2*a+3)*x - 3*a) / (x**2-2) = (x**2 - 2*x - 3) / (x-a)")

@XFAIL
def test_H21():
    raise NotImplementedError("let b**3==2, c**2==3; evaluate (b+c)**4")

@XFAIL
def test_H22():
    raise NotImplementedError("factor x**4 - 3*x**2 + 1 mod 5")

@XFAIL
def test_H23():
    raise NotImplementedError("factor x**11 + x + 1 mod 65537")

@XFAIL
def test_H24():
    raise NotImplementedError("factor x**4 - 3*x**2 + 1, GoldenRatio")

@XFAIL
def test_H25():
    e = (x - 2*y**2 + 3*z**3) ** 20
    assert factor(expand(e)) == e

@XFAIL
def test_H26():
    raise NotImplementedError("factor(expand((sin(x) - 2*cos(y)**2 + 3*tan(z)**3)**20))")

# ... I'm bored with XFAILs. Goes up to H33.

# I. Trigonometry

@XFAIL
def test_I1():
    # I'm not sure that nsimplify is the correct thing to try. There ought to be
    # a way to get to the target expression exactly.
    assert tan(7*pi/10) == -sqrt(1+2/sqrt(5))

@XFAIL
def test_I2():
    assert sqrt((1+cos(6))/2) == -cos(3)

@XFAIL
def test_I3():
    assert cos(n*pi) + sin((4*n-1)*pi/2) == (-1)**n - 1

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

@XFAIL
def test_I11():
    assert limit((tan(x)**2 + 1 - cos(x)**-2) / (sin(x)**2 + cos(x)**2 - 1), x, 0) != 0

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
    assert gamma(R(-1,2)) == -2*sqrt(pi)

@XFAIL
def test_J5():
    assert polygamma(0, R(1,3)) == -EulerGamma - pi/2*sqrt(R(1,3)) - R(3,2)*log(3)

def test_J6():
    assert mpmath.jn(2, 1+1j).ae(mpc('0.04157988694396212', '0.24739764151330632'))

@XFAIL
def test_J7():
    raise NotImplementedError("jv(R(-5,2), pi/2) == 12/(pi**2)")

@XFAIL
def test_J8():
    raise NotImplementedError("jv(R(3,2), z) == sqrt(2/(pi*z))*(sin(z)/z - cos(z))")

@XFAIL
def test_J9():
    raise NotImplementedError("diff(j0(z), z) == -j1(z)")

@XFAIL
def test_J10():
    mu, nu = symbols('mu', 'nu', integer=True)
    assert assoc_legendre(nu, mu, 0) == 2**mu*sqrt(pi)/gamma((nu-mu)/2+1)/gamma((-nu-mu+1)/2)

def test_J11():
    assert assoc_legendre(3,1,x) == sqrt(1 - x**2)*(R(3,2) - R(15,2)*x**2)

def test_J12():
    skip('takes too much time')
    assert simplify(chebyshevt(1008,x) - 2*x*chebyshevt(1007,x) + chebyshevt(1006,x)) == 0

@XFAIL
def test_J13():
    a = Symbol("a", integer=True, negative=False)
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
    assert deltaintegrate(f((x+2)/5)*DiracDelta((x-2)/3) - g(x)*diff(DiracDelta(x-1),x), (x,0,3))

@XFAIL
def test_J18():
    raise NotImplementedError("define an antisymmetric function")


# K. The Complex Domain

def test_K1():
    z1, z2 = symbols('z1', 'z2', complex=True)
    assert re(z1+I*z2) == -im(z2) + re(z1)
    assert im(z1+I*z2) ==  im(z1) + re(z2)

def test_K2():
    assert abs(3 - sqrt(7) + I*sqrt(6*sqrt(7)-15)) == 1

@XFAIL
def test_K3():
    a, b = symbols('a', 'b', real=True)
    assert simplify(abs(1/(a+I/a+I*b))) == 1/sqrt(a**2 + (I/a+b)**2)

def test_K4():
    assert log(3+4*I).expand(complex=True) == log(5) + I*atan(R(4,3))

def test_K5():
    x, y = symbols('x', 'y', real=True)
    assert tan(x+I*y).expand(complex=True) == sin(x)*cos(x) / (cos(x)**2 + sinh(y)**2) + I*sinh(y)*cosh(y) / (cos(x)**2 + sinh(y)**2)

@XFAIL
def test_K6():
    expr = sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z))
    sexpr = simplify(expr)
    assert sexpr == sqrt(x*y)/sqrt(x)
    assert sexpr != sqrt(y)

def test_K7():
    y = Symbol('y', negative=False)
    expr = sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z))
    sexpr = simplify(expr)
    assert sexpr == sqrt(y)

def test_K8():
    z = Symbol('z', complex=True)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) != 0
    z = Symbol('z', complex=True, negative=False)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0

def test_K9():
    z = Symbol('z', real=True, positive=True)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0

def test_K10():
    z = Symbol('z', real=True, negative=True)
    assert simplify(sqrt(1/z) + 1/sqrt(z)) == 0


# L. Determining Zero Equivalence

def test_L1():
    assert sqrt(997)-(997**3)**R(1,6) == 0

@XFAIL
def test_L2():
    assert sqrt(999983)-(999983**3)**R(1,6) == 0

def test_L3():
    assert simplify((2**R(1,3)+4**R(1,3))**3-6*(2**R(1,3)+4**R(1,3))-6) == 0

def test_L4():
    assert trigsimp(cos(x)**3+cos(x)*sin(x)**2-cos(x)) == 0

@XFAIL
def test_L5():
    assert log(tan(R(1,2)*x+pi/4))-asinh(tan(x)) == 0

def test_L6():
    assert (log(tan(x/2+pi/4))-asinh(tan(x))).diff(x).subs({x:0}) == 0

@XFAIL
def test_L7():
    assert simplify(log((2*sqrt(x)+1)/(sqrt(4*x+4*sqrt(x)+1)))) == 0

@XFAIL
def test_L8():
    assert (4*x+4*sqrt(x)+1)**(sqrt(x)/(2*sqrt(x)+1))*(2*sqrt(x)+1)**(1/(2*sqrt(x)+1))-2*sqrt(x)-1 == 0

@XFAIL
def test_L9():
    z = symbols('z', complex=True)
    assert 2**(1-z)*gamma(z)*zeta(z)*cos(z*pi/2)-pi**2*zeta(1-z) == 0

# M. Equations

@XFAIL
def test_M1():
    assert Equality(x,2)/2 + Equality(1,1) == Equality(x/2+1,2)

def test_M2():
    # This takes a bit or work, but SymPy is capable of recognizing that all the
    # roots of this equation are real.
    sol = solve(3*x**3-18*x**2+33*x-19,x)
    for i in sol:
        assert re(i.as_real_imag()) == i.as_real.imag()

