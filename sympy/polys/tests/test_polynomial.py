
from sympy import symbols, expand, sin, sqrt, re, im, I, Rational

from sympy.polys.monomial import *
from sympy.polys.polynomial import *
from sympy.polys.algorithms import *

import py

a,b,c,x,y,z,u,v,t = symbols('abcxyzuvt')

def test_monomial_cmp():
    assert monomial_lex_cmp((3,2,1), (1,2,4)) == 1
    assert monomial_grlex_cmp((2,4,1), (1,6,0)) == 1
    assert monomial_grevlex_cmp((1,3,1), (1,2,2)) == 1
    assert monomial_1_el_cmp((2,0,1), (1,2,0)) == 1

def test_poly_basics():
    f = x*y**2*z + x*y + x + 1

    assert Poly(f, x).as_basic().expand() == f
    assert Poly(f, x, y).as_basic().expand() == f
    assert Poly(f, x, y, z).as_basic().expand() == f

    assert Poly((), x, y, z).as_basic() == 0
    assert Poly([], x, y, z).as_basic() == 0
    assert Poly({}, x, y, z).as_basic() == 0

    assert Poly(((), ()), x, y, z).as_basic() == 0
    assert Poly(([], []), x, y, z).as_basic() == 0

    assert Poly((Integer(17), (2,5,4)),
        x, y, z).as_basic() == 17*x**2*y**5*z**4
    assert Poly(((Integer(17),), ((2,5,4),)),
        x, y, z).as_basic() == 17*x**2*y**5*z**4

    p = Poly(x**2 + x + 1, x)

    assert Poly(p, y) == \
        Poly(((1 + x + x**2,), ((0,),)), y, order='grlex')

    assert Poly(p, x, y) == \
        Poly(((1, 1, 1), ((2, 0), (1, 0), (0, 0))), x, y, order='grlex')

    assert Poly(p, y, x) == \
        Poly(((1, 1, 1), ((0, 2), (0, 1), (0, 0))), y, x, order='grlex')

    assert Poly(p, z, x, y) == \
        Poly(((1, 1, 1), ((0, 2, 0), (0, 1, 0), (0, 0, 0))), z, x, y, order='grlex')

    p = Poly(x*y**3 + x**2*y + 1, x)

    assert Poly(p, y) == \
        Poly(((x, x**2, 1), ((3,), (1,), (0,))), y, order='grlex')

    assert Poly(p, x, y) == \
        Poly(((1, 1, 1), ((1, 3), (2, 1), (0, 0))), x, y, order='grlex')

    assert Poly(p, y, x) == \
        Poly(((1, 1, 1), ((3, 1), (1, 2), (0, 0))), y, x, order='grlex')

    assert Poly(p, x, z) == \
        Poly(((y, y**3, 1), ((2, 0), (1, 0), (0, 0))), x, z, order='grlex')

    assert Poly(p, z, x) == \
        Poly(((y, y**3, 1), ((0, 2), (0, 1), (0, 0))), z, x, order='grlex')

    assert Poly(p, z, x, y) == \
        Poly(((1, 1, 1), ((0, 1, 3), (0, 2, 1), (0, 0, 0))), z, x, y, order='grlex')

    assert Poly(p, x, z, y) == \
        Poly(((1, 1, 1), ((1, 0, 3), (2, 0, 1), (0, 0, 0))), x, z, y, order='grlex')

    f = x*t**4 + y*t**2 + z

    assert Poly([(x, (4,)), (y, (2,)), (z, (0,))], t) == Poly(f, t)
    assert Poly([(x, 4), (y, 2), (z, 0)], t) == Poly(f, t)

    f = x*t**4 + z

    assert Poly([(x, 4), (0, 2), (z, 0)], t) == Poly(f, t)

    f = x*t**4 + y*t**2 + z

    assert Poly({(4,): x, (2,): y, (0,): z}, t) == Poly(f, t)
    assert Poly({4: x, 2: y, 0: z}, t) == Poly(f, t)

def test_poly_internals():
    p = Poly(x**2*y*z + x*y*z**3 + x*y + y*z, x, y, z)

    assert p.as_dict() == \
        {(1, 1, 3): 1, (1, 1, 0): 1, (2, 1, 1): 1, (0, 1, 1): 1}

    assert Poly._permute(p, x) == \
        {(2,): y*z, (0,): y*z, (1,): y + y*z**3}

    assert Poly._permute(p, y) == \
        {(1,): x + z + x*z**3 + z*x**2}

    assert Poly._permute(p, z) == \
        {(0,): x*y, (3,): x*y, (1,): y + y*x**2}

    assert Poly._permute(p, x, y) == \
        {(0, 1): z, (1, 1): 1 + z**3, (2, 1): z}

    assert Poly._permute(p, y, x) == \
        {(1, 2): z, (1, 0): z, (1, 1): 1 + z**3}

    assert Poly._permute(p, x, z) == \
        {(0, 1): y, (1, 0): y, (1, 3): y, (2, 1): y}

    assert Poly._permute(p, z, x) == \
        {(1, 2): y, (0, 1): y, (1, 0): y, (3, 1): y}

    assert Poly._permute(p, y, z) == \
        {(1, 0): x, (1, 3): x, (1, 1): 1 + x**2}

    assert Poly._permute(p, z, y) == \
        {(0, 1): x, (3, 1): x, (1, 1): 1 + x**2}

    q = Poly(x**2*y*z + 2*x*y*z**3 + 3*x*y + 4*y*z, x, y, z)

    assert q.as_dict() == \
        {(1, 1, 3): 2, (1, 1, 0): 3, (2, 1, 1): 1, (0, 1, 1): 4}

    assert Poly._permute(q, z, y, x) == \
        {(0, 1, 1): 3, (1, 1, 0): 4, (3, 1, 1): 2, (1, 1, 2): 1}

def test_poly_cancel():
    assert Poly._cancel(x) == x
    assert Poly._cancel(x+1) == x+1
    assert Poly._cancel((x+1)/(x-1)) == (x+1)/(x-1)

    assert Poly._cancel((x**2-1)/(x-1)) == x+1
    assert Poly._cancel((x**2-y**2)/(x-y)) == x+y

    assert Poly._cancel((x**2-y)/(x-y)) == (x**2 - y)/(x - y)
    assert Poly._cancel((x**2-2)/(x+sqrt(2))) == x - sqrt(2)

def test_poly_characteristics():
    f = -3*x**5*y*z**4 + 2*x**2*y**8 - x*y**4 + x*y*z**3

    p = Poly(f, x, y, z, order='lex')

    assert p.total_degree == 10
    assert p.degree == 10
    assert p.length == 4
    assert p.norm == 3

    assert p.LT == (-3, (5, 1, 4))
    assert p.LM == (5, 1, 4)
    assert p.LC == -3

    assert p.coeffs == (-3,2,-1,1)
    assert p.monoms == ((5,1,4), (2,8,0), (1,4,0), (1,1,3))

    q = Poly(f, x, y, z, order='grlex')

    assert q.total_degree == 10
    assert q.degree == 10
    assert q.length == 4
    assert q.norm == 3

    assert q.LT == (-3, (5, 1, 4))
    assert q.LM == (5, 1, 4)
    assert q.LC == -3

    assert q.coeffs == (-3,2,-1,1)
    assert q.monoms == ((5,1,4), (2,8,0), (1,4,0), (1,1,3))

    r = Poly(f, x, y, z, order='grevlex')

    assert r.total_degree == 10
    assert r.degree == 10
    assert r.length == 4
    assert r.norm == 3

    assert r.LT == (2, (2, 8, 0))
    assert r.LM == (2, 8, 0)
    assert r.LC == 2

    assert r.coeffs == (2,-3,-1,1)
    assert r.monoms == ((2,8,0), (5,1,4), (1,4,0), (1,1,3))

def test_poly_properties():
    f = Poly(0, x)

    assert f.is_zero == True
    assert f.is_one == False
    assert f.is_number == True
    assert f.is_constant == True
    assert f.is_monomial == True
    assert f.is_univariate == True
    assert f.is_multivariate == False
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

    f = Poly(0, x, y, z)

    assert f.is_zero == True
    assert f.is_one == False
    assert f.is_number == True
    assert f.is_constant == True
    assert f.is_monomial == True
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

    f = Poly(1, x)

    assert f.is_zero == False
    assert f.is_one == True
    assert f.is_number == True
    assert f.is_constant == True
    assert f.is_monomial == True
    assert f.is_univariate == True
    assert f.is_multivariate == False
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == True

    f = Poly(43242, x, y, z)

    assert f.is_zero == False
    assert f.is_one == False
    assert f.is_number == True
    assert f.is_constant == True
    assert f.is_monomial == True
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

    f = Poly(x, x, y, z)

    assert f.is_zero == False
    assert f.is_one == False
    assert f.is_number == False
    assert f.is_constant == False
    assert f.is_monomial == True
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == True
    assert f.is_inhomogeneous == False
    assert f.is_monic == True

    f = Poly(t, x, y, z)

    assert f.is_zero == False
    assert f.is_one == False
    assert f.is_number == False
    assert f.is_constant == True
    assert f.is_monomial == True
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

    f = Poly(2*x*y*z + x*y + 1, x, y, z)

    assert f.is_zero == False
    assert f.is_one == False
    assert f.is_number == False
    assert f.is_constant == False
    assert f.is_monomial == False
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

def test_as_monic():
    assert Poly(0, x,y,z).as_monic() == Poly(0, x,y,z)
    assert Poly(1, x,y,z).as_monic() == Poly(1, x,y,z)
    assert Poly(t, x,y,z).as_monic() == Poly(1, x,y,z)

    f = x**3*y*z**2 + 7*x*y*z + 14*x*y + 3

    assert Poly(f, x,y,z).as_monic() == Poly(f, x,y,z)

    f = 7*x**3*y*z**2 + 7*x*y*z + 14*x*y + 3
    g = x**3*y*z**2 + x*y*z + 2*x*y + Rational(3,7)

    assert Poly(f, x,y,z).as_monic() == Poly(g, x,y,z)

    f = (z-1)*x**3*y + (z**2-1)*x*y + (z-1)**2*x + 7
    g = x**3*y + (z+1)*x*y + (z-1)*x + 7/(z-1)

    assert Poly(f, x,y).as_monic() == Poly(g, x,y)

    f = y*x**2 + y**2 + y
    g = x**2 + y + 1

    assert Poly(f, x).as_monic() == Poly(g, x)

def test_poly_add():
    f = -Rational(1,6)*x**2-Rational(5,36)+Rational(17,18)
    g = -Rational(1,6)*x**2-Rational(5,36)+Rational(35,18)

    assert Poly(1, x) + Poly(f, x) == Poly(g, x)

def test_poly_sub():
    f = -Rational(1,6)*x**2-Rational(5,36)+Rational(17,18)
    g = Rational(1,6)*x**2+Rational(5,36)+Rational(1,18)

    assert Poly(1, x) - Poly(f, x) == Poly(g, x)

def test_poly_mul():
    f = x**3-12*x**2-42
    g = x**2+x-3

    assert Poly(f, x)*Poly(g, x) == f*g
    assert Poly(f, x)*g == f*g

    f = x**3*y-2*x*y**2-3*z+1
    g = x**2*y*z+x*y**3-3

    assert Poly(f, x, y, z)*Poly(g, x, y, z) == f*g
    assert Poly(f, x, y, z)*g == f*g

    p = Poly(f, x, z)*Poly(g, y, z)

    assert p.symbols == (x, y, z)

def test_poly_div():
    f = Poly(x**3-12*x**2-42, x)
    g = Poly(x**2+x-3, x)

    h = poly_div(f, g)

    assert h[0] == Poly(x-13, x)
    assert h[1] == Poly(16*x-81, x)

    assert h[0].as_basic() == x-13
    assert h[1].as_basic() == 16*x-81

    assert f / g == Poly(x-13, x)
    assert f % g == Poly(16*x-81, x)

    assert divmod(f, g) == (Poly(x-13, x), Poly(16*x-81, x))

    assert f / x == Poly(x**2-12*x, x)
    assert f % x == Poly(-42, x)

    assert divmod(f, x) == (Poly(x**2-12*x, x), Poly(-42, x))

    assert f / sin(x) == f.as_basic() / sin(x)
    assert f % sin(x) == 0

    assert divmod(f, sin(x)) == (f.as_basic() / sin(x), 0)

    assert poly_div(4*x**2*y-2*x*y-2*y+4*x+8, 2, x, y) == \
        (Poly(2*x**2*y-x*y-y+2*x+4, x, y), Poly((), x, y))

    assert poly_div(4*x**2*y-2*x*y-2*y+4*x+8, 2*y, x, y) == \
        (Poly(2*x**2-x-1, x, y), Poly(4*x+8, x, y))

    assert poly_div(x-1, y-1, x, y) == (Poly((), x, y), Poly(x-1, x, y))

    assert poly_div(x**3-12*x**2-42, x-3, x) == \
        (Poly(x**2-9*x-27, x), Poly(-123, x))

    assert poly_div(2+2*x+x**2, 1, x) == (Poly(2+2*x+x**2, x), Poly(0, x))
    assert poly_div(2+2*x+x**2, 2, x) == (Poly(1+x+x**2/2, x), Poly(0, x))

    assert poly_div(3*x**3, x**2, x) == (Poly(3*x, x), Poly(0, x))

    assert poly_div(1, x, x) == (Poly(0, x), Poly(1, x))

    assert poly_div(x*y+2*x+y,x,x) == (Poly(2+y, x), Poly(y, x))

    assert poly_div(x*y**2 + 1, [x*y+1, y+1], x, y) == \
        ((Poly(y, x, y), Poly(-1, x, y)), Poly(2, x, y))

    assert poly_div(x**2*y+x*y**2+y**2, [x*y-1, y**2-1], x, y) == \
        ((Poly(x+y, x, y), Poly(1, x, y)), Poly(1+x+y, x, y))
    assert poly_div(x**2*y+x*y**2+y**2, [y**2-1, x*y-1], x, y) == \
        ((Poly(1+x, x, y), Poly(x, x, y)), Poly(1+2*x, x, y))

    f, g = 3*x**3 + x**2 + x + 5, 5*x**2 - 3*x + 1

    q = Poly(Rational(3,5)*x + Rational(14, 25), x)
    r = Poly(Rational(52, 25)*x + Rational(111, 25), x)

    assert poly_div(f, g, x) == (q, r)
    assert poly_div(Poly(f, x), Poly(g, x)) == (q, r)

    q = Poly(15*x + 14, x)
    r = Poly(52*x + 111, x)

    assert poly_pdiv(f, g, x) == (q, r)
    assert poly_pdiv(Poly(f, x), Poly(g, x)) == (q, r)

def test_poly_lcm():
    assert poly_lcm(2, 6, x) == Poly(6, x)
    assert poly_lcm(2, 6, x, y) == Poly(6, x, y)

    assert poly_lcm(x, y, x, y) == Poly(x*y, x, y)

    assert poly_lcm(2*x**3, 6*x, x) == Poly(6*x**3, x)
    assert poly_lcm(2*x**3, 3*x, x) == Poly(6*x**3, x)

    assert poly_lcm(2*x**3, 6*x*y**2, x, y) == Poly(6*x**3*y**2, x, y)
    assert poly_lcm(2*x**3, 3*x*y**2, x, y) == Poly(6*x**3*y**2, x, y)

    assert poly_lcm(x**2+x, x, x) == Poly(x**2+x, x)
    assert poly_lcm(x**2+x, 2*x, x) == Poly(2*x**2+2*x, x)
    assert poly_lcm(x**2+2*x, x, x) == Poly(x**2+2*x, x)
    assert poly_lcm(2*x**2+x, x, x) == Poly(x**2+x/2, x)
    assert poly_lcm(2*x**2+x, 2*x, x) == Poly(2*x**2+x, x)

    assert poly_lcm(x**2*y, x*y**2, x, y) == Poly(x**2*y**2, x, y)

    f, g = (x+y)**2*(x-5*y), (x+y)**3*(x+3*y)

    assert poly_lcm(f, g, x, y) == \
        Poly((x+y)**3*(x-5*y)*(x+3*y), x, y)

    f = 3*x*y**2 - 2*x*y**3 - 3*x*y**4 + 2*x*y**5
    g = y**2 - 2*y**4 + y**6

    assert poly_lcm(f, g, x, y) == \
        Poly(-3*x*y**2/2+x*y**3+3*x*y**4-2*x*y**5-3*x*y**6/2+x*y**7, x, y)

def test_poly_gcd():
    assert poly_gcd(2, 6, x) == Poly(2, x)
    assert poly_gcd(2, 6, x, y) == Poly(2, x, y)

    assert poly_gcd(x, y, x, y) == Poly(1, x, y)

    assert poly_gcd(2*x**3, 6*x, x) == Poly(2*x, x)
    assert poly_gcd(2*x**3, 3*x, x) == Poly(x, x)

    assert poly_gcd(2*x**3*y, 6*x*y**2, x, y) == Poly(2*x*y, x, y)
    assert poly_gcd(2*x**3*y, 3*x*y**2, x, y) == Poly(x*y, x, y)

    assert poly_gcd(x**2+2*x+1, x+1, x) == Poly(x+1, x)
    assert poly_gcd(x**2+2*x+2, x+1, x) == Poly(1, x)

    assert poly_gcd(x**2+2*x+1, 2+2*x, x) == Poly(x+1, x)
    assert poly_gcd(x**2+2*x+2, 2+2*x, x) == Poly(1, x)

    assert poly_gcd(sin(z)*(x+y), x**2+2*x*y+y**2,
        x, y) == Poly(x+y, x, y)

    f = x**8+x**6-3*x**4-3*x**3+8*x**2+2*x-5
    g = 3*x**6+5*x**4-4*x**2-9*x+21

    assert poly_gcd(f, g, x) == Poly(1, x)

def test_poly_gcdex():
    f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    g = x**3 + x**2 - 4*x - 4

    assert poly_half_gcdex(f, g, x) == \
        (Poly(-x+3, x), Poly(5*x+5, x))
    assert poly_half_gcdex(Poly(f, x), Poly(g, x)) == \
        (Poly(-x+3, x), Poly(5*x+5, x))

    assert poly_gcdex(f, g, x) == \
        (Poly(-x+3, x), Poly(x**2-6*x+10, x), Poly(5*x+5, x))
    assert poly_gcdex(Poly(f, x), Poly(g, x)) == \
        (Poly(-x+3, x), Poly(x**2-6*x+10, x), Poly(5*x+5, x))

    f = x**4 + 4*x**3 - x + 1
    g = x**3 - x + 1

    s, t, h = poly_gcdex(f, g, x)
    S, T, H = poly_gcdex(g, f, x)

    assert s*f + t*g == h
    assert S*g + T*f == H

def test_poly_resultant():
    assert poly_resultant(x**2-1, x**3-x**2+2, x) == 0
    assert poly_resultant(3*x**3-x, 5*x**2+1, x) == 64
    assert poly_resultant(x**2-2*x+7, x**3-x+5, x) == 265
    assert poly_resultant((x-a)**2-2, a**2-3, a) == 1 - 10*x**2 + x**4
    assert poly_resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-6), x) == -8640
    assert poly_resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-1), x) == 0
    assert poly_resultant(x**3-1, x**3+2*x**2+2*x-1, x) == 16
    assert poly_resultant(x**8-2, x-1, x) == -1
    assert poly_resultant(3*x**2+2*a*x+3*a**2-2,
        3*x**2-2*a*x+3*a**2-2, x) == 144*a**4 - 96*a**2
    assert poly_resultant((x-a)*(x-b), x-c, x) == a*b-a*c-b*c+c**2

def test_poly_subresultants():
    f = x**8+x**6-3*x**4-3*x**3+8*x**2+2*x-5
    g = 3*x**6+5*x**4-4*x**2-9*x+21

    assert poly_subresultants(f, g, x) == \
        (Poly(f, x), Poly(g, x),
         Poly(15*x**4 - 3*x**2 + 9,  x),
         Poly(65*x**2 + 125*x - 245, x),
         Poly(9326*x - 12300, x),
         Poly(260708, x))

    assert poly_subresultants((x-1)**2, x**2-1, x) == \
        (Poly((x-1)**2, x), Poly(x**2-1, x), Poly(2*x - 2, x))

def test_poly_groebner():
    assert poly_groebner(0, x) == (Poly((), x),)

    assert poly_groebner(x*y, x) == (Poly(x, x),)
    assert poly_groebner(x*y, z) == (Poly(1, z),)

    assert poly_groebner((x**2 + 2*x*y**2, x*y + 2*y**3 - 1), y, x, order='lex') == \
        (Poly(y**3 - Rational(1,2), y, x, order='lex'),
         Poly(x, y, x, order='lex'))

    assert poly_groebner((y-x**2, z-x**3), y, z, x, order='lex') == \
        (Poly(-x**2+y, y, z, x, order='lex'),
         Poly(z-x**3, y, z, x, order='lex'))

    assert poly_groebner((x**3-2*x*y, x**2*y-2*y**2+x), x, y, order='grlex') == \
        (Poly(x**2, x, y, order='grlex'),
         Poly(x*y, x, y, order='grlex'),
         Poly(y**2-x/2, x, y, order='grlex'))

def test_map_coeffs():
    p = Poly(x**2 + 2*x*y, x, y)
    q = p.map_coeffs(lambda c: 2*c)

    assert q.as_basic() == 2*x**2 + 4*x*y

    p = Poly(u*x**2 + v*x*y, x, y)
    q = p.map_coeffs(expand, complex=True)

    assert q.as_basic() == x**2*(I*im(u) + re(u)) + x*y*(I*im(v) + re(v))

    py.test.raises(PolynomialError, "p.map_coeffs(lambda c: x*c)")

def test_coeff():
    p = Poly(3*x**2*y + 4*x*y**2 + 1, x, y)

    assert p.coeff(0, 0) == p.coeff() == 1

    assert p.coeff(2, 1) == 3
    assert p.coeff(1, 2) == 4
    assert p.coeff(1, 1) == 0

def test_add_sub_term():
    f = Poly((), x, y, z)

    assert Poly(2, x) - Poly(1, x) == Poly(1, x)
    assert Poly(1, x) - Poly(1, x) == Poly(0, x)
    assert Poly(1, x) - Poly(2, x) == Poly(-1, x)

    assert f.add_term(12, (1,2,3)) == Poly(((12,), ((1,2,3),)), x,y,z)
    assert f.sub_term(12, (1,2,3)) == Poly(((-12,), ((1,2,3),)), x,y,z)

def test_mul_div_term():
    f = Poly(x*y**2 + 2*y, x, y)

    assert f.mul_term(0, (0, 1)) == Poly(0, x, y)

    assert f.mul_term(1) == Poly(x*y**2 + 2*y, x, y)
    assert f.mul_term(1, (0, 1)) == Poly(x*y**3 + 2*y**2, x, y)

    py.test.raises(ZeroDivisionError, "f.div_term(0)")

    assert f.div_term(1) == Poly(x*y**2 + 2*y, x, y)
    assert f.div_term(1, (0, 1)) == Poly(x*y + 2, x, y)

    assert f.mul_term(2) == Poly(2*x*y**2 + 4*y, x, y)
    assert f.mul_term(2, (0, 1)) == Poly(2*x*y**3 + 4*y**2, x, y)

    assert f.div_term(2) == Poly(x*y**2 / 2 + y, x, y)
    assert f.div_term(2, (0, 1)) == Poly(x*y / 2 + 1, x, y)

def test_kill_term():
    f = Poly(2*x**17*y + 3*x*y**15, x, y)

    assert f.kill_lead_term() == Poly(3*x*y**15, x, y)
    assert f.kill_last_term() == Poly(2*x**17*y, x, y)

def test_call():
    assert Poly(0, x)(2) == 0
    assert Poly(0, x, y, z)(1, 2, 3) == 0

    assert Poly(1, x)(2) == 1
    assert Poly(1, x, y, z)(1, 2, 3) == 1

    assert Poly(x, x)(2) == 2
    assert Poly(x, x, y, z)(1, 2, 3) == 1

    assert Poly(x, x)(t) == t
    assert Poly(x, x, y, z)(t, 2, 3) == t

    assert Poly(x*y, x)(2) == 2*y
    assert Poly(x*y, x, y, z)(1, 2, 3) == 2

    f = 2*x**17 + 3*x**15 + x**10 + 17*x**3 + 1

    assert Poly(f, x)(2) == f.subs(x, 2)
    assert Poly(f, x)(y) == 1 + y**3*(17 + y**7*(1 + y**5*(3 + 2*y**2)))

    f = x**3*y + x**2*z + x**2*y*z

    assert Poly(f, x, y, z)(1, 2, 3) == f.subs({x: 1, y: 2, z: 3})
    assert Poly(f, x, y, z)(u, v, t) == u**2*(v*(u + t) + t)

def test_content():
    assert Poly(0, x, y, z).content == 0
    assert Poly(1, x, y, z).content == 1

    assert Poly(2*x + 5*x*y, x, y).content == 1
    assert Poly(6*x + 4*x*y, x, y).content == 2
    assert Poly(2*x + z*x*y, x, y).content == 1

def test_primitive():
    assert Poly(0, x, y, z).as_primitive() == (0, Poly(0, x, y, z))
    assert Poly(1, x, y, z).as_primitive() == (1, Poly(1, x, y, z))

    assert Poly(2*x + 5*x*y, x, y).as_primitive() == (1, Poly(2*x + 5*x*y, x, y))
    assert Poly(6*x + 4*x*y, x, y).as_primitive() == (2, Poly(3*x + 2*x*y, x, y))
    assert Poly(2*x + z*x*y, x, y).as_primitive() == (1, Poly(2*x + z*x*y, x, y))

def test_reduced():
    assert Poly(x**2+1, x).as_reduced() == ((0,), Poly(x**2+1, x))
    assert Poly(x**2*y+1, x, y).as_reduced() == ((0,0), Poly(x**2*y+1, x, y))

    assert Poly(x**3+x, x).as_reduced() == ((1,), Poly(x**2+1, x))
    assert Poly(x**3*y+x**2*y**2, x, y).as_reduced() == ((2, 1), Poly(x+y, x, y))

def test_squarefree():
    assert Poly(x-1, x).is_squarefree == True
    assert Poly((x-1)**2, x).is_squarefree == False

    assert Poly(3*x**2, x).as_squarefree() == Poly(3*x, x)
    assert Poly(x**2+2*x+1, x).as_squarefree() == Poly(x+1, x)
    assert Poly(x**5-x**4-x+1, x).as_squarefree() == Poly(x**4-1, x)

    assert poly_sqf(1, x) == (Poly(1, x),)
    assert poly_sqf(x, x) == (Poly(x, x),)

    assert poly_sqf(3*x**2, x) == (Poly(3, x), Poly(x, x))
    assert poly_sqf(x**2+2*x+1, x) == (Poly(1, x), Poly(x+1, x))

    assert poly_sqf(x**5-x**4-x+1, x) == \
        (Poly(x**3 + x**2 + x + 1, x), Poly(x-1, x))
    assert poly_sqf(x**8+6*x**6+12*x**4+8*x**2, x) == \
        (Poly(1, x), Poly(x, x), Poly(x**2+2, x))

def test_decompose():
    assert poly_decompose(1, x) == (Poly(1, x),)
    assert poly_decompose(x, x) == (Poly(x, x),)

    assert poly_decompose(z*x**3, x) == (Poly(z*x**3, x),)
    assert poly_decompose(z*x**5, x) == (Poly(z*x**5, x),)

    assert poly_decompose(x**4, x) == (Poly(x**2, x), Poly(x**2, x))
    assert poly_decompose(z*x**4+1, x) == (Poly(z*x**2+1, x), Poly(x**2, x))

    assert poly_decompose(x**4+2*x**2+z, x) == (Poly(x**2+2*x+z, x), Poly(x**2, x))

    f, g = x**4 - 2*x + z, x**3 + 5*x

    assert poly_decompose(f.subs(x, g), x) == (Poly(f, x), Poly(g, x))
    assert poly_decompose(2*f.subs(x, g), x) == (Poly(2*f, x), Poly(g, x))
    assert poly_decompose(f.subs(x, g-2), x) == (Poly(f.subs(x, x-2), x), Poly(g, x))

def test_evaluate():
    f = x**2*y*z + 2*x*y*z**3 + 3*x*y + 4*y*z

    p = Poly(f, x, y, z)

    assert p.evaluate({x: 7}) == Poly(f.subs({x: 7}), y, z)
    assert p.evaluate({x: 7, y: 5}) == Poly(f.subs({x: 7, y: 5}), z)
    assert p.evaluate({x: 7, y: 5, z: 4}) == f.subs({x: 7, y: 5, z: 4})

    py.test.raises(PolynomialError, "Poly(x + y, x, y).evaluate({x: y})")

def test_subs():
    p = Poly(t*x*y**2 + x*y + t**2, x, y)

    assert p.subs(x, 2) == \
        Poly(((2*t, 2, t**2), ((2,), (1,), (0,))), y, 'grlex')
    assert p.subs(y, 2) == \
        Poly(((2 + 4*t, t**2), ((1,), (0,))), x)

    assert p.subs(x, y) == \
        Poly(((t, 1, t**2), ((3,), (2,), (0,))), y)
    assert p.subs(y, x) == \
        Poly(((t, 1, t**2), ((3,), (2,), (0,))), x)

    assert p.subs(x, z) == \
        Poly(((t, 1, t**2), ((1, 2), (1, 1), (0, 0))), z, y)
    assert p.subs(y, z) == \
        Poly(((t, 1, t**2), ((1, 2), (1, 1), (0, 0))), x, z)
    assert p.subs(t, z) == \
        Poly(((z, 1, z**2), ((1, 2), (1, 1), (0, 0))), x, y)

    assert p.subs(z, t) == \
        Poly(((t, 1, t**2), ((1, 2), (1, 1), (0, 0))), x, y)

    assert p.subs(t, x) == \
        Poly(((1, 1, 1), ((2, 2), (2, 0), (1, 1))), x, y)
    assert p.subs(t, y) == \
        Poly(((1, 1, 1), ((1, 3), (1, 1), (0, 2))), x, y)

    assert p.subs(t, sin(1)) == \
        Poly(((sin(1), 1, sin(1)**2), ((1, 2), (1, 1), (0, 0))), x, y)

    assert p.subs(t, sin(x)) == \
        x*y**2*sin(x) + x*y + sin(x)**2

    assert p.subs(x, sin(x)) == \
        t*y**2*sin(x) + y*sin(x) + t**2

def test_unify():
    p = Poly(x**2+x*y, x, y)
    q = Poly(x**2+2*x*y+1, x)
    r = Poly(x*z+y*z+1, x,y,z)

    assert p.unify_with(q) == p.unify_with(x**2+2*x*y+1) == \
        (Poly(x**2+x*y, x, y), Poly(x**2+2*x*y+1, x, y))

    assert p.unify_with([q, r]) == p.unify_with([x**2+2*x*y+1, r]) == \
        (Poly(x**2+x*y, x, y, z), [Poly(x**2+2*x*y+1, x,y,z), Poly(x*z+y*z+1, x,y,z)])

    assert p.unify_with((q, r)) == p.unify_with((x**2+2*x*y+1, r)) == \
        (Poly(x**2+x*y, x, y, z), (Poly(x**2+2*x*y+1, x,y,z), Poly(x*z+y*z+1, x,y,z)))

def test_eq_ne():
    p = Poly(x**2+x*y, x, y)
    q = Poly(x**2+x*y+1, x, y)

    assert (p == q) == False
    assert (p != q) == True

    assert (p+1 == q) == True
    assert (p+1 != q) == False

    assert (Poly(1, x) == Poly(0, x)) == False

    assert (Poly(1, x) == Poly(1, x)) == True
    assert (Poly(1, x) == Poly(1, x, y)) == True
    assert (Poly(1, x, y, z) == Poly(1, x, y)) == True

    assert (Poly(x, x, y, z) == Poly(x, x, y)) == True

    assert (Poly(x*z, x, y, z) == Poly(x, x, y)) == False
    assert (Poly(x*z, x, y, z) == Poly(x*z, x, y)) == True

    assert (Poly(x, x, y, z) == Poly(x*z, x, y)) == False
    assert (Poly(x, x, y) == Poly(x*z, x, y, z)) == False
    assert (Poly(x, x, y, z) != Poly(x*z, x, y)) == True
    assert (Poly(x, x, y) != Poly(x*z, x, y, z)) == True

    assert (Poly({ (2, 3) : 4 }, x, y) == Poly(4*x**2*y**3, x, y)) == True
    assert (Poly({ (2, 3) : 4 }, x, y) != Poly(4*x**2*y**3, x, y)) == False

    p = Poly(x*y**2 + x*y + x**2, x, y, order='lex')
    q = Poly(x*y**2 + x*y + x**2, x, y, order='grlex')

    assert (p == q) == True
    assert (p != q) == False

def test_nonzero():
    assert bool(Poly((), x, y, z)) == False
    assert bool(Poly(102, x, y, z)) == True
    assert bool(Poly(x+21*y, x, y)) == True

def test_atoms():
    assert Poly(0, x).atoms() == set([S.Zero])
    assert Poly(1, x).atoms() == set([S.One])
    assert Poly(x, x).atoms() == set([S.One, x])
    assert Poly(x, x, y).atoms() == set([S.One, x])
    assert Poly(x + y, x, y).atoms() == set([S.One, x, y])
    assert Poly(x + y, x, y, z).atoms() == set([S.One, x, y])
    assert Poly(x + y*t, x, y, z).atoms() == set([S.One, x, t, y])

    assert Poly(x+1, x).atoms(type=Symbol) == set([x])
    assert Poly(x+1, x).atoms(type=(Symbol,)) == set([x])
