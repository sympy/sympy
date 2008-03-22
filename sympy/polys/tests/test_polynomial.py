
from sympy import symbols, expand, re, im, I

from sympy.polys.monomial import *
from sympy.polys.polynomial import *
from sympy.polys.algorithms import *

import py

x,y,z,u,v,t = symbols('xyzuvt')

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
    assert f.is_constant == False
    assert f.is_monomial == True
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == True
    assert f.is_inhomogeneous == False
    assert f.is_monic == True

    f = Poly(2*x*y*z + x*y + 1, x, y, z)

    assert f.is_zero == False
    assert f.is_one == False
    assert f.is_constant == False
    assert f.is_monomial == False
    assert f.is_univariate == False
    assert f.is_multivariate == True
    assert f.is_homogeneous == False
    assert f.is_inhomogeneous == True
    assert f.is_monic == False

def test_poly_mul():
    f = x**3-12*x**2-42
    g = x**2+x-3

    assert Poly(f, x)*Poly(g, x) == f*g
    assert Poly(f, x)*g == f*g

    f = x**3*y-2*x*y**2-3*z+1
    g = x**2*y*z+x*y**3-3

    assert Poly(f, x, y, z)*Poly(g, x, y, z) == f*g
    assert Poly(f, x, y, z)*g == f*g

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

def test_map_coeffs():
    p = Poly(x**2 + 2*x*y, x, y)
    q = p.map_coeffs(lambda c: 2*c)

    assert q.as_basic() == 2*x**2 + 4*x*y

    p = Poly(u*x**2 + v*x*y, x, y)
    q = p.map_coeffs(expand, complex=True)

    assert q.as_basic() == x**2*(I*im(u) + re(u)) + x*y*(I*im(v) + re(v))

    code = \
    """
    p = Poly(x**2 + 2*x*y, x, y)
    q = p.map_coeffs(lambda c: x*c)

    """

    py.test.raises(PolynomialError, code)

def test_calculus():
    pass

def test_add_sub_term():
    pass

def test_kill_term():
    pass

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

def test_subs():
    pass

def test_unify():
    p = Poly(x**2+x*y, x, y)
    q = Poly(x**2+x*y+1, x)

    assert p.unify(q) == \
        (Poly(x**2+x*y, x, y), Poly(x**2+x*y+1, x, y))

    assert p.unify(x**2+x*y+1) == \
        (Poly(x**2+x*y, x, y), Poly(x**2+x*y+1, x, y))

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
