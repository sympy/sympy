from sympy import symbols, expand, sin, sqrt, re, im, I, Rational, Lambda, \
        separate, raises, Integer, Symbol, Poly, RootsOf, RootSum, RootOf, S

from sympy.polys.monomial import monomial_lex_cmp, monomial_grlex_cmp, \
        monomial_grevlex_cmp, monomial_1_el_cmp

from sympy.polys.algorithms import poly_groebner, poly_subresultants,    \
        poly_resultant, poly_half_gcdex, poly_gcdex, poly_gcd, poly_lcm, \
        poly_div, poly_pdiv, poly_decompose, poly_sqf, poly_reduce, \
        poly_discriminant

from sympy.polys.rootfinding import poly_root_factors, roots_linear,  \
        roots_quadratic, roots_cubic, roots_quartic, roots_binomial, \
        roots_rational, roots, number_of_real_roots, poly_sturm

from sympy.polys.factortools import poly_factors, factors, factor

from sympy.polys.polynomial import Poly, PolynomialError, \
        CoefficientError, SymbolsError

from sympy.utilities.pytest import skip

a,b,c,d,x,y,z,u,v,t = symbols('abcdxyzuvt')

def test_monomial_cmp():
    assert monomial_lex_cmp((3,2,1), (1,2,4)) == 1
    assert monomial_grlex_cmp((2,4,1), (1,6,0)) == 1
    assert monomial_grevlex_cmp((1,3,1), (1,2,2)) == 1
    assert monomial_1_el_cmp((2,0,1), (1,2,0)) == 1

def test_poly_has():
    f = x*y**2*z + I*x*y + x + 1
    assert Poly(f, x).has(x) == True
    assert Poly(f, x, y).has(x) == True
    assert Poly(f, y).has(x) == True
    assert Poly(f, x, y).has(z) == True
    assert Poly(f, x, y).has(t) == False
    assert Poly(f, x, y).has(I) == True
    assert Poly(f, y).has(I) == True
    assert Poly(f, x).has(I) == True

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

    assert Poly(Poly(x, x)) == Poly(x, x)

    assert Poly(Poly(x*y, x, y), order='lex') == \
        Poly(x*y, x, y, order='lex')

    raises(SymbolsError, "Poly(x, 2)")
    raises(SymbolsError, "Poly(x, 2*x)")

    raises(SymbolsError, "Poly(x, 2, x)")
    raises(SymbolsError, "Poly(x, 2*x, x)")

    raises(SymbolsError, "Poly(x, x, 2)")
    raises(SymbolsError, "Poly(x, x, 2*x)")

    A, B = symbols('AB', commutative=False)

    raises(SymbolsError, "Poly(x + A**2 + B**2, x, A, B)")

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
    assert Poly.cancel(x) == x
    assert Poly.cancel(x+1) == x+1
    assert Poly.cancel((x+1)/(1-x)) == (x+1)/(1-x)

    assert Poly.cancel((x**2-1)/(x-1)) == x+1
    assert Poly.cancel((x**2-y**2)/(x-y)) == x+y

    assert Poly.cancel((x**2-y)/(x-y)) == (x**2 - y)/(x - y)
    assert Poly.cancel((x**2-2)/(x+sqrt(2))) == x - sqrt(2)

    assert Poly.cancel((x**2-y**2)/(x-y), x) == x+y

    assert Poly.cancel((x, S.One), x) == x
    assert Poly.cancel((x+1, S.One), x) == x+1
    assert Poly.cancel((x+1, x-1), x) == (x+1)/(x-1)

    assert Poly.cancel((x**2-1, x-1), x) == x+1
    assert Poly.cancel((x**2-y**2, x-y), x, y) == x+y

    assert Poly.cancel((x**2-y, x-y), x, y) == (x**2 - y)/(x - y)
    assert Poly.cancel((x**2-2, x+sqrt(2)), x) == x - sqrt(2)

    assert Poly.cancel((x**2-y**2, x-y), x) == x+y

    assert Poly.cancel(((x**2-y**2).as_poly(x), (x-y).as_poly(x))) == x+y

    f = -1/(3 + 2*sqrt(2))*(1 + 1/(3 + 2*sqrt(2))*(7 + 5*sqrt(2)))

    assert Poly.cancel(f) == -2 + sqrt(2)

    a, b = x, 1/(1/y + 1/(x+y))

    assert Poly.cancel(y/(x+y) * b/(a+b), x, y) == y**2/(x**2 + 3*x*y + y**2 )

    raises(SymbolsError, "Poly.cancel((x**2-y**2, x-y))")

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

    assert p.TT == (1, (1, 1, 3))
    assert p.TM == (1, 1, 3)
    assert p.TC == 1

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

    assert p.TT == (1, (1, 1, 3))
    assert p.TM == (1, 1, 3)
    assert p.TC == 1

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

    assert p.TT == (1, (1, 1, 3))
    assert p.TM == (1, 1, 3)
    assert p.TC == 1

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

    assert Poly(x + y - 1, x).is_linear == True
    assert Poly(x + y - 1, y).is_linear == True

    assert Poly(x + y - 1, x, y).is_linear == True

    assert Poly(x*x + y - 1, x, y).is_linear == False
    assert Poly(x*y + y - 1, x, y).is_linear == False

def test_iterators():
    f = Poly(x**5 + x, x)

    assert list(f.iter_all_coeffs()) == \
        [1, 0, 0, 0, 1, 0]

    assert list(f.iter_all_monoms()) == \
        [(5,), (4,), (3,), (2,), (1,), (0,)]

    assert list(f.iter_all_terms()) == \
        [(1, (5,)), (0, (4,)), (0, (3,)), (0, (2,)), (1, (1,)), (0, (0,))]

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

def test_as_integer():
    assert Poly(3*x**2 + x, x).as_integer() == \
        (Integer(1), Poly(3*x**2 + x, x))

    assert Poly(x**2 + x/2, x).as_integer() == \
        (Integer(2), Poly(2*x**2 + x, x))

    assert Poly(25.0*x**7 + 5.0*x**2, x).as_integer() == \
        (Integer(1), Poly(25*x**7 + 5*x**2, x))

    raises(CoefficientError, "Poly(x**2 + t*x, x).as_integer()")
    raises(CoefficientError, "Poly(x**2 + 0.1, x).as_integer()")

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

def test_poly_pow():
    assert Poly(x**3-12*x**2-42, x)**2 == \
        Poly(x**6 - 24*x**5 + 144*x**4 - 84*x**3 + 1008*x**2 + 1764, x)

    assert Poly(x**3*y-2*x*y**2-3*z+1, x)**2 == \
        Poly(y**2*x**6 - 4*y**3*x**4 + (2*y - 6*y*z)*x**3 + \
        4*y**4*x**2 + (-4*y**2 + 12*z*y**2)*x + 1 - 6*z + 9*z**2, x)

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
        ([Poly(y, x, y), Poly(-1, x, y)], Poly(2, x, y))

    assert poly_div(x**2*y+x*y**2+y**2, [x*y-1, y**2-1], x, y) == \
        ([Poly(x+y, x, y), Poly(1, x, y)], Poly(1+x+y, x, y))
    assert poly_div(x**2*y+x*y**2+y**2, [y**2-1, x*y-1], x, y) == \
        ([Poly(1+x, x, y), Poly(x, x, y)], Poly(1+2*x, x, y))

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
    assert poly_gcd(0, 0, x) == Poly(0, x)
    assert poly_gcd(0, 1, x) == Poly(1, x)
    assert poly_gcd(1, 0, x) == Poly(1, x)
    assert poly_gcd(1, 1, x) == Poly(1, x)

    assert poly_gcd(x-1, 0, x) == Poly(x-1, x)
    assert poly_gcd(0, x-1, x) == Poly(x-1, x)

    assert poly_gcd(-x-1, 0, x) == Poly(x+1, x)
    assert poly_gcd(0, -x-1, x) == Poly(x+1, x)

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
        (Poly(-x/5+Rational(3,5), x), Poly(x+1, x))
    assert poly_half_gcdex(Poly(f, x), Poly(g, x)) == \
        (Poly(-x/5+Rational(3,5), x), Poly(x+1, x))

    assert poly_gcdex(f, g, x) == \
        (Poly(-x/5+Rational(3,5), x), Poly(x**2/5-6*x/5+2, x), Poly(x+1, x))
    assert poly_gcdex(Poly(f, x), Poly(g, x)) == \
        (Poly(-x/5+Rational(3,5), x), Poly(x**2/5-6*x/5+2, x), Poly(x+1, x))

    f = x**4 + 4*x**3 - x + 1
    g = x**3 - x + 1

    s, t, h = poly_gcdex(f, g, x)
    S, T, H = poly_gcdex(g, f, x)

    assert s*f + t*g == h
    assert S*g + T*f == H

    assert poly_gcdex(2*x, x**2-16, x) == \
        (Poly(x/32, x), Poly(-Rational(1,16), x), Poly(1, x))

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
        (Poly(260708, x), [Poly(f, x), Poly(g, x),
                           Poly(15*x**4 - 3*x**2 + 9,  x),
                           Poly(65*x**2 + 125*x - 245, x),
                           Poly(9326*x - 12300, x),
                           Poly(260708, x)])

    assert poly_subresultants((x-1)**2, x**2-1, x) == \
        (Poly(0, x), [Poly((x-1)**2, x),
                      Poly(x**2-1, x),
                      Poly(2*x - 2, x)])

    f = Poly(-x**3 + 5, x)
    g = Poly((1 + 3*t)*x**2, x)

    assert poly_subresultants(f, g) == \
        (Poly(25 + 225*t + 675*t**2 + 675*t**3, x), [Poly(-x**3 + 5, x),
                                                     Poly((1 + 3*t)*x**2, x),
                                                     Poly(5 + 30*t + 45*t**2, x)])

def test_poly_groebner():
    assert poly_groebner(0, x) == [Poly((), x)]

    assert poly_groebner(x*y, x) == [Poly(x, x)]
    assert poly_groebner(x*y, z) == [Poly(1, z)]

    assert poly_groebner((x**2 + 2*x*y**2, x*y + 2*y**3 - 1), y, x, order='lex') == \
        [Poly(y**3 - Rational(1,2), y, x, order='lex'),
         Poly(x, y, x, order='lex')]

    assert poly_groebner((y-x**2, z-x**3), y, z, x, order='lex') == \
        [Poly(-x**2+y, y, z, x, order='lex'),
         Poly(z-x**3, y, z, x, order='lex')]

    assert poly_groebner((x**3-2*x*y, x**2*y-2*y**2+x), x, y, order='grlex') == \
        [Poly(x**2, x, y, order='grlex'),
         Poly(x*y, x, y, order='grlex'),
         Poly(y**2-x/2, x, y, order='grlex')]

def test_map_coeffs():
    p = Poly(x**2 + 2*x*y, x, y)
    q = p.map_coeffs(lambda c: 2*c)

    assert q.as_basic() == 2*x**2 + 4*x*y

    p = Poly(u*x**2 + v*x*y, x, y)
    q = p.map_coeffs(expand, complex=True)

    assert q.as_basic() == x**2*(I*im(u) + re(u)) + x*y*(I*im(v) + re(v))

    raises(PolynomialError, "p.map_coeffs(lambda c: x*c)")

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

    raises(ZeroDivisionError, "f.div_term(0)")

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
    assert Poly(0, x, y, z).content == 1
    assert Poly(1, x, y, z).content == 1

    assert Poly(2*x + 5*x*y, x, y).content == 1
    assert Poly(6*x + 4*x*y, x, y).content == 2
    assert Poly(2*x + z*x*y, x, y).content == 1

def test_primitive():
    assert Poly(0, x, y, z).as_primitive() == (1, Poly(0, x, y, z))
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

    assert poly_sqf(1, x) == [Poly(1, x)]
    assert poly_sqf(x, x) == [Poly(x, x)]

    assert poly_sqf(3*x**2, x) == [Poly(3, x), Poly(x, x)]
    assert poly_sqf(x**2+2*x+1, x) == [Poly(1, x), Poly(x+1, x)]

    assert poly_sqf(x**5-x**4-x+1, x) == \
        [Poly(x**3 + x**2 + x + 1, x), Poly(x-1, x)]
    assert poly_sqf(x**8+6*x**6+12*x**4+8*x**2, x) == \
        [Poly(1, x), Poly(x, x), Poly(x**2+2, x)]

    # Bronstein, Symbolic Integration, pp. 52

    A = Poly(x**4 - 3*x**2 + 6, x)
    D = Poly(x**6 - 5*x**4 + 5*x**2 + 4, x)

    f, g = D, A - D.diff(x).mul_term(t)

    res, R = poly_subresultants(f, g)
    S = poly_sqf(Poly(res, t))

    assert S == [Poly(45796, t), Poly(1, t), Poly(4*t**2 + 1, t)]

def test_decompose():
    assert poly_decompose(1, x) == [Poly(1, x)]
    assert poly_decompose(x, x) == [Poly(x, x)]

    assert poly_decompose(z*x**3, x) == [Poly(z*x**3, x)]
    assert poly_decompose(z*x**5, x) == [Poly(z*x**5, x)]

    assert poly_decompose(x**4, x) == [Poly(x**2, x), Poly(x**2, x)]
    assert poly_decompose(z*x**4+1, x) == [Poly(z*x**2+1, x), Poly(x**2, x)]

    assert poly_decompose(x**4+2*x**2+z, x) == [Poly(x**2+2*x+z, x), Poly(x**2, x)]

    f, g = x**4 - 2*x + z, x**3 + 5*x

    assert poly_decompose(f.subs(x, g), x) == [Poly(f, x), Poly(g, x)]
    assert poly_decompose(2*f.subs(x, g), x) == [Poly(2*f, x), Poly(g, x)]
    assert poly_decompose(f.subs(x, g-2), x) == [Poly(f.subs(x, x-2), x), Poly(g, x)]

def test_reduce():
    f = Poly(2930944*x**6 + 2198208*x**4 + 549552*x**2 + 45796, x)
    g = Poly(17585664*x**5 + 8792832*x**3 + 1099104*x, x)

    F, G = poly_reduce(f, g)

    assert F == Poly(64*x**6 + 48*x**4 + 12*x**2 + 1, x)
    assert G == Poly(384*x**5 + 192*x**3 + 24*x, x)

def test_evaluate():
    f = x**2*y*z + 2*x*y*z**3 + 3*x*y + 4*y*z

    p = Poly(f, x, y, z)

    assert p.evaluate({x: 7}) == Poly(f.subs({x: 7}), y, z)
    assert p.evaluate({x: 7, y: 5}) == Poly(f.subs({x: 7, y: 5}), z)
    assert p.evaluate({x: 7, y: 5, z: 4}) == f.subs({x: 7, y: 5, z: 4})

    p = Poly(x**2 + x*y*sin(z), x, y, order='lex')

    assert p.evaluate({x: 3, y: 0}) == 9
    assert p.evaluate({x: 0, y: 0}) == 0

    q = p.evaluate({x: 0})

    assert q == Poly(0, y, order='lex')
    assert q.is_zero == True

    assert p.evaluate({y: 0}) == \
        Poly(x**2, x, order='lex')

    raises(PolynomialError, "Poly(x + y, x, y).evaluate({x: y})")

def test_subs():
    p = Poly(t*x*y**2 + x*y + t**2, x, y)

    assert p.subs(x, 2) == \
        Poly(((2*t, 2, t**2), ((2,), (1,), (0,))), y)
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

    coeffs = [ Symbol('a' + str(i)) for i in [3,2,1,0] ]
    values = [ 4, 0, 0, 1 ]

    f = Poly(zip(coeffs, [3,2,1,0]), x)

    assert f.subs(dict(zip(coeffs, values))) == Poly(4*x**3+1, x)

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

def test_diff():
    f = Poly(a*x**2 + b*x + 2, x)

    assert f.diff(x) == Poly(2*a*x + b, x)
    assert f.diff(y) == Poly(0, x)
    assert f.diff(a) == Poly(x**2, x)
    assert f.diff(b) == Poly(x, x)

    assert f.diff() == Poly(2*a*x + b, x)

    g = Poly(a*x**2 + b*x*y + 2, x, y)

    assert g.diff(x) == Poly(2*a*x + b*y, x, y)
    assert g.diff(y) == Poly(b*x, x, y)
    assert g.diff(a) == Poly(x**2, x, y)
    assert g.diff(b) == Poly(x*y, x, y)

    assert g.diff() == g

def test_invert():
    assert Poly(2*x, x).invert(x**2-16) == Poly(x/32, x)

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

def test_sturm():
    assert poly_sturm(5, x) == [Poly(5, x)]
    assert poly_sturm(x, x) == [Poly(x, x), Poly(1, x)]

    assert poly_sturm(x**3-2*x**2+3*x-5, x) == \
        [Poly(x**3-2*x**2+3*x-5, x), Poly(3*x**2-4*x+3, x),
         Poly(-Rational(10,9)*x+Rational(13,3), x),
         Poly(Rational(-3303,100), x)]

def test_number_of_real_roots():
    assert number_of_real_roots(0, x) == 0
    assert number_of_real_roots(7, x) == 0

    f = Poly(x - 1, x)

    assert number_of_real_roots(f) == 1
    assert number_of_real_roots(f, sup=0) == 0
    assert number_of_real_roots(f, inf=1) == 0
    assert number_of_real_roots(f, sup=0, inf=1) == 1

    assert number_of_real_roots(f, sup=1, inf=0) == 1

    f = x**2 - 4

    assert number_of_real_roots(f, x) == 2
    assert number_of_real_roots(f, x, sup=0) == 1
    assert number_of_real_roots(f, x, inf=-1, sup=1) == 0

    raises(ValueError, "number_of_real_roots(f, x, inf=t)")
    raises(ValueError, "number_of_real_roots(f, x, sup=t)")

def test_roots_linear():
    assert roots_linear(Poly(2*x+1, x)) == [-Rational(1, 2)]

def test_roots_quadratic():
    assert roots_quadratic(Poly(2*x**2, x)) == [0, 0]
    assert roots_quadratic(Poly(2*x**2+3*x, x)) == [0, -Rational(3, 2)]
    assert roots_quadratic(Poly(2*x**2+3, x)) == [I*sqrt(6)/2, -I*sqrt(6)/2]

    assert roots_quadratic(Poly(2*x**2+4*x+3, x)) == \
        [-1 + I*sqrt(2)/2, -1 - I*sqrt(2)/2]

def test_roots_cubic():
    assert roots_cubic(Poly(2*x**3, x)) == [0, 0, 0]
    assert roots_cubic(Poly(x**3-3*x**2+3*x-1, x)) == [1, 1, 1]

    assert roots_cubic(Poly(x**3+1, x)) == \
        [-1, S.Half - I*sqrt(3)/2, S.Half + I*sqrt(3)/2]

def test_roots_quartic():
    assert roots_quartic(Poly(x**4, x)) == [0, 0, 0, 0]
    assert roots_quartic(Poly(x**4 + x**3, x)) in [
        [-1,0,0,0],
        [0,-1,0,0],
        [0,0,-1,0],
        [0,0,0,-1]
    ]
    assert roots_quartic(Poly(x**4 - x**3, x)) in [
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
    ]
    assert roots_quartic(Poly(x**4 + x, x)) == [S.Half + I*sqrt(3)/2, S.Half - I*sqrt(3)/2, 0, -1]

def test_roots_binomial():
    assert roots_binomial(Poly(5*x, x)) == [0]
    assert roots_binomial(Poly(5*x**4, x)) == [0, 0, 0, 0]
    assert roots_binomial(Poly(5*x+2, x)) == [-Rational(2, 5)]

    A = 10**Rational(3, 4)/10

    assert roots_binomial(Poly(5*x**4+2, x)) == \
        [A+I*A, -A+I*A, -A-I*A, A-I*A]

    a1 = Symbol('a1', nonnegative=True)
    b1 = Symbol('b1', nonnegative=True)

    r0 = roots_quadratic(Poly(a1*x**2 + b1, x))
    r1 = roots_binomial(Poly(a1*x**2 + b1, x))

    assert separate(r0[0]) == separate(r1[0])
    assert separate(r0[1]) == separate(r1[1])

def test_roots_rational():
    assert roots_rational(Poly(x**2-1, x)) == [S.One, -S.One]
    assert roots_rational(Poly(x**2-x, x)) == [S.Zero, S.One]

    assert roots_rational(Poly(x**2-x/2, x)) == [S.Zero]
    assert roots_rational(Poly(2*x**2-x, x)) == [S.Zero]

    assert roots_rational(Poly(t*x**2-x, x)) == []

def test_roots():
    assert roots(1, x) == {}
    assert roots(x, x) == {S.Zero: 1}
    assert roots(x**9, x) == {S.Zero: 9}
    assert roots(((x-2)*(x+3)*(x-4)).expand(), x) == {-S(3): 1, S(2): 1, S(4): 1}

    assert roots(2*x+1, x) == {-S.Half: 1}
    assert roots((2*x+1)**2, x) == {-S.Half: 2}
    assert roots((2*x+1)**5, x) == {-S.Half: 5}
    assert roots((2*x+1)**10, x) == {-S.Half: 10}

    assert roots(x**4 - 1, x) == {I: 1, S.One: 1, -S.One: 1, -I: 1}
    assert roots((x**4 - 1)**2, x) == {I: 2, S.One: 2, -S.One: 2, -I: 2}

    assert roots(((2*x-3)**2).expand(), x) == { Rational(3,2): 2}
    assert roots(((2*x+3)**2).expand(), x) == {-Rational(3,2): 2}

    assert roots(((2*x-3)**3).expand(), x) == { Rational(3,2): 3}
    assert roots(((2*x+3)**3).expand(), x) == {-Rational(3,2): 3}

    assert roots(((2*x-3)**5).expand(), x) == { Rational(3,2): 5}
    assert roots(((2*x+3)**5).expand(), x) == {-Rational(3,2): 5}

    assert roots(((a*x-b)**5).expand(), x) == { b/a: 5}
    assert roots(((a*x+b)**5).expand(), x) == {-b/a: 5}

    assert roots(x**4-2*x**2+1, x) == {S.One: 2, -S.One: 2}

    assert roots(x**6-4*x**4+4*x**3-x**2, x) == \
        {S.One: 2, -1 - sqrt(2): 1, S.Zero: 2, -1 + sqrt(2): 1}

    assert roots(x**8-1, x) == {
         2**S.Half/2 + I*2**S.Half/2: 1,
         2**S.Half/2 - I*2**S.Half/2: 1,
        -2**S.Half/2 + I*2**S.Half/2: 1,
        -2**S.Half/2 - I*2**S.Half/2: 1,
        S.One: 1, -S.One: 1, I: 1, -I: 1
    }

    assert roots(-2016*x**2 - 5616*x**3 - 2056*x**4 + 3324*x**5 + 2176*x**6 \
        - 224*x**7 - 384*x**8 - 64*x**9, x) == {S(0): 2, -S(2): 2, S(2): 1, -S(7)/2: 1,\
                                            -S(3)/2: 1, -S(1)/2: 1, S(3)/2: 1}

    assert roots((a+b+c)*x + a+b+c+d, x) == \
        { (-a-b-c-d) / (a+b+c) : 1 }

    assert roots(x**3+x**2-x+1, x, cubics=False) == {}
    assert roots(((x-2)*(x+3)*(x-4)).expand(), x, cubics=False) == {-S(3): 1, S(2): 1, S(4): 1}
    assert roots(((x-2)*(x+3)*(x-4)*(x-5)).expand(), x, cubics=False) == \
            {-S(3): 1, S(2): 1, S(4): 1, S(5): 1}
    assert roots(x**3 + 2*x**2 + 4*x + 8, x) == {-S(2): 1, -2*I: 1, 2*I: 1}
    assert roots(x**3 + 2*x**2 + 4*x + 8, x, cubics=True) == \
                {-2*I: 1, 2*I: 1, -S(2): 1}
    assert roots((x**2 - x)*(x**3 + 2*x**2 + 4*x + 8), x ) == \
                {S(1): 1, S(0): 1, -S(2): 1, -2*I: 1, 2*I: 1}

    r1_2, r1_3, r1_9, r4_9, r19_27 = [ Rational(*r) \
        for r in ((1,2), (1,3), (1,9), (4,9), (19,27)) ]

    assert roots(x**3+x**2-x+1, x, cubics=True) in [
            {
        -r1_3 - (r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 - \
        r4_9*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3): 1,

        -r1_3 + r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 + \
        r4_9/(r1_2 + r1_2*I*3**r1_2)*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3) + \
        r1_2*I*3**r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3: 1,

        -r1_3 + r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 + \
        r4_9/(r1_2 - r1_2*I*3**r1_2)*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3) - \
        r1_2*I*3**r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3: 1,
            },
            {
        -r1_3 - (r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 - \
        r4_9*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3): 1,

        -r1_3 + r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 - \
        r4_9/(-r1_2 - r1_2*I*3**r1_2)*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3) + \
        r1_2*I*3**r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3: 1,

        -r1_3 + r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3 + \
        r4_9/(r1_2 - r1_2*I*3**r1_2)*(r19_27 + r1_9*3**r1_2*11**r1_2)**(-r1_3) - \
        r1_2*I*3**r1_2*(r19_27 + r1_9*3**r1_2*11**r1_2)**r1_3: 1,
            },
            ]

    f = (x**2+2*x+3).subs(x, 2*x**2 + 3*x).subs(x, 5*x-4)

    r1_2, r13_20, r1_100 = [ Rational(*r) \
        for r in ((1,2), (13,20), (1,100)) ]

    assert roots(f, x) == {
        r13_20 + r1_100*(25 - 200*I*2**r1_2)**r1_2: 1,
        r13_20 - r1_100*(25 - 200*I*2**r1_2)**r1_2: 1,
        r13_20 + r1_100*(25 + 200*I*2**r1_2)**r1_2: 1,
        r13_20 - r1_100*(25 + 200*I*2**r1_2)**r1_2: 1,
    }

    p = Poly(z**3 + (-2 - y)*z**2 + (1 + 2*y - 2*x**2)*z - y + 2*x**2, z)

    assert roots(p) == {
        S.One: 1,
        S.Half + S.Half*y + S.Half*(1 - 2*y + y**2 + 8*x**2)**S.Half: 1,
        S.Half + S.Half*y - S.Half*(1 - 2*y + y**2 + 8*x**2)**S.Half: 1,
    }

    assert roots(a*b*c*x**3 + 2*x**2 + 4*x + 8, x, cubics=False) == {}
    assert roots(a*b*c*x**3 + 2*x**2 + 4*x + 8, x, cubics=True) != {}

    assert roots(x**4-1, x, domain='Z') == {S.One: 1, -S.One: 1}
    assert roots(x**4-1, x, domain='I') == {I: 1, -I: 1}

    assert roots((x-1)*(x+1), x) == {S.One: 1, -S.One: 1}
    assert roots((x-1)*(x+1), x, predicate=lambda r: r.is_positive) == {S.One: 1}

    assert roots(x**4-1, x, domain='Z', multiple=True) == [S.One, -S.One]
    assert roots(x**4-1, x, domain='I', multiple=True) in ([I, -I], [-I, I])

    assert roots(x**3, x, multiple=True) == [S.Zero, S.Zero, S.Zero]
    assert roots(1234, x, multiple=True) == []

def test_roots2():
    """Just test that calculating these roots does not hang
    (final result is not checked)
    """
    a, b, c, d, x = symbols("a b c d x")
    f1 = x**2*c + (a/b) + x*c*d - a
    f2 = x**2*(a + b*(c-d)*a) + x*a*b*c/(b*d-d) + (a*d-c/d)
    assert roots(f1, x).values() == [1, 1]
    assert roots(f2, x).values() == [1, 1]

def test_root_factors():
    assert poly_root_factors(1, x) == [Poly(1, x)]
    assert poly_root_factors(x, x) == [Poly(x, x)]

    assert poly_root_factors(x**2-1, x) == [Poly(x-1, x), Poly(x+1, x)]

    assert set( poly_root_factors((x**4 - 1)**2, x) ) == \
    set([Poly(((S.One, -I),     ((1,), (0,))), x),
         Poly(((S.One, -I),     ((1,), (0,))), x),
         Poly(((S.One, -S.One), ((1,), (0,))), x),
         Poly(((S.One, -S.One), ((1,), (0,))), x),
         Poly(((S.One, S.One),  ((1,), (0,))), x),
         Poly(((S.One, S.One),  ((1,), (0,))), x),
         Poly(((S.One, I),      ((1,), (0,))), x),
         Poly(((S.One, I),      ((1,), (0,))), x)])

    assert poly_root_factors(x**4-1, x, domain='Z') == \
        [Poly(x-1, x), Poly(x+1, x), Poly(x**2+1, x)]

def test_RootsOf():
    f = Poly((x-4)**4, x)

    roots = RootsOf(f)

    assert roots.count == 4

    assert list(roots.roots()) == [ Integer(4),
        Integer(4), Integer(4), Integer(4) ]

    assert RootSum(lambda r: r**2, f) == 64

    roots = RootsOf(x**5+x+1, x)

    assert roots.count == 5

    f = Poly(x**5+x+1, x)

    assert list(roots.roots()) == [ RootOf(f, 0), RootOf(f, 1),
        RootOf(f, 2), RootOf(f, 3), RootOf(f, 4) ]

    assert RootSum(lambda r: r**2, f).doit() == RootOf(f, 0)**2 + \
        RootOf(f, 1)**2 + RootOf(f, 2)**2 + RootOf(f, 3)**2 + RootOf(f, 4)**2

    assert RootSum(Lambda(x, x), Poly(0, x), evaluate=True)  == S.Zero
    assert RootSum(Lambda(x, x), Poly(0, x), evaluate=False) != S.Zero

    assert RootSum(Lambda(x, x), Poly(x-1, x), evaluate=False).doit() == S.One

def test_factor():
    assert factor(-2) == -2
    assert factor(-x) == -x

    assert factor(Rational(3, 8)*x**2 - Rational(3, 2)) == \
        Rational(3, 8)*((x - 2)*(x + 2))

    assert factor(x**3 - 1) == (x - 1)*(x**2 + x + 1)
    assert factor(x**3 - x) == x*(x - 1)*(x + 1)

    assert factor(x**2 + 2*x + 1) == (x + 1)**2
    assert factor(x**2 +   x - 2) == (x - 1)*(x + 2)

    assert factor(2*x**2 + 5*x + 2) == (x + 2)*(2*x + 1)

    assert factor(x**2 + y**2) == x**2 + y**2
    assert factor(x**5 - y**2) == x**5 - y**2

    assert factor(x*y + x*z + y*z) == x*y + x*z + y*z
    assert factor(x*y + x*z + x**2) == x*(x + y + z)

    assert poly_factors((a*x - b)**5, x) == \
        (1, [(Poly(a*x - b, x), 5)])

    assert poly_factors((a*x - b)**5, x, a) == \
        (1, [(Poly(x*a - b, x, a), 5)])

    assert poly_factors((a*x - b)**5, x, a, b) == \
        (1, [(Poly(x*a - b, x, a, b), 5)])

    assert poly_factors(-2*x**2 + x, x) == \
        (-1, [(Poly(x, x), 1),
              (Poly(2*x - 1, x), 1)])

    assert poly_factors(x**3 - 3*x**2 + 3*x - 1, x) == \
        (1, [(Poly(x - 1, x), 3)])

    assert poly_factors(x**6 - 1, x) == \
        (1, [(Poly(x - 1, x), 1),
             (Poly(x + 1, x), 1),
             (Poly(x**2 - x + 1, x), 1),
             (Poly(x**2 + x + 1, x), 1)])

    assert poly_factors(2*x**3*y - 2*a**2*x*y - 3*a**2*x**2 + 3*a**4, a, x, y) == \
        (1, [(Poly(a - x, a, x, y), 1),
             (Poly(a + x, a, x, y), 1),
             (Poly(3*a**2 - 2*x*y, a, x, y), 1)])

    assert poly_factors(x**20 - z**5*y**20, x, y, z) == \
        (1, [(Poly(-y**4*z + x**4, x, y, z), 1),
             (Poly(y**16*z**4 + x**4*y**12*z**3 + x**8*y**8*z**2 + x**12*y**4*z + x**16, x, y, z) , 1)])

def test_discriminant():
    e = Symbol('e')
    assert poly_discriminant(Poly(a, x), x) == S(0)
    assert poly_discriminant(Poly(a*x + b, x), x) == S(1)
    assert poly_discriminant(Poly(a*x**2 + b*x + c, x), x) == -4*a*c + b**2
    assert poly_discriminant(Poly(a*x**2 + b*x + c, x), y) == S.Zero
    assert poly_discriminant(Poly(a*x**3 + b*x**2 + c*x + d, x), x) == \
    18*a*b*c*d + b**2*c**2 - 27*a**2*d**2 - 4*a*c**3 - 4*d*b**3
    assert poly_discriminant(Poly(a*x**4 + b*x**3 + c*x**2 + d*x + e, x), x) == \
    -27*b**4*e**2 + b**2*c**2*d**2 - 128*a**2*c**2*e**2 - 4*a*c**3*d**2 - \
    192*b*d*a**2*e**2 - 6*a*e*b**2*d**2 + 144*a*c*b**2*e**2 + 144*c*e*a**2*d**2 \
    - 80*a*b*d*e*c**2 - 4*b**3*d**3 + 256*a**3*e**3 - 4*e*b**2*c**3 + \
    18*a*b*c*d**3 + 18*c*d*e*b**3 - 27*a**2*d**4 + 16*a*e*c**4
    assert poly_discriminant(Poly(5*x**5 + x**3 + 2, x), x) == S(31252160)
    assert poly_discriminant(5*x**5 + x**3 + 2, x) == S(31252160)
    assert poly_discriminant(Poly(12*x**7 + 15*x**4 + 30*x**3 + x**2 + 1, x), \
    x) == S(-220289699947514112)
    # (x - 1)**2*(x + 2 - 3*I)*(x + 2 + 3*I)
    assert poly_discriminant(Poly(13 - 22*x + 6*x**2 + 2*x**3 + x**4, x), x) == S(0)
    # (x - 1)*(x + 2 - 3*I)*(x + 2 + 3*I)
    assert poly_discriminant(Poly(-13 + 9*x + 3*x**2 + x**3, x), x) == S(-11664)
    assert poly_discriminant(Poly(x**2*y + 2*y, x, y), x) == -8*y**2
    assert poly_discriminant(Poly(x**2*y + 2*y, x, y), y) == 1

def test_discriminant_5th():
    skip('takes too much time')
    # but the result is correct nonetheless
    e, f = symbols('ef')
    assert poly_discriminant(Poly(a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x + f, x), x) \
    == -128*b**4*d**2*f**2 - 27*a**2*d**4*e**2 - 27*b**2*c**4*f**2 - \
    4*b**3*d**3*e**2 + 108*a*c**5*f**2 + b**2*c**2*d**2*e**2 - \
    900*b*a**2*d**3*f**2 - 900*e*a**2*c**3*f**2 - 192*c*e*b**4*f**2 - \
    50*a**2*b**2*e**2*f**2 - 6*f*b**3*c**2*e**2 - 4*a*c**3*d**2*e**2 + \
    144*d*f*b**4*e**2 + 144*d*b**3*c**2*f**2 + 825*a**2*c**2*d**2*f**2 + \
    2000*c*a**3*e**2*f**2 + 2250*e*a**3*d**2*f**2 - 630*a*b*d*c**3*f**2 - \
    80*c*e*f*b**3*d**2 + 18*a*b*c*d**3*e**2 + 24*a*b*f*c**3*e**2 + \
    160*a*d*e*b**3*f**2 + 560*a*c*b**2*d**2*f**2 + 560*d*f*a**2*c**2*e**2 + \
    1020*a*e*b**2*c**2*f**2 + 1020*b*f*a**2*d**2*e**2 - 2050*b*c*d*e*a**2*f**2\
    - 746*a*c*d*f*b**2*e**2 + 356*a*b*e*f*c**2*d**2 + 256*b**5*f**3 - \
    4*b**2*c**3*e**3 + 16*a*c**4*e**3 - 3750*c*d*a**3*f**3 - 2500*b*e*a**3*f**3\
    - 1600*a*c*b**3*f**3 - 1600*d*f*a**3*e**3 - 36*a*f*b**3*e**3 - \
    6*a*b**2*d**2*e**3 - 4*f*b**2*c**2*d**3 + 16*a*f*c**3*d**3 + \
    18*c*d*b**3*e**3 + 144*c*a**2*d**2*e**3 + 2000*d*a**2*b**2*f**3 + \
    2250*b*a**2*c**2*f**3 - 630*c*e*f*a**2*d**3 - 80*a*b*d*c**2*e**3 + \
    18*d*e*f*b**2*c**3 + 24*a*e*f*b**2*d**3 + 160*b*c*f*a**2*e**3 - \
    27*b**4*e**4 + 3125*a**4*f**4 - 128*a**2*c**2*e**4 + 16*f*b**3*d**4 - \
    192*b*d*a**2*e**4 + 144*a*c*b**2*e**4 - 72*a*b*c*f*d**4 - 72*a*d*e*f*c**4 + \
    256*a**3*e**5 + 108*f*a**2*d**5

