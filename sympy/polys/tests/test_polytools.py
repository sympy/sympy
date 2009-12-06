"""Tests for user-friendly public interface to polynomial functions. """

from sympy.polys.polytools import (
    Poly, _polify_basic,
    pdiv, prem, pquo, pexquo,
    div, rem, quo, exquo,
    half_gcdex, gcdex, invert,
    subresultants,
    resultant, discriminant,
    cofactors, gcd, lcm, terms_gcd,
    monic, content, primitive,
    compose, decompose,
    sqf_part, sqf_list, sqf,
    factor_list, factor,
    cancel, sturm,
    groebner,
)

from sympy.polys.polytools import (
    GeneratorsNeeded,
    CoercionFailed,
)

from sympy.polys.algebratools import ZZ, QQ, EX

from sympy import S, Integer, Rational, symbols, sqrt, exp, sin, expand, raises

x,y,z,p,q,r,s,t,u,v,w,a,b,c = symbols('x,y,z,p,q,r,s,t,u,v,w,a,b,c')

def test_Poly__gens():
    assert Poly((x-p)*(x-q), x).gens == (x,)
    assert Poly((x-p)*(x-q), p).gens == (p,)
    assert Poly((x-p)*(x-q), q).gens == (q,)

    assert Poly((x-p)*(x-q), x, p).gens == (x, p)
    assert Poly((x-p)*(x-q), x, q).gens == (x, q)

    assert Poly((x-p)*(x-q), x, p, q).gens == (x, p, q)
    assert Poly((x-p)*(x-q), p, x, q).gens == (p, x, q)
    assert Poly((x-p)*(x-q), p, q, x).gens == (p, q, x)

    assert Poly((x-p)*(x-q)).gens == (x, p, q)

    assert Poly((x-p)*(x-q), sort='x < p < q').gens == (x, p, q)
    assert Poly((x-p)*(x-q), sort='p < x < q').gens == (p, x, q)
    assert Poly((x-p)*(x-q), sort='p < q < x').gens == (p, q, x)

    assert Poly((x-p)*(x-q), x, p, q, sort='p < q < x').gens == (x, p, q)

    assert Poly((x-p)*(x-q), wrt='x').gens == (x, p, q)
    assert Poly((x-p)*(x-q), wrt='p').gens == (p, x, q)
    assert Poly((x-p)*(x-q), wrt='q').gens == (q, x, p)

    assert Poly((x-p)*(x-q), wrt=x).gens == (x, p, q)
    assert Poly((x-p)*(x-q), wrt=p).gens == (p, x, q)
    assert Poly((x-p)*(x-q), wrt=q).gens == (q, x, p)

    assert Poly((x-p)*(x-q), x, p, q, wrt='p').gens == (x, p, q)

    assert Poly((x-p)*(x-q), wrt='p', sort='q < x').gens == (p, q, x)
    assert Poly((x-p)*(x-q), wrt='q', sort='p < x').gens == (q, p, x)

def test_Poly_get_domain():
    assert Poly(2*x).get_domain() == ZZ

    assert Poly(2*x, domain='ZZ').get_domain() == ZZ
    assert Poly(2*x, domain='QQ').get_domain() == QQ

    assert Poly(x/2).get_domain() == QQ

    raises(CoercionFailed, "Poly(x/2, domain='ZZ')")
    assert Poly(x/2, domain='QQ').get_domain() == QQ

def test__init_poly_from_poly():
    assert Poly(x*y, x, y).gens == (x, y)
    assert Poly(x*y, y, x).gens == (y, x)

    assert Poly(Poly(x*y, x, y), y, x).gens == (y, x)

def test__init_poly_from_basic():
    assert Poly(0).is_Poly == False
    assert Poly(1).is_Poly == False

    assert Poly(0, x).is_Poly == True
    assert Poly(1, x).is_Poly == True

def test_Poly_abs():
    assert Poly(-x+1, x).abs() == abs(Poly(-x+1, x)) == Poly(x+1, x)

def test_Poly_neg():
    assert Poly(-x+1, x).neg() == -Poly(-x+1, x) == Poly(x-1, x)

def test_Poly_add():
    assert Poly(0, x).add(Poly(0, x)) == Poly(0, x)
    assert Poly(0, x) + Poly(0, x) == Poly(0, x)

    assert Poly(1, x).add(Poly(0, x)) == Poly(1, x)
    assert Poly(1, x, y) + Poly(0, x) == Poly(1, x, y)
    assert Poly(0, x).add(Poly(1, x, y)) == Poly(1, x, y)
    assert Poly(0, x, y) + Poly(1, x, y) == Poly(1, x, y)

    assert Poly(1, x) + x == Poly(x+1, x)
    assert Poly(1, x) + sin(x) == 1+sin(x)

def test_Poly_sub():
    assert Poly(0, x).sub(Poly(0, x)) == Poly(0, x)
    assert Poly(0, x) - Poly(0, x) == Poly(0, x)

    assert Poly(1, x).sub(Poly(0, x)) == Poly(1, x)
    assert Poly(1, x, y) - Poly(0, x) == Poly(1, x, y)
    assert Poly(0, x).sub(Poly(1, x, y)) == Poly(-1, x, y)
    assert Poly(0, x, y) - Poly(1, x, y) == Poly(-1, x, y)

    assert Poly(1, x) - x == Poly(1-x, x)
    assert Poly(1, x) - sin(x) == 1-sin(x)

def test_Poly_mul():
    assert Poly(0, x).mul(Poly(0, x)) == Poly(0, x)
    assert Poly(0, x) * Poly(0, x) == Poly(0, x)

    assert Poly(2, x).mul(Poly(4, x)) == Poly(8, x)
    assert Poly(2, x, y) * Poly(4, x) == Poly(8, x, y)
    assert Poly(4, x).mul(Poly(2, x, y)) == Poly(8, x, y)
    assert Poly(4, x, y) * Poly(2, x, y) == Poly(8, x, y)

    assert Poly(1, x) * x == Poly(x, x)
    assert Poly(1, x) * sin(x) == sin(x)

def test_Poly_sqr():
    assert Poly(x*y, x, y).sqr() == Poly(x**2*y**2, x, y)

def test_Poly_pow():
    assert Poly(x, x).pow(10) == Poly(x**10, x)
    assert Poly(x, x).pow(Integer(10)) == Poly(x**10, x)

    assert Poly(2*y, x, y).pow(4) == Poly(16*y**4, x, y)
    assert Poly(2*y, x, y).pow(Integer(4)) == Poly(16*y**4, x, y)

    assert Poly(7*x*y, x, y)**3 == Poly(343*x**3*y**3, x, y)

    assert Poly(x*y+1, x, y)**(-1) == (x*y+1)**(-1)
    assert Poly(x*y+1, x, y)**x == (x*y+1)**x

def test_Poly_divmod():
    f, g = Poly(x**2), Poly(x)
    q, r = g, Poly(0, x)

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    assert divmod(f, x) == (q, r)
    assert f // x == q
    assert f % x == r

def test_Poly_eq_ne():
    assert (Poly(x+y, x, y) == Poly(x+y, x, y)) == True
    assert (Poly(x+y, x) == Poly(x+y, x, y)) == False
    assert (Poly(x+y, x, y) == Poly(x+y, x)) == False
    assert (Poly(x+y, x) == Poly(x+y, x)) == True
    assert (Poly(x+y, y) == Poly(x+y, y)) == True

    assert (Poly(x+y, x, y) == x+y) == True
    assert (Poly(x+y, x) == x+y) == True
    assert (Poly(x+y, x, y) == x+y) == True
    assert (Poly(x+y, x) == x+y) == True
    assert (Poly(x+y, y) == x+y) == True

    assert (Poly(x+y, x, y) != Poly(x+y, x, y)) == False
    assert (Poly(x+y, x) != Poly(x+y, x, y)) == True
    assert (Poly(x+y, x, y) != Poly(x+y, x)) == True
    assert (Poly(x+y, x) != Poly(x+y, x)) == False
    assert (Poly(x+y, y) != Poly(x+y, y)) == False

    assert (Poly(x+y, x, y) != x+y) == False
    assert (Poly(x+y, x) != x+y) == False
    assert (Poly(x+y, x, y) != x+y) == False
    assert (Poly(x+y, x) != x+y) == False
    assert (Poly(x+y, y) != x+y) == False

    assert (Poly(x, x) == sin(x)) == False
    assert (Poly(x, x) != sin(x)) == True

def test_Poly_nonzero():
    assert not bool(Poly(0, x)) == True
    assert not bool(Poly(1, x)) == False

def test_Poly_properties():
    assert Poly(0, x).is_zero == True
    assert Poly(1, x).is_zero == False

    assert Poly(1, x).is_one == True
    assert Poly(2, x).is_one == False

    assert Poly(x-1, x).is_sqf == True
    assert Poly((x-1)**2, x).is_sqf == False

    assert Poly(x-1, x).is_monic == True
    assert Poly(2*x-1, x).is_monic == False

    assert Poly(3*x+2, x).is_primitive == True
    assert Poly(4*x+2, x).is_primitive == False

    assert Poly(1, x).is_ground == True
    assert Poly(x, x).is_ground == False

    assert Poly(x*y*z+1).is_linear == True
    assert Poly(x**2*y*z+1).is_linear == False

    assert Poly(x*y).is_monomial == True
    assert Poly(x*y+1).is_monomial == False

    assert Poly(x*y+x).is_homogeneous == True
    assert Poly(x*y+x+1).is_homogeneous == False

    assert Poly(x).is_univariate == True
    assert Poly(x*y).is_univariate == False

    assert Poly(x*y).is_multivariate == True
    assert Poly(x).is_multivariate == False

def test__polify_basic():
    assert _polify_basic(x-1, x**2-1, x) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(Poly(x-1, x), x**2-1, x) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(x-1, Poly(x**2-1, x), x) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(Poly(x-1, x), Poly(x**2-1, x), x) == (Poly(x-1, x), Poly(x**2-1, x))

    assert _polify_basic(x-1, x**2-1, x, y) == (Poly(x-1, x, y), Poly(x**2-1, x, y))
    assert _polify_basic(Poly(x-1, x), x**2-1, x, y) == (Poly(x-1, x, y), Poly(x**2-1, x, y))
    assert _polify_basic(x-1, Poly(x**2-1, x), x, y) == (Poly(x-1, x, y), Poly(x**2-1, x, y))
    assert _polify_basic(Poly(x-1, x), Poly(x**2-1, x), x, y) == (Poly(x-1, x, y), Poly(x**2-1, x, y))

    assert _polify_basic(x-1, x**2-1) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(Poly(x-1, x), x**2-1) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(x-1, Poly(x**2-1, x)) == (Poly(x-1, x), Poly(x**2-1, x))
    assert _polify_basic(Poly(x-1, x), Poly(x**2-1, x)) == (Poly(x-1, x), Poly(x**2-1, x))

    assert _polify_basic(1, x**2-1) == (Poly(1, x), Poly(x**2-1, x))
    assert _polify_basic(1, x**2-1) == (Poly(1, x), Poly(x**2-1, x))
    assert _polify_basic(1, Poly(x**2-1, x)) == (Poly(1, x), Poly(x**2-1, x))
    assert _polify_basic(1, Poly(x**2-1, x)) == (Poly(1, x), Poly(x**2-1, x))

    assert _polify_basic(x**2-1, 1) == (Poly(x**2-1, x), Poly(1, x))
    assert _polify_basic(x**2-1, 1) == (Poly(x**2-1, x), Poly(1, x))
    assert _polify_basic(Poly(x**2-1, x), 1) == (Poly(x**2-1, x), Poly(1, x))
    assert _polify_basic(Poly(x**2-1, x), 1) == (Poly(x**2-1, x), Poly(1, x))

    raises(CoercionFailed, "_polify_basic(1, 2)")

def test_pdiv():
    f, g = x**2 - y**2, x - y
    q, r = x + y, 0

    F, G, Q, R = [ Poly(h, x, y) for h in (f, g, q, r) ]

    assert F.pdiv(G) == (Q, R)
    assert F.pexquo(G) == Q
    assert F.pquo(G) == Q
    assert F.prem(G) == R

    assert pdiv(f, g) == (q, r)
    assert pexquo(f, g) == q
    assert pquo(f, g) == q
    assert prem(f, g) == r

    assert pdiv(f, g, x, y) == (q, r)
    assert pexquo(f, g, x, y) == q
    assert pquo(f, g, x, y) == q
    assert prem(f, g, x, y) == r

    assert pdiv(f, g, (x,y)) == (q, r)
    assert pexquo(f, g, (x,y)) == q
    assert pquo(f, g, (x,y)) == q
    assert prem(f, g, (x,y)) == r

    assert pdiv(F, G) == (Q, R)
    assert pexquo(F, G) == Q
    assert pquo(F, G) == Q
    assert prem(F, G) == R

    assert pdiv(f, g, polys=True) == (Q, R)
    assert pexquo(f, g, polys=True) == Q
    assert pquo(f, g, polys=True) == Q
    assert prem(f, g, polys=True) == R

    assert pdiv(F, G, polys=False) == (q, r)
    assert pexquo(F, G, polys=False) == q
    assert pquo(F, G, polys=False) == q
    assert prem(F, G, polys=False) == r

    raises(GeneratorsNeeded, "pdiv(4, 2)")
    raises(GeneratorsNeeded, "pexquo(4, 2)")
    raises(GeneratorsNeeded, "pquo(4, 2)")
    raises(GeneratorsNeeded, "prem(4, 2)")

def test_div():
    f, g = x**2 - y**2, x - y
    q, r = x + y, 0

    F, G, Q, R = [ Poly(h, x, y) for h in (f, g, q, r) ]

    assert F.div(G) == (Q, R)
    assert F.exquo(G) == Q
    assert F.quo(G) == Q
    assert F.rem(G) == R

    assert div(f, g) == (q, r)
    assert exquo(f, g) == q
    assert quo(f, g) == q
    assert rem(f, g) == r

    assert div(f, g, x, y) == (q, r)
    assert exquo(f, g, x, y) == q
    assert quo(f, g, x, y) == q
    assert rem(f, g, x, y) == r

    assert div(f, g, (x,y)) == (q, r)
    assert exquo(f, g, (x,y)) == q
    assert quo(f, g, (x,y)) == q
    assert rem(f, g, (x,y)) == r

    assert div(F, G) == (Q, R)
    assert exquo(F, G) == Q
    assert quo(F, G) == Q
    assert rem(F, G) == R

    assert div(f, g, polys=True) == (Q, R)
    assert exquo(f, g, polys=True) == Q
    assert quo(f, g, polys=True) == Q
    assert rem(f, g, polys=True) == R

    assert div(F, G, polys=False) == (q, r)
    assert exquo(F, G, polys=False) == q
    assert quo(F, G, polys=False) == q
    assert rem(F, G, polys=False) == r

    raises(GeneratorsNeeded, "div(4, 2)")
    raises(GeneratorsNeeded, "exquo(4, 2)")
    raises(GeneratorsNeeded, "quo(4, 2)")
    raises(GeneratorsNeeded, "rem(4, 2)")

def test_gcdex():
    f, g = 2*x, x**2 - 16
    s, t, h = x/32, -Rational(1,16), 1

    F, G, S, T, H = [ Poly(u, domain='QQ') for u in (f, g, s, t, h) ]

    assert F.half_gcdex(G) == (S, H)
    assert F.gcdex(G) == (S, T, H)
    assert F.invert(G) == S

    assert half_gcdex(f, g) == (s, h)
    assert gcdex(f, g) == (s, t, h)
    assert invert(f, g) == s

    assert half_gcdex(f, g, x) == (s, h)
    assert gcdex(f, g, x) == (s, t, h)
    assert invert(f, g, x) == s

    assert half_gcdex(f, g, (x,)) == (s, h)
    assert gcdex(f, g, (x,)) == (s, t, h)
    assert invert(f, g, (x,)) == s

    assert half_gcdex(F, G) == (S, H)
    assert gcdex(F, G) == (S, T, H)
    assert invert(F, G) == S

    assert half_gcdex(f, g, polys=True) == (S, H)
    assert gcdex(f, g, polys=True) == (S, T, H)
    assert invert(f, g, polys=True) == S

    assert half_gcdex(F, G, polys=False) == (s, h)
    assert gcdex(F, G, polys=False) == (s, t, h)
    assert invert(F, G, polys=False) == s

    assert half_gcdex(100, 2004) == (-20, 4)
    assert gcdex(100, 2004) == (-20, 1, 4)
    assert invert(3, 7) == 5

def test_subresultants():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 2*x - 2
    F, G, H = Poly(f), Poly(g), Poly(h)

    assert F.subresultants(G) == [F, G, H]
    assert subresultants(f, g) == [f, g, h]
    assert subresultants(f, g, x) == [f, g, h]
    assert subresultants(f, g, (x,)) == [f, g, h]
    assert subresultants(F, G) == [F, G, H]
    assert subresultants(f, g, polys=True) == [F, G, H]
    assert subresultants(F, G, polys=False) == [f, g, h]

    raises(GeneratorsNeeded, "subresultants(4, 2)")

def test_resultant():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 0
    F, G = Poly(f), Poly(g)

    assert F.resultant(G) == h
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(f, g, (x,)) == h
    assert resultant(F, G) == h
    assert resultant(f, g, polys=True) == h
    assert resultant(F, G, polys=False) == h

    f, g, h = x - a, x - b, a - b
    F, G, H = Poly(f), Poly(g), Poly(h)

    assert F.resultant(G) == H
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(f, g, (x,)) == h
    assert resultant(F, G) == H
    assert resultant(f, g, polys=True) == H
    assert resultant(F, G, polys=False) == h

    raises(GeneratorsNeeded, "resultant(4, 2)")

def test_discriminant():
    f, g = x**3 + 3*x**2 + 9*x - 13, -11664
    F = Poly(f)

    assert F.discriminant() == g
    assert discriminant(f) == g
    assert discriminant(f, x) == g
    assert discriminant(f, (x,)) == g
    assert discriminant(F) == g
    assert discriminant(f, polys=True) == g
    assert discriminant(F, polys=False) == g

    f, g = a*x**2 + b*x + c, b**2 - 4*a*c
    F, G = Poly(f), Poly(g)

    assert F.discriminant() == G
    assert discriminant(f) == g
    assert discriminant(f, x, a, b, c) == g
    assert discriminant(f, (x, a, b, c)) == g
    assert discriminant(F) == G
    assert discriminant(f, polys=True) == G
    assert discriminant(F, polys=False) == g

    raises(GeneratorsNeeded, "discriminant(4)")

def test_gcd():
    f, g = x**3 - 1, x**2 - 1
    s, t = x**2 + x + 1, x + 1
    h, r = x - 1, x**4 + x**3 - x - 1

    F, G, S, T, H, R = [ Poly(u) for u in (f, g, s, t, h, r) ]

    assert F.cofactors(G) == (H, S, T)
    assert F.gcd(G) == H
    assert F.lcm(G) == R

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == r

    assert cofactors(f, g, x) == (h, s, t)
    assert gcd(f, g, x) == h
    assert lcm(f, g, x) == r

    assert cofactors(f, g, (x,)) == (h, s, t)
    assert gcd(f, g, (x,)) == h
    assert lcm(f, g, (x,)) == r

    assert cofactors(F, G) == (H, S, T)
    assert gcd(F, G) == H
    assert lcm(F, G) == R

    assert cofactors(f, g, polys=True) == (H, S, T)
    assert gcd(f, g, polys=True) == H
    assert lcm(f, g, polys=True) == R

    assert cofactors(F, G, polys=False) == (h, s, t)
    assert gcd(F, G, polys=False) == h
    assert lcm(F, G, polys=False) == r

    assert cofactors(8, 6) == (2, 4, 3)
    assert gcd(8, 6) == 2
    assert lcm(8, 6) == 24

def test_terms_gcd():
    assert terms_gcd(1) == 1
    assert terms_gcd(1, x) == 1

    assert terms_gcd(x**3*y - x*y**3) == x*y*(x**2 - y**2)
    assert terms_gcd(2*x**3*y - 2*x*y**3) == 2*x*y*(x**2 - y**2)
    assert terms_gcd(x**3*y/2 - x*y**3/2) == x*y/2*(x**2 - y**2)

def test_monic():
    f, g = 2*x - 1, x - S(1)/2
    F, G = Poly(f, domain='QQ'), Poly(g)

    assert F.monic() == G
    assert monic(f) == g
    assert monic(f, x) == g
    assert monic(f, (x,)) == g
    assert monic(F) == G
    assert monic(f, polys=True) == G
    assert monic(F, polys=False) == g

    raises(GeneratorsNeeded, "monic(4)")

def test_content():
    f, F = 4*x + 2, Poly(4*x + 2)

    F.content() == 2
    content(f) == 2

    raises(GeneratorsNeeded, "content(4)")

def test_primitive():
    f, g = 4*x + 2, 2*x + 1
    F, G = Poly(f), Poly(g)

    assert F.primitive() == (2, G)
    assert primitive(f) == (2, g)
    assert primitive(f, x) == (2, g)
    assert primitive(f, (x,)) == (2, g)
    assert primitive(F) == (2, G)
    assert primitive(f, polys=True) == (2, G)
    assert primitive(F, polys=False) == (2, g)

    raises(GeneratorsNeeded, "primitive(4)")

def test_compose():
    f = x**12+20*x**10+150*x**8+500*x**6+625*x**4-2*x**3-10*x+9
    g = x**4 - 2*x + 9
    h = x**3 + 5*x

    F, G, H = map(Poly, (f, g, h))

    assert G.compose(H) == F
    assert compose(g, h) == f
    assert compose(g, h, x) == f
    assert compose(g, h, (x,)) == f
    assert compose(G, H) == F
    assert compose(g, h, polys=True) == F
    assert compose(G, H, polys=False) == f

    assert F.decompose() == [G, H]
    assert decompose(f) == [g, h]
    assert decompose(f, x) == [g, h]
    assert decompose(f, (x,)) == [g, h]
    assert decompose(F) == [G, H]
    assert decompose(f, polys=True) == [G, H]
    assert decompose(F, polys=False) == [g, h]

    raises(GeneratorsNeeded, "compose(4, 2)")
    raises(GeneratorsNeeded, "decompose(4)")

def test_sqf():
    f = x**5 - x**3 - x**2 + 1
    g = x**3 + 2*x**2 + 2*x + 1
    h = x - 1

    p = x**4 + x**3 - x - 1

    F, G, H, P = map(Poly, (f, g, h, p))

    assert F.sqf_part() == P
    assert sqf_part(f) == p
    assert sqf_part(f, x) == p
    assert sqf_part(f, (x,)) == p
    assert sqf_part(F) == P
    assert sqf_part(f, polys=True) == P
    assert sqf_part(F, polys=False) == p

    assert F.sqf_list() == (1, [(G, 1), (H, 2)])
    assert sqf_list(f) == (1, [(g, 1), (h, 2)])
    assert sqf_list(f, x) == (1, [(g, 1), (h, 2)])
    assert sqf_list(f, (x,)) == (1, [(g, 1), (h, 2)])
    assert sqf_list(F) == (1, [(G, 1), (H, 2)])
    assert sqf_list(f, polys=True) == (1, [(G, 1), (H, 2)])
    assert sqf_list(F, polys=False) == (1, [(g, 1), (h, 2)])

    assert sqf_list(f, include=True) == [(g, 1), (h, 2)]

    raises(GeneratorsNeeded, "sqf_part(4)")
    raises(GeneratorsNeeded, "sqf_list(4)")

    assert sqf(1) == 1
    assert sqf(1, frac=True) == 1

    assert sqf(f) == g*h**2
    assert sqf(f, x) == g*h**2
    assert sqf(f, (x,)) == g*h**2

    d = x**2 + y**2

    assert sqf(f/d, frac=True) == (g*h**2)/d
    assert sqf(f/d, x, frac=True) == (g*h**2)/d
    assert sqf(f/d, (x,), frac=True) == (g*h**2)/d

def test_factor():
    f = x**5 - x**3 - x**2 + 1

    u = x + 1
    v = x - 1
    w = x**2 + x + 1

    F, U, V, W = map(Poly, (f, u, v, w))

    assert F.factor_list() == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(f) == (1, [(u, 1), (v, 2), (w, 1)])
    assert factor_list(f, x) == (1, [(u, 1), (v, 2), (w, 1)])
    assert factor_list(f, (x,)) == (1, [(u, 1), (v, 2), (w, 1)])
    assert factor_list(F) == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(f, polys=True) == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(F, polys=False) == (1, [(u, 1), (v, 2), (w, 1)])

    # TODO: assert factor_list(f, include=True) == [(u, 1), (v, 2), (w, 1)]

    raises(GeneratorsNeeded, "factor_list(4)")

    assert factor(1) == 1
    assert factor(1, frac=True) == 1

    assert factor(f) == u*v**2*w
    assert factor(f, x) == u*v**2*w
    assert factor(f, (x,)) == u*v**2*w

    g, p, q = x**2 - y**2, x - y, x + y

    assert factor(f/g, frac=True) == (u*v**2*w)/(p*q)
    assert factor(f/g, x, frac=True) == (u*v**2*w)/(p*q)
    assert factor(f/g, (x,), frac=True) == (u*v**2*w)/(p*q)

def test_sturm():
    f, F = x, Poly(x, domain='QQ')
    g, G = 1, Poly(1, x, domain='QQ')

    assert F.sturm() == [F, G]
    assert sturm(f) == [f, g]
    assert sturm(f, x) == [f, g]
    assert sturm(f, (x,)) == [f, g]
    assert sturm(F) == [F, G]
    assert sturm(f, polys=True) == [F, G]
    assert sturm(F, polys=False) == [f, g]

    raises(GeneratorsNeeded, "sturm(4)")

def test_cancel():
    assert cancel(0) == 0
    assert cancel(7) == 7
    assert cancel(x) == x

    f, g, p, q = 4*x**2-4, 2*x-2, 2*x+2, 1
    F, G, P, Q = [ Poly(u, x) for u in (f, g, p, q) ]

    assert F.cancel(G) == (1, P, Q)
    assert cancel((f, g)) == (1, p, q)
    assert cancel((f, g), x) == (1, p, q)
    assert cancel((f, g), (x,)) == (1, p, q)
    assert cancel((F, G)) == (1, P, Q)
    assert cancel((f, g), polys=True) == (1, P, Q)
    assert cancel((F, G), polys=False) == (1, p, q)

    f = (x**2 - 2)/(x + sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x - sqrt(2)

    f = (x**2 - 2)/(x - sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x + sqrt(2)

    assert cancel((x**2-y)/(x-y)) == 1/(x - y)*(x**2 - y)

    assert cancel((x**2-y**2)/(x-y), x) == x + y
    assert cancel((x**2-y**2)/(x-y), y) == x + y
    assert cancel((x**2-y**2)/(x-y)) == x + y

    assert cancel((x**3-1)/(x**2-1)) == (x**2+x+1)/(x+1)
    assert cancel((x**3/2-S(1)/2)/(x**2-1)) == (x**2+x+1)/(2*x+2)

    assert cancel((exp(2*x) + 2*exp(x) + 1)/(exp(x) + 1)) == exp(x) + 1

def test_groebner():
    pass

