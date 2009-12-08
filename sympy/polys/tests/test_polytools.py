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

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
    PolynomialError,
    CoercionFailed,
    DomainError,
)

from sympy.polys.algebratools import ZZ, QQ, EX

from sympy import S, Integer, Rational, symbols, sqrt, exp, sin, expand, raises, oo

x,y,z,p,q,r,s,t,u,v,w,a,b,c = symbols('x,y,z,p,q,r,s,t,u,v,w,a,b,c')

def test__init_poly_from_poly():
    assert Poly(x*y, x, y).gens == (x, y)
    assert Poly(x*y, y, x).gens == (y, x)

    assert Poly(Poly(x*y, x, y), y, x).gens == (y, x)

def test__init_poly_from_basic():
    assert Poly(0).is_Poly == False
    assert Poly(1).is_Poly == False

    assert Poly(0, x).is_Poly == True
    assert Poly(1, x).is_Poly == True

def test_Poly__args():
    assert Poly(x**2 + 1).args == [x**2 + 1]

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

def test_Poly__analyze_domain():
    assert Poly._analyze_domain({}) is None
    assert Poly._analyze_domain({'domain': ZZ}) == ZZ
    assert Poly._analyze_domain({'domain': 'ZZ'}) == ZZ

def test_Poly__parse_domain():
    assert Poly._parse_domain(ZZ) == ZZ
    assert Poly._parse_domain(QQ) == QQ
    assert Poly._parse_domain(EX) == EX
    assert Poly._parse_domain(ZZ[x,y]) == ZZ[x,y]

    assert Poly._parse_domain('Z') == ZZ
    assert Poly._parse_domain('Q') == QQ

    assert Poly._parse_domain('ZZ') == ZZ
    assert Poly._parse_domain('QQ') == QQ

    assert Poly._parse_domain('EX') == EX

    raises(ValueError, "Poly._parse_domain('Z[]')")

    assert Poly._parse_domain('Z[x]') == ZZ[x]
    assert Poly._parse_domain('Q[x]') == QQ[x]

    assert Poly._parse_domain('ZZ[x]') == ZZ[x]
    assert Poly._parse_domain('QQ[x]') == QQ[x]

    assert Poly._parse_domain('Z[x,y]') == ZZ[x,y]
    assert Poly._parse_domain('Q[x,y]') == QQ[x,y]

    assert Poly._parse_domain('ZZ[x,y]') == ZZ[x,y]
    assert Poly._parse_domain('QQ[x,y]') == QQ[x,y]

    raises(ValueError, "Poly._parse_domain('Z()')")

    assert Poly._parse_domain('Z(x)') == ZZ.frac_field(x)
    assert Poly._parse_domain('Q(x)') == QQ.frac_field(x)

    assert Poly._parse_domain('ZZ(x)') == ZZ.frac_field(x)
    assert Poly._parse_domain('QQ(x)') == QQ.frac_field(x)

    assert Poly._parse_domain('Z(x,y)') == ZZ.frac_field(x,y)
    assert Poly._parse_domain('Q(x,y)') == QQ.frac_field(x,y)

    assert Poly._parse_domain('ZZ(x,y)') == ZZ.frac_field(x,y)
    assert Poly._parse_domain('QQ(x,y)') == QQ.frac_field(x,y)

def test_Poly_get_domain():
    assert Poly(2*x).get_domain() == ZZ

    assert Poly(2*x, domain='ZZ').get_domain() == ZZ
    assert Poly(2*x, domain='QQ').get_domain() == QQ

    assert Poly(x/2).get_domain() == QQ

    raises(CoercionFailed, "Poly(x/2, domain='ZZ')")
    assert Poly(x/2, domain='QQ').get_domain() == QQ

def test_Poly_set_domain():
    assert Poly(2*x + 1).set_domain(ZZ) == Poly(2*x + 1)
    assert Poly(2*x + 1).set_domain('ZZ') == Poly(2*x + 1)

    assert Poly(2*x + 1).set_domain(QQ) == Poly(2*x + 1, domain='QQ')
    assert Poly(2*x + 1).set_domain('QQ') == Poly(2*x + 1, domain='QQ')

    raises(CoercionFailed, "Poly(x/2 + 1).set_domain(ZZ)")
    raises(DomainError, "Poly(x + 1, modulus=2).set_domain(QQ)")

def test_Poly__analyze_modulus():
    assert Poly._analyze_modulus({}) is None
    assert Poly._analyze_modulus({'modulus': 2}) == 2
    assert Poly._analyze_modulus({'modulus': Integer(2)}) == 2

def test_Poly__parse_modulus():
    assert Poly._parse_modulus(5) == 5
    assert Poly._parse_modulus(Integer(5)) == 5

    raises(ValueError, "Poly._parse_modulus(1)")
    raises(ValueError, "Poly._parse_modulus(x)")

def test_Poly_get_modulus():
    Poly(x**2 + 1, modulus=2).get_modulus() == 2
    raises(PolynomialError, "Poly(x**2 + 1).get_modulus()")

def test_Poly_set_modulus():
    Poly(x**2 + 1, modulus=2).set_modulus(7) == Poly(x**2 + 1, modulus=7)
    Poly(x**2 + 5, modulus=7).set_modulus(2) == Poly(x**2 + 1, modulus=2)

    Poly(x**2 + 1).set_modulus(2) == Poly(x**2 + 1, modulus=2)

    raises(PolynomialError, "Poly(x/2 + 1).set_modulus(2)")

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

def test_Poly_to_ring():
    assert Poly(2*x+1, domain='ZZ').to_ring() == Poly(2*x+1, domain='ZZ')
    assert Poly(2*x+1, domain='QQ').to_ring() == Poly(2*x+1, domain='ZZ')

    raises(CoercionFailed, "Poly(x/2+1).to_ring()")

def test_Poly_to_field():
    assert Poly(2*x+1, domain='ZZ').to_field() == Poly(2*x+1, domain='QQ')
    assert Poly(2*x+1, domain='QQ').to_field() == Poly(2*x+1, domain='QQ')

    assert Poly(x/2+1, domain='QQ').to_field() == Poly(x/2+1, domain='QQ')

def test_Poly_coeffs():
    assert Poly(0, x).coeffs() == [0]
    assert Poly(1, x).coeffs() == [1]

    assert Poly(2*x+1, x).coeffs() == [2,1]

    assert Poly(7*x**2+2*x+1, x).coeffs() == [7,2,1]
    assert Poly(7*x**4+2*x+1, x).coeffs() == [7,2,1]

def test_Poly_monoms():
    assert Poly(0, x).monoms() == [(0,)]
    assert Poly(1, x).monoms() == [(0,)]

    assert Poly(2*x+1, x).monoms() == [(1,),(0,)]

    assert Poly(7*x**2+2*x+1, x).monoms() == [(2,),(1,),(0,)]
    assert Poly(7*x**4+2*x+1, x).monoms() == [(4,),(1,),(0,)]

def test_Poly_terms():
    assert Poly(0, x).terms() == [((0,), 0)]
    assert Poly(1, x).terms() == [((0,), 1)]

    assert Poly(2*x+1, x).terms() == [((1,), 2),((0,), 1)]

    assert Poly(7*x**2+2*x+1, x).terms() == [((2,), 7),((1,), 2),((0,), 1)]
    assert Poly(7*x**4+2*x+1, x).terms() == [((4,), 7),((1,), 2),((0,), 1)]

def test_Poly_all_coeffs():
    assert Poly(0, x).all_coeffs() == [0]
    assert Poly(1, x).all_coeffs() == [1]

    assert Poly(2*x+1, x).all_coeffs() == [2,1]

    assert Poly(7*x**2+2*x+1, x).all_coeffs() == [7,2,1]
    assert Poly(7*x**4+2*x+1, x).all_coeffs() == [7,0,0,2,1]

def test_Poly_all_monoms():
    assert Poly(0, x).all_monoms() == [(0,)]
    assert Poly(1, x).all_monoms() == [(0,)]

    assert Poly(2*x+1, x).all_monoms() == [(1,),(0,)]

    assert Poly(7*x**2+2*x+1, x).all_monoms() == [(2,),(1,),(0,)]
    assert Poly(7*x**4+2*x+1, x).all_monoms() == [(4,),(3,),(2,),(1,),(0,)]

def test_Poly_all_terms():
    assert Poly(0, x).all_terms() == [((0,), 0)]
    assert Poly(1, x).all_terms() == [((0,), 1)]

    assert Poly(2*x+1, x).all_terms() == [((1,), 2),((0,), 1)]

    assert Poly(7*x**2+2*x+1, x).all_terms() == [((2,), 7),((1,), 2),((0,), 1)]
    assert Poly(7*x**4+2*x+1, x).all_terms() == [((4,), 7),((3,),0),((2,),0),((1,), 2),((0,), 1)]

def test_Poly_length():
    assert Poly(0, x).length() == 0
    assert Poly(1, x).length() == 1
    assert Poly(x, x).length() == 1

    assert Poly(x+1, x).length() == 2
    assert Poly(x**2+1, x).length() == 2
    assert Poly(x**2+x+1, x).length() == 3

def test_Poly_as_dict():
    assert Poly(0, x).as_dict() == {}
    assert Poly(0, x, y, z).as_dict() == {}

    assert Poly(1, x).as_dict() == {(0,): 1}
    assert Poly(1, x, y, z).as_dict() == {(0,0,0): 1}

    assert Poly(x**2+3, x).as_dict() == {(2,): 1, (0,): 3}
    assert Poly(x**2+3, x, y, z).as_dict() == {(2,0,0): 1, (0,0,0): 3}

    assert Poly(3*x**2*y*z**3+4*x*y+5*x*z).as_dict() == {(2,1,3): 3, (1,1,0): 4, (1,0,1): 5}

def test_Poly_as_basic():
    assert Poly(0, x).as_basic() == 0
    assert Poly(0, x, y, z).as_basic() == 0

    assert Poly(1, x).as_basic() == 1
    assert Poly(1, x, y, z).as_basic() == 1

    assert Poly(x**2+3, x).as_basic() == x**2+3
    assert Poly(x**2+3, x, y, z).as_basic() == x**2+3

    assert Poly(3*x**2*y*z**3+4*x*y+5*x*z).as_basic() == 3*x**2*y*z**3+4*x*y+5*x*z

def test_Poly_deflate():
    assert Poly(0, x).deflate() == ((1,), Poly(0, x))
    assert Poly(1, x).deflate() == ((1,), Poly(1, x))
    assert Poly(x, x).deflate() == ((1,), Poly(x, x))

    assert Poly(x**2, x).deflate() == ((2,), Poly(x, x))
    assert Poly(x**17, x).deflate() == ((17,), Poly(x, x))

    assert Poly(x**2*y*z**11+x**4*z**11).deflate() == ((2,1,11), Poly(x*y*z+x**2*z))

def test_Poly__gen_to_level():
    assert Poly(1, x, y)._gen_to_level(-2) == 0
    assert Poly(1, x, y)._gen_to_level(-1) == 1
    assert Poly(1, x, y)._gen_to_level( 0) == 0
    assert Poly(1, x, y)._gen_to_level( 1) == 1

    raises(PolynomialError, "Poly(1, x, y)._gen_to_level(-3)")
    raises(PolynomialError, "Poly(1, x, y)._gen_to_level( 2)")

    assert Poly(1, x, y)._gen_to_level(x) == 0
    assert Poly(1, x, y)._gen_to_level(y) == 1

    assert Poly(1, x, y)._gen_to_level('x') == 0
    assert Poly(1, x, y)._gen_to_level('y') == 1

    raises(PolynomialError, "Poly(1, x, y)._gen_to_level(z)")
    raises(PolynomialError, "Poly(1, x, y)._gen_to_level('z')")

def test_Poly_degree():
    assert Poly(0, x).degree() ==-1
    assert Poly(1, x).degree() == 0
    assert Poly(x, x).degree() == 1

    assert Poly(0, x).degree(gen=0) ==-1
    assert Poly(1, x).degree(gen=0) == 0
    assert Poly(x, x).degree(gen=0) == 1

    assert Poly(0, x).degree(gen=x) ==-1
    assert Poly(1, x).degree(gen=x) == 0
    assert Poly(x, x).degree(gen=x) == 1

    assert Poly(0, x).degree(gen='x') ==-1
    assert Poly(1, x).degree(gen='x') == 0
    assert Poly(x, x).degree(gen='x') == 1

    raises(PolynomialError, "Poly(1, x).degree(gen=1)")
    raises(PolynomialError, "Poly(1, x).degree(gen=y)")
    raises(PolynomialError, "Poly(1, x).degree(gen='y')")

    assert Poly(1, x, y).degree() == 0
    assert Poly(2*y, x, y).degree() == 0
    assert Poly(x*y, x, y).degree() == 1

    assert Poly(1, x, y).degree(gen=x) == 0
    assert Poly(2*y, x, y).degree(gen=x) == 0
    assert Poly(x*y, x, y).degree(gen=x) == 1

    assert Poly(1, x, y).degree(gen=y) == 0
    assert Poly(2*y, x, y).degree(gen=y) == 1
    assert Poly(x*y, x, y).degree(gen=y) == 1

def test_Poly_degree_list():
    assert Poly(0, x).degree_list() == (-1,)
    assert Poly(0, x, y).degree_list() == (-1,-1)
    assert Poly(0, x, y, z).degree_list() == (-1,-1,-1)

    assert Poly(1, x).degree_list() == (0,)
    assert Poly(1, x, y).degree_list() == (0,0)
    assert Poly(1, x, y, z).degree_list() == (0,0,0)

    assert Poly(x**2*y+x**3*z**2+1).degree_list() == (3,1,2)

def test_Poly_total_degree():
    assert Poly(x**2*y+x**3*z**2+1).total_degree() == 6

def test_Poly_LC():
    assert Poly(0, x).LC() == 0
    assert Poly(1, x).LC() == 1
    assert Poly(2*x**2+x, x).LC() == 2

def test_Poly_TC():
    assert Poly(0, x).TC() == 0
    assert Poly(1, x).TC() == 1
    assert Poly(2*x**2+x, x).TC() == 0

def test_Poly_EC():
    assert Poly(0, x).EC() == 0
    assert Poly(1, x).EC() == 1
    assert Poly(2*x**2+x, x).EC() == 1

def test_Poly_nth():
    assert Poly(0, x).nth(0) == 0
    assert Poly(0, x).nth(1) == 0

    assert Poly(1, x).nth(0) == 1
    assert Poly(1, x).nth(1) == 0

    assert Poly(x**8, x).nth(0) == 0
    assert Poly(x**8, x).nth(7) == 0
    assert Poly(x**8, x).nth(8) == 1
    assert Poly(x**8, x).nth(9) == 0

    assert Poly(3*x*y**2 + 1).nth(0, 0) == 1
    assert Poly(3*x*y**2 + 1).nth(1, 2) == 3

def test_Poly_LM():
    assert Poly(0, x).LM() == (0,)
    assert Poly(1, x).LM() == (0,)
    assert Poly(2*x**2+x, x).LM() == (2,)

def test_Poly_EM():
    assert Poly(0, x).EM() == (0,)
    assert Poly(1, x).EM() == (0,)
    assert Poly(2*x**2+x, x).EM() == (1,)

def test_Poly_LT():
    assert Poly(0, x).LT() == ((0,), 0)
    assert Poly(1, x).LT() == ((0,), 1)
    assert Poly(2*x**2+x, x).LT() == ((2,), 2)

def test_Poly_ET():
    assert Poly(0, x).ET() == ((0,), 0)
    assert Poly(1, x).ET() == ((0,), 1)
    assert Poly(2*x**2+x, x).ET() == ((1,), 1)

def test_Poly_max_norm():
    assert Poly(-1, x).max_norm() == 1
    assert Poly( 0, x).max_norm() == 0
    assert Poly( 1, x).max_norm() == 1

def test_Poly_l1_norm():
    assert Poly(-1, x).l1_norm() == 1
    assert Poly( 0, x).l1_norm() == 0
    assert Poly( 1, x).l1_norm() == 1

def test_Poly_ground_to_ring():
    assert Poly(2*x + 1).ground_to_ring() == (1, Poly(2*x + 1, domain='ZZ'))
    assert Poly(x/2 + 1).ground_to_ring() == (2, Poly(x + 2, domain='QQ'))

def test_Poly_integrate():
    assert Poly(x + 1).integrate() == Poly(x**2/2 + x)
    assert Poly(x + 1).integrate(x) == Poly(x**2/2 + x)
    assert Poly(x + 1).integrate((x, 1)) == Poly(x**2/2 + x)

    assert Poly(x*y + 1).integrate(x) == Poly(x**2*y/2 + x)
    assert Poly(x*y + 1).integrate(y) == Poly(x*y**2/2 + y)

    assert Poly(x*y + 1).integrate(x, x) == Poly(x**3*y/6 + x**2/2)
    assert Poly(x*y + 1).integrate(y, y) == Poly(x*y**3/6 + y**2/2)

    assert Poly(x*y + 1).integrate((x, 2)) == Poly(x**3*y/6 + x**2/2)
    assert Poly(x*y + 1).integrate((y, 2)) == Poly(x*y**3/6 + y**2/2)

    assert Poly(x*y + 1).integrate(x, y) == Poly(x**2*y**2/4 + x*y)
    assert Poly(x*y + 1).integrate(y, x) == Poly(x**2*y**2/4 + x*y)

def test_Poly_diff():
    assert Poly(x**2 + x).diff() == Poly(2*x + 1)
    assert Poly(x**2 + x).diff(x) == Poly(2*x + 1)
    assert Poly(x**2 + x).diff((x, 1)) == Poly(2*x + 1)

    assert Poly(x**2*y**2 + x*y).diff(x) == Poly(2*x*y**2 + y)
    assert Poly(x**2*y**2 + x*y).diff(y) == Poly(2*x**2*y + x)

    assert Poly(x**2*y**2 + x*y).diff(x, x) == Poly(2*y**2, x, y)
    assert Poly(x**2*y**2 + x*y).diff(y, y) == Poly(2*x**2, x, y)

    assert Poly(x**2*y**2 + x*y).diff((x, 2)) == Poly(2*y**2, x, y)
    assert Poly(x**2*y**2 + x*y).diff((y, 2)) == Poly(2*x**2, x, y)

    assert Poly(x**2*y**2 + x*y).diff(x, y) == Poly(4*x*y + 1)
    assert Poly(x**2*y**2 + x*y).diff(y, x) == Poly(4*x*y + 1)

def test_Poly_eval():
    assert Poly(0, x).eval(7) == 0
    assert Poly(1, x).eval(7) == 1
    assert Poly(x, x).eval(7) == 7

    assert Poly(0, x).eval(7, gen=0) == 0
    assert Poly(1, x).eval(7, gen=0) == 1
    assert Poly(x, x).eval(7, gen=0) == 7

    assert Poly(0, x).eval(7, gen=x) == 0
    assert Poly(1, x).eval(7, gen=x) == 1
    assert Poly(x, x).eval(7, gen=x) == 7

    assert Poly(0, x).eval(7, gen='x') == 0
    assert Poly(1, x).eval(7, gen='x') == 1
    assert Poly(x, x).eval(7, gen='x') == 7

    raises(PolynomialError, "Poly(1, x).eval(7, gen=1)")
    raises(PolynomialError, "Poly(1, x).eval(7, gen=y)")
    raises(PolynomialError, "Poly(1, x).eval(7, gen='y')")

    assert Poly(1, x, y).eval(7) == Poly(1, y)
    assert Poly(2*y, x, y).eval(7) == Poly(2*y, y)
    assert Poly(x*y, x, y).eval(7) == Poly(7*y, y)

    assert Poly(1, x, y).eval(7, gen=x) == Poly(1, y)
    assert Poly(2*y, x, y).eval(7, gen=x) == Poly(2*y, y)
    assert Poly(x*y, x, y).eval(7, gen=x) == Poly(7*y, y)

    assert Poly(1, x, y).eval(7, gen=y) == Poly(1, x)
    assert Poly(2*y, x, y).eval(7, gen=y) == Poly(14, x)
    assert Poly(x*y, x, y).eval(7, gen=y) == Poly(7*x, x)

    raises(CoercionFailed, "Poly(x+1, domain='ZZ').eval(S(1)/2)")

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

    f = Poly(sin(1)*x + 1, x, domain=EX)

    assert f.factor_list() == (1, [(f, 1)])

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

    assert cancel(oo) == oo

    assert cancel((1, 0)) == (1, 1, 0)
    assert cancel((0, 1)) == (1, 0, 1)

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

    assert cancel((x**2/4 - 1, x/2 - 1)) == (S(1)/2, x + 2, 1)

    assert cancel((x**2-y)/(x-y)) == 1/(x - y)*(x**2 - y)

    assert cancel((x**2-y**2)/(x-y), x) == x + y
    assert cancel((x**2-y**2)/(x-y), y) == x + y
    assert cancel((x**2-y**2)/(x-y)) == x + y

    assert cancel((x**3-1)/(x**2-1)) == (x**2+x+1)/(x+1)
    assert cancel((x**3/2-S(1)/2)/(x**2-1)) == (x**2+x+1)/(2*x+2)

    assert cancel((exp(2*x) + 2*exp(x) + 1)/(exp(x) + 1)) == exp(x) + 1

def test_groebner():
    pass

