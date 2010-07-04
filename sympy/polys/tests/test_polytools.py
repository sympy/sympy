"""Tests for user-friendly public interface to polynomial functions. """

from sympy.polys.polytools import (
    Poly, poly, _polify_basic,
    _construct_domain,
    _init_poly_from_dict,
    _init_poly_from_list,
    _init_poly_from_poly,
    _init_poly_from_basic,
    degree, degree_list,
    LC, LM, LT,
    pdiv, prem, pquo, pexquo,
    div, rem, quo, exquo,
    half_gcdex, gcdex, invert,
    subresultants,
    resultant, discriminant,
    cofactors, gcd, lcm, terms_gcd,
    trunc,
    monic, content, primitive,
    compose, decompose,
    sturm,
    sqf_norm, sqf_part, sqf_list, sqf,
    factor_list, factor,
    intervals, refine_root, count_roots,
    real_roots, nroots,
    cancel,
    reduced, groebner,
    symmetrize,
    horner,
)

from sympy.polys.polyerrors import (
    OperationNotSupported,
    ExactQuotientFailed,
    UnificationFailed,
    RefinementFailed,
    GeneratorsNeeded,
    GeneratorsError,
    PolynomialError,
    CoercionFailed,
    NotAlgebraic,
    DomainError,
    OptionError,
)

from sympy.polys.monomialtools import (
    monomial_lex_key,
)

from sympy.polys.polyclasses import GFP, DMP, DMF

from sympy.polys.algebratools import ZZ, QQ, RR, EX

from sympy import (
    S, Integer, Rational, Mul, symbols, sqrt, exp,
    sin, expand, raises, oo, I, pi, re, im, RootOf,
)

from sympy.utilities import all

x,y,z,p,q,r,s,t,u,v,w,a,b,c,d,e = symbols('x,y,z,p,q,r,s,t,u,v,w,a,b,c,d,e')

def _eq(a, b):
    for x, y in zip(a, b):
        if abs(x-y) > 1e-10:
            return False
    return True

def test__construct_domain():
    assert _construct_domain({(0,): 1, (1,): 2}) == \
        (ZZ, {(0,): ZZ(1), (1,): ZZ(2)})
    assert _construct_domain({(0,): 1, (1,): 2}, field=True) == \
        (QQ, {(0,): QQ(1), (1,): QQ(2)})

    assert _construct_domain({(0,): S(1), (1,): S(2)}) == \
        (ZZ, {(0,): ZZ(1), (1,): ZZ(2)})
    assert _construct_domain({(0,): S(1), (1,): S(2)}, field=True) == \
        (QQ, {(0,): QQ(1), (1,): QQ(2)})

    assert _construct_domain({(0,): S(1)/2, (1,): S(2)}) == \
        (QQ, {(0,): QQ(1,2), (1,): QQ(2)})

    assert _construct_domain({(0,): 3.14, (1,): 1, (2,): S(1)/2}) == \
        (RR, {(0,): RR(3.14), (1,): RR(1.0), (2,): RR(0.5)})

    assert _construct_domain({(0,): 3.14, (1,): sqrt(2)}, extension=None) == \
        (EX, {(0,): EX(3.14), (1,): EX(sqrt(2))})
    assert _construct_domain({(0,): 3.14, (1,): sqrt(2)}, extension=True) == \
        (EX, {(0,): EX(3.14), (1,): EX(sqrt(2))})

    assert _construct_domain({(0,): 1, (1,): sqrt(2)}, extension=None) == \
        (EX, {(0,): EX(1), (1,): EX(sqrt(2))})

    ALG = QQ.algebraic_field(sqrt(2))

    assert _construct_domain({(0,): 7, (1,): S(1)/2, (2,): sqrt(2)}, extension=True) == \
        (ALG, {(0,): ALG.convert(7), (1,): ALG.convert(S(1)/2), (2,): ALG.convert(sqrt(2))})

    ALG = QQ.algebraic_field(sqrt(2)+sqrt(3))

    assert _construct_domain({(0,): 7, (1,): sqrt(2), (2,): sqrt(3)}, extension=True) == \
        (ALG, {(0,): ALG.convert(7), (1,): ALG.convert(sqrt(2)), (2,): ALG.convert(sqrt(3))})

    assert _construct_domain({(0,): 2*x, (1,): 3}) == \
        (ZZ[x], {(0,): DMP([2,0], ZZ), (1,): DMP([3], ZZ)})
    assert _construct_domain({(0,): 2*x, (1,): 3*y}) == \
        (ZZ[x,y], {(0,): DMP([[2],[]], ZZ), (1,): DMP([[3,0]], ZZ)})

    assert _construct_domain({(0,): x/2, (1,): 3}) == \
        (QQ[x], {(0,): DMP([QQ(1,2),QQ(0)], QQ), (1,): DMP([QQ(3)], QQ)})
    assert _construct_domain({(0,): x/2, (1,): 3*y}) == \
        (QQ[x,y], {(0,): DMP([[QQ(1,2)],[]], QQ), (1,): DMP([[QQ(3),QQ(0)]], QQ)})

    assert _construct_domain({(0,): 2/x, (1,): 3}) == \
        (ZZ.frac_field(x), {(0,): DMF(([2], [1,0]), ZZ), (1,): DMF(([3], [1]), ZZ)})
    assert _construct_domain({(0,): 2/x, (1,): 3*y}) == \
        (ZZ.frac_field(x,y), {(0,): DMF(([[2]], [[1],[]]), ZZ), (1,): DMF(([[3,0]], [[1]]), ZZ)})

    assert _construct_domain({(1,): sin(y)}, composite=False) == \
        (EX, {(1,): EX(sin(y))})
    assert _construct_domain({(1,): y}, composite=False) == \
        (EX, {(1,): EX(y)})
    assert _construct_domain({(1, 1): 1}, composite=False) == \
        (ZZ, {(1, 1): 1})
    assert _construct_domain({(1, 0): y}, composite=False) == \
        (EX, {(1, 0): EX(y)})


def test__init_poly_from_dict():
    raises(PolynomialError, "_init_poly_from_dict({0: 1, 1: 2}, x, y, modulus=3, domain=ZZ)")

    assert _init_poly_from_dict({0: 1, 1: 2}, x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)
    assert _init_poly_from_dict({0: 1, 1: 5}, x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)

    assert _init_poly_from_dict({(0,): 1, (1,): 2}, x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)
    assert _init_poly_from_dict({(0,): 1, (1,): 5}, x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)

    raises(DomainError, "_init_poly_from_dict({0: 1, 1: 2}, x, modulus=3, domain=QQ)")

    assert _init_poly_from_dict({0: 1, 1: 2}, x) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_dict({0: 1, 1: 2}, x, field=True) == DMP([QQ(2),QQ(1)], QQ)

    assert _init_poly_from_dict({0: 1, 1: 2}, x, domain=ZZ) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_dict({0: 1, 1: 2}, x, domain=QQ) == DMP([QQ(2),QQ(1)], QQ)

    assert _init_poly_from_dict({(0,): 1, (1,): 2}, x) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_dict({(0,): 1, (1,): 2}, x, field=True) == DMP([QQ(2),QQ(1)], QQ)

    assert _init_poly_from_dict({(0,): 1, (1,): 2}, x, domain=ZZ) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_dict({(0,): 1, (1,): 2}, x, domain=QQ) == DMP([QQ(2),QQ(1)], QQ)

def test__init_poly_from_list():
    raises(PolynomialError, "_init_poly_from_list([[]], x, y)")

    assert _init_poly_from_list([2,1], x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)
    assert _init_poly_from_list([5,1], x, modulus=3, domain=ZZ) == GFP([2,1], 3, ZZ)

    assert _init_poly_from_list([2,1], x) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_list([2,1], x, field=True) == DMP([QQ(2),QQ(1)], QQ)

    assert _init_poly_from_list([2,1], x, domain=ZZ) == DMP([ZZ(2),ZZ(1)], ZZ)
    assert _init_poly_from_list([2,1], x, domain=QQ) == DMP([QQ(2),QQ(1)], QQ)

def test__init_poly_from_poly():
    f = Poly(x+7, x, domain=ZZ)
    g = Poly(x+2, x, modulus=3)
    h = Poly(x+y, x, y, domain=ZZ)

    assert _init_poly_from_poly(f) == f
    assert _init_poly_from_poly(f, domain=ZZ) == (DMP([1,7], ZZ), (x,))
    assert _init_poly_from_poly(f, domain=QQ) == (DMP([1,7], QQ), (x,))
    assert _init_poly_from_poly(f, modulus=5) == (GFP([1,2], 5, ZZ), (x,))

    assert _init_poly_from_poly(f, x) == f
    assert _init_poly_from_poly(f, x, domain=ZZ) == (DMP([1,7], ZZ), (x,))
    assert _init_poly_from_poly(f, x, domain=QQ) == (DMP([1,7], QQ), (x,))
    assert _init_poly_from_poly(f, x, modulus=5) == (GFP([1,2], 5, ZZ), (x,))

    assert _init_poly_from_poly(f, y,) == Poly(x + 7, y, domain='ZZ[x]')
    raises(CoercionFailed, "_init_poly_from_poly(f, y, domain=ZZ)")
    raises(CoercionFailed, "_init_poly_from_poly(f, y, domain=QQ)")
    raises(CoercionFailed, "_init_poly_from_poly(f, y, modulus=5)")

    assert _init_poly_from_poly(f, x, y) == Poly(x + 7, x, y, domain='ZZ')
    assert _init_poly_from_poly(f, x, y, domain=ZZ) == Poly(x + 7, x, y, domain='ZZ')
    assert _init_poly_from_poly(f, x, y, domain=QQ) == Poly(x + 7, x, y, domain='QQ')
    raises(PolynomialError, "_init_poly_from_poly(f, x, y, modulus=3)")

    assert _init_poly_from_poly(g) == g
    assert _init_poly_from_poly(g, domain=ZZ) == (GFP([1,2], 3, ZZ), (x,))
    raises(DomainError, "_init_poly_from_poly(g, domain=QQ)")
    assert _init_poly_from_poly(g, modulus=2) == (GFP([1,0], 2, ZZ), (x,))

    assert _init_poly_from_poly(g, x) == g
    assert _init_poly_from_poly(g, x, domain=ZZ) == (GFP([1,2], 3, ZZ), (x,))
    raises(DomainError, "_init_poly_from_poly(g, x, domain=QQ)")
    assert _init_poly_from_poly(g, x, modulus=2) == (GFP([1,0], 2, ZZ), (x,))

    raises(PolynomialError, "_init_poly_from_poly(g, y)")
    raises(PolynomialError, "_init_poly_from_poly(g, y, domain=ZZ)")
    raises(PolynomialError, "_init_poly_from_poly(g, y, domain=QQ)")
    raises(PolynomialError, "_init_poly_from_poly(g, y, modulus=3)")

    raises(PolynomialError, "_init_poly_from_poly(g, x, y)")
    raises(PolynomialError, "_init_poly_from_poly(g, x, y, domain=ZZ)")
    raises(PolynomialError, "_init_poly_from_poly(g, x, y, domain=QQ)")
    raises(PolynomialError, "_init_poly_from_poly(g, x, y, modulus=3)")

    assert _init_poly_from_poly(h) == h
    assert _init_poly_from_poly(h, domain=ZZ) == (DMP([[ZZ(1)],[ZZ(1),ZZ(0)]], ZZ), (x,y))
    assert _init_poly_from_poly(h, domain=QQ) == (DMP([[QQ(1)],[QQ(1),QQ(0)]], QQ), (x,y))
    raises(PolynomialError, "_init_poly_from_poly(h, modulus=3)")

    assert _init_poly_from_poly(h, x) == Poly(x+y, x, domain=ZZ[y])
    raises(CoercionFailed, "_init_poly_from_poly(h, x, domain=ZZ)")
    assert _init_poly_from_poly(h, x, domain=ZZ[y]) == Poly(x+y, x, domain=ZZ[y])
    raises(CoercionFailed, "_init_poly_from_poly(h, x, domain=QQ)")
    assert _init_poly_from_poly(h, x, domain=QQ[y]) == Poly(x+y, x, domain=QQ[y])
    raises(CoercionFailed, "_init_poly_from_poly(h, x, modulus=3)")

    assert _init_poly_from_poly(h, y) == Poly(x+y, y, domain=ZZ[x])
    raises(CoercionFailed, "_init_poly_from_poly(h, y, domain=ZZ)")
    assert _init_poly_from_poly(h, y, domain=ZZ[x]) == Poly(x+y, y, domain=ZZ[x])
    raises(CoercionFailed, "_init_poly_from_poly(h, y, domain=QQ)")
    assert _init_poly_from_poly(h, y, domain=QQ[x]) == Poly(x+y, y, domain=QQ[x])
    raises(CoercionFailed, "_init_poly_from_poly(h, y, modulus=3)")

    assert _init_poly_from_poly(h, x, y) == h
    assert _init_poly_from_poly(h, x, y, domain=ZZ) == (DMP([[ZZ(1)],[ZZ(1),ZZ(0)]], ZZ), (x,y))
    assert _init_poly_from_poly(h, x, y, domain=QQ) == (DMP([[QQ(1)],[QQ(1),QQ(0)]], QQ), (x,y))
    raises(PolynomialError, "_init_poly_from_poly(h, x, y, modulus=3)")

    assert _init_poly_from_poly(h, y, x) == (DMP([[ZZ(1)],[ZZ(1),ZZ(0)]], ZZ), (y, x))
    assert _init_poly_from_poly(h, y, x, domain=ZZ) == (DMP([[ZZ(1)],[ZZ(1),ZZ(0)]], ZZ), (y, x))
    assert _init_poly_from_poly(h, y, x, domain=QQ) == (DMP([[QQ(1)],[QQ(1),QQ(0)]], QQ), (y, x))
    raises(PolynomialError, "_init_poly_from_poly(h, y, x, modulus=3)")

    assert _init_poly_from_poly(h, x, y, field=True) == (DMP([[QQ(1)],[QQ(1),QQ(0)]], QQ), (x, y))
    assert _init_poly_from_poly(h, x, y, field=True) == (DMP([[QQ(1)],[QQ(1),QQ(0)]], QQ), (x, y))

def test__init_poly_from_basic():
    assert _init_poly_from_basic(S(0)) == 0
    assert _init_poly_from_basic(S(7)) == 7

    assert _init_poly_from_basic(x + 5, modulus=3, domain=ZZ) == (GFP([1,2], 3, ZZ), (x,))
    assert _init_poly_from_basic(y + 5, modulus=3, domain=ZZ) == (GFP([1,2], 3, ZZ), (y,))

    assert _init_poly_from_basic(x + 5, x, modulus=3, domain=ZZ) == (GFP([1,2], 3, ZZ), (x,))
    assert _init_poly_from_basic(y + 5, y, modulus=3, domain=ZZ) == (GFP([1,2], 3, ZZ), (y,))

    raises(PolynomialError, "_init_poly_from_basic(x + y, modulus=3, domain=ZZ)")
    raises(PolynomialError, "_init_poly_from_basic(x + y, x, y, modulus=3, domain=ZZ)")

    assert _init_poly_from_basic(x + 5) == (DMP([1,5], ZZ), (x,))
    assert _init_poly_from_basic(y + 5) == (DMP([1,5], ZZ), (y,))

    assert _init_poly_from_basic(x + 5, x) == (DMP([1,5], ZZ), (x,))
    assert _init_poly_from_basic(y + 5, y) == (DMP([1,5], ZZ), (y,))

    assert _init_poly_from_basic(x + 5, domain=ZZ) == (DMP([1,5], ZZ), (x,))
    assert _init_poly_from_basic(y + 5, domain=ZZ) == (DMP([1,5], ZZ), (y,))

    assert _init_poly_from_basic(x + 5, x, domain=ZZ) == (DMP([1,5], ZZ), (x,))
    assert _init_poly_from_basic(y + 5, y, domain=ZZ) == (DMP([1,5], ZZ), (y,))

    assert _init_poly_from_basic(x + 5, x, y, domain=ZZ) == (DMP([[1],[5]], ZZ), (x,y))
    assert _init_poly_from_basic(y + 5, x, y, domain=ZZ) == (DMP([[1,5]], ZZ), (x,y))

def test_Poly__new__():
    raises(GeneratorsError, "Poly(x+1, x, x)")

    raises(PolynomialError, "Poly(DMP([1,2], ZZ), x, y)")
    raises(PolynomialError, "Poly(GFP([1,2], 3, ZZ), x, y)")

    raises(PolynomialError, "Poly(DMP([1,2], ZZ), x, domain=ZZ)")
    raises(PolynomialError, "Poly(GFP([1,2], 3, ZZ), x, domain=ZZ)")

    raises(PolynomialError, "Poly(DMP([1,2], ZZ), x, modulus=3)")
    raises(PolynomialError, "Poly(GFP([1,2], 3, ZZ), x, modulus=3)")

    raises(OptionError, "Poly(x, x, symmetric=True)")

    raises(PolynomialError, "Poly(x+y, x, y, domain=ZZ[x])")
    raises(PolynomialError, "Poly(x+y, x, y, domain=ZZ[y])")

    raises(PolynomialError, "Poly(x+2, x, modulus=3, domain=QQ)")

    raises(OptionError, "Poly(x+2, x, domain=ZZ, gaussian=True)")
    raises(OptionError, "Poly(x+2, x, modulus=3, gaussian=True)")

    raises(OptionError, "Poly(x+2, x, domain=ZZ, extension=[sqrt(3)])")
    raises(OptionError, "Poly(x+2, x, modulus=3, extension=[sqrt(3)])")

    raises(OptionError, "Poly(x+2, x, domain=ZZ, extension=True)")
    raises(OptionError, "Poly(x+2, x, modulus=3, extension=True)")

    raises(OptionError, "Poly(x+2, x, domain=ZZ, greedy=True)")
    raises(OptionError, "Poly(x+2, x, domain=QQ, field=True)")

    raises(OptionError, "Poly(x+2, x, domain=ZZ, greedy=False)")
    raises(OptionError, "Poly(x+2, x, domain=QQ, field=False)")

    raises(NotImplementedError, "Poly(x+1, x, modulus=3, order='grlex')")
    raises(NotImplementedError, "Poly(x+1, x, order='grlex')")

    raises(GeneratorsNeeded, "Poly({1: 2, 0: 1})")
    raises(GeneratorsNeeded, "Poly([2, 1])")

    raises(GeneratorsNeeded, "Poly(1)")
    assert Poly(1, strict=False) == 1

    assert Poly(Poly(a*x + b*y, x, y), x) == Poly(a*x + b*y, x)

    assert Poly(3*x**2 + 2*x + 1, domain='ZZ').all_coeffs() == [3, 2, 1]
    assert Poly(3*x**2 + 2*x + 1, domain='QQ').all_coeffs() == [3, 2, 1]
    assert Poly(3*x**2 + 2*x + 1, domain='RR').all_coeffs() == [3.0, 2.0, 1.0]

    raises(CoercionFailed, "Poly(3*x**2/5 + 2*x/5 + 1, domain='ZZ')")
    assert Poly(3*x**2/5 + 2*x/5 + 1, domain='QQ').all_coeffs() == [S(3)/5, S(2)/5, 1]
    assert _eq(Poly(3*x**2/5 + 2*x/5 + 1, domain='RR').all_coeffs(),
            [0.6, 0.4, 1.0])

    assert Poly(3.0*x**2 + 2.0*x + 1, domain='ZZ').all_coeffs() == [3, 2, 1]
    assert Poly(3.0*x**2 + 2.0*x + 1, domain='QQ').all_coeffs() == [3, 2, 1]
    assert Poly(3.0*x**2 + 2.0*x + 1, domain='RR').all_coeffs() == [3.0, 2.0, 1.0]

    raises(CoercionFailed, "Poly(3.1*x**2 + 2.1*x + 1, domain='ZZ')")
    assert Poly(3.1*x**2 + 2.1*x + 1, domain='QQ').all_coeffs() == [S(31)/10, S(21)/10, 1]
    assert Poly(3.1*x**2 + 2.1*x + 1, domain='RR').all_coeffs() == [3.1, 2.1, 1.0]

    assert Poly({(2,1): 1, (1,2): 2, (1,1): 3}, x, y) == \
        Poly(x**2*y + 2*x*y**2 + 3*x*y, x, y)

    assert Poly(x**2 + 1, extension=I).get_domain() == QQ.algebraic_field(I)

    f = 3*x**5 - x**4 + x**3 - x** 2 + 65538

    assert Poly(f, x, modulus=65537, symmetric=True) == \
        Poly(3*x**5 - x**4 + x**3 - x** 2 + 1, x, modulus=65537, symmetric=True)
    assert Poly(f, x, modulus=65537, symmetric=False) == \
        Poly(3*x**5 + 65536*x**4 + x**3 + 65536*x** 2 + 1, x, modulus=65537, symmetric=False)

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

    assert Poly((x-p)*(x-q), sort='x > p > q').gens == (x, p, q)
    assert Poly((x-p)*(x-q), sort='p > x > q').gens == (p, x, q)
    assert Poly((x-p)*(x-q), sort='p > q > x').gens == (p, q, x)

    assert Poly((x-p)*(x-q), x, p, q, sort='p > q > x').gens == (x, p, q)

    assert Poly((x-p)*(x-q), wrt='x').gens == (x, p, q)
    assert Poly((x-p)*(x-q), wrt='p').gens == (p, x, q)
    assert Poly((x-p)*(x-q), wrt='q').gens == (q, x, p)

    assert Poly((x-p)*(x-q), wrt=x).gens == (x, p, q)
    assert Poly((x-p)*(x-q), wrt=p).gens == (p, x, q)
    assert Poly((x-p)*(x-q), wrt=q).gens == (q, x, p)

    assert Poly((x-p)*(x-q), x, p, q, wrt='p').gens == (x, p, q)

    assert Poly((x-p)*(x-q), wrt='p', sort='q > x').gens == (p, q, x)
    assert Poly((x-p)*(x-q), wrt='q', sort='p > x').gens == (q, p, x)

def test_Poly_unify():
    raises(UnificationFailed, "Poly(x).unify(y)")

    raises(UnificationFailed, "Poly(x, x, modulus=3).unify(Poly(x, x, modulus=5))")
    raises(UnificationFailed, "Poly(x, x, modulus=3).unify(Poly(y, y, modulus=3))")

    raises(UnificationFailed, "Poly(x, x, y).unify(Poly(x, x, modulus=3))")
    raises(UnificationFailed, "Poly(x, x, y).unify(Poly(x, x, modulus=3))")

    raises(UnificationFailed, "Poly(x, x, modulus=3).unify(Poly(x, x, y))")
    raises(UnificationFailed, "Poly(x, x, modulus=3).unify(Poly(x, x, y))")

    assert Poly(x+1, x).unify(Poly(x+2, x))[2:] == (DMP([1, 1], ZZ), DMP([1, 2], ZZ))
    assert Poly(x+1, x, domain='QQ').unify(Poly(x+2, x))[2:] == (DMP([1, 1], QQ), DMP([1, 2], QQ))
    assert Poly(x+1, x).unify(Poly(x+2, x, domain='QQ'))[2:] == (DMP([1, 1], QQ), DMP([1, 2], QQ))

    assert Poly(x+1, x).unify(Poly(x+2, x, y))[2:] == (DMP([[1], [1]], ZZ), DMP([[1], [2]], ZZ))
    assert Poly(x+1, x, domain='QQ').unify(Poly(x+2, x, y))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))
    assert Poly(x+1, x).unify(Poly(x+2, x, y, domain='QQ'))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))

    assert Poly(x+1, x, y).unify(Poly(x+2, x))[2:] == (DMP([[1], [1]], ZZ), DMP([[1], [2]], ZZ))
    assert Poly(x+1, x, y, domain='QQ').unify(Poly(x+2, x))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))
    assert Poly(x+1, x, y).unify(Poly(x+2, x, domain='QQ'))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))

    assert Poly(x+1, x, y).unify(Poly(x+2, x, y))[2:] == (DMP([[1], [1]], ZZ), DMP([[1], [2]], ZZ))
    assert Poly(x+1, x, y, domain='QQ').unify(Poly(x+2, x, y))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))
    assert Poly(x+1, x, y).unify(Poly(x+2, x, y, domain='QQ'))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))

    assert Poly(x+1, x).unify(Poly(x+2, y, x))[2:] == (DMP([[1, 1]], ZZ), DMP([[1, 2]], ZZ))
    assert Poly(x+1, x, domain='QQ').unify(Poly(x+2, y, x))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))
    assert Poly(x+1, x).unify(Poly(x+2, y, x, domain='QQ'))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))

    assert Poly(x+1, y, x).unify(Poly(x+2, x))[2:] == (DMP([[1, 1]], ZZ), DMP([[1, 2]], ZZ))
    assert Poly(x+1, y, x, domain='QQ').unify(Poly(x+2, x))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))
    assert Poly(x+1, y, x).unify(Poly(x+2, x, domain='QQ'))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))

    assert Poly(x+1, x, y).unify(Poly(x+2, y, x))[2:] == (DMP([[1], [1]], ZZ), DMP([[1], [2]], ZZ))
    assert Poly(x+1, x, y, domain='QQ').unify(Poly(x+2, y, x))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))
    assert Poly(x+1, x, y).unify(Poly(x+2, y, x, domain='QQ'))[2:] == (DMP([[1], [1]], QQ), DMP([[1], [2]], QQ))

    assert Poly(x+1, y, x).unify(Poly(x+2, x, y))[2:] == (DMP([[1, 1]], ZZ), DMP([[1, 2]], ZZ))
    assert Poly(x+1, y, x, domain='QQ').unify(Poly(x+2, x, y))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))
    assert Poly(x+1, y, x).unify(Poly(x+2, x, y, domain='QQ'))[2:] == (DMP([[1, 1]], QQ), DMP([[1, 2]], QQ))

    assert Poly(a*x, x, domain='ZZ[a]').unify(Poly(a*b*x, x, domain='ZZ(a,b)'))[2:] == \
        (DMP([DMF(([[1], []], [[1]]), ZZ), DMF(([[]], [[1]]), ZZ)], ZZ.frac_field(a,b)),
         DMP([DMF(([[1, 0], []], [[1]]), ZZ), DMF(([[]], [[1]]), ZZ)], ZZ.frac_field(a,b)))
    assert Poly(a*x, x, domain='ZZ(a)').unify(Poly(a*b*x, x, domain='ZZ(a,b)'))[2:] == \
        (DMP([DMF(([[1], []], [[1]]), ZZ), DMF(([[]], [[1]]), ZZ)], ZZ.frac_field(a,b)),
         DMP([DMF(([[1, 0], []], [[1]]), ZZ), DMF(([[]], [[1]]), ZZ)], ZZ.frac_field(a,b)))
    raises(CoercionFailed, "Poly(Poly(x**2 + x**2*z, y), domain='ZZ(x)')")
    raises(CoercionFailed, "Poly(Poly(x**2 + x**2*z, y, field=True), domain='ZZ(x)')")

    assert Poly(2*x+5, x).unify(Poly(x+2, x, modulus=3))[2:] == (GFP([2, 2], 3, ZZ), GFP([1, 2], 3, ZZ))
    assert Poly(x+2, x, modulus=3).unify(Poly(2*x+5, x))[2:] == (GFP([1, 2], 3, ZZ), GFP([2, 2], 3, ZZ))

    assert Poly(x+5, x, modulus=3).unify(Poly(x+7, x, modulus=3))[2:] == (GFP([1, 2], 3, ZZ), GFP([1, 1], 3, ZZ))

    assert Poly(x+5, x, modulus=3, symmetric=True).unify(Poly(x+7, x, modulus=3, symmetric=False))[2:] == \
        (GFP([1, 2], 3, ZZ, symmetric=True), GFP([1, 1], 3, ZZ, symmetric=True))

def test_Poly__eq__():
    assert (Poly(x, x) == Poly(x, x)) == True
    assert (Poly(x, x, domain=QQ) == Poly(x, x)) == True
    assert (Poly(x, x) == Poly(x, x, domain=QQ)) == True

    assert (Poly(x, x, domain=ZZ[a]) == Poly(x, x)) == True
    assert (Poly(x, x) == Poly(x, x, domain=ZZ[a])) == True

    assert (Poly(x*y, x, y) == Poly(x, x)) == False

def test_Poly__analyze_order():
    assert Poly._analyze_order({}) is None
    assert Poly._analyze_order({'order': 'lex'}) == monomial_lex_key

    raises(ValueError, "Poly._analyze_order({'order': 'foo'})")
    raises(ValueError, "Poly._analyze_order({'order': 1})")

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

    assert Poly._parse_domain('Q<I>') == QQ.algebraic_field(I)
    assert Poly._parse_domain('QQ<I>') == QQ.algebraic_field(I)

    assert Poly._parse_domain('Q<sqrt(2), I>') == QQ.algebraic_field(sqrt(2), I)
    assert Poly._parse_domain('QQ<sqrt(2), I>') == QQ.algebraic_field(sqrt(2), I)

def test_Poly_get_domain():
    assert Poly(2*x).get_domain() == ZZ

    assert Poly(2*x, domain='ZZ').get_domain() == ZZ
    assert Poly(2*x, domain='QQ').get_domain() == QQ

    assert Poly(x/2).get_domain() == QQ

    raises(CoercionFailed, "Poly(x/2, domain='ZZ')")
    assert Poly(x/2, domain='QQ').get_domain() == QQ

    assert Poly(0.2*x).get_domain() == RR

def test_Poly_set_domain():
    assert Poly(2*x + 1).set_domain(ZZ) == Poly(2*x + 1)
    assert Poly(2*x + 1).set_domain('ZZ') == Poly(2*x + 1)

    assert Poly(2*x + 1).set_domain(QQ) == Poly(2*x + 1, domain='QQ')
    assert Poly(2*x + 1).set_domain('QQ') == Poly(2*x + 1, domain='QQ')

    assert Poly(S(2)/10*x + S(1)/10).set_domain('RR') == Poly(0.2*x + 0.1)
    assert Poly(0.2*x + 0.1).set_domain('QQ') == Poly(S(2)/10*x + S(1)/10)

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

def test_Poly__analyze_extension():
    assert Poly._analyze_extension({}) is None
    assert Poly._analyze_extension({'extension': []}) is None
    assert Poly._analyze_extension({'extension': sqrt(2)}) == set([sqrt(2)])
    assert Poly._analyze_extension({'extension': [sqrt(2),sqrt(3)]}) == set([sqrt(2),sqrt(3)])

    assert Poly._analyze_extension({'extension': True}) is True
    assert Poly._analyze_extension({'extension': False}) is None

    assert Poly._analyze_extension({'extension': I}) == set([I])
    assert Poly._analyze_extension({'gaussian': True}) == set([I])

    raises(PolynomialError, "Poly._analyze_extension({'gaussian': True, 'extension': I})")
    raises(PolynomialError, "Poly._analyze_extension({'gaussian': True, 'split': True})")
    raises(PolynomialError, "Poly._analyze_extension({'extension': I, 'split': True})")

    raises(NotImplementedError, "Poly._analyze_extension({'split': True})")

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

    assert Poly(x+y+z+1).is_linear == True
    assert Poly(x*y*z+1).is_linear == False

    assert Poly(x*y).is_monomial == True
    assert Poly(x*y+1).is_monomial == False

    assert Poly(x*y+x).is_homogeneous == True
    assert Poly(x*y+x+1).is_homogeneous == False

    assert Poly(x).is_univariate == True
    assert Poly(x*y).is_univariate == False

    assert Poly(x*y).is_multivariate == True
    assert Poly(x).is_multivariate == False

def test_Poly_is_irreducible():
    assert Poly(7*x + 3, modulus=11).is_irreducible == True
    assert Poly(7*x**2 + 3*x + 1, modulus=11).is_irreducible == False

def test_Poly_replace():
    assert Poly(x+1).replace(x) == Poly(x+1)
    assert Poly(x+1).replace(y) == Poly(y+1)

    raises(PolynomialError, "Poly(x+y).replace(z)")

    assert Poly(x+1).replace(x, x) == Poly(x+1)
    assert Poly(x+1).replace(x, y) == Poly(y+1)

    assert Poly(x+y).replace(x, x) == Poly(x+y)
    assert Poly(x+y).replace(x, z) == Poly(z+y, z, y)

    assert Poly(x+y).replace(y, y) == Poly(x+y)
    assert Poly(x+y).replace(y, z) == Poly(x+z, x, z)

    raises(PolynomialError, "Poly(x+y).replace(x, y)")
    raises(PolynomialError, "Poly(x+y).replace(z, t)")

    assert Poly(x+y, x).replace(x, z) == Poly(z+y, z)
    assert Poly(x+y, y).replace(y, z) == Poly(x+z, z)

    raises(PolynomialError, "Poly(x+y, x).replace(x, y)")
    raises(PolynomialError, "Poly(x+y, y).replace(y, x)")

def test_Poly_reorder():
    raises(PolynomialError, "Poly(x+y).reorder(x, z)")

    assert Poly(x + y, x, y).reorder(x, y) == Poly(x + y, x, y)
    assert Poly(x + y, x, y).reorder(y, x) == Poly(x + y, y, x)

    assert Poly(x + y, y, x).reorder(x, y) == Poly(x + y, x, y)
    assert Poly(x + y, y, x).reorder(y, x) == Poly(x + y, y, x)

    assert Poly(x + y, x, y).reorder(wrt=x) == Poly(x + y, x, y)
    assert Poly(x + y, x, y).reorder(wrt=y) == Poly(x + y, y, x)

def test_Poly_to_ring():
    assert Poly(2*x+1, domain='ZZ').to_ring() == Poly(2*x+1, domain='ZZ')
    assert Poly(2*x+1, domain='QQ').to_ring() == Poly(2*x+1, domain='ZZ')

    raises(CoercionFailed, "Poly(x/2+1).to_ring()")
    raises(OperationNotSupported, "Poly(2*x+1, modulus=3).to_ring()")

def test_Poly_to_field():
    assert Poly(2*x+1, domain='ZZ').to_field() == Poly(2*x+1, domain='QQ')
    assert Poly(2*x+1, domain='QQ').to_field() == Poly(2*x+1, domain='QQ')

    assert Poly(x/2+1, domain='QQ').to_field() == Poly(x/2+1, domain='QQ')
    assert Poly(2*x+1, modulus=3).to_field() == Poly(2*x+1, modulus=3)

    raises(DomainError, "Poly(2.0*x + 1.0).to_field()")

def test_Poly_to_exact():
    assert Poly(2*x).to_exact() == Poly(2*x)
    assert Poly(x/2).to_exact() == Poly(x/2)

    assert Poly(0.1*x).to_exact() == Poly(x/10)

    raises(OperationNotSupported, "Poly(x, modulus=2).to_exact()")

def test_Poly_coeffs():
    assert Poly(0, x).coeffs() == [0]
    assert Poly(1, x).coeffs() == [1]

    assert Poly(2*x+1, x).coeffs() == [2,1]

    assert Poly(7*x**2+2*x+1, x).coeffs() == [7,2,1]
    assert Poly(7*x**4+2*x+1, x).coeffs() == [7,2,1]

    assert Poly(x*y**7 + 2*x**2*y**3).coeffs('lex') == [2, 1]
    assert Poly(x*y**7 + 2*x**2*y**3).coeffs('grlex') == [1, 2]

def test_Poly_monoms():
    assert Poly(0, x).monoms() == [(0,)]
    assert Poly(1, x).monoms() == [(0,)]

    assert Poly(2*x+1, x).monoms() == [(1,),(0,)]

    assert Poly(7*x**2+2*x+1, x).monoms() == [(2,),(1,),(0,)]
    assert Poly(7*x**4+2*x+1, x).monoms() == [(4,),(1,),(0,)]

    assert Poly(x*y**7 + 2*x**2*y**3).monoms('lex') == [(2, 3), (1, 7)]
    assert Poly(x*y**7 + 2*x**2*y**3).monoms('grlex') == [(1, 7), (2, 3)]

def test_Poly_terms():
    assert Poly(0, x).terms() == [((0,), 0)]
    assert Poly(1, x).terms() == [((0,), 1)]

    assert Poly(2*x+1, x).terms() == [((1,), 2),((0,), 1)]

    assert Poly(7*x**2+2*x+1, x).terms() == [((2,), 7),((1,), 2),((0,), 1)]
    assert Poly(7*x**4+2*x+1, x).terms() == [((4,), 7),((1,), 2),((0,), 1)]

    assert Poly(x*y**7 + 2*x**2*y**3).terms('lex') == [((2, 3), 2), ((1, 7), 1)]
    assert Poly(x*y**7 + 2*x**2*y**3).terms('grlex') == [((1, 7), 1), ((2, 3), 2)]

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

def test_Poly_lift():
    assert Poly(x**4 - I*x + 17*I, x, gaussian=True).lift() == \
        Poly(x**16 + 2*x**10 + 578*x**8 + x**4 - 578*x**2 + 83521, x, domain='QQ')

def test_Poly_deflate():
    assert Poly(0, x).deflate() == ((1,), Poly(0, x))
    assert Poly(1, x).deflate() == ((1,), Poly(1, x))
    assert Poly(x, x).deflate() == ((1,), Poly(x, x))

    assert Poly(x**2, x).deflate() == ((2,), Poly(x, x))
    assert Poly(x**17, x).deflate() == ((17,), Poly(x, x))

    assert Poly(x**2*y*z**11+x**4*z**11).deflate() == ((2,1,11), Poly(x*y*z+x**2*z))

def test_Poly_exclude():
    assert Poly(x, x, y).exclude() == Poly(x, x)
    assert Poly(x*y, x, y).exclude() == Poly(x*y, x, y)
    assert Poly(1, x, y).exclude() == Poly(1, x, y)

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

    assert degree(1, x) == 0
    assert degree(x, x) == 1

    assert degree(x*y**2, gen=x) == 1
    assert degree(x*y**2, gen=y) == 2

    assert degree(x*y**2, x, y) == 1
    assert degree(x*y**2, y, x) == 2

    raises(GeneratorsNeeded, "degree(1)")

def test_Poly_degree_list():
    assert Poly(0, x).degree_list() == (-1,)
    assert Poly(0, x, y).degree_list() == (-1,-1)
    assert Poly(0, x, y, z).degree_list() == (-1,-1,-1)

    assert Poly(1, x).degree_list() == (0,)
    assert Poly(1, x, y).degree_list() == (0,0)
    assert Poly(1, x, y, z).degree_list() == (0,0,0)

    assert Poly(x**2*y+x**3*z**2+1).degree_list() == (3,1,2)

    assert degree_list(1, x) == (0,)
    assert degree_list(x, x) == (1,)

    assert degree_list(x*y**2) == (1,2)

    raises(GeneratorsNeeded, "degree_list(1)")

def test_Poly_total_degree():
    assert Poly(x**2*y+x**3*z**2+1).total_degree() == 6

def test_Poly_LC():
    assert Poly(0, x).LC() == 0
    assert Poly(1, x).LC() == 1
    assert Poly(2*x**2+x, x).LC() == 2

    assert Poly(x*y**7 + 2*x**2*y**3).LC('lex') == 2
    assert Poly(x*y**7 + 2*x**2*y**3).LC('grlex') == 1

    assert LC(x*y**7 + 2*x**2*y**3, order='lex') == 2
    assert LC(x*y**7 + 2*x**2*y**3, order='grlex') == 1

def test_Poly_TC():
    assert Poly(0, x).TC() == 0
    assert Poly(1, x).TC() == 1
    assert Poly(2*x**2+x, x).TC() == 0

def test_Poly_EC():
    assert Poly(0, x).EC() == 0
    assert Poly(1, x).EC() == 1
    assert Poly(2*x**2+x, x).EC() == 1

    assert Poly(x*y**7 + 2*x**2*y**3).EC('lex') == 1
    assert Poly(x*y**7 + 2*x**2*y**3).EC('grlex') == 2

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

    assert Poly(x*y**7 + 2*x**2*y**3).LM('lex') == (2, 3)
    assert Poly(x*y**7 + 2*x**2*y**3).LM('grlex') == (1, 7)

    assert LM(x*y**7 + 2*x**2*y**3, order='lex') == x**2*y**3
    assert LM(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

def test_Poly_EM():
    assert Poly(0, x).EM() == (0,)
    assert Poly(1, x).EM() == (0,)
    assert Poly(2*x**2+x, x).EM() == (1,)

    assert Poly(x*y**7 + 2*x**2*y**3).EM('lex') == (1, 7)
    assert Poly(x*y**7 + 2*x**2*y**3).EM('grlex') == (2, 3)

def test_Poly_LT():
    assert Poly(0, x).LT() == ((0,), 0)
    assert Poly(1, x).LT() == ((0,), 1)
    assert Poly(2*x**2+x, x).LT() == ((2,), 2)

    assert Poly(x*y**7 + 2*x**2*y**3).LT('lex') == ((2, 3), 2)
    assert Poly(x*y**7 + 2*x**2*y**3).LT('grlex') == ((1, 7), 1)

    assert LT(x*y**7 + 2*x**2*y**3, order='lex') == 2*x**2*y**3
    assert LT(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

def test_Poly_ET():
    assert Poly(0, x).ET() == ((0,), 0)
    assert Poly(1, x).ET() == ((0,), 1)
    assert Poly(2*x**2+x, x).ET() == ((1,), 1)

    assert Poly(x*y**7 + 2*x**2*y**3).ET('lex') == ((1, 7), 1)
    assert Poly(x*y**7 + 2*x**2*y**3).ET('grlex') == ((2, 3), 2)

def test_Poly_max_norm():
    assert Poly(-1, x).max_norm() == 1
    assert Poly( 0, x).max_norm() == 0
    assert Poly( 1, x).max_norm() == 1

def test_Poly_l1_norm():
    assert Poly(-1, x).l1_norm() == 1
    assert Poly( 0, x).l1_norm() == 0
    assert Poly( 1, x).l1_norm() == 1

def test_Poly_clear_denoms():
    coeff, poly = Poly(x + 2).clear_denoms()

    assert coeff == 1 and poly == Poly(x + 2, domain='ZZ') and poly.get_domain() == ZZ

    coeff, poly = Poly(x/2 + 1).clear_denoms()

    assert coeff == 2 and poly == Poly(x + 2, domain='QQ') and poly.get_domain() == QQ

    coeff, poly = Poly(x/2 + 1).clear_denoms(convert=True)

    assert coeff == 2 and poly == Poly(x + 2, domain='ZZ') and poly.get_domain() == ZZ

    coeff, poly = Poly(1/x*y, y, domain='ZZ(x)').clear_denoms()

    assert coeff == x and poly == Poly(y, y) and poly.get_domain() == ZZ.frac_field(x)

    coeff, poly = Poly(sin(x)/x*y, y, domain='EX').clear_denoms()
    assert coeff == x and poly == Poly(sin(x)*y, y) and poly.get_domain() == EX


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

def test_poly_cancel():
    a = Poly(y, y, domain='ZZ(x)')
    b = Poly(1, y, domain='ZZ[x]')
    assert a.cancel(b) == (1, Poly(y, y, domain='ZZ(x)'), Poly(1, y, domain='ZZ(x)'))
    assert a.cancel(b, include=True) == (Poly(y, y, domain='ZZ(x)'), Poly(1, y, domain='ZZ(x)'))
    a = Poly(5*x*y + x, y, domain='ZZ(x)')
    b = Poly(2*x**2*y, y, domain='ZZ(x)')
    assert a.cancel(b, include=True) == (Poly(5*y + 1, y, domain='ZZ(x)'),
        Poly(2*x*y, y, domain='ZZ(x)'))


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

    F, G, S, T, H = [ Poly(u, x, domain='QQ') for u in (f, g, s, t, h) ]

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

    raises(DomainError, "half_gcdex(x + 1, 2*x + 1, auto=False)")
    raises(DomainError, "gcdex(x + 1, 2*x + 1, auto=False)")
    raises(DomainError, "invert(x + 1, 2*x + 1, auto=False)")

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
    assert resultant(f, g, includePRS=True) == (h, [f, g, 2*x - 2])

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

    f, g = x**2 - 1, x - 1.0
    h, s, t = g, x + 1.0, 1.0

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == f

    f, g = x**2 - 1.0, x - 1
    h, s, t = g, x + 1.0, 1.0

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == f

    assert cofactors(8, 6) == (2, 4, 3)
    assert gcd(8, 6) == 2
    assert lcm(8, 6) == 24

    f, g = x**2 - 3*x - 4, x**3 - 4*x**2 + x - 4
    l = x**4 - 3*x**3 - 3*x**2 - 3*x - 4
    h, s, t = x - 4, x + 1, x**2 + 1

    assert cofactors(f, g, modulus=11) == (h, s, t)
    assert gcd(f, g, modulus=11) == h
    assert lcm(f, g, modulus=11) == l

    f, g = x**2 + 8*x + 7, x**3 + 7*x**2 + x + 7
    l = x**4 + 8*x**3 + 8*x**2 + 8*x + 7
    h, s, t = x + 7, x + 1, x**2 + 1

    assert cofactors(f, g, modulus=11, symmetric=False) == (h, s, t)
    assert gcd(f, g, modulus=11, symmetric=False) == h
    assert lcm(f, g, modulus=11, symmetric=False) == l

def test_terms_gcd():
    assert terms_gcd(1) == 1
    assert terms_gcd(1, x) == 1

    assert terms_gcd(x - 1) == x - 1
    assert terms_gcd(-x - 1) == -x - 1

    assert terms_gcd(2*x + 3) != Mul(1, 2*x + 3, evaluate=False)
    assert terms_gcd(6*x + 4) == Mul(2, 3*x + 2, evaluate=False)

    assert terms_gcd(x**3*y + x*y**3) == x*y*(x**2 + y**2)
    assert terms_gcd(2*x**3*y + 2*x*y**3) == 2*x*y*(x**2 + y**2)
    assert terms_gcd(x**3*y/2 + x*y**3/2) == x*y/2*(x**2 + y**2)

    assert terms_gcd(x**3*y + 2*x*y**3) == x*y*(x**2 + 2*y**2)
    assert terms_gcd(2*x**3*y + 4*x*y**3) == 2*x*y*(x**2 + 2*y**2)
    assert terms_gcd(2*x**3*y/3 + 4*x*y**3/5) == 2*x*y/15*(5*x**2 + 6*y**2)

    assert terms_gcd(2.0*x**3*y + 4.1*x*y**3) == x*y*(2.0*x**2 + 4.1*y**2)

def test_trunc():
    f, g = x**5 + 2*x**4 + 3*x**3 + 4*x**2 + 5*x + 6, x**5 - x**4 + x**2 - x
    F, G = Poly(f), Poly(g)

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(f, 3, (x,)) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f, g = 6*x**5 + 5*x**4 + 4*x**3 + 3*x**2 + 2*x + 1, -x**4 + x**3 - x + 1
    F, G = Poly(f), Poly(g)

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(f, 3, (x,)) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f = Poly(x**2 + 2*x + 3, modulus=5)

    assert f.trunc(2) == Poly(x**2 + 1, modulus=2)

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

    assert monic(2*x**2 + 6*x + 4, auto=False) == x**2 + 3*x + 2
    raises(ExactQuotientFailed, "monic(2*x + 6*x + 1, auto=False)")

    assert monic(2.0*x**2 + 6.0*x + 4.0) == x**2 + 3.0*x + 2.0
    assert monic(2*x**2 + 3*x + 4, modulus=5) == x**2 - x + 2

def test_content():
    f, F = 4*x + 2, Poly(4*x + 2)

    F.content() == 2
    content(f) == 2

    raises(GeneratorsNeeded, "content(4)")
    raises(OperationNotSupported, "Poly(2*x, modulus=3).content()")

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
    raises(OperationNotSupported, "Poly(2*x, modulus=3).primitive()")

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

    assert compose(x**2 - y**2, x - y, x, y) ==  x**2 - 2*x*y
    assert compose(x**2 - y**2, x - y, y, x) == -y**2 + 2*x*y

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
    raises(DomainError, "sturm(f, auto=False)")

    f = Poly(S(1024)/(15625*pi**8)*x**5   \
           - S(4096)/(625*pi**8)*x**4     \
           + S(32)/(15625*pi**4)*x**3     \
           - S(128)/(625*pi**4)*x**2      \
           + S(1)/62500*x                 \
           - S(1)/625, x, domain='ZZ(pi)')

    assert sturm(f) == \
        [Poly(x**3 - 100*x**2 + pi**4/64*x - 25*pi**4/16, x, domain='ZZ(pi)'),
         Poly(3*x**2 - 200*x + pi**4/64, x, domain='ZZ(pi)'),
         Poly((S(20000)/9 - pi**4/96)*x + 25*pi**4/18, x, domain='ZZ(pi)'),
         Poly((-3686400000000*pi**4 - 11520000*pi**8 - 9*pi**12)/(26214400000000 - 245760000*pi**4 + 576*pi**8), x, domain='ZZ(pi)')]

def test_sqf_norm():
    assert sqf_norm(x**2-2, extension=sqrt(3)) == \
        (1, x**2 - 2*sqrt(3)*x + 1, x**4 - 10*x**2 + 1)
    assert sqf_norm(x**2-3, extension=sqrt(2)) == \
        (1, x**2 - 2*sqrt(2)*x - 1, x**4 - 10*x**2 + 1)

    assert Poly(x**2-2, extension=sqrt(3)).sqf_norm() == \
        (1, Poly(x**2 - 2*sqrt(3)*x + 1, x, extension=sqrt(3)),
            Poly(x**4 - 10*x**2 + 1, x, domain='QQ'))

    assert Poly(x**2-3, extension=sqrt(2)).sqf_norm() == \
        (1, Poly(x**2 - 2*sqrt(2)*x - 1, x, extension=sqrt(2)),
            Poly(x**4 - 10*x**2 + 1, x, domain='QQ'))

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

    assert F.sqf_list_include() == [(G, 1), (H, 2)]
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

    assert sqf(x - 1) == x - 1
    assert sqf(-x - 1) == -x - 1

    assert sqf(x - 1) != Mul(1, x - 1, evaluate=False)
    assert sqf(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert sqf((6*x - 10)/(3*x - 6), frac=True) == S(2)/3*((3*x - 5)/(x - 2))

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

    assert F.factor_list_include() == [(U, 1), (V, 2), (W, 1)]
    assert factor_list(f, include=True) == [(u, 1), (v, 2), (w, 1)]

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

    f = x**4 + 1

    assert factor(f) == f
    assert factor(f, extension=I) == (x**2 - I)*(x**2 + I)
    assert factor(f, gaussian=True) == (x**2 - I)*(x**2 + I)
    assert factor(f, extension=sqrt(2)) == (x**2 + sqrt(2)*x + 1)*(x**2 - sqrt(2)*x + 1)

    f = x**2 + 2*sqrt(2)*x + 2

    assert factor(f, extension=sqrt(2)) == (x + sqrt(2))**2
    assert factor(f**3, extension=sqrt(2)) == (x + sqrt(2))**6

    assert factor(x**2 - 2*y**2, extension=sqrt(2)) == \
        (x + sqrt(2)*y)*(x - sqrt(2)*y)
    assert factor(2*x**2 - 4*y**2, extension=sqrt(2)) == \
        2*((x + sqrt(2)*y)*(x - sqrt(2)*y))

    assert factor(x - 1) == x - 1
    assert factor(-x - 1) == -x - 1

    assert factor(x - 1) != Mul(1, x - 1, evaluate=False)
    assert factor(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert factor((6*x - 10)/(3*x - 6), frac=True) == S(2)/3*((3*x - 5)/(x - 2))

    assert factor(x**11 + x + 1, modulus=65537, symmetric=True) == \
        (x**2 + x + 1)*(x**9 - x**8 + x**6 - x**5 + x**3 - x** 2 + 1)
    assert factor(x**11 + x + 1, modulus=65537, symmetric=False) == \
        (x**2 + x + 1)*(x**9 + 65536*x**8 + x**6 + 65536*x**5 + x**3 + 65536*x** 2 + 1)

    assert factor(x/pi + x*sin(x)/pi) == x*(sin(x) + 1)/pi

    f = y/(pi**2 + 2*pi + 1) + y*sin(x)/(pi**2 + 2*pi + 1)

    assert factor(f) == y*(sin(x) + 1)/(pi**2 + 2*pi + 1)
    assert factor(f, frac=True) == y*(sin(x) + 1)/(pi + 1)**2

def test_intervals():
    assert intervals(0) == []
    assert intervals(1) == []

    assert intervals(x, sqf=True) == [(0, 0)]
    assert intervals(x) == [((0, 0), 1)]

    assert intervals(x**128) == [((0, 0), 128)]
    assert intervals([x**2, x**4]) == [((0, 0), {0: 2, 1: 4})]

    f = Poly((2*x/5 - S(17)/3)*(4*x + S(1)/257))

    assert f.intervals(sqf=True) == [(-1, 0), (14, 15)]
    assert f.intervals() == [((-1, 0), 1), ((14, 15), 1)]

    assert f.intervals(fast=True, sqf=True) == [(-1, 0), (14, 15)]
    assert f.intervals(fast=True) == [((-1, 0), 1), ((14, 15), 1)]

    assert f.intervals(eps=S(1)/10) == f.intervals(eps=0.1) == \
        [((-S(1)/258, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert f.intervals(eps=S(1)/100) == f.intervals(eps=0.01) == \
        [((-S(1)/258, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert f.intervals(eps=S(1)/1000) == f.intervals(eps=0.001) == \
        [((-S(1)/1005, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert f.intervals(eps=S(1)/10000) == f.intervals(eps=0.0001) == \
        [((-S(1)/1028, -S(1)/1028), 1), ((S(85)/6, S(85)/6), 1)]

    f = (2*x/5 - S(17)/3)*(4*x + S(1)/257)

    assert intervals(f, sqf=True) == [(-1, 0), (14, 15)]
    assert intervals(f) == [((-1, 0), 1), ((14, 15), 1)]

    assert intervals(f, eps=S(1)/10) == intervals(f, eps=0.1) == \
        [((-S(1)/258, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert intervals(f, eps=S(1)/100) == intervals(f, eps=0.01) == \
        [((-S(1)/258, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert intervals(f, eps=S(1)/1000) == intervals(f, eps=0.001) == \
        [((-S(1)/1005, 0), 1), ((S(85)/6, S(85)/6), 1)]
    assert intervals(f, eps=S(1)/10000) == intervals(f, eps=0.0001) == \
        [((-S(1)/1028, -S(1)/1028), 1), ((S(85)/6, S(85)/6), 1)]

    f = Poly((x**2 - 2)*(x**2-3)**7*(x+1)*(7*x+3)**3)

    assert f.intervals() == \
        [((-2, -S(3)/2), 7), ((-S(3)/2, -1), 1),
         ((-1, -1), 1), ((-1, 0), 3),
         ((1, S(3)/2), 1), ((S(3)/2, 2), 7)]

    assert intervals([x**5 - 200, x**5 - 201]) == \
        [((S(75)/26, S(101)/35), {0: 1}), ((S(283)/98, S(26)/9), {1: 1})]

    assert intervals([x**5 - 200, x**5 - 201], fast=True) == \
        [((S(75)/26, S(101)/35), {0: 1}), ((S(283)/98, S(26)/9), {1: 1})]

    assert intervals([x**2 - 200, x**2 - 201]) == \
        [((-S(71)/5, -S(85)/6), {1: 1}), ((-S(85)/6, -14), {0: 1}), ((14, S(85)/6), {0: 1}), ((S(85)/6, S(71)/5), {1: 1})]

    assert intervals([x+1, x+2, x-1, x+1, 1, x-1, x-1, (x-2)**2]) == \
        [((-2, -2), {1: 1}), ((-1, -1), {0: 1, 3: 1}), ((1, 1), {2: 1, 5: 1, 6: 1}), ((2, 2), {7: 2})]

    f, g, h = x**2 - 2, x**4 - 4*x**2 + 4, x - 1

    assert intervals(f, inf=S(7)/4, sqf=True) == []
    assert intervals(f, inf=S(7)/5, sqf=True) == [(S(7)/5, S(3)/2)]
    assert intervals(f, sup=S(7)/4, sqf=True) == [(-2, -1), (1, S(3)/2)]
    assert intervals(f, sup=S(7)/5, sqf=True) == [(-2, -1)]

    assert intervals(g, inf=S(7)/4) == []
    assert intervals(g, inf=S(7)/5) == [((S(7)/5, S(3)/2), 2)]
    assert intervals(g, sup=S(7)/4) == [((-2, -1), 2), ((1, S(3)/2), 2)]
    assert intervals(g, sup=S(7)/5) == [((-2, -1), 2)]

    assert intervals([g, h], inf=S(7)/4) == []
    assert intervals([g, h], inf=S(7)/5) == [((S(7)/5, S(3)/2), {0: 2})]
    assert intervals([g, h], sup=S(7)/4) == [((-2, -1), {0: 2}), ((1, 1), {1: 1}), ((1, S(3)/2), {0: 2})]
    assert intervals([g, h], sup=S(7)/5) == [((-2, -1), {0: 2}), ((1, 1), {1: 1})]

    assert intervals([x+2, x**2 - 2]) == \
        [((-2, -2), {0: 1}), ((-2, -1), {1: 1}), ((1, 2), {1: 1})]
    assert intervals([x+2, x**2 - 2], strict=True) == \
        [((-2, -2), {0: 1}), ((-S(3)/2, -1), {1: 1}), ((1, 2), {1: 1})]

    f = 7*z**4 - 19*z**3 + 20*z**2 + 17*z + 20

    assert intervals(f) == []

    roots = sorted(nroots(f), key=lambda r: (re(r), -im(r)))

    real_part, complex_part = intervals(f, all=True, sqf=True)

    assert real_part == []
    assert all([ re(a) < re(r) < re(b) and im(a) < im(r) < im(b) for (a, b), r in zip(complex_part, roots) ])

    assert complex_part == [(-S(40)/7, 40*I/7), (-S(40)/7 - 40*I/7, 0),
                            (0, S(40)/7 + 40*I/7), (-40*I/7, S(40)/7)]

    real_part, complex_part = intervals(f, all=True, sqf=True, eps=S(1)/10)

    assert real_part == []
    assert all([ re(a) < re(r) < re(b) and im(a) < im(r) < im(b) for (a, b), r in zip(complex_part, roots) ])

def test_refine_root():
    f = Poly(x**2 - 2)

    assert f.refine_root(1, 2, steps=0) == (1, 2)
    assert f.refine_root(-2, -1, steps=0) == (-2, -1)

    assert f.refine_root(1, 2, steps=None) == (1, S(3)/2)
    assert f.refine_root(-2, -1, steps=None) == (-S(3)/2, -1)

    assert f.refine_root(1, 2, steps=1) == (1, S(3)/2)
    assert f.refine_root(-2, -1, steps=1) == (-S(3)/2, -1)

    assert f.refine_root(1, 2, steps=1, fast=True) == (1, S(3)/2)
    assert f.refine_root(-2, -1, steps=1, fast=True) == (-S(3)/2, -1)

    assert f.refine_root(1, 2, eps=S(1)/100) == (S(24)/17, S(17)/12)
    assert f.refine_root(1, 2, eps=1e-2) == (S(24)/17, S(17)/12)

    raises(PolynomialError, "(f**2).refine_root(1, 2, check_sqf=True)")

    raises(RefinementFailed, "(f**2).refine_root(1, 2)")
    raises(RefinementFailed, "(f**2).refine_root(2, 3)")

    f = x**2 - 2

    assert refine_root(f, 1, 2, steps=1) == (1, S(3)/2)
    assert refine_root(f, -2, -1, steps=1) == (-S(3)/2, -1)

    assert refine_root(f, 1, 2, steps=1, fast=True) == (1, S(3)/2)
    assert refine_root(f, -2, -1, steps=1, fast=True) == (-S(3)/2, -1)

    assert refine_root(f, 1, 2, eps=S(1)/100) == (S(24)/17, S(17)/12)
    assert refine_root(f, 1, 2, eps=1e-2) == (S(24)/17, S(17)/12)

    raises(PolynomialError, "refine_root(1, 7, 8, eps=S(1)/100)")

def test_count_roots():
    assert count_roots(x**2 - 2) == 2

    assert count_roots(x**2 - 2, inf=-oo) == 2
    assert count_roots(x**2 - 2, sup=+oo) == 2
    assert count_roots(x**2 - 2, inf=-oo, sup=+oo) == 2

    assert count_roots(x**2 - 2, inf=-2) == 2
    assert count_roots(x**2 - 2, inf=-1) == 1

    assert count_roots(x**2 - 2, sup=1) == 1
    assert count_roots(x**2 - 2, sup=2) == 2

    assert count_roots(x**2 - 2, inf=-1, sup=1) == 0
    assert count_roots(x**2 - 2, inf=-2, sup=2) == 2

    assert count_roots(x**2 - 2, inf=-1, sup=1) == 0
    assert count_roots(x**2 - 2, inf=-2, sup=2) == 2

    assert count_roots(x**2 + 2) == 0
    assert count_roots(x**2 + 2, inf=-2*I) == 2
    assert count_roots(x**2 + 2, sup=+2*I) == 2
    assert count_roots(x**2 + 2, inf=-2*I, sup=+2*I) == 2

    assert count_roots(x**2 + 2, inf=0) == 0
    assert count_roots(x**2 + 2, sup=0) == 0

    assert count_roots(x**2 + 2, inf=-I) == 1
    assert count_roots(x**2 + 2, sup=+I) == 1

    assert count_roots(x**2 + 2, inf=+I/2, sup=+I) == 0
    assert count_roots(x**2 + 2, inf=-I, sup=-I/2) == 0

    raises(PolynomialError, "count_roots(1)")

def test_real_roots():
    assert real_roots(x) == [0]
    assert real_roots(x, multiple=False) == [(0, 1)]

    assert real_roots(x**3) == [0, 0, 0]
    assert real_roots(x**3, multiple=False) == [(0, 3)]

    assert real_roots(x*(x**3 + x + 3)) == [RootOf(x**3 + x + 3, 0), 0]
    assert real_roots(x*(x**3 + x + 3), multiple=False) == [(RootOf(x**3 + x + 3, 0), 1), (0, 1)]

    assert real_roots(x**3*(x**3 + x + 3)) == [RootOf(x**3 + x + 3, 0), 0, 0, 0]
    assert real_roots(x**3*(x**3 + x + 3), multiple=False) == [(RootOf(x**3 + x + 3, 0), 1), (0, 3)]

def test_nroots():
    assert Poly(0, x).nroots() == []
    assert Poly(1, x).nroots() == []

    assert Poly(x**2 - 1, x).nroots() == [-1.0, 1.0]
    assert Poly(x**2 + 1, x).nroots() == [-I, I]

    roots, error = Poly(x**2 - 1, x).nroots(error=True)
    assert roots == [-1.0, 1.0] and error < 1e25;

    roots, error = Poly(x**2 + 1, x).nroots(error=True)
    assert roots == [-I, I] and error < 1e25;

    roots, error = Poly(x**2/3 - S(1)/3, x).nroots(error=True)
    assert roots == [-1.0, 1.0] and error < 1e25;

    roots, error = Poly(x**2/3 + S(1)/3, x).nroots(error=True)
    assert roots == [-I, I] and error < 1e25;

    assert Poly(x**2 + 2*I, x).nroots() == [-1.0 + I, 1.0 - I]
    assert Poly(x**2 + 2*I, x, extension=I).nroots() == [-1.0 + I, 1.0 - I]

    assert Poly(0.2*x + 0.1).nroots() == [-0.5]

    raises(DomainError, "Poly(x+y, x).nroots()")
    raises(PolynomialError, "Poly(x+y).nroots()")

    assert nroots(x**2 - 1) == [-1.0, 1.0]

    roots, error = nroots(x**2 - 1, error=True)
    assert roots == [-1.0, 1.0] and error < 1e25;

    assert nroots(x + I) == [-I]
    assert nroots(x + 2*I) == [-2*I]

    raises(PolynomialError, "nroots(0)")

def test_cancel():
    assert cancel(0) == 0
    assert cancel(7) == 7
    assert cancel(x) == x

    assert cancel(oo) == oo

    assert cancel((2, 3)) == (1, 2, 3)

    assert cancel((1, 0), x) == (1, 1, 0)
    assert cancel((0, 1), x) == (1, 0, 1)

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

    f = Poly(x**2 - a**2, x)
    g = Poly(x - a, x)

    F = Poly(x + a, x)
    G = Poly(1, x)

    assert cancel((f, g)) == (1, F, G)

    f = x**3 + (sqrt(2) - 2)*x**2 - (2*sqrt(2) + 3)*x - 3*sqrt(2)
    g = x**2 - 2

    assert cancel((f, g), extension=True) == (1, x**2 - 2*x - 3, x - sqrt(2))

def test_reduced():
    raises(PolynomialError, "reduced(x, [x], x, modulus=3)")

    f = 2*x**4 + y**2 - x**2 + y**3
    G = [x**3 - x, y**3 - y]

    Q = [2*x, 1]
    r = x**2 + y**2 + y

    assert reduced(f, G) == (Q, r)
    assert reduced(f, G, x, y) == (Q, r)

    Q = [Poly(2*x, x, y), Poly(1, x, y)]
    r = Poly(x**2 + y**2 + y, x, y)

    assert reduced(f, G, polys=True) == (Q, r)
    assert reduced(f, G, x, y, polys=True) == (Q, r)

def test_groebner():
    raises(PolynomialError, "groebner([x], x, modulus=3)")

    assert groebner([], x, y, z) == []

    assert groebner([x**2 + 1, y**4*x + x**3],
        x, y, order='lex') == [1 + x**2, -1 + y**4]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3],
        x, y, z, order='grevlex') == [-1 + y**4, z**3, 1 + x**2]

    assert groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex', polys=True) == \
        [Poly(1 + x**2, x, y), Poly(-1 + y**4, x, y)]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex', polys=True) == \
        [Poly(-1 + y**4, x, y, z), Poly(z**3, x, y, z), Poly(1 + x**2, x, y, z)]

def test_symmetrize():
    assert symmetrize(0) == (0, 0)
    assert symmetrize(1) == (1, 0)

    assert symmetrize(0, x, y, z) == (0, 0)
    assert symmetrize(1, x, y, z) == (1, 0)

    assert symmetrize(0, formal=True) == (0, 0, {})
    assert symmetrize(1, formal=True) == (1, 0, {})

    s1 = x + y + z
    s2 = x*y + x*z + y*z
    s3 = x*y*z

    assert symmetrize(x) == (x, 0)
    assert symmetrize(x + 1) == (x + 1, 0)

    assert symmetrize(x, x, y) == (x + y, -y)
    assert symmetrize(x + 1, x, y) == (x + y + 1, -y)

    assert symmetrize(x, x, y, z) == (s1, -y - z)
    assert symmetrize(x + 1, x, y, z) == (s1 + 1, -y - z)

    assert symmetrize(x**2, x, y, z) == (s1**2 - 2*s2, -y**2 - z**2)

    assert symmetrize(x**2 + y**2) == (-2*x*y + (x + y)**2, 0)
    assert symmetrize(x**2 - y**2) == (-2*x*y + (x + y)**2, -2*y**2)

    assert symmetrize(x**3 + y**2 + a*x**2 + b*y**3, x, y) == \
        (-3*x*y*(x + y) - 2*a*x*y + a*(x + y)**2 + (x + y)**3, y**2*(1 - a) - y**3*(1 - b))

def test_horner():
    assert horner(0) == 0
    assert horner(1) == 1
    assert horner(x) == x

    assert horner(x + 1) == x + 1
    assert horner(x**2 + 1) == x**2 + 1
    assert horner(x**2 + x) == (x + 1)*x
    assert horner(x**2 + x + 1) == (x + 1)*x + 1

    assert horner(9*x**4 + 8*x**3 + 7*x**2 + 6*x + 5) == (((9*x + 8)*x + 7)*x + 6)*x + 5
    assert horner(a*x**4 + b*x**3 + c*x**2 + d*x + e) == (((a*x + b)*x + c)*x + d)*x + e

    assert horner(4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y, wrt=x) == ((4*y + 2)*x*y + (2*y + 1)*y)*x
    assert horner(4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y, wrt=y) == ((4*x + 2)*y*x + (2*x + 1)*x)*y

def test_poly():
    assert poly(x) == Poly(x, x)
    assert poly(y) == Poly(y, y)

    assert poly(x + y) == Poly(x + y, x, y)
    assert poly(x + sin(x)) == Poly(x + sin(x), x, sin(x))

    assert poly(x + y, wrt=y) == Poly(x + y, y, x)
    assert poly(x + sin(x), wrt=sin(x)) == Poly(x + sin(x), sin(x), x)

    assert poly(x*y + 2*x*z**2 + 17) == Poly(x*y + 2*x*z**2 + 17, x, y, z)

    assert poly(2*(y + z)**2 - 1) == Poly(2*y**2 + 4*y*z + 2*z**2 - 1, y, z)
    assert poly(x*(y + z)**2 - 1) == Poly(x*y**2 + 2*x*y*z + x*z**2 - 1, x, y, z)
    assert poly(2*x*(y + z)**2 - 1) == Poly(2*x*y**2 + 4*x*y*z + 2*x*z**2 - 1, x, y, z)

    assert poly(2*(y + z)**2 - x - 1) == Poly(2*y**2 + 4*y*z + 2*z**2 - x - 1, x, y, z)
    assert poly(x*(y + z)**2 - x - 1) == Poly(x*y**2 + 2*x*y*z + x*z**2 - x - 1, x, y, z)
    assert poly(2*x*(y + z)**2 - x - 1) == Poly(2*x*y**2 + 4*x*y*z + 2*x*z**2 - x - 1, x, y, z)

    assert poly(x*y + (x + y)**2 + (x + z)**2) == \
        Poly(2*x*z + 3*x*y + y**2 + z**2 + 2*x**2, x, y, z)
    assert poly(x*y*(x + y)*(x + z)**2) == \
        Poly(x**3*y**2 + x*y**2*z**2 + y*x**2*z**2 + 2*z*x**2*y**2 + 2*y*z*x**3 + y*x**4, x, y, z)

    assert poly(Poly(x + y + z, y, x, z)) == Poly(x + y + z, x, y, z)
    assert poly(Poly(x + y + z, y, x, z), wrt=z) == Poly(x + y + z, z, x, y)

    raises(GeneratorsNeeded, "poly(1)")

