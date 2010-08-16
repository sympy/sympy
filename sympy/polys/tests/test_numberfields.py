"""Tests for computational algebraic number field theory. """

from sympy import S, Rational, Symbol, Poly, raises, sin, sqrt, I, oo, pure

from sympy.polys.numberfields import (
    minimal_polynomial,
    primitive_element,
    is_isomorphism_possible,
    field_isomorphism_pslq,
    field_isomorphism,
    to_number_field,
    AlgebraicNumber,
    isolate,
)

from sympy.polys.polyerrors import (
    IsomorphismFailed,
    NotAlgebraic,
)

from sympy.polys.polyclasses import DMP
from sympy.polys.domains import QQ

from sympy.abc import x, y

Q = Rational

def test_minimal_polynomial():
    assert minimal_polynomial(-7, x) == x + 7
    assert minimal_polynomial(-1, x) == x + 1
    assert minimal_polynomial( 0, x) == x
    assert minimal_polynomial( 1, x) == x - 1
    assert minimal_polynomial( 7, x) == x - 7

    assert minimal_polynomial(sqrt(2), x) == x**2 - 2
    assert minimal_polynomial(sqrt(5), x) == x**2 - 5
    assert minimal_polynomial(sqrt(6), x) == x**2 - 6

    assert minimal_polynomial(2*sqrt(2), x) == x**2 - 8
    assert minimal_polynomial(3*sqrt(5), x) == x**2 - 45
    assert minimal_polynomial(4*sqrt(6), x) == x**2 - 96

    assert minimal_polynomial(2*sqrt(2) + 3, x) == x**2 -  6*x +  1
    assert minimal_polynomial(3*sqrt(5) + 6, x) == x**2 - 12*x -  9
    assert minimal_polynomial(4*sqrt(6) + 7, x) == x**2 - 14*x - 47

    assert minimal_polynomial(2*sqrt(2) - 3, x) == x**2 +  6*x +  1
    assert minimal_polynomial(3*sqrt(5) - 6, x) == x**2 + 12*x -  9
    assert minimal_polynomial(4*sqrt(6) - 7, x) == x**2 + 14*x - 47

    assert minimal_polynomial(sqrt(1 + sqrt(6)), x) == x**4 -  2*x**2 -  5
    assert minimal_polynomial(sqrt(I + sqrt(6)), x) == x**8 - 10*x**4 + 49

    assert minimal_polynomial(2*I + sqrt(2 + I), x) == x**4 + 4*x**2 + 8*x + 37

    assert minimal_polynomial(sqrt(2) + sqrt(3), x) == x**4 - 10*x**2 + 1
    assert minimal_polynomial(sqrt(2) + sqrt(3) + sqrt(6), x) == x**4 - 22*x**2 - 48*x - 23

    a = 1 - 9*sqrt(2) + 7*sqrt(3)

    assert minimal_polynomial(1/a, x) == 392*x**4 - 1232*x**3 + 612*x**2 + 4*x - 1
    assert minimal_polynomial(1/sqrt(a), x) == 392*x**8 - 1232*x**6 + 612*x**4 + 4*x**2 - 1

    raises(NotAlgebraic, "minimal_polynomial(y, x)")
    raises(NotAlgebraic, "minimal_polynomial(oo, x)")
    raises(NotAlgebraic, "minimal_polynomial(2**y, x)")
    raises(NotAlgebraic, "minimal_polynomial(sin(1), x)")

    assert minimal_polynomial(sqrt(2)) == pure**2 - 2
    assert minimal_polynomial(sqrt(2), x) == x**2 - 2

    assert minimal_polynomial(sqrt(2), polys=True) == Poly(pure**2 - 2)
    assert minimal_polynomial(sqrt(2), x, polys=True) == Poly(x**2 - 2)

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(3))

    assert minimal_polynomial(a, x) == x**2 - 2
    assert minimal_polynomial(b, x) == x**2 - 3

    assert minimal_polynomial(a, x, polys=True) == Poly(x**2 - 2)
    assert minimal_polynomial(b, x, polys=True) == Poly(x**2 - 3)

    assert minimal_polynomial(sqrt(a/2 + 17), x) == 2*x**4 -  68*x**2 +  577
    assert minimal_polynomial(sqrt(b/2 + 17), x) == 4*x**4 - 136*x**2 + 1153

    a, b = sqrt(2)/3 + 7, AlgebraicNumber(sqrt(2)/3 + 7)

    f = 81*x**8 - 2268*x**6 - 4536*x**5 + 22644*x**4 + 63216*x**3 - 31608*x**2 - 189648*x + 141358

    assert minimal_polynomial(sqrt(a) + sqrt(sqrt(a)), x) == f
    assert minimal_polynomial(sqrt(b) + sqrt(sqrt(b)), x) == f

    assert minimal_polynomial(a**Rational(3, 2), x) == 729*x**4 - 506898*x**2 + 84604519

def test_primitive_element():
    assert primitive_element([sqrt(2)], x) == (x**2 - 2, [1])
    assert primitive_element([sqrt(2), sqrt(3)], x) == (x**4 - 10*x**2 + 1, [1, 1])

    assert primitive_element([sqrt(2)], x, polys=True) == (Poly(x**2 - 2), [1])
    assert primitive_element([sqrt(2), sqrt(3)], x, polys=True) == (Poly(x**4 - 10*x**2 + 1), [1, 1])

    assert primitive_element([sqrt(2)], x, ex=True) == (x**2 - 2, [1], [[1, 0]])
    assert primitive_element([sqrt(2), sqrt(3)], x, ex=True) == \
        (x**4 - 10*x**2 + 1, [1, 1], [[Q(1,2), 0, -Q(9,2), 0], [-Q(1,2), 0, Q(11,2), 0]])

    assert primitive_element([sqrt(2)], x, ex=True, polys=True) == (Poly(x**2 - 2), [1], [[1, 0]])
    assert primitive_element([sqrt(2), sqrt(3)], x, ex=True, polys=True) == \
        (Poly(x**4 - 10*x**2 + 1), [1, 1], [[Q(1,2), 0, -Q(9,2), 0], [-Q(1,2), 0, Q(11,2), 0]])

    assert primitive_element([sqrt(2)]) == (pure**2 - 2, [1])

    raises(ValueError, "primitive_element([], x, ex=False)")
    raises(ValueError, "primitive_element([], x, ex=True)")

def test_field_isomorphism_pslq():
    a = AlgebraicNumber(I)
    b = AlgebraicNumber(I*sqrt(3))

    raises(NotImplementedError, "field_isomorphism_pslq(a, b)")

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(3))
    c = AlgebraicNumber(sqrt(7))
    d = AlgebraicNumber(sqrt(2)+sqrt(3))
    e = AlgebraicNumber(sqrt(2)+sqrt(3)+sqrt(7))

    assert field_isomorphism_pslq(a, a) == [1, 0]
    assert field_isomorphism_pslq(a, b) == None
    assert field_isomorphism_pslq(a, c) == None
    assert field_isomorphism_pslq(a, d) == [Q(1,2), 0, -Q(9,2), 0]
    assert field_isomorphism_pslq(a, e) == [Q(1,80), 0, -Q(1,2), 0, Q(59,20), 0]

    assert field_isomorphism_pslq(b, a) == None
    assert field_isomorphism_pslq(b, b) == [1, 0]
    assert field_isomorphism_pslq(b, c) == None
    assert field_isomorphism_pslq(b, d) == [-Q(1,2), 0, Q(11,2), 0]
    assert field_isomorphism_pslq(b, e) == [-Q(3,640), 0, Q(67,320), 0, -Q(297,160), 0, Q(313,80), 0]

    assert field_isomorphism_pslq(c, a) == None
    assert field_isomorphism_pslq(c, b) == None
    assert field_isomorphism_pslq(c, c) == [1, 0]
    assert field_isomorphism_pslq(c, d) == None
    assert field_isomorphism_pslq(c, e) == [Q(3,640), 0, -Q(71,320), 0, Q(377,160), 0, -Q(469,80), 0]

    assert field_isomorphism_pslq(d, a) == None
    assert field_isomorphism_pslq(d, b) == None
    assert field_isomorphism_pslq(d, c) == None
    assert field_isomorphism_pslq(d, d) == [1, 0]
    assert field_isomorphism_pslq(d, e) == [-Q(3,640), 0, Q(71,320), 0, -Q(377,160), 0, Q(549,80), 0]

    assert field_isomorphism_pslq(e, a) == None
    assert field_isomorphism_pslq(e, b) == None
    assert field_isomorphism_pslq(e, c) == None
    assert field_isomorphism_pslq(e, d) == None
    assert field_isomorphism_pslq(e, e) == [1, 0]

    f = AlgebraicNumber(3*sqrt(2)+8*sqrt(7)-5)

    assert field_isomorphism_pslq(f, e) == [Q(3,80), 0, -Q(139,80), 0, Q(347,20), 0, -Q(761,20), -5]

def test_field_isomorphism():
    assert field_isomorphism(3, sqrt(2)) == [3]

    assert field_isomorphism( I*sqrt(3), I*sqrt(3)/2) == [ 2, 0]
    assert field_isomorphism(-I*sqrt(3), I*sqrt(3)/2) == [-2, 0]

    assert field_isomorphism( I*sqrt(3),-I*sqrt(3)/2) == [-2, 0]
    assert field_isomorphism(-I*sqrt(3),-I*sqrt(3)/2) == [ 2, 0]

    assert field_isomorphism( 2*I*sqrt(3)/7, 5*I*sqrt(3)/3) == [ S(6)/35, 0]
    assert field_isomorphism(-2*I*sqrt(3)/7, 5*I*sqrt(3)/3) == [-S(6)/35, 0]

    assert field_isomorphism( 2*I*sqrt(3)/7,-5*I*sqrt(3)/3) == [-S(6)/35, 0]
    assert field_isomorphism(-2*I*sqrt(3)/7,-5*I*sqrt(3)/3) == [ S(6)/35, 0]

    assert field_isomorphism( 2*I*sqrt(3)/7+27, 5*I*sqrt(3)/3) == [ S(6)/35, 27]
    assert field_isomorphism(-2*I*sqrt(3)/7+27, 5*I*sqrt(3)/3) == [-S(6)/35, 27]

    assert field_isomorphism( 2*I*sqrt(3)/7+27,-5*I*sqrt(3)/3) == [-S(6)/35, 27]
    assert field_isomorphism(-2*I*sqrt(3)/7+27,-5*I*sqrt(3)/3) == [ S(6)/35, 27]

    p = AlgebraicNumber( sqrt(2) + sqrt(3))
    q = AlgebraicNumber(-sqrt(2) + sqrt(3))
    r = AlgebraicNumber( sqrt(2) - sqrt(3))
    s = AlgebraicNumber(-sqrt(2) - sqrt(3))

    pos_coeffs = [ S(1)/2, S(0), -S(9)/2, S(0)]
    neg_coeffs = [-S(1)/2, S(0),  S(9)/2, S(0)]

    a = AlgebraicNumber(sqrt(2))

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == pos_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_coeffs

    a = AlgebraicNumber(-sqrt(2))

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == pos_coeffs
    assert field_isomorphism(a, r, fast=True) == neg_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == pos_coeffs
    assert field_isomorphism(a, r, fast=False) == neg_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    pos_coeffs = [ S(1)/2, S(0), -S(11)/2, S(0)]
    neg_coeffs = [-S(1)/2, S(0),  S(11)/2, S(0)]

    a = AlgebraicNumber(sqrt(3))

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    a = AlgebraicNumber(-sqrt(3))

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == pos_coeffs
    assert field_isomorphism(a, q, fast=True) == pos_coeffs
    assert field_isomorphism(a, r, fast=True) == neg_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_coeffs
    assert field_isomorphism(a, q, fast=False) == pos_coeffs
    assert field_isomorphism(a, r, fast=False) == neg_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_coeffs

    pos_coeffs = [ S(3)/2, S(0), -S(33)/2, -S(8)]
    neg_coeffs = [-S(3)/2, S(0),  S(33)/2, -S(8)]

    a = AlgebraicNumber(3*sqrt(3)-8)

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    a = AlgebraicNumber(3*sqrt(2)+2*sqrt(3)+1)

    pos_1_coeffs = [ S(1)/2, S(0), -S(5)/2,  S(1)]
    neg_5_coeffs = [-S(5)/2, S(0),  S(49)/2, S(1)]
    pos_5_coeffs = [ S(5)/2, S(0), -S(49)/2, S(1)]
    neg_1_coeffs = [-S(1)/2, S(0),  S(5)/2,  S(1)]

    assert is_isomorphism_possible(a, p) == True
    assert is_isomorphism_possible(a, q) == True
    assert is_isomorphism_possible(a, r) == True
    assert is_isomorphism_possible(a, s) == True

    assert field_isomorphism(a, p, fast=True) == pos_1_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_5_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_5_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_1_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_1_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_5_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_5_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_1_coeffs

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(3))
    c = AlgebraicNumber(sqrt(7))

    assert is_isomorphism_possible(a, b) == True
    assert is_isomorphism_possible(b, a) == True

    assert is_isomorphism_possible(c, p) == False

    assert field_isomorphism(sqrt(2), sqrt(3), fast=True) is None
    assert field_isomorphism(sqrt(3), sqrt(2), fast=True) is None

    assert field_isomorphism(sqrt(2), sqrt(3), fast=False) is None
    assert field_isomorphism(sqrt(3), sqrt(2), fast=False) is None

def test_to_number_field():
    assert to_number_field(sqrt(2)) == AlgebraicNumber(sqrt(2))
    assert to_number_field([sqrt(2), sqrt(3)]) == AlgebraicNumber(sqrt(2)+sqrt(3))

    a = AlgebraicNumber(sqrt(2)+sqrt(3), [S(1)/2, S(0), -S(9)/2, S(0)])

    assert to_number_field(sqrt(2), sqrt(2)+sqrt(3)) == a
    assert to_number_field(sqrt(2), AlgebraicNumber(sqrt(2)+sqrt(3))) == a

    raises(IsomorphismFailed, "to_number_field(sqrt(2), sqrt(3))")

def test_AlgebraicNumber():
    minpoly, root = x**2 - 2, sqrt(2)

    a = AlgebraicNumber(root, gen=x)

    assert a.rep == DMP([QQ(1),QQ(0)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly

    assert a.is_aliased == False

    assert a.coeffs() == [S(1), S(0)]
    assert a.native_coeffs() == [QQ(1), QQ(0)]

    a = AlgebraicNumber(root, gen=x, alias='y')

    assert a.rep == DMP([QQ(1),QQ(0)], QQ)
    assert a.root == root
    assert a.alias == Symbol('y')
    assert a.minpoly == minpoly

    assert a.is_aliased == True

    a = AlgebraicNumber(root, gen=x, alias=Symbol('y'))

    assert a.rep == DMP([QQ(1),QQ(0)], QQ)
    assert a.root == root
    assert a.alias == Symbol('y')
    assert a.minpoly == minpoly

    assert a.is_aliased == True

    assert AlgebraicNumber(sqrt(2), []).rep == DMP([], QQ)

    assert AlgebraicNumber(sqrt(2), [8]).rep == DMP([QQ(8)], QQ)
    assert AlgebraicNumber(sqrt(2), [S(8)/3]).rep == DMP([QQ(8,3)], QQ)

    assert AlgebraicNumber(sqrt(2), [7, 3]).rep == DMP([QQ(7),QQ(3)], QQ)
    assert AlgebraicNumber(sqrt(2), [S(7)/9, S(3)/2]).rep == DMP([QQ(7,9),QQ(3,2)], QQ)

    assert AlgebraicNumber(sqrt(2), [1, 2, 3]).rep == DMP([QQ(2),QQ(5)], QQ)

    a = AlgebraicNumber(AlgebraicNumber(root, gen=x), [1,2])

    assert a.rep == DMP([QQ(1),QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly

    assert a.is_aliased == False

    assert a.coeffs() == [S(1), S(2)]
    assert a.native_coeffs() == [QQ(1), QQ(2)]

    a = AlgebraicNumber((minpoly, root), [1,2])

    assert a.rep == DMP([QQ(1),QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly

    assert a.is_aliased == False

    a = AlgebraicNumber((Poly(minpoly), root), [1,2])

    assert a.rep == DMP([QQ(1),QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly

    assert a.is_aliased == False

    assert AlgebraicNumber( sqrt(3)).rep == DMP([ QQ(1),QQ(0)], QQ)
    assert AlgebraicNumber(-sqrt(3)).rep == DMP([-QQ(1),QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(2))

    assert a == b and a == sqrt(2)

    a = AlgebraicNumber(sqrt(2), gen=x)
    b = AlgebraicNumber(sqrt(2), gen=x)

    assert a == b and a == sqrt(2)

    a = AlgebraicNumber(sqrt(2), [1,2])
    b = AlgebraicNumber(sqrt(2), [1,3])

    assert a != b and a != sqrt(2)+3

    assert (a == x) == False and (a != x) == True

    a = AlgebraicNumber(sqrt(2), [1,0])
    b = AlgebraicNumber(sqrt(2), [1,0], alias=y)

    assert a.as_poly(x) == Poly(x)
    assert b.as_poly()  == Poly(y)

    assert a.as_basic()  == sqrt(2)
    assert a.as_basic(x) == x
    assert b.as_basic()  == sqrt(2)
    assert b.as_basic(x) == x

    a = AlgebraicNumber(sqrt(2), [2,3])
    b = AlgebraicNumber(sqrt(2), [2,3], alias=y)

    p = a.as_poly()

    assert p == Poly(2*p.gen+3)

    assert a.as_poly(x) == Poly(2*x+3)
    assert b.as_poly()  == Poly(2*y+3)

    assert a.as_basic()  == 2*sqrt(2)+3
    assert a.as_basic(x) == 2*x+3
    assert b.as_basic()  == 2*sqrt(2)+3
    assert b.as_basic(x) == 2*x+3

def test_to_algebraic_integer():
    a = AlgebraicNumber(sqrt(3), gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 3
    assert a.root    == sqrt(3)
    assert a.rep     == DMP([QQ(1),QQ(0)], QQ)

    a = AlgebraicNumber(2*sqrt(3), gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root    == 2*sqrt(3)
    assert a.rep     == DMP([QQ(1),QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root    == 2*sqrt(3)
    assert a.rep     == DMP([QQ(1),QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, [S(7)/19, 3], gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root    == 2*sqrt(3)
    assert a.rep     == DMP([QQ(7,19),QQ(3)], QQ)

def test_isolate():
    assert isolate(1) == (1, 1)
    assert isolate(S(1)/2) == (S(1)/2, S(1)/2)

    assert isolate(sqrt(2)) == (1, 2)
    assert isolate(-sqrt(2)) == (-2, -1)

    assert isolate(sqrt(2), eps=S(1)/100) == (S(24)/17, S(17)/12)
    assert isolate(-sqrt(2), eps=S(1)/100) == (-S(17)/12, -S(24)/17)

    raises(NotImplementedError, "isolate(I)")

