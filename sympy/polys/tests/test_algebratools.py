"""Tests for classes defining properties of ground domains, e.g. ZZ, QQ, ZZ[x] ... """

from sympy.polys.algebratools import (
    ZZ, QQ, RR, CC, PolynomialRing, FractionField, EX, ZZ_sympy, QQ_sympy
)

from sympy.polys.polyerrors import (
    UnificationFailed,
    GeneratorsNeeded,
    DomainError,
)

from sympy.polys.polyclasses import DMF, DMP

from sympy import S, sqrt, sin, oo, raises, Integer, Rational, I, Real

from sympy.abc import x, y

ALG = QQ.algebraic_field(sqrt(2)+sqrt(3))

def test_Algebra__unify():
    assert ZZ.unify(ZZ) == ZZ
    assert QQ.unify(QQ) == QQ

    assert ZZ.unify(QQ) == QQ
    assert QQ.unify(ZZ) == QQ

    assert EX.unify(EX) == EX

    assert ZZ.unify(EX) == EX
    assert QQ.unify(EX) == EX
    assert EX.unify(ZZ) == EX
    assert EX.unify(QQ) == EX

    assert ZZ.poly_ring('x').unify(EX) == EX
    assert ZZ.frac_field('x').unify(EX) == EX
    assert EX.unify(ZZ.poly_ring('x')) == EX
    assert EX.unify(ZZ.frac_field('x')) == EX

    assert ZZ.poly_ring('x','y').unify(EX) == EX
    assert ZZ.frac_field('x','y').unify(EX) == EX
    assert EX.unify(ZZ.poly_ring('x','y')) == EX
    assert EX.unify(ZZ.frac_field('x','y')) == EX

    assert QQ.poly_ring('x').unify(EX) == EX
    assert QQ.frac_field('x').unify(EX) == EX
    assert EX.unify(QQ.poly_ring('x')) == EX
    assert EX.unify(QQ.frac_field('x')) == EX

    assert QQ.poly_ring('x','y').unify(EX) == EX
    assert QQ.frac_field('x','y').unify(EX) == EX
    assert EX.unify(QQ.poly_ring('x','y')) == EX
    assert EX.unify(QQ.frac_field('x','y')) == EX

    assert ZZ.poly_ring('x').unify(ZZ) == ZZ.poly_ring('x')
    assert ZZ.poly_ring('x').unify(QQ) == QQ.poly_ring('x')
    assert QQ.poly_ring('x').unify(ZZ) == QQ.poly_ring('x')
    assert QQ.poly_ring('x').unify(QQ) == QQ.poly_ring('x')

    assert ZZ.unify(ZZ.poly_ring('x')) == ZZ.poly_ring('x')
    assert QQ.unify(ZZ.poly_ring('x')) == QQ.poly_ring('x')
    assert ZZ.unify(QQ.poly_ring('x')) == QQ.poly_ring('x')
    assert QQ.unify(QQ.poly_ring('x')) == QQ.poly_ring('x')

    assert ZZ.poly_ring('x','y').unify(ZZ) == ZZ.poly_ring('x','y')
    assert ZZ.poly_ring('x','y').unify(QQ) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x','y').unify(ZZ) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x','y').unify(QQ) == QQ.poly_ring('x','y')

    assert ZZ.unify(ZZ.poly_ring('x','y')) == ZZ.poly_ring('x','y')
    assert QQ.unify(ZZ.poly_ring('x','y')) == QQ.poly_ring('x','y')
    assert ZZ.unify(QQ.poly_ring('x','y')) == QQ.poly_ring('x','y')
    assert QQ.unify(QQ.poly_ring('x','y')) == QQ.poly_ring('x','y')

    assert ZZ.frac_field('x').unify(ZZ) == ZZ.frac_field('x')
    assert ZZ.frac_field('x').unify(QQ) == EX # QQ.frac_field('x')
    assert QQ.frac_field('x').unify(ZZ) == EX # QQ.frac_field('x')
    assert QQ.frac_field('x').unify(QQ) == QQ.frac_field('x')

    assert ZZ.unify(ZZ.frac_field('x')) == ZZ.frac_field('x')
    assert QQ.unify(ZZ.frac_field('x')) == EX # QQ.frac_field('x')
    assert ZZ.unify(QQ.frac_field('x')) == EX # QQ.frac_field('x')
    assert QQ.unify(QQ.frac_field('x')) == QQ.frac_field('x')

    assert ZZ.frac_field('x','y').unify(ZZ) == ZZ.frac_field('x','y')
    assert ZZ.frac_field('x','y').unify(QQ) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(ZZ) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(QQ) == QQ.frac_field('x','y')

    assert ZZ.unify(ZZ.frac_field('x','y')) == ZZ.frac_field('x','y')
    assert QQ.unify(ZZ.frac_field('x','y')) == EX # QQ.frac_field('x','y')
    assert ZZ.unify(QQ.frac_field('x','y')) == EX # QQ.frac_field('x','y')
    assert QQ.unify(QQ.frac_field('x','y')) == QQ.frac_field('x','y')

    assert ZZ.poly_ring('x').unify(ZZ.poly_ring('x')) == ZZ.poly_ring('x')
    assert ZZ.poly_ring('x').unify(QQ.poly_ring('x')) == QQ.poly_ring('x')
    assert QQ.poly_ring('x').unify(ZZ.poly_ring('x')) == QQ.poly_ring('x')
    assert QQ.poly_ring('x').unify(QQ.poly_ring('x')) == QQ.poly_ring('x')

    assert ZZ.poly_ring('x','y').unify(ZZ.poly_ring('x')) == ZZ.poly_ring('x','y')
    assert ZZ.poly_ring('x','y').unify(QQ.poly_ring('x')) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x','y').unify(ZZ.poly_ring('x')) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x','y').unify(QQ.poly_ring('x')) == QQ.poly_ring('x','y')

    assert ZZ.poly_ring('x').unify(ZZ.poly_ring('x','y')) == ZZ.poly_ring('x','y')
    assert ZZ.poly_ring('x').unify(QQ.poly_ring('x','y')) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x').unify(ZZ.poly_ring('x','y')) == QQ.poly_ring('x','y')
    assert QQ.poly_ring('x').unify(QQ.poly_ring('x','y')) == QQ.poly_ring('x','y')

    assert ZZ.poly_ring('x','y').unify(ZZ.poly_ring('x','z')) == ZZ.poly_ring('x','y','z')
    assert ZZ.poly_ring('x','y').unify(QQ.poly_ring('x','z')) == QQ.poly_ring('x','y','z')
    assert QQ.poly_ring('x','y').unify(ZZ.poly_ring('x','z')) == QQ.poly_ring('x','y','z')
    assert QQ.poly_ring('x','y').unify(QQ.poly_ring('x','z')) == QQ.poly_ring('x','y','z')

    assert ZZ.frac_field('x').unify(ZZ.frac_field('x')) == ZZ.frac_field('x')
    assert ZZ.frac_field('x').unify(QQ.frac_field('x')) == QQ.frac_field('x')
    assert QQ.frac_field('x').unify(ZZ.frac_field('x')) == QQ.frac_field('x')
    assert QQ.frac_field('x').unify(QQ.frac_field('x')) == QQ.frac_field('x')

    assert ZZ.frac_field('x','y').unify(ZZ.frac_field('x')) == ZZ.frac_field('x','y')
    assert ZZ.frac_field('x','y').unify(QQ.frac_field('x')) == QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(ZZ.frac_field('x')) == QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(QQ.frac_field('x')) == QQ.frac_field('x','y')

    assert ZZ.frac_field('x').unify(ZZ.frac_field('x','y')) == ZZ.frac_field('x','y')
    assert ZZ.frac_field('x').unify(QQ.frac_field('x','y')) == QQ.frac_field('x','y')
    assert QQ.frac_field('x').unify(ZZ.frac_field('x','y')) == QQ.frac_field('x','y')
    assert QQ.frac_field('x').unify(QQ.frac_field('x','y')) == QQ.frac_field('x','y')

    assert ZZ.frac_field('x','y').unify(ZZ.frac_field('x','z')) == ZZ.frac_field('x','y','z')
    assert ZZ.frac_field('x','y').unify(QQ.frac_field('x','z')) == QQ.frac_field('x','y','z')
    assert QQ.frac_field('x','y').unify(ZZ.frac_field('x','z')) == QQ.frac_field('x','y','z')
    assert QQ.frac_field('x','y').unify(QQ.frac_field('x','z')) == QQ.frac_field('x','y','z')

    assert ZZ.poly_ring('x').unify(ZZ.frac_field('x')) == ZZ.frac_field('x')
    assert ZZ.poly_ring('x').unify(QQ.frac_field('x')) == EX # QQ.frac_field('x')
    assert QQ.poly_ring('x').unify(ZZ.frac_field('x')) == EX # QQ.frac_field('x')
    assert QQ.poly_ring('x').unify(QQ.frac_field('x')) == QQ.frac_field('x')

    assert ZZ.poly_ring('x','y').unify(ZZ.frac_field('x')) == ZZ.frac_field('x','y')
    assert ZZ.poly_ring('x','y').unify(QQ.frac_field('x')) == EX # QQ.frac_field('x','y')
    assert QQ.poly_ring('x','y').unify(ZZ.frac_field('x')) == EX # QQ.frac_field('x','y')
    assert QQ.poly_ring('x','y').unify(QQ.frac_field('x')) == QQ.frac_field('x','y')

    assert ZZ.poly_ring('x').unify(ZZ.frac_field('x','y')) == ZZ.frac_field('x','y')
    assert ZZ.poly_ring('x').unify(QQ.frac_field('x','y')) == EX # QQ.frac_field('x','y')
    assert QQ.poly_ring('x').unify(ZZ.frac_field('x','y')) == EX # QQ.frac_field('x','y')
    assert QQ.poly_ring('x').unify(QQ.frac_field('x','y')) == QQ.frac_field('x','y')

    assert ZZ.poly_ring('x','y').unify(ZZ.frac_field('x','z')) == ZZ.frac_field('x','y','z')
    assert ZZ.poly_ring('x','y').unify(QQ.frac_field('x','z')) == EX # QQ.frac_field('x','y','z')
    assert QQ.poly_ring('x','y').unify(ZZ.frac_field('x','z')) == EX # QQ.frac_field('x','y','z')
    assert QQ.poly_ring('x','y').unify(QQ.frac_field('x','z')) == QQ.frac_field('x','y','z')

    assert ZZ.frac_field('x').unify(ZZ.poly_ring('x')) == ZZ.frac_field('x')
    assert ZZ.frac_field('x').unify(QQ.poly_ring('x')) == EX # QQ.frac_field('x')
    assert QQ.frac_field('x').unify(ZZ.poly_ring('x')) == EX # QQ.frac_field('x')
    assert QQ.frac_field('x').unify(QQ.poly_ring('x')) == QQ.frac_field('x')

    assert ZZ.frac_field('x','y').unify(ZZ.poly_ring('x')) == ZZ.frac_field('x','y')
    assert ZZ.frac_field('x','y').unify(QQ.poly_ring('x')) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(ZZ.poly_ring('x')) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x','y').unify(QQ.poly_ring('x')) == QQ.frac_field('x','y')

    assert ZZ.frac_field('x').unify(ZZ.poly_ring('x','y')) == ZZ.frac_field('x','y')
    assert ZZ.frac_field('x').unify(QQ.poly_ring('x','y')) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x').unify(ZZ.poly_ring('x','y')) == EX # QQ.frac_field('x','y')
    assert QQ.frac_field('x').unify(QQ.poly_ring('x','y')) == QQ.frac_field('x','y')

    assert ZZ.frac_field('x','y').unify(ZZ.poly_ring('x','z')) == ZZ.frac_field('x','y','z')
    assert ZZ.frac_field('x','y').unify(QQ.poly_ring('x','z')) == EX # QQ.frac_field('x','y','z')
    assert QQ.frac_field('x','y').unify(ZZ.poly_ring('x','z')) == EX # QQ.frac_field('x','y','z')
    assert QQ.frac_field('x','y').unify(QQ.poly_ring('x','z')) == QQ.frac_field('x','y','z')

    raises(UnificationFailed, "ZZ.poly_ring('x','y').unify(ZZ, gens=('y', 'z'))")
    raises(UnificationFailed, "ZZ.unify(ZZ.poly_ring('x','y'), gens=('y', 'z'))")

def test_Algebra__contains__():
    assert (0 in EX) == True
    assert (0 in ZZ) == True
    assert (0 in QQ) == True
    assert (0 in RR) == True
    assert (0 in ALG) == True
    assert (0 in ZZ[x,y]) == True
    assert (0 in QQ[x,y]) == True
    assert (0 in RR[x,y]) == True

    assert (-7 in EX) == True
    assert (-7 in ZZ) == True
    assert (-7 in QQ) == True
    assert (-7 in RR) == True
    assert (-7 in ALG) == True
    assert (-7 in ZZ[x,y]) == True
    assert (-7 in QQ[x,y]) == True
    assert (-7 in RR[x,y]) == True

    assert (17 in EX) == True
    assert (17 in ZZ) == True
    assert (17 in QQ) == True
    assert (17 in RR) == True
    assert (17 in ALG) == True
    assert (17 in ZZ[x,y]) == True
    assert (17 in QQ[x,y]) == True
    assert (17 in RR[x,y]) == True

    assert (-S(1)/7 in EX) == True
    assert (-S(1)/7 in ZZ) == False
    assert (-S(1)/7 in QQ) == True
    assert (-S(1)/7 in RR) == True
    assert (-S(1)/7 in ALG) == True
    assert (-S(1)/7 in ZZ[x,y]) == False
    assert (-S(1)/7 in QQ[x,y]) == True
    assert (-S(1)/7 in RR[x,y]) == True

    assert (S(3)/5 in EX) == True
    assert (S(3)/5 in ZZ) == False
    assert (S(3)/5 in QQ) == True
    assert (S(3)/5 in RR) == True
    assert (S(3)/5 in ALG) == True
    assert (S(3)/5 in ZZ[x,y]) == False
    assert (S(3)/5 in QQ[x,y]) == True
    assert (S(3)/5 in RR[x,y]) == True

    assert (3.0 in EX) == True
    assert (3.0 in ZZ) == True
    assert (3.0 in QQ) == True
    assert (3.0 in RR) == True
    assert (3.0 in ALG) == True
    assert (3.0 in ZZ[x,y]) == True
    assert (3.0 in QQ[x,y]) == True
    assert (3.0 in RR[x,y]) == True

    assert (3.14 in EX) == True
    assert (3.14 in ZZ) == False
    assert (3.14 in QQ) == True
    assert (3.14 in RR) == True
    assert (3.14 in ALG) == True
    assert (3.14 in ZZ[x,y]) == False
    assert (3.14 in QQ[x,y]) == True
    assert (3.14 in RR[x,y]) == True

    assert (oo in EX) == True
    assert (oo in ZZ) == False
    assert (oo in QQ) == False
    assert (oo in RR) == False
    assert (oo in ALG) == False
    assert (oo in ZZ[x,y]) == False
    assert (oo in QQ[x,y]) == False
    assert (oo in RR[x,y]) == False

    assert (-oo in EX) == True
    assert (-oo in ZZ) == False
    assert (-oo in QQ) == False
    assert (-oo in RR) == False
    assert (-oo in ALG) == False
    assert (-oo in ZZ[x,y]) == False
    assert (-oo in QQ[x,y]) == False
    assert (-oo in RR[x,y]) == False

    assert (sqrt(7) in EX) == True
    assert (sqrt(7) in ZZ) == False
    assert (sqrt(7) in QQ) == False
    assert (sqrt(7) in RR) == True
    assert (sqrt(7) in ALG) == False
    assert (sqrt(7) in ZZ[x,y]) == False
    assert (sqrt(7) in QQ[x,y]) == False
    assert (sqrt(7) in RR[x,y]) == True

    assert (2*sqrt(3)+1 in EX) == True
    assert (2*sqrt(3)+1 in ZZ) == False
    assert (2*sqrt(3)+1 in QQ) == False
    assert (2*sqrt(3)+1 in RR) == True
    assert (2*sqrt(3)+1 in ALG) == True
    assert (2*sqrt(3)+1 in ZZ[x,y]) == False
    assert (2*sqrt(3)+1 in QQ[x,y]) == False
    assert (2*sqrt(3)+1 in RR[x,y]) == True

    assert (sin(1) in EX) == True
    assert (sin(1) in ZZ) == False
    assert (sin(1) in QQ) == False
    assert (sin(1) in RR) == True
    assert (sin(1) in ALG) == False
    assert (sin(1) in ZZ[x,y]) == False
    assert (sin(1) in QQ[x,y]) == False
    assert (sin(1) in RR[x,y]) == True

    assert (x**2 + 1 in EX) == True
    assert (x**2 + 1 in ZZ) == False
    assert (x**2 + 1 in QQ) == False
    assert (x**2 + 1 in RR) == False
    assert (x**2 + 1 in ALG) == False
    assert (x**2 + 1 in ZZ[x]) == True
    assert (x**2 + 1 in QQ[x]) == True
    assert (x**2 + 1 in RR[x]) == True
    assert (x**2 + 1 in ZZ[x,y]) == True
    assert (x**2 + 1 in QQ[x,y]) == True
    assert (x**2 + 1 in RR[x,y]) == True

    assert (x**2 + y**2 in EX) == True
    assert (x**2 + y**2 in ZZ) == False
    assert (x**2 + y**2 in QQ) == False
    assert (x**2 + y**2 in RR) == False
    assert (x**2 + y**2 in ALG) == False
    assert (x**2 + y**2 in ZZ[x]) == False
    assert (x**2 + y**2 in QQ[x]) == False
    assert (x**2 + y**2 in RR[x]) == False
    assert (x**2 + y**2 in ZZ[x,y]) == True
    assert (x**2 + y**2 in QQ[x,y]) == True
    assert (x**2 + y**2 in RR[x,y]) == True

def test_Algebra_get_ring():
    assert ZZ.has_assoc_Ring == True
    assert QQ.has_assoc_Ring == True
    assert ZZ[x].has_assoc_Ring == True
    assert QQ[x].has_assoc_Ring == True
    assert ZZ[x,y].has_assoc_Ring == True
    assert QQ[x,y].has_assoc_Ring == True
    assert ZZ.frac_field(x).has_assoc_Ring == True
    assert QQ.frac_field(x).has_assoc_Ring == True
    assert ZZ.frac_field(x,y).has_assoc_Ring == True
    assert QQ.frac_field(x,y).has_assoc_Ring == True

    assert EX.has_assoc_Ring == False
    assert RR.has_assoc_Ring == False
    assert ALG.has_assoc_Ring == False

    assert ZZ.get_ring() == ZZ
    assert QQ.get_ring() == ZZ
    assert ZZ[x].get_ring() == ZZ[x]
    assert QQ[x].get_ring() == QQ[x]
    assert ZZ[x,y].get_ring() == ZZ[x,y]
    assert QQ[x,y].get_ring() == QQ[x,y]
    assert ZZ.frac_field(x).get_ring() == ZZ[x]
    assert QQ.frac_field(x).get_ring() == QQ[x]
    assert ZZ.frac_field(x,y).get_ring() == ZZ[x,y]
    assert QQ.frac_field(x,y).get_ring() == QQ[x,y]

    raises(DomainError, "EX.get_ring()")
    raises(DomainError, "RR.get_ring()")
    raises(DomainError, "ALG.get_ring()")

def test_Algebra_get_field():
    assert EX.has_assoc_Field == True
    assert ZZ.has_assoc_Field == True
    assert QQ.has_assoc_Field == True
    assert ALG.has_assoc_Field == True
    assert ZZ[x].has_assoc_Field == True
    assert QQ[x].has_assoc_Field == True
    assert ZZ[x,y].has_assoc_Field == True
    assert QQ[x,y].has_assoc_Field == True

    assert RR.has_assoc_Field == False

    assert EX.get_field() == EX
    assert ZZ.get_field() == QQ
    assert QQ.get_field() == QQ
    assert ALG.get_field() == ALG
    assert ZZ[x].get_field() == ZZ.frac_field(x)
    assert QQ[x].get_field() == QQ.frac_field(x)
    assert ZZ[x,y].get_field() == ZZ.frac_field(x,y)
    assert QQ[x,y].get_field() == QQ.frac_field(x,y)

    raises(DomainError, "RR.get_field()")

def test_Algebra_get_exact():
    assert EX.get_exact() == EX
    assert ZZ.get_exact() == ZZ
    assert QQ.get_exact() == QQ
    assert RR.get_exact() == QQ
    assert ALG.get_exact() == ALG
    assert ZZ[x].get_exact() == ZZ[x]
    assert QQ[x].get_exact() == QQ[x]
    assert ZZ[x,y].get_exact() == ZZ[x,y]
    assert QQ[x,y].get_exact() == QQ[x,y]
    assert ZZ.frac_field(x).get_exact() == ZZ.frac_field(x)
    assert QQ.frac_field(x).get_exact() == QQ.frac_field(x)
    assert ZZ.frac_field(x,y).get_exact() == ZZ.frac_field(x,y)
    assert QQ.frac_field(x,y).get_exact() == QQ.frac_field(x,y)

def test_Algebra_convert():
    assert QQ.convert(10e-52) != QQ(0)

def test_PolynomialRing__init():
    raises(GeneratorsNeeded, "ZZ.poly_ring()")

def test_PolynomialRing_from_FractionField():
    x = DMF(([1, 0, 1], [1, 1]), ZZ)
    y = DMF(([1, 0, 1], [1]), ZZ)
    assert ZZ['x'].from_FractionField(x, ZZ['x']) is None
    assert ZZ['x'].from_FractionField(y, ZZ['x']) == \
        DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)

def test_FractionField__init():
    raises(GeneratorsNeeded, "ZZ.frac_field()")

def test_sympy_of_type():
    assert ZZ_sympy().of_type(Integer(1))
    assert ZZ_sympy().of_type(Integer(0))
    assert ZZ_sympy().of_type(Integer(-1))
    assert ZZ_sympy().of_type(Integer(2))
    assert not ZZ_sympy().of_type(Rational(1, 2))
    assert QQ_sympy().of_type(Rational(1))
    assert QQ_sympy().of_type(Rational(-1))
    assert QQ_sympy().of_type(Rational(0))
    assert QQ_sympy().of_type(Rational(2))
    assert QQ_sympy().of_type(Rational(1, 2))
    assert QQ_sympy().of_type(Rational(3, 2))

def test_CC_to_from_sympy():
    assert CC.from_sympy(Real(1.0)) == 1+0j
    assert CC.from_sympy(I) == 1j
    assert CC.from_sympy(2*I) == 2j
    assert CC.from_sympy(1 + 2*I) == 1+2j
    assert CC.to_sympy(1+0j) == S.One
    assert CC.to_sympy(1+2j) == S.One + S(2)*I
    assert CC.to_sympy(2j) == S(2)*I

def test___eq__():
    assert not QQ['x'] == ZZ['x']
    assert not QQ.frac_field(x) == ZZ.frac_field(x)
