"""Tests for classes defining properties of ground domains, e.g. ZZ, QQ,
ZZ[x] ... """

from sympy import S, sqrt, sin, oo, nan, Poly, Integer, Rational
from sympy.abc import x, y, z

from sympy.polys.domains import (
    ZZ, QQ, RR, FF, PythonRationalType as Q, ZZ_sympy, QQ_sympy,
    RR_mpmath, RR_sympy, PolynomialRing, FractionField, EX)

from sympy.polys.domains.modularinteger import ModularIntegerFactory

from sympy.polys.polyerrors import (
    UnificationFailed,
    GeneratorsNeeded,
    GeneratorsError,
    CoercionFailed,
    NotInvertible,
    DomainError)

from sympy.polys.polyclasses import DMP, DMF
from sympy.utilities.pytest import raises

ALG = QQ.algebraic_field(sqrt(2) + sqrt(3))

def test_Domain__unify():
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

    alg = QQ.algebraic_field(sqrt(5))

    assert alg.unify(alg['x','y']) == alg['x','y']
    assert alg['x','y'].unify(alg) == alg['x','y']

    assert alg.unify(alg.frac_field('x','y')) == alg.frac_field('x','y')
    assert alg.frac_field('x','y').unify(alg) == alg.frac_field('x','y')

    ext = QQ.algebraic_field(sqrt(7))

    raises(NotImplementedError, lambda: alg.unify(ext))

    raises(UnificationFailed, lambda: ZZ.poly_ring('x','y').unify(ZZ, gens=('y', 'z')))
    raises(UnificationFailed, lambda: ZZ.unify(ZZ.poly_ring('x','y'), gens=('y', 'z')))

def test_Domain__contains__():
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

    assert (S(3)/2*x/(y + 1) - z in QQ[x, y, z]) == False

def test_Domain_get_ring():
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

    raises(DomainError, lambda: EX.get_ring())
    raises(DomainError, lambda: RR.get_ring())
    raises(DomainError, lambda: ALG.get_ring())

def test_Domain_get_field():
    assert EX.has_assoc_Field == True
    assert ZZ.has_assoc_Field == True
    assert QQ.has_assoc_Field == True
    assert RR.has_assoc_Field == False
    assert ALG.has_assoc_Field == True
    assert ZZ[x].has_assoc_Field == True
    assert QQ[x].has_assoc_Field == True
    assert ZZ[x,y].has_assoc_Field == True
    assert QQ[x,y].has_assoc_Field == True

    assert EX.get_field() == EX
    assert ZZ.get_field() == QQ
    assert QQ.get_field() == QQ
    raises(DomainError, lambda: RR.get_field())
    assert ALG.get_field() == ALG
    assert ZZ[x].get_field() == ZZ.frac_field(x)
    assert QQ[x].get_field() == QQ.frac_field(x)
    assert ZZ[x,y].get_field() == ZZ.frac_field(x,y)
    assert QQ[x,y].get_field() == QQ.frac_field(x,y)

def test_Domain_get_exact():
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

def test_Domain_convert():
    assert QQ.convert(10e-52) != QQ(0)
    assert ZZ.convert(DMP([[ZZ(1)]], ZZ)) == ZZ(1)

def test_PolynomialRing__init():
    raises(GeneratorsNeeded, lambda: ZZ.poly_ring())

def test_PolynomialRing_from_FractionField():
    x = DMF(([1, 0, 1], [1, 1]), ZZ)
    y = DMF(([1, 0, 1], [1]), ZZ)

    assert ZZ['x'].from_FractionField(x, ZZ['x']) is None
    assert ZZ['x'].from_FractionField(y, ZZ['x']) == DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)

def test_FractionField__init():
    raises(GeneratorsNeeded, lambda: ZZ.frac_field())

def test_inject():
    assert ZZ.inject(x, y, z) == ZZ[x, y, z]
    assert ZZ[x].inject(y, z) == ZZ[x, y, z]
    assert ZZ.frac_field(x).inject(y, z) == ZZ.frac_field(x, y, z)
    raises(GeneratorsError, lambda: ZZ[x].inject(x))

def test_Domain_map():
    seq = ZZ.map([1, 2, 3, 4])

    assert all(ZZ.of_type(elt) for elt in seq)

    seq = ZZ.map([[1, 2, 3, 4]])

    assert all(ZZ.of_type(elt) for elt in seq[0]) and len(seq) == 1

def test_Domain___eq__():
    assert (ZZ[x,y] == ZZ[x,y]) == True
    assert (QQ[x,y] == QQ[x,y]) == True

    assert (ZZ[x,y] == QQ[x,y]) == False
    assert (QQ[x,y] == ZZ[x,y]) == False

    assert (ZZ.frac_field(x,y) == ZZ.frac_field(x,y)) == True
    assert (QQ.frac_field(x,y) == QQ.frac_field(x,y)) == True

    assert (ZZ.frac_field(x,y) == QQ.frac_field(x,y)) == False
    assert (QQ.frac_field(x,y) == ZZ.frac_field(x,y)) == False

def test_Domain__algebraic_field():
    alg = ZZ.algebraic_field(sqrt(2))
    assert alg.ext.minpoly == Poly(x**2 - 2)
    assert alg.dom == QQ

    alg = QQ.algebraic_field(sqrt(2))
    assert alg.ext.minpoly == Poly(x**2 - 2)
    assert alg.dom == QQ

    alg = alg.algebraic_field(sqrt(3))
    assert alg.ext.minpoly == Poly(x**4 - 10*x**2 + 1)
    assert alg.dom == QQ

def test_PolynomialRing__from_FractionField():
    f = DMF(([1, 0, 1], [1, 1]), ZZ)
    g = DMF(([1, 0, 1], [1]), ZZ)

    assert ZZ[x].from_FractionField(f, ZZ[x]) is None
    assert ZZ[x].from_FractionField(g, ZZ[x]) == DMP([ZZ(1), ZZ(0), ZZ(1)], ZZ)

def test_PythonRationalType__init__():
    assert Q(0).p == 0
    assert Q(0).q == 1
    assert Q(0, 1).p == 0
    assert Q(0, 1).q == 1
    assert Q(0,-1).p == 0
    assert Q(0,-1).q == 1

    assert Q(1).p == 1
    assert Q(1).q == 1
    assert Q(1, 1).p == 1
    assert Q(1, 1).q == 1
    assert Q(-1,-1).p == 1
    assert Q(-1,-1).q == 1

    assert Q(-1).p == -1
    assert Q(-1).q ==  1
    assert Q(-1, 1).p == -1
    assert Q(-1, 1).q ==  1
    assert Q( 1,-1).p == -1
    assert Q( 1,-1).q ==  1

    assert Q(1, 2).p == 1
    assert Q(1, 2).q == 2
    assert Q(3, 4).p == 3
    assert Q(3, 4).q == 4

    assert Q(2, 2).p == 1
    assert Q(2, 2).q == 1
    assert Q(2, 4).p == 1
    assert Q(2, 4).q == 2

def test_PythonRationalType__hash__():
    assert hash(Q(0)) == hash(0)
    assert hash(Q(1)) == hash(1)
    assert hash(Q(117)) == hash(117)

def test_PythonRationalType__int__():
    assert int(Q(-1, 4)) == 0
    assert int(Q( 1, 4)) == 0
    assert int(Q(-5, 4)) == -1
    assert int(Q( 5, 4)) ==  1

def test_PythonRationalType__float__():
    assert float(Q(-1, 2)) == -0.5
    assert float(Q( 1, 2)) ==  0.5

def test_PythonRationalType__abs__():
    assert abs(Q(-1, 2)) == Q(1, 2)
    assert abs(Q( 1, 2)) == Q(1, 2)

def test_PythonRationalType__pos__():
    assert +Q(-1, 2) == Q(-1, 2)
    assert +Q( 1, 2) == Q( 1, 2)

def test_PythonRationalType__neg__():
    assert -Q(-1, 2) == Q( 1, 2)
    assert -Q( 1, 2) == Q(-1, 2)

def test_PythonRationalType__add__():
    assert Q(-1, 2) + Q( 1, 2) == Q(0)
    assert Q( 1, 2) + Q(-1, 2) == Q(0)

    assert Q(1, 2) + Q(1, 2) == Q(1)
    assert Q(1, 2) + Q(3, 2) == Q(2)
    assert Q(3, 2) + Q(1, 2) == Q(2)
    assert Q(3, 2) + Q(3, 2) == Q(3)

    assert 1 + Q(1, 2) == Q(3, 2)
    assert Q(1, 2) + 1 == Q(3, 2)

def test_PythonRationalType__sub__():
    assert Q(-1, 2) - Q( 1, 2) == Q(-1)
    assert Q( 1, 2) - Q(-1, 2) == Q( 1)

    assert Q(1, 2) - Q(1, 2) == Q( 0)
    assert Q(1, 2) - Q(3, 2) == Q(-1)
    assert Q(3, 2) - Q(1, 2) == Q( 1)
    assert Q(3, 2) - Q(3, 2) == Q( 0)

    assert 1 - Q(1, 2) == Q( 1, 2)
    assert Q(1, 2) - 1 == Q(-1, 2)

def test_PythonRationalType__mul__():
    assert Q(-1, 2) * Q( 1, 2) == Q(-1, 4)
    assert Q( 1, 2) * Q(-1, 2) == Q(-1, 4)

    assert Q(1, 2) * Q(1, 2) == Q(1, 4)
    assert Q(1, 2) * Q(3, 2) == Q(3, 4)
    assert Q(3, 2) * Q(1, 2) == Q(3, 4)
    assert Q(3, 2) * Q(3, 2) == Q(9, 4)

    assert 2 * Q(1, 2) == Q(1)
    assert Q(1, 2) * 2 == Q(1)

def test_PythonRationalType__div__():
    assert Q(-1, 2) / Q( 1, 2) == Q(-1)
    assert Q( 1, 2) / Q(-1, 2) == Q(-1)

    assert Q(1, 2) / Q(1, 2) == Q(1)
    assert Q(1, 2) / Q(3, 2) == Q(1, 3)
    assert Q(3, 2) / Q(1, 2) == Q(3)
    assert Q(3, 2) / Q(3, 2) == Q(1)

    assert 2 / Q(1, 2) == Q(4)
    assert Q(1, 2) / 2 == Q(1, 4)

    raises(ZeroDivisionError, lambda: Q(1, 2) / Q(0))
    raises(ZeroDivisionError, lambda: Q(1, 2) / 0)

def test_PythonRationalType__pow__():
    assert Q(1)**10 == Q(1)
    assert Q(2)**10 == Q(1024)

    assert Q(1)**(-10) == Q(1)
    assert Q(2)**(-10) == Q(1, 1024)

def test_PythonRationalType__eq__():
    assert (Q(1, 2) == Q(1, 2)) is True
    assert (Q(1, 2) != Q(1, 2)) is False

    assert (Q(1, 2) == Q(1, 3)) is False
    assert (Q(1, 2) != Q(1, 3)) is True

def test_PythonRationalType__lt_le_gt_ge__():
    assert (Q(1, 2) <  Q(1, 4)) is False
    assert (Q(1, 2) <= Q(1, 4)) is False
    assert (Q(1, 2) >  Q(1, 4)) is True
    assert (Q(1, 2) >= Q(1, 4)) is True

    assert (Q(1, 4) <  Q(1, 2)) is True
    assert (Q(1, 4) <= Q(1, 2)) is True
    assert (Q(1, 4) >  Q(1, 2)) is False
    assert (Q(1, 4) >= Q(1, 2)) is False

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

def test_FF_of_type():
    assert FF(3).of_type(FF(3)(1)) is True
    assert FF(5).of_type(FF(5)(3)) is True
    assert FF(5).of_type(FF(7)(3)) is False

def test___eq__():
    assert not QQ['x'] == ZZ['x']
    assert not QQ.frac_field(x) == ZZ.frac_field(x)

def test_RealDomain_from_sympy():
    RR = RR_mpmath()

    assert RR.convert(S(0)) == RR.dtype(0)
    assert RR.convert(S(0.0)) == RR.dtype(0.0)
    assert RR.convert(S(1)) == RR.dtype(1)
    assert RR.convert(S(1.0)) == RR.dtype(1.0)
    assert RR.convert(sin(1)) == RR.dtype(sin(1).evalf())
    raises(CoercionFailed, lambda: RR.convert(x))
    raises(CoercionFailed, lambda: RR.convert(oo))
    raises(CoercionFailed, lambda: RR.convert(-oo))

    RR = RR_sympy()

    assert RR.convert(S(0)) == RR.dtype(0)
    assert RR.convert(S(0.0)) == RR.dtype(0.0)
    assert RR.convert(S(1)) == RR.dtype(1)
    assert RR.convert(S(1.0)) == RR.dtype(1.0)
    assert RR.convert(sin(1)) == RR.dtype(sin(1).evalf())
    assert RR.n(3, 2) == RR.evalf(3, 2) == Rational(3).n(2)
    raises(CoercionFailed, lambda: RR.convert(x))
    raises(CoercionFailed, lambda: RR.convert(oo))
    raises(CoercionFailed, lambda: RR.convert(-oo))

def test_ModularInteger():
    GF = ModularIntegerFactory(3)

    a = GF(0)
    assert isinstance(a, GF) and a == 0
    a = GF(1)
    assert isinstance(a, GF) and a == 1
    a = GF(2)
    assert isinstance(a, GF) and a == 2
    a = GF(3)
    assert isinstance(a, GF) and a == 0
    a = GF(4)
    assert isinstance(a, GF) and a == 1

    a = GF(GF(0))
    assert isinstance(a, GF) and a == 0
    a = GF(GF(1))
    assert isinstance(a, GF) and a == 1
    a = GF(GF(2))
    assert isinstance(a, GF) and a == 2
    a = GF(GF(3))
    assert isinstance(a, GF) and a == 0
    a = GF(GF(4))
    assert isinstance(a, GF) and a == 1

    a = -GF(1)
    assert isinstance(a, GF) and a == 2
    a = -GF(2)
    assert isinstance(a, GF) and a == 1

    a = 2 + GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2) + 2
    assert isinstance(a, GF) and a == 1
    a = GF(2) + GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2) + GF(2)
    assert isinstance(a, GF) and a == 1

    a = 3 - GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(3) - 2
    assert isinstance(a, GF) and a == 1
    a = GF(3) - GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(3) - GF(2)
    assert isinstance(a, GF) and a == 1

    a = 2*GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2)*2
    assert isinstance(a, GF) and a == 1
    a = GF(2)*GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2)*GF(2)
    assert isinstance(a, GF) and a == 1

    a = 2/GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2)/2
    assert isinstance(a, GF) and a == 1
    a = GF(2)/GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(2)/GF(2)
    assert isinstance(a, GF) and a == 1

    a = 1 % GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(1) % 2
    assert isinstance(a, GF) and a == 1
    a = GF(1) % GF(2)
    assert isinstance(a, GF) and a == 1
    a = GF(1) % GF(2)
    assert isinstance(a, GF) and a == 1

    a = GF(2)**0
    assert isinstance(a, GF) and a == 1
    a = GF(2)**1
    assert isinstance(a, GF) and a == 2
    a = GF(2)**2
    assert isinstance(a, GF) and a == 1

    assert bool(GF(3)) is False
    assert bool(GF(4)) is True

    GF = ModularIntegerFactory(5)

    a = GF(1)**(-1)
    assert isinstance(a, GF) and a == 1
    a = GF(2)**(-1)
    assert isinstance(a, GF) and a == 3
    a = GF(3)**(-1)
    assert isinstance(a, GF) and a == 2
    a = GF(4)**(-1)
    assert isinstance(a, GF) and a == 4

    assert (GF(1) <  GF(2)) == True
    assert (GF(1) <= GF(2)) == True
    assert (GF(1) >  GF(2)) == False
    assert (GF(1) >= GF(2)) == False

    assert (GF(3) <  GF(2)) == False
    assert (GF(3) <= GF(2)) == False
    assert (GF(3) >  GF(2)) == True
    assert (GF(3) >= GF(2)) == True

    assert (GF(1) <  GF(7)) == True
    assert (GF(1) <= GF(7)) == True
    assert (GF(1) >  GF(7)) == False
    assert (GF(1) >= GF(7)) == False

    assert (GF(3) <  GF(7)) == False
    assert (GF(3) <= GF(7)) == False
    assert (GF(3) >  GF(7)) == True
    assert (GF(3) >= GF(7)) == True

    assert (GF(1) <  2) == True
    assert (GF(1) <= 2) == True
    assert (GF(1) >  2) == False
    assert (GF(1) >= 2) == False

    assert (GF(3) <  2) == False
    assert (GF(3) <= 2) == False
    assert (GF(3) >  2) == True
    assert (GF(3) >= 2) == True

    assert (GF(1) <  7) == True
    assert (GF(1) <= 7) == True
    assert (GF(1) >  7) == False
    assert (GF(1) >= 7) == False

    assert (GF(3) <  7) == False
    assert (GF(3) <= 7) == False
    assert (GF(3) >  7) == True
    assert (GF(3) >= 7) == True

    raises(NotInvertible, lambda: GF(0)**(-1))
    raises(NotInvertible, lambda: GF(5)**(-1))

    raises(ValueError, lambda: ModularIntegerFactory(0))
    raises(ValueError, lambda: ModularIntegerFactory(2.1))
    raises(TypeError, lambda: ModularIntegerFactory(3, QQ))
