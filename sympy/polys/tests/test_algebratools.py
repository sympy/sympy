"""Tests for classes defining properties of ground domains, e.g. ZZ, QQ, ZZ[x] ... """

from sympy.polys.algebratools import (
    ZZ, QQ, PolynomialRing, FractionField, EX
)

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
)

from sympy.utilities import raises

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

def test_PolynomialRing__init():
    raises(GeneratorsNeeded, "ZZ.poly_ring()")

def test_FractionField__init():
    raises(GeneratorsNeeded, "ZZ.frac_field()")

