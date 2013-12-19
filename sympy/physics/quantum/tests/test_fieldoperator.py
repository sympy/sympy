from sympy.physics.quantum import Dagger, Commutator, AntiCommutator
from sympy.physics.quantum.fieldoperator import BosonOperator, FermionOperator
from sympy.physics.quantum.fieldoperator import (normal_order,
                                                 normal_ordered_form)


def test_bosonoperator():
    a = BosonOperator('a')
    b = BosonOperator('b')

    assert isinstance(a, BosonOperator)
    assert isinstance(Dagger(a), BosonOperator)

    assert a.is_annihilation
    assert not Dagger(a).is_annihilation

    assert BosonOperator("a") == BosonOperator("a")
    assert BosonOperator("a") != BosonOperator("c")
    assert BosonOperator("a", True) != BosonOperator("a", False)

    assert Commutator(a, Dagger(a)).doit() == 1

    assert Commutator(a, Dagger(b)).doit() == a * Dagger(b) - Dagger(b) * a


def test_fermionoperator():
    c = FermionOperator('c')
    d = FermionOperator('d')

    assert isinstance(c, FermionOperator)
    assert isinstance(Dagger(c), FermionOperator)

    assert c.is_annihilation
    assert not Dagger(c).is_annihilation

    assert FermionOperator("c") == FermionOperator("c")
    assert FermionOperator("c") != FermionOperator("d")
    assert FermionOperator("c", True) != FermionOperator("c", False)

    assert AntiCommutator(c, Dagger(c)).doit() == 1

    assert AntiCommutator(c, Dagger(d)).doit() == c * Dagger(d) + Dagger(d) * c


def test_normal_order():
    a = BosonOperator('a')
    b = BosonOperator('b')

    c = FermionOperator('c')
    d = FermionOperator('d')

    assert normal_order(a * Dagger(a)) == Dagger(a) * a
    assert normal_order(Dagger(a) * a) == Dagger(a) * a
    assert normal_order(a * Dagger(a) ** 2) == Dagger(a) ** 2 * a

    assert normal_order(c * Dagger(c)) == - Dagger(c) * c
    assert normal_order(Dagger(c) * c) == Dagger(c) * c
    assert normal_order(c * Dagger(c) ** 2) == Dagger(c) ** 2 * c


def test_normal_ordered_form():
    a = BosonOperator('a')
    b = BosonOperator('b')

    c = FermionOperator('c')
    d = FermionOperator('d')

    assert normal_ordered_form(Dagger(a) * a) == Dagger(a) * a
    assert normal_ordered_form(a * Dagger(a)) == 1 + Dagger(a) * a
    assert normal_ordered_form(a ** 2 * Dagger(a)) == \
        2 * a + Dagger(a) * a ** 2
    assert normal_ordered_form(a ** 3 * Dagger(a)) == \
        3 * a ** 2 + Dagger(a) * a ** 3

    assert normal_ordered_form(Dagger(c) * c) == Dagger(c) * c
    assert normal_ordered_form(c * Dagger(c)) == 1 - Dagger(c) * c
    assert normal_ordered_form(c ** 2 * Dagger(c)) == Dagger(c) * c ** 2
    assert normal_ordered_form(c ** 3 * Dagger(c)) == \
        c ** 2 - Dagger(c) * c ** 3
