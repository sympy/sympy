from __future__ import annotations
from sympy import Integer
from sympy.sets.ordinals import Ordinal, OrdinalOmega, OrdinalZero, OmegaPower, ord0, omega
from sympy.testing.pytest import raises


def test_hash_finite():
    for value in (0, 1, 2, 5, 2**70):
        ordinal = Ordinal.convert(value)
        sympy_integer = Integer(value)

        assert ordinal == value
        assert ordinal == sympy_integer
        assert hash(ordinal) == hash(value)
        assert hash(ordinal) == hash(sympy_integer)

        assert len({value, ordinal}) == 1
        assert {value: 'value'}[ordinal] == 'value'
        assert {ordinal: 'ordinal'}[value] == 'ordinal'

        if value:
            term = OmegaPower(0, value)
            assert term == value
            assert hash(term) == hash(value)
            assert term == sympy_integer
            assert hash(term) == hash(sympy_integer)
            assert len({value, term}) == 1
            assert {value: 'value'}[term] == 'value'
            assert {term: 'term'}[value] == 'term'


def test_hash_omega():
    zero = OrdinalZero()
    assert zero == ord0
    assert hash(zero) == hash(ord0)
    assert hash(zero) == hash(0)
    assert len({zero, ord0}) == 1
    assert {zero: 'zero'}[ord0] == 'zero'

    ordinal_omega = OrdinalOmega()
    omega_representation = Ordinal(OmegaPower(1, 1))
    assert ordinal_omega == omega_representation
    assert hash(ordinal_omega) == hash(omega_representation)
    assert len({ordinal_omega, omega_representation}) == 1


def test_hash_nested():
    ordinal = Ordinal(OmegaPower(omega, 1))
    equivalent_ordinal = Ordinal(OmegaPower(Ordinal(OmegaPower(1, 1)), 1))
    assert ordinal == equivalent_ordinal
    assert hash(ordinal) == hash(equivalent_ordinal)

    multi_term_exponent = omega + 1
    equivalent_multi_term_exponent = Ordinal(
        OmegaPower(1, 1), OmegaPower(0, 1))
    ordinal = Ordinal(OmegaPower(multi_term_exponent, 1), OmegaPower(0, 2))
    equivalent_ordinal = Ordinal(
        OmegaPower(equivalent_multi_term_exponent, 1), OmegaPower(0, 2))
    assert ordinal == equivalent_ordinal
    assert hash(ordinal) == hash(equivalent_ordinal)

def test_addition_with_integers():
    assert 3 + Ordinal(OmegaPower(5, 3)) == Ordinal(OmegaPower(5, 3))
    assert Ordinal(OmegaPower(5, 3))+3 == Ordinal(OmegaPower(5, 3), OmegaPower(0, 3))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(0, 2))+3 == \
        Ordinal(OmegaPower(5, 3), OmegaPower(0, 5))


def test_addition_with_ordinals():
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) + Ordinal(OmegaPower(3, 3)) == \
        Ordinal(OmegaPower(5, 3), OmegaPower(3, 5))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) + Ordinal(OmegaPower(4, 2)) == \
        Ordinal(OmegaPower(5, 3), OmegaPower(4, 2))
    assert Ordinal(OmegaPower(omega, 2), OmegaPower(3, 2)) + Ordinal(OmegaPower(4, 2)) == \
        Ordinal(OmegaPower(omega, 2), OmegaPower(4, 2))

def test_comparison():
    assert Ordinal(OmegaPower(5, 3)) > Ordinal(OmegaPower(4, 3), OmegaPower(2, 1))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) < Ordinal(OmegaPower(5, 4))
    assert Ordinal(OmegaPower(5, 4)) < Ordinal(OmegaPower(5, 5), OmegaPower(4, 1))

    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == \
        Ordinal(OmegaPower(5, 3), OmegaPower(3, 2))
    assert not Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == Ordinal(OmegaPower(5, 3))
    assert Ordinal(OmegaPower(omega, 3)) > Ordinal(OmegaPower(5, 3))

def test_multiplication_with_integers():
    w = omega
    assert 3*w == w
    assert w*9 == Ordinal(OmegaPower(1, 9))

def test_multiplication():
    w = omega
    assert w*(w + 1) == w*w + w
    assert (w + 1)*(w + 1) ==  w*w + w + 1
    assert w*1 == w
    assert 1*w == w
    assert w*ord0 == ord0
    assert ord0*w == ord0
    assert w**w == w * w**w
    assert (w**w)*w*w == w**(w + 2)

def test_exponentiation():
    w = omega
    assert ord0**ord0 == Ordinal(OmegaPower(0, 1))
    assert 1**w == Ordinal(OmegaPower(0, 1))
    assert 2**(w**w+w**2+w*3+4) == w**w**w * w**w * w**3 * 16
    assert (w**w*2+1)**w == w**(w**2)
    assert (w**2*3+2)**3==(w**2*3+2)*(w**2*3+2)*(w**2*3+2)
    assert w**2 == w*w
    assert w**3 == w*w*w
    assert w**(w + 1) == Ordinal(OmegaPower(omega + 1, 1))
    assert (w**w)*(w**w) == w**(w*2)

def test_comapre_not_instance():
    w = OmegaPower(omega + 1, 1)
    assert(not (w is None))
    assert(not (w < 5))
    raises(TypeError, lambda: w < 6.66)

def test_is_successort():
    w = Ordinal(OmegaPower(5, 1))
    assert not w.is_successor_ordinal

def test_atomic_classes():
    raises(TypeError, lambda: OrdinalZero(OmegaPower(0,1)))
    raises(TypeError, lambda: OrdinalOmega(OmegaPower(0,1)))

def test_ordinal():
    assert Ordinal() is ord0

def test_ordinal_canonical():
    assert raises(ValueError, lambda: Ordinal(OmegaPower(0, 3), OmegaPower(0, 2)))
    assert raises(ValueError, lambda: Ordinal(OmegaPower(0, 1), OmegaPower(1, 1)))
