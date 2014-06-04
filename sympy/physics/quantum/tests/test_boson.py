from sympy import sqrt, exp, S, prod
from sympy.physics.quantum import Dagger, Commutator, qapply
from sympy.physics.quantum.boson import BosonOperator
from sympy.physics.quantum.boson import (
    BosonFockKet, BosonFockBra, BosonCoherentKet, BosonCoherentBra)


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


def test_boson_states():
    a = BosonOperator("a")

    # Fock states
    n = 3
    assert (BosonFockBra(0) * BosonFockKet(1)).doit() == 0
    assert (BosonFockBra(1) * BosonFockKet(1)).doit() == 1
    assert qapply(BosonFockBra(n) * Dagger(a)**n * BosonFockKet(0)) \
        == sqrt(prod(range(1, n+1)))

    # Coherent states
    alpha1, alpha2 = 1.2, 4.3
    assert (BosonCoherentBra(alpha1) * BosonCoherentKet(alpha1)).doit() == 1
    assert (BosonCoherentBra(alpha2) * BosonCoherentKet(alpha2)).doit() == 1
    assert abs((BosonCoherentBra(alpha1) * BosonCoherentKet(alpha2)).doit() -
               exp(-S(1) / 2 * (alpha1 - alpha2) ** 2)) < 1e-12
    assert qapply(a * BosonCoherentKet(alpha1)) == \
        alpha1 * BosonCoherentKet(alpha1)
