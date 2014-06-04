from sympy.physics.quantum import Dagger, AntiCommutator, qapply
from sympy.physics.quantum.fermion import FermionOperator
from sympy.physics.quantum.fermion import FermionFockKet, FermionFockBra


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


def test_fermion_states():
    c = FermionOperator("c")

    # Fock states
    assert (FermionFockBra(0) * FermionFockKet(1)).doit() == 0
    assert (FermionFockBra(1) * FermionFockKet(1)).doit() == 1

    assert qapply(c * FermionFockKet(1)) == FermionFockKet(0)
    assert qapply(c * FermionFockKet(0)) == 0

    assert qapply(Dagger(c) * FermionFockKet(0)) == FermionFockKet(1)
    assert qapply(Dagger(c) * FermionFockKet(1)) == 0