from sympy import I
from sympy.physics.quantum import Dagger, Commutator, qapply
from sympy.physics.quantum.pauli import SigmaOpBase, SigmaX, SigmaY, SigmaZ
from sympy.physics.quantum.pauli import SigmaZKet, SigmaZBra


def test_pauli_operators():
    sx, sy, sz = SigmaX(), SigmaY(), SigmaZ()

    assert isinstance(sx, SigmaOpBase)
    assert isinstance(sy, SigmaOpBase)
    assert isinstance(sz, SigmaOpBase)

    assert Commutator(sx, sy).doit() == 2 * I * sz
    assert Commutator(sy, sz).doit() == 2 * I * sx
    assert Commutator(sz, sx).doit() == 2 * I * sy

    assert Dagger(sx) == sx
    assert Dagger(sy) == sy
    assert Dagger(sz) == sz


def test_pauli_states():
    sx, sz = SigmaX(), SigmaZ()

    up = SigmaZKet(0)
    down = SigmaZKet(1)

    assert qapply(sx * up) == down
    assert qapply(sx * down) == up
    assert qapply(sz * up) == up
    assert qapply(sz * down) == - down

    up = SigmaZBra(0)
    down = SigmaZBra(1)

    assert qapply(up * sx, dagger=True) == down
    assert qapply(down * sx, dagger=True) == up
    assert qapply(up * sz, dagger=True) == up
    assert qapply(down * sz, dagger=True) == - down

    assert Dagger(SigmaZKet(0)) == SigmaZBra(0)
    assert Dagger(SigmaZBra(1)) == SigmaZKet(1)