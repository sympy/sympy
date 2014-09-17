from sympy.physics.quantum import Dagger, AntiCommutator, Commutator, qapply
from sympy.physics.quantum.pauli import SigmaX, SigmaY, SigmaZ
from sympy.physics.quantum.pauli import SigmaZKet, SigmaZBra


def test_pauli_operators():
    sx, sy, sz = SigmaX(), SigmaY(), SigmaZ()

    assert isinstance(sx, SigmaX)
    assert isinstance(sy, SigmaY)
    assert isinstance(sz, SigmaZ)

    assert Commutator(sx, sy).doit() == 2 * I * sz
    assert Commutator(sy, sz).doit() == 2 * I * sx
    assert Commutator(sz, sx).doit() == 2 * I * sy


def test_pauli_states():
    sx, sy, sz = SigmaX(), SigmaY(), SigmaZ()

    up = SigmaZKet(0)
    down = SigmaZKet(1)

    assert qapply(sx * up) == down
    assert qapply(sx * down) == up
    assert qapply(sz * up) == up
    assert qapply(sz * down) == - down
