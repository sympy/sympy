import random

from sympy import Integer, Matrix, Rational, sqrt, symbols
from sympy.physics.quantum.qubit import (measure_all, measure_partial,
                                         matrix_to_qubit, matrix_to_density,
                                         qubit_to_matrix, IntQubit,
                                         IntQubitBra, QubitBra)
from sympy.physics.quantum.gate import (HadamardGate, CNOT, XGate, YGate,
                                        ZGate, PhaseGate)
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.shor import Qubit
from sympy.utilities.pytest import raises
from sympy.physics.quantum.density import Density
from sympy.core.trace import Tr

x, y = symbols('x,y')

epsilon = .000001


def test_Qubit():
    array = [0, 0, 1, 1, 0]
    qb = Qubit('00110')
    assert qb.flip(0) == Qubit('00111')
    assert qb.flip(1) == Qubit('00100')
    assert qb.flip(4) == Qubit('10110')
    assert qb.dimension == 5
    for i in range(5):
        assert qb[i] == array[4 - i]
    assert len(qb) == 5
    qb = Qubit('110')


def test_QubitBra():
    assert Qubit(0).dual_class() == QubitBra
    assert QubitBra(0).dual_class() == Qubit
    assert represent(Qubit(1, 1, 0), nqubits=3).H == \
        represent(QubitBra(1, 1, 0), nqubits=3)
    assert Qubit(
        0, 1)._eval_innerproduct_QubitBra(QubitBra(1, 0)) == Integer(0)
    assert Qubit(
        0, 1)._eval_innerproduct_QubitBra(QubitBra(0, 1)) == Integer(1)


def test_IntQubit():
    assert IntQubit(8).as_int() == 8
    assert IntQubit(8).qubit_values == (1, 0, 0, 0)
    assert IntQubit(7, 4).qubit_values == (0, 1, 1, 1)
    assert IntQubit(3) == IntQubit(3, 2)

    #test Dual Classes
    assert IntQubit(3).dual_class() == IntQubitBra
    assert IntQubitBra(3).dual_class() == IntQubit

    assert IntQubit(5)._eval_innerproduct_IntQubitBra(IntQubitBra(5)) == Integer(1)
    assert IntQubit(4)._eval_innerproduct_IntQubitBra(IntQubitBra(5)) == Integer(0)
    raises(ValueError, lambda: IntQubit(4, 1))


def test_superposition_of_states():
    assert qapply(CNOT(0, 1)*HadamardGate(0)*(1/sqrt(2)*Qubit('01') + 1/sqrt(2)*Qubit('10'))).expand() == (Qubit('01')/2 + Qubit('00')/2 - Qubit('11')/2 +
     Qubit('10')/2)

    assert matrix_to_qubit(represent(CNOT(0, 1)*HadamardGate(0)
        *(1/sqrt(2)*Qubit('01') + 1/sqrt(2)*Qubit('10')), nqubits=2)) == \
        (Qubit('01')/2 + Qubit('00')/2 - Qubit('11')/2 + Qubit('10')/2)


#test apply methods
def test_apply_represent_equality():
    gates = [HadamardGate(int(3*random.random())),
     XGate(int(3*random.random())), ZGate(int(3*random.random())),
        YGate(int(3*random.random())), ZGate(int(3*random.random())),
        PhaseGate(int(3*random.random()))]

    circuit = Qubit(int(random.random()*2), int(random.random()*2),
    int(random.random()*2), int(random.random()*2), int(random.random()*2),
        int(random.random()*2))
    for i in range(int(random.random()*6)):
        circuit = gates[int(random.random()*6)]*circuit

    mat = represent(circuit, nqubits=6)
    states = qapply(circuit)
    state_rep = matrix_to_qubit(mat)
    states = states.expand()
    state_rep = state_rep.expand()
    assert state_rep == states


def test_matrix_to_qubits():
    assert matrix_to_qubit(
        Matrix([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])) == \
        Qubit(0, 0, 0, 0)
    assert qubit_to_matrix(Qubit(0, 0, 0, 0)) == \
        Matrix([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert matrix_to_qubit(sqrt(2)*2*Matrix([1, 1, 1, 1, 1, 1, 1, 1])) == \
        (2*sqrt(2)*(Qubit(0, 0, 0) + Qubit(0, 0, 1) + Qubit(0, 1, 0) +
        Qubit(0, 1, 1) + Qubit(1, 0, 0) + Qubit(1, 0, 1) + Qubit(1, 1, 0) +
         Qubit(1, 1, 1))).expand()
    assert qubit_to_matrix(2*sqrt(2)*(Qubit(0, 0, 0) + Qubit(0, 0, 1) +
        Qubit(0, 1, 0) + Qubit(0, 1, 1) + Qubit(1, 0, 0) + Qubit(1, 0, 1) +
        Qubit(1, 1, 0) + Qubit(1, 1, 1))) == \
        sqrt(2)*2*Matrix([1, 1, 1, 1, 1, 1, 1, 1])


def test_measure_normalize():
    a, b = symbols('a b')
    state = a*Qubit('110') + b*Qubit('111')
    assert measure_partial(state, (0,), normalize=False) == \
        [(a*Qubit('110'), a*a.conjugate()), (b*Qubit('111'), b*b.conjugate())]
    assert measure_all(state, normalize=False) == \
        [(Qubit('110'), a*a.conjugate()), (Qubit('111'), b*b.conjugate())]


def test_measure_partial():
    #Basic test of collapse of entangled two qubits (Bell States)
    state = Qubit('01') + Qubit('10')
    assert measure_partial(state, (0,)) == \
        [(Qubit('10'), Rational(1, 2)), (Qubit('01'), Rational(1, 2))]
    assert measure_partial(state, (0,)) == \
        measure_partial(state, (1,))[::-1]

    #Test of more complex collapse and probability calculation
    state1 = sqrt(2)/sqrt(3)*Qubit('00001') + 1/sqrt(3)*Qubit('11111')
    assert measure_partial(state1, (0,)) == \
        [(sqrt(2)/sqrt(3)*Qubit('00001') + 1/sqrt(3)*Qubit('11111'), 1)]
    assert measure_partial(state1, (1, 2)) == measure_partial(state1, (3, 4))
    assert measure_partial(state1, (1, 2, 3)) == \
        [(Qubit('00001'), Rational(2, 3)), (Qubit('11111'), Rational(1, 3))]

    #test of measuring multiple bits at once
    state2 = Qubit('1111') + Qubit('1101') + Qubit('1011') + Qubit('1000')
    assert measure_partial(state2, (0, 1, 3)) == \
        [(Qubit('1000'), Rational(1, 4)), (Qubit('1101'), Rational(1, 4)),
         (Qubit('1011')/sqrt(2) + Qubit('1111')/sqrt(2), Rational(1, 2))]
    assert measure_partial(state2, (0,)) == \
        [(Qubit('1000'), Rational(1, 4)),
         (Qubit('1111')/sqrt(3) + Qubit('1101')/sqrt(3) +
          Qubit('1011')/sqrt(3), Rational(3, 4))]


def test_measure_all():
    assert measure_all(Qubit('11')) == [(Qubit('11'), 1)]
    state = Qubit('11') + Qubit('10')
    assert measure_all(state) == [(Qubit('10'), Rational(1, 2)),
           (Qubit('11'), Rational(1, 2))]
    state2 = Qubit('11')/sqrt(5) + 2*Qubit('00')/sqrt(5)
    assert measure_all(state2) == \
        [(Qubit('00'), Rational(4, 5)), (Qubit('11'), Rational(1, 5))]


def test_eval_trace():
    q1 = Qubit('10110')
    q2 = Qubit('01010')
    d = Density([q1, 0.6], [q2, 0.4])

    t = Tr(d)
    assert t.doit() == 1

    # extreme bits
    t = Tr(d, 0)
    assert t.doit() == (0.4*Density([Qubit('0101'), 1]) +
                        0.6*Density([Qubit('1011'), 1]))
    t = Tr(d, 4)
    assert t.doit() == (0.4*Density([Qubit('1010'), 1]) +
                        0.6*Density([Qubit('0110'), 1]))
    # index somewhere in between
    t = Tr(d, 2)
    assert t.doit() == (0.4*Density([Qubit('0110'), 1]) +
                        0.6*Density([Qubit('1010'), 1]))
    #trace all indices
    t = Tr(d, [0, 1, 2, 3, 4])
    assert t.doit() == 1

    # trace some indices, initialized in
    # non-canonical order
    t = Tr(d, [2, 1, 3])
    assert t.doit() == (0.4*Density([Qubit('00'), 1]) +
                        0.6*Density([Qubit('10'), 1]))

    # mixed states
    q = (1/sqrt(2)) * (Qubit('00') + Qubit('11'))
    d = Density( [q, 1.0] )
    t = Tr(d, 0)
    assert t.doit() == (0.5*Density([Qubit('0'), 1]) +
                        0.5*Density([Qubit('1'), 1]))


def test_matrix_to_density():
    mat = Matrix([[0, 0], [0, 1]])
    assert matrix_to_density(mat) == Density([Qubit('1'), 1])

    mat = Matrix([[1, 0], [0, 0]])
    assert matrix_to_density(mat) == Density([Qubit('0'), 1])

    mat = Matrix([[0, 0], [0, 0]])
    assert matrix_to_density(mat) == 0

    mat = Matrix([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 0]])

    assert matrix_to_density(mat) == Density([Qubit('10'), 1])

    mat = Matrix([[1, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]])

    assert matrix_to_density(mat) == Density([Qubit('00'), 1])
