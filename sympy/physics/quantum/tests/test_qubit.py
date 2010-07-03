from sympy.physics.quantum.qubit import *
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.qft import *
from sympy.physics.quantum.represent import *
from sympy.physics.quantum.qapply import *
from sympy import symbols, Rational
from sympy.core.numbers import *
from sympy.functions.elementary import *
from sympy.physics.quantum.shor import *
from sympy.core.containers import Tuple
from sympy.matrices.matrices import Matrix
import random
x, y = symbols('x,y')

epsilon = .000001

def test_Qubit():
    array = [0,0,1,1,0]
    qb = Qubit('00110')
    assert qb.flip(0) == Qubit('00111')
    assert qb.flip(1) == Qubit('00100')
    assert qb.flip(4) == Qubit('10110')
    assert qb.dimension == 5
    for i in range(5):
        assert qb[i] == array[4-i]
    assert len(qb) == 5
    qb = Qubit('110')

def test_QubitBra():
    assert Qubit(0).dual_class == QubitBra
    assert QubitBra(0).dual_class == Qubit
    assert represent(Qubit(1,1,0), nqubits=3).H ==\
           represent(QubitBra(1,1,0), nqubits=3)
    assert Qubit(0,1)._eval_innerproduct_QubitBra(QubitBra(1,0)) == Integer(0)
    assert Qubit(0,1)._eval_innerproduct_QubitBra(QubitBra(0,1)) == Integer(1)

def test_IntQubit():
    assert IntQubit(8).as_int() == 8
    assert IntQubit(8).qubit_values == (1,0,0,0)
    assert IntQubit(7, 4).qubit_values == (0,1,1,1)
    assert IntQubit(3) == IntQubit(3,2)

    #test Dual Classes
    assert IntQubit(3).dual_class == IntQubitBra
    assert IntQubitBra(3).dual_class == IntQubit

def test_superposition_of_states():
    assert qapply(CNOT(0,1)*HadamardGate(0)*(1/sqrt(2)*Qubit('01') + 1/sqrt(2)*Qubit('10'))).expand() == (Qubit('01')/2 + Qubit('00')/2 - Qubit('11')/2 +\
     Qubit('10')/2)

    assert matrix_to_qubit(represent(CNOT(0,1)*HadamardGate(0)\
    *(1/sqrt(2)*Qubit('01') + 1/sqrt(2)*Qubit('10')), nqubits=2))\
     == (Qubit('01')/2 + Qubit('00')/2 - Qubit('11')/2 + Qubit('10')/2)


#test apply methods
def test_apply_represent_equality():
    gates = [HadamardGate(int(3*random.random())),\
     XGate(int(3*random.random())), ZGate(int(3*random.random())),\
      YGate(int(3*random.random())), ZGate(int(3*random.random())),\
       PhaseGate(int(3*random.random()))]

    circuit = Qubit(int(random.random()*2),int(random.random()*2),\
    int(random.random()*2),int(random.random()*2),int(random.random()*2),\
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
    assert matrix_to_qubit(Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))\
    == Qubit(0,0,0,0)
    assert qubit_to_matrix(Qubit(0,0,0,0)) ==\
    Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    assert matrix_to_qubit(sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])) ==\
    (2*sqrt(2)*(Qubit(0,0,0) + Qubit(0,0,1) + Qubit(0,1,0) + Qubit(0,1,1)\
    + Qubit(1,0,0) + Qubit(1,0,1) + Qubit(1,1,0) + Qubit(1,1,1))).expand()
    assert qubit_to_matrix(2*sqrt(2)*(Qubit(0,0,0) + Qubit(0,0,1) + Qubit(0,1,0)\
    + Qubit(0,1,1) + Qubit(1,0,0) + Qubit(1,0,1) + Qubit(1,1,0) + Qubit(1,1,1)))\
    == sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])

def test_measure_partial():
    #Basic test of collapse of entangled two qubits (Bell States)
    state = Qubit('01') + Qubit('10')
    assert sorted(measure_partial(state, (0,))) ==\
           [(Qubit('01'), Rational(1,2)), (Qubit('10'), Rational(1,2))]
    assert sorted(measure_partial(state, (0,))) ==\
           sorted(measure_partial(state, (1,)))

    #Test of more complex collapse and probability calculation
    state1 = sqrt(2)/sqrt(3)*Qubit('00001') + 1/sqrt(3)*Qubit('11111')
    assert measure_partial(state1, (0,)) ==\
           [(sqrt(2)/sqrt(3)*Qubit('00001') + 1/sqrt(3)*Qubit('11111'), 1)]
    assert measure_partial(state1, (1,2)) == measure_partial(state1, (3,4))
    assert measure_partial(state1, (1,2,3)) ==\
           [(Qubit('00001'), Rational(2,3)), (Qubit('11111'), Rational(1,3))]

    #test of measuring multiple bits at once
    state2 = Qubit('1111') + Qubit('1101') + Qubit('1011') + Qubit('1000')
    assert sorted(measure_partial(state2, (0,1,3))) ==\
           sorted([(Qubit('1011')/sqrt(2) + Qubit('1111')/sqrt(2), Rational(1,2)), \
           (Qubit('1101'), Rational(1,4)), (Qubit('1000'), Rational(1,4))])
    assert sorted(measure_partial(state2, (0,))) ==\
           sorted([(Qubit('1111')/sqrt(3) + Qubit('1101')/sqrt(3) + Qubit('1011')/sqrt(3), Rational(3,4)),\
           (Qubit('1000'), Rational(1,4))])


def test_measure_all():
    assert measure_all(Qubit('11')) == [(Qubit('11'), 1)]
    state = Qubit('11') + Qubit('10')
    assert sorted(measure_all(state)) == sorted([(Qubit('11'), Rational(1,2)),\
           (Qubit('10'), Rational(1,2))])
    state2 = Qubit('11')/sqrt(5) + 2*Qubit('00')/sqrt(5)
    assert sorted(measure_all(state2)) == sorted([(Qubit('11'), Rational(1,5)), \
           (Qubit('00'), Rational(4,5))])
