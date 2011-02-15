from sympy import sqrt
from sympy.physics.quantum.applyops import apply_operators
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.grover import *
from sympy.physics.quantum import grover

def return_one_on_two(qubits):
    return True if qubits == IntQubit(2, qubits.nqubits) else False

def return_one_on_one(qubits):
    return True if qubits == IntQubit(1, qubits.nqubits) else False

def test_create_computational_basis():
    nbits = 2
    first_half_state = IntQubit(0, nbits)/2 + IntQubit(1, nbits)/2
    second_half_state = IntQubit(2, nbits)/2 + IntQubit(3, nbits)/2
    super_state = first_half_state + second_half_state
    assert super_state == grover._create_computational_basis(nbits)

    nbits = 3
    first_q = (1/sqrt(8))*IntQubit(0, nbits) + (1/sqrt(8))*IntQubit(1, nbits)
    second_q = (1/sqrt(8))*IntQubit(2, nbits) + (1/sqrt(8))*IntQubit(3, nbits)  
    third_q = (1/sqrt(8))*IntQubit(4, nbits) + (1/sqrt(8))*IntQubit(5, nbits)  
    fourth_q = (1/sqrt(8))*IntQubit(6, nbits) + (1/sqrt(8))*IntQubit(7, nbits)  
    super_state = first_q + second_q + third_q + fourth_q
    assert super_state == grover._create_computational_basis(nbits)

def test_OracleGate():
    v = OracleGate(1, lambda qubits: True if qubits == IntQubit(0) else False)
    assert apply_operators(v*IntQubit(0)) == -IntQubit(0)
    assert apply_operators(v*IntQubit(1)) == IntQubit(1)

    nbits = 2
    v = OracleGate(2, return_one_on_two)
    assert apply_operators(v*IntQubit(0, nbits)) == IntQubit(0, nbits)
    assert apply_operators(v*IntQubit(1, nbits)) == IntQubit(1, nbits)  
    assert apply_operators(v*IntQubit(2, nbits)) == -IntQubit(2, nbits)
    assert apply_operators(v*IntQubit(3, nbits)) == IntQubit(3, nbits)

def test_WGate():
    numqubits = 2
    basis_states = grover._create_computational_basis(numqubits)
    w = WGate(numqubits)
    assert apply_operators(w*basis_states) == basis_states

    qubit_one = IntQubit(1, numqubits)
    expected = ((2/sqrt(pow(2, numqubits)))*basis_states) - qubit_one
    assert apply_operators(w*qubit_one) == expected

def test_grover_iteration_1():
    numqubits = 2
    basis_states = grover._create_computational_basis(numqubits)
    v = OracleGate(numqubits, return_one_on_one)
    iterated = grover.grover_iteration(basis_states, v)
    expected = IntQubit(1, numqubits)
    assert apply_operators(iterated) == expected

def test_grover_iteration_2():
    numqubits = 4
    basis_states = grover._create_computational_basis(numqubits)
    v = OracleGate(numqubits, return_one_on_two)
    # After (pi/4)sqrt(pow(2, n)), IntQubit(2) should have highest prob
    # In this case, after around pi times (3 or 4)
    # print ''
    # print basis_states
    iterated = grover.grover_iteration(basis_states, v) 
    iterated = apply_operators(iterated)
    # print iterated
    iterated = grover.grover_iteration(iterated, v) 
    iterated = apply_operators(iterated)
    # print iterated
    iterated = grover.grover_iteration(iterated, v) 
    iterated = apply_operators(iterated)
    # print iterated
    # In this case, probability was highest after 3 iterations
    # Probability of Qubit('0010') was 251/256 (3) vs 781/1024 (4)
    # Ask about measurement
    expected = (-13*basis_states)/64
    expected = expected + 264*IntQubit(2, numqubits)/256
    assert apply_operators(expected) == iterated

def test_grover():
    nqubits = 2
    actual = grover.apply_grover(return_one_on_one, nqubits)
    assert actual == IntQubit(1, nqubits)

    nqubits = 4
    basis_states = grover._create_computational_basis(nqubits)
    expected = (-13*basis_states)/64
    expected = expected + 264*IntQubit(2, nqubits)/256
    actual = grover.apply_grover(return_one_on_two, 4)
    assert actual == apply_operators(expected)

