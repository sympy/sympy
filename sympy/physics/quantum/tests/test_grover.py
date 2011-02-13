from sympy import sqrt, power
from sympy.physics.quantum.applyops import apply_operators
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.grover import *
import sympy.physics.quantum.grover

gpath = sympy.physics.quantum.grover

def return_one_on_two(qubits):
    qubitTwo = gpath._create_zeroes(qubits.nqubits - 2) + '10'
    return 1 if qubits == Qubit(qubitTwo) else 0

def return_one_on_one(qubits):
    qubitOne = gpath._create_zeroes(qubits.nqubits - 2) + '01'
    return 1 if qubits == Qubit(qubitOne) else 0

def test_create_basis_states():
    first_half_state = Qubit('00')/2 + Qubit('01')/2
    second_half_state = Qubit('10')/2 + Qubit('11')/2
    super_state = first_half_state + second_half_state
    assert super_state == gpath._create_basis_states(2)

    first_q = (1/sqrt(8))*Qubit('000') + (1/sqrt(8))*Qubit('001')
    second_q = (1/sqrt(8))*Qubit('010') + (1/sqrt(8))*Qubit('011')  
    third_q = (1/sqrt(8))*Qubit('100') + (1/sqrt(8))*Qubit('101')  
    fourth_q = (1/sqrt(8))*Qubit('110') + (1/sqrt(8))*Qubit('111')  
    super_state = first_q + second_q + third_q + fourth_q
    assert super_state == gpath._create_basis_states(3)

def test_OracleGate():
    v = OracleGate(1, lambda qubits: 1 if qubits == Qubit('0') else 0)
    assert apply_operators(v*Qubit('0')) == -Qubit('0')
    assert apply_operators(v*Qubit('1')) == Qubit('1')

    v = OracleGate(2, return_one_on_two)
    assert apply_operators(v*Qubit('00')) == Qubit('00')
    assert apply_operators(v*Qubit('01')) == Qubit('01')  
    assert apply_operators(v*Qubit('10')) == -Qubit('10')
    assert apply_operators(v*Qubit('11')) == Qubit('11')

def test_WGate():
    numqubits = 2
    basis_states = gpath._create_basis_states(numqubits)
    w = WGate(numqubits)
    assert apply_operators(w*basis_states) == basis_states

    qubit_one = Qubit('01')
    expected = ((2/sqrt(pow(2, numqubits)))*basis_states) - qubit_one
    assert apply_operators(w*qubit_one) == expected

def test_grover_iteration_1():
    numqubits = 2
    basis_states = gpath._create_basis_states(numqubits)
    v = OracleGate(numqubits, return_one_on_one)
    iterated = gpath.grover_iteration(basis_states, v)
    expected = Qubit('01')
    assert apply_operators(iterated) == expected

def test_grover_iteration_2():
    numqubits = 4
    basis_states = gpath._create_basis_states(numqubits)
    v = OracleGate(numqubits, return_one_on_two)
    # Should return Qubit('0010') after (pi/4)sqrt(pow(2, n))
    # In this case, after around pi times (3 or 4)
    #print ''
    iterated = gpath.grover_iteration(basis_states, v) 
    iterated = apply_operators(iterated)
    #print iterated
    iterated = gpath.grover_iteration(iterated, v) 
    iterated = apply_operators(iterated)
    #print iterated
    iterated = gpath.grover_iteration(iterated, v) 
    iterated = apply_operators(iterated)
    #print iterated
    # In this case, probability was highest after 3 iterations
    # Probability of Qubit('0010') was 251/256 (3) vs 781/1024 (4)
    # Ask about measurement
    assert True

def test_grover():
    # Needs implementation
    # Ask about measurement
    assert True
