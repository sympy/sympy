from sympy import sqrt, pprint
from sympy.physics.quantum import apply_operators
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.circuitplot import circuit_plot
from sympy.physics.quantum.grover import superposition_basis

"""
"""
psi = superposition_basis(2); psi

# Dense coding demo:

# Assume Alice has the left QBit in psi
# pprint(psi)

# The corresponding gates applied to Alice's QBit are:
# Identity Gate (1), Not Gate (X), Z Gate (Z), Z Gate and Not Gate (ZX)
# Then there's the controlled not gate (with Alice's as control):CNOT(1, 0)
# And the Hadamard gate applied to Alice's Qbit: H(1)

# To Send Bob the message |0>|0>
circuit = H(1)*CNOT(1,0)
result = apply_operators(circuit*psi); result
# pprint(result)

# To send Bob the message |0>|1>
circuit = H(1)*CNOT(1,0)*X(1)
result = apply_operators(circuit*psi); result
# pprint(result)

# To send Bob the message |1>|0>
circuit = H(1)*CNOT(1,0)*Z(1)
result = apply_operators(circuit*psi); result
# pprint(result)

# To send Bob the message |1>|1>
circuit = H(1)*CNOT(1,0)*Z(1)*X(1)
result = apply_operators(circuit*psi); result
# pprint(result)

# Using circuit diagrams, the following two are equivalent
circ_classic = CNOT(3,1)*CNOT(2,0)
close('all'); circuit_plot(circ_classic, nqubits=4)
circ_quantum = H(1)*CNOT(1,0)*CGate(3,Z(1))*CNOT(2,1)*CNOT(1,0)*H(1)
close('all'); circuit_plot(circ_quantum, nqubits=4)

# With a 6 qubit example
circ_classic = CNOT(5,2)*CNOT(4,1)*CNOT(3,0)
#circuit_plot(circ_classic, nqubits=6)

"""
Demonstrates dense coding with 2 qubits
Create circuit plots for visual representation of dense coding
"""
def conference_material():
    psi = superposition_basis(2); psi

    # The following sends Bob the message |00>, |01>, |10>, |11> respectively
    circuit = H(1)*CNOT(1,0)
    result = apply_operators(circuit*psi); result

    circuit = H(1)*CNOT(1,0)*X(1)
    result = apply_operators(circuit*psi); result

    circuit = H(1)*CNOT(1,0)*Z(1)
    result = apply_operators(circuit*psi); result

    circuit = H(1)*CNOT(1,0)*Z(1)*X(1)
    result = apply_operators(circuit*psi); result

    # Using circuit diagrams, the following two are equivalent
    circ_classic = CNOT(3,1)*CNOT(2,0)
    close('all'); circuit_plot(circ_classic, nqubits=4)
    circ_quantum = H(1)*CNOT(1,0)*CGate(3,Z(1))*CNOT(2,1)*CNOT(1,0)*H(1)
    close('all'); circuit_plot(circ_quantum, nqubits=4)
   
