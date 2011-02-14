from sympy import sqrt, pprint
from sympy.physics.quantum import *
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.qubit import *

psi = (1/(sqrt(2)))*(Qubit('00')+Qubit('11'))

# Dense coding demo:

# Assume Alice has the left QBit in psi
# pprint(psi)

# The corresponding gates applied to Alice's QBit are:
# Identity Gate (1), Not Gate (X), Z Gate (Z), Z Gate and Not Gate (ZX)
# Then there's the controlled not gate (with Alice's as control):CNOT(1, 0)
# And the Hadamard gate applied to Alice's Qbit: H(1)

# To Send Bob the message |0>|0>
circuit = H(1)*CNOT(1,0)*psi
result = apply_operators(circuit)
pprint(result)

# To send Bob the message |0>|1>
circuit = H(1)*CNOT(1,0)*X(1)*psi
result = apply_operators(circuit)
pprint(result)

# To send Bob the message |1>|0>
circuit = H(1)*CNOT(1,0)*Z(1)*psi
result = apply_operators(circuit)
pprint(result)

# To send Bob the message |1>|1>
circuit = H(1)*CNOT(1,0)*Z(1)*X(1)*psi
result = apply_operators(circuit)
pprint(result)
