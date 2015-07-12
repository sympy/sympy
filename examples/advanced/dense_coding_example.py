#!/usr/bin/env python

"""Demonstration of quantum dense coding."""

from __future__ import division, print_function

from sympy import sqrt, pprint
from sympy.physics.quantum import qapply
from sympy.physics.quantum.gate import H, X, Z, CNOT
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.circuitplot import circuit_plot
from sympy.physics.quantum.grover import superposition_basis


def main():
    psi = superposition_basis(2)
    psi

    # Dense coding demo:

    # Assume Alice has the left QBit in psi
    print("An even superposition of 2 qubits.  Assume Alice has the left QBit.")
    pprint(psi)

    # The corresponding gates applied to Alice's QBit are:
    # Identity Gate (1), Not Gate (X), Z Gate (Z), Z Gate and Not Gate (ZX)
    # Then there's the controlled not gate (with Alice's as control):CNOT(1, 0)
    # And the Hadamard gate applied to Alice's Qbit: H(1)

    # To Send Bob the message |0>|0>
    print("To Send Bob the message |00>.")
    circuit = H(1)*CNOT(1, 0)
    result = qapply(circuit*psi)
    result
    pprint(result)

    # To send Bob the message |0>|1>
    print("To Send Bob the message |01>.")
    circuit = H(1)*CNOT(1, 0)*X(1)
    result = qapply(circuit*psi)
    result
    pprint(result)

    # To send Bob the message |1>|0>
    print("To Send Bob the message |10>.")
    circuit = H(1)*CNOT(1, 0)*Z(1)
    result = qapply(circuit*psi)
    result
    pprint(result)

    # To send Bob the message |1>|1>
    print("To Send Bob the message |11>.")
    circuit = H(1)*CNOT(1, 0)*Z(1)*X(1)
    result = qapply(circuit*psi)
    result
    pprint(result)

if __name__ == "__main__":
    main()
