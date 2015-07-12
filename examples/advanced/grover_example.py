#!/usr/bin/env python

"""Grover's quantum search algorithm example."""

from __future__ import division, print_function

from sympy import pprint
from sympy.physics.quantum import qapply
from sympy.physics.quantum.qubit import IntQubit
from sympy.physics.quantum.grover import (OracleGate, superposition_basis,
        WGate, grover_iteration)


def demo_vgate_app(v):
    for i in range(2**v.nqubits):
        print('qapply(v*IntQubit(%i, %r))' % (i, v.nqubits))
        pprint(qapply(v*IntQubit(i, v.nqubits)))
        qapply(v*IntQubit(i, v.nqubits))


def black_box(qubits):
    return True if qubits == IntQubit(1, qubits.nqubits) else False


def main():
    print()
    print('Demonstration of Grover\'s Algorithm')
    print('The OracleGate or V Gate carries the unknown function f(x)')
    print('> V|x> = ((-1)^f(x))|x> where f(x) = 1 when x = a (True in our case)')
    print('> and 0 (False in our case) otherwise')
    print()

    nqubits = 2
    print('nqubits = ', nqubits)

    v = OracleGate(nqubits, black_box)
    print('Oracle or v = OracleGate(%r, black_box)' % nqubits)
    print()

    psi = superposition_basis(nqubits)
    print('psi:')
    pprint(psi)
    demo_vgate_app(v)
    print('qapply(v*psi)')
    pprint(qapply(v*psi))
    print()

    w = WGate(nqubits)
    print('WGate or w = WGate(%r)' % nqubits)
    print('On a 2 Qubit system like psi, 1 iteration is enough to yield |1>')
    print('qapply(w*v*psi)')
    pprint(qapply(w*v*psi))
    print()

    nqubits = 3
    print('On a 3 Qubit system, it requires 2 iterations to achieve')
    print('|1> with high enough probability')
    psi = superposition_basis(nqubits)
    print('psi:')
    pprint(psi)

    v = OracleGate(nqubits, black_box)
    print('Oracle or v = OracleGate(%r, black_box)' % nqubits)
    print()

    print('iter1 = grover.grover_iteration(psi, v)')
    iter1 = qapply(grover_iteration(psi, v))
    pprint(iter1)
    print()

    print('iter2 = grover.grover_iteration(iter1, v)')
    iter2 = qapply(grover_iteration(iter1, v))
    pprint(iter2)
    print()

if __name__ == "__main__":
    main()
