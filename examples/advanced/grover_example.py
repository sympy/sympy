from sympy import pprint
from sympy.physics.quantum import apply_operators
from sympy.physics.quantum.qubit import IntQubit
from sympy.physics.quantum.grover import *

def demo_vgate_app(v):
    for i in range(2**v.nqubits):
        # print 'apply_operators(v*IntQubit({0}, {1}))'.format(i, v.nqubits)
        # pprint(apply_operators(v*IntQubit(i, v.nqubits)))
        apply_operators(v*IntQubit(i, v.nqubits))

def black_box(qubits):
    return True if qubits == IntQubit(1, qubits.nqubits) else False

def main():
    # print ''
    # print 'Demonstration of Grover\'s Algorithm'
    # print 'The OracleGate or V Gate carries the unknown function f(x)'
    # print '> V|x> = ((-1)^f(x))|x> where f(x) = 1 (True in our case)'
    # print '> when x = a, 0 (False in our case) otherwise'
    # print ''

    nqubits = 2
    # print 'nqubits = ', nqubits

    v = OracleGate(nqubits, black_box)
    # print 'Oracle or v = OracleGate(%r, black_box)' % nqubits
    # print ''

    psi = superposition_basis(nqubits)
    # print 'psi:'
    # pprint(psi)
    # demo_vgate_app(v)
    # print 'apply_operators(v*psi)'
    # pprint(apply_operators(v*psi))
    # print ''

    w = WGate(nqubits)
    # print 'WGate or w = WGate(%r)' % nqubits
    # print 'On a 2 Qubit system like psi, 1 iteration is enough to yield |1>'
    # print 'apply_operators(w*v*psi)'
    # pprint(apply_operators(w*v*psi))
    # print ''

    nqubits = 3
    # print 'On a 3 Qubit system, it requires 2 iterations to achieve'
    # print '|1> with high enough probability'
    psi = superposition_basis(nqubits)
    # print 'psi:'
    # pprint(psi)
  
    v = OracleGate(nqubits, black_box)
    # print 'Oracle or v = OracleGate(%r, black_box)' % nqubits
    # print ''

    # print 'iter1 = grover.grover_iteration(psi, v)'
    iter1 = apply_operators(grover_iteration(psi, v))
    # pprint(iter1)
    # print '' 

    # print 'iter2 = grover.grover_iteration(iter1, v)'
    iter2 = apply_operators(grover_iteration(iter1, v))
    # pprint(iter2)
    # print '' 

"""
Allows application of Grover's algorithm to arbitrary number of qubits
Allows specification of number of iterations
Includes W gate and V (Oracle) gate
V gate accepts a "black box" function that returns true or false on a 
 set of qubits
"""
def conference_material():
    nqubits = 3
    psi = superposition_basis(nqubits); psi
    v = OracleGate(nqubits, black_box)
    # With a 3 qubit system, 2 iterations achieve high enough probability
    iter1 = apply_operators(grover_iteration(psi, v))
    iter2 = apply_operators(grover_iteration(iter1, v)); iter2

main()
conference_material()
