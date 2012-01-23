from sympy.physics.quantum.qcevolve import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, CGate, gate_simp)

def test_find_subcircuit():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    x0 = X(0)
    x1 = X(1)

    circuit = (x, y, z)

    assert find_subcircuit(circuit, (x0,)) == 0
    assert find_subcircuit(circuit, (x1,)) == -1
    assert find_subcircuit(circuit, (y,)) == 1
    assert find_subcircuit(circuit, (h,)) == -1
    assert find_subcircuit(circuit, (x, h)) == -1
    assert find_subcircuit(circuit, (x, y, z)) == 0
    assert find_subcircuit(circuit, (y, z)) == 1
    assert find_subcircuit(circuit, (x, y, z, h)) == -1
    assert find_subcircuit(circuit, (z, y, x)) == -1

