from sympy import Symbol, Matrix
from sympy.circuit.circuit import (Circuit, Resistor, Capacitor, Inductor,
    VoltageSource, CurrentSource, OpAmp, parse_netlist,
    _g_matrix, _b_matrix, _c_matrix, _d_matrix, _e_matrix, _j_matrix,
    _v_matrix, _i_matrix, _A_matrix, _z_matrix)


def test_solve():
    elements = []
    elements.append(VoltageSource('V1', '1', '0','12'))
    elements.append(Resistor('R1', '1', '2','1000'))
    elements.append(Resistor('R2', '2', '0','2000'))
    elements.append(Capacitor('C1', '2', '0','0.0002'))

    s = Symbol('s')
    V1 = Symbol('V1')
    R1 = Symbol('R1')
    R2 = Symbol('R2')
    C1 = Symbol('C1')

    circuit = Circuit(elements)
    solution = circuit.solve()
    assert solution.g == Matrix([[ 1/R1, -1/R1],[-1/R1, C1*s + 1/R2 + 1/R1]])
    assert solution.b == Matrix([[1], [0]])
    assert solution.c == Matrix([[1, 0]])
    assert solution.d == Matrix([[0]])
    assert solution.e == Matrix([[V1]])
    assert solution.i == Matrix([[0],[0]])
    assert solution.A == Matrix([[ 1/R1, -1/R1, 1],[-1/R1, C1*s + 1/R2 + 1/R1, 0],[1, 0, 0]])
    assert solution.Z == Matrix([[0], [0], [V1]])
    assert solution.x == Matrix([[V1],[R2*V1/(C1*R1*R2*s + R1 + R2)],[-V1*(-C1*R2*s - 1)/(-C1*R1*R2*s - R1 - R2)]])
    assert solution.I == {C1: C1*R2*V1*s/(C1*R1*R2*s + R1 + R2), R2: V1/(C1*R1*R2*s + R1 + R2), R1: (R2*V1/(C1*R1*R2*s + R1 + R2) - V1)/R1}
    assert solution.V == [0, V1, R2*V1/(C1*R1*R2*s + R1 + R2)]
