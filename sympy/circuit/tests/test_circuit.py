from sympy import Symbol, Matrix
from sympy.circuit.circuit import (Circuit, Resistor, Capacitor, Inductor,
    VoltageSource, CurrentSource, OpAmp, solve, parse_netlist,
    g_matrix, b_matrix, c_matrix, d_matrix, e_matrix, j_matrix,
    v_matrix, i_matrix, A_matrix, z_matrix)


def test_solve():
    cir = Circuit()
    cir.add_vsource('V1', '1', '0','12')
    cir.add_resistor('R1', '1', '2','1000')
    cir.add_resistor('R2', '2', '0','2000')
    cir.add_capacitor('C1', '2', '0','0.0002')

    s = Symbol('s')
    V1 = Symbol('V1')
    R1 = Symbol('R1')
    R2 = Symbol('R2')
    C1 = Symbol('C1')

    solve(cir)
    assert cir.G == Matrix([[ 1/R1, -1/R1],[-1/R1, C1*s + 1/R2 + 1/R1]])
    assert cir.B == Matrix([[1], [0]])
    assert cir.C == Matrix([[1, 0]])
    assert cir.D == Matrix([[0]])
    assert cir.E == Matrix([[V1]])
    assert cir.I == Matrix([[0],[0]])
    assert cir.A == Matrix([[ 1/R1, -1/R1, 1],[-1/R1, C1*s + 1/R2 + 1/R1, 0],[1, 0, 0]])
    assert cir.Z == Matrix([[0], [0], [V1]])
    assert cir.x == Matrix([[V1],[R2*V1/(C1*R1*R2*s + R1 + R2)],[-V1*(-C1*R2*s - 1)/(-C1*R1*R2*s - R1 - R2)]])
