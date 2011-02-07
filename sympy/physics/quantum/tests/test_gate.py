from sympy.physics.quantum.gate import *
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.applyops import apply_operators
from sympy.physics.quantum.qubit import Qubit, IntQubit, qubit_to_matrix,\
     matrix_to_qubit

from sympy.matrices.matrices import Matrix

from sympy import exp
from sympy import symbols
from sympy import sqrt
from sympy.core.numbers import *

I = ImaginaryUnit()

def test_Gate():
    h = HadamardGate(1)
    assert h.min_qubits == 2
    assert h.nqubits == 1

def test_UGate():
    a,b,c,d = symbols('abcd')
    uMat = Matrix([[a,b],[c,d]])

    #test basic case where gate exists in 1-qubit space
    u1 = UGate((0,), uMat)
    assert represent(u1, ZGate(0), nqubits = 1) == uMat
    assert apply_operators(u1*Qubit('0')) == a*Qubit('0') + c*Qubit('1')
    assert apply_operators(u1*Qubit('1')) == b*Qubit('0') + d*Qubit('1')

    #test case where gate exists in a larger space
    u2 = UGate((1,), uMat)
    u2Rep = represent(u2, ZGate(0), nqubits=2)
    for i in range(4):
        assert u2Rep*qubit_to_matrix(IntQubit(i,2)) ==\
            qubit_to_matrix(apply_operators(u2*IntQubit(i,2)))

def test_CGate():
    #test single control functionality
    CNOTMatrix = Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    assert represent(CGate(1, XGate(0)), ZGate(0), nqubits=2) == CNOTMatrix

    #test multiple control bit functionality
    ToffoliGate = CGate((1,2), XGate(0))
    assert represent(ToffoliGate, ZGate(0), nqubits=3) == \
    Matrix([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],\
    [0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],\
    [0,0,0,0,0,0,1,0]])

    ToffoliGate = CGate((3,0), XGate(1))
    assert apply_operators(ToffoliGate*Qubit('1001')) == \
    matrix_to_qubit(represent(ToffoliGate*Qubit('1001'), ZGate(0), nqubits=4))
    assert apply_operators(ToffoliGate*Qubit('0000')) == \
    matrix_to_qubit(represent(ToffoliGate*Qubit('0000'), ZGate(0), nqubits=4))

    CZGate = CGate(0, ZGate(1))
    assert apply_operators(CZGate*Qubit('11')) == -Qubit('11')
    assert matrix_to_qubit(represent(CZGate*Qubit('11'),\
     ZGate(0), nqubits=2)) == -Qubit('11')

    CPhaseGate = CGate(0, PhaseGate(1))
    assert apply_operators(CPhaseGate*Qubit('11')) == ImaginaryUnit()*Qubit('11')
    assert matrix_to_qubit(represent(CPhaseGate*Qubit('11'),\
     ZGate(0), nqubits=2)) == ImaginaryUnit()*Qubit('11')

def test_UGate_CGate_combo():
    a,b,c,d = symbols('abcd')
    uMat = Matrix([[a,b],[c,d]])
    cMat = Matrix([[1,0,0,0],[0,1,0,0],[0,0,a,b],[0,0,c,d]])

    #test basic case where gate exists in 1-qubit space
    u1 = UGate((0,), uMat)
    cu1 = CGate(1, u1)
    assert represent(cu1, ZGate(0), nqubits = 2) == cMat
    assert apply_operators(cu1*Qubit('10')) == a*Qubit('10') + c*Qubit('11')
    assert apply_operators(cu1*Qubit('11')) == b*Qubit('10') + d*Qubit('11')
    assert apply_operators(cu1*Qubit('01')) == Qubit('01')
    assert apply_operators(cu1*Qubit('00')) == Qubit('00')

    #test case where gate exists in a larger space
    u2 = UGate((1,), uMat)
    cu2 = CGate(0, u2)
    u2Rep = represent(u2, ZGate(0), nqubits=2)
    for i in range(4):
        assert u2Rep*qubit_to_matrix(IntQubit(i,2)) ==\
            qubit_to_matrix(apply_operators(u2*IntQubit(i,2)))

def test_represent_Hadamard_():
    circuit = HadamardGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0), nqubits=2)
    # check that the answers are same to within an epsilon
    assert answer == Matrix([1/sqrt(2),1/sqrt(2), 0, 0])

def test_represent_XGate_():
    circuit = XGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0),nqubits=2)
    assert Matrix([0, 1, 0, 0]) == answer

def test_represent_YGate_():
    circuit = YGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0), nqubits=2)
    assert answer[0] == 0 and answer[1] == I and \
    answer[2] == 0 and answer[3] == 0

def test_represent_ZGate_():
    circuit = ZGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0), nqubits=2)
    assert Matrix([1, 0, 0, 0]) == answer

def test_represent_PhaseGate_():
    circuit = PhaseGate(0)*Qubit('01')
    answer = represent(circuit, ZGate(0), nqubits=2)
    assert Matrix([0, ImaginaryUnit(),0,0]) == answer

def test_represent_TGate_():
    circuit = TGate(0)*Qubit('01')
    assert Matrix([0, exp(I*Pi()/4), 0, 0]) == represent(circuit, ZGate(0), nqubits=2)

def test_CompoundGates_():
    circuit = YGate(0)*ZGate(0)*XGate(0)*HadamardGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0), nqubits=2)
    assert Matrix([ImaginaryUnit()/sqrt(2),ImaginaryUnit()/sqrt(2), 0, 0])\
    == answer

def test_CNOTGate():
    circuit = CNotGate(1,0)
    assert represent(circuit, ZGate(0), nqubits=2) == \
    Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    circuit = circuit*Qubit('111')
    assert matrix_to_qubit(represent(circuit, ZGate(0), nqubits=3)) == \
    apply_operators(circuit)

"""def test_gateSort():
    assert gate_sort(XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)*XGate(1)*XGate(0))\
     == HadamardGate(0)**2*XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)

def test_gate_simp():
     assert gate_simp(HadamardGate(0)*XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)\
     *XGate(1)**3*XGate(0)*ZGate(3)**2*PhaseGate(4)**3) == HadamardGate(0)*\
     XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)*ZGate(4)*PhaseGate(4)
"""

def test_SwapGate():
    SWAP_gate_matrix = Matrix(((1,0,0,0),(0,0,1,0),(0,1,0,0),(0,0,0,1)))
    #test SWAP gate decompose method
    assert represent(SwapGate(1,0).decompose(), ZGate(0), nqubits=2) == SWAP_gate_matrix
