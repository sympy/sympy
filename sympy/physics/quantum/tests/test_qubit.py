from sympy.physics.quantum.qubit import *
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.qft import *
from sympy.physics.quantum.represent import *
from sympy.physics.quantum.applyops import *
from sympy import symbols, Rational
from sympy.core.numbers import *
from sympy.functions.elementary import *
from sympy.physics.quantum.shor import *
from sympy.core.containers import Tuple
import random
x, y = symbols('xy')

epsilon = .000001

def test_Qubit():
    array = [0,0,1,1,0]
    qb = Qubit('00110')
    assert qb.flip(0) == Qubit('00111')
    assert qb.flip(1) == Qubit('00100')
    assert qb.flip(4) == Qubit('10110')
    assert qb.dimension == 5
    for i in range(5):
        assert qb[i] == array[4-i]
    assert len(qb) == 5
    qb = Qubit('110')

def test_Fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)\
    *RkGate(1,0,2)*HadamardGate(1)*RkGate(2,0,3)*RkGate(2,1,2)*HadamardGate(2)
    assert QFT(0,3).input_number == 2
    assert IQFT(0,3).decompose() == HadamardGate(2)*IRkGate(2,1,2)\
    *IRkGate(2,0,3)*HadamardGate(1)*IRkGate(1,0,2)*HadamardGate(0)*SwapGate(0,2)

def test_represent_HilbertSpace():
    import numpy as np
    a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = symbols('abcdefghijklmnop')
    gateMat = Matrix([[a,b,c,d],[e,f,g,h],[i,j,k,l],[m,n,o,p]])
    assert represent_hilbert_space(gateMat, 3, (0,1)) == \
    Matrix([[a,c,b,d,0,0,0,0],[i,k,j,l,0,0,0,0],[e,g,f,h,0,0,0,0],\
    [m,o,n,p,0,0,0,0],[0,0,0,0,a,c,b,d],[0,0,0,0,i,k,j,l],\
    [0,0,0,0,e,g,f,h],[0,0,0,0,m,o,n,p]])
    assert type(represent_hilbert_space(gateMat, 2, \
    (0,1), format = 'numpy')) == type(np.matrix(1))

def test_gateSort():
    assert gate_sort(XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)*XGate(1)*XGate(0))\
     == HadamardGate(0)**2*XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)

def test_gate_simp():
     assert gate_simp(HadamardGate(0)*XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)\
     *XGate(1)**3*XGate(0)*ZGate(3)**2*PhaseGate(4)**3) == HadamardGate(0)*\
     XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)*ZGate(4)*PhaseGate(4)

def test_ArbMat4_Equality():

    class Arb(Gate):
        @property
        def matrix(self):
            a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = symbols('abcdefghijklmnop')
            return Matrix([[a,b,c,d],[e,f,g,h],[i,j,k,l],[m,n,o,p]])

    for i in range(4):
        for j in range(4):
            if j != i:
                assert apply_operators(Arb(i,j)*(Qubit('10110'))) ==\
                matrix_to_qubits(represent(Arb(i,j)*Qubit('10110'),\
                ZGate(0)**5))

def test_Arb8_Matrix_Equality():
    class Arb(Gate):
        @property
        def matrix(self):
            a,b,c,d,e,f,g,h = symbols('abcdefgh')
            symlist = [a,b,c,d,e,f,g,h]
            lout = []
            for i in range(8):
                lin = []
                for j in range(8):
                    lin.append(symlist[i]**j)
                lout.append(lin)
            return Matrix(lout)

    for i in range(1):
        for j in range(4):
            for k in range(4):
                if j != i and k != i and k != j:
                    assert apply_operators(Arb(i,j,k)*(Qubit((0,1,1,1,0)))) ==\
                     matrix_to_qubits(represent(Arb(i,j,k)*Qubit((0,1,1,1,0)),\
                     ZGate(0)**5))

def test_superposition_of_states():
    assert apply_operators(CNOTGate(0,1)*HadamardGate(0)*(1/sqrt(2)*Qubit((0,1))\
     + 1/sqrt(2)*Qubit((1,0)))).expand() == (Qubit((0,1))/2 + Qubit((0,0))/2 - Qubit((1,1))/2 +\
     Qubit((1,0))/2)

    assert matrix_to_qubits(represent(CNOTGate(0,1)*HadamardGate(0)\
    *(1/sqrt(2)*Qubit((0,1)) + 1/sqrt(2)*Qubit((1,0))), ZGate(0)**(2)))\
     == (Qubit((0,1))/2 + Qubit((0,0))/2 - Qubit((1,1))/2 + Qubit((1,0))/2)


#test apply methods
def test_apply_represent_equality():
    gates = [HadamardGate(int(3*random.random())),\
     XGate(int(3*random.random())), ZGate(int(3*random.random())),\
      YGate(int(3*random.random())), ZGate(int(3*random.random())),\
       PhaseGate(int(3*random.random()))]

    circuit = Qubit((int(random.random()*2),int(random.random()*2),\
    int(random.random()*2),int(random.random()*2),int(random.random()*2),\
    int(random.random()*2)))
    for i in range(int(random.random()*6)):
        circuit = gates[int(random.random()*6)]*circuit


    mat = represent(circuit, ZGate(0)**(6))
    states = apply_operators(circuit)
    state_rep = matrix_to_qubits(mat)
    states = states.expand()
    state_rep = state_rep.expand()
    assert state_rep == states

def test_matrix_to_qubits():
    assert matrix_to_qubits(Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))\
    == Qubit((0,0,0,0))
    assert qubits_to_matrix(Qubit((0,0,0,0))) ==\
    Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    assert matrix_to_qubits(sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])) ==\
    (2*sqrt(2)*(Qubit((0,0,0)) + Qubit((0,0,1)) + Qubit((0,1,0)) + Qubit((0,1,1))\
    + Qubit((1,0,0)) + Qubit((1,0,1)) + Qubit((1,1,0)) + Qubit((1,1,1)))).expand()
    assert qubits_to_matrix(2*sqrt(2)*(Qubit((0,0,0)) + Qubit((0,0,1)) + Qubit((0,1,0))\
    + Qubit((0,1,1)) + Qubit((1,0,0)) + Qubit((1,0,1)) + Qubit((1,1,0)) + Qubit((1,1,1))))\
    == sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])

def test_RkGate_and_inverse():
    assert RkGate(1,2,x).k == x
    assert RkGate(1,2,x).args[0] == Tuple(1,2)
    assert IRkGate(1,2,x).k == x
    assert IRkGate(1,2,x).args[0] == Tuple(1,2)

    assert represent(RkGate(0,1,2), ZGate(0)**2) ==\
    Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,\
    exp(2*ImaginaryUnit()*Pi()/2**2)]])

    assert represent(IRkGate(0,1,3), ZGate(0)**(2)) ==\
    Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,\
    exp(-2*ImaginaryUnit()*Pi()/2**3)]])

def test_quantum_fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)*RkGate(1,0,2)\
    *HadamardGate(1)*RkGate(2,0,3)*RkGate(2,1,2)*HadamardGate(2)
    assert IQFT(0,3).decompose() == HadamardGate(2)*IRkGate(2,1,2)*IRkGate(2,0,3)\
    * HadamardGate(1)*IRkGate(1,0,2)*HadamardGate(0)*SwapGate(0,2)
    assert represent(QFT(0,3).decompose()*IQFT(0,3).decompose(), ZGate(0)**(3))\
     == eye(8)
    assert apply_operators(QFT(0,3).decompose()*Qubit((0,0,0))) ==\
     apply_operators(HadamardGate(0)*HadamardGate(1)*HadamardGate(2)*Qubit((0,0,0))).expand()
