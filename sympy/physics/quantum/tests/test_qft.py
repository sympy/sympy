from sympy.physics.quantum.qft import *
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.represent import represent
from sympy.matrices.matrices import Matrix, eye
from sympy.core.symbol import Symbol
from sympy import exp, I, pi


def test_RkGate():
    x = Symbol('x')
    assert RkGate(1,x).k == x
    assert RkGate(1,x).targets == (1,)
    assert RkGate(1,1) == ZGate(1)
    assert RkGate(2,2) == PhaseGate(2)
    assert RkGate(3,3) == TGate(3)

    assert represent(RkGate(0,x), ZGate(0), nqubits =1) ==\
    Matrix([[1,0],[0,exp(2*I*pi/2**x)]])

def test_RkGate_controled():
    pass

def test_quantum_fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)*CGate((0,), PhaseGate(1))\
    *HadamardGate(1)*CGate((0,), TGate(2))*CGate((1,), PhaseGate(2))*HadamardGate(2)
    
    assert IQFT(0,3).decompose() == HadamardGate(2)*CGate((1,), RkGate(2,-2))*CGate((0,),RkGate(2,-3))\
    *HadamardGate(1)*CGate((0,), RkGate(1,-2))*HadamardGate(0)*SwapGate(0,2)
    
    assert represent(QFT(0,3)*IQFT(0,3), ZGate(0), nqubits=3)\
     == eye(8)
     
    assert QFT(0,4).decompose() #non-trivial decomposition
    assert apply_operators(QFT(0,3).decompose()*Qubit((0,0,0))).expand() ==\
     apply_operators(HadamardGate(0)*HadamardGate(1)*HadamardGate(2)*Qubit((0,0,0))).expand()
     
     
def test_Fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)\
    *RkGate(1,0,2)*HadamardGate(1)*RkGate(2,0,3)*RkGate(2,1,2)*HadamardGate(2)
    assert QFT(0,3).input_number == 2
    assert IQFT(0,3).decompose() == HadamardGate(2)*IRkGate(2,1,2)\
    *IRkGate(2,0,3)*HadamardGate(1)*IRkGate(1,0,2)*HadamardGate(0)*SwapGate(0,2)
