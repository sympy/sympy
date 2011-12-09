from sympy.physics.quantum.densityOp import Density, matrix_to_density
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.gate import HadamardGate
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.dagger import Dagger
from sympy.functions import sqrt


    
def test_representDensity():
    mat = represent(Qubit('01'), nqubits=2)
    assert represent(Density([Qubit('01'),1]), nqubits=2) == mat*Dagger(mat)    
    
def test_apply_operators_density():    
    assert apply_operators(Density([Qubit('00'),1]).operate_on(HadamardGate(0)\
           *HadamardGate(1)))\
        == Density([apply_operators(HadamardGate(0)*HadamardGate(1)*Qubit('00')),1])

