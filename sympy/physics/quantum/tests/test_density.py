from sympy.physics.quantum.densityOp import Density, matrix_to_density
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.gate import HadamardGate
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.dagger import Dagger
from sympy.functions import sqrt

def test_matrix_to_density():
    assert matrix_to_density(represent(Density([Qubit('00'),1]), nqubits=2)) ==\
     Density([Qubit('00'),1])
    
    den = Density([Qubit('00'),1], [(Qubit('00')+Qubit('11'))/sqrt(2),1])
    assert represent(matrix_to_density(represent(den), nqubits=2))\
      == represent(den, nqubits=2)

def test_densityOp():
    assert Density([Qubit(0,0),1]) == Density([Qubit(0,0),1], 1, 1)
    
def test_representDensity():
    mat = represent(Qubit('01'), nqubits=2)
    assert represent(Density([Qubit('01'),1]), nqubits=2) == mat*Dagger(mat)    
    
def test_apply_operators_density():    
    assert apply_operators(Density([Qubit('00'),1]).operate_on(HadamardGate(0)\
           *HadamardGate(1)))\
        == Density([apply_operators(HadamardGate(0)*HadamardGate(1)*Qubit('00')),1])

