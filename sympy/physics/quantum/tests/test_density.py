from sympy.physics.quantum.densityOp import Density
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.applyops import apply_operators

def test_densityOp():
    assert Density([1, Qubit(0,0)]) == Density([1, Qubit(0,0)], 1, 1)
    
def test_apply_operators_density():    
    assert apply_operators(Density([1, Qubit('00')]).operate_on(HadamardGate(0)*HadamardGate(1)))\
        == Density([1, apply_operators(HadamardGate(0)*HadamardGate(1)*Qubit('00'))])
    
