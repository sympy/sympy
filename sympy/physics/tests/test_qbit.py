from sympy.physics.qbit import *
from sympy import symbols, Rational
from sympy.core.numbers import *
from sympy.functions.elementary import *
import random
x, y = symbols('xy')

epsilon = .000001

def test_represent_Hadamard_Z():
    circuit = HadamardGate(0)*Qbit(0)
    answer = represent(circuit, ZBasisSet())
    # check that the answers are same to within an epsilon
    assert answer == Matrix([1/sqrt(2),1/sqrt(2)])

def test_represent_XGate_Z():
    circuit = XGate(0)*Qbit(0)
    answer = represent(circuit, ZBasisSet())
    assert Matrix([0, 1]) == answer

def test_represent_YGate_Z():
    circuit = YGate(0)*Qbit(0)
    answer = represent(circuit, ZBasisSet())
    assert answer[0] == 0 and answer[1] == ImaginaryUnit()

def test_represent_ZGate_Z():
    circuit = ZGate(0)*Qbit(0)
    answer = represent(circuit, ZBasisSet())
    assert Matrix([1, 0]) == answer

def test_represent_PhaseGate_Z():
    circuit = PhaseGate(0)*Qbit(1)
    answer = represent(circuit, ZBasisSet())
    assert Matrix([0, complex(0,1)]) == answer 

def test_represent_TGate_Z():
    pass

def test_CompoundGates_Z():
    circuit = YGate(0)*ZGate(0)*XGate(0)*HadamardGate(0)*Qbit(0)
    answer = represent(circuit)
    assert Matrix([.5*ImaginaryUnit()*sqrt(2),ImaginaryUnit()/sqrt(2)]) == answer

def test_tensor_product():
    try:
        import numpy as np
    except ImportError:
        return
    l1 = zeros(4)
    for i in range(16):
        l1[i] = 2**i
    l2 = zeros(4)
    for i in range(16):
        l2[i] = i
    l3 = zeros(2)
    for i in range(4):
        l3[i] = i
    vec = Matrix([1,2,3])

    #test for Matrix known 4x4 matricies
    numpyl1 = np.matrix(l1.tolist())
    numpyl2 = np.matrix(l2.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, l2]
    sympy_product = TensorProduct(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l2, l1]
    sympy_product = TensorProduct(*args)    
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for other known matrix of different dimensions
    numpyl2 = np.matrix(l3.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, l3]
    sympy_product = TensorProduct(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l3, l1]
    sympy_product = TensorProduct(*args)    
    assert numpy_product.tolist() == sympy_product.tolist()    

    #test for non square matrix
    numpyl2 = np.matrix(vec.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, vec]
    sympy_product = TensorProduct(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [vec, l1]
    sympy_product = TensorProduct(*args)    
    assert numpy_product.tolist() == sympy_product.tolist()   

    #test for random matrix with random values that are floats    
    random_matrix1 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    random_matrix2 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    numpy_product = np.kron(random_matrix1,random_matrix2)
    args = [Matrix(random_matrix1.tolist()),Matrix(random_matrix2.tolist())]
    sympy_product = TensorProduct(*args)
    assert not (sympy_product - Matrix(numpy_product.tolist())).tolist() > (ones((sympy_product.rows,sympy_product.cols))*epsilon).tolist()

    #test for three matrix kronecker
    sympy_product = TensorProduct(l1,vec,l2)
    npl1 = np.matrix(l1.tolist())
    npl2 = np.matrix(l2.tolist())
    npvec = np.matrix(vec.tolist())
    
    numpy_product = np.kron(l1,np.kron(vec,l2)) 
    assert numpy_product.tolist() == sympy_product.tolist()

#test apply methods
def test_apply_represent_equality():
    gates = [HadamardGate(int(5*random.random())), XGate(int(5*random.random())), ZGate(int(5*random.random())), YGate(int(5*random.random())), ZGate(int(5*random.random())), PhaseGate(int(5*random.random()))]
    
    circuit = Qbit(int(random.random()*2),int(random.random()*2),int(random.random()*2),int(random.random()*2),int(random.random()*2),int(random.random()*2))
    circuit = HadamardGate(2)*HadamardGate(4)*HadamardGate(5)*HadamardGate(1)*HadamardGate(0)*HadamardGate(3)*circuit
    for i in range(int(random.random()*6)):
        circuit = gates[int(random.random()*6)]*circuit


    mat = represent(circuit)
    states = apply_gates(circuit)
    state_rep = matrix_to_qbits(mat)
    states = states.expand()
    state_rep = state_rep.expand()
    print "\nthis is circuit and results",circuit
    print states  
    print state_rep 
    assert state_rep == states

