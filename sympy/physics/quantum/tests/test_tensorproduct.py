from sympy import Matrix, I, symbols
from sympy.matrices.matrices import *

from sympy.physics.quantum.commutator import Commutator as Comm
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.tensorproduct import TensorProduct as TP
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.tensorproduct import tensor_product_simp

from sympy.physics.quantum.dagger import Dagger

epsilon = .000001

A,B,C = symbols('ABC', commutative=False)
x,y,z = symbols('xyz')

def test_matrix_tensor_product():
    try:
        import numpy as np
    except ImportError:
        print 'import error numpy'
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
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l2, l1]
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for other known matrix of different dimensions
    numpyl2 = np.matrix(l3.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, l3]
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l3, l1]
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for non square matrix
    numpyl2 = np.matrix(vec.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, vec]
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [vec, l1]
    sympy_product = matrix_tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for random matrix with random values that are floats
    random_matrix1 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    random_matrix2 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    numpy_product = np.kron(random_matrix1,random_matrix2)
    args = [Matrix(random_matrix1.tolist()),Matrix(random_matrix2.tolist())]
    sympy_product = matrix_tensor_product(*args)
    assert not (sympy_product - Matrix(numpy_product.tolist())).tolist() > \
    (ones((sympy_product.rows,sympy_product.cols))*epsilon).tolist()

    #test for three matrix kronecker
    sympy_product = matrix_tensor_product(l1,vec,l2)

    numpy_product = np.kron(l1,np.kron(vec,l2))
    assert numpy_product.tolist() == sympy_product.tolist()


def test_tensor_product_dagger():
    assert Dagger(TensorProduct(I*A, B)) ==\
           -I*TensorProduct(Dagger(A),Dagger(B))


def test_tensor_product_abstract():

    assert TP(x*A,2*B) == x*2*TP(A,B)
    assert TP(A,B) != TP(B,A)
    assert TP(A,B).is_commutative == False
    assert isinstance(TP(A,B), TP)
    assert TP(A,B).subs(A,C) == TP(C,B)


def test_tensor_product_expand():
    assert TP(A+B,B+C).expand(tensorproduct=True) ==\
        TP(A,B) + TP(A,C) + TP(B,B) + TP(B,C)


def test_tensor_product_commutator():
    assert TP(Comm(A,B),C).doit().expand(tensorproduct=True) ==\
        TP(A*B,C) - TP(B*A,C)
    assert Comm(TP(A,B),TP(B,C)).doit() ==\
        TP(A,B)*TP(B,C) - TP(B,C)*TP(A,B)


def test_tensor_product_simp():
    assert tensor_product_simp(TP(A,B)*TP(B,C)) == TP(A*B,B*C)
