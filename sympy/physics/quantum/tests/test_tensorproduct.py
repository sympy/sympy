from sympy import I, symbols, Matrix

from sympy.physics.quantum.commutator import Commutator as Comm
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.tensorproduct import TensorProduct as TP
from sympy.physics.quantum.tensorproduct import tensor_product_simp
from sympy.physics.quantum.dagger import Dagger


A,B,C = symbols('A,B,C', commutative=False)
x = symbols('x')

mat1 = Matrix([[1,2*I],[1+I,3]])
mat2 = Matrix([[2*I,3],[4*I,2]])

def test_tensor_product_dagger():
    assert Dagger(TensorProduct(I*A, B)) ==\
           -I*TensorProduct(Dagger(A),Dagger(B))
    assert Dagger(TensorProduct(mat1,mat2)) ==\
        TensorProduct(Dagger(mat1),Dagger(mat2))

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
