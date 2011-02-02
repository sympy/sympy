from sympy import I, Matrix, symbols, conjugate

from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.operator import (
    Operator, HermitianOperator, UnitaryOperator
)
from sympy.physics.quantum.tensorproduct import TensorProduct

def test_scalars():
    x,y,z = symbols('xyz')
    i,j,k = symbols('ijk',real=True)
    assert Dagger(x) == conjugate(x)
    assert Dagger(i) == i
    assert Dagger(I*x) == -I*conjugate(x)


def test_matrix():
    x = symbols('x')
    m = Matrix([[I,x*I],[2,4]])
    assert Dagger(m) == m.H


def test_bra_ket():
    x = symbols('x')
    k = Ket('k')
    b = Bra('b')
    assert Dagger(k) == Bra('k')
    assert Dagger(b) == Ket('b')

    k2 = Ket('k2')
    e = 2*I*k + x*k2
    assert Dagger(e) == conjugate(x)*Dagger(k2) - 2*I*Dagger(k)


def test_operator():
    A = Operator('A')
    B = Operator('B')
    assert isinstance(Dagger(A), Dagger)

    H = HermitianOperator('H')
    assert Dagger(H) == H

    U = UnitaryOperator('U')
    assert Dagger(U)*U == 1

    assert Dagger(A*B) == Dagger(B)*Dagger(A)
    assert Dagger(A+B) == Dagger(A) + Dagger(B)
    assert Dagger(A**2) == Dagger(A)**2


def test_outer_product():
    k = Ket('k')
    b = Bra('b')
    op = k*b
    assert Dagger(op) == Dagger(b)*Dagger(k)


def test_inner_product():
    k = Ket('k')
    b = Bra('b')
    ip = b*k
    assert Dagger(ip) == Dagger(k)*Dagger(b)


def test_commutator():
    A = Operator('A')
    B = Operator('B')
    C = Operator('C')
    comm = Commutator(A*B,C)
    assert Dagger(comm).expand(commutator=True) ==\
        - Commutator(Dagger(B),Dagger(C))*Dagger(A) -\
        Dagger(B)*Commutator(Dagger(A),Dagger(C))


def test_anticommutator():
    A = Operator('A')
    B = Operator('B')
    assert Dagger(AntiCommutator(A, B)) == AntiCommutator(Dagger(A),Dagger(B))


def test_tensor_product():
    A = Operator('A')
    B = Operator('B')
    assert Dagger(TensorProduct(I*A, B)) ==\
           -I*TensorProduct(Dagger(A),Dagger(B))
