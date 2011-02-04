from sympy import Matrix, symbols, I

from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.state import Bra, Ket
from sympy.physics.quantum.operator import Operator, OuterProduct
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.tensorproduct import matrix_tensor_product
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.innerproduct import InnerProduct

Amat = Matrix([[1,I],[-I,1]])
Bmat = Matrix([[1,2],[3,4]])
Avec = Matrix([[1],[I]])

class AKet(Ket):

    @property
    def dual_class(self):
        return ABra

    def _represent_AOp(self, basis, **options):
        return Avec


class ABra(Bra):

    @property
    def dual_class(self):
        return AKet


class AOp(Operator):

    def _represent_AOp(self, basis, **options):
        return Amat


class BOp(Operator):

    def _represent_AOp(self, basis, **options):
        return Bmat


k = AKet('a')
b = ABra('a')
A = AOp('A')
B = BOp('B')


def test_bra():
    assert represent(b, A) == Dagger(Avec)
    assert represent(Dagger(b), A) == Avec


def test_ket():
    assert represent(k, A) == Avec
    assert represent(Dagger(k), A) == Dagger(Avec)


def test_op():
    assert represent(A, A) == Amat
    assert represent(Dagger(A), A) == Dagger(Amat)


def test_outerproduct():
    op = k*b
    assert represent(op, A) == Avec*Avec.H


def test_tensor_product():
    assert represent(TensorProduct(A,B),A) == matrix_tensor_product(Amat,Bmat)


def test_pow():
    assert represent(A**2,A) == Amat**2


def test_add_mul():
    assert represent(A*B + 2*A, A) == Amat*Bmat + 2*Amat


def test_commutator():
    assert represent(Commutator(A,B), A) == Amat*Bmat - Bmat*Amat


def test_anticommutator():
    assert represent(AntiCommutator(A,B), A) == Amat*Bmat + Bmat*Amat


def test_innerproduct():
    assert represent(InnerProduct(b,k), A) == (Avec.H*Avec)[0]
