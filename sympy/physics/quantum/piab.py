from sympy import Symbol, pi, sqrt, sin, conjugate, Interval, S

from sympy.physics.quantum.operator import HermitianOperator
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.kronecker import KroneckerDelta
from sympy.physics.quantum.hilbert import L2

m = Symbol('m')
L = Symbol('L')


class PIABHamiltonian(HermitianOperator):

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    def _apply_operator_PIABKet(self, ket, **options):
        n = ket.label[0]
        return (n**2*pi**2*hbar**2)/(2*m*L**2)*ket


class PIABKet(Ket):

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    @property
    def dual_class(self):
        return PIABBra

    def _represent_XOp(self, basis, **options):
        x = Symbol('x')
        n = Symbol('n')
        subs_info = options.get('subs',{})
        return sqrt(2/L)*sin(n*pi*x/L).subs(subs_info)

    def _eval_innerproduct_PIABBra(self, bra):
        return KroneckerDelta(bra.label[0], self.label[0])


class PIABBra(Bra):

    @classmethod
    def _eval_hilbert_space(cls, label):
        return L2(Interval(S.NegativeInfinity,S.Infinity))

    @property
    def dual_class(self):
        return PIABBra

    def _represent_XOp(self, basis, **options):
        return conjugate(self.dual._represent_XOp(basis, **options))
