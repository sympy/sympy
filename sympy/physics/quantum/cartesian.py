from sympy import I, S, sqrt, pi
from sympy import exp, conjugate
from sympy import Interval, DiracDelta

from sympy.physics.quantum.operator import HermitianOperator
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.hilbert import L2


class XOp(HermitianOperator):

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _eval_commutator_PxOp(self, other):
        return I*hbar

    def _apply_operator_XKet(self, ket):
        return ket.position*ket


class PxOp(HermitianOperator):

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _apply_operator_PxKet(self, ket):
        return ket.momentum*ket


X = XOp('X')
Px = PxOp('Px')


class XKet(Ket):

    @property
    def dual_class(self):
        return XBra

    @property
    def position(self):
        return self.label[0]

    def _eval_innerproduct_XBra(self, bra, **hints):
        return DiracDelta(self.position-bra.position)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return exp(-I*self.position*bra.momentum/hbar)/sqrt(2*pi*hbar)

    def _represent_XOp(self, basis, **options):
        return self.position


class XBra(Bra):

    @property
    def dual_class(self):
        return XKet

    @property
    def position(self):
        return self.label[0]

    def _represent_XOp(self, basis, **options):
        # TODO: really conjugate here?
        return conjugate(self.position)


class PxKet(Ket):

    @property
    def dual_class(self):
        return PxBra

    @property
    def momentum(self):
        return self.label[0]

    def _eval_innerproduct_XBra(self, bra, **hints):
        return exp(I*self.position*bra.momentum/hbar)/sqrt(2*pi*hbar)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return DiracDelta(self.momentum-bra.momentum)

    def _represent_PxOp(self, basis, **options):
        return self.momentum


class PxBra(Bra):

    @property
    def dual_class(self):
        return PxKet

    @property
    def momentum(self):
        return self.label[0]

    