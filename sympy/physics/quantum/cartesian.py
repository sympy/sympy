"""Operators and states for 1D cartesian position and momentum."""

from sympy import I, S, sqrt, pi
from sympy import exp
from sympy import Interval, DiracDelta

from sympy.physics.quantum.operator import HermitianOperator
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.hilbert import L2


__all__ = [
    'XOp',
    'PxOp',
    'X',
    'Px',
    'XKet',
    'XBra',
    'PxKet',
    'PxBra'
]

class XOp(HermitianOperator):
    """1D cartesian position operator."""

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _eval_commutator_PxOp(self, other):
        return I*hbar

    def _apply_operator_XKet(self, ket):
        return ket.position*ket


class PxOp(HermitianOperator):
    """1D cartesian momentum operator."""

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _apply_operator_PxKet(self, ket):
        return ket.momentum*ket


X = XOp('X')
Px = PxOp('Px')


class XKet(Ket):
    """1D cartesian position eigenket."""

    @property
    def dual_class(self):
        return XBra

    @property
    def position(self):
        """The position of the state."""
        return self.label[0]

    def _eval_innerproduct_XBra(self, bra, **hints):
        return DiracDelta(self.position-bra.position)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return exp(-I*self.position*bra.momentum/hbar)/sqrt(2*pi*hbar)

    def _represent_default_basis(self, **options):
        return self._represent_XOp(None, **options)

    def _represent_XOp(self, basis, **options):
        return self.position


class XBra(Bra):
    """1D cartesian position eigenbra."""

    @property
    def dual_class(self):
        return XKet

    @property
    def position(self):
        """The position of the state."""
        return self.label[0]


class PxKet(Ket):
    """1D cartesian momentum eigenket."""

    @property
    def dual_class(self):
        return PxBra

    @property
    def momentum(self):
        """The momentum of the state."""
        return self.label[0]

    def _eval_innerproduct_XBra(self, bra, **hints):
        return exp(I*self.momentum*bra.position/hbar)/sqrt(2*pi*hbar)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return DiracDelta(self.momentum-bra.momentum)

    def _represent_default_basis(self, **options):
        return self._represent_PxOp(None, **options)

    def _represent_PxOp(self, basis, **options):
        return self.momentum


class PxBra(Bra):
    """1D cartesian momentum eigenbra."""

    @property
    def dual_class(self):
        return PxKet

    @property
    def momentum(self):
        """The momentum of the state."""
        return self.label[0]
