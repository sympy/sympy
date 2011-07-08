"""Operators and states for 1D cartesian position and momentum."""

from sympy import I, S, sqrt, pi
from sympy import exp
from sympy import Interval, DiracDelta
from sympy import Symbol

from sympy.physics.quantum.operator import HermitianOperator, DifferentialOperator
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
    def default_args(self):
        return ("X",)

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _eval_commutator_PxOp(self, other):
        return I*hbar

    def _apply_operator_XKet(self, ket):
        return ket.position*ket

    def _represent_PxKet(self, basis, **options):
        index = options.pop("index", 1)

        states = basis._enumerate_state(2, start_index = index)
        coord1 = states[0].momentum
        coord2 = states[1].momentum
        d = DifferentialOperator(coord1)
        delta = DiracDelta(coord1 - coord2)

        return I*hbar*(d*delta)


class PxOp(HermitianOperator):
    """1D cartesian momentum operator."""

    @classmethod
    def default_args(self):
        return ("Px",)

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _apply_operator_PxKet(self, ket):
        return ket.momentum*ket

    def _represent_XKet(self, basis, **options):
        index = options.pop("index", 1)

        states = basis._enumerate_state(2, start_index = index)
        coord1 = states[0].position
        coord2 = states[1].position
        d = DifferentialOperator(coord1)
        delta = DiracDelta(coord1 - coord2)

        return -I*hbar*(d*delta)

X = XOp('X')
Px = PxOp('Px')


class XKet(Ket):
    """1D cartesian position eigenket."""

    @classmethod
    def _operators_to_state(self, op, **options):
        return self.__new__(self, str(op.label[0]).lower(), **options)

    def _state_to_operators(self, op_class, **options):
        return op_class.__new__(op_class, str(self.label[0]).upper(), **options)

    @classmethod
    def default_args(self):
        return ("x",)

    @classmethod
    def dual_class(self):
        return XBra

    @property
    def position(self):
        """The position of the state."""
        return self.label[0]

    def _enumerate_state(self, num_states, **options):
        return _enumerate_continuous_1D(self, num_states, **options)

    def _eval_innerproduct_XBra(self, bra, **hints):
        return DiracDelta(self.position-bra.position)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return exp(-I*self.position*bra.momentum/hbar)/sqrt(2*pi*hbar)

class XBra(Bra):
    """1D cartesian position eigenbra."""

    @classmethod
    def default_args(self):
        return ("x",)

    @classmethod
    def dual_class(self):
        return XKet

    @property
    def position(self):
        """The position of the state."""
        return self.label[0]

class PxKet(Ket):
    """1D cartesian momentum eigenket."""

    @classmethod
    def _operators_to_state(self, op, **options):
        return self.__new__(self, str(op.label[0]).lower(), **options)

    def _state_to_operators(self, op_class, **options):
        lab = str(self.label[0])
        lab = lab[0].upper() + lab[1:]
        return op_class.__new__(op_class, lab, **options)

    @classmethod
    def default_args(self):
        return ("px",)

    @classmethod
    def dual_class(self):
        return PxBra

    @property
    def momentum(self):
        """The momentum of the state."""
        return self.label[0]

    def _enumerate_state(self, *args, **options):
        return _enumerate_continuous_1D(self, *args, **options)

    def _eval_innerproduct_XBra(self, bra, **hints):
        return exp(I*self.momentum*bra.position/hbar)/sqrt(2*pi*hbar)

    def _eval_innerproduct_PxBra(self, bra, **hints):
        return DiracDelta(self.momentum-bra.momentum)

class PxBra(Bra):
    """1D cartesian momentum eigenbra."""

    @classmethod
    def default_args(self):
        return ("px",)

    @classmethod
    def dual_class(self):
        return PxKet

    @property
    def momentum(self):
        """The momentum of the state."""
        return self.label[0]

def _enumerate_continuous_1D(*args, **options):
    state = args[0]
    num_states = args[1]
    state_class = state.__class__
    index_list = options.pop('index_list', [])

    if len(index_list) == 0:
        start_index = options.pop('start_index', 1)
        index_list = range(start_index, start_index + num_states)

    enum_states = [0 for i in range(len(index_list))]

    for i, ind in enumerate(index_list):
        label = state.args[0]
        enum_states[i] = state_class(str(label) + "_" + str(ind), **options)

    return enum_states
