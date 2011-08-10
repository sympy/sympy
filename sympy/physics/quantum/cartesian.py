"""Operators and states for 1D cartesian position and momentum.

TODO:
1) Add 3D classes to mappings in operatorset.py

"""

from sympy import (DiracDelta, Derivative, exp, Function, I, Interval,
                   pi, S, sqrt, Symbol)
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.hilbert import L2
from sympy.physics.quantum.operator import DifferentialOperator, HermitianOperator
from sympy.physics.quantum.state import Ket, Bra, State

__all__ = [
    'XOp',
    'YOp',
    'ZOp',
    'PxOp',
    'X',
    'Y',
    'Z',
    'Px',
    'XKet',
    'XBra',
    'PxKet',
    'PxBra',
    'PositionState3D',
    'PositionKet3D',
    'PositionBra3D'
]

#-------------------------------------------------------------------------
# Position operators
#-------------------------------------------------------------------------

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

    def _apply_operator_PositionKet3D(self, ket):
        return ket.position_x*ket

    def _represent_PxKet(self, basis, **options):
        index = options.pop("index", 1)

        states = basis._enumerate_state(2, start_index = index)
        coord1 = states[0].momentum
        coord2 = states[1].momentum
        f = Function('f')
        d = DifferentialOperator(Derivative(f(coord1), coord1), f(coord1))
        delta = DiracDelta(coord1 - coord2)

        return I*hbar*delta*d

    def _represent_XKet(self, basis, **options):
        from sympy.physics.quantum.represent import expectation_helper
        options['basis'] = basis
        return expectation_helper(self, wrap_wavefunction=True, **options)

    def _get_default_basis(self, **options):
        return XKet

class YOp(HermitianOperator):
    """ Y cartesian coordinate operator (for 2D or 3D systems) """

    @classmethod
    def default_args(self):
        return ("Y",)

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _apply_operator_PositionKet3D(self, ket):
        return ket.position_y*ket

class ZOp(HermitianOperator):
    """ Z cartesian coordinate operator (for 3D systems) """

    @classmethod
    def default_args(self):
        return ("Z",)

    @classmethod
    def _eval_hilbert_space(self, args):
        return L2(Interval(S.NegativeInfinity, S.Infinity))

    def _apply_operator_PositionKet3D(self, ket):
        return ket.position_z*ket

#-------------------------------------------------------------------------
# Momentum operators
#-------------------------------------------------------------------------

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
        f = Function('f')
        d = DifferentialOperator(Derivative(f(coord1), coord1), f(coord1))
        delta = DiracDelta(coord1 - coord2)

        return -I*hbar*delta*d

    def _represent_PxKet(self, basis, **options):
        from sympy.physics.quantum.represent import expectation_helper
        options['basis'] = basis
        return expectation_helper(self, wrap_wavefunction=True, **options)

    def _get_default_basis(self, **options):
        return PxKet

X = XOp('X')
Y = YOp('Y')
Z = ZOp('Z')
Px = PxOp('Px')

#-------------------------------------------------------------------------
# Position eigenstates
#-------------------------------------------------------------------------

class XKet(Ket):
    """1D cartesian position eigenket."""

    @classmethod
    def _operators_to_state(self, op, **options):
        return self.__new__(self, *_lowercase_labels(op), **options)

    def _state_to_operators(self, op_class, **options):
        return op_class.__new__(op_class, \
                                *_uppercase_labels(self), **options)

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

    def _represent_XKet(self, basis, **options):
        from sympy.physics.quantum.represent import innerproduct_helper
        options['basis'] = basis
        return innerproduct_helper(self, wrap_wavefunction=True, **options)

    def _represent_PxKet(self, basis, **options):
        from sympy.physics.quantum.represent import innerproduct_helper
        options['basis'] = basis
        return innerproduct_helper(self, wrap_wavefunction=True, **options)

    def _get_default_basis(self, **options):
        return XKet

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

class PositionState3D(State):
    """ Base class for 3D cartesian position eigenstates """

    @classmethod
    def _operators_to_state(self, ops, **options):
        if not isinstance(ops, set) or len(ops) != 3:
            raise TypeError("This is not the set of commuting operators expected!")

        op_list = list(ops)
        order = [XOp, YOp, ZOp]
        op_list.sort(lambda a, b: _class_cmp(order, a, b))
        return self.__new__(self, *_lowercase_labels(op_list), **options)

    def _state_to_operators(self, op_classes, **options):
        if not isinstance(op_classes, set) or len(op_classes) != 3:
            raise TypeError("This is not the set of commuting operators expected!")

        op_list = list(op_classes)
        order = [XOp, YOp, ZOp]
        op_list.sort(lambda a, b: _class_cmp(order, a, b, is_class=True))

        labs = _uppercase_labels(self)
        op_insts = [None for op in op_list]

        for i, cls in enumerate(op_list):
            op_insts[i] = cls.__new__(cls, labs[i], **options)

        return set(op_insts)

    @classmethod
    def default_args(self):
        return ("x", "y", "z")

    @property
    def position_x(self):
        """ The x coordinate of the state """
        return self.label[0]

    @property
    def position_y(self):
        """ The y coordinate of the state """
        return self.label[1]

    @property
    def position_z(self):
        """ The z coordinate of the state """
        return self.label[2]

class PositionKet3D(Ket, PositionState3D):
    """ 3D cartesian position eigenket """

    def _eval_innerproduct_PositionBra3D(self, bra, **options):
        x_diff = self.position_x - bra.position_x
        y_diff = self.position_y - bra.position_y
        z_diff = self.position_z - bra.position_z

        return DiracDelta(x_diff)*DiracDelta(y_diff)*DiracDelta(z_diff)

    @classmethod
    def dual_class(self):
        return PositionBra3D

class PositionBra3D(Bra, PositionState3D):
    """ 3D cartesian position eigenbra """

    @classmethod
    def dual_class(self):
        return PositionKet3D

#-------------------------------------------------------------------------
# Momentum eigenstates
#-------------------------------------------------------------------------

class PxKet(Ket):
    """1D cartesian momentum eigenket."""

    @classmethod
    def _operators_to_state(self, op, **options):
        return self.__new__(self, *_lowercase_labels(op), **options)

    def _state_to_operators(self, op_class, **options):
        return op_class.__new__(op_class, \
                                *_uppercase_labels(self), **options)

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

    def _represent_PxKet(self, basis, **options):
        from sympy.physics.quantum.represent import innerproduct_helper
        options['basis'] = basis
        return innerproduct_helper(self, wrap_wavefunction=True, **options)

    def _represent_XKet(self, basis, **options):
        from sympy.physics.quantum.represent import innerproduct_helper
        options['basis'] = basis
        return innerproduct_helper(self, wrap_wavefunction=True, **options)

    def _get_default_basis(self, **options):
        return PxKet


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

#-------------------------------------------------------------------------
# Global helper functions
#-------------------------------------------------------------------------

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

def _lowercase_labels(ops):
    if not isinstance(ops, list):
        ops = [ops]

    return [str(arg.label[0]).lower() for arg in ops]

def _uppercase_labels(state):
    new_args = [str(lab)[0].upper() + \
                str(lab)[1:] for lab in state.label]

    return new_args

def _class_cmp(ordering, a, b, **options):
    ordering = list(ordering)

    is_class = options.pop("is_class", False)

    if not is_class:
        a = a.__class__
        b = b.__class__

    try:
        ind1 = ordering.index(a)
        ind2 = ordering.index(b)
    except ValueError:
        raise TypeError("Expected Operators for this eigenstate not found!")

    if ind1 < ind2:
        return -1
    else:
        return 1
