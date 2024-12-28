"""Kinds for Operators, Bras, and Kets."""

from sympy.core.mul import Mul
from sympy.core.kind import Kind, _NumberKind, NumberKind

class _KetKind(Kind):

    def __new__(cls):
        obj = super().__new__(cls)
        return obj

    def __repr__(self):
        return "KetKind"

KetKind = _KetKind()


class _BraKind(Kind):

    def __new__(cls):
        obj = super().__new__(cls)
        return obj

    def __repr__(self):
        return "BraKind"

BraKind = _BraKind()


from sympy.core.kind import Kind

class _OperatorKind(Kind):

    def __new__(cls):
        obj = super().__new__(cls)
        return obj

    def __repr__(self):
        return "OperatorKind"

OperatorKind = _OperatorKind()

#-----------------------------------------------------------------------------
# Kind resolution.
#-----------------------------------------------------------------------------

# Note: We can't currently add kind dispatchers for the following combinations
#       as the Mul._kind_dispatcher is set to commutative and will also
#       register the opposite order, which isn't correct for these pairs:
# 
# 1. (_OperatorKind, _KetKind)
# 2. (_BraKind, _OperatorKind)
# 3. (_BraKind, _KetKind)


@Mul._kind_dispatcher.register(_NumberKind, _KetKind)
def _mul_number_ket_kind(lhs, rhs):
    return KetKind


@Mul._kind_dispatcher.register(_NumberKind, _BraKind)
def _mul_number_bra_kind(lhs, rhs):
    return BraKind


@Mul._kind_dispatcher.register(_NumberKind, _OperatorKind)
def _mul_operator_kind(lhs, rhs):
    return OperatorKind
