"""A kind for Operators."""

from sympy.core.kind import Kind

class _OperatorKind(Kind):

    def __new__(cls):
        obj = super().__new__(cls)
        return obj

    def __repr__(self):
        return "OperatorKind"

OperatorKind = _OperatorKind()
