"""Kinds for Operators, Bras, and Kets."""

from sympy.core.kind import Kind

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