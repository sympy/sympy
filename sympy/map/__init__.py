"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'AppliedMap',
    'CompositeMap', 'IteratedMap',
    'BinaryOperator', 'LeftDivision', 'RightDivision',
    'InverseOperator', 'ExponentOperator',
    'AppliedBinaryOperator', 'InverseElement', 'ExponentElement',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    AppliedMap,
)
from .operator import (
    BinaryOperator, LeftDivision, RightDivision,
    InverseOperator, ExponentOperator,
    AppliedBinaryOperator, InverseElement, ExponentElement,
)
from .composite import (
    CompositeMap, IteratedMap,
)

from .add import (
    AdditionOperator, Addition,
    scalar_add,
)

from .mul import (
    MultiplicationOperator, Multiplication,
    scalar_mul, scalar_pow,
)
