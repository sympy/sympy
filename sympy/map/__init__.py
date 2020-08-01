"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'AppliedMap',
    'isapplied',
    'BinaryOperator', 'LeftDivision', 'RightDivision',
    'InverseOperator', 'ExponentOperator',
    'InverseElement', 'ExponentElement',
    'CompositeMap', 'IteratedMap',
    'AdditionOperator', 'Addition',
    'scalar_add',
    'MultiplicationOperator', 'Multiplication',
    'scalar_mul', 'scalar_pow', 'scalar_divide',
    'FunctionSet',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    AppliedMap,
    isapplied,
)
from .operator import (
    BinaryOperator, LeftDivision, RightDivision,
    InverseOperator, ExponentOperator,
    InverseElement, ExponentElement,
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
    scalar_mul, scalar_pow, scalar_divide,
)

from .functionspace import (
    FunctionSet,
)
