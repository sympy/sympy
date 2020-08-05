"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'AppliedMap',
    'isappliedmap',
    'BinaryOperator', 'LeftDivision', 'RightDivision',
    'InverseOperator', 'ExponentOperator',
    'InverseElement', 'ExponentElement',
    'FunctionSet', 'function_set',
    'CompositionOperator', 'composite_op', 'CompositeMap',
    'IterationOperator', 'IteratedMap',
    'AdditionOperator', 'Addition',
    'scalar_add',
    'MultiplicationOperator', 'Multiplication',
    'scalar_mul', 'scalar_pow', 'scalar_divide',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    AppliedMap,
    isappliedmap,
)
from .operator import (
    BinaryOperator, LeftDivision, RightDivision,
    InverseOperator, ExponentOperator,
    InverseElement, ExponentElement,
)
from .composite import (
    FunctionSet, function_set,
    CompositionOperator, composite_op, CompositeMap,
    IterationOperator, IteratedMap,
)

from .add import (
    AdditionOperator, Addition,
    scalar_add,
)

from .mul import (
    MultiplicationOperator, Multiplication,
    scalar_mul, scalar_pow, scalar_divide,
)
