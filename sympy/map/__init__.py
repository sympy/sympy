"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'AppliedMap',
    'isappliedmap',
    'BinaryOperator', 'LeftDivisionOperator', 'RightDivisionOperator',
    'InverseOperator', 'ExponentOperator',
    'InverseElement', 'ExponentElement',
    'FunctionSet', 'function_set',
    'CompositionOperator', 'composite_op', 'CompositeMap',
    'IterationOperator', 'IteratedMap',
    'AdditionOperator', 'NumericAdditionOperator',
    'VectorAdditionOperator',
    'Addition',
    'MultiplicationOperator', 'NumericMultiplicationOperator',
    'ScalarMultiplicationOperator',
    'Multiplication',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    AppliedMap,
    isappliedmap,
)
from .operator import (
    BinaryOperator, LeftDivisionOperator, RightDivisionOperator,
    InverseOperator, ExponentOperator,
    InverseElement, ExponentElement,
)
from .composite import (
    FunctionSet, function_set,
    CompositionOperator, composite_op, CompositeMap,
    IterationOperator, IteratedMap,
)

from .add import (
    AdditionOperator, NumericAdditionOperator,
    VectorAdditionOperator,
    Addition,
)

from .mul import (
    MultiplicationOperator, NumericMultiplicationOperator,
    ScalarMultiplicationOperator,
    Multiplication,
)
