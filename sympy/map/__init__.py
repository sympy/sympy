"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'FunctionSet', 'function_set',
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'ConstantMap',
    'AppliedMap',
    'isappliedmap',
    'BinaryOperator', 'LeftDivisionOperator', 'RightDivisionOperator',
    'InverseOperator', 'ExponentOperator',
    'InverseElement', 'ExponentElement',
    'CompositionOperator', 'composite_op', 'CompositeMap',
    'IterationOperator', 'IteratedMap',
    'AdditionOperator', 'NumericAdditionOperator',
    'VectorAdditionOperator',
    'Addition',
    'MultiplicationOperator', 'NumericMultiplicationOperator',
    'ScalarMultiplicationOperator', 'VectorMultiplicationOperator',
    'Multiplication',

    'Sin', 'Cos', 'Tan',
]

from .map import (
    FunctionSet, function_set,
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    ConstantMap,
    AppliedMap,
    isappliedmap,
)
from .operator import (
    BinaryOperator, LeftDivisionOperator, RightDivisionOperator,
    InverseOperator, ExponentOperator,
    InverseElement, ExponentElement,
)
from .composite import (
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
    ScalarMultiplicationOperator, VectorMultiplicationOperator,
    Multiplication,
)

from .elementary import (
    Sin, Cos, Tan,
)
