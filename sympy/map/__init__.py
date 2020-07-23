"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'RestrictedMap',
    'AppliedMap',
    'CompositeMap', 'IteratedMap',
    'BinaryOperator', 'AssociativeOperator', 'AppliedBinaryOperator', 'AppliedAssociativeOperator',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, RestrictedMap,
    AppliedMap,
)
from .operator import (
    BinaryOperator, AssociativeOperator, AppliedBinaryOperator, AppliedAssociativeOperator,
)
from .composite import (
    CompositeMap, IteratedMap,
)
