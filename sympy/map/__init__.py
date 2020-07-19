"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'AppliedMap',
    'CompositeMap', 'IteratedMap',
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, AppliedMap
)
from .binary import (
    BinaryOperator, AssociativeOperator
)
from .composite import (
    CompositeMap, IteratedMap,
)
