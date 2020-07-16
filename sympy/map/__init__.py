"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'AppliedMap',
    'CompositeMap'
]

from .map import (
    Map, UndefinedMap, InverseMap, IdentityMap, AppliedMap
)
from .composite import (
    CompositeMap,
)
