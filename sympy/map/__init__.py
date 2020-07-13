"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'InverseMap', 'IdentityMap', 'AppliedMap'
]

from .map import (
    Map, InverseMap, IdentityMap, AppliedMap
)
