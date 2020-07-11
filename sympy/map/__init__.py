"""
A module that implements maps as Expr.

Includes any maps such as function, differential operator, etc.
"""

__all__ = [
    'Map', 'AppliedMap'
]

from .map import Map, AppliedMap
