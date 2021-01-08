"""
A module to implement finitary relations [1] as predicate.

This module expands assumption module to provide predicates for
finitary relations. Relations can be symbolically manipulated, evaluated
to boolean, and assumed.

References
==========

 .. [1] https://en.wikipedia.org/wiki/Finitary_relation
"""

__all__ = ['BinaryRelation', 'AppliedBinaryRelation',
    'Equal',
    'GreaterThan', 'GreaterEq', 'LessThan', 'LessEq',
    'rearrange', 'eqnsimp', 'solveeqn']

from .binrel import BinaryRelation, AppliedBinaryRelation
from .equality import Equal
from .inequality import GreaterThan, GreaterEq, LessThan, LessEq
from .reltools import rearrange, eqnsimp, solveeqn
