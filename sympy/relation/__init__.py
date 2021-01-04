"""
A module to implement finitary relations [1] as predicate.
This module expands assumption module to provide predicates for
finitary relations.

References
==========

 .. [1] https://en.wikipedia.org/wiki/Finitary_relation
"""

__all__ = ['BinaryRelation', 'AppliedBinaryRelation',
    'Equal',
    'eqnsimp']

from .binrel import BinaryRelation, AppliedBinaryRelation
from .equality import Equal
from .reltools import eqnsimp
