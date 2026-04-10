"""
A module to implement finitary relations [1] as predicate.

References
==========

.. [1] https://en.wikipedia.org/wiki/Finitary_relation

"""
from __future__ import annotations

__all__ = ['BinaryRelation', 'AppliedBinaryRelation']

from .binrel import BinaryRelation, AppliedBinaryRelation
