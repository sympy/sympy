"""
A module to implement logical predicates and assumption system.
"""
from __future__ import annotations

from .assume import (
    AppliedPredicate, Predicate, AssumptionsContext, assuming,
    global_assumptions
)
from .ask import Q, ask
from .refine import refine
from .relation import BinaryRelation, AppliedBinaryRelation

__all__ = [
    'AppliedPredicate', 'Predicate', 'AssumptionsContext', 'assuming',
    'global_assumptions', 'Q', 'ask',
    'refine',
    'BinaryRelation', 'AppliedBinaryRelation'
]
