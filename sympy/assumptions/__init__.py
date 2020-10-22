"""
This module implements logical predicate and assumption system.

Predicate is a function which returns boolean value. Currently, only
unary predicate which takes a single argument is supported. Polyadic
predicate, which takes multiple arguments, is a work in progress.

Every predicate can be accessed via ``Q``, such as ``Q.even`` or
``Q.prime``. It is also possible to define custome predicates.

Applied predicate can be evaluated to boolean value by ``ask`` function.
The result is ``True`` if the proposition is true with respect to given
context, ``False`` if it is false, and ``None`` if it cannot be
determined.

"""

from .assume import (
    AppliedPredicate, Predicate, AssumptionsContext, assuming,
    global_assumptions
)
from .ask import Q, ask, register_handler
from .refine import refine

__all__ = [
    'AppliedPredicate', 'Predicate', 'AssumptionsContext', 'assuming',
    'global_assumptions', 'Q', 'ask', 'register_handler', 'refine',
]
