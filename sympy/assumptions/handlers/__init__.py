"""
Multipledispatch handlers for ``Predicate`` are implemented here.
Handlers in this module are not directly imported to other modules in
order to avoid circular import problem.
"""

from .common import test_closed_group

__all__ = [
    'test_closed_group'
]
