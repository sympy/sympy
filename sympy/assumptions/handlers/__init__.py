"""
Multipledispatch handlers for ``Predicate`` are implemented in this
module. Handlers in this module are not directly imported to other
module in order to avoid circular import problem.
"""

from .common import (AskHandlerClass, CommonHandler, AskCommutativeHandler,
    TautologicalHandler, test_closed_group)

__all__ = [
    'AskHandlerClass', 'CommonHandler', 'AskCommutativeHandler',
    'TautologicalHandler', 'test_closed_group'
]
