"""
A module to implement the operators as instances of ``Expr``.

SymPy's ``sin``, ``cos``, etc are classes, not instances. Hence, combining operators
such as ``sin + cos`` is impossible - It is only possible by defining a ``Function``
class or ``Lambda`` instance.

This modlue provides a wrapper of such classes, so that the classes themselves
can behave as if they are instances of ``Expr``.
"""

from .operator import Op

__all__ = ['Op']
