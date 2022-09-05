"""Helper module for setting up interactive SymPy sessions. """

from sympy.core.cache import lazy_function
from .traversal import interactive_traversal

init_printing = lazy_function('sympy.interactive.printing', 'init_printing')
init_session = lazy_function('sympy.interactive.session', 'init_session')


__all__ = ['init_printing', 'init_session', 'interactive_traversal']
