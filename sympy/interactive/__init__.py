"""Helper module for setting up interactive SymPy sessions. """

__all__ = []

from .printing import init_printing
__all__ += ["init_printing"]

from .session import init_session
__all__ += ["init_session"]
