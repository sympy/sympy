"""Definitions of common exceptions for :mod:`sympy.core` module. """

from __future__ import division, print_function


class BaseCoreError(Exception):
    """Base class for core related exceptions. """


class NonCommutativeExpression(BaseCoreError):
    """Raised when expression didn't have commutative property. """
