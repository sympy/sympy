"""Implementations of pathways for use by actuators.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from __future__ import annotations

from abc import ABC, abstractmethod

from sympy.core.backend import S, Symbol
from sympy.core.expr import Expr
from sympy.physics.mechanics import Point
from sympy.physics.vector import Vector, dynamicsymbols


__all__ = ['LinearPathway']


class PathwayBase(ABC):
    """Abstract base class for all pathway classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom pathway types through subclassing.

    """

    def __init__(
        self,
        *,
        origin: Point,
        insertion: Point,
    ) -> None:
        """Initializer for ``PathwayBase``."""
        self.origin = origin
        self.insertion = insertion

    @property
    def origin(self) -> Point:
        """The ``Point`` from which the pathway originates."""
        return self._origin

    @origin.setter
    def origin(self, origin: Point) -> None:
        if not isinstance(origin, Point):
            msg = (
                f'Value {repr(origin)} passed to `origin` was of type '
                f'{type(origin)}, must be {type(Point)}.'
            )
            raise TypeError(msg)
        self._origin = origin

    @property
    def insertion(self) -> Point:
        """The ``Point`` into which the pathway inserts."""
        return self._insertion

    @insertion.setter
    def insertion(self, insertion: Point) -> None:
        if not isinstance(insertion, Point):
            msg = (
                f'Value {repr(insertion)} passed to `insertion` was of type '
                f'{type(insertion)}, must be {type(Point)}.'
            )
            raise TypeError(msg)
        self._insertion = insertion

    @abstractmethod
    def _true_length(self) -> Expr:
        """Exact analytical expression for the pathway's length."""
        pass

    @abstractmethod
    def _true_shortening_velocity(self) -> Expr:
        """Exact analytical expression for the pathway's shortening velocity."""
        pass

    @property
    def length(self) -> Expr:
        """An expression representing the pathway's length."""
        return self._true_length()

    @property
    def shortening_velocity(self) -> Expr:
        """An expression representing the pathway's shortening velocity."""
        return self._true_shortening_velocity()

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}(origin={self.origin}, '
            f'insertion={self.insertion})'
        )
