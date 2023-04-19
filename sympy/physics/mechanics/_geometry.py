"""Implementations of geometry objects for use by wrapping pathways.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.core.expr import Expr
    from sympy.physics.mechanics import Point, Vector


__all__ = []


class GeometryBase(ABC):
    """Abstract base class for all geometry classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom geometry types through subclassing.

    """

    @abstractmethod
    def _point_is_on_surface(self, point: Point) -> bool:
        """Determine if a point is on the geometry's surface.

        Parameters
        ==========
        point : Point
            The point for which it's to be ascertained if it's on the
            geometry's surface or not.

        """

    @abstractmethod
    def geodesic_length(self, point_1: Point, point_2: Point) -> Expr:
        """The shortest distance between two points on a geometry's surface.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic length should be calculated.
        point_2 : Point
            The point to which the geodesic length should be calculated.

        """
        pass

    def __repr__(self) -> str:
        """Default representation of a geometry object."""
        return f'{self.__class__.__name__}()'
