"""
This module contains various types of mirrors.

Following sign convention is followed:

    For a spherically curved mirror in air, the magnitude of the focal length
    is equal to the radius of curvature of the mirror divided by two. The focal
    length is positive for a concave mirror, and negative for a convex mirror.
    In the sign convention used in optical design, a concave mirror has negative
    radius of curvature.

**Contains**

* Mirror
* SphericalMirror
* PlaneMirror
* ConvexMirror
* ConcaveMirror

"""

from __future__ import division

__all__ = ['Mirror', 'SphericalMirror', 'ConcaveMirror', 'ConvexMirror']

from sympy import sympify, sqrt, symbols, oo
from sympy.core.basic import Basic
from sympy.core.compatibility import is_sequence
from sympy.geometry.point import Point
from sympy.geometry.point3d import Point3D
from sympy.geometry.line3d import Ray3D, Line3D
from sympy.geometry.plane import Plane
from sympy.matrices import Matrix


class Mirror(Basic):

    """
    This is the base class of all kinds of mirrors. This class is not for
    public usage.

    """
    pass


class SphericalMirror(Mirror):

    """
    This is the base class of spherical mirrors.
    Plane mirror is also derived from this class with a limiting condition
    on its radius as it becomes infinite.

    Notes
    =====

    This assumes that the aperture is very small and all the rays are
    paraxial rays.

    Arguments
    =========

    center : Point3D or a sequence of length 3
        This is the center of curvature.
    radius : Sympifiable
        Radius of the curvature.

    """
    def __new__(cls, center, radius):
        obj = super(SphericalMirror, cls).__new__(cls)
        if not isinstance(center, Point3D):
            if not is_sequence(center):
                raise TypeError("'center' can only be Point3D or any sequence")
            else:
                # This will automatically take care of the
                # length of the sequence
                obj._center = Point3D(center)
        else:
            obj._center = center
            obj._radius = sympify(radius)
            return obj

    @property
    def center(self):
        """
        Returns the center of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.center
        Point3D(0, 0, 0)
        """
        return self._center

    @property
    def pole(self):
        """
        Returns the pole of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.pole
        Point3D(1, 2, 3)
        """
        return self._pole

    @property
    def focus(self):
        """
        Returns focus of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.focus
        Point3D(1/2, 1, 1)
        """
        return self._pole.midpoint()


class ConcaveMirror(SphericalMirror):

    """
    A concave mirror.

    Arguments
    =========

    center : Point3D or a sequence of length 3
        This is the center of curvature.
    pole : Point3D or a sequence of length 3
        This is the pole of the mirror.
    aperture_size : Sympifiable
        Linear aperture size of the mirror.

    Examples
    ========

    >>> from sympy import Point3D
    >>> from sympy.physics.optics import ConcaveMirror
    >>> m = ConcaveMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))

    """
    def __new__(cls, center, pole, aperture_size=symbols('a')):
        if not isinstance(center, Point3D):
            if not is_sequence(center):
                raise TypeError("'center' can only be Point3D or any sequence")
            else:
                center = Point3D(center)
        if not isinstance(pole, Point3D):
            if not is_sequence(pole):
                raise TypeError("'pole' can only be Point3D or any sequence")
            else:
                pole = Point3D(pole)
        radius = pole.distance(center)
        tbase = sqrt(radius**2 - (aperture_size/2)**2)
        int_x = (pole.x*tbase + center.x*(radius - tbase))/radius
        int_y = (pole.y*tbase + center.y*(radius - tbase))/radius
        int_z = (pole.z*tbase + center.z*(radius - tbase))/radius
        # These are the coordinates of the center of the circle which
        # cuts the the sphere.
        pt = Point3D(int_x, int_y, int_z)
        aprt_plane = Plane(pt, normal_vector=Ray3D(pt, pole).direction_ratio)

        obj = super(ConcaveMirror, cls).__new__(cls, center, radius)

        obj._radius = radius
        obj._center = center
        obj._pole = pole
        obj._aperture_size = aperture_size
        return obj

    @property
    def radius(self):
        """
        Returns the radius of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConcaveMirror
        >>> m = ConcaveMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.radius
        -sqrt(14)
        """
        return -self._radius

    @property
    def focal_length(self):
        """
        Returns focal length of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConcaveMirror
        >>> m = ConcaveMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.focal_length
        sqrt(14)/2
        """
        return -self.radius/2

    @property
    def power(self):
        """
        Returns power of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConcaveMirror
        >>> m = ConcaveMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.power
        sqrt(14)/7
        """
        return 1/self.focal_length

    def principal_axis(self, eq=False):
        """
        Returns the principal axis of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConcaveMirror
        >>> m = ConcaveMirror(Point3D(1, 1, 1), Point3D(-1, 0, 0))
        >>> m.principal_axis()
        Line3D(Point3D(1, 1, 1), Point3D(-1, 0, 0))
        >>> m.principal_axis(eq=True)
        (-x/2 + 1/2, -y + 1, -z + 1, k)

        """
        if eq:
            return Line3D(self._center, self._pole).equation()
        else:
            return Line3D(self._center, self._pole)

    def __repr__(self):
        from sympy.printing import sstr
        return type(self).__name__ + sstr(
            (self._center, self._pole, self._aperture_size)
        )

    __str__ = __repr__


class ConvexMirror(SphericalMirror):

    """
    A convex mirror.

    Arguments
    =========

    center : Point3D or a sequence of length 3
        This is the center of curvature.
    pole : Point3D or a sequence of length 3
        This is the pole of the mirror.
    aperture_size : Sympifiable
        Linear aperture size of the mirror.

    Examples
    ========

    >>> from sympy import Point3D
    >>> from sympy.physics.optics import ConvexMirror
    >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))

    """
    def __new__(cls, center, pole, aperture_size=symbols('a')):
        if not isinstance(center, Point3D):
            if not is_sequence(center):
                raise TypeError("'center' can only be Point3D or any sequence")
            else:
                center = Point3D(center)
        if not isinstance(pole, Point3D):
            if not is_sequence(pole):
                raise TypeError("'pole' can only be Point3D or any sequence")
            else:
                pole = Point3D(pole)
        radius = pole.distance(center)
        tbase = sqrt(radius**2 - (aperture_size/2)**2)
        int_x = (pole.x*tbase + center.x*(radius - tbase))/radius
        int_y = (pole.y*tbase + center.y*(radius - tbase))/radius
        int_z = (pole.z*tbase + center.z*(radius - tbase))/radius
        # These are the coordinates of the center of the circle which
        # cuts the the sphere.
        # WIP: Calculating extreme points of the curvature
        obj = super(ConvexMirror, cls).__new__(cls, center, radius)

        obj._radius = radius
        obj._center = center
        obj._pole = pole
        obj._aperture_size = aperture_size
        return obj

    @property
    def radius(self):
        """
        Returns the radius of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.radius
        sqrt(14)
        """
        return self._radius

    @property
    def focal_length(self):
        """
        Returns focal length of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.focal_length
        -sqrt(14)/2
        """
        return -self.radius/2

    @property
    def power(self):
        """
        Returns power of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(0, 0, 0), Point3D(1, 2, 3))
        >>> m.power
        -sqrt(14)/7
        """
        return 1/self.focal_length

    def principal_axis(self, eq=False):
        """
        Returns the principal axis of the mirror.

        Examples
        ========

        >>> from sympy import Point3D
        >>> from sympy.physics.optics import ConvexMirror
        >>> m = ConvexMirror(Point3D(1, 1, 1), Point3D(2, 0, 0))
        >>> m.principal_axis()
        Line3D(Point3D(1, 1, 1), Point3D(2, 0, 0))
        >>> m.principal_axis(eq=True)
        (x - 1, -y + 1, -z + 1, k)

        """
        if eq:
            return Line3D(self._center, self._pole).equation()
        else:
            return Line3D(self._center, self._pole)

    def __repr__(self):
        from sympy.printing import sstr
        return type(self).__name__ + sstr(
            (self._center, self._pole, self._aperture_size)
        )

    __str__ = __repr__
