"""Geometry objects for use by wrapping pathways."""

from abc import ABC, abstractmethod

from sympy import Integer, acos, pi, sqrt, sympify, tan, cos, sin, Piecewise
from sympy.core.relational import Eq
from sympy.functions.elementary.trigonometric import atan2
from sympy.polys.polytools import cancel
from sympy.physics.vector import Vector, dot
from sympy.simplify.simplify import trigsimp
from sympy.physics.mechanics import Point


__all__ = [
    'WrappingGeometryBase',
    'WrappingCylinder',
    'WrappingSphere',
    'WrappingCone'
]


class WrappingGeometryBase(ABC):
    """Abstract base class for all geometry classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom geometry types through subclassing.

    """

    @property
    @abstractmethod
    def point(cls):
        """The point with which the geometry is associated."""
        pass

    @abstractmethod
    def geodesic_length(self, point_1, point_2):
        """Returns the shortest distance between two points on a geometry's
        surface.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic length should be calculated.
        point_2 : Point
            The point to which the geodesic length should be calculated.

        """
        pass

    @abstractmethod
    def geodesic_end_vectors(self, point_1, point_2):
        """The vectors parallel to the geodesic at the two end points.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic originates.
        point_2 : Point
            The point at which the geodesic terminates.

        """
        pass

    def __repr__(self):
        """Default representation of a geometry object."""
        return f'{self.__class__.__name__}()'


class WrappingSphere(WrappingGeometryBase):
    """A solid spherical object.

    Explanation
    ===========

    A wrapping geometry that allows for circular arcs to be defined between
    pairs of points. These paths are always geodetic (the shortest possible).

    Examples
    ========

    To create a ``WrappingSphere`` instance, a ``Symbol`` denoting its radius
    and ``Point`` at which its center will be located are needed:

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import Point, WrappingSphere
    >>> r = symbols('r')
    >>> pO = Point('pO')

    A sphere with radius ``r`` centered on ``pO`` can be instantiated with:

    >>> WrappingSphere(r, pO)
    WrappingSphere(radius=r, point=pO)

    Parameters
    ==========

    radius : Symbol
        Radius of the sphere. This symbol must represent a value that is
        positive and constant, i.e. it cannot be a dynamic symbol, nor can it
        be an expression.
    point : Point
        A point at which the sphere is centered.

    See Also
    ========

    WrappingCylinder: Cylindrical geometry where the wrapping direction can be
        defined.

    """

    def __init__(self, radius, point):
        """Initializer for ``WrappingSphere``.

        Parameters
        ==========

        radius : Symbol
            The radius of the sphere.
        point : Point
            A point on which the sphere is centered.

        """
        self.radius = radius
        self.point = point

    @property
    def radius(self):
        """Radius of the sphere."""
        return self._radius

    @radius.setter
    def radius(self, radius):
        self._radius = radius

    @property
    def point(self):
        """A point on which the sphere is centered."""
        return self._point

    @point.setter
    def point(self, point):
        self._point = point

    def point_on_surface(self, point):
        """Returns a symbolic equality for whether a point is on the sphere's
        surface.

        Parameters
        ==========

        point : Point
            The point for which the expression is to be generated.

        """
        point_vector = point.pos_from(self.point)
        if isinstance(point_vector, Vector):
            point_radius_squared = dot(point_vector, point_vector)
        else:
            point_radius_squared = point_vector**2
        return Eq(point_radius_squared, self.radius**2, evaluate=False)

    def geodesic_length(self, point_1, point_2):
        r"""Returns the shortest distance between two points on the sphere's
        surface.

        Explanation
        ===========

        The geodesic length, i.e. the shortest arc along the surface of a
        sphere, connecting two points can be calculated using the formula:

        .. math::

           l = \arccos\left(\mathbf{v}_1 \cdot \mathbf{v}_2\right)

        where $\mathbf{v}_1$ and $\mathbf{v}_2$ are the unit vectors from the
        sphere's center to the first and second points on the sphere's surface
        respectively. Note that the actual path that the geodesic will take is
        undefined when the two points are directly opposite one another.

        Examples
        ========

        A geodesic length can only be calculated between two points on the
        sphere's surface. Firstly, a ``WrappingSphere`` instance must be
        created along with two points that will lie on its surface:

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import (Point, ReferenceFrame,
        ...     WrappingSphere)
        >>> N = ReferenceFrame('N')
        >>> r = symbols('r')
        >>> pO = Point('pO')
        >>> pO.set_vel(N, 0)
        >>> sphere = WrappingSphere(r, pO)
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')

        Let's assume that ``p1`` lies at a distance of ``r`` in the ``N.x``
        direction from ``pO`` and that ``p2`` is located on the sphere's
        surface in the ``N.y + N.z`` direction from ``pO``. These positions can
        be set with:

        >>> p1.set_pos(pO, r*N.x)
        >>> p1.pos_from(pO)
        r*N.x
        >>> p2.set_pos(pO, r*(N.y + N.z).normalize())
        >>> p2.pos_from(pO)
        sqrt(2)*r/2*N.y + sqrt(2)*r/2*N.z

        The geodesic length, which is in this case is a quarter of the sphere's
        circumference, can be calculated using the ``geodesic_length`` method:

        >>> sphere.geodesic_length(p1, p2)
        pi*r/2

        If the ``geodesic_length`` method is passed an argument that doesn't
        lie on the sphere's surface then unexpected results may be obtained.
        It is the user's responsibility to ensure that the points lie on the
        sphere's surface.

        Parameters
        ==========

        point_1 : Point
            Point from which the geodesic length should be calculated.
        point_2 : Point
            Point to which the geodesic length should be calculated.

        """
        point_1_vector = point_1.pos_from(self.point).normalize()
        point_2_vector = point_2.pos_from(self.point).normalize()
        central_angle = acos(point_2_vector.dot(point_1_vector))
        geodesic_length = self.radius*central_angle
        return geodesic_length

    def geodesic_end_vectors(self, point_1, point_2):
        """The vectors parallel to the geodesic at the two end points.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic originates.
        point_2 : Point
            The point at which the geodesic terminates.

        """
        pA, pB = point_1, point_2
        pO = self.point
        pA_vec = pA.pos_from(pO)
        pB_vec = pB.pos_from(pO)

        if pA_vec.cross(pB_vec) == 0:
            msg = (
                f'Can\'t compute geodesic end vectors for the pair of points '
                f'{pA} and {pB} on a sphere {self} as they are diametrically '
                f'opposed, thus the geodesic is not defined.'
            )
            raise ValueError(msg)

        return (
            pA_vec.cross(pB.pos_from(pA)).cross(pA_vec).normalize(),
            pB_vec.cross(pA.pos_from(pB)).cross(pB_vec).normalize(),
        )

    def __repr__(self):
        """Representation of a ``WrappingSphere``."""
        return (
            f'{self.__class__.__name__}(radius={self.radius}, '
            f'point={self.point})'
        )


class WrappingCylinder(WrappingGeometryBase):
    """A solid (infinite) cylindrical object.

    Explanation
    ===========

    A wrapping geometry that allows for circular arcs to be defined between
    pairs of points. These paths are always geodetic (the shortest possible) in
    the sense that they will be a straight line on the unwrapped cylinder's
    surface. However, it is also possible for a direction to be specified, i.e.
    paths can be influenced such that they either wrap along the shortest side
    or the longest side of the cylinder. To define these directions, rotations
    are in the positive direction following the right-hand rule.

    Examples
    ========

    To create a ``WrappingCylinder`` instance, a ``Symbol`` denoting its
    radius, a ``Vector`` defining its axis, and a ``Point`` through which its
    axis passes are needed:

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import (Point, ReferenceFrame,
    ...     WrappingCylinder)
    >>> N = ReferenceFrame('N')
    >>> r = symbols('r')
    >>> pO = Point('pO')
    >>> ax = N.x

    A cylinder with radius ``r``, and axis parallel to ``N.x`` passing through
    ``pO`` can be instantiated with:

    >>> WrappingCylinder(r, pO, ax)
    WrappingCylinder(radius=r, point=pO, axis=N.x)

    Parameters
    ==========

    radius : Symbol
        The radius of the cylinder.
    point : Point
        A point through which the cylinder's axis passes.
    axis : Vector
        The axis along which the cylinder is aligned.

    See Also
    ========

    WrappingSphere: Spherical geometry where the wrapping direction is always
        geodetic.

    """

    def __init__(self, radius, point, axis):
        """Initializer for ``WrappingCylinder``.

        Parameters
        ==========

        radius : Symbol
            The radius of the cylinder. This symbol must represent a value that
            is positive and constant, i.e. it cannot be a dynamic symbol.
        point : Point
            A point through which the cylinder's axis passes.
        axis : Vector
            The axis along which the cylinder is aligned.

        """
        self.radius = radius
        self.point = point
        self.axis = axis

    @property
    def radius(self):
        """Radius of the cylinder."""
        return self._radius

    @radius.setter
    def radius(self, radius):
        self._radius = radius

    @property
    def point(self):
        """A point through which the cylinder's axis passes."""
        return self._point

    @point.setter
    def point(self, point):
        self._point = point

    @property
    def axis(self):
        """Axis along which the cylinder is aligned."""
        return self._axis

    @axis.setter
    def axis(self, axis):
        self._axis = axis.normalize()

    def point_on_surface(self, point):
        """Returns a symbolic equality for whether a point is on the cylinder's
        surface.

        Parameters
        ==========

        point : Point
            The point for which the expression is to be generated.

        """
        relative_position = point.pos_from(self.point)
        parallel = relative_position.dot(self.axis) * self.axis
        point_vector = relative_position - parallel
        if isinstance(point_vector, Vector):
            point_radius_squared = dot(point_vector, point_vector)
        else:
            point_radius_squared = point_vector**2
        return Eq(trigsimp(point_radius_squared), self.radius**2, evaluate=False)

    def geodesic_length(self, point_1, point_2):
        """The shortest distance between two points on a geometry's surface.

        Explanation
        ===========

        The geodesic length, i.e. the shortest arc along the surface of a
        cylinder, connecting two points. It can be calculated using Pythagoras'
        theorem. The first short side is the distance between the two points on
        the cylinder's surface parallel to the cylinder's axis. The second
        short side is the arc of a circle between the two points of the
        cylinder's surface perpendicular to the cylinder's axis. The resulting
        hypotenuse is the geodesic length.

        Examples
        ========

        A geodesic length can only be calculated between two points on the
        cylinder's surface. Firstly, a ``WrappingCylinder`` instance must be
        created along with two points that will lie on its surface:

        >>> from sympy import symbols, cos, sin
        >>> from sympy.physics.mechanics import (Point, ReferenceFrame,
        ...     WrappingCylinder, dynamicsymbols)
        >>> N = ReferenceFrame('N')
        >>> r = symbols('r')
        >>> pO = Point('pO')
        >>> pO.set_vel(N, 0)
        >>> cylinder = WrappingCylinder(r, pO, N.x)
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')

        Let's assume that ``p1`` is located at ``N.x + r*N.y`` relative to
        ``pO`` and that ``p2`` is located at ``r*(cos(q)*N.y + sin(q)*N.z)``
        relative to ``pO``, where ``q(t)`` is a generalized coordinate
        specifying the angle rotated around the ``N.x`` axis according to the
        right-hand rule where ``N.y`` is zero. These positions can be set with:

        >>> q = dynamicsymbols('q')
        >>> p1.set_pos(pO, N.x + r*N.y)
        >>> p1.pos_from(pO)
        N.x + r*N.y
        >>> p2.set_pos(pO, r*(cos(q)*N.y + sin(q)*N.z).normalize())
        >>> p2.pos_from(pO).simplify()
        r*cos(q(t))*N.y + r*sin(q(t))*N.z

        The geodesic length, which is in this case a is the hypotenuse of a
        right triangle where the other two side lengths are ``1`` (parallel to
        the cylinder's axis) and ``r*q(t)`` (parallel to the cylinder's cross
        section), can be calculated using the ``geodesic_length`` method:

        >>> cylinder.geodesic_length(p1, p2).simplify()
        sqrt(r**2*q(t)**2 + 1)

        If the ``geodesic_length`` method is passed an argument ``Point`` that
        doesn't lie on the cylinder's surface then unexpected results may occur.
        It is the user's responsibility to ensure that the points lie on the
        cylinder's surface.

        Parameters
        ==========

        point_1 : Point
            Point from which the geodesic length should be calculated.
        point_2 : Point
            Point to which the geodesic length should be calculated.

        """
        relative_position = point_2.pos_from(point_1)
        parallel_length = relative_position.dot(self.axis)

        point_1_relative_position = point_1.pos_from(self.point)
        point_1_perpendicular_vector = (
            point_1_relative_position
            - point_1_relative_position.dot(self.axis)*self.axis
        ).normalize()

        point_2_relative_position = point_2.pos_from(self.point)
        point_2_perpendicular_vector = (
            point_2_relative_position
            - point_2_relative_position.dot(self.axis)*self.axis
        ).normalize()

        central_angle = _directional_atan(
            cancel(point_1_perpendicular_vector
                .cross(point_2_perpendicular_vector)
                .dot(self.axis)),
            cancel(point_1_perpendicular_vector.dot(point_2_perpendicular_vector)),
        )

        planar_arc_length = self.radius*central_angle
        geodesic_length = sqrt(parallel_length**2 + planar_arc_length**2)
        return geodesic_length

    def geodesic_end_vectors(self, point_1, point_2):
        """The vectors parallel to the geodesic at the two end points.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic originates.
        point_2 : Point
            The point at which the geodesic terminates.

        """
        point_1_from_origin_point = point_1.pos_from(self.point)
        point_2_from_origin_point = point_2.pos_from(self.point)

        if point_1_from_origin_point == point_2_from_origin_point:
            msg = (
                f'Cannot compute geodesic end vectors for coincident points '
                f'{point_1} and {point_2} as no geodesic exists.'
            )
            raise ValueError(msg)

        point_1_parallel = point_1_from_origin_point.dot(self.axis) * self.axis
        point_2_parallel = point_2_from_origin_point.dot(self.axis) * self.axis
        point_1_normal = (point_1_from_origin_point - point_1_parallel)
        point_2_normal = (point_2_from_origin_point - point_2_parallel)

        if point_1_normal == point_2_normal:
            point_1_perpendicular = Vector(0)
            point_2_perpendicular = Vector(0)
        else:
            point_1_perpendicular = self.axis.cross(point_1_normal).normalize()
            point_2_perpendicular = -self.axis.cross(point_2_normal).normalize()

        geodesic_length = self.geodesic_length(point_1, point_2)
        relative_position = point_2.pos_from(point_1)
        parallel_length = relative_position.dot(self.axis)
        planar_arc_length = sqrt(geodesic_length**2 - parallel_length**2)

        point_1_vector = (
            planar_arc_length * point_1_perpendicular
            + parallel_length * self.axis
        ).normalize()
        point_2_vector = (
            planar_arc_length * point_2_perpendicular
            - parallel_length * self.axis
        ).normalize()

        return (point_1_vector, point_2_vector)

    def __repr__(self):
        """Representation of a ``WrappingCylinder``."""
        return (
            f'{self.__class__.__name__}(radius={self.radius}, '
            f'point={self.point}, axis={self.axis})'
        )


class WrappingCone(WrappingGeometryBase):
    """A solid (infinite) conical object.

    Explanation
    ===========

    A wrapping geometry that allows for circular arcs to be defined between
    pairs of points on the surface of a cone. These paths are always geodetic
    (the shortest possible) in the sense that they become straight lines on
    the unwrapped conical surface.

    Examples
    ========

    To create a ``WrappingCone`` instance, a ``Symbol`` denoting its semi-vertical
    angle, a ``Point`` defining its apex, and a ``Vector`` specifying its axis
    are needed:

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import Point, ReferenceFrame, WrappingCone
    >>> N = ReferenceFrame('N')
    >>> alpha = symbols('alpha')
    >>> pO = Point('pO')
    >>> ax = N.z

    A cone with semi-vertical angle ``alpha``, apex at ``pO``, and axis aligned
    with ``N.z`` can be instantiated with:

    >>> WrappingCone(alpha, pO, ax)
    WrappingCone(alpha=alpha, apex=pO, axis=N.z)

    Parameters
    ==========

    alpha : Symbol
        The semi-vertical angle of the cone.
    apex : Point
        The tip of the cone where the curved surface meets.
    axis : Vector
        The axis along which the cone is aligned.

    See Also
    ========

    WrappingCylinder
        Cylindrical geometry where wrapping arcs are geodetic on the cylinder.
    WrappingSphere
        Spherical geometry where the wrapping direction is always geodetic.

    """
    def __init__(self, alpha, apex, axis):
        """
        Initializer for ``WrappingCone``.

        Parameters
        ==========

        alpha: Symbol
            The semi vertical angle of the cone.

        apex: Point
            The tip of the cone where the curved surface meets.

        axis: Vector
            The axis along which the cone is aligned.

        """
        if alpha.is_number:
            if alpha == 0:
                raise ValueError(
                    "Cone angle alpha must be positive."
                )
            if alpha == pi / 2:
                raise ValueError(
                    "Cone angle alpha must be less than pi/2."
                )
        elif alpha.is_real is False or alpha.is_positive is False:
            raise ValueError(
                    "Cone angle alpha must be real and positive."
                )
        self._alpha = alpha

        if not isinstance(apex, Point):
            raise TypeError("The 'apex' must be a Point object.")
        self._apex = apex

        if not isinstance(axis, Vector):
            raise TypeError("The 'axis' must be a Vector object.")
        self._axis = axis.normalize()

    @property
    def point(self):
        """This method is implemented as required by WrappingGeometryBase, use WrappingCone.apex instead."""
        return self._apex

    @property
    def alpha(self):
        """The semi vertical angle of the cone."""
        return self._alpha

    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha

    @property
    def apex(self):
        """The tip of the cone where the curved surface meets."""
        return self._apex

    @apex.setter
    def apex(self, apex):
        self._apex = apex

    @property
    def axis(self):
        """The axis along which the cone is aligned."""
        return self._axis

    @axis.setter
    def axis(self, axis):
        self._axis = axis.normalize()

    def point_on_surface(self, point):
        """
        Returns a symbolic equality for whether a point is on the cone's
        surface.

        Parameters
        ==========

        point : Point
            The point for which the expression is to be generated.
        """
        position = point.pos_from(self.apex)
        axis_component = position.dot(self.axis) * self.axis
        radial_component = position - axis_component
        lhs = radial_component.dot(radial_component)
        rhs = axis_component.dot(axis_component) * tan(self.alpha) ** 2
        return Eq(lhs, rhs, evaluate=False)

    def geodesic_length(self, point_1, point_2):
        """
        The shortest distance between two points on a conical surface.

        Explanation
        ===========
        Computes the geodesic by "unwrapping" the cone into a planar sector
        and measuring the straight line distance between the corresponding
        points in that sector.

        Examples
        ========
        >>> from sympy import pi, sqrt
        >>> from sympy.physics.mechanics import Point, ReferenceFrame, WrappingCone
        >>> N = ReferenceFrame('N')
        >>> alpha = pi/6
        >>> apex = Point('O')
        >>> cone = WrappingCone(alpha, apex, N.z)
        >>> p1 = Point('A')
        >>> p1.set_pos(apex, N.x / sqrt(3) + N.z)
        >>> p2 = Point('B')
        >>> p2.set_pos(apex, N.y / sqrt(3) + N.z)
        >>> cone.geodesic_length(p1, p2)
        sqrt(8/3 - 4*sqrt(2)/3)

        Parameters
        ==========
        point_1 : Point
            Starting point on the cone's surface.
        point_2 : Point
            Ending point on the cone's surface.
        """
        pos1 = point_1.pos_from(self.apex)
        pos2 = point_2.pos_from(self.apex)
        z1 = pos1.dot(self.axis)
        z2 = pos2.dot(self.axis)
        s1 = z1 / cos(self.alpha)
        s2 = z2 / cos(self.alpha)
        n1 = (pos1 - z1 * self.axis).normalize()
        n2 = (pos2 - z2 * self.axis).normalize()
        central = _directional_atan(
            cancel((n1.cross(n2)).dot(self.axis)),
            cancel(n1.dot(n2))
        )

        central = Piecewise((central, central <= pi), (2 * pi - central, True))

        delta_u = central * sin(self.alpha)
        return sqrt(s1**2 + s2**2 - 2*s1*s2*cos(delta_u))


    def geodesic_end_vectors(self, point_1, point_2):
        """
        Computes the unit tangent vectors for the geodesic path at its
        endpoints.

        The tangent vector of the geodesic is found by considering the
        straight-line path in the unrolled planar sector representation of the
        cone. This vector is then mapped back to the 3D space at each
        endpoint.
        """
        # Get the position vectors of the points relative to the cone's apex.
        # All subsequent calculations are performed in the cone's reference frame.
        pos1 = point_1.pos_from(self.apex)
        pos2 = point_2.pos_from(self.apex)

        # A unique geodesic cannot be defined between two identical points.
        if pos1 == pos2:
            raise ValueError(
                f'No unique geodesic exists for coincident points {point_1} and {point_2}.'
            )

        # If one point is the apex, the geodesic is the straight line along the
        # cone's surface to the other point. The tangent at the apex points
        # towards the other point, and the tangent at the other point points
        # away from the apex.
        if pos1.magnitude() == 0:
            v = pos2.normalize()
            return (v, -v)
        if pos2.magnitude() == 0:
            v = pos1.normalize()
            return (-v, v)

        # Project the position vectors onto the cone's axis to get the z-height.
        z1 = pos1.dot(self.axis)
        z2 = pos2.dot(self.axis)
        # Calculate the slant height (distance from apex to point along the
        # cone's surface). This is R in polar coordinates for the unrolled cone.
        s1 = z1 / cos(self.alpha)
        s2 = z2 / cos(self.alpha)

        # Calculate the geodesic length (shortest distance on the surface).
        # This will be the denominator when normalizing the tangent vectors.
        L = self.geodesic_length(point_1, point_2)

        # If the points are the same, the tangent vectors are zero vectors.
        if L == 0:
            return (Vector(0), Vector(0))

        # Unrolling the cone
        # Find the unit vectors perpendicular to the cone's axis that point
        # towards each point's projection on the xy-plane.
        n1 = (pos1 - z1 * self.axis).normalize()
        n2 = (pos2 - z2 * self.axis).normalize()

        # Calculate the central angle (theta) between the two points in the
        # plane perpendicular to the cone's axis.
        central = _directional_atan(
            cancel((n1.cross(n2)).dot(self.axis)),
            cancel(n1.dot(n2))
        )

        # The shortest path is chosen by ensuring the angle is not reflex.
        central = Piecewise((central, central <= pi), (2 * pi - central, True))

        # Convert the 3D central angle to the corresponding angle (phi) in the
        # unrolled 2D planar sector.
        delta_u = central * sin(self.alpha)

        # At each point, define an orthogonal basis on the tangent plane.
        # All vectors are expressed in the cone's main reference frame.
        # g: The generator vector, pointing radially away from the apex.
        # c: The circumferential vector, tangential to the circular base.
        g1 = pos1.normalize()
        c1 = self.axis.cross(n1)
        g2 = pos2.normalize()
        c2 = self.axis.cross(n2)

        # In the unrolled 2D plane, the tangent vector is constant. We find its
        # components in the local polar basis at each point.
        # v_radial: Component along the generator vector 'g'.
        # v_circ: Component along the circumferential vector 'c'.

        # Components for v1 (tangent vector at point_1)
        v1_radial_comp = (s2 * cos(delta_u) - s1) / L
        v1_circ_comp = (s2 * sin(delta_u)) / L
        v1 = v1_radial_comp * g1 + v1_circ_comp * c1

        # Components for v2 (tangent vector at point_2)
        v2_radial_comp = (s1 * cos(delta_u) - s2) / L
        v2_circ_comp = (-s1 * sin(delta_u)) / L
        v2 = v2_radial_comp * g2 + v2_circ_comp * c2

        return (v1, v2)

    def __repr__(self):
        """Representation of a ``WrappingCone``."""
        return (
            f'{self.__class__.__name__}(alpha={self.alpha}, '
            f'apex={self.apex}, axis={self.axis})'
        )


def _directional_atan(numerator, denominator):
    """Compute atan in a directional sense as required for geodesics.

    Explanation
    ===========

    To be able to control the direction of the geodesic length along the
    surface of a cylinder a dedicated arctangent function is needed that
    properly handles the directionality of different case. This function
    ensures that the central angle is always positive but shifting the case
    where ``atan2`` would return a negative angle to be centered around
    ``2*pi``.

    Notes
    =====

    This function only handles very specific cases, i.e. the ones that are
    expected to be encountered when calculating symbolic geodesics on uniformly
    curved surfaces. As such, ``NotImplemented`` errors can be raised in many
    cases. This function is named with a leader underscore to indicate that it
    only aims to provide very specific functionality within the private scope
    of this module.

    """

    if numerator.is_number and denominator.is_number:
        angle = atan2(numerator, denominator)
        if angle < 0:
            angle += 2 * pi
    elif numerator.is_number:
        msg = (
            f'Cannot compute a directional atan when the numerator {numerator} '
            f'is numeric and the denominator {denominator} is symbolic.'
        )
        raise NotImplementedError(msg)
    elif denominator.is_number:
        msg = (
            f'Cannot compute a directional atan when the numerator {numerator} '
            f'is symbolic and the denominator {denominator} is numeric.'
        )
        raise NotImplementedError(msg)
    else:
        ratio = sympify(trigsimp(numerator / denominator))
        if isinstance(ratio, tan):
            angle = ratio.args[0]
        elif (
            ratio.is_Mul
            and ratio.args[0] == Integer(-1)
            and isinstance(ratio.args[1], tan)
        ):
            angle = 2 * pi - ratio.args[1].args[0]
        else:
            msg = f'Cannot compute a directional atan for the value {ratio}.'
            raise NotImplementedError(msg)

    return angle
