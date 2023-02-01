from abc import ABC, abstractmethod
from collections import namedtuple
from sympy.physics.mechanics.body_base import BodyBase
from sympy.physics.vector import Vector, ReferenceFrame, Point

__all__ = ['LoadBase', 'Force', 'Torque']


class LoadBase(ABC, namedtuple('LoadBase', ['location', 'vector'])):
    """Abstract base class for the various loading types."""

    def __new__(cls, location, vector):
        return super().__new__(cls, cls._generate_location(location),
                               cls._generate_vector(vector))

    @staticmethod
    @abstractmethod
    def _generate_location(location):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def _generate_vector(vector):
        raise NotImplementedError

    def __add__(self, other):
        raise TypeError(f"unsupported operand type(s) for +: "
                        f"'{self.__class__.__name__}' and "
                        f"'{other.__class__.__name__}'")

    def __mul__(self, other):
        raise TypeError(f"unsupported operand type(s) for *: "
                        f"'{self.__class__.__name__}' and "
                        f"'{other.__class__.__name__}'")

    __radd__ = __add__
    __rmul__ = __mul__


class Force(LoadBase):
    """Force acting upon a point.

    Explanation
    ===========

    A force is a vector that is bound to a line of action. This class stores
    both a point, which lies on the line of action, and the vector. A tuple can
    also be used, with the location as the first entry and the vector as second
    entry.

    Examples
    ========

    A force of magnitude 2 along N.x acting on a point Po can be created as
    follows:

    >>> from sympy.physics.mechanics import Point, ReferenceFrame, Force
    >>> N = ReferenceFrame('N')
    >>> Po = Point('Po')
    >>> Force(Po, 2 * N.x)
    Force(point=Po, force=2*N.x)

    If a body is supplied, then the center of mass of that body is used.

    >>> from sympy.physics.mechanics import Particle
    >>> P = Particle('P', point=Po)
    >>> Force(P, 2 * N.x)
    Force(point=Po, force=2*N.x)

    """

    def __new__(cls, point, force):
        return super().__new__(cls, point, force)

    def __repr__(self):
        return f'Force(point={self.point}, force={self.force})'

    @staticmethod
    def _generate_location(point):
        if isinstance(point, BodyBase):
            return point.masscenter
        if isinstance(point, Point):
            return point
        raise TypeError('Force location should of type Point.')

    @staticmethod
    def _generate_vector(force):
        if isinstance(force, Vector):
            return force
        raise TypeError('Force vector should be of type Vector.')

    @property
    def point(self):
        return self.location

    @property
    def force(self):
        return self.vector


class Torque(LoadBase):
    """Force acting upon a point.

    Explanation
    ===========

    A torque is a free vector that is acting on a reference frame, which is
    associated with a rigid body. This class stores both the frame and the
    vector. A tuple can also be used, with the location as the first entry and
    the vector as second entry.

    Examples
    ========

    A torque of magnitude 2 along N.x acting on a frame N can be created as
    follows:

    >>> from sympy.physics.mechanics import ReferenceFrame, Torque
    >>> N = ReferenceFrame('N')
    >>> Torque(N, 2 * N.x)
    Torque(frame=N, torque=2*N.x)

    If a body is supplied, then the frame fixed to that body is used.

    >>> from sympy.physics.mechanics import RigidBody
    >>> rb = RigidBody('rb', frame=N)
    >>> Torque(rb, 2 * N.x)
    Torque(frame=N, torque=2*N.x)

    """

    def __new__(cls, frame, torque):
        return super().__new__(cls, frame, torque)

    def __repr__(self):
        return f'Torque(frame={self.frame}, torque={self.torque})'

    @staticmethod
    def _generate_location(frame):
        if isinstance(frame, BodyBase):
            return frame.frame
        if isinstance(frame, ReferenceFrame):
            return frame
        raise TypeError('Torque location should of type ReferenceFrame.')

    @staticmethod
    def _generate_vector(torque):
        if isinstance(torque, Vector):
            return torque
        raise TypeError('Torque vector should be of type Vector.')

    @property
    def frame(self):
        return self.location

    @property
    def torque(self):
        return self.vector


def gravity(acceleration, *bodies):
    """
    Returns a list of gravity forces given the acceleration
    due to gravity and any number of particles or rigidbodies.

    Example
    =======

    >>> from sympy.physics.mechanics import ReferenceFrame, Particle, RigidBody
    >>> from sympy.physics.mechanics.loads import gravity
    >>> from sympy import symbols
    >>> N = ReferenceFrame('N')
    >>> g = symbols('g')
    >>> P = Particle('P')
    >>> B = RigidBody('B')
    >>> gravity(g*N.y, P, B)
    [Force(point=P_masscenter, force=P_mass*g*N.y),
     Force(point=B_masscenter, force=B_mass*g*N.y)]

    """

    gravity_force = []
    for body in bodies:
        if not isinstance(body, BodyBase):
            raise TypeError(f'{type(body)} is not a body type')
        gravity_force.append(Force(body.masscenter, body.mass * acceleration))
    return gravity_force
