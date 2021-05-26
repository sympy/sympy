from sympy.physics.mechanics import Particle, RigidBody
from abc import ABC, abstractmethod

class Joint(ABC):
    """ A basic Joint.

    Explanation
    ===========

    Abstract base class for all joints.

    Attributes
    ==========

    name : string
        The Joint's name.
    parent : Particle/RigidBody/Body
        The parent body of joint.
    child : Particle/RigidBody/Body
        The child body of joint.
    coordinates: List
        Coordinates of joint.
    speeds : List
        Speeds of joint.

    Examples
    ========

    >>> from sympy.physics.vector import dynamicsymbols
    >>> from sympy.physics.mechanics import Body
    >>> from sympy.physics.mechanics.joints import Joint

    >>> class J(Joint):
    ...     def _generate_coordinates(self, coordinates):
    ...         if coordinates is None:
    ...             a = dynamicsymbols('a')
    ...             return [a]
    ...         if not isinstance(coordinates, list):
    ...             raise TypeError()
    ...         return coordinates
    ...
    ...     def _generate_speeds(self, speeds):
    ...         if speeds is None:
    ...             b = dynamicsymbols('b')
    ...             return [b]
    ...         if not isinstance(speeds, list):
    ...             raise TypeError()
    ...         return speeds
    ...
    ...     def apply_joint(self):
    ...         t = dynamicsymbols._t
    ...         return [self._coordinates[0].diff(t) - self._speeds[0]]

    >>> P = Body('P')
    >>> C = Body('B')
    >>> J1 = J('J1', P, C)
    >>> J1.kde()
    [-b(t) + Derivative(a(t), t)]

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name

        if not isinstance(parent, (Particle, RigidBody)):
            raise TypeError('Parent must be an instance of Particle/RigidBody.')
        self._parent = parent

        if not isinstance(child, (Particle, RigidBody)):
            raise TypeError('Parent must be an instance of Particle/RigidBody.')
        self._child = child

        self._coordinates = self._generate_coordinates(coordinates)
        self._speeds = self._generate_speeds(speeds)
        self._kde = self.apply_joint()

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.__str__()

    @property
    def parent(self):
        """ Parent body of Joint."""
        return self._parent

    @property
    def child(self):
        """ Child body of Joint."""
        return self._child

    def coordinates(self):
        """ List of Coordinates of Joint."""
        return self._coordinates

    def speeds(self):
        """ List of speeds of Joint."""
        return self._speeds

    def kde(self):
        """ KDE of Joint."""
        return self._kde

    @abstractmethod
    def _generate_coordinates(self, coordinates):
        """ Generate list of coordinates."""
        pass

    @abstractmethod
    def _generate_speeds(self, speeds):
        """ Generate list of speeds. """
        pass

    @abstractmethod
    def apply_joint(self):
        """ Apply the particular Joint."""
        pass
