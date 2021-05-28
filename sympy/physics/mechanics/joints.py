from sympy.physics.vector.functions import dynamicsymbols
from sympy.physics.mechanics import Body
from sympy.physics.vector import Vector, Point
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
    parent_load_point : Point
        Point about which load is applied(masscenter)
    child_load_point : Point
        Point about which load is applied(masscenter)
    child_axis : Vector
        Axis of Child frame.
    parent_axis : Vector
        Axis of parent frame about which child would be oriented.

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
    >>> C = Body('C')
    >>> J1 = J('J1', P, C)
    >>> J1.kdes()
    [-b(t) + Derivative(a(t), t)]
    >>> J1.parent()
    P
    >>> J1.child()
    C

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_load_point=None, 
        child_load_point=None, parent_axis = None, child_axis=None):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name

        if not isinstance(parent, Body):
            raise TypeError('Parent must be an instance of Particle/RigidBody.')
        self._parent = parent

        if not isinstance(child, Body):
            raise TypeError('Parent must be an instance of Particle/RigidBody.')
        self._child = child

        self._coordinates = self._generate_coordinates(coordinates)
        self._speeds = self._generate_speeds(speeds)
        self._kdes = self._generate_kdes()

        self.child_axis = self._axis(child, child_axis)
        self.parent_axis = self._axis(parent, parent_axis)

        self.parent_load_pos = self._generate_load_point(parent,parent_load_point)
        self.child_load_pos = self._generate_load_point(parent, child_load_point)

        self._orient_frames()
        self._set_angular_velocity()
        self._set_linear_velocity()

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.__str__()

    def parent(self):
        """ Parent body of Joint."""
        return self._parent

    def child(self):
        """ Child body of Joint."""
        return self._child

    def coordinates(self):
        """ List of Coordinates of Joint."""
        return self._coordinates

    def speeds(self):
        """ List of speeds of Joint."""
        return self._speeds

    def kdes(self):
        """ KDE of Joint."""
        return self._kdes

    @abstractmethod
    def _generate_coordinates(self, coordinates):
        """ Generate list of coordinates."""
        pass

    @abstractmethod
    def _generate_speeds(self, speeds):
        """ Generate list of speeds. """
        pass

    @abstractmethod
    def _orient_frames(self):
        """Orient frames as per the joint"""
        pass
    
    @abstractmethod
    def _set_angular_velocity(self):
        pass

    @abstractmethod
    def _set_linear_velocity(self):
        pass

    def _generate_kdes(self):
        kdes = []
        t = dynamicsymbols._t
        for i in range(len(self._coordinates)):
            kdes.append(-self._coordinates[0].diff(t) + self._speeds[0])
    
    def _axis(self, body, ax):
        if ax is None:
            ax = body.frame.x
            return ax
        if not isinstance(ax, Vector):
            raise TypeError("Axis must be of type Vector. Example-> A.x wehere 'A' is ReferenceFrame.")
        return ax
        
    def _generate_load_point(self, body, load_point):
        if load_point is None:
            return body.masscenter
        if not isinstance(load_point, Point):
            raise ValueError
        return load_point
        
