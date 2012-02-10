__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.essential import ReferenceFrame, Dyadic

class RigidBody(object):
    """An idealized rigid body.

    This is essentially a container which holds the various components which
    describe a rigid body: a name, mass, point, reference frame, and inertia.

    All of these need to be supplied on creation, but can be changed
    afterwards.

    Attributes
    ==========
    name : str
        The body's name
    mass : Sympifyable
        The body's mass
    inertia : (Dyadic, Point)
        The body's inertia about a point; stored in a tuple as shown above
    mc : Point
        The point which represents the mass center of the rigid body
    frame : ReferenceFrame
        The ReferenceFrame which the rigid body is fixed in

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
    >>> from sympy.physics.mechanics import outer
    >>> m = Symbol('m')
    >>> A = ReferenceFrame('A')
    >>> P = Point('P')
    >>> I = outer (A.x, A.x)
    >>> Inertia_tuple = (I, P)
    >>> B = RigidBody('B', P, A, m, Inertia_tuple)
    >>> # Or you could change them afterwards
    >>> m2 = Symbol('m2')
    >>> B.mass = m2

    """

    def __init__(self, name, mc, frame, mass, inertia):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name
        self.set_mc(mc)
        self.set_mass(mass)
        self.set_frame(frame)
        self.set_inertia(inertia)

    def __str__(self):
        return self._name

    __repr__ = __str__


    def get_frame(self):
        return self._frame

    def set_frame(self, F):
        if not isinstance(F, ReferenceFrame):
            raise TypeError("RigdBody frame must be a ReferenceFrame object.")
        self._frame = F

    frame = property(get_frame, set_frame)

    def get_mc(self):
        return self._mc

    def set_mc(self, p):
        if not isinstance(p, Point):
            raise TypeError("RigidBody mass center must be a Point object.")
        self._mc = p

    mc = property(get_mc, set_mc)

    def get_mass(self):
        return self._mass

    def set_mass(self, m):
        self._mass = sympify(m)

    mass = property(get_mass, set_mass)

    def get_inertia(self):
        return (self._inertia, self._inertia_point)

    def set_inertia(self, I):
        if not isinstance(I[0], Dyadic):
            raise TypeError("RigidBody inertia must be a Dyadic object.")
        if not isinstance(I[1], Point):
            raise TypeError("RigidBody inertia must be about a Point.")
        self._inertia = I[0]
        self._inertia_point = I[1]

    inertia = property(get_inertia, set_inertia)
