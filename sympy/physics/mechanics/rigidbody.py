__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.essential import ReferenceFrame, Dyadic

class RigidBody(object):
    """An idealized rigid body.

    This is essentially a container which holds the various components which
    describe a rigid body: mass, point, reference frame, and inertia.

    Attributes
    ==========

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
    >>> B = RigidBody()
    >>> B.mass = m
    >>> B.frame = A
    >>> B.mc = P
    >>> B.inertia = (I, B.mc)

    """

    def __init__(self):
        self._frame = None
        self._mc = None
        self._mass = None
        self._inertia = None
        self._inertia_point = None

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
            raise TypeError("RigidBody mass center must be a Point object")
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
