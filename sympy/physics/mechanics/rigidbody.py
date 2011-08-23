__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.essential import ReferenceFrame, Dyad

class RigidBody(object):
    """An idealized rigid body.

    This is essentially a container which holds the various components which
    describe a rigid body: mass, point, reference frame, and inertia.

    Attributes
    ==========
    mass : Sympifyable
        The body's mass
    inertia : (Dyad, Point)
        The body's inertia about a point; stored in a tuple as shown above
    mc : Point
        The point which represents the mass center of the rigid body
    frame : ReferenceFrame
        The ReferenceFrame which the rigid body is fixed in

    Example
    =======

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

    @property
    def frame(self):
        return self._frame

    @frame.setter
    def frame(self, F):
        if not isinstance(F, ReferenceFrame):
            raise TypeError("RigdBody frame must be a ReferenceFrame object.")
        self._frame = F

    @property
    def mc(self):
        return self._mc

    @mc.setter
    def mc(self, p):
        if not isinstance(p, Point):
            raise TypeError("RigidBody mass center must be a Point object")
        self._mc = p

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, m):
        self._mass = sympify(m)

    @property
    def inertia(self):
        return (self._inertia, self._inertia_point)

    @inertia.setter
    def inertia(self, I):
        if not isinstance(I[0], Dyad):
            raise TypeError("RigidBody inertia must be a Dyad object.")
        if not isinstance(I[1], Point):
            raise TypeError("RigidBody inertia must be about a Point.")
        self._inertia = I[0]
        self._inertia_point = I[1]

