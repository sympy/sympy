__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.essential import ReferenceFrame, Dyadic

class RigidBody(object):
    """An idealized rigid body.

    This is essentially a container which holds the various components which
    describe a rigid body: a name, mass, center of mass, reference frame, and
    inertia.

    All of these need to be supplied on creation, but can be changed
    afterwards.

    Attributes
    ==========
    name : string
        The body's name.
    masscenter : Point
        The point which represents the center of mass of the rigid body.
    frame : ReferenceFrame
        The ReferenceFrame which the rigid body is fixed in.
    mass : Sympifyable
        The body's mass.
    inertia : (Dyadic, Point)
        The body's inertia about a point; stored in a tuple as shown above.

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
    >>> from sympy.physics.mechanics import outer
    >>> m = Symbol('m')
    >>> A = ReferenceFrame('A')
    >>> P = Point('P')
    >>> I = outer (A.x, A.x)
    >>> inertia_tuple = (I, P)
    >>> B = RigidBody('B', P, A, m, inertia_tuple)
    >>> # Or you could change them afterwards
    >>> m2 = Symbol('m2')
    >>> B.mass = m2

    """

    def __init__(self, name, masscenter, frame, mass, inertia):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name
        self.set_masscenter(masscenter)
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

    def get_masscenter(self):
        return self._masscenter

    def set_masscenter(self, p):
        if not isinstance(p, Point):
            raise TypeError("RigidBody center of mass must be a Point object.")
        self._masscenter = p

    masscenter = property(get_masscenter, set_masscenter)

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

    def linmom(self, frame):
        """ Linear momentum of a rigid body.

        The linear mometum, L, of a rigid body, B, of mass, M, whose mass
        center, B*, is translating with an inertial velocity, ^N v^B* is given
        by:

        L = M * ^N v^B*

        Parameters
        ==========

        frame : ReferenceFrame
            Linear momentum is defined in an inertial reference frame; the
            user is responsible for entering the correct reference frame.

        Examples
        ========

        >>> from sympy.physics.mechanics import Point, ReferenceFrame, outer
        >>> from sympy.physics.mechanics import RigidBody, dynamicsymbols
        >>> M, v = dynamicsymbols('M v')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> P.set_vel(N, v * N.x)
        >>> I = outer (N.x, N.x)
        >>> Inertia_tuple = (I, P)
        >>> B = RigidBody('B', P, N, M, Inertia_tuple)
        >>> B.linmom(N)
        M*v*N.x

        """

        return self.mass * self.masscenter.vel(frame)

    inertia = property(get_inertia, set_inertia)
