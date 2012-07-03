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

    inertia = property(get_inertia, set_inertia)

    def linearmomentum(self, frame):
        """ Linear momentum of the rigid body.

        The linear momentum L, of a rigid body B, with respect to frame N is
        given by

        L = M * v*

        where M is the mass of the rigid body and v* is the velocity of
        the mass center of B in the frame, N.

        Parameters
        ==========

        frame : ReferenceFrame
            The frame in which linear momentum is desired.

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
        >>> B.linearmomentum(N)
        M*v*N.x

        """

        return self.mass * self.masscenter.vel(frame)

    def angularmomentum(self, point, frame):
        """ Angular momentum of the rigid body.

        The angular momentum H, about some point O, of a rigid body B, in a
        frame N is given by

        H = I* . omega + r* x (M * v)

        where I* is the central inertia dyadic of B, omega is the angular
        velocity of body B in the frame, N, r* is the position vector from
        point O to the mass center of B, and v is the velocity of point O in
        the frame, N.

        Parameters
        ==========

        point : Point
            The point about which angular momentum is desired.

        frame : ReferenceFrame
            The frame in which angular momentum is desired.

        Examples
        ========

        >>> from sympy.physics.mechanics import Point, ReferenceFrame, outer
        >>> from sympy.physics.mechanics import RigidBody, dynamicsymbols
        >>> M, v, r, omega = dynamicsymbols('M v r omega')
        >>> N = ReferenceFrame('N')
        >>> b = ReferenceFrame('b')
        >>> b.set_ang_vel(N, omega * b.x)
        >>> P = Point('P')
        >>> P.set_vel(N, 1 * N.x)
        >>> I = outer (b.x, b.x)
        >>> Inertia_tuple = (I, P)
        >>> B = RigidBody('B', P, b, M, Inertia_tuple)
        >>> B.angularmomentum(P, N)
        omega*b.x

        """

        return ((self.inertia[0] & self.frame.ang_vel_in(frame)) +
                (point.vel(frame) ^ -self.masscenter.pos_from(point)) *
                self.mass)
