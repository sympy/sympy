from __future__ import print_function, division

__all__ = ['RigidBody']

from sympy import sympify
from sympy.physics.vector import Point, ReferenceFrame, Dyadic


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
        self._pe = sympify(0)

    def __str__(self):
        return self._name

    __repr__ = __str__

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        if isinstance(parent, RigidBody):
            self._parent = parent
        else:
            raise TypeError("parent of {0} should be a Rigidbody".format(self._name))

    @property
    def child(self):
        return self.child

    @child.setter
    def child(self, child):
        if isinstance(child, RigidBody):
            self._child = child
        else:
            raise TypeError("child of {0} should be a Rigidbody".format(self._name))

    @property
    def system(self):
        return self.system

    @system.setter
    def system(self, system):
        # TODO Checking if this is an instance of pydy.System will include
        # pydy as an dependency. How to deal with this?
        self.coordinates = []
        self.diff_coordinates = []
        self.velocities = []
        self.kdequations = []
        self.force_list = []
        self._system = system

    def add_coordinate(self, coordinate):
        self.coordinates.append(coordinate)

    def add_diff_velocity(self, diff_velocity):
        self.diff_coordinates.append(diff_velocity)

    def add_velocity(self, velocity):
        self.velocities.append(velocity)

    def add_kdequation(self, kdequation):
        self.kdequations.append(kdequation)

    def add_force(self, force):
        self.force_list.append(force)

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
        # have I S/O, want I S/S*
        # I S/O = I S/S* + I S*/O; I S/S* = I S/O - I S*/O
        # I_S/S* = I_S/O - I_S*/O
        from sympy.physics.mechanics.functions import inertia_of_point_mass
        I_Ss_O = inertia_of_point_mass(self.mass,
                                       self.masscenter.pos_from(I[1]),
                                       self.frame)
        self._central_inertia = I[0] - I_Ss_O

    inertia = property(get_inertia, set_inertia)

    @property
    def central_inertia(self):
        """The body's central inertia dyadic."""
        return self._central_inertia

    def linear_momentum(self, frame):
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
        >>> B.linear_momentum(N)
        M*v*N.x

        """

        return self.mass * self.masscenter.vel(frame)

    def angular_momentum(self, point, frame):
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
        >>> B.angular_momentum(P, N)
        omega*b.x

        """

        return ((self.central_inertia & self.frame.ang_vel_in(frame)) +
                (point.vel(frame) ^ -self.masscenter.pos_from(point)) *
                self.mass)

    def kinetic_energy(self, frame):
        """Kinetic energy of the rigid body

        The kinetic energy, T, of a rigid body, B, is given by

        'T = 1/2 (I omega^2 + m v^2)'

        where I and m are the central inertia dyadic and mass of rigid body B,
        respectively, omega is the body's angular velocity and v is the
        velocity of the body's mass center in the supplied ReferenceFrame.

        Parameters
        ==========

        frame : ReferenceFrame
            The RigidBody's angular velocity and the velocity of it's mass
            center are typically defined with respect to an inertial frame but
            any relevant frame in which the velocities are known can be supplied.

        Examples
        ========

        >>> from sympy.physics.mechanics import Point, ReferenceFrame, outer
        >>> from sympy.physics.mechanics import RigidBody
        >>> from sympy import symbols
        >>> M, v, r, omega = symbols('M v r omega')
        >>> N = ReferenceFrame('N')
        >>> b = ReferenceFrame('b')
        >>> b.set_ang_vel(N, omega * b.x)
        >>> P = Point('P')
        >>> P.set_vel(N, v * N.x)
        >>> I = outer (b.x, b.x)
        >>> inertia_tuple = (I, P)
        >>> B = RigidBody('B', P, b, M, inertia_tuple)
        >>> B.kinetic_energy(N)
        M*v**2/2 + omega**2/2

        """

        rotational_KE = (self.frame.ang_vel_in(frame) & (self.central_inertia &
                self.frame.ang_vel_in(frame)) / sympify(2))

        translational_KE = (self.mass * (self.masscenter.vel(frame) &
            self.masscenter.vel(frame)) / sympify(2))

        return rotational_KE + translational_KE

    def set_potential_energy(self, scalar):
        """Used to set the potential energy of this RigidBody.

        Parameters
        ==========

        scalar: Sympifyable
            The potential energy (a scalar) of the RigidBody.

        Examples
        ========

        >>> from sympy.physics.mechanics import Particle, Point, outer
        >>> from sympy.physics.mechanics import RigidBody, ReferenceFrame
        >>> from sympy import symbols
        >>> b = ReferenceFrame('b')
        >>> M, g, h = symbols('M g h')
        >>> P = Point('P')
        >>> I = outer (b.x, b.x)
        >>> Inertia_tuple = (I, P)
        >>> B = RigidBody('B', P, b, M, Inertia_tuple)
        >>> B.set_potential_energy(M * g * h)

        """

        self._pe = sympify(scalar)

    @property
    def potential_energy(self):
        """The potential energy of the RigidBody.

        Examples
        ========

        >>> from sympy.physics.mechanics import RigidBody, Point, outer, ReferenceFrame
        >>> from sympy import symbols
        >>> M, g, h = symbols('M g h')
        >>> b = ReferenceFrame('b')
        >>> P = Point('P')
        >>> I = outer (b.x, b.x)
        >>> Inertia_tuple = (I, P)
        >>> B = RigidBody('B', P, b, M, Inertia_tuple)
        >>> B.set_potential_energy(M * g * h)
        >>> B.potential_energy
        M*g*h

        """

        return self._pe
