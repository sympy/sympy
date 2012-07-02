__all__ = ['Particle']

from sympy import sympify
from sympy.physics.mechanics.point import Point

class Particle(object):
    """A particle.

    Particles have a non-zero mass and lack spatial extension; they take up no
    space.

    Values need to be supplied on initialization, but can be changed later.

    Parameters
    ==========
    name : str
        Name of particle
    mass : sympifyable
        A SymPy expression representing the Particle's mass
    point : Point
        A physics/mechanics Point which represents the position, velocity, and
        acceleration of this Particle

    Examples
    ========

    >>> from sympy.physics.mechanics import Particle, Point
    >>> from sympy import Symbol
    >>> po = Point('po')
    >>> m = Symbol('m')
    >>> pa = Particle('pa', po, m)
    >>> # Or you could change these later
    >>> pa.mass = m
    >>> pa.point = po

    """

    def __init__(self, name, point, mass):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name
        self.set_mass(mass)
        self.set_point(point)

    def __str__(self):
        return self._name

    __repr__ = __str__

    def get_mass(self):
        """Mass of the particle."""
        return self._mass

    def set_mass(self, mass):
        self._mass = sympify(mass)

    mass = property(get_mass, set_mass)

    def get_point(self):
        """Point of the particle."""
        return self._point

    def set_point(self, p):
        if not isinstance(p, Point):
            raise TypeError("Particle point attribute must be a Point object.")
        self._point = p

    point = property(get_point, set_point)

    def linearmomentum(self, frame):
        """Linear momentum of the particle.

        The linear momentum L, of a particle P, with respect to frame N is
        given by

        L = m * v

        where m is the mass of the particle, and v is the velocity of the
        particle in the frame N.

        Parameters
        ==========

        frame : ReferenceFrame
            The frame in which linear momentum is desired.

        Examples
        ========

        >>> from sympy.physics.mechanics import Particle, Point, ReferenceFrame
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> m, v = dynamicsymbols('m v')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> A = Particle('A', P, m)
        >>> P.set_vel(N, v * N.x)
        >>> print A.linearmomentum(N)
        m*v*N.x

        """

        return self.mass * self.point.vel(frame)

    def angularmomentum(self, point, frame):
        """Angular momentum of the particle about the point.

        The angular momentum H, about some point O of a particle, P, is given
        by:

        H = r x m * v

        where r is the position vector from point O to the particle P, m is
        the mass of the particle, and v is the velocity of the particle in
        the inertial frame, N.

        Parameters
        ==========

        point : Point
            The point about which angular momentum of the particle is desired.

        frame : ReferenceFrame
            The frame in which angular momentum is desired.

        Examples
        ========

        >>> from sympy.physics.mechanics import Particle, Point, ReferenceFrame
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> m, v, r = dynamicsymbols('m v r')
        >>> N = ReferenceFrame('N')
        >>> O = Point('O')
        >>> A = O.locatenew('A', r * N.x)
        >>> P = Particle('P', A, m)
        >>> P.point.set_vel(N, v * N.y)
        >>> print P.angularmomentum(O, N)
        m*r*v*N.z

        """

        return self.point.pos_from(point) ^ (self.mass * self.point.vel(frame))
