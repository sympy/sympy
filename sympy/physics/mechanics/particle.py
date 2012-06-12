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

    def linmom(self, frame):
        """Linear momentum of a particle.

        The linear momentum, L, of a particle, P, of mass, m, with an inertial
        velocity, ^N v^P is given by:

        L = m * ^N v^P

        Parameters
        ==========

        frame : ReferenceFrame
            Linear momentum is defined in an inertial reference frame; the
            user is responsible for entering the correct reference frame.

        Examples
        ========

        >>> from sympy.physics.mechanics import Particle, Point, ReferenceFrame
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> m, v = dynamicsymbols('m v')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> A = Particle('A', P, m)
        >>> P.set_vel(N, v * N.x)
        >>> print A.linmom(N)
        m*v*N.x

        """

        return self.mass * self.point.vel(frame)

    def angmom(self, point, frame):
        """Angular momentum of a particle.

        The angular momentum, H, about some point O of a particle, P, of mass,
        m, with an inertial velocity, ^N v^P is given by:

        H = r^OP x m * ^N v^P

        Parameters
        ==========

        point : Point
            The point about which angular momentum of the particle is being
            determined.

        frame : ReferenceFrame
            Linear momentum is defined in an inertial reference frame; the
            user is responsible for entering the correct reference frame.

        Examples
        ========

        >>> from sympy.physics.mechanics import Particle, Point, ReferenceFrame
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> m, v, r = dynamicsymbols('m v r')
        >>> N = ReferenceFrame('N')
        >>> O = Point('O')
        >>> A = O.locatenew('A', r * N.x)
        >>> P = Particle('P', A, m)
        >>> P.set_vel(N, v * N.y)
        >>> print P.angmom(O, N)
        m*r*v*N.z

        """

        return self.point.pos_from(point) ^ self.mass * self.point.vel(frame)


    point = property(get_point, set_point)
