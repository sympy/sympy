__all__ = ['Particle']

from sympy import sympify
from sympy.physics.mechanics.point import Point

class Particle(object):
    """A particle.

    Particles have a non-zero mass and lack spatial extension; they take up no
    space.

    Parameters
    ==========
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
    >>> pa = Particle()
    >>> m = Symbol('m')
    >>> pa.mass = m
    >>> pa.point = po

    """

    def __init__(self):
        self._mass = None
        self._point = None

    def get_mass(self):
        """Mass of the particle.

        >>> from sympy.physics.mechanics import Particle, Point
        >>> from sympy import Symbol
        >>> pa = Particle()
        >>> m = Symbol('m')
        >>> pa.mass = m
        >>> pa.get_mass()
        m

        """
        return self._mass

    def set_mass(self, mass):
        """Set the mass of a particle

        >>> from sympy.physics.mechanics import Particle, Point
        >>> from sympy import Symbol
        >>> pa = Particle()
        >>> m = Symbol('m')
        >>> pa.set_mass(m)
        >>> pa.get_mass()
        m

        """
        self._mass = sympify(mass)

    mass = property(get_mass, set_mass)

    def get_point(self):
        """Point of the particle.

        >>> from sympy.physics.mechanics import Particle, Point
        >>> from sympy import Symbol
        >>> pa = Particle()
        >>> po = Point('po')
        >>> pa.point = po
        >>> pa.get_point()
        po

        """
        return self._point

    def set_point(self, p):
        """Set the point of a particle

        >>> from sympy.physics.mechanics import Particle, Point
        >>> from sympy import Symbol
        >>> pa = Particle()
        >>> po = Point('po')
        >>> pa.set_point(po)
        >>> pa.get_point()
        po

        """
        if not isinstance(p, Point):
            raise TypeError("Particle point attribute must be a Point object.")
        self._point = p

    point = property(get_point, set_point)
