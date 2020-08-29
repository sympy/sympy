from sympy.core.backend import sympify
from sympy.physics.vector import Point
from sympy import pi
from sympy.physics.units.definitions.unit_definitions import e0

from sympy.utilities.exceptions import SymPyDeprecationWarning

from sympy import symbols

__all__ = ['Point_charge']

class Point_charge(object):
    """A particle.

    Particles have a non-zero mass and lack spatial extension; they take up no
    space.

    Values need to be supplied on initialization, but can be changed later.

    Parameters
    ==========
    name : str
        Name of particle
    point : Point
        A physics/mechanics Point which represents the position, velocity, and
        acceleration of this Particle
    charge : sympifyable
        A SymPy expression representing the Particle's charge

    Examples
    ========

    >>> from sympy.physics.electrostatics import Particle, Point
    >>> from sympy import Symbol
    >>> po = Point('po')
    >>> q = Symbol('q')
    >>> pa = Particle('pa', po, q)
    >>> # Or you could change these later
    >>> pa.charge = q
    >>> pa.point = po

    """

    def __init__(self, name, point, charge):
        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name
        self.charge = charge
        self.point = point

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.__str__()

    @property
    def charge(self):
        """Charge of the particle."""
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = sympify(value)

    @property
    def point(self):
        """Point of the particle."""
        return self._point

    @point.setter
    def point(self, p):
        if not isinstance(p, Point):
            raise TypeError("Particle point attribute must be a Point object.")
        self._point = p

    def value(self, val):
        """Gives numerical value or expression as per user's choice."""
        if val:
            return 9 * 10**9
        else:
            return 1/(4 * pi * e0)

    def is_vector(vector):
        #Converts scalar value into vector
        
        pass

    def electric_potential(self, r, Er = 1, val = False):
        """Electric Potential due to a point charge.

        Electric potential V due to point charge q of Particle P at a distance r
        is given by

        V = k * 1/Er * q / r

        where k = 1/4*pi*e0 , q is charge of particle , r is distance at which
        potential is required.

        Parameters
        ==========

        r : Distance
            Distance at which potential is required.

        Er(optional) : Permittivity of a medium
            Permitivity of the medium in which Electric Potential is required.

        val(optional) : True/False
            Choice to get value or expression of k.

        Examples
        ========

        >>> from sympy.physics.electrostatics import Particle, Point
        >>> from sympy import Symbol
        >>> po = Point('po')
        >>> q, r, E = Symbol('q r E')
        >>> pa = Particle('pa', po, q)
        >>> # Or you could change these later
        >>> pa.electric_potential(r,Er=E)
        q/(4*pi*E*r*Eo)

    """

        return self.value(val) * 1/Er * self.charge /sympify(r)

    def electrostatic_force(self, P, r, Er = 1, val = False):
        """Electrostatic Force on one particle due other.

        Electric Force F due to point charge q of Particle P at a distance r
        is given by

        V = k * q1*q2/ r**2

        where k = 1/4*pi*e0*Er , q is charge of particle , r is distance at which
        potential is required.

        Parameters
        ==========

        P : Insatnce of class Particle
            The point charge witch wich Electrostatics force is required.

        r : Distance
            Distance at which potential is required.

        Er(optional) : Permittivity of a medium
            Permitivity of the medium in which Electric Potential is required.

        val(optional) : True/False
            Choice to get value or expression of k.

        Examples
        ========

        >>> from sympy.physics.electrostatics import Particle, Point
        >>> from sympy import Symbol
        >>> po = Point('po')
        >>> q1 q2 r E = Symbol('q r E')
        >>> pa = Particle('pa', po, q)
        >>> S = Point('S')
        >>> s = Particle('pa2', S, q2)
        >>> pa.electric_potential(s, r, Er=E)
        q1*q2/(4*pi*E*Eo*r**2)

    """

        if not isinstance(P, Point_charge):
            raise TypeError('Supply a Particle instance.')
        return self.value(val) * 1/Er * ( self.charge * P.charge)/r**2

    def electric_field_intensity(self, r, Er = 1, val = False):
        """Electric Field Intensity due to a point charge.

        Electric Field Intensity E due to point charge q of Particle P at a distance r
        is given by

        E = k * q/ r**2

        where k = 1/4*pi*e0*Er , q is charge of particle , r is distance at which
        electric field intensity is required.

        Parameters
        ==========

        P : Insatnce of class Particle
            The point charge witch wich Electrostatics force is required.

        r : Distance
            Distance at which potential is required.

        Er(optional) : Permittivity of a medium
            Permitivity of the medium in which Electric Potential is required.

        val(optional) : True/False
            Choice to get value or expression of k.

        Examples
        ========

        >>> from sympy.physics.electrostatics import Particle, Point
        >>> from sympy import Symbol
        >>> po = Point('po')
        >>> q, r, E = Symbol('q r E')
        >>> pa = Particle('pa', po, q)
        >>> pa.electric_field_intensity(r, Er=E)
        q/(4*pi*E*Eo*r**2)

    """

        return self.value(val) * 1/Er * self.charge/r**2
