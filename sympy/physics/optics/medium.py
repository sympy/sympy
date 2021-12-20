"""
**Contains**

* Medium
"""
from sympy.physics.units import second, meter, kilogram, ampere

__all__ = ['Medium']

from sympy.core.basic import Basic
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.units import speed_of_light, u0, e0
from sympy.printing import sstr


c = speed_of_light.convert_to(meter/second)
_e0mksa = e0.convert_to(ampere**2*second**4/(kilogram*meter**3))
_u0mksa = u0.convert_to(meter*kilogram/(ampere**2*second**2))


class Medium(Basic):

    """
    This class represents an optical medium. The prime reason to implement this is
    to facilitate refraction, Fermat's principle, etc.

    Explanation
    ===========

    An optical medium is a material through which electromagnetic waves propagate.
    The permittivity and permeability of the medium define how electromagnetic
    waves propagate in it.


    Parameters
    ==========

    name: string
        The display name of the Medium.

    permittivity: Sympifyable
        Electric permittivity of the space.

    permeability: Sympifyable
        Magnetic permeability of the space.

    n: Sympifyable
        Index of refraction of the medium.


    Examples
    ========

    >>> from sympy.abc import epsilon, mu
    >>> from sympy.physics.optics import Medium
    >>> m1 = Medium('m1')
    >>> m2 = Medium('m2', epsilon, mu)
    >>> m1.intrinsic_impedance
    149896229*pi*kilogram*meter**2/(1250000*ampere**2*second**3)
    >>> m2.refractive_index
    299792458*meter*sqrt(epsilon*mu)/second


    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Optical_medium

    """

    def __new__(cls, name, permittivity=None, permeability=None, n=None):
        if not isinstance(name, Str):
            name = Str(name)

        if n is not None:
            if permittivity is not None and permeability is None:
                permeability = n**2/(c**2*permittivity)
                return MediumPP(name, _sympify(permittivity), _sympify(permeability))
            elif permeability is not None and permittivity is None:
                permittivity = n**2/(c**2*permeability)
                return MediumPP(name, _sympify(permittivity), _sympify(permeability))
            elif permittivity is not None and permittivity is not None:
                expr = abs(n - c * sqrt(permittivity * permeability))
                expr = expr.subs({meter: 1, second: 1})
                if len(expr.free_symbols) == 0 and expr > 1e-6:
                    raise ValueError("Values are not consistent.")
                return MediumPP(name, _sympify(permittivity), _sympify(permeability))
            else:
                return MediumN(name, _sympify(n))
        elif permittivity is not None and permeability is not None:
            return MediumPP(name, _sympify(permittivity), _sympify(permeability))
        elif permittivity is not None and permeability is None:
            raise ValueError("At least one of n and permeability must not be None if permittivity is not None")
        elif permittivity is None and permeability is not None:
            raise ValueError("At least one of n and permittivity must not be None if permeability is not None")
        elif permittivity is None and permeability is None:
            permittivity = _e0mksa
            permeability = _u0mksa
            return MediumPP(name, _sympify(permittivity), _sympify(permeability))

    @property
    def speed(self):
        """
        Returns speed of the electromagnetic wave travelling in the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.speed
        299792458*meter/second
        >>> m2 = Medium('m2', n=1)
        >>> m.speed == m2.speed
        True

        """
        return c / self.n

    @property
    def refractive_index(self):
        """
        Returns refractive index of the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.refractive_index
        1

        """
        return (c/self.speed)

    def __lt__(self, other):
        """
        Compares based on refractive index of the medium.
        """
        return self.refractive_index < other.refractive_index

    def __gt__(self, other):
        return not self < other

    def __eq__(self, other):
        return self.refractive_index == other.refractive_index

    def __ne__(self, other):
        return not self == other


class MediumN(Medium):

    def __new__(cls, name, n):
        obj = super(Medium, cls).__new__(cls, name, n)
        obj.name = name
        obj._n = n
        return obj

    @property
    def n(self):
        return self._n

    def __str__(self):
        return type(self).__name__ + ': ' + sstr([self.n])


class MediumPP(Medium):

    def __new__(cls, name, permittivity, permeability):
        obj = super(Medium, cls).__new__(cls, name, permittivity, permeability)
        obj.name = name
        obj._permittivity = permittivity
        obj._permeability = permeability
        return obj

    @property
    def intrinsic_impedance(self):
        """
        Returns intrinsic impedance of the medium.

        Explanation
        ===========

        The intrinsic impedance of a medium is the ratio of the
        transverse components of the electric and magnetic fields
        of the electromagnetic wave travelling in the medium.
        In a region with no electrical conductivity it simplifies
        to the square root of ratio of magnetic permeability to
        electric permittivity.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.intrinsic_impedance
        149896229*pi*kilogram*meter**2/(1250000*ampere**2*second**3)

        """
        return sqrt(self._permeability / self._permittivity)

    @property
    def permittivity(self):
        """
        Returns electric permittivity of the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.permittivity
        625000*ampere**2*second**4/(22468879468420441*pi*kilogram*meter**3)

        """
        return self._permittivity

    @property
    def permeability(self):
        """
        Returns magnetic permeability of the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.permeability
        pi*kilogram*meter/(2500000*ampere**2*second**2)

        """
        return self._permeability

    @property
    def n(self):
        return c*sqrt(self.permittivity*self.permeability)

    def __str__(self):
        return type(self).__name__ + ': ' + sstr([self._permittivity,
                self._permeability, self.n])

