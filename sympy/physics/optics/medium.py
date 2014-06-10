"""
**Contains**

* Medium
"""

from __future__ import division

__all__ = ['Medium']

from sympy import Symbol, sympify, sqrt
from sympy.physics.units import c, u0, e0

class Medium(Symbol):

    """
    This class represents an optical medium. The prime reason to implement this is
    to facilitate refraction, Fermat's priciple, etc.

    An optical medium is a material through which electromagnetic waves propagate.
    The permittivity and permeability of the medium define how electromagnetic
    waves propagate in it.


    Prameters
    =========

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
    149896229*pi*kg*m**2/(1250000*A**2*s**3)
    >>> m2.refractive_index
    299792458*m*sqrt(epsilon*mu)/s


    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Optical_medium

    """

    def __new__(cls, name, permittivity=e0, permeability=u0, n=None):
        obj = super(Medium, cls).__new__(cls, name)
        obj._permittivity = sympify(permittivity)
        obj._permeability = sympify(permeability)
        obj._n = sympify(n)
        if n is not None:
            if permittivity != e0 and permeability == u0:
                obj._permeability = n**2/(c**2*obj._permittivity)
            if permeability != u0 and permittivity == e0:
                obj._permittivity = n**2/(c**2*obj._permeability)
            # XXX: There's issue with precision. Values may be
            # different slightly.
            #if permittivity != u0 and permittivity != e0:
                # if n != c*sqrt(permittivity*permeability):
                #    raise ValueError("Values are not consistent.")
        else:
            obj._n = c*sqrt(permittivity*permeability)
        if n is None and permittivity == e0 and permeability == u0:
            obj._flag = False
        else:
            obj._flag = True
        return obj

    @property
    def intrinsic_impedance(self):
        """
        Returns intrinsic impedance of the medium.

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
        149896229*pi*kg*m**2/(1250000*A**2*s**3)

        """
        return sqrt(self._permeability/self._permittivity)

    @property
    def speed(self):
        """
        Returns speed of the electromagnetic wave travelling in the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.speed
        299792458*m/s

        """
        return 1/sqrt(self._permittivity*self._permeability)

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
        return c/self.speed

    @property
    def permittivity(self):
        """
        Returns electric permittivity of the medium.

        Examples
        ========

        >>> from sympy.physics.optics import Medium
        >>> m = Medium('m')
        >>> m.permittivity
        625000*A**2*s**4/(22468879468420441*pi*kg*m**3)

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
        pi*kg*m/(2500000*A**2*s**2)

        """
        return self._permeability

    def __str__(self):
        from sympy.printing import sstr
        return type(self).__name__ + sstr(self.args)

    def __lt__(self, other):
        """
        Compares based on refractive index of the medium.
        """
        return self.refractive_index < other.refractive_index

    def __gt__(self, other):
        return not self.__lt__(other)

    def __eq__(self, other):
        return self.refractive_index == other.refractive_index

    def __ne__(self, other):
        return not self.__eq__(other)
