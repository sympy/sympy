"""
This module has all the classes and functions related to waves in optics.

**Contains**

* TWave
"""

from __future__ import print_function, division

__all__ = ['TWave']

from sympy import sympify, pi, cos, sqrt, simplify, Symbol, S
from sympy.core.expr import Expr


class TWave(Expr):

    r"""
    This is a simple transverse wave travelling in a two dimensional space.
    Basic properties are required at the time of creation of the object but
    they can be changed later with respective methods provided.

    It has been represented as :math:`A \times cos(\omega \times t + \phi )`
    where :math:`A` is amplitude, :math:`\omega` is angular velocity and
    :math:`\phi` is phase angle of the wave.


    Arguments
    =========

    amplitude : Sympifyable
        Amplitude of the wave.
    frequency : Sympifyable
        Frequency of the wave.
    phase : Sympifyable
        Phase angle of the wave.
    time_period : Sympifyable
        Time period of the wave.
    n : Sympifyable
        Refractive index of the medium.

    Raises
    =======

    ValueError : When niether frequency nor time period is provided.
    TypeError : When anyting other than TWave objects is added.


    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.optics import TWave
    >>> A1, phi1, A2, phi2, f = symbols('A1, phi1, A2, phi2, f')
    >>> w1 = TWave(A1, f, phi1)
    >>> w2 = TWave(A2, f, phi2)
    >>> w3 = w1 + w2  # Superposition of two waves
    >>> w3
    TWave(sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2), f, phi1 + phi2)
    >>> w3.amplitude
    sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2)
    >>> w3.phase
    phi1 + phi2
    >>> w3.speed
    c/n
    >>> w3.angular_velocity
    2*pi*f

    """

    def __init__(
            self,
            amplitude,
            frequency=None,
            phase=S.Zero,
            time_period=None,
            n=Symbol('n')):
        frequency = sympify(frequency)
        amplitude = sympify(amplitude)
        phase = sympify(phase)
        time_period = sympify(time_period)
        n = sympify(n)
        self._frequency = frequency
        self._amplitude = amplitude
        self._phase = phase
        self._time_period = time_period
        self._n = n
        self.c = Symbol('c')  # Speed of light in vacuum
        if time_period is not None:
            self._frequency = 1/self._time_period
        if frequency is not None:
            self._time_period = 1/self._frequency
            if time_period is not None:
                if frequency != 1/time_period:
                    raise ValueError("frequency and time_period should be consistent.")
        if frequency is None and time_period is None:
            raise ValueError("Either frequency or time period is needed.")

    @property
    def frequency(self):
        """
        Returns the frequency of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.frequency
        f
        """
        return self._frequency

    @property
    def time_period(self):
        """
        Returns the time period of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.time_period
        1/f
        """
        return self._time_period

    @property
    def wavelength(self):
        """
        Returns wavelength of the wave.
        It depends on the medium of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.wavelength
        c/(f*n)
        """
        return self.c/(self._frequency*self._n)

    @property
    def amplitude(self):
        """
        Returns the amplitude of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.amplitude
        A
        """
        return self._amplitude

    @property
    def phase(self):
        """
        Returns the phase angle of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.phase
        phi
        """
        return self._phase

    @property
    def speed(self):
        """
        Returns the speed of travelling wave.
        It is medium dependent.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.speed
        c/n
        """
        return self.wavelength*self._frequency

    @property
    def angular_velocity(self):
        """
        Returns angular velocity of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.angular_velocity
        2*pi*f
        """
        return 2*pi*self._frequency

    def equation(self, type='cosine'):
        """
        Returns equation of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.equation('cosine')
        A*cos(2*pi*f*t + phi)
        """
        if not isinstance(type, str):
            raise TypeError("type can only be a string.")
        if type == 'cosine':
            return self._amplitude*cos(self.angular_velocity*Symbol('t') + self._phase)

    def __str__(self):
        """String representation of a TWave."""
        from sympy.printing import sstr
        return type(self).__name__ + sstr(self.args)

    __repr__ = __str__

    def __add__(self, other):
        """
        Addition of two waves will result in their superposition.
        The type of interference will depend on their phase angles.
        """
        if isinstance(other, TWave):
            if self._frequency == other._frequency and self.wavelength == other.wavelength:
                return TWave(sqrt(self._amplitude**2 + other._amplitude**2 + 2 *
                                  self.amplitude*other.amplitude*cos(
                                      self._phase - other.phase)),
                             self.frequency,
                             self._phase + other._phase
                             )
        else:
            raise TypeError(type(other).__name__ + " and TWave objects can't be added.")
