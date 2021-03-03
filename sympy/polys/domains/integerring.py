"""Implementation of :class:`IntegerRing` class. """


from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.domains.ring import Ring
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.utilities import public

import math

@public
class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    r"""The domain ``ZZ`` representing the integers `\mathbb{Z}`.

    The :py:class:`IntegerRing` class represents the ring of integers as a
    :py:class:`~.Domain` in the domain system. :py:class:`IntegerRing` is a
    super class of :py:class:`PythonIntegerRing` and
    :py:class:`GMPYIntegerRing` one of which will be the implementation for
    :ref:`ZZ` depending on whether or not ``gmpy`` or ``gmpy2`` is installed.

    See also
    ========

    Domain
    """

    rep = 'ZZ'

    is_IntegerRing = is_ZZ = True
    is_Numerical = True
    is_PID = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def get_field(self):
        r"""Return the associated field of fractions :ref:`QQ`

        Returns
        =======

        :ref:`QQ`:
            The associated field of fractions :ref:`QQ`, a
            :py:class:`~.Domain` representing the rational numbers
            `\mathbb{Q}`.

        Examples
        ========

        >>> from sympy import ZZ
        >>> ZZ.get_field()
        QQ
        """
        from sympy.polys.domains import QQ
        return QQ

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`.

        Parameters
        ==========

        *extension: One or more Expr.
            Generators of the extension. These should be expressions that are
            algebraic over `\mathbb{Q}`.

        Returns
        =======

        :py:class:`~.AlgebraicField`
            A :py:class:`~.Domain` representing the algebraic field extension.

        Examples
        ========

        >>> from sympy import ZZ, sqrt
        >>> ZZ.algebraic_field(sqrt(2))
        QQ<sqrt(2)>
        """
        return self.get_field().algebraic_field(*extension)

    def from_AlgebraicField(K1, a, K0):
        """Convert a :py:class:`~.ANP` object to :ref:`ZZ`.

        See :py:meth:`~.Domain.convert`.
        """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)

    def log(self, a, b):
        r"""logarithm of *a* to the base *b*

        Parameters
        ==========

        a: number
        b: number

        Returns
        =======

        $\\lfloor\log(a, b)\\rfloor$:
            Floor of the logarithm of *a* to the base *b*

        Examples
        ========

        >>> from sympy import ZZ
        >>> ZZ.log(ZZ(8), ZZ(2))
        3
        >>> ZZ.log(ZZ(9), ZZ(2))
        3

        Notes
        =====

        This function uses ``math.log`` which is based on ``float`` so it will
        fail for large integer arguments.
        """
        return self.dtype(math.log(int(a), b))
