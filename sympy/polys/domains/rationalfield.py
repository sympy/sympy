"""Implementation of :class:`RationalField` class. """


from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.utilities import public

@public
class RationalField(Field, CharacteristicZero, SimpleDomain):
    """The domain ``QQ`` representing the integers `\mathbb{Q}`.

    The :py:class:`RationalField` class represents the ring of integers as a
    :py:class:`~.Domain` in the domain system. :py:class:`RationalField` is a
    super class of :py:class:`PythonRationalField` and
    :py:class:`GMPYRationalField` one of which will be the implementation for
    :ref:`QQ` depending on whether or not ``gmpy`` or ``gmpy2`` is installed.

    See also
    ========

    Domain
    """

    rep = 'QQ'

    is_RationalField = is_QQ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`.

        Parameters
        ==========

        *extension: One or more Expr
            Generators of the extension. These should be expressions that are
            algebraic over `\mathbb{Q}`.

        Returns
        =======

        :py:class:`~.AlgebraicField`
            A :py:class:`~.Domain` representing the algebraic field extension.

        Examples
        ========

        >>> from sympy import QQ, sqrt
        >>> QQ.algebraic_field(sqrt(2))
        QQ<sqrt(2)>
        """
        from sympy.polys.domains import AlgebraicField
        return AlgebraicField(self, *extension)

    def from_AlgebraicField(K1, a, K0):
        """Convert a :py:class:`~.ANP` object to :ref:`QQ`.

        See :py:meth:`~.Domain.convert`
        """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)
