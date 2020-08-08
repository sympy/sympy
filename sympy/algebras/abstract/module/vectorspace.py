from ..ring import Field
from .module import Module

__all__ = [
    "VectorSpace",
]

class VectorSpace(Module):
    """
    A base class for vector space.

    Explanation
    ===========

    Vector space is a module whose ring is field.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of two structures
        The first one is scalar field, and the
        second one is abelian group.

    operators : tuple of one Map
        This map is scalar multiplication operator.

    Examples
    ========

    >>> from sympy import (
    ... Set, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, S, VectorSpace
    ... )

    >>> F = S.RealsField

    >>> X = Set('X')
    >>> x, e = [X.element(i) for i in 'xe']
    >>> vadd = VectorAdditionOperator(X**2, X, e)
    >>> G = AbelianGroup('G', (X,), (vadd,))

    >>> smul = ScalarMultiplicationOperator(F*G, G)

    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> V.div(x, 2, evaluate=True)
    (2**(-1))*x
    >>> V.div(x, 1, evaluate=True)
    x

    """
    def __new__(cls, name, sets, operators, **kwargs):
        scalar_field, abelian_group = sets

        if not isinstance(scalar_field, Field):
            raise TypeError("%s is not field." % scalar_field)

        return super().__new__(cls, name, sets, operators)

    @property
    def field(self):
        return self.args[1].args[0]

    def divide(self, a, b, **kwargs):

        if not self.check_scalar(b) == True:
            raise TypeError("Divisor must be scalar.")

        inv_b = self.field.mul_op.inverse_element(b, evaluate=True)
        return self.mul(a, inv_b, evaluate=True)
    div = divide
