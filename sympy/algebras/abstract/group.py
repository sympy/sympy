"""Group-like algebraic structures"""

from .structure import AlgebraicStructure

class Magma(AlgebraicStructure):
    """
    Magma is algebraic structure consists of one set and one binary operation.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Magma

    >>> S = Set('S')
    >>> a,b = S.element('a'), S.element('b')

    >>> class MagmaOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    >>> op = MagmaOp()

    >>> M = Magma('M', (S,), (op,))

    >>> a in M and b in M
    True

    >>> op(a, b)
    a * b
    >>> op(a, b) in M
    True

    """
    def __new__(cls, name, sets, operators, **kwargs):

        if not len(sets) == 1:
            raise TypeError("%s consists of one set." % cls)
        if not (len(operators) == 1 and operators[0].arity == 2):
            raise TypeError("%s consists of one binary operator." % cls)

        obj = super().__new__(cls, name, sets, operators)
        return obj

