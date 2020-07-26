from .structure import AlgebraicStructure

__all__ = [
    "Magma",
]

class Magma(AlgebraicStructure):
    """
    A base class for algebraic magma.

    Explanation
    ===========

    Magma is algebraic structure consists of one set and one binary operation.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Magma

    >>> S = Set('S')
    >>> a, b, c = [S.element(n) for n in 'abc']

    Define operator for magma.

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

    Operation of magma does not need to be associative.

    >>> op(a, op(b, c), evaluate=True)
    a * (b * c)

    """
    def __new__(cls, name, sets, operators, **kwargs):

        if not len(sets) == 1:
            raise TypeError("%s consists of one set." % cls)
        if not (len(operators) == 1 and operators[0].arity == 2):
            raise TypeError("%s consists of one binary operator." % cls)

        obj = super().__new__(cls, name, sets, operators)
        return obj