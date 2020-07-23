"""Group-like algebraic structures"""

from .structure import AlgebraicStructure

__all__ = [
    "Magma", "Semigroup", "Monoid",
]

class Magma(AlgebraicStructure):
    """
    Magma is algebraic structure consists of one set and one binary operation.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Magma

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

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

class Semigroup(Magma):
    """
    Semigroup is magma whose operation is associative.

    Examples
    ========

    >>> from sympy import Set, AssociativeOperator, Semigroup

    >>> S = Set('S')
    >>> a, b, c = S.element('a'), S.element('b'), S.element('c')

    >>> class SemigroupOp(AssociativeOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    >>> op = SemigroupOp()

    >>> G = Semigroup('G', (S,), (op,))

    >>> op(a, b)
    a * b
    >>> op(a, b) in G
    True

    Operation of semigroup is associative.

    >>> op(a, op(b, c), evaluate=True)
    a * b * c

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not op.is_associative:
            raise TypeError("%s is not associative." % op)

        return super().__new__(cls, name, sets, operators)

class Monoid(Semigroup):
    """
    Monoid is semigroup whose operation has two-sided identity.

    Examples
    ========

    >>> from sympy import Set, AssociativeOperator, Monoid

    >>> S = Set('S')
    >>> a, b, e = S.element('a'), S.element('b'), S.element('e')

    >>> class MonoidOp(AssociativeOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     identity = e
    >>> op = MonoidOp()

    >>> M = Monoid('M', (S,), (op,))

    >>> op(a, b)
    a * b
    >>> op(a, b) in M
    True

    Operation of monoid is associative.

    >>> op(a, op(b, b), evaluate=True)
    a * b * b

    Operation of monoid has identity.

    >>> op(a, op(e, b), evaluate=True)
    a * b

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if getattr(op, 'identity', None) is None:
            raise TypeError("%s does not have identity." % op)

        return super().__new__(cls, name, sets, operators)
