from sympy.assumptions import ask, Q
from .monoid import Monoid
from .loop import Loop

__all__ = [
    "Group", "AbelianGroup",
]

class Group(Monoid, Loop):
    """
    A base class for algebraic group.

    Explanation
    ===========

    Group is both monoid and loop. It consists of a set and an
    operator which is closed, associative, has a two-sided identity
    element, and inverse element for every element.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Group

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class GroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = is_right_divisible = True
    ...     is_associative = True
    ...     identity = e
    >>> op = GroupOp()

    >>> G = Group('G', (S,), (op,))

    Operation of group is associative.

    >>> op(a, op(b, c), evaluate=True)
    a * b * c

    Exponent is defined.

    >>> op(a, a, evaluate=True)
    a**2

    Operation of group has identity.

    >>> op(a, op(e, b), evaluate=True)
    a * b

    Operation of group has inverse element.

    >>> op(a, G.inverse(a), evaluate=True)
    e
    >>> op(G.inverse(b), b, evaluate=True)
    e

    """
    @property
    def inverse_operator(self):
        op = self.operators[0]
        return op.inverse_operator()

    def inverse(self, element):
        return self.inverse_operator(element, evaluate=True)

class AbelianGroup(Group):
    """
    A base class for abelian group.

    Explanation
    ===========

    Abelian group is a group whose operation is commutative.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, AbelianGroup

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class AbelianGroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = is_right_divisible = True
    ...     is_associative = True
    ...     is_commutative = True
    ...     identity = e
    >>> op = AbelianGroupOp()

    >>> G = AbelianGroup('G', (S,), (op,))

    Operation of abelian group is similar to natural scalar operations.

    >>> op(a, b, c, G.inverse(a), b, evaluate=True)
    b**2 * c

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.commutative(op)):
            raise TypeError("%s is not commutative." % op)

        return super().__new__(cls, name, sets, operators)
