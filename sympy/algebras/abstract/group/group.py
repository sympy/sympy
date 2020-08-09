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

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of one Map
        See sympy.map module.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, Group

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class GroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     left_divisible = right_divisible = True
    ...     associative = True
    ...     identity = e
    >>> op = GroupOp()

    >>> G = Group('G', (S,), (op,))
    >>> G_op = G.operator

    Operation of group is associative.

    >>> G_op(a, op(b, c), evaluate=True)
    a*b*c

    Exponent is defined.

    >>> G_op(a, a, evaluate=True)
    a**2

    Operation of group has identity.

    >>> G_op(a, G_op(e, b), evaluate=True)
    a*b

    Operation of group has inverse element.

    >>> G_op(a, G.inverse(a), evaluate=True)
    e
    >>> G_op(G.inverse(b), b, evaluate=True)
    e

    """

    @property
    def group(self):
        # for compatibility with other structures
        return self

    @property
    def inverse(self):
        op = self.operator
        return op.inverse_operator()

    def negate(self, a, evaluate=False):
        return self.inverse(a, evaluate=evaluate)
    neg = negate

class AbelianGroup(Group):
    """
    A base class for abelian group.

    Explanation
    ===========

    Abelian group is a group whose operation is commutative.

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of one Map
        See sympy.map module.

    Examples
    ========

    >>> from sympy import Set, BinaryOperator, AbelianGroup

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class AbelianGroupOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     left_divisible = right_divisible = True
    ...     associative = True
    ...     commutative = True
    ...     identity = e
    >>> op = AbelianGroupOp()

    >>> G = AbelianGroup('G', (S,), (op,))
    >>> G_op = G.operator

    Operation of abelian group is similar to natural scalar operations.

    >>> G_op(a, b, c, G.inverse(a), b, evaluate=True)
    (b**2)*c

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if not ask(Q.commutative(op)):
            raise TypeError("%s is not commutative." % op)

        return super().__new__(cls, name, sets, operators)
