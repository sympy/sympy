from .semigroup import Semigroup

__all__ = [
    "Monoid",
]

class Monoid(Semigroup):
    """
    A base class for monoid.

    Explanation
    ===========

    Monoid is semigroup whose operation has two-sided identity.

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

    >>> from sympy import Set, BinaryOperator, Monoid

    >>> S = Set('S')
    >>> a, b, c, e = [S.element(n) for n in 'abce']

    >>> class MonoidOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     associative = True
    ...     identity = e
    >>> op = MonoidOp()

    >>> M = Monoid('M', (S,), (op,))
    >>> M_op = M.operator

    Operation of monoid is associative.

    >>> M_op(a, op(b, c), evaluate=True)
    a*b*c

    Operation of monoid has identity.

    >>> M_op(a, op(e, b), evaluate=True)
    a*b

    """
    def __new__(cls, name, sets, operators, **kwargs):
        op = operators[0]
        if getattr(op, 'identity', None) is None:
            raise TypeError("%s does not have identity." % op)

        return super().__new__(cls, name, sets, operators)

    @property
    def identity(self):
        op = self.operators[0]
        return op.identity
