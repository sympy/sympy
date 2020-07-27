from .quasigroup import Quasigroup

__all__ = [
    "Loop"
]

class Loop(Quasigroup):
    """
    A base class for algebraic loop.

    Explanation
    ===========

    Loop is quasigroup whose operation has two-sided identity.

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

    >>> from sympy import Set, BinaryOperator, Loop

    >>> S = Set('S')
    >>> a, b, e = [S.element(n) for n in 'abe']

    >>> class LoopOp(BinaryOperator):
    ...     name = '*'
    ...     domain = S*S
    ...     codomain = S
    ...     is_left_divisible = is_right_divisible = True
    ...     identity = e
    >>> op = LoopOp()

    >>> L = Loop('L', (S,), (op,))

    Left division and right division exist.

    >>> op_ld, op_rd = op.left_division(), op.right_division()
    >>> op_ld(a, b)
    a \ b
    >>> op_rd(a, b)
    a / b

    Identity exists.

    >>> op(a, e, evaluate=True)
    a
    >>> op(e, b, evaluate=True)
    b

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
