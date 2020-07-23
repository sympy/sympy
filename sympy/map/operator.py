from sympy.core.sympify import _sympify
from .map import Map, AppliedMap

__all__ = [
    'BinaryOperator', 'AssociativeOperator', 'AppliedBinaryOperator', 'AppliedAssociativeOperator'
]

class BinaryOperator(Map):
    """
    An abstract class for general binary operator.

    Explanation
    ===========

    In mathematics, a binary operation is a calculation that combines
    two elements (called operands) to produce another element [1].
    Unlike ordinary mappings, binary operators are often written using
    infix notation [1].

    Examples
    ========

    >>> from sympy.map import BinaryOperator
    >>> from sympy.abc import x, y

    >>> class AddOp(BinaryOperator):
    ...     name = '+'
    >>> addop = AddOp()

    >>> addop(x, y)
    x + y

    BinaryOperator does not flatten nested arguments.

    >>> addop(addop(x, y), y)
    (x + y) + y

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Binary_operation

    """

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedBinaryOperator(self, args, evaluate=evaluate)

class AssociativeOperator(BinaryOperator):
    r"""
    An abstract class for associative binary operator

    Explanation
    ===========

    Associative operator $*$ satisfies that $\left( f * g \right) * h$
    equals $f * \left( g * h \right)$.

    .. note::
       Although associative operators are mathematically binary, nested
       application of them can be designed to be n-ary in SymPy for
       performance reason.

    Examples
    ========

    >>> from sympy.map import AssociativeOperator
    >>> from sympy.abc import x, y

    >>> class AddOp(AssociativeOperator):
    ...     name = '+'
    >>> addop = AddOp()

    AssociativeOperator flattens nested arguments.

    >>> addop(addop(x, y), y, evaluate=True)
    x + y + y

    """

    is_associative = True

    def flatten(self, seq):
        new_seq = []
        while seq:
            o = seq.pop()
            if isinstance(o, AppliedMap) and o.map.is_restriction(self):
                seq.extend(o.arguments)
            else:
                new_seq.append(o)
        new_seq.reverse()
        return new_seq

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedAssociativeOperator(self, args, evaluate=evaluate)

class AppliedBinaryOperator(AppliedMap):
    r"""
    Unevaluated result of BinaryOperator applied to arguments.

    Examples
    ========

    >>> from sympy.map import BinaryOperator
    >>> from sympy.abc import x, y

    >>> class AddOp(BinaryOperator):
    ...     name = '+'
    >>> addop = AddOp()
    >>> expr = addop(addop(x, y), y, evaluate=True)

    >>> expr.arguments # Not flattened
    (x + y, y)

    """

class AppliedAssociativeOperator(AppliedBinaryOperator):
    """
    Unevaluated result of AssociativeOperator applied to arguments.

    Examples
    ========

    >>> from sympy.map import AssociativeOperator
    >>> from sympy.abc import x, y

    >>> class AddOp(AssociativeOperator):
    ...     name = '+'
    >>> addop = AddOp()
    >>> expr = addop(addop(x, y), y, evaluate=True)

    >>> expr.arguments  # Flattened
    (x, y, y)

    """
    def __new__(cls, mapping, args, evaluate=False, **kwargs):
        args = [_sympify(a) for a in args]

        if evaluate:
            args = mapping.flatten(args)

        return super().__new__(cls, mapping, args, evaluate=evaluate, **kwargs)
