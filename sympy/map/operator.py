from .map import Map, AppliedMap

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

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Binary_operation

    """

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

    """

    def flatten(self, seq):
        new_seq = []
        while seq:
            o = seq.pop()
            if isinstance(o, AppliedMap) and o.map == self:
                seq.extend(o.arguments)
            else:
                new_seq.append(o)
        new_seq.reverse()
        return new_seq
