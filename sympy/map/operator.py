from sympy.assumptions import ask, Q
from sympy.core.sympify import _sympify
from .map import Map, AppliedMap

__all__ = [
    'BinaryOperator', 'AppliedBinaryOperator',
]

class BinaryOperator(Map):
    """
    An abstract class for general algebraic binary operator.

    Explanation
    ===========

    In mathematics, a binary operation is a calculation that combines
    two elements (called operands) to produce another element [1].

    Binary operator can multiple left identities or right identities. However,
    if it has both left and right identity, there is a single two-sided identity
    and no other left or right identity. For this case, you must give ``identity``
    attribute to the operator.

    Examples
    ========

    >>> from sympy.map import BinaryOperator
    >>> from sympy.abc import a, b, e

    >>> class Op(BinaryOperator):
    ...     name = '.'
    ...     identity = e
    >>> op = Op()

    >>> op(a, b)
    a . b

    Binary operation with identity is evaluated.

    >>> op(a, e)
    a . e
    >>> op(a, e, evaluate=True)
    a
    >>> op(e, b, evaluate=True)
    b

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Binary_operation

    """

    # to be overridden
    @property
    def identity(self):
        return

    def flatten(self, seq):
        new_seq = []
        for o in seq:
            if isinstance(o, AppliedMap) and o.map.is_restriction(self):
                new_seq.extend(o.arguments)
            else:
                new_seq.append(o)
        return new_seq

    def left_identity(self, element):
        r"""
        Returns ``True`` if *element* is left identity of *self*.

        Explanation
        ===========

        For a set $S$ and binary operation $*:S \times S \rightarrow S$, $e \in S$ is
        left identity of $*$ if $ e * x = x$ for any $x \in S$.
        There can be multiple left identities. However, if left identity and right identity
        both exist, there is only a single two-sided identity.

        Examples
        ========

        >>> from sympy.map import BinaryOperator
        >>> from sympy.abc import a, b, c

        >>> class Op(BinaryOperator):
        ...     name = '.'
        ...     def _eval_left_identity(self, element):
        ...         return element in (a, b)
        >>> op = Op()

        >>> op(c, a, evaluate=True)
        c . a
        >>> op(a, c, evaluate=True)
        c
        >>> op(b, c, evaluate=True)
        c

        """
        if self.identity is not None:
            if element == self.identity:
                return True
            return False
        result = self._eval_left_identity(element)
        if result is not None:
            return result
        return

    def _eval_left_identity(self, element):
        return

    def right_identity(self, element):
        r"""
        Returns ``True`` *element* is right identity of *self*.

        Explanation
        ===========

        For a set $S$ and binary operation $*:S \times S \rightarrow S$, $e \in S$ is
        left identity of $*$ if $ x * e = x$ for any $x \in S$.
        There can be multiple right identities. However, if left identity and right identity
        both exist, there is only a single two-sided identity.

        Examples
        ========

        >>> from sympy.map import BinaryOperator
        >>> from sympy.abc import a, b, c

        >>> class Op(BinaryOperator):
        ...     name = '.'
        ...     def _eval_right_identity(self, element):
        ...         return element in (a, b)
        >>> op = Op()

        >>> op(a, c, evaluate=True)
        a . c
        >>> op(c, a, evaluate=True)
        c
        >>> op(c, a, evaluate=True)
        c

        """
        if self.identity is not None:
            if element == self.identity:
                return True
            return False
        result = self._eval_right_identity(element)
        if result is not None:
            return result
        return

    def _eval_right_identity(self, element):
        return

    def remove_identity(self, seq):
        # this method is deliberately designed to work for len(seq) > 2 case
        # to be compatible if self is associative
        if self.identity is not None:
            return [o for o in seq if o != self.identity]

        newseq = []
        for i,o in enumerate(seq):
            if self.left_identity(o):
                # a*e != a if e is left identity
                if i != len(seq)-1:
                    continue
            elif self.right_identity(o):
                # e*a != a if e is right identity
                if i != 0:
                    continue
            newseq.append(o)
        return newseq

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedBinaryOperator(self, args, evaluate=evaluate)

class AppliedBinaryOperator(AppliedMap):
    """
    Unevaluated result of BinaryOperator applied to arguments.

    """

    def __new__(cls, mapping, args, evaluate=False, **kwargs):

        associative = ask(Q.associative(mapping))

        if len(args) < 2:
            raise TypeError(
        "%s argument given to binary operator." % len(args)
        )
        if not associative and len(args) > 2:
            raise TypeError(
        "Multiple aruments given to non-associative binary operator."
        )
        args = [_sympify(a) for a in args]

        if evaluate:

            if associative:
                args = mapping.flatten(args)

            args = mapping.remove_identity(args)

            if not args:
                return mapping.identity
            if len(args) == 1:
                return args[0]

        return super().__new__(cls, mapping, args, evaluate=evaluate, **kwargs)
