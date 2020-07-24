from sympy.assumptions import ask, Q
from sympy.core.sympify import _sympify
from .map import Map, AppliedMap

__all__ = [
    'BinaryOperator', 'LeftDivision', 'RightDivision',
    'AppliedBinaryOperator',
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
        Returns ``True`` if *element* is right identity of *self*.

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

    def left_division(self):
        return LeftDivision(self)

    def right_division(self):
        return RightDivision(self)

    def is_left_division(self, other):
        """
        Return ``True`` if *self* is left division operator of *other*.
        """
        if isinstance(self, LeftDivision) and self.base == other:
            return True
        return False

    def is_right_division(self, other):
        """
        Return ``True`` if *self* is right division operator of *other*.
        """
        if isinstance(self, RightDivision) and self.base == other:
            return True
        return False

    def cancel_division(self, seq):
        # this method is deliberately designed to work for len(seq) > 2 case
        # to be compatible if self is associative
        newseq = []
        skip = False
        for i,o in enumerate(seq):
            if skip:
                skip = False
                continue
            if isinstance(o, AppliedBinaryOperator):
                if o.map.is_left_division(self) or self.is_left_division(o.map):
                    if i>0 and newseq[-1] == o.arguments[0]:
                        newseq.pop()
                        newseq.extend(o.arguments[1:])
                    else:
                        newseq.append(o)
                    continue
                if o.map.is_right_division(self) or self.is_right_division(o.map):
                    if i<len(seq)-1 and o.arguments[-1] == seq[i+1]:
                        newseq.extend(o.arguments[:-1])
                        skip = True
                        continue
                    else:
                        newseq.append(o)
                    continue
            newseq.append(o)
        return newseq

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedBinaryOperator(self, args, evaluate=evaluate)

class LeftDivision(BinaryOperator):
    r"""
    Left division operator, derived from a binary operation.

    Explanation
    ===========

    For a set $S$ and binary operator $*$ defined on it, a left division $\backslash$
    is an operator that satisfies $a * (a \backslash b) = b$ for any $a \in S$ and $ b \in S$.

    Parameters
    ==========

    base_op : BinaryOperator
        Base operator from where left division is derived

    Examples
    ========

    >>> from sympy import BinaryOperator
    >>> from sympy.abc import a, b
    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     is_left_divisible = True
    >>> op = Op()
    >>> op_ld = op.left_division()

    >>> op(b, op_ld(b, a))
    b * (b \ a)
    >>> op(b, op_ld(b, a), evaluate=True)
    a

    """

    latex_name = r'\backslash'
    str_name = '\\'

    def __new__(cls, base_op, **kwargs):

        if not ask(Q.left_divisible(base_op)):
            raise TypeError("Left division of %s does not exist." % base_op)

        return super().__new__(cls, base_op)

    @property
    def base(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

class RightDivision(BinaryOperator):
    r"""
    Right division operator, derived from a binary operation.

    Explanation
    ===========

    For a set $S$ and binary operator $*$ defined on it, a right division $/$
    is an operator that satisfies $(a / b) * b = a$ for any $a \in S$ and $ b \in S$.

    Parameters
    ==========

    base_op : BinaryOperator
        Base operator from where left division is derived

    Examples
    ========

    >>> from sympy import BinaryOperator
    >>> from sympy.abc import a, b
    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     is_right_divisible = True
    >>> op = Op()
    >>> op_rd = op.right_division()

    >>> op(op_rd(a, b), b)
    (a / b) * b
    >>> op(op_rd(a, b), b, evaluate=True)
    a

    """
    name = '/'

    def __new__(cls, base_op, **kwargs):

        if not ask(Q.right_divisible(base_op)):
            raise TypeError("Right division of %s does not exist." % base_op)

        return super().__new__(cls, base_op)

    @property
    def base(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base.domain

    @property
    def codomain(self):
        return self.base.codomain

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

        if associative:
            s = mapping.domain.args[0]
            for a in args:
                if s.contains(a) == False:
                    raise TypeError(
                "%s is not in %s's set %s." % (a, mapping, s)
                )
        elif not associative and mapping.domain.contains(tuple(args)) == False:
            raise TypeError(
        "%s is not in %s's domain %s." % (tuple(args), mapping, mapping.domain)
        )

        if evaluate:

            if associative:
                args = mapping.flatten(args)

            args = mapping.remove_identity(args)
            args = mapping.cancel_division(args)

            if not args:
                result = mapping.identity
            elif len(args) == 1:
                result = args[0]
            else:
                result = super().__new__(cls, mapping, args, evaluate=True, **kwargs)

            if mapping.codomain.contains(result) == False:
                raise TypeError(
            "%s is not in %s's codomain %s." % (result, mapping, mapping.codomain)
            )

            return result

        return super().__new__(cls, mapping, args, evaluate=False, **kwargs)
