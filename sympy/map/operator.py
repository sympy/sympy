from functools import cmp_to_key
from sympy.assumptions import ask, Q
from sympy.core import Basic, S
from sympy.core.sympify import _sympify
from .map import Map, AppliedMap, IdentityMap

__all__ = [
    'BinaryOperator', 'LeftDivision', 'RightDivision',
    'InverseOperator', 'ExponentOperator',
    'AppliedBinaryOperator', 'InverseElement', 'ExponentElement',
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
    domain = S.UniversalSet**2
    codomain = S.UniversalSet

    # to be overridden
    @property
    def identity(self):
        return

    def flatten(self, seq):
        #Flatten nested structure for associative operator.
        #Other procedures such as inverse cancelling are implemented in other methods

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
        """
        Return left division operator with respect to *self*
        """
        return LeftDivision(self)

    def right_division(self):
        """
        Return right division operator with respect to *self*
        """
        return RightDivision(self)

    def is_left_division(self, other):
        """
        Return ``True`` if *self* is left division operator of *other*.
        """
        if isinstance(self, LeftDivision) and self.base_op == other:
            return True
        return False

    def is_right_division(self, other):
        """
        Return ``True`` if *self* is right division operator of *other*.
        """
        if isinstance(self, RightDivision) and self.base_op == other:
            return True
        return False

    def inverse_operator(self):
        """
        Return unary operator which returns two-sided inverse element
        with repect to *self*. Not to be confused with ``inverse`` method.
        """
        return InverseOperator(self)

    def are_inverse(self, a, b):
        """
        Return ``True`` if *a* and *b* are in inverse relation with respect to *self*
        """
        if getattr(self, 'identity', None) is None:
            return False
        invop = self.inverse_operator()
        return invop(a) == b or invop(b) == a

    def exponent_operator(self):
        return ExponentOperator(self)

    def _binary_cancel(self, a, b):
        # attempt to cancel a and b using inverse relation or division relation.
        # 1. inverse element
        if self.are_inverse(a, b):
            return []
        # 2. right division
        if isinstance(a, AppliedMap):
            if self.is_right_division(a.map) or a.map.is_right_division(self):
                if a.arguments[-1] == b:
                    return [*a.arguments[:-1]]
        # 3. left division
        if isinstance(b, AppliedMap):
            if self.is_left_division(b.map) or b.map.is_left_division(self):
                if a == b.arguments[0]:
                    return [*b.arguments[1:]]
        return [a, b]
        
    def cancel(self, seq):
        # cancel division and inverse
        # this method is deliberately designed to work for len(seq) > 2 case
        # to be compatible if self is associative
        newseq = []
        for o in seq:
            if not newseq:
                newseq.append(o)
                continue
            a = newseq.pop()
            b = o
            newseq.extend(self._binary_cancel(a, b))
        return newseq

    def collect_iterated(self, seq):
        # convert iterated argument to exponentiation.
        # run only when *self* is associative
        newseq = []
        exp_op = self.exponent_operator()
        b, e = None, None
        for o in seq:
            o_b, o_e = o.as_base_exp(self)
            if b is None:
                b, e = o_b, o_e
            elif b == o_b:
                e += o_e
            else:
                newseq.append(exp_op(b, e, evaluate=True))
                b, e = o_b, o_e
        if not (b is None or e is None):
            newseq.append(exp_op(b, e, evaluate=True))
        return newseq

    def assoc_comm_process(self, seq):
        # Special process for associative and commutative operator
        # Common operators such as addition or multiplication use this method.
        seq = self.remove_identity(seq)

        base_exp_dict = {}
        for o in seq:
            b, e = o.as_base_exp(self)
            if b not in base_exp_dict:
                base_exp_dict[b] = e
            else:
                base_exp_dict[b] += e
        result = []
        exp_op = self.exponent_operator()
        for b,e in base_exp_dict.items():
            if e != 0:
                result.append(exp_op(b, e, evaluate=True))
        return result

    def process_args(self, seq):
        associative = ask(Q.associative(self))
        commutative = ask(Q.commutative(self))

        if associative:
            seq = self.flatten(seq)
        if associative and commutative:
            return self.assoc_comm_process(seq)
        seq = self.remove_identity(seq)
        seq = self.cancel(seq)
        if associative:
            seq = self.collect_iterated(seq)
        return seq

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
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain

    @property
    def codomain(self):
        return self.base_op.codomain

    def eval(self, a, b):
        base_op = self.base_op
        if base_op.identity is not None and ask(Q.associative(base_op)):
            if b == base_op.identity:
                return InverseOperator(base_op)(a)

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
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain

    @property
    def codomain(self):
        return self.base_op.codomain

    def eval(self, a, b):
        base_op = self.base_op
        if base_op.identity is not None and ask(Q.associative(base_op)):
            if a == base_op.identity:
                return InverseOperator(base_op)(b)

class InverseOperator(Map):
    """
    Unary map which returns two-sided inverse element.

    .. note::
       Not to be confused with ``InverseMap``.

    Parameters
    ==========

    base_op : BinaryOperator
        Base operator from where inverse operator is derived

    Examples
    ========

    >>> from sympy import BinaryOperator, InverseOperator, S
    >>> from sympy.abc import x
    >>> class Op(BinaryOperator):
    ...     identity = S.One
    >>> op = Op()
    >>> op_invop = InverseOperator(op)

    >>> op(x, op_invop(x), evaluate=True)
    1
    >>> op(op_invop(x), x, evaluate=True)
    1

    """
    def __new__(cls, base_op, **kwargs):

        if base_op.identity is None:
            raise TypeError("%s does not have identity." % base_op)

        return super().__new__(cls, base_op)

    @property
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain.args[0]

    @property
    def codomain(self):
        return self.base_op.codomain

    def is_left_division(self, other):
        return False

    def is_right_division(self, other):
        return False

    def __call__(self, x, evaluate=False, **kwargs):
        return InverseElement(self, (x,), evaluate=evaluate)

    def eval(self, x):
        if self.base_op.identity is not None and x == self.base_op.identity:
            return x

    def _eval_as_base_exp(self, a):
        return a, S.NegativeOne

    def _eval_iterate(self, n):
        if ask(Q.even(n)):
            return IdentityMap(self.domain)
        if ask(Q.odd(n)):
            return self

class ExponentOperator(Map):
    """
    Binary map which returns exponentiation.

    Parameters
    ==========

    base_op : BinaryOperator
        Base operator from where exponent operator is derived

    Examples
    ========

    >>> from sympy import BinaryOperator, ExponentOperator, S
    >>> from sympy.abc import x
    >>> class Op(BinaryOperator):
    ...     identity = S.One
    ...     is_associative = True
    >>> op = Op()
    >>> op_expop = ExponentOperator(op)

    >>> op(x, x, evaluate=True) == op_expop(x, 2)
    True

    >>> op_expop(x, 1, evaluate=True)
    x
    >>> op_expop(x, 0, evaluate=True)
    1
    >>> op_expop(x, -1, evaluate=True) == op.inverse_operator()(x)
    True

    """
    def __new__(cls, base_op, **kwargs):

        if ask(Q.associative(base_op)) is False:
            raise TypeError("%s is not associative." % base_op)

        return super().__new__(cls, base_op)

    @property
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain.args[0] * S.Integers

    @property
    def codomain(self):
        return self.base_op.domain.args[0]

    def is_left_division(self, other):
        return False

    def is_right_division(self, other):
        return False

    def __call__(self, x, n, evaluate=False, **kwargs):
        return ExponentElement(self, (x, n), evaluate=evaluate)

    def eval(self, x, n):

        if n <= 0:
            if getattr(self.base_op, "identity", None) is None:
                raise TypeError(
            "Negative exponent cannot be defined since %s does not have identity." % self
            )

        if n == 1:
            return x
        if n == 0:
            return self.base_op.identity
        if n == -1:
            return self.base_op.inverse_operator()(x)
        if n < -1:
            return self(self.base_op.inverse_operator()(x), -n)

    def _eval_as_base_exp(self, x, n):
        b, e = x.as_base_exp(self.base_op)
        return b, n*e

class AppliedBinaryOperator(AppliedMap):
    """
    Unevaluated result of BinaryOperator applied to arguments.

    """

    def __new__(cls, mapping, args, evaluate=False, **kwargs):

        associative = ask(Q.associative(mapping))
        commutative = ask(Q.commutative(mapping))

    # check argument
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

            args = mapping.process_args(args)

            if commutative:
                args.sort(key=cmp_to_key(Basic.compare))

            if not args:
                result = mapping.identity
            elif len(args) == 1:
                result = args[0]
            else:
                result = super().__new__(cls, mapping, args, evaluate=True, **kwargs)

            # check result
            if mapping.codomain.contains(result) == False:
                raise TypeError(
            "%s is not in %s's codomain %s." % (result, mapping, mapping.codomain)
            )

            return result

        return super().__new__(cls, mapping, args, evaluate=False, **kwargs)

class InverseElement(AppliedMap):
    """
    Result of InverseOperator applied to element.

    """
    def __new__(cls, mapping, args, evaluate=False, **kwargs):
        x = _sympify(args[0])
        if evaluate:
            base, exp = x.as_base_exp(mapping.base_op)
            if exp == -1:
                return base
            elif exp != 1:
                exp_op = mapping.base_op.exponent_operator()
                inv_elem = mapping(base, evaluate=True)
                return exp_op(inv_elem, exp, evaluate=True)
        return super().__new__(cls, mapping, (x,), evaluate=evaluate, **kwargs)

    def as_base_exp(self, operator):
        if self.map.base_op.is_restriction(operator):
            return self.map._eval_as_base_exp(*self.arguments)
        return self, S.One

class ExponentElement(AppliedMap):
    """
    Result of ExponentOperator applied to element.

    """
    def __new__(cls, mapping, args, evaluate=False, **kwargs):
        x, n = map(_sympify, args)
        if evaluate:
            base, exp = x.as_base_exp(mapping.base_op)
            x, n = base, exp*n
        return super().__new__(cls, mapping, (x, n), evaluate=evaluate, **kwargs)

    def as_base_exp(self, operator):
        if self.map.base_op.is_restriction(operator):
            return self.map._eval_as_base_exp(*self.arguments)
        return self, S.One
