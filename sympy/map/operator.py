from functools import cmp_to_key
from sympy.assumptions import ask, Q
from sympy.core import Basic, S, Expr, Tuple
from sympy.core.sympify import _sympify
from .map import Map, AppliedMap, IdentityMap, isappliedmap

__all__ = [
    'BinaryOperator', 'LeftDivisionOperator', 'RightDivisionOperator',
    'InverseOperator', 'ExponentOperator',
    'InverseElement', 'ExponentElement',
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

    Binary operation with identity element

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

    def apply(self, *args, **kwargs):
        # This method completely overrides Map.apply

        evaluate = kwargs.get('evaluate', False)
        associative = ask(Q.associative(self))
        commutative = ask(Q.commutative(self))

        # check number of arguments
        if len(args) < 2:
            raise TypeError(
        "%s argument given to binary operator." % len(args)
        )
        if not associative and len(args) > 2:
            raise TypeError(
        "Multiple aruments given to non-associative binary operator."
        )

        # sympify the arguments
        args = [_sympify(a) for a in args]

        if evaluate:
            # flattening, identity removing, etc
            args = self.process_args(args)
            if commutative:
                args.sort(key=cmp_to_key(Basic.compare))

            if not args:
                return self.identity
            elif len(args) == 1:
                return args[0]

            result = self.eval(*args)
            if result is not None:
                return result

        args = Tuple(*args)
        return super(Expr, AppliedMap).__new__(AppliedMap, self, args)

    def flatten(self, seq):
        #Flatten nested structure for associative operator.
        #Other procedures such as inverse cancelling are implemented in other methods

        new_seq = []
        for o in seq:
            if isappliedmap(o, self):
                new_seq.extend(o.arguments)
            else:
                new_seq.append(o)
        return new_seq

    def check_left_identity(self, element):
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
        ...     def _eval_check_left_identity(self, element):
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
        result = self._eval_check_left_identity(element)
        if result is not None:
            return result
        return

    def _eval_check_left_identity(self, element):
        return

    def check_right_identity(self, element):
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
        ...     def _eval_check_right_identity(self, element):
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
        result = self._eval_check_right_identity(element)
        if result is not None:
            return result
        return

    def _eval_check_right_identity(self, element):
        return

    def remove_identity(self, seq):
        # This method is deliberately designed to work for len(seq) > 2 case
        # to be compatible if self is associative
        if self.identity is not None:
            return [o for o in seq if o != self.identity]

        newseq = []
        for i,o in enumerate(seq):
            if self.check_left_identity(o):
                if i != len(seq)-1:
                    # a*e != a if e is left identity
                    continue
            elif self.check_right_identity(o):
                if i != 0:
                    # e*a != a if e is right identity
                    continue
            newseq.append(o)
        return newseq

    def left_division_operator(self):
        """
        Return left division operator with respect to *self*
        """
        return LeftDivisionOperator(self)

    def right_division_operator(self):
        """
        Return right division operator with respect to *self*
        """
        return RightDivisionOperator(self)

    def is_left_division(self, other):
        """
        Return ``True`` if *self* is left division operator of *other*.
        """
        if isinstance(self, LeftDivisionOperator) and self.base_op == other:
            return True
        return False

    def is_right_division(self, other):
        """
        Return ``True`` if *self* is right division operator of *other*.
        """
        if isinstance(self, RightDivisionOperator) and self.base_op == other:
            return True
        return False

    def inverse_operator(self):
        """
        Return unary operator which returns two-sided inverse element
        with repect to *self*.
        """
        return InverseOperator(self)

    def are_inverse(self, a, b):
        """
        Return ``True`` if *a* and *b* are in inverse relation with respect to *self*

        """
        if getattr(self, 'identity', None) is None:
            return False

        a_base, a_exp = a.as_base_exp(self)
        b_base, b_exp = b.as_base_exp(self)
        if a_base == b_base and a_exp == -b_exp:
            return True

    def exponent_operator(self):
        """
        Return binary operator which represents repetitive operation
        of same elements by **self**.

        """
        return ExponentOperator(self)

    def _binary_cancel(self, a, b):
        # attempt to cancel a and b using inverse relation or division relation.
        # used by ``cancel`` method

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
        # method run to evaluate the operation of *self* on *weq*
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

class LeftDivisionOperator(BinaryOperator):
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
    >>> op_ld = op.left_division_operator()

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

        return super().__new__(cls, base_op, **kwargs)

    @property
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain

    @property
    def codomain(self):
        return self.base_op.codomain

    def apply(self, divisor, dividend, **kwargs):
        # if base_op is associative and identity exists,
        # convert the result to operation between dividend and
        # inverse element of divisor.
        base_op = self.base_op
        if base_op.identity is not None and ask(Q.associative(base_op)):
            kwargs["evaluate"] = True
            inv_divisor = base_op.exponent_operator()(divisor, -1, **kwargs)
            return base_op(inv_divisor, dividend, **kwargs)
        return super().apply(divisor, dividend, **kwargs)

class RightDivisionOperator(BinaryOperator):
    r"""
    Right division operator, derived from a binary operation.

    Explanation
    ===========

    For a set $S$ and binary operator $*$ defined on it, a right division $/$
    is an operator that satisfies $(a / b) * b = a$ for any $a \in S$ and $ b \in S$.

    Parameters
    ==========

    base_op : BinaryOperator
        Base operator from where right division is derived

    Examples
    ========

    >>> from sympy import BinaryOperator
    >>> from sympy.abc import a, b
    >>> class Op(BinaryOperator):
    ...     name = '*'
    ...     is_right_divisible = True
    >>> op = Op()
    >>> op_rd = op.right_division_operator()

    >>> op(op_rd(a, b), b)
    (a / b) * b
    >>> op(op_rd(a, b), b, evaluate=True)
    a

    """
    name = '/'

    def __new__(cls, base_op, **kwargs):

        if not ask(Q.right_divisible(base_op)):
            raise TypeError("Right division of %s does not exist." % base_op)

        return super().__new__(cls, base_op, **kwargs)

    @property
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain

    @property
    def codomain(self):
        return self.base_op.codomain

    def apply(self, dividend, divisor, **kwargs):
        # if base_op is associative and identity exists,
        # convert the result to operation between dividend and
        # inverse element of divisor.
        base_op = self.base_op
        if base_op.identity is not None and ask(Q.associative(base_op)):
            kwargs["evaluate"] = True
            inv_divisor = base_op.exponent_operator()(divisor, -1, **kwargs)
            return base_op(dividend, inv_divisor, **kwargs)
        return super().apply(dividend, divisor, **kwargs)

class InverseOperator(Map):
    """
    Unary map which returns two-sided inverse element.

    .. note::
       If the base operator is associative, applying this
       returns ``ExponentElement`` instead of ``InverseElement``.

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

        return super().__new__(cls, base_op, **kwargs)

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
        return InverseElement(self, (x,), evaluate=evaluate, **kwargs)

    def apply(self, x, **kwargs):
        base_op = self.base_op
        if ask(Q.associative(base_op)):
            return base_op.exponent_operator()(x, -1, **kwargs)
        return super().apply(x, **kwargs)

    def eval(self, x):
        if x == self.base_op.identity:
            return x

        base, exp = x.as_base_exp(self.base_op)
        if exp == -1:
            return base

    def _eval_as_base_exp(self, a):
        return a, S.NegativeOne

    def _eval_iterate(self, n):
        if ask(Q.even(n)):
            return IdentityMap(self.domain)
        if ask(Q.odd(n)):
            return self

class ExponentOperator(Map):
    """
    Binary map which represents the repetitive operation with same elements.

    Explanations
    ============

    If an operation is associative, its exponentiation can be defined as
    repetitive application of same elements.
    If identity element exists, zeroth exponentiation is defined and it
    returns the identity element.
    If inverse element exists, negative exponentiation is defined and it
    returns the positive exponentiation of inverse element.

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

    """
    def __new__(cls, base_op, **kwargs):

        if ask(Q.associative(base_op)) is False:
            raise TypeError("%s is not associative." % base_op)

        return super().__new__(cls, base_op, **kwargs)

    @property
    def base_op(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base_op.domain.args[0] * S.Complexes

    @property
    def codomain(self):
        return self.base_op.domain.args[0]

    def is_left_division(self, other):
        return False

    def is_right_division(self, other):
        return False

    def __call__(self, x, n, evaluate=False, **kwargs):
        return ExponentElement(self, (x, n), evaluate=evaluate)

    def apply(self, x, n, **kwargs):
        if not ask(Q.positive(n)):
            if getattr(self.base_op, "identity", None) is None:
                raise TypeError("%s does not have identity." % self.base_op)

        if kwargs.get("evaluate", False):
            x, n = _sympify(x), _sympify(n)
            base, exp = x.as_base_exp(self.base_op)
            if exp != 1:
                # collect x**2**3 to x**6
                x, n = base, exp*n
                return self(x, n, **kwargs)

        return super().apply(x, n, **kwargs)

    def eval(self, x, n):

        if x == self.base_op.identity:
            return x
        if n == 1:
            return x
        if n == 0:
            return self.base_op.identity

    def _eval_as_base_exp(self, x, n):
        return x, n

class InverseElement(AppliedMap):
    """
    Result of InverseOperator applied to element.

    """
    def as_base_exp(self, operator):
        if self.map.base_op.is_restriction(operator):
            return self.map._eval_as_base_exp(*self.arguments)
        return self, S.One

class ExponentElement(AppliedMap):
    """
    Result of ExponentOperator applied to element.

    """

    def as_base_exp(self, operator):
        if self.map.base_op.is_restriction(operator):
            return self.map._eval_as_base_exp(*self.arguments)
        return self, S.One
