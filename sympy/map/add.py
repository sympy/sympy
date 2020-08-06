"""Common addition operator"""

from functools import cmp_to_key
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap
from .operator import BinaryOperator

__all__ = [
    "AdditionOperator", "NumericAdditionOperator", "Addition"
]

class AdditionOperator(BinaryOperator):
    """
    Class for abstract addition operator defined in algebraic ring.

    Explanation
    ===========

    Addition operater forms an abelian group over ring structure, and can be subject to
    distribution by multiplication.
    When ``AdditionOperator`` instance is called, multiplication operator must be passed
    as *mul_op* parameter for this relation. By using ``Ring`` class from algebras module,
    this is done automatically.

    .. note::
        Since generalized addition does not have to satisfy the common addition
        of numbers, repetitive addition or inverse element of addition does not
        return the result of multiplication.

    Examples
    ========

    >>> from sympy import AdditionOperator, Set, BinaryOperator

    >>> A = Set('A')
    >>> a, b, e1, e2 = [A.element(n) for n in ('a', 'b', 'e1', 'e2')]

    >>> add = AdditionOperator(A**2, A, e1)
    >>> class MonoidOp(BinaryOperator):
    ...     name = '*'
    ...     domain = A*A
    ...     codomain = A
    ...     is_associative = True
    ...     identity = e2
    >>> mul = MonoidOp()

    >>> add(a, b, mul_op=mul)
    a + b
    >>> add(a, a, mul_op=mul, evaluate=True)
    (e2 + e2)*a

    Constructing a ring makes applying the operator easier

    >>> from sympy import Ring
    >>> R = Ring('R', (A,), (add, mul))

    >>> R.add(a, a, evaluate=True)
    (e2 + e2)*a

    See Also
    ========

    algebras.Ring

    """
    name = '+'
    is_associative = True
    is_commutative = True
    is_left_divisible = is_right_divisible = True

    def __new__(cls, domain, codomain, identity, **kwargs):
        domain, codomain, identity = _sympify(domain), _sympify(codomain), _sympify(identity)
        return super().__new__(cls, domain, codomain, identity, **kwargs)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    @property
    def identity(self):
        return self.args[2]

    def subtraction_operator(self):
        return self.right_division_operator()

    def __call__(self, *args, mul_op, evaluate=False, **kwargs):
        return Addition(self, args, mul_op, evaluate=evaluate, **kwargs)

    def apply(self, *args, **kwargs):
        mul_op = kwargs.get('mul_op')
        evaluate = kwargs.get('evaluate', False)

        # sympify the arguments
        args = [_sympify(a) for a in args]

        if evaluate:
            args = self.addition_process(args, mul_op)

            if not args:
                return self.identity
            elif len(args) == 1:
                return args[0]

        # return Addition class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])

        result = super(AppliedMap, Addition).__new__(Addition, self, args, mul_op)
        return result

    def undistribute(self, seq, mul_op, evaluate=False):
        # collect coefficients, e.g. 2*x + 3*x -> 5*x
        # also, remove the identity element from seq
        terms = {}
        for o in seq:
            coeff, term = o.as_coeff_Mul(mul_op=mul_op)
            if term not in terms:
                terms[term] = mul_op.identity
            else:
                terms[term] = self(terms[term], coeff, mul_op=mul_op, evaluate=False)

        result = []
        for term, coeff in terms.items():
            if coeff == self.identity:
                continue
            elif coeff == mul_op.identity:
                result.append(term)
            else:
                result.append(mul_op(coeff, term, add_op=self, evaluate=evaluate))
        return result

    def cancel(self, seq, evaluate=False):
        terms = {}
        for o in seq:
            base, exp = o.as_base_exp(self)
            if base not in terms:
                terms[base] = exp
            else:
                terms[base] += exp

        result = []
        expop = self.exponent_operator()
        for base, exp in terms.items():
            result.append(expop(base, exp, evaluate=evaluate))
        return result

    def addition_process(self, seq, mul_op):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        seq = self.undistribute(seq, mul_op, evaluate=True)
        seq = self.cancel(seq, evaluate=True)
        seq.sort(key=cmp_to_key(Basic.compare))
        return seq

    def _ldiv_apply(self, ldivop, divisor, dividend, **kwargs):
        inv_divisor = self.exponent_operator()(divisor, -1, **kwargs)
        return self(inv_divisor, dividend, **kwargs)

    def _rdiv_apply(self, rdivop, dividend, divisor, **kwargs):
        inv_divisor = self.exponent_operator()(divisor, -1, **kwargs)
        return self(dividend, inv_divisor, **kwargs)

class NumericAdditionOperator(AdditionOperator):
    """
    Class for common numeric addition operator.

    Explanation
    ===========

    Unlike ``AdditionOperator`` which is purely abstract operation, this class is
    common addition operator with its identity element fixed to 0. Also, multiplication
    operator does not need to be passed when applying the arguments.

    Examples
    ========

    >>> from sympy import NumericAdditionOperator, S, Symbol
    >>> a = Symbol('a', real=True)

    >>> add = NumericAdditionOperator(S.Reals**2, S.Reals)

    >>> add(1, 2, evaluate=True)
    3
    >>> add(a, 0, evaluate=True)
    a
    >>> add(a, a, evaluate=True)
    2*a

    """
    def __new__(cls, domain, codomain, **kwargs):
        domain, codomain = _sympify(domain), _sympify(codomain)
        return super(BinaryOperator, cls).__new__(cls, domain, codomain, **kwargs)

    @property
    def identity(self):
        return S.Zero

    def multiplication_operator(self):
        from .mul import NumericMultiplicationOperator
        return NumericMultiplicationOperator(self.domain, self.codomain)
    mul_op = multiplication_operator

    def __call__(self, *args, evaluate=False, **kwargs):
        kwargs.pop('mul_op', None)
        mul_op = self.mul_op()
        return Addition(self, args, mul_op, evaluate=evaluate, **kwargs)

    def undistribute(self, seq, mul_op, evaluate=False):
        terms = {}
        for o in seq:
            coeff, term = o.as_coeff_Mul(mul_op=mul_op)
            if term not in terms:
                terms[term] = mul_op.identity
            else:
                terms[term] += coeff

        result = []
        for term, coeff in terms.items():
            if coeff == self.identity:
                continue
            elif coeff == mul_op.identity:
                result.append(term)
            else:
                result.append(mul_op(coeff, term, evaluate=evaluate))
        return result

    def gather_num(self, seq):
        num = self.identity
        newseq = []
        for o in seq:
            if o.is_Number:
                num += o
            else:
                newseq.append(o)
        if num != self.identity:
            newseq.insert(0, num)
        return newseq

    def addition_process(self, seq, mul_op):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        seq = self.undistribute(seq, mul_op, evaluate=True)
        seq = self.gather_num(seq)
        seq.sort(key=cmp_to_key(Basic.compare))
        return seq

    def _expop_apply(self, x, n, **kwargs):
        kwargs['evaluate'] = True
        mul_op = self.mul_op()
        return mul_op(n, x, **kwargs)

class Addition(AppliedMap):
    """
    Class for the unevaluated result of general addition.

    """
    def __new__(cls, mapping, args, mul_op, evaluate=False, **kwargs):
        kwargs.update(
            evaluate=evaluate,
            mul_op=mul_op
        )

        return super().__new__(cls, mapping, args, **kwargs)

    @property
    def mul_op(self):
        return self.args[2]

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args, self.mul_op)

    def as_coeff_Add(self, rational=False, add_op=None, **kwargs):
        if add_op is None or isappliedmap(self, add_op):
            coeff, args = self.arguments[0], self.arguments[1:]
            if coeff.is_Number and not rational or coeff.is_Rational:
                if len(args) == 1:
                    return coeff, args[0]
                else:
                    return coeff, self._new_rawargs(*args)
        return self.map.identity, self

    def undistribute(self, mul_op, evaluate=False):
        # collect coefficients, e.g. 2*x + 3*x -> 5*x
        seq = [*self.arguments]
        return self.map.undistribute(seq, mul_op=mul_op, evaluate=evaluate)
