"""Common addition operator"""

from functools import cmp_to_key
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap
from .operator import (
    BinaryOperator, AppliedBinaryOperator,
)

__all__ = [
    "AdditionOperator", "Addition",
    "scalar_add",
]

class AdditionOperator(BinaryOperator):
    """
    Class for general addition operator.

    Explanation
    ===========

    When called, *mul_op* can be passed to transform the repetitive arguments.

    Examples
    ========

    >>> from sympy import scalar_add, scalar_mul
    >>> from sympy.abc import x, y

    >>> scalar_add(x, y)
    x + y

    Multiplication operator not given:

    >>> scalar_add(x, x, evaluate=True)
    x + x

    Multiplication operator given:

    >>> scalar_add(x, x, mul_op=scalar_mul, evaluate=True)
    2 * x

    """
    name = '+'
    is_associative = True
    is_commutative = True
    is_left_divisible = is_right_divisible = True

    def __new__(cls, domain, codomain, identity, **kwargs):
        domain, codomain, identity = _sympify(domain), _sympify(codomain), _sympify(identity)
        return super().__new__(cls, domain, codomain, identity)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    @property
    def identity(self):
        return self.args[2]

    def __call__(self, *args, mul_op=None, evaluate=False):
        return Addition(self, args, mul_op, evaluate=evaluate)

    def undestribute(self, seq, mul_op):
        # collect coefficients, e.g. 2*x + 3*x -> 5*x
        terms = {}
        for o in seq:
            coeff, term = o.as_coeff_Mul(mul_op=mul_op)
            if term not in terms:
                terms[term] = mul_op.identity
            else:
                terms[term] += coeff

        result = []
        for term, coeff in terms.items():
            if coeff == S.Zero:
                continue
            elif coeff == mul_op.identity:
                result.append(term)
            else:
                result.append(mul_op(coeff, term, evaluate=True))
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
        if mul_op is not None:
            result = self.undestribute(seq, mul_op)
        else:
            result = self.cancel(seq)
        result = self.gather_num(result)
        return result

class Addition(AppliedBinaryOperator):
    """
    Class for the unevaluated result of general addition.

    """
    def __new__(cls, mapping, args, mul_op=None, evaluate=False, **kwargs):
        args = [_sympify(a) for a in args]

        if evaluate:
            args = mapping.addition_process(args, mul_op)

        if not args:
            return mapping.identity
        elif len(args) == 1:
            return args[0]

        args.sort(key=cmp_to_key(Basic.compare))
        args = Tuple(*args)

        if mul_op is None:
            return super(AppliedMap, cls).__new__(cls, mapping, args)

        return super(AppliedMap, cls).__new__(cls, mapping, args, mul_op)

    @property
    def mul_op(self):
        if len(self.args) < 3:
            return None
        return self.args[2]

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args, self.mul_op)

scalar_add = AdditionOperator(S.Complexes**2, S.Complexes, S.Zero)
