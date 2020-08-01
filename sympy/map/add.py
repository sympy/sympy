"""Common addition operator"""

from functools import cmp_to_key
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap
from .operator import BinaryOperator

__all__ = [
    "AdditionOperator", "Addition",
    "scalar_add",
]

class AdditionOperator(BinaryOperator):
    """
    Class for general addition operator.

    Explanation
    ===========

    When defined as group operator, addition operator has no special feature. However,
    when defined as ring operator, addition is paired with multiplication and subject to
    distribution by multiplication.
    When ``AdditionOperator`` is called, multiplication operator can be passed as *mul_op*
    parameter for this relation. This implies that two operators form a ring over their
    domain.
    If you use ``Ring`` class from algebras module, two operators are automatically related.

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

    def __call__(self, *args, mul_op=None, evaluate=False):
        return Addition(self, args, mul_op, evaluate=evaluate)

    def apply(self, *args, **kwargs):
        evaluate = kwargs.get('evaluate', False)
        mul_op = kwargs.get('mul_op', None)

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

        if mul_op is None:
            result = super(AppliedMap, Addition).__new__(Addition, self, args)
        else:
            result = super(AppliedMap, Addition).__new__(Addition, self, args, mul_op)
        return result

    def undistribute(self, seq, mul_op, evaluate=False):
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
        if mul_op is not None:
            result = self.undistribute(seq, mul_op, evaluate=True)
        else:
            result = self.cancel(seq)
        result = self.gather_num(result)
        result.sort(key=cmp_to_key(Basic.compare))
        return result

class Addition(AppliedMap):
    """
    Class for the unevaluated result of general addition.

    """
    def __new__(cls, mapping, args, mul_op=None, evaluate=False, **kwargs):
        kwargs.update(
            evaluate=evaluate,
            mul_op=mul_op
        )

        return super().__new__(cls, mapping, args, **kwargs)

    @property
    def mul_op(self):
        if len(self.args) < 3:
            return None
        return self.args[2]

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args, self.mul_op)

    def undistribute(self, mul_op=None, evaluate=False):
        # collect coefficients, e.g. 2*x + 3*x -> 5*x
        if mul_op is None:
            mul_op = self.mul_op

        if mul_op is not None:
            seq = [*self.arguments]
            return self.map.undistribute(seq, mul_op=mul_op, evaluate=evaluate)
        else:
            return self

scalar_add = AdditionOperator(S.Complexes**2, S.Complexes, S.Zero)
