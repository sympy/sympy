"""Common multiplication operator"""

from functools import cmp_to_key
import itertools
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap, isappliedmap
from .operator import BinaryOperator

__all__ = [
    "MultiplicationOperator", "Multiplication",
    "scalar_mul", "scalar_pow", "scalar_divide",
]

class MultiplicationOperator(BinaryOperator):
    """
    Class for general scalar multiplication operator.

    Explanation
    ===========

    When defined as group operator, multiplication operator has no special feature. However,
    when defined as ring operator, multiplication is paired with addition and distributive
    with respect to it.
    When ``MultiplicationOperator`` is called, addition operator can be passed as *add_op*
    parameter for this relation. This implies that two operators form a ring over their
    domain.
    If you use ``Ring`` class from algebras module, two operators are automatically related.

    Examples
    ========

    >>> from sympy import scalar_mul
    >>> from sympy.abc import x, y

    >>> scalar_mul(x, y)
    x * y

    >>> scalar_mul(x, x)
    x * x
    >>> scalar_mul(x, x, evaluate=True)
    x**2

    """
    str_name = '*'
    pretty_name = '⋅'
    latex_name = r'\times'
    is_left_divisible = is_right_divisible = True
    is_associative = True
    is_commutative = True

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

    def __call__(self, *args, add_op=None, evaluate=False):
        return Multiplication(self, args, add_op, evaluate=evaluate)

    def apply(self, *args, **kwargs):
        evaluate = kwargs.get('evaluate', False)
        add_op = kwargs.get('add_op', None)

        # sympify the arguments
        args = [_sympify(a) for a in args]

        if evaluate:
            args = self.multiplication_process(args)

            if not args:
                return self.identity
            elif len(args) == 1:
                return args[0]

        # return Multiplication class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])

        if add_op is None:
            result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args)
        else:
            result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args, add_op)
        return result

    def gather_num(self, seq):
        num = self.identity
        newseq = []
        for o in seq:
            if o.is_Number:
                num *= o
            else:
                newseq.append(o)
        if num != self.identity:
            newseq.insert(0, num)
        return newseq

    def multiplication_process(self, seq):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        seq = self.assoc_comm_process(seq)
        seq = self.gather_num(seq)
        return seq

    def distribute(self, seq, add_op, evaluate=False):
        # x*(y+z) -> x*y + x*z
        # not used in construction algorithm because it is
        # discouraged in core/parameters
        terms = []
        for o in seq:
            if isappliedmap(o, add_op):
                terms.append(o.arguments)
            else:
                terms.append((o,))
        terms = [self(*l, evaluate=evaluate) for l in itertools.product(*terms)]
        return add_op(*terms, evaluate=evaluate)


class Multiplication(AppliedMap):
    """
    Class for the unevaluated result of general multiplication.

    """
    def __new__(cls, mapping, args, add_op=None, evaluate=False, **kwargs):
        kwargs.update(
            evaluate=evaluate,
            add_op=add_op
        )

        return super().__new__(cls, mapping, args, **kwargs)

    @property
    def add_op(self):
        if len(self.args) < 3:
            return None
        return self.args[2]

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args)

    def as_coeff_Mul(self, rational=False, mul_op=None):
        if mul_op is None or isappliedmap(self, mul_op):
            coeff, args = self.arguments[0], self.arguments[1:]
            if coeff.is_Number:
                if not rational or coeff.is_Rational:
                    if len(args) == 1:
                        return coeff, args[0]
                    else:
                        return coeff, self._new_rawargs(*args)
                elif coeff.is_extended_negative:
                    return -self.map.identity, self._new_rawargs(*((-coeff,) + args))
        return self.map.identity, self

    def distribute(self, add_op=None, evaluate=False):
        # x*(y+z) -> x*y + x*z
        # not used in construction algorithm because it is
        # discouraged according to core/parameters
        if add_op is None:
            add_op = self.add_op

        if add_op is not None:
            seq = [*self.arguments]
            return self.map.distribute(seq, add_op=add_op, evaluate=evaluate)
        else:
            return self

scalar_mul = MultiplicationOperator(S.Complexes**2, S.Complexes, S.One)
scalar_pow = scalar_mul.exponent_operator()
scalar_divide = scalar_mul.right_division()
