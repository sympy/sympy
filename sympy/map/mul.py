"""Common multiplication operator"""

from functools import cmp_to_key
import itertools
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap, isappliedmap
from .operator import BinaryOperator

__all__ = [
    "MultiplicationOperator", "NumericMultiplicationOperator",
    "ScalarMultiplicationOperator",
    "Multiplication",
]

class MultiplicationOperator(BinaryOperator):
    """
    Class for general scalar multiplication operator for algebraic field, which multiplies
    a scalar and a vector.

    Explanation
    ===========

    Multiplication operater forms an abelian group over field structure, and is distributive
    over addition.
    When ``MultiplicationOperator`` instance is called, addition operator must be passed
    as *add_op* parameter for this relation. By using ``Field`` class from algebras module,
    this is done automatically.

    Examples
    ========

    >>> from sympy import AdditionOperator, MultiplicationOperator, Set

    >>> A = Set('A')
    >>> a, b, e1, e2 = [A.element(n) for n in ('a', 'b', 'e1', 'e2')]

    >>> add = AdditionOperator(A**2, A, e1)
    >>> mul = MultiplicationOperator(A**2, A, e2)

    >>> mul(a, b, add_op=add)
    a*b
    >>> mul(a, a, add_op=add, evaluate=True)
    a**2

    Constructing a Field makes applying the operator easier

    >>> from sympy import Field
    >>> F = Field('F', (A,), (add, mul))

    >>> F.mul(a, a, evaluate=True)
    a**2

    See Also
    ========

    algebras.Field

    """
    str_name = '*'
    pretty_name = '⋅'
    latex_name = r' \times '
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

    def __call__(self, *args, add_op, evaluate=False):
        return Multiplication(self, args, (add_op,), evaluate=evaluate)

    def apply(self, *args, aux, **kwargs):
        evaluate = kwargs.get('evaluate', False)

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
        aux = Tuple(*[_sympify(a) for a in aux])
        result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args, aux)
        return result

    def multiplication_process(self, seq):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        seq = self.assoc_comm_process(seq)
        return seq

    def distribute(self, seq, add_op, evaluate=False):
        # x*(y+z) -> x*y + x*z
        # not used in construction algorithm because it is
        # discouraged according to core/parameters
        terms = []
        for o in seq:
            if isappliedmap(o, add_op):
                terms.append(o.arguments)
            else:
                terms.append((o,))
        terms = [self(*l, add_op=add_op, evaluate=evaluate) for l in itertools.product(*terms)]
        return add_op(*terms, mul_op=self, evaluate=evaluate)

    def _ldiv_apply(self, ldivop, divisor, dividend, **kwargs):
        inv_divisor = self.exponent_operator()(divisor, -1, **kwargs)
        return self(inv_divisor, dividend, **kwargs)

    def _rdiv_apply(self, rdivop, dividend, divisor, **kwargs):
        inv_divisor = self.exponent_operator()(divisor, -1, **kwargs)
        return self(dividend, inv_divisor, **kwargs)

class NumericMultiplicationOperator(MultiplicationOperator):
    """
    Class for common numeric multiplication operator

    Explanation
    ===========

    Unlike ``MultiplicationOperator`` which is purely abstract operation, this class is
    common multiplication operator with its identity element fixed to 1. Also, addition
    operator does not need to be passed when applying the arguments.

    Examples
    ========

    >>> from sympy import NumericMultiplicationOperator, S, symbols
    >>> a, b = symbols('a b', real=True)

    >>> mul = NumericMultiplicationOperator(S.Reals**2, S.Reals)

    >>> mul(2, 3, evaluate=True)
    6
    >>> mul(a, b, evaluate=True)
    a*b
    >>> mul(a, 1, evaluate=True)
    a
    >>> mul(a, a, evaluate=True)
    a**2

    """
    def __new__(cls, domain, codomain, **kwargs):
        domain, codomain = _sympify(domain), _sympify(codomain)
        return super(BinaryOperator, cls).__new__(cls, domain, codomain, **kwargs)

    @property
    def identity(self):
        return S.One

    def addition_operator(self):
        from .add import NumericAdditionOperator
        return NumericAdditionOperator(self.domain, self.codomain)
    add_op = addition_operator

    def __call__(self, *args, evaluate=False, **kwargs):
        kwargs.pop('add_op', None)
        add_op = self.add_op()
        return Multiplication(self, args, (add_op,), evaluate=evaluate, **kwargs)

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

class ScalarMultiplicationOperator(MultiplicationOperator):
    r"""
    Class for abstract scalar multiplication operator defined in algebraic module.

    Explanation
    ===========

    To be called, scalar multiplication operator needs keyword arguments *vv_add*
    (vector-vector addition operator), *ss_add* (scalar-scalar addition operator)
    and *ss_mul(scalar-scalar multiplication operator). This is to perform the
    distribution of vectors, i.e. $\left(1+1 \right) \mathbf(v) = \mathbf(v) + \mathbf(v)$ and
    $a \left(\mathbf(v) + \mathbf(w) \right) = a \mathbf(v) + a \mathbf(w)$.

    Examples
    ========

    >>> from sympy import (
    ... Set, VectorAdditionOperator, AdditionOperator, MultiplicationOperator,
    ... ScalarMultiplicationOperator
    ... )

    >>> S, V = Set('S'), Set('V') # scalar set and vector set
    >>> S_e1, S_e2 = [S.element(i) for i in ('e1', 'e2')]
    >>> v, V_e = [V.element(i) for i in ('v', 'e')]

    >>> vv_add = VectorAdditionOperator(V**2, V, V_e)
    >>> ss_add = AdditionOperator(S**2, S, S_e1)
    >>> ss_mul = MultiplicationOperator(S**2, S, S_e2)

    >>> mul = ScalarMultiplicationOperator(S*V, V)

    >>> mul(S_e2, v, vv_add=vv_add, ss_add=ss_add, ss_mul=ss_mul)
    e2*v
    >>> _.doit()
    v

    Constructing a Module makes applying the operator easier

    >>> from sympy import Ring, AbelianGroup, Module
    >>> S_ring = Ring('S', (S,), (ss_add, ss_mul))
    >>> V_group = AbelianGroup('V', (V,), (vv_add,))
    >>> M = Module('M', (S_ring, V_group), (mul,))

    >>> M.mul(S_e2, v, evaluate=True)
    v

    """
    str_name = '*'
    pretty_name = ' '
    latex_name = ' '
    is_left_divisible = is_right_divisible = False
    is_associative = False

    def __new__(cls, domain, codomain, **kwargs):
        domain, codomain = _sympify(domain), _sympify(codomain)
        return super(MultiplicationOperator, cls).__new__(cls, domain, codomain, **kwargs)

    @property
    def scalars(self):
        return self.domain.args[0]

    @property
    def vectors(self):
        return self.domain.args[1]

    def __call__(self, a, b, vv_add, ss_add, ss_mul, evaluate=False, **kwargs):
        # vv_add : vector-vector add
        # ss_add : scalar-scalar add
        # ss_mul : scalar-scalar mul
        return Multiplication(
            self, (a, b), (vv_add, ss_add, ss_mul), evaluate=evaluate, **kwargs
        )

    def apply(self, a, b, aux, **kwargs):
        vv_add, ss_add, ss_mul = aux
        evaluate = kwargs.get('evaluate', False)

        # sympify the arguments
        args = [_sympify(i) for i in (a, b)]

        a, = [i for i in args if self.scalars.contains(i) == True]
        b, = [i for i in args if self.vectors.contains(i) == True]

        if evaluate:
            args = self.multiplication_process(args, vv_add, ss_add, ss_mul)

            if len(args) == 1:
                return args[0]

        # return Multiplication class with processed arguments
        args = Tuple(*[_sympify(a) for a in args])
        aux = Tuple(*[_sympify(a) for a in aux])
        result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args, aux)
        return result

    def distribute(self, seq, vv_add, ss_add, ss_mul, evaluate=False):
        # not used in construction algorithm because it is
        # discouraged according to core/parameters
        terms = []
        for o in seq:
            if isappliedmap(o, (vv_add, ss_add)):
                terms.append(o.arguments)
            else:
                terms.append((o,))
        terms = [
            self(*l, vv_add=vv_add, ss_add=ss_add, ss_mul=ss_mul, evaluate=evaluate)
            for l in itertools.product(*terms)
        ]
        return vv_add(*terms, sv_mul=self, ss_add=ss_add, ss_mul=ss_mul, evaluate=evaluate)

    def multiplication_process(self, seq, vv_add, ss_add, ss_mul):
        scalar, vector = seq
        if scalar == ss_mul.identity:
            return [vector]
        return seq

    def _eval_as_coeff_Mul(self, expr, **kwargs):
        if isappliedmap(expr, self):
            return expr.arguments
        ss_mul = kwargs.get('ss_mul')
        return ss_mul.identity, expr

class Multiplication(AppliedMap):
    """
    Class for the unevaluated result of general multiplication.

    Parameters
    ==========

    map : Map

    args : tuple of arguments
        Arguments applied to *map*

    aux : tuple of Maps
        Auxillary operators

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    """
    def __new__(cls, mapping, args, aux, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)

        # consult mapping.apply
        result = mapping.apply(*args, aux=aux, **kwargs)

        # check codomain
        if mapping.codomain.contains(result) == False:
                    raise TypeError(
                "%s is not in %s's codomain %s." % (result, mapping, mapping.codomain)
                )

        return result

    @property
    def aux(self):
        return self.args[2]

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args, self.aux, **kwargs)

    def as_coeff_Mul(self, rational=False, mul_op=None, **kwargs):
        if mul_op is None:
            mul_op = self.map
        kwargs.update(rational=rational, mul_op=mul_op)
        return mul_op._eval_as_coeff_Mul(self, **kwargs)

def _coeff_isneg(a):
    """
    Return True if the leading number is negative.

    Examples
    ========

    >>> from sympy import NumericMultiplicationOperator, S, symbols
    >>> from sympy.map.mul import _coeff_isneg
    >>> a, b = symbols('a b', real=True)
    >>> mul = NumericMultiplicationOperator(S.Reals**2, S.Reals)

    >>> _coeff_isneg(S(-3))
    True
    >>> _coeff_isneg(mul(a, b))
    False
    >>> _coeff_isneg(mul(-3, a))
    True

    """
    if isinstance(a, Multiplication):
        a = a.arguments[0]
        return _coeff_isneg(a)
    if a.is_Mul:
        a = a.args[0]
        return _coeff_isneg(a)
    return a.is_Number and a.is_extended_negative
