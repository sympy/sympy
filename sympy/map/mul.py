"""Common multiplication operator"""

from functools import cmp_to_key
import itertools
from sympy.assumptions import ask, Q
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap, isappliedmap
from .operator import BinaryOperator

__all__ = [
    "MultiplicationOperator", "NumericMultiplicationOperator",
    "ScalarMultiplicationOperator", "VectorMultiplicationOperator"
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

    Parameters
    ==========

    domain : ProductSet of Set

    codomain : Set

    identity : identity element

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
    latex_name = r' \times '
    left_divisible = right_divisible = True
    associative = True
    commutative = True

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

        if evaluate:
            args = self.multiplication_process(args)

            if not args:
                return self.identity
            elif len(args) == 1:
                return args[0]

            if ask(Q.commutative(self)):
                args.sort(key=cmp_to_key(Basic.compare))

        # return Multiplication class with processed arguments
        args = Tuple(*args)
        aux = Tuple(*aux)
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
    Class for common numeric multiplication operator.

    Explanation
    ===========

    Unlike ``MultiplicationOperator`` which is purely abstract operation, this class is
    common multiplication operator with its identity element fixed to 1. Also, addition
    operator does not need to be passed when applying the arguments.

    Parameters
    ==========

    domain : ProductSet of Set

    codomain : Set

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
    Class for abstract scalar multiplication operator defined between
    scalar and vector of algebraic module.

    Parameters
    ==========

    domain : ProductSet of Ring and AbelianGroup

    codomain : AbelianGroup

    Examples
    ========

    >>> from sympy import (
    ... Set, AdditionOperator, MultiplicationOperator, Ring,
    ... VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator
    ... )

    Build underlying structures

    >>> S = Set('S') # scalar set
    >>> S_e1, S_e2 = [S.element(i) for i in ('e1', 'e2')]
    >>> ss_add = AdditionOperator(S**2, S, S_e1)
    >>> ss_mul = MultiplicationOperator(S**2, S, S_e2)
    >>> S_ring = Ring('S', (S,), (ss_add, ss_mul))

    >>> V = Set('V') # vector set
    >>> v, V_e = [V.element(i) for i in ('v', 'e')]
    >>> vv_add = VectorAdditionOperator(V**2, V, V_e)
    >>> V_group = AbelianGroup('V', (V,), (vv_add,))

    >>> mul = ScalarMultiplicationOperator(S_ring*V_group, V_group)

    >>> mul(S_e2, v)
    e2*v
    >>> _.doit()
    v

    See Also
    ========

    algebras.Module

    """
    latex_name = ' '
    left_divisible = right_divisible = False
    associative = False

    def __new__(cls, domain, codomain, **kwargs):
        return super(MultiplicationOperator, cls).__new__(cls, domain, codomain, **kwargs)

    @property
    def scalar_ring(self):
        return self.domain.args[0]

    @property
    def vector_group(self):
        return self.domain.args[1]

    @property
    def vv_add(self):
        return self.vector_group.operator

    @property
    def ss_add(self):
        return self.scalar_ring.add_op

    @property
    def ss_mul(self):
        return self.scalar_ring.mul_op

    def __call__(self, a, b, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)
        return Multiplication(self, (a, b), (), **kwargs)

    def apply(self, a, b, aux, **kwargs):
        evaluate = kwargs.get('evaluate', False)
        args = (a, b)

        # order a and b so that a is scalar and b is vector
        a, = [i for i in args if self.scalar_ring.contains(i) == True]
        b, = [i for i in args if self.vector_group.contains(i) == True]

        if evaluate:
            args = self.multiplication_process(args)

            if len(args) == 1:
                return args[0]

            if ask(Q.commutative(self)):
                args.sort(key=cmp_to_key(Basic.compare))

        # return Multiplication class with processed arguments
        args = Tuple(*args)
        aux = Tuple(*aux)
        result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args, aux)
        return result

    def multiplication_process(self, seq):
        scalar, vector = seq
        v_s, v_v = vector.as_coeff_Mul(mul_op=self)
        scalar, vector = self.ss_mul(scalar, v_s, add_op=self.ss_add, evaluate=True), v_v
        if scalar == self.ss_mul.identity:
            return [vector]
        return [scalar, vector]

    def distribute(self, seq, evaluate=False):
        # not used in construction algorithm because it is
        # discouraged according to core/parameters
        vv_add, ss_add, ss_mul = self.vv_add, self.ss_add, self.ss_mul

        terms = []
        for o in seq:
            if isappliedmap(o, (vv_add, ss_add)):
                terms.append(o.arguments)
            else:
                terms.append((o,))
        terms = [
            self(*l, evaluate=evaluate)
            for l in itertools.product(*terms)
        ]
        return vv_add(
            *terms,
            sv_mul=self, ss_add=ss_add, ss_mul=ss_mul,
            evaluate=evaluate
        )

    def _eval_as_coeff_Mul(self, expr, **kwargs):
        if isappliedmap(expr, self):
            return expr.arguments
        ss_mul = self.ss_mul
        return ss_mul.identity, expr

class VectorMultiplicationOperator(MultiplicationOperator):
    r"""
    Class for bilinear product operator between two vectors which is
    defined in an algebra over a vector space [1].

    Parameters
    ==========

    domain : ProductSet of two VectorSpace

    codomain : VectorSpace

    Examples
    ========

    >>> from sympy import (
    ... Set, AdditionOperator, MultiplicationOperator, Field,
    ... VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, VectorSpace,
    ... VectorMultiplicationOperator
    ... )

    Build underlying vector space

    >>> S = Set('S') # scalar set
    >>> S_e1, S_e2 = [S.element(i) for i in ('e1', 'e2')]
    >>> ss_add = AdditionOperator(S**2, S, S_e1)
    >>> ss_mul = MultiplicationOperator(S**2, S, S_e2)
    >>> S_field = Field('S', (S,), (ss_add, ss_mul))

    >>> V = Set('V') # vector set
    >>> v, V_e = [V.element(i) for i in ('v', 'e')]
    >>> vv_add = VectorAdditionOperator(V**2, V, V_e)
    >>> V_group = AbelianGroup('V', (V,), (vv_add,))

    >>> sv_mul = ScalarMultiplicationOperator(S_field*V_group, V_group)
    >>> VS = VectorSpace('VS', (S_field, V_group), (sv_mul,))

    >>> vv_mul = VectorMultiplicationOperator(VS*VS, V_group)

    >>> vv_mul(v, v)
    v*v

    See Also
    ========

    algebras.Algebra

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Algebra_over_a_field

    """

    left_divisible = right_divisible = False
    associative = None
    commutative = None
    identity = None # can be overridden to construct unital algebra

    def __new__(cls, domain, codomain, **kwargs):
        return super(MultiplicationOperator, cls).__new__(cls, domain, codomain, **kwargs)

    @property
    def vectorspace(self):
        return self.domain.args[0]

    @property
    def vv_add(self):
        return self.vectorspace.group.operator

    @property
    def ss_add(self):
        return self.vectorspace.field.add_op

    @property
    def sv_mul(self):
        return self.vectorspace.smul_op

    @property
    def ss_mul(self):
        return self.vectorspace.field.mul_op

    def __call__(self, *args, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)
        return Multiplication(self, args, (), **kwargs)

    def apply(self, *args, aux, **kwargs):
        evaluate = kwargs.get('evaluate', False)

        if evaluate:
            args = self.multiplication_process(args)

            if len(args) == 1:
                return args[0]

            if ask(Q.commutative(self)):
                args.sort(key=cmp_to_key(Basic.compare))

        # return Multiplication class with processed arguments
        args = Tuple(*args)
        aux = Tuple(*aux)
        result = super(AppliedMap, Multiplication).__new__(Multiplication, self, args, aux)
        return result

    def multiplication_process(self, seq):
        seq = self.remove_identity(seq)
        return seq

    def distribute(self, seq, evaluate=False):
        # not used in construction algorithm because it is
        # discouraged according to core/parameters
        vv_add, ss_add, sv_add, ss_mul = self.vv_add, self.ss_add, self.sv_mul, self.ss_mul

        coeffs = []
        vectors = []
        for o in seq:
            if isappliedmap(o, sv_mul):
                c, v = o.as_coeff_Mul(mul_op=sv_mul)
                coeffs.append(c)
                vectors.append((v,))
            elif isappliedmap(o, vv_add):
                vectors.append(o.arguments)
            else:
                vectors.append((o,))
        terms = [
            self(*l, evaluate=evaluate)
            for l in itertools.product(*vectors)
        ]
        vector = vv_add(
            *terms,
            sv_mul=sv_mul, ss_add=ss_add, ss_mul=ss_mul,
            evaluate=evaluate
        )

        if len(coeffs) == 0:
            return vector
        elif len(coeffs) == 1:
            coeff, = coeffs
            return sv_mul(
                coeff, vector,
                vv_add=vv_add, ss_add=ss_add, ss_mul=ss_mul,
                evaluate=evaluate
            )

    def _eval_as_coeff_Mul(self, expr, **kwargs):
        ss_mul = self.ss_mul
        return ss_mul.identity, expr

class Multiplication(AppliedMap):
    """
    Class for the unevaluated result of general multiplication.

    Parameters
    ==========

    mapping : Map

    args : tuple of arguments
        Arguments applied to *mapping*

    aux : tuple of Maps
        Auxillary operators

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    """
    def __new__(cls, mapping, args, aux, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)

        args = [_sympify(a) for a in args]

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
