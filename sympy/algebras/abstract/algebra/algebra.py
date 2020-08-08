from ..structure import AlgebraicStructure
from ..module import Module

__all__ = [
    'Algebra',
]

class Algebra(AlgebraicStructure):
    """
    Class for algebra over a ring.

    Explanation
    ===========

    Algebra (over a ring) is algebraic structure over a module,
    where an additional operator for vector-vector multiplication
    is introduced [1].

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of one Module

    operators : tuple of one Map
        This map is vector multiplication operator.

    Examples
    ========

    >>> from sympy import (
    ... Set, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, S, VectorSpace,
    ... VectorMultiplicationOperator, Algebra, 
    ... )

    >>> F = S.RealsField

    >>> X = Set('X')
    >>> x, v, e = [X.element(i) for i in 'xve']
    >>> vadd = VectorAdditionOperator(X**2, X, e)
    >>> G = AbelianGroup('G', (X,), (vadd,))

    >>> smul = ScalarMultiplicationOperator(F*G, G)
    >>> V = VectorSpace('V', (F, G), (smul,))

    >>> vv_mul = VectorMultiplicationOperator(V*V, G)
    >>> A = Algebra('A', (V,), (vv_mul,))

    >>> A.mul(v, v)
    v*v

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Algebra_over_a_field

    """
    def __new__(cls, name, sets, operators, **kwargs):
        if not len(sets) == 1 and isinstance(sets[0], Module):
            raise TypeError("%s must consist of one underlying module." % cls)

        if not (len(operators) == 1 and all(o.arity == 2 for o in operators)):
            raise TypeError("%s must consist of one binary vector multiplication operator." % cls)

        return super().__new__(cls, name, sets, operators)

    @property
    def module(self):
        return self.args[1].args[0]

    @property
    def vector_multiplication_operator(self):
        return self.args[2].args[0]
    vmul_op = vector_multiplication_operator

    def check_vector(self, a):
        return self.module.check_vector(a)

    def check_scalar(self, a):
        return self.module.check_scalar(a)

    def add(self, *args, evaluate=False):
        return self.module.add(*args, evaluate=evaluate)

    def subtract(self, a, b, evaluate=False):
        return self.module.subtract(a, b, evaluate=evaluate)
    sub = subtract

    def multiply(self, *args, evaluate=False):
        vmul_op = self.vmul_op
        vectors = [i for i in args if self.check_vector(i)]

        # vector multiplication
        if len(vectors) > 1:
            vmul_op = self.vmul_op
            scalars = [i for i in args if self.check_scalar(i)]

            if len(scalars) == 0:
                coeff = self.module.ring.mul_op.identity
            elif len(scalars) == 1:
                coeff, = scalars
            else:
                coeff = self.module.mul(*scalars, evaluate=evaluate)

            vector = vmul_op(*vectors, evaluate=evaluate)
            return self.module.mul(coeff, vector, evaluate=True)

        else:
            return self.module.mul(*args, evaluate=evaluate)
    mul = multiply

    def power(self, x, n, evaluate=False):
        if self.check_vector(x):
            powop = self.vmul_op.exponent_operator()
            return powop(x, n, evaluate=evaluate)

        else:
            return self.module.power(x, n, evaluate=evaluate)
    pow = power

    def divide(self, a, b, **kwargs):
        if self.check_vector(b):
            inv_b = self.vmul_op.inverse_element(b, evaluate=True)
            return self.mul(a, inv_b, evaluate=True)

        else:
            return self.module.div(a, b, **kwargs)
    div = divide
