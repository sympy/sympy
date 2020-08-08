from ..structure import AlgebraicStructure
from ..group import AbelianGroup
from ..ring import Ring

__all__ = [
    "Module",
]

### TODO: Implement LeftModule and RightModule

class Module(AlgebraicStructure):
    """
    A base class for algebraic module.

    Explanation
    ===========

    Module is composite algebraic structure which consists of a ring $R$ that is
    called scalar, and an abelian group $M$. Scalar multiplication operation
    between $R$ and $M$ is defined [1].

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of two structures
        The first one is scalar ring, and the
        second one is abelian group.

    operators : tuple of one Map
        This map is scalar multiplication operator.

    Examples
    ========

    >>> from sympy import (
    ... Set, VectorAdditionOperator, AbelianGroup,
    ... ScalarMultiplicationOperator, S, Module
    ... )

    >>> R = S.IntegersRing

    >>> X = Set('X')
    >>> x, y, e = [X.element(i) for i in 'xye']
    >>> vadd = VectorAdditionOperator(X**2, X, e)
    >>> G = AbelianGroup('G', (X,), (vadd,))

    >>> smul = ScalarMultiplicationOperator(R*G, G)

    >>> M = Module('M', (R, G), (smul,))

    >>> M.add(2, 3, evaluate=True)
    5
    >>> M.add(x, x, evaluate=True)
    2*x

    >>> M.sub(2, 3, evaluate=True)
    -1
    >>> M.sub(x, y, evaluate=True)
    (-y) + x

    >>> M.mul(2, 3, evaluate=True)
    6
    >>> M.mul(3, x, evaluate=True)
    3*x

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Module_(mathematics)

    """
    def __new__(cls, name, sets, operators, **kwargs):
        if not len(sets) == 2:
            raise TypeError("%s must consist of two underlying structures." % cls)
        if not (len(operators) == 1 and all(o.arity == 2 for o in operators)):
            raise TypeError("%s must consist of one binary scalar multiplication operator." % cls)

        scalar_ring, abelian_group = sets

        if not isinstance(scalar_ring, Ring):
            raise TypeError("%s is not ring." % scalar_ring)
        if not isinstance(abelian_group, AbelianGroup):
            raise TypeError("%s is not abelian group." % abelian_group)

        return super().__new__(cls, name, sets, operators)

    @property
    def ring(self):
        return self.args[1].args[0]

    @property
    def group(self):
        return self.args[1].args[1]

    @property
    def scalar_multiplication_operator(self):
        return self.args[2].args[0]
    smul_op = scalar_multiplication_operator

    def check_vector(self, a):
        """
        Return ``True`` if a is vector in *self*.
        """
        return self.group.contains(a) == True

    def check_scalar(self, a):
        """
        Return ``True`` if a is scalar in *self*.
        """
        return self.ring.contains(a) == True

    def add(self, *args, evaluate=False):
        # scalar addition
        if all(self.check_scalar(a) == True for a in args):
            return self.ring.add(*args, evaluate=evaluate)

        # vector addition
        if all(self.check_vector(a) == True for a in args):
            op = self.group.operator
            return op(
                *args, sv_mul=self.smul_op, ss_add=self.ring.add_op, ss_mul=self.ring.mul_op,
                evaluate=evaluate
            )
        raise TypeError("Mismatching argument for module addition")

    def subtract(self, a, b, evaluate=False):
        # scalar subtraction
        if all(self.check_scalar(i) == True for i in [a, b]):
            return self.ring.sub(a, b, evaluate=evaluate)

        # vector subtraction
        if all(self.check_vector(i) == True for i in [a, b]):
            sv_mul, ss_add, ss_mul = self.smul_op, self.ring.add_op, self.ring.mul_op
            op = self.group.operator
            inv_b = op.inverse_element(b, sv_mul=sv_mul, ss_add=ss_add, ss_mul=ss_mul)
            return op(
                a, inv_b,
                sv_mul=sv_mul, ss_add=ss_add, ss_mul=ss_mul,
                evaluate=evaluate
            )
        raise TypeError("Mismatching argument for module subtraction")
    sub = subtract

    def multiply(self, *args, evaluate=False):
        # scalar multiplication
        if all(self.check_scalar(a) == True for a in args):
            return self.ring.mul(*args, evaluate=evaluate)

        # vector multiplication
        vectors = [a for a in args if self.check_vector(a) == True]
        if len(vectors) == 1:

            coeffs = [a for a in args if not self.check_vector(a) == True]
            if len(coeffs) == 0:
                coeff = self.ring.mul_op.identity
            elif len(coeffs) == 1:
                coeff, = coeffs
            else:
                coeff = self.mul(*coeffs, evaluate=evaluate)

            vector, = vectors
            return self.smul_op(coeff, vector, evaluate=evaluate)
        raise TypeError("Mismatching argument for module multiplication")
    mul = multiply

    def power(self, x, n, evaluate=False):
        if self.check_scalar(x) == True:
            return self.ring.pow(x, n, evaluate=evaluate)
        raise TypeError("Mismatching argument for module power")
    pow = power
