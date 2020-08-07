from ..structure import AlgebraicStructure
from ..group import AbelianGroup
from ..ring import Ring

__all__ = [
    "Module",
]

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
    >>> op = VectorAdditionOperator(X**2, X, e)
    >>> G = AbelianGroup('G', (X,), (op,))

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

    def add(self, *args, evaluate=False):
        # scalar addition
        if all(self.ring.contains(a) == True for a in args):
            return self.ring.add(*args, evaluate=evaluate)
        # vector addition
        if all(self.group.contains(a) == True for a in args):
            op = self.group.operator
            return op(
                *args, sv_mul=self.smul_op, ss_add=self.ring.add_op, ss_mul=self.ring.mul_op,
                evaluate=evaluate
            )
        raise TypeError("Mismatching argument for module addition")

    def subtract(self, a, b, evaluate=False):
        # scalar subtraction
        if all(self.ring.contains(i) == True for i in [a, b]):
            return self.ring.sub(a, b, evaluate=evaluate)
        # vector subtraction
        if all(self.group.contains(i) == True for i in [a, b]):
            op = self.group.operator
            inv_b = op.inverse_element(b)
            return op(
                a, inv_b,
                sv_mul=self.smul_op, ss_add=self.ring.add_op, ss_mul=self.ring.mul_op,
                evaluate=evaluate
            )
        raise TypeError("Mismatching argument for module subtraction")
    sub = subtract

    def multiply(self, *args, evaluate=False):
        # scalar multiplication
        if all(self.ring.contains(a) == True for a in args):
            return self.ring.mul(*args, evaluate=evaluate)
        # vector addition
        vectors = [a for a in args if self.group.contains(a) == True]
        if len(vectors) == 1:
            coeffs = [a for a in args if not self.group.contains(a) == True]
            if len(coeffs) == 0:
                coeff = self.ring.mul_op.identity
            elif len(coeffs) == 1:
                coeff, = coeffs
            else:
                coeff = self.mul(*coeffs, evaluate=evaluate)
            vector, = vectors
            return self.smul_op(
                coeff, vector,
                vv_add=self.group.operator, ss_add=self.ring.add_op, ss_mul=self.ring.mul_op,
                evaluate=evaluate
            )
        raise TypeError("Mismatching argument for module multiplication")
    mul = multiply

    def power(self, x, n, evaluate=False):
        if self.ring.contains(x) == True:
            return self.ring.pow(x, n, evaluate=evaluate)
        raise TypeError("Mismatching argument for module power")
    pow = power
