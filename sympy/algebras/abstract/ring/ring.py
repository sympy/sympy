from sympy.assumptions import ask, Q
from ..structure import AlgebraicStructure
from ..group import Monoid, AbelianGroup

__all__ = [
    "Ring", "CommutativeRing",
]

class Ring(AlgebraicStructure):
    """
    A base class for algebraic ring.

    Explanation
    ===========

    Ring is algebraic structure that consists of two binary operation;
    addition and multiplication. Ring must be and abelian group under
    addition, and monoid under multiplication. Also, multiplication must
    be distributive with respect to addition [1].

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets

    operators : tuple of two Maps
        The first one is addition operator, and the
        second is multiplication operator.

    Examples
    ========

    >>> from sympy import Ring, Set, AdditionOperator, BinaryOperator

    Build a purely abstract ring with no number involved.

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

    >>> R = Ring('R', (A,), (add, mul))

    >>> R.add(a, b)
    a + b
    >>> R.add(a, e1, evaluate=True)
    a

    >>> R.sub(a, b)
    a + -b
    >>> R.sub(a, a, evaluate=True)
    e1

    >>> R.mul(a, b)
    a*b
    >>> R.mul(b, e2, evaluate=True)
    b

    Ring knows how to collect the added arguments to multiplication.

    >>> R.add(a, a, evaluate=True)
    (e2 + e2)*a

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Ring_(algebra)

    """
    def __new__(cls, name, sets, operators, **kwargs):
        if not len(sets) == 1:
            raise TypeError("%s must consist of one set." % cls)
        if not (len(operators) == 2 and all(o.arity == 2 for o in operators)):
            raise TypeError("%s must consist of two binary operators." % cls)

        add, mul = operators

        obj = super().__new__(cls, name, sets, operators)

        # this is in __new__ to automatically checks the operators
        obj._add_group = AbelianGroup(name, sets, (add,))
        obj._mul_monoid = Monoid(name, sets, (mul,))
        return obj

    @property
    def add_group(self):
        return self._add_group

    @property
    def addition_operator(self):
        return self.operators[0]
    add_op = addition_operator

    def add(self, *args, evaluate=False):
        return self.add_op(*args, mul_op=self.mul_op, evaluate=evaluate)

    @property
    def subtraction_operator(self):
        return self.add_op.subtraction_operator()
    sub_op = subtraction_operator

    def subtract(self, a, b, evaluate=False):
        return self.sub_op(a, b, mul_op=self.mul_op, evaluate=evaluate)
    sub = subtract

    @property
    def mul_monoid(self):
        return self._mul_monoid

    @property
    def multiplication_operator(self):
        return self.operators[1]
    mul_op = multiplication_operator

    def multiply(self, *args, evaluate=False):
        return self.mul_op(*args, add_op=self.add_op, evaluate=evaluate)
    mul = multiply

    @property
    def power_operator(self):
        return self.mul_op.exponent_operator()
    pow_op = power_operator

    def pow(self, x, n, evaluate=False):
        return self.pow_op(x, n, evaluate=evaluate)

    def distribute(self, mul, evaluate=False):
        return mul.distribute(add_op=self.add_op, evaluate=evaluate)

    def undistribute(self, add, evaluate=False):
        return add.undistribute(mul_op=self.mul_op, evaluate=evaluate)

class CommutativeRing(Ring):
    """
    A base class for commutative ring.

    Explanation
    ===========

    Commutative ring is a ring in which the multiplication operation
    is commutative [1]. 

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets

    operators : tuple of two Maps
        The first one is addition operator, and the
        second is multiplication operator.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Commutative_ring

    """

    def __new__(cls, name, sets, operators, **kwargs):
        add, mul = operators

        if not ask(Q.commutative(mul)):
            raise TypeError("Multiplication operator must be commutative for commutative ring.")

        obj = super().__new__(cls, name, sets, operators)
        return obj
