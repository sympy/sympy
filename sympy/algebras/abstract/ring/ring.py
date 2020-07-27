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

    >>> from sympy import Ring, S, scalar_add, scalar_mul
    >>> from sympy.abc import x, y
    >>> R = Ring('R', (S.Complexes,), (scalar_add, scalar_mul))

    >>> R.add(x, y)
    x + y
    >>> R.mul(x, y)
    x * y

    Ring knows how to collect the added arguments to multiplication.

    >>> R.add(x, x, evaluate=True)
    2 * x

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

        # this automatically checks the operators
        obj._add_group = AbelianGroup(name, sets, (add,))
        obj._mul_monoid = Monoid(name, sets, (mul,))
        return obj

    @property
    def add_op(self):
        return self.operators[0]

    @property
    def mul_op(self):
        return self.operators[1]

    @property
    def add_group(self):
        return self._add_group

    @property
    def mul_monoid(self):
        return self._mul_monoid

    def addition(self, *args, evaluate=False):
        return self.add_op(*args, mul_op=self.mul_op, evaluate=evaluate)
    add = addition

    def multiplication(self, *args, evaluate=False):
        return self.mul_op(*args, evaluate=evaluate)
    mul = multiplication

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
