from sympy import ceiling, floor, solve, Dummy, S
from sympy.core.sets import *
oo = S.Infinity

class Naturals(CountableSet):
    """
    Represents the Natural Numbers. The Naturals are available as a singleton
    as S.Naturals

    Examples
    ========

        >>> from sympy import S, Interval

    """

    __metaclass__ = Singleton

    def _intersect(self, other):
        if other.is_Interval:
            other = other.intersect(Interval(0,oo))
            return FiniteSet(range(ceiling(other.left), floor(other.right)))
        return None

    def _contains(self, other):
        if other<0:
            return False
        other = sympify(other)
        return sympify(other).is_integer or (other-int(other)) == 0

    def __iter__(self):
        def all_naturals():
            i = 1
            while True:
                yield i
                i = i + 1
        return all_naturals()

    @property
    def _inf(self):
        return 1

    @property
    def _sup(self):
        return oo

class Integers(CountableSet):
    """
    Represents the Integers. The Integers are available as a singleton
    as S.Integers

    Examples
    ========

        >>> from sympy import S, Interval

    """

    __metaclass__ = Singleton

    def _intersect(self, other):
        if other.is_Interval:
            return FiniteSet(range(ceiling(other.left), floor(other.right)))
        return None

    def _contains(self, other):
        other = sympify(other)
        return sympify(other).is_integer or (other-int(other)) == 0

    def __iter__(self):
        def all_ints():
            yield 0
            i = 1
            while True:
                yield i
                yield -i
                i = i + 1
        return all_ints()

    @property
    def _inf(self):
        return -oo

    @property
    def _sup(self):
        return oo

class Isomorphic(Set):
    """
    A set that is isomorphic to another through some algebraic expressions

    Examples
    --------
    >>> from sympy import Isomorphic, S, FiniteSet
    >>> squares = Isomorphic(x**2, x, S.Naturals) # {x**2 for x in Naturals}
    >>> 4 in squares
    True
    >>> 5 in squares
    False
    >>> FiniteSet(0,1,2,3,4,5,6,7,9,10).intersect(squares)
    {1, 4}
    >>> square_iterable = iter(squares)
    >>> for i in range(5):
    ...     square_iterable.next()
    0
    1
    4
    9
    16
    """
    def __new__(cls, expr, variable, base_set):
        return Basic.__new__(cls, expr, variable, base_set)

    expr     = property(lambda self: self.args[0])
    variable = property(lambda self: self.args[1])
    base_set = property(lambda self: self.args[2])

    def __iter__(self):
        return (self.expr.subs(self.variable, i) for i in self.base_set)

    def _contains(self, other):
        solns = solve(self.expr - other, self.variable)
        for soln in solns:
            try:
                if soln in self.base_set:           return True
            except TypeError:
                if soln.evalf() in self.base_set:   return True
        return False

    @property
    def is_iterable(self):
        return self.base_set.is_iterable

class IsomorphicToN(Isomorphic):
    """
    A set that is isomorphic to the natural numbers through an algebraic
    expressions. A countable set.

    Examples
    --------
    >>> from sympy import IsomorphicToN, S, FiniteSet
    >>> squares = IsomorphicToN(x**2, x) # {x**2 for x in Naturals}
    >>> 4 in squares
    True
    >>> 5 in squares
    False

    >>> FiniteSet(0,1,2,3,4,5,6,7,9,10).intersect(squares)
    {1, 4}
    >>> square_iterable = iter(squares)
    >>> for i in range(5):
    ...     square_iterable.next()
    0
    1
    4
    9
    16
    """
    def __new__(cls, expr, variable):
        return Isomorphic.__new__(cls, expr, variable, S.Naturals)

x = Dummy('x')
harmonics = IsomorphicToN(1/x, x)
squares = IsomorphicToN(x**2, x)
