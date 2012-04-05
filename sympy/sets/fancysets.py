from sympy import Dummy, S, symbols, Lambda, pi, Basic, sympify, ask, Q
from sympy.functions.elementary.integers import floor, ceiling
from sympy.core.compatibility import iterable
from sympy.core.sets import Set, Interval, FiniteSet, Intersection
from sympy.core.singleton import Singleton, S
from sympy.solvers import solve
oo = S.Infinity

class Naturals(Set):
    """
    Represents the Natural Numbers. The Naturals are available as a singleton
    as S.Naturals

    Examples
    ========

        >>> from sympy import S, Interval
        >>> 5 in S.Naturals
        True
        >>> iterable = iter(S.Naturals)
        >>> print iterable.next()
        1
        >>> print iterable.next()
        2
        >>> print iterable.next()
        3
        >>> S.Naturals.intersect(Interval(0, 10))
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
    """

    __metaclass__ = Singleton
    is_iterable = True

    def _intersect(self, other):
        if other.is_Interval:
            return Intersection(S.Integers, other, Interval(1, oo))
        return None

    def _contains(self, other):
        if ask(Q.positive(other)) and ask(Q.integer(other)):
            return True
        return False

    def __iter__(self):
        i = S(1)
        while True:
            yield i
            i = i + 1

    @property
    def _inf(self):
        return S.One

    @property
    def _sup(self):
        return oo

class Integers(Set):
    """
    Represents the Integers. The Integers are available as a singleton
    as S.Integers

    Examples
    ========

        >>> from sympy import S, Interval
        >>> 5 in S.Naturals
        True
        >>> iterable = iter(S.Integers)
        >>> print iterable.next()
        0
        >>> print iterable.next()
        1
        >>> print iterable.next()
        -1
        >>> print iterable.next()
        2

        >>> S.Integers.intersect(Interval(-4, 4))
        {-4, -3, -2, -1, 0, 1, 2, 3, 4}
    """

    __metaclass__ = Singleton
    is_iterable = True

    def _intersect(self, other):
        if other.is_Interval:
            s = FiniteSet(range(ceiling(other.left), floor(other.right) + 1))
            return s.intersect(other) # take out endpoints if open interval
        return None

    def _contains(self, other):
        if ask(Q.integer(other)):
            return True
        return False

    def __iter__(self):
        yield S.Zero
        i = S(1)
        while True:
            yield i
            yield -i
            i = i + 1

    @property
    def _inf(self):
        return -oo

    @property
    def _sup(self):
        return oo

class TransformationSet(Set):
    """
    A set that is a transformation of another through some algebraic expression

    Examples
    --------
    >>> from sympy import Symbol, S, TransformationSet, FiniteSet, Lambda

    >>> x = Symbol('x')
    >>> N = S.Naturals
    >>> squares = TransformationSet(Lambda(x, x**2), N) # {x**2 for x in N}
    >>> 4 in squares
    True
    >>> 5 in squares
    False

    >>> FiniteSet(0,1,2,3,4,5,6,7,9,10).intersect(squares)
    {1, 4, 9}

    >>> square_iterable = iter(squares)
    >>> for i in range(4):
    ...     square_iterable.next()
    1
    4
    9
    16
    """
    def __new__(cls, lamda, base_set):
        return Basic.__new__(cls, lamda, base_set)

    lamda    = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def __iter__(self):
        already_seen = set()
        def iterator():
            for i in self.base_set:
                val = self.lamda(i)
                if val in already_seen:
                    continue
                else:
                    already_seen.add(val)
                    yield val

        return iterator()

    def _is_multivariate(self):
        return len(self.lamda.variables) > 1

    def _contains(self, other):
        L = self.lamda
        if self._is_multivariate():
            solns = solve([expr - val for val, expr in zip(other, L.expr)],
                    L.variables)
        else:
            solns = solve(L.expr - other, L.variables[0])

        for soln in solns:
            try:
                if soln in self.base_set:           return True
            except TypeError:
                if soln.evalf() in self.base_set:   return True
        return False

    @property
    def is_iterable(self):
        return self.base_set.is_iterable
