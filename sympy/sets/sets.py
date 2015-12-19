from __future__ import print_function, division

from itertools import product

from sympy.core.sympify import _sympify, sympify
from sympy.core.basic import Basic
from sympy.core.singleton import Singleton, S
from sympy.core.evalf import EvalfMixin
from sympy.core.numbers import Float
from sympy.core.compatibility import iterable, with_metaclass, ordered, range
from sympy.core.evaluate import global_evaluate
from sympy.core.decorators import deprecated
from sympy.core.mul import Mul
from sympy.core.relational import Eq
from sympy.core.symbol import Symbol
from sympy.sets.contains import Contains
from sympy.utilities.misc import func_name

from mpmath import mpi, mpf
from sympy.logic.boolalg import And, Or, Not, true, false
from sympy.utilities import subsets


class Set(Basic):
    """
    The base class for any kind of set.

    This is not meant to be used directly as a container of items. It does not
    behave like the builtin ``set``; see :class:`FiniteSet` for that.

    Real intervals are represented by the :class:`Interval` class and unions of
    sets by the :class:`Union` class. The empty set is represented by the
    :class:`EmptySet` class and available as a singleton as ``S.EmptySet``.
    """
    is_number = False
    is_iterable = False
    is_interval = False

    is_FiniteSet = False
    is_Interval = False
    is_ProductSet = False
    is_Union = False
    is_Intersection = None
    is_EmptySet = None
    is_UniversalSet = None
    is_Complement = None
    is_ComplexRegion = False

    @staticmethod
    def _infimum_key(expr):
        """
        Return infimum (if possible) else S.Infinity.
        """
        try:
            infimum = expr.inf
            assert infimum.is_comparable
        except (NotImplementedError,
                AttributeError, AssertionError, ValueError):
            infimum = S.Infinity
        return infimum

    def union(self, other):
        """
        Returns the union of 'self' and 'other'.

        Examples
        ========

        As a shortcut it is possible to use the '+' operator:

        >>> from sympy import Interval, FiniteSet
        >>> Interval(0, 1).union(Interval(2, 3))
        [0, 1] U [2, 3]
        >>> Interval(0, 1) + Interval(2, 3)
        [0, 1] U [2, 3]
        >>> Interval(1, 2, True, True) + FiniteSet(2, 3)
        (1, 2] U {3}

        Similarly it is possible to use the '-' operator for set differences:

        >>> Interval(0, 2) - Interval(0, 1)
        (1, 2]
        >>> Interval(1, 3) - FiniteSet(2)
        [1, 2) U (2, 3]

        """
        return Union(self, other)

    def intersect(self, other):
        """
        Returns the intersection of 'self' and 'other'.

        >>> from sympy import Interval

        >>> Interval(1, 3).intersect(Interval(1, 2))
        [1, 2]

        >>> from sympy import imageset, Lambda, symbols, S
        >>> n, m = symbols('n m')
        >>> a = imageset(Lambda(n, 2*n), S.Integers)
        >>> a.intersect(imageset(Lambda(m, 2*m + 1), S.Integers))
        EmptySet()

        """
        return Intersection(self, other)

    def intersection(self, other):
        """
        Alias for :meth:`intersect()`
        """
        return self.intersect(other)

    def _intersect(self, other):
        """
        This function should only be used internally

        self._intersect(other) returns a new, intersected set if self knows how
        to intersect itself with other, otherwise it returns ``None``

        When making a new set class you can be assured that other will not
        be a :class:`Union`, :class:`FiniteSet`, or :class:`EmptySet`

        Used within the :class:`Intersection` class
        """
        return None

    def is_disjoint(self, other):
        """
        Returns True if 'self' and 'other' are disjoint

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 2).is_disjoint(Interval(1, 2))
        False
        >>> Interval(0, 2).is_disjoint(Interval(3, 4))
        True

        References
        ==========

        .. [1] http://en.wikipedia.org/wiki/Disjoint_sets
        """
        return self.intersect(other) == S.EmptySet

    def isdisjoint(self, other):
        """
        Alias for :meth:`is_disjoint()`
        """
        return self.is_disjoint(other)

    def _union(self, other):
        """
        This function should only be used internally

        self._union(other) returns a new, joined set if self knows how
        to join itself with other, otherwise it returns ``None``.
        It may also return a python set of SymPy Sets if they are somehow
        simpler. If it does this it must be idempotent i.e. the sets returned
        must return ``None`` with _union'ed with each other

        Used within the :class:`Union` class
        """
        return None

    def complement(self, universe):
        """
        The complement of 'self' w.r.t the given the universe.

        Examples
        ========

        >>> from sympy import Interval, S
        >>> Interval(0, 1).complement(S.Reals)
        (-oo, 0) U (1, oo)

        >>> Interval(0, 1).complement(S.UniversalSet)
        UniversalSet() \ [0, 1]

        """
        return Complement(universe, self)

    def _complement(self, other):
        # this behaves as other - self
        if isinstance(other, ProductSet):
            # For each set consider it or it's complement
            # We need at least one of the sets to be complemented
            # Consider all 2^n combinations.
            # We can conveniently represent these options easily using a
            # ProductSet

            # XXX: this doesn't work if the dimentions of the sets isn't same.
            # A - B is essentially same as A if B has a different
            # dimentionality than A
            switch_sets = ProductSet(FiniteSet(o, o - s) for s, o in
                                     zip(self.sets, other.sets))
            product_sets = (ProductSet(*set) for set in switch_sets)
            # Union of all combinations but this one
            return Union(p for p in product_sets if p != other)

        elif isinstance(other, Interval):
            if isinstance(self, Interval) or isinstance(self, FiniteSet):
                return Intersection(other, self.complement(S.Reals))

        elif isinstance(other, Union):
            return Union(o - self for o in other.args)

        elif isinstance(other, Complement):
            return Complement(other.args[0], Union(other.args[1], self), evaluate=False)

        elif isinstance(other, EmptySet):
            return S.EmptySet

        elif isinstance(other, FiniteSet):
            return FiniteSet(*[el for el in other if self.contains(el) != True])

    def symmetric_difference(self, other):
        return SymmetricDifference(self, other)

    def _symmetric_difference(self, other):
        return Union(Complement(self, other), Complement(other, self))

    @property
    def inf(self):
        """
        The infimum of 'self'

        Examples
        ========

        >>> from sympy import Interval, Union
        >>> Interval(0, 1).inf
        0
        >>> Union(Interval(0, 1), Interval(2, 3)).inf
        0

        """
        return self._inf

    @property
    def _inf(self):
        raise NotImplementedError("(%s)._inf" % self)

    @property
    def sup(self):
        """
        The supremum of 'self'

        Examples
        ========

        >>> from sympy import Interval, Union
        >>> Interval(0, 1).sup
        1
        >>> Union(Interval(0, 1), Interval(2, 3)).sup
        3

        """
        return self._sup

    @property
    def _sup(self):
        raise NotImplementedError("(%s)._sup" % self)

    def contains(self, other):
        """
        Returns True if 'other' is contained in 'self' as an element.

        As a shortcut it is possible to use the 'in' operator:

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1).contains(0.5)
        True
        >>> 0.5 in Interval(0, 1)
        True

        """
        other = sympify(other, strict=True)
        ret = sympify(self._contains(other))
        if ret is None:
            ret = Contains(other, self, evaluate=False)
        return ret

    def _contains(self, other):
        raise NotImplementedError("(%s)._contains(%s)" % (self, other))

    @deprecated(useinstead="is_subset", issue=7460, deprecated_since_version="0.7.6")
    def subset(self, other):
        """
        Returns True if 'other' is a subset of 'self'.
        """
        return other.is_subset(self)

    def is_subset(self, other):
        """
        Returns True if 'self' is a subset of 'other'.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 0.5).is_subset(Interval(0, 1))
        True
        >>> Interval(0, 1).is_subset(Interval(0, 1, left_open=True))
        False

        """
        if isinstance(other, Set):
            return self.intersect(other) == self
        else:
            raise ValueError("Unknown argument '%s'" % other)

    def issubset(self, other):
        """
        Alias for :meth:`is_subset()`
        """
        return self.is_subset(other)

    def is_proper_subset(self, other):
        """
        Returns True if 'self' is a proper subset of 'other'.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 0.5).is_proper_subset(Interval(0, 1))
        True
        >>> Interval(0, 1).is_proper_subset(Interval(0, 1))
        False

        """
        if isinstance(other, Set):
            return self != other and self.is_subset(other)
        else:
            raise ValueError("Unknown argument '%s'" % other)

    def is_superset(self, other):
        """
        Returns True if 'self' is a superset of 'other'.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 0.5).is_superset(Interval(0, 1))
        False
        >>> Interval(0, 1).is_superset(Interval(0, 1, left_open=True))
        True

        """
        if isinstance(other, Set):
            return other.is_subset(self)
        else:
            raise ValueError("Unknown argument '%s'" % other)

    def issuperset(self, other):
        """
        Alias for :meth:`is_superset()`
        """
        return self.is_superset(other)

    def is_proper_superset(self, other):
        """
        Returns True if 'self' is a proper superset of 'other'.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1).is_proper_superset(Interval(0, 0.5))
        True
        >>> Interval(0, 1).is_proper_superset(Interval(0, 1))
        False

        """
        if isinstance(other, Set):
            return self != other and self.is_superset(other)
        else:
            raise ValueError("Unknown argument '%s'" % other)

    def _eval_powerset(self):
        raise NotImplementedError('Power set not defined for: %s' % self.func)

    def powerset(self):
        """
        Find the Power set of 'self'.

        Examples
        ========

        >>> from sympy import FiniteSet, EmptySet
        >>> A = EmptySet()
        >>> A.powerset()
        {EmptySet()}
        >>> A = FiniteSet(1, 2)
        >>> a, b, c = FiniteSet(1), FiniteSet(2), FiniteSet(1, 2)
        >>> A.powerset() == FiniteSet(a, b, c, EmptySet())
        True

        References
        ==========

        .. [1] http://en.wikipedia.org/wiki/Power_set

        """
        return self._eval_powerset()

    @property
    def measure(self):
        """
        The (Lebesgue) measure of 'self'

        Examples
        ========

        >>> from sympy import Interval, Union
        >>> Interval(0, 1).measure
        1
        >>> Union(Interval(0, 1), Interval(2, 3)).measure
        2

        """
        return self._measure

    @property
    def boundary(self):
        """
        The boundary or frontier of a set

        A point x is on the boundary of a set S if

        1.  x is in the closure of S.
            I.e. Every neighborhood of x contains a point in S.
        2.  x is not in the interior of S.
            I.e. There does not exist an open set centered on x contained
            entirely within S.

        There are the points on the outer rim of S.  If S is open then these
        points need not actually be contained within S.

        For example, the boundary of an interval is its start and end points.
        This is true regardless of whether or not the interval is open.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1).boundary
        {0, 1}
        >>> Interval(0, 1, True, False).boundary
        {0, 1}
        """
        return self._boundary

    @property
    def is_open(self):
        if not Intersection(self, self.boundary):
            return True
        # We can't confidently claim that an intersection exists
        return None

    @property
    def is_closed(self):
        return self.boundary.is_subset(self)

    @property
    def closure(self):
        return self + self.boundary

    @property
    def interior(self):
        return self - self.boundary

    @property
    def _boundary(self):
        raise NotImplementedError()

    def _eval_imageset(self, f):
        from sympy.sets.fancysets import ImageSet
        return ImageSet(f, self)

    @property
    def _measure(self):
        raise NotImplementedError("(%s)._measure" % self)

    def __add__(self, other):
        return self.union(other)

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersect(other)

    def __mul__(self, other):
        return ProductSet(self, other)

    def __xor__(self, other):
        return SymmetricDifference(self, other)

    def __pow__(self, exp):
        if not sympify(exp).is_Integer and exp >= 0:
            raise ValueError("%s: Exponent must be a positive Integer" % exp)
        return ProductSet([self]*exp)

    def __sub__(self, other):
        return Complement(self, other)

    def __contains__(self, other):
        symb = sympify(self.contains(other))
        if not (symb is S.true or symb is S.false):
            raise TypeError('contains did not evaluate to a bool: %r' % symb)
        return bool(symb)

    @property
    @deprecated(useinstead="is_subset(S.Reals)", issue=6212, deprecated_since_version="0.7.6")
    def is_real(self):
        return None


class ProductSet(Set):
    """
    Represents a Cartesian Product of Sets.

    Returns a Cartesian product given several sets as either an iterable
    or individual arguments.

    Can use '*' operator on any sets for convenient shorthand.

    Examples
    ========

    >>> from sympy import Interval, FiniteSet, ProductSet
    >>> I = Interval(0, 5); S = FiniteSet(1, 2, 3)
    >>> ProductSet(I, S)
    [0, 5] x {1, 2, 3}

    >>> (2, 2) in ProductSet(I, S)
    True

    >>> Interval(0, 1) * Interval(0, 1) # The unit square
    [0, 1] x [0, 1]

    >>> coin = FiniteSet('H', 'T')
    >>> set(coin**2)
    set([(H, H), (H, T), (T, H), (T, T)])


    Notes
    =====

    - Passes most operations down to the argument sets
    - Flattens Products of ProductSets

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Cartesian_product
    """
    is_ProductSet = True

    def __new__(cls, *sets, **assumptions):
        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_ProductSet:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            elif iterable(arg):
                return sum(map(flatten, arg), [])
            raise TypeError("Input must be Sets or iterables of Sets")
        sets = flatten(list(sets))

        if EmptySet() in sets or len(sets) == 0:
            return EmptySet()

        if len(sets) == 1:
            return sets[0]

        return Basic.__new__(cls, *sets, **assumptions)

    def _eval_Eq(self, other):
        if not other.is_ProductSet:
            return

        if len(self.args) != len(other.args):
            return false

        return And(*(Eq(x, y) for x, y in zip(self.args, other.args)))

    def _contains(self, element):
        """
        'in' operator for ProductSets

        Examples
        ========

        >>> from sympy import Interval
        >>> (2, 3) in Interval(0, 5) * Interval(0, 5)
        True

        >>> (10, 10) in Interval(0, 5) * Interval(0, 5)
        False

        Passes operation on to constituent sets
        """
        try:
            if len(element) != len(self.args):
                return false
        except TypeError:  # maybe element isn't an iterable
            return false
        return And(*
            [set.contains(item) for set, item in zip(self.sets, element)])

    def _intersect(self, other):
        """
        This function should only be used internally

        See Set._intersect for docstring
        """
        if not other.is_ProductSet:
            return None
        if len(other.args) != len(self.args):
            return S.EmptySet
        return ProductSet(a.intersect(b)
                for a, b in zip(self.sets, other.sets))

    def _union(self, other):
        if not other.is_ProductSet:
            return None
        if len(other.args) != len(self.args):
            return None
        if self.args[0] == other.args[0]:
            return self.args[0] * Union(ProductSet(self.args[1:]),
                                        ProductSet(other.args[1:]))
        if self.args[-1] == other.args[-1]:
            return Union(ProductSet(self.args[:-1]),
                         ProductSet(other.args[:-1])) * self.args[-1]
        return None

    @property
    def sets(self):
        return self.args

    @property
    def _boundary(self):
        return Union(ProductSet(b + b.boundary if i != j else b.boundary
                                for j, b in enumerate(self.sets))
                                for i, a in enumerate(self.sets))


    @property
    @deprecated(useinstead="is_subset(S.Reals)", issue=6212, deprecated_since_version="0.7.6")
    def is_real(self):
        return all(set.is_real for set in self.sets)

    @property
    def is_iterable(self):
        return all(set.is_iterable for set in self.sets)

    def __iter__(self):
        if self.is_iterable:
            return product(*self.sets)
        else:
            raise TypeError("Not all constituent sets are iterable")

    @property
    def _measure(self):
        measure = 1
        for set in self.sets:
            measure *= set.measure
        return measure

    def __len__(self):
        return Mul(*[len(s) for s in self.args])


class Interval(Set, EvalfMixin):
    """
    Represents a real interval as a Set.

    Usage:
        Returns an interval with end points "start" and "end".

        For left_open=True (default left_open is False) the interval
        will be open on the left. Similarly, for right_open=True the interval
        will be open on the right.

    Examples
    ========

    >>> from sympy import Symbol, Interval
    >>> Interval(0, 1)
    [0, 1]
    >>> Interval(0, 1, False, True)
    [0, 1)
    >>> Interval.Ropen(0, 1)
    [0, 1)
    >>> Interval.Lopen(0, 1)
    (0, 1]
    >>> Interval.open(0, 1)
    (0, 1)

    >>> a = Symbol('a', real=True)
    >>> Interval(0, a)
    [0, a]

    Notes
    =====
    - Only real end points are supported
    - Interval(a, b) with a > b will return the empty set
    - Use the evalf() method to turn an Interval into an mpmath
      'mpi' interval instance

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Interval_%28mathematics%29
    """
    is_Interval = True

    @property
    @deprecated(useinstead="is_subset(S.Reals)", issue=6212, deprecated_since_version="0.7.6")
    def is_real(self):
        return True

    def __new__(cls, start, end, left_open=False, right_open=False):

        start = _sympify(start)
        end = _sympify(end)
        left_open = _sympify(left_open)
        right_open = _sympify(right_open)

        if not all(isinstance(a, (type(true), type(false)))
            for a in [left_open, right_open]):
            raise NotImplementedError(
                "left_open and right_open can have only true/false values, "
                "got %s and %s" % (left_open, right_open))

        inftys = [S.Infinity, S.NegativeInfinity]
        # Only allow real intervals (use symbols with 'is_real=True').
        if not all(i.is_real is not False or i in inftys for i in (start, end)):
            raise ValueError("Non-real intervals are not supported")

        # evaluate if possible
        if (end < start) == True:
            return S.EmptySet
        elif (end - start).is_negative:
            return S.EmptySet

        if end == start and (left_open or right_open):
            return S.EmptySet
        if end == start and not (left_open or right_open):
            return FiniteSet(end)

        # Make sure infinite interval end points are open.
        if start == S.NegativeInfinity:
            left_open = true
        if end == S.Infinity:
            right_open = true

        return Basic.__new__(cls, start, end, left_open, right_open)

    @property
    def start(self):
        """
        The left end point of 'self'.

        This property takes the same value as the 'inf' property.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1).start
        0

        """
        return self._args[0]

    _inf = left = start

    @classmethod
    def open(cls, a, b):
        """Return an interval including neither boundary."""
        return cls(a, b, True, True)

    @classmethod
    def Lopen(cls, a, b):
        """Return an interval not including the left boundary."""
        return cls(a, b, True, False)

    @classmethod
    def Ropen(cls, a, b):
        """Return an interval not including the right boundary."""
        return cls(a, b, False, True)

    @property
    def end(self):
        """
        The right end point of 'self'.

        This property takes the same value as the 'sup' property.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1).end
        1

        """
        return self._args[1]

    _sup = right = end

    @property
    def left_open(self):
        """
        True if 'self' is left-open.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1, left_open=True).left_open
        True
        >>> Interval(0, 1, left_open=False).left_open
        False

        """
        return self._args[2]

    @property
    def right_open(self):
        """
        True if 'self' is right-open.

        Examples
        ========

        >>> from sympy import Interval
        >>> Interval(0, 1, right_open=True).right_open
        True
        >>> Interval(0, 1, right_open=False).right_open
        False

        """
        return self._args[3]

    def _intersect(self, other):
        """
        This function should only be used internally

        See Set._intersect for docstring
        """
        # We only know how to intersect with other intervals
        if not other.is_Interval:
            return None

        # handle (-oo, oo)
        infty = S.NegativeInfinity, S.Infinity
        if self == Interval(*infty):
            l, r = self.left, self.right
            if l.is_real or l in infty or r.is_real or r in infty:
                return other

        # We can't intersect [0,3] with [x,6] -- we don't know if x>0 or x<0
        if not self._is_comparable(other):
            return None

        empty = False

        if self.start <= other.end and other.start <= self.end:
            # Get topology right.
            if self.start < other.start:
                start = other.start
                left_open = other.left_open
            elif self.start > other.start:
                start = self.start
                left_open = self.left_open
            else:
                start = self.start
                left_open = self.left_open or other.left_open

            if self.end < other.end:
                end = self.end
                right_open = self.right_open
            elif self.end > other.end:
                end = other.end
                right_open = other.right_open
            else:
                end = self.end
                right_open = self.right_open or other.right_open

            if end - start == 0 and (left_open or right_open):
                empty = True
        else:
            empty = True

        if empty:
            return S.EmptySet

        return Interval(start, end, left_open, right_open)


    def _complement(self, other):
        if other == S.Reals:
            a = Interval(S.NegativeInfinity, self.start,
                         True, not self.left_open)
            b = Interval(self.end, S.Infinity, not self.right_open, True)
            return Union(a, b)

        if isinstance(other, FiniteSet):
            nums = [m for m in other.args if m.is_number]
            if nums == []:
                return None

        return Set._complement(self, other)


    def _union(self, other):
        """
        This function should only be used internally

        See Set._union for docstring
        """
        if other.is_Interval and self._is_comparable(other):
            from sympy.functions.elementary.miscellaneous import Min, Max
            # Non-overlapping intervals
            end = Min(self.end, other.end)
            start = Max(self.start, other.start)
            if (end < start or
               (end == start and (end not in self and end not in other))):
                return None
            else:
                start = Min(self.start, other.start)
                end = Max(self.end, other.end)

                left_open = ((self.start != start or self.left_open) and
                             (other.start != start or other.left_open))
                right_open = ((self.end != end or self.right_open) and
                              (other.end != end or other.right_open))

                return Interval(start, end, left_open, right_open)

        # If I have open end points and these endpoints are contained in other
        if ((self.left_open and sympify(other.contains(self.start)) is S.true) or
                (self.right_open and sympify(other.contains(self.end)) is S.true)):
            # Fill in my end points and return
            open_left = self.left_open and self.start not in other
            open_right = self.right_open and self.end not in other
            new_self = Interval(self.start, self.end, open_left, open_right)
            return set((new_self, other))

        return None

    @property
    def _boundary(self):
        return FiniteSet(self.start, self.end)

    def _contains(self, other):
        if other.is_real is False or other is S.NegativeInfinity or other is S.Infinity:
            return false


        if self.start is S.NegativeInfinity and self.end is S.Infinity:
            if not other.is_real is None:
                return other.is_real

        if self.left_open:
            expr = other > self.start
        else:
            expr = other >= self.start

        if self.right_open:
            expr = And(expr, other < self.end)
        else:
            expr = And(expr, other <= self.end)

        return _sympify(expr)

    def _eval_imageset(self, f):
        from sympy.functions.elementary.miscellaneous import Min, Max
        from sympy.solvers.solveset import solveset
        from sympy.core.function import diff, Lambda
        from sympy.series import limit
        from sympy.calculus.singularities import singularities
        # TODO: handle functions with infinitely many solutions (eg, sin, tan)
        # TODO: handle multivariate functions

        expr = f.expr
        if len(expr.free_symbols) > 1 or len(f.variables) != 1:
            return
        var = f.variables[0]

        if expr.is_Piecewise:
            result = S.EmptySet
            domain_set = self
            for (p_expr, p_cond) in expr.args:
                if p_cond is true:
                    intrvl = domain_set
                else:
                    intrvl = p_cond.as_set()
                    intrvl = Intersection(domain_set, intrvl)

                if p_expr.is_Number:
                    image = FiniteSet(p_expr)
                else:
                    image = imageset(Lambda(var, p_expr), intrvl)
                result = Union(result, image)

                # remove the part which has been `imaged`
                domain_set = Complement(domain_set, intrvl)
                if domain_set.is_EmptySet:
                    break
            return result

        if not self.start.is_comparable or not self.end.is_comparable:
            return

        try:
            sing = [x for x in singularities(expr, var)
                if x.is_real and x in self]
        except NotImplementedError:
            return

        if self.left_open:
            _start = limit(expr, var, self.start, dir="+")
        elif self.start not in sing:
            _start = f(self.start)
        if self.right_open:
            _end = limit(expr, var, self.end, dir="-")
        elif self.end not in sing:
            _end = f(self.end)

        if len(sing) == 0:
            solns = list(solveset(diff(expr, var), var))

            extr = [_start, _end] + [f(x) for x in solns
                                     if x.is_real and x in self]
            start, end = Min(*extr), Max(*extr)

            left_open, right_open = False, False
            if _start <= _end:
                # the minimum or maximum value can occur simultaneously
                # on both the edge of the interval and in some interior
                # point
                if start == _start and start not in solns:
                    left_open = self.left_open
                if end == _end and end not in solns:
                    right_open = self.right_open
            else:
                if start == _end and start not in solns:
                    left_open = self.right_open
                if end == _start and end not in solns:
                    right_open = self.left_open

            return Interval(start, end, left_open, right_open)
        else:
            return imageset(f, Interval(self.start, sing[0],
                                        self.left_open, True)) + \
                Union(*[imageset(f, Interval(sing[i], sing[i + 1], True, True))
                        for i in range(0, len(sing) - 1)]) + \
                imageset(f, Interval(sing[-1], self.end, True, self.right_open))

    @property
    def _measure(self):
        return self.end - self.start

    def to_mpi(self, prec=53):
        return mpi(mpf(self.start._eval_evalf(prec)),
            mpf(self.end._eval_evalf(prec)))

    def _eval_evalf(self, prec):
        return Interval(self.left._eval_evalf(prec),
            self.right._eval_evalf(prec),
                        left_open=self.left_open, right_open=self.right_open)

    def _is_comparable(self, other):
        is_comparable = self.start.is_comparable
        is_comparable &= self.end.is_comparable
        is_comparable &= other.start.is_comparable
        is_comparable &= other.end.is_comparable

        return is_comparable

    @property
    def is_left_unbounded(self):
        """Return ``True`` if the left endpoint is negative infinity. """
        return self.left is S.NegativeInfinity or self.left == Float("-inf")

    @property
    def is_right_unbounded(self):
        """Return ``True`` if the right endpoint is positive infinity. """
        return self.right is S.Infinity or self.right == Float("+inf")

    def as_relational(self, x):
        """Rewrite an interval in terms of inequalities and logic operators."""
        x = sympify(x)
        if self.right_open:
            right = x < self.end
        else:
            right = x <= self.end
        if self.left_open:
            left = self.start < x
        else:
            left = self.start <= x
        return And(left, right)

    def _eval_Eq(self, other):
        if not other.is_Interval:
            if (other.is_Union or other.is_Complement or
                other.is_Intersection or other.is_ProductSet):
                return

            return false

        return And(Eq(self.left, other.left),
                   Eq(self.right, other.right),
                   self.left_open == other.left_open,
                   self.right_open == other.right_open)


class Union(Set, EvalfMixin):
    """
    Represents a union of sets as a :class:`Set`.

    Examples
    ========

    >>> from sympy import Union, Interval
    >>> Union(Interval(1, 2), Interval(3, 4))
    [1, 2] U [3, 4]

    The Union constructor will always try to merge overlapping intervals,
    if possible. For example:

    >>> Union(Interval(1, 2), Interval(2, 3))
    [1, 3]

    See Also
    ========

    Intersection

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Union_%28set_theory%29
    """
    is_Union = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs to merge intersections and iterables
        args = list(args)

        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_Union:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):  # and not isinstance(arg, Set) (implicit)
                return sum(map(flatten, arg), [])
            raise TypeError("Input must be Sets or iterables of Sets")
        args = flatten(args)

        # Union of no sets is EmptySet
        if len(args) == 0:
            return S.EmptySet

        # Reduce sets using known rules
        if evaluate:
            return Union.reduce(args)

        args = list(ordered(args, Set._infimum_key))

        return Basic.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a :class:`Union` using known rules

        We first start with global rules like
        'Merge all FiniteSets'

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent
        """

        # ===== Global Rules =====
        # Merge all finite sets
        finite_sets = [x for x in args if x.is_FiniteSet]
        if len(finite_sets) > 1:
            a = (x for set in finite_sets for x in set)
            finite_set = FiniteSet(*a)
            args = [finite_set] + [x for x in args if not x.is_FiniteSet]

        # ===== Pair-wise Rules =====
        # Here we depend on rules built into the constituent sets
        args = set(args)
        new_args = True
        while(new_args):
            for s in args:
                new_args = False
                for t in args - set((s,)):
                    new_set = s._union(t)
                    # This returns None if s does not know how to intersect
                    # with t. Returns the newly intersected set otherwise
                    if new_set is not None:
                        if not isinstance(new_set, set):
                            new_set = set((new_set, ))
                        new_args = (args - set((s, t))).union(new_set)
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return Union(args, evaluate=False)

    def _complement(self, universe):
        # DeMorgan's Law
        return Intersection(s.complement(universe) for s in self.args)

    @property
    def _inf(self):
        # We use Min so that sup is meaningful in combination with symbolic
        # interval end points.
        from sympy.functions.elementary.miscellaneous import Min
        return Min(*[set.inf for set in self.args])

    @property
    def _sup(self):
        # We use Max so that sup is meaningful in combination with symbolic
        # end points.
        from sympy.functions.elementary.miscellaneous import Max
        return Max(*[set.sup for set in self.args])

    def _contains(self, other):
        return Or(*[set.contains(other) for set in self.args])

    @property
    def _measure(self):
        # Measure of a union is the sum of the measures of the sets minus
        # the sum of their pairwise intersections plus the sum of their
        # triple-wise intersections minus ... etc...

        # Sets is a collection of intersections and a set of elementary
        # sets which made up those intersections (called "sos" for set of sets)
        # An example element might of this list might be:
        #    ( {A,B,C}, A.intersect(B).intersect(C) )

        # Start with just elementary sets (  ({A}, A), ({B}, B), ... )
        # Then get and subtract (  ({A,B}, (A int B), ... ) while non-zero
        sets = [(FiniteSet(s), s) for s in self.args]
        measure = 0
        parity = 1
        while sets:
            # Add up the measure of these sets and add or subtract it to total
            measure += parity * sum(inter.measure for sos, inter in sets)

            # For each intersection in sets, compute the intersection with every
            # other set not already part of the intersection.
            sets = ((sos + FiniteSet(newset), newset.intersect(intersection))
                    for sos, intersection in sets for newset in self.args
                    if newset not in sos)

            # Clear out sets with no measure
            sets = [(sos, inter) for sos, inter in sets if inter.measure != 0]

            # Clear out duplicates
            sos_list = []
            sets_list = []
            for set in sets:
                if set[0] in sos_list:
                    continue
                else:
                    sos_list.append(set[0])
                    sets_list.append(set)
            sets = sets_list

            # Flip Parity - next time subtract/add if we added/subtracted here
            parity *= -1
        return measure

    @property
    def _boundary(self):
        def boundary_of_set(i):
            """ The boundary of set i minus interior of all other sets """
            b = self.args[i].boundary
            for j, a in enumerate(self.args):
                if j != i:
                    b = b - a.interior
            return b
        return Union(map(boundary_of_set, range(len(self.args))))

    def _eval_imageset(self, f):
        return Union(imageset(f, arg) for arg in self.args)

    def as_relational(self, symbol):
        """Rewrite a Union in terms of equalities and logic operators. """
        return Or(*[set.as_relational(symbol) for set in self.args])

    @property
    def is_iterable(self):
        return all(arg.is_iterable for arg in self.args)

    def _eval_evalf(self, prec):
        try:
            return Union(set._eval_evalf(prec) for set in self.args)
        except Exception:
            raise TypeError("Not all sets are evalf-able")

    def __iter__(self):
        import itertools

        # roundrobin recipe taken from itertools documentation:
        # https://docs.python.org/2/library/itertools.html#recipes
        def roundrobin(*iterables):
            "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
            # Recipe credited to George Sakkis
            pending = len(iterables)
            nexts = itertools.cycle(iter(it).next for it in iterables)
            while pending:
                try:
                    for next in nexts:
                        yield next()
                except StopIteration:
                    pending -= 1
                    nexts = itertools.cycle(itertools.islice(nexts, pending))

        if all(set.is_iterable for set in self.args):
            return roundrobin(*(iter(arg) for arg in self.args))
        else:
            raise TypeError("Not all constituent sets are iterable")

    @property
    @deprecated(useinstead="is_subset(S.Reals)", issue=6212, deprecated_since_version="0.7.6")
    def is_real(self):
        return all(set.is_real for set in self.args)


class Intersection(Set):
    """
    Represents an intersection of sets as a :class:`Set`.

    Examples
    ========

    >>> from sympy import Intersection, Interval
    >>> Intersection(Interval(1, 3), Interval(2, 4))
    [2, 3]

    We often use the .intersect method

    >>> Interval(1,3).intersect(Interval(2,4))
    [2, 3]

    See Also
    ========

    Union

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Intersection_%28set_theory%29
    """
    is_Intersection = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs to merge intersections and iterables
        args = list(args)

        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_Intersection:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):  # and not isinstance(arg, Set) (implicit)
                return sum(map(flatten, arg), [])
            raise TypeError("Input must be Sets or iterables of Sets")
        args = flatten(args)

        if len(args) == 0:
            return S.EmptySet

        # args can't be ordered for Partition see issue #9608
        if 'Partition' not in [type(a).__name__ for a in args]:
            args = list(ordered(args, Set._infimum_key))

        # Reduce sets using known rules
        if evaluate:
            return Intersection.reduce(args)

        return Basic.__new__(cls, *args)

    @property
    def is_iterable(self):
        return any(arg.is_iterable for arg in self.args)

    @property
    def _inf(self):
        raise NotImplementedError()

    @property
    def _sup(self):
        raise NotImplementedError()

    def _eval_imageset(self, f):
        return Intersection(imageset(f, arg) for arg in self.args)

    def _contains(self, other):
        return And(*[set.contains(other) for set in self.args])

    def __iter__(self):
        no_iter = True
        for s in self.args:
            if s.is_iterable:
                no_iter = False
                other_sets = set(self.args) - set((s,))
                other = Intersection(other_sets, evaluate=False)
                for x in s:
                    c = sympify(other.contains(x))
                    if c is S.true:
                        yield x
                    elif c is S.false:
                        pass
                    else:
                        yield c

        if no_iter:
            raise ValueError("None of the constituent sets are iterable")

    @staticmethod
    def _handle_finite_sets(args):
        from sympy.core.logic import fuzzy_and, fuzzy_bool
        from sympy.core.compatibility import zip_longest

        new_args = []
        fs_args = []
        for s in args:
            if s.is_FiniteSet:
                fs_args.append(s)
            else:
                new_args.append(s)
        if not fs_args:
            return
        s = fs_args[0]
        fs_args = fs_args[1:]
        res = []
        unk = []
        for x in s:
            c = fuzzy_and(fuzzy_bool(o.contains(x))
                for o in fs_args + new_args)
            if c:
                res.append(x)
            elif c is None:
                unk.append(x)
            else:
                pass  # drop arg
        res = FiniteSet(
            *res, evaluate=False) if res else S.EmptySet
        if unk:
            symbolic_s_list = [x for x in s if x.has(Symbol)]
            non_symbolic_s = s - FiniteSet(
                *symbolic_s_list, evaluate=False)
            while fs_args:
                v = fs_args.pop()
                if all(i == j for i, j in zip_longest(
                        symbolic_s_list,
                        (x for x in v if x.has(Symbol)))):
                    # all the symbolic elements of `v` are the same
                    # as in `s` so remove the non-symbol containing
                    # expressions from `unk`, since they cannot be
                    # contained
                    for x in non_symbolic_s:
                        if x in unk:
                            unk.remove(x)
                else:
                    # if only a subset of elements in `s` are
                    # contained in `v` then remove them from `v`
                    # and add this as a new arg
                    contained = [x for x in symbolic_s_list
                        if sympify(v.contains(x)) is S.true]
                    if contained != symbolic_s_list:
                        new_args.append(
                            v - FiniteSet(
                            *contained, evaluate=False))
                    else:
                        pass  # for coverage

            other_sets = Intersection(*new_args)
            if not other_sets:
                return S.EmptySet  # b/c we use evaluate=False below
            res += Intersection(
                FiniteSet(*unk),
                other_sets, evaluate=False)
        return res

    @staticmethod
    def reduce(args):
        """
        Simplify an intersection using known rules

        We first start with global rules like
        'if any empty sets return empty set' and 'distribute any unions'

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent
        """

        # ===== Global Rules =====
        # If any EmptySets return EmptySet
        if any(s.is_EmptySet for s in args):
            return S.EmptySet

        # Handle Finite sets
        rv = Intersection._handle_finite_sets(args)
        if rv is not None:
            return rv

        # If any of the sets are unions, return a Union of Intersections
        for s in args:
            if s.is_Union:
                other_sets = set(args) - set((s,))
                if len(other_sets) > 0:
                    other = Intersection(other_sets)
                    return Union(Intersection(arg, other) for arg in s.args)
                else:
                    return Union(arg for arg in s.args)

        for s in args:
            if s.is_Complement:
                args.remove(s)
                other_sets = args + [s.args[0]]
                return Complement(Intersection(*other_sets), s.args[1])

        # At this stage we are guaranteed not to have any
        # EmptySets, FiniteSets, or Unions in the intersection

        # ===== Pair-wise Rules =====
        # Here we depend on rules built into the constituent sets
        args = set(args)
        new_args = True
        while(new_args):
            for s in args:
                new_args = False
                for t in args - set((s,)):
                    new_set = s._intersect(t)
                    # This returns None if s does not know how to intersect
                    # with t. Returns the newly intersected set otherwise
                    if new_set is not None:
                        new_args = (args - set((s, t))).union(set((new_set, )))
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return Intersection(args, evaluate=False)

    def as_relational(self, symbol):
        """Rewrite an Intersection in terms of equalities and logic operators"""
        return And(*[set.as_relational(symbol) for set in self.args])


class Complement(Set, EvalfMixin):
    """Represents the set difference or relative complement of a set with
    another set.

    `A - B = \{x \in A| x \\notin B\}`


    Examples
    ========

    >>> from sympy import Complement, FiniteSet
    >>> Complement(FiniteSet(0, 1, 2), FiniteSet(1))
    {0, 2}

    See Also
    =========

    Intersection, Union

    References
    ==========

    .. [1] http://mathworld.wolfram.com/ComplementSet.html
    """

    is_Complement = True

    def __new__(cls, a, b, evaluate=True):
        if evaluate:
            return Complement.reduce(a, b)

        return Basic.__new__(cls, a, b)

    @staticmethod
    def reduce(A, B):
        """
        Simplify a :class:`Complement`.

        """
        if B == S.UniversalSet:
            return EmptySet()

        if isinstance(B, Union):
            return Intersection(s.complement(A) for s in B.args)

        result = B._complement(A)
        if result != None:
            return result
        else:
            return Complement(A, B, evaluate=False)

    def _contains(self, other):
        A = self.args[0]
        B = self.args[1]
        return And(A.contains(other), Not(B.contains(other)))


class EmptySet(with_metaclass(Singleton, Set)):
    """
    Represents the empty set. The empty set is available as a singleton
    as S.EmptySet.

    Examples
    ========

    >>> from sympy import S, Interval
    >>> S.EmptySet
    EmptySet()

    >>> Interval(1, 2).intersect(S.EmptySet)
    EmptySet()

    See Also
    ========

    UniversalSet

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Empty_set
    """
    is_EmptySet = True
    is_FiniteSet = True

    def _intersect(self, other):
        return S.EmptySet

    @property
    def _measure(self):
        return 0

    def _contains(self, other):
        return false

    def as_relational(self, symbol):
        return False

    def __len__(self):
        return 0

    def _union(self, other):
        return other

    def __iter__(self):
        return iter([])

    def _eval_imageset(self, f):
        return self

    def _eval_powerset(self):
        return FiniteSet(self)

    @property
    def _boundary(self):
        return self

    def _complement(self, other):
        return other

    def _symmetric_difference(self, other):
        return other


class UniversalSet(with_metaclass(Singleton, Set)):
    """
    Represents the set of all things.
    The universal set is available as a singleton as S.UniversalSet

    Examples
    ========

    >>> from sympy import S, Interval
    >>> S.UniversalSet
    UniversalSet()

    >>> Interval(1, 2).intersect(S.UniversalSet)
    [1, 2]

    See Also
    ========

    EmptySet

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Universal_set
    """

    is_UniversalSet = True

    def _intersect(self, other):
        return other

    def _complement(self, other):
        return S.EmptySet

    def _symmetric_difference(self, other):
        return other

    @property
    def _measure(self):
        return S.Infinity

    def _contains(self, other):
        return true

    def as_relational(self, symbol):
        return True

    def _union(self, other):
        return self

    @property
    def _boundary(self):
        return EmptySet()


class FiniteSet(Set, EvalfMixin):
    """
    Represents a finite set of discrete numbers

    Examples
    ========

    >>> from sympy import FiniteSet
    >>> FiniteSet(1, 2, 3, 4)
    {1, 2, 3, 4}
    >>> 3 in FiniteSet(1, 2, 3, 4)
    True

    >>> members = [1, 2, 3, 4]
    >>> FiniteSet(*members)
    {1, 2, 3, 4}

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Finite_set
    """
    is_FiniteSet = True
    is_iterable = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])
        if evaluate:
            args = list(map(sympify, args))

            if len(args) == 0:
                return EmptySet()
        else:
            args = list(map(sympify, args))

        args = list(ordered(frozenset(tuple(args)), Set._infimum_key))
        obj = Basic.__new__(cls, *args)
        obj._elements = frozenset(args)
        return obj

    def _eval_Eq(self, other):
        if not other.is_FiniteSet:
            if (other.is_Union or other.is_Complement or
                other.is_Intersection or other.is_ProductSet):
                return

            return false

        if len(self) != len(other):
            return false

        return And(*(Eq(x, y) for x, y in zip(self.args, other.args)))

    def __iter__(self):
        return iter(self.args)

    def _intersect(self, other):
        """
        This function should only be used internally

        See Set._intersect for docstring
        """
        if isinstance(other, self.__class__):
            return self.__class__(*(self._elements & other._elements))
        return self.__class__(*[el for el in self if el in other])

    def _complement(self, other):
        if isinstance(other, Interval):
            nums = sorted(m for m in self.args if m.is_number)
            if other == S.Reals and nums != []:
                syms = [m for m in self.args if m.is_Symbol]
                # Reals cannot contain elements other than numbers and symbols.

                intervals = []  # Build up a list of intervals between the elements
                intervals += [Interval(S.NegativeInfinity, nums[0], True, True)]
                for a, b in zip(nums[:-1], nums[1:]):
                    intervals.append(Interval(a, b, True, True))  # both open
                intervals.append(Interval(nums[-1], S.Infinity, True, True))

                if syms != []:
                    return Complement(Union(intervals, evaluate=False),
                            FiniteSet(*syms), evaluate=False)
                else:
                    return Union(intervals, evaluate=False)
            elif nums == []:
                return None

        elif isinstance(other, FiniteSet):
            unk = []
            for i in self:
                c = sympify(other.contains(i))
                if c is not S.true and c is not S.false:
                    unk.append(i)
            unk = FiniteSet(*unk)
            if unk == self:
                return
            not_true = []
            for i in other:
                c = sympify(self.contains(i))
                if c is not S.true:
                    not_true.append(i)
            return Complement(FiniteSet(*not_true), unk)

        return Set._complement(self, other)


    def _union(self, other):
        """
        This function should only be used internally

        See Set._union for docstring
        """
        if other.is_FiniteSet:
            return FiniteSet(*(self._elements | other._elements))

        # If other set contains one of my elements, remove it from myself
        if any(sympify(other.contains(x)) is S.true for x in self):
            return set((
                FiniteSet(*[x for x in self
                    if other.contains(x) != True]), other))

        return None


    def _contains(self, other):
        """
        Tests whether an element, other, is in the set.

        Relies on Python's set class. This tests for object equality
        All inputs are sympified

        Examples
        ========

        >>> from sympy import FiniteSet
        >>> 1 in FiniteSet(1, 2)
        True
        >>> 5 in FiniteSet(1, 2)
        False

        """
        r = false
        for e in self._elements:
            t = Eq(e, other, evaluate=True)
            if isinstance(t, Eq):
                t = t.simplify()
            if t == true:
                return t
            elif t != false:
                r = None
        return r

    def _eval_imageset(self, f):
        return FiniteSet(*map(f, self))

    @property
    def _boundary(self):
        return self

    @property
    def _inf(self):
        from sympy.functions.elementary.miscellaneous import Min
        return Min(*self)

    @property
    def _sup(self):
        from sympy.functions.elementary.miscellaneous import Max
        return Max(*self)

    @property
    def measure(self):
        return 0

    def __len__(self):
        return len(self.args)

    def as_relational(self, symbol):
        """Rewrite a FiniteSet in terms of equalities and logic operators. """
        from sympy.core.relational import Eq
        return Or(*[Eq(symbol, elem) for elem in self])

    @property
    @deprecated(useinstead="is_subset(S.Reals)", issue=6212, deprecated_since_version="0.7.6")
    def is_real(self):
        return all(el.is_real for el in self)

    def compare(self, other):
        return (hash(self) - hash(other))

    def _eval_evalf(self, prec):
        return FiniteSet(*[elem._eval_evalf(prec) for elem in self])

    def _hashable_content(self):
        return (self._elements,)

    @property
    def _sorted_args(self):
        return tuple(ordered(self.args, Set._infimum_key))

    def _eval_powerset(self):
        return self.func(*[self.func(*s) for s in subsets(self.args)])

    def __ge__(self, other):
        if not isinstance(other, Set):
            raise TypeError("Invalid comparison of set with %s" % func_name(other))
        return other.is_subset(self)

    def __gt__(self, other):
        if not isinstance(other, Set):
            raise TypeError("Invalid comparison of set with %s" % func_name(other))
        return self.is_proper_superset(other)

    def __le__(self, other):
        if not isinstance(other, Set):
            raise TypeError("Invalid comparison of set with %s" % func_name(other))
        return self.is_subset(other)

    def __lt__(self, other):
        if not isinstance(other, Set):
            raise TypeError("Invalid comparison of set with %s" % func_name(other))
        return self.is_proper_subset(other)


class SymmetricDifference(Set):
    """Represents the set of elements which are in either of the
    sets and not in their intersection.

    Examples
    ========

    >>> from sympy import SymmetricDifference, FiniteSet
    >>> SymmetricDifference(FiniteSet(1, 2, 3), FiniteSet(3, 4, 5))
    {1, 2, 4, 5}

    See Also
    ========

    Complement, Union

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Symmetric_difference
    """

    is_SymmetricDifference = True

    def __new__(cls, a, b, evaluate=True):
        if evaluate:
            return SymmetricDifference.reduce(a, b)

        return Basic.__new__(cls, a, b)

    @staticmethod
    def reduce(A, B):
        result = B._symmetric_difference(A)
        if result is not None:
            return result
        else:
            return SymmetricDifference(A, B, evaluate=False)


def imageset(*args):
    r"""
    Image of set under transformation ``f``.

    If this function can't compute the image, it returns an
    unevaluated ImageSet object.

    .. math::
        { f(x) | x \in self }

    Examples
    ========

    >>> from sympy import Interval, Symbol, imageset, sin, Lambda
    >>> x = Symbol('x')

    >>> imageset(x, 2*x, Interval(0, 2))
    [0, 4]

    >>> imageset(lambda x: 2*x, Interval(0, 2))
    [0, 4]

    >>> imageset(Lambda(x, sin(x)), Interval(-2, 1))
    ImageSet(Lambda(x, sin(x)), [-2, 1])

    See Also
    ========

    sympy.sets.fancysets.ImageSet

    """
    from sympy.core import Dummy, Lambda
    from sympy.sets.fancysets import ImageSet
    if len(args) == 3:
        f = Lambda(*args[:2])
    else:
        # var and expr are being defined this way to
        # support Python lambda and not just sympy Lambda
        f = args[0]
        if not isinstance(f, Lambda):
            var = Dummy()
            expr = args[0](var)
            f = Lambda(var, expr)
    set = args[-1]

    r = set._eval_imageset(f)
    if isinstance(r, ImageSet):
        f, set = r.args

    if f.variables[0] == f.expr:
        return set

    if isinstance(set, ImageSet):
        if len(set.lamda.variables) == 1 and len(f.variables) == 1:
            return imageset(Lambda(set.lamda.variables[0],
                                   f.expr.subs(f.variables[0], set.lamda.expr)),
                            set.base_set)

    if r is not None:
        return r

    return ImageSet(f, set)
