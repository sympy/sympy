from basic import Basic, SingletonMeta, S
from sympify import _sympify
from sympy.mpmath import mpi

class Set(Basic):
    """
    Represents any kind of set.

    Real intervals are represented by the Interval class and unions of sets
    by the Union class. The empty set is represented by the EmptySet class
    and available as a singleton as S.EmptySet.
    """

    def union(self, other):
        """
        Returns the union of 'self' and 'other'. As a shortcut it is possible
        to use the '+' operator:

        >>> from sympy import *

        >>> Interval(0, 1).union(Interval(2, 3))
        Union([0, 1], [2, 3])
        >>> Interval(0, 1) + Interval(2, 3)
        Union([0, 1], [2, 3])

        Similarly it is possible to use the '-' operator for set
        differences:

        >>> Interval(0, 2) - Interval(0, 1)
        (1, 2]

        """
        return Union(self, other)

    def intersect(self, other):
        """
        Returns the intersection of 'self' and 'other'.

        >>> from sympy import *

        >>> Interval(1, 3).intersect(Interval(1, 2))
        [1, 2]

        """
        return self._intersect(other)

    def _intersect(self, other):
        raise NotImplementedError("(%s)._intersect(%s)" % (self, other))

    @property
    def complement(self):
        """
        The complement of 'self'.

        As a shortcut it is possible to use the '~' or '-' operators:

        >>> from sympy import *

        >>> Interval(0, 1).complement
        Union((-oo, 0), (1, oo))
        >>> ~Interval(0, 1)
        Union((-oo, 0), (1, oo))
        >>> -Interval(0, 1)
        Union((-oo, 0), (1, oo))

        """
        return self._complement

    @property
    def _complement(self):
        raise NotImplementedError("(%s)._complement" % self)

    @property
    def inf(self):
        """
        The infimum of 'self'.

        >>> from sympy import *

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
        """ The supremum of 'self'.

        >>> from sympy import *

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

        >>> from sympy import *

        >>> Interval(0, 1).contains(0.5)
        True
        >>> 0.5 in Interval(0, 1)
        True

        """
        return self._contains(other)

    def _contains(self, other):
        raise NotImplementedError("(%s)._contains(%s)" % (self, other))

    def subset(self, other):
        """
        Returns True if 'other' is a subset of 'self'.

        >>> from sympy import *

        >>> Interval(0, 1).contains(0)
        True
        >>> Interval(0, 1, left_open=True).contains(0)
        False

        """
        if isinstance(other, Set):
            return self.intersect(other) == other
        else:
            raise ValueError, "Unknown argument '%s'" % other

    @property
    def measure(self):
        """
        The (Lebesgue) measure of 'self'.

        >>> from sympy import *

        >>> Interval(0, 1).measure
        1
        >>> Union(Interval(0, 1), Interval(2, 3)).measure
        2

        """
        return self._measure

    @property
    def _measure(self):
        raise NotImplementedError("(%s)._measure" % self)

    def __add__(self, other):
        return self.union(other)

    def __sub__(self, other):
        return self.intersect(other.complement)

    def __neg__(self):
        return self.complement

    def __invert__(self):
        return self.complement

    def __contains__(self, other):
        return self.contains(other)

    def _eval_subs(self, old, new):
        if self == old: return new
        new_args = []
        for arg in self.args:
            if arg == old:
                new_args.append(new)
            elif isinstance(arg, Basic):
                new_args.append(arg._eval_subs(old, new))
            else:
                new_args.append(arg)
        return self.__class__(*new_args)

class Interval(Set):
    """
    Represents a real interval as a Set.

    Usage:
        Returns an interval with end points "start" and "end".

        For left_open=True (default left_open is False) the interval
        will be open on the left. Similarly, for right_open=True the interval
        will be open on the right.

    Examples:
        >>> from sympy import Symbol, Interval, sets

        >>> Interval(0, 1)
        [0, 1]
        >>> Interval(0, 1, False, True)
        [0, 1)

        >>> a = Symbol('a', real=True)
        >>> Interval(0, a)
        [0, a]

    Notes:
        - Only real end points are supported
        - Interval(a, b) with a > b will return the empty set
        - Use the evalf() method to turn an Interval into an mpmath
          'mpi' interval instance
    """

    def __new__(cls, start, end,
                left_open=False, right_open=False, **assumptions):

        start = _sympify(start)
        end = _sympify(end)

        # Only allow real intervals (use symbols with 'is_real=True').
        if not start.is_real or not end.is_real:
            raise ValueError, "Only real intervals are supported"

        # Make sure that the created interval will be valid.
        if end.is_comparable and start.is_comparable:
            if end < start:
                return S.EmptySet

        if end == start and (left_open or right_open):
            return S.EmptySet

        # Make sure infinite interval end points are open.
        if start == S.NegativeInfinity:
            left_open = True
        if end == S.Infinity:
            right_open = True

        return Basic.__new__(cls, start, end,
                             left_open, right_open, **assumptions)

    @property
    def start(self):
        """
        The left end point of 'self'. This property takes the same value as the
        'inf' property.

        >>> from sympy import *

        >>> Interval(0, 1).start
        0

        """
        return self._args[0]

    _inf = start

    @property
    def end(self):
        """
        The right end point of 'self'. This property takes the same value as the
        'sup' property.

        >>> from sympy import *

        >>> Interval(0, 1).end
        1

        """
        return self._args[1]

    _sup = end

    @property
    def left_open(self):
        """
        True if 'self' is left-open.

        >>> from sympy import *

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

        >>> from sympy import *

        >>> Interval(0, 1, right_open=True).right_open
        True
        >>> Interval(0, 1, right_open=False).right_open
        False

        """
        return self._args[3]

    def _intersect(self, other):
        if not isinstance(other, Interval):
            return other.intersect(self)

        if not self._is_comparable(other):
            raise NotImplementedError("Intersection of intervals with symbolic "
                                      "end points is not yet implemented")

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

        return self.__class__(start, end, left_open, right_open)

    @property
    def _complement(self):
        a = Interval(S.NegativeInfinity, self.start, True, not self.left_open)
        b = Interval(self.end, S.Infinity, not self.right_open, True)
        return Union(a, b)

    def _contains(self, other):
        # We use the logic module here so that this method is meaningful
        # when used with symbolic end points.
        from sympy.logic.boolalg import And

        other = _sympify(other)

        if self.left_open:
            expr = other > self.start
        else:
            expr = other >= self.start

        if self.right_open:
            expr = And(expr, other < self.end)
        else:
            expr = And(expr, other <= self.end)

        return expr

    @property
    def _measure(self):
        return self.end - self.start

    def _eval_evalf(self, prec):
        return mpi(self.start.evalf(prec), self.end.evalf(prec))

    def _is_comparable(self, other):
        is_comparable = self.start.is_comparable
        is_comparable &= self.end.is_comparable
        is_comparable &= other.start.is_comparable
        is_comparable &= other.end.is_comparable

        return is_comparable

class Union(Set):
    """
    Represents a union of sets as a Set.

    Examples:
        >>> from sympy import *

        >>> Union(Interval(1, 2), Interval(3, 4))
        Union([1, 2], [3, 4])

        The Union constructor will always try to merge overlapping intervals,
        if possible. For example:

        >>> Union(Interval(1, 2), Interval(2, 3))
        [1, 3]

    """

    def __new__(cls, *args, **assumptions):
        intervals, other_sets = [], []
        for arg in args:
            if isinstance(arg, EmptySet):
                continue

            elif isinstance(arg, Interval):
                intervals.append(arg)

            elif isinstance(arg, Union):
                intervals += arg.args

            elif isinstance(arg, Set):
                other_sets.append(arg)

            else:
                raise ValueError, "Unknown argument '%s'" % arg

        # Any non-empty sets at all?
        if len(intervals) == 0 and len(other_sets) == 0:
            return S.EmptySet

        # Sort intervals according to their infimum
        intervals.sort(lambda i, j: cmp(i.start, j.start))

        # Merge comparable overlapping intervals
        i = 0
        while i < len(intervals) - 1:
            cur = intervals[i]
            next = intervals[i + 1]

            merge = False
            if cur._is_comparable(next):
                if next.start < cur.end:
                    merge = True
                elif next.start == cur.end:
                    # Must be careful with boundaries.
                    merge = not(next.left_open and cur.right_open)

            if merge:
                if cur.start == next.start:
                    left_open = cur.left_open and next.left_open
                else:
                    left_open = cur.left_open

                if cur.end < next.end:
                    right_open = next.right_open
                    end = next.end
                elif cur.end > next.end:
                    right_open = cur.right_open
                    end = cur.end
                else:
                    right_open = cur.right_open and next.right_open
                    end = cur.end

                intervals[i] = Interval(cur.start, end, left_open, right_open)
                del intervals[i + 1]
            else:
                i += 1

        # If a single set is left over, don't create a new Union object but
        # rather return the single set.
        if len(intervals) == 1 and len(other_sets) == 0:
            return intervals[0]
        elif len(intervals) == 0 and len(other_sets) == 1:
            return other_sets[0]

        return Basic.__new__(cls, *(intervals + other_sets), **assumptions)

    @property
    def _inf(self):
        # We use min_ so that sup is meaningful in combination with symbolic
        # interval end points.
        from sympy.functions.elementary.miscellaneous import min_

        inf = self.args[0].inf
        for set in self.args[1:]:
            inf = min_(inf, set.inf)
        return inf

    @property
    def _sup(self):
        # We use max_ so that sup is meaningful in combination with symbolic
        # end points.
        from sympy.functions.elementary.miscellaneous import max_

        sup = self.args[0].sup
        for set in self.args[1:]:
            sup = max_(sup, set.sup)
        return sup

    def _intersect(self, other):
        # Distributivity.
        if isinstance(other, Interval):
            intersections = []
            for interval in self.args:
                intersections.append(interval.intersect(other))
            return self.__class__(*intersections)

        elif isinstance(other, Union):
            intersections = []
            for interval in other.args:
                intersections.append(self.intersect(interval))
            return self.__class__(*intersections)

        else:
            return other.intersect(self)

    @property
    def _complement(self):
        # De Morgan's formula.
        complement = self.args[0].complement
        for set in self.args[1:]:
            complement = complement.intersect(set.complement)
        return complement

    def _contains(self, other):
        for set in self.args:
            if other in set:
                return True
        return False

    @property
    def _measure(self):
        measure = 0
        for set in self.args:
            measure += set.measure
        return measure

class EmptySet(Set):
    """
    Represents the empty set. The empty set is available as a singleton
    as S.EmptySet.

    Examples:
        >>> from sympy import *

        >>> S.EmptySet
        EmptySet()

        >>> Interval(1, 2).intersect(S.EmptySet)
        EmptySet()

    """

    __metaclass__ = SingletonMeta

    def _intersect(self, other):
        return S.EmptySet

    @property
    def _complement(self):
        return Interval(S.NegativeInfinity, S.Infinity)

    @property
    def _measure(self):
        return 0

    def _contains(self, other):
        return False
