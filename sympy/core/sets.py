from __future__ import print_function, division

from itertools import product

from sympy.core.sympify import _sympify, sympify
from sympy.core.basic import Basic
from sympy.core.singleton import Singleton, S
from sympy.core.evalf import EvalfMixin
from sympy.core.numbers import Float
from sympy.core.compatibility import iterable, with_metaclass

from sympy.mpmath import mpi, mpf
from sympy.logic.boolalg import And, Or, true, false

from sympy.utilities import default_sort_key


class Set(Basic):
    """
    The base class for any kind of set.

    This is not meant to be used directly as a container of items.
    It does not behave like the builtin set; see FiniteSet for that.

    Real intervals are represented by the Interval class and unions of sets
    by the Union class. The empty set is represented by the EmptySet class
    and available as a singleton as S.EmptySet.
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

    def sort_key(self, order=None):
        """
        Give sort_key of infimum (if possible) else sort_key of the set.
        """
        try:
            infimum = self.inf
            if infimum.is_comparable:
                return default_sort_key(infimum, order)
        except (NotImplementedError, ValueError):
            pass
        args = tuple([default_sort_key(a, order) for a in self._sorted_args])
        return self.class_key(), (len(args), args), S.One.class_key(), S.One

    def union(self, other):
        """
        Returns the union of 'self' and 'other'.

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

        """
        return Intersection(self, other)

    def _intersect(self, other):
        """
        This function should only be used internally

        self._intersect(other) returns a new, intersected set if self knows how
        to intersect itself with other, otherwise it returns None

        When making a new set class you can be assured that other will not
        be a Union, FiniteSet, or EmptySet

        Used within the Intersection class
        """
        return None

    def _union(self, other):
        """
        This function should only be used internally

        self._union(other) returns a new, joined set if self knows how
        to join itself with other, otherwise it returns None.
        It may also return a python set of SymPy Sets if they are somehow
        simpler. If it does this it must be idempotent i.e. the sets returned
        must return None with _union'ed with each other

        Used within the Union class
        """
        return None

    @property
    def complement(self):
        """
        The complement of 'self'.

        As a shortcut it is possible to use the '~' or '-' operators:

        >>> from sympy import Interval

        >>> Interval(0, 1).complement
        (-oo, 0) U (1, oo)
        >>> ~Interval(0, 1)
        (-oo, 0) U (1, oo)
        >>> -Interval(0, 1)
        (-oo, 0) U (1, oo)

        """
        return self._complement

    @property
    def _complement(self):
        raise NotImplementedError("(%s)._complement" % self)

    @property
    def inf(self):
        """
        The infimum of 'self'

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

        >>> from sympy import Interval

        >>> Interval(0, 1).contains(0.5)
        True
        >>> 0.5 in Interval(0, 1)
        True

        """
        c = self._contains(sympify(other, strict=True))
        if c in (true, false):
            # TODO: would we want to return the Basic type here?
            return bool(c)
        return c

    def _contains(self, other):
        raise NotImplementedError("(%s)._contains(%s)" % (self, other))

    def subset(self, other):
        """
        Returns True if 'other' is a subset of 'self'.

        >>> from sympy import Interval

        >>> Interval(0, 1).subset(Interval(0, 0.5))
        True
        >>> Interval(0, 1, left_open=True).subset(Interval(0, 1))
        False

        """
        if isinstance(other, Set):
            return self.intersect(other) == other
        else:
            raise ValueError("Unknown argument '%s'" % other)

    @property
    def measure(self):
        """
        The (Lebesgue) measure of 'self'

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
        return self.subset(self.boundary)

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

    def __pow__(self, exp):
        if not sympify(exp).is_Integer and exp >= 0:
            raise ValueError("%s: Exponent must be a positive Integer" % exp)
        return ProductSet([self]*exp)

    def __sub__(self, other):
        return self.intersect(other.complement)

    def __neg__(self):
        return self.complement

    def __invert__(self):
        return self.complement

    def __contains__(self, other):
        from sympy.assumptions import ask
        symb = self.contains(other)
        result = ask(symb)
        if result is None:
            raise TypeError('contains did not evaluate to a bool: %r' % symb)
        return result

    @property
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
    http://en.wikipedia.org/wiki/Cartesian_product
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

    def _contains(self, element):
        """
        'in' operator for ProductSets

        >>> from sympy import Interval

        >>> (2, 3) in Interval(0, 5) * Interval(0, 5)
        True

        >>> (10, 10) in Interval(0, 5) * Interval(0, 5)
        False

        Passes operation on to constituent sets
        """
        try:
            if len(element) != len(self.args):
                return False
        except TypeError:  # maybe element isn't an iterable
            return False
        return And(*[set.contains(item) for set, item in zip(self.sets, element)])

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
    def _complement(self):
        # For each set consider it or it's complement
        # We need at least one of the sets to be complemented
        # Consider all 2^n combinations.
        # We can conveniently represent these options easily using a ProductSet
        switch_sets = ProductSet(FiniteSet(s, s.complement) for s in self.sets)
        product_sets = (ProductSet(*set) for set in switch_sets)
        # Union of all combinations but this one
        return Union(p for p in product_sets if p != self)

    @property
    def _boundary(self):
        return Union(ProductSet(b + b.boundary if i != j else b.boundary
                                for j, b in enumerate(self.sets))
                                for i, a in enumerate(self.sets))


    @property
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

    <http://en.wikipedia.org/wiki/Interval_(mathematics)>
    """
    is_Interval = True
    is_real = True

    def __new__(cls, start, end, left_open=False, right_open=False):

        start = _sympify(start)
        end = _sympify(end)

        inftys = [S.Infinity, S.NegativeInfinity]
        # Only allow real intervals (use symbols with 'is_real=True').
        if not (start.is_real or start in inftys) or not (end.is_real or end in inftys):
            raise ValueError("Only real intervals are supported")

        # Make sure that the created interval will be valid.
        if end.is_comparable and start.is_comparable:
            if end < start:
                return S.EmptySet

        if end == start and (left_open or right_open):
            return S.EmptySet
        if end == start and not (left_open or right_open):
            return FiniteSet(end)

        # Make sure infinite interval end points are open.
        if start == S.NegativeInfinity:
            left_open = True
        if end == S.Infinity:
            right_open = True

        return Basic.__new__(cls, start, end, left_open, right_open)

    @property
    def start(self):
        """
        The left end point of 'self'.

        This property takes the same value as the 'inf' property.

        >>> from sympy import Interval

        >>> Interval(0, 1).start
        0

        """
        return self._args[0]

    _inf = left = start

    @property
    def end(self):
        """
        The right end point of 'self'.

        This property takes the same value as the 'sup' property.

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
        if ((self.left_open and other.contains(self.start) is True) or
                (self.right_open and other.contains(self.end) is True)):
            # Fill in my end points and return
            open_left = self.left_open and self.start not in other
            open_right = self.right_open and self.end not in other
            new_self = Interval(self.start, self.end, open_left, open_right)
            return set((new_self, other))

        return None

    @property
    def _complement(self):
        a = Interval(S.NegativeInfinity, self.start, True, not self.left_open)
        b = Interval(self.end, S.Infinity, not self.right_open, True)
        return Union(a, b)

    @property
    def _boundary(self):
        return FiniteSet(self.start, self.end)

    def _contains(self, other):
        if self.left_open:
            expr = other > self.start
        else:
            expr = other >= self.start

        if self.right_open:
            expr = And(expr, other < self.end)
        else:
            expr = And(expr, other <= self.end)

        return expr

    def _eval_imageset(self, f):
        from sympy import Dummy
        from sympy.functions.elementary.miscellaneous import Min, Max
        from sympy.solvers import solve
        from sympy.core.function import diff
        from sympy.series import limit
        from sympy.calculus.singularities import singularities
        # TODO: handle piecewise defined functions
        # TODO: handle functions with infinitely many solutions (eg, sin, tan)
        # TODO: handle multivariate functions

        # var and expr are being defined this way to
        # support Python lambda and not just sympy Lambda
        try:
            var = Dummy()
            expr = f(var)
            if len(expr.free_symbols) > 1:
                raise TypeError
        except TypeError:
            raise NotImplementedError("Sorry, Multivariate imagesets are"
                                      " not yet implemented, you are welcome"
                                      " to add this feature in Sympy")

        if not self.start.is_comparable or not self.end.is_comparable:
            raise NotImplementedError("Sets with non comparable/variable"
                                      " arguments are not supported")

        sing = [x for x in singularities(expr, var) if x.is_real and x in self]

        if self.left_open:
            _start = limit(expr, var, self.start, dir="+")
        elif self.start not in sing:
            _start = f(self.start)
        if self.right_open:
            _end = limit(expr, var, self.end, dir="-")
        elif self.end not in sing:
            _end = f(self.end)

        if len(sing) == 0:
            solns = solve(diff(expr, var), var)

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
                Union(*[imageset(f, Interval(sing[i], sing[i + 1]), True, True)
                        for i in range(1, len(sing) - 1)]) + \
                imageset(f, Interval(sing[-1], self.end, True, self.right_open))

    @property
    def _measure(self):
        return self.end - self.start

    def to_mpi(self, prec=53):
        return mpi(mpf(self.start.evalf(prec)), mpf(self.end.evalf(prec)))

    def _eval_evalf(self, prec):
        return Interval(self.left.evalf(), self.right.evalf(),
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

    def as_relational(self, symbol):
        """Rewrite an interval in terms of inequalities and logic operators. """
        other = sympify(symbol)
        if self.right_open:
            right = other < self.end
        else:
            right = other <= self.end
        if right is True:
            if self.left_open:
                return other > self.start
            else:
                return other >= self.start
        if self.left_open:
            left = self.start < other
        else:
            left = self.start <= other
        return And(left, right)

    @property
    def free_symbols(self):
        return self.start.free_symbols | self.end.free_symbols

class Union(Set, EvalfMixin):
    """
    Represents a union of sets as a Set.

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
    <http://en.wikipedia.org/wiki/Union_(set_theory)>
    """
    is_Union = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', True)

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

        args = sorted(args, key=default_sort_key)

        # Reduce sets using known rules
        if evaluate:
            return Union.reduce(args)

        return Basic.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a Union using known rules

        We first start with global rules like
        'Merge all FiniteSets'

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent
        """

        # ===== Global Rules =====
        # Merge all finite sets
        finite_sets = [x for x in args if x.is_FiniteSet]
        if len(finite_sets) > 1:
            finite_set = FiniteSet(x for set in finite_sets for x in set)
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

    @property
    def _complement(self):
        # De Morgan's formula.
        complement = self.args[0].complement
        for set in self.args[1:]:
            complement = complement.intersect(set.complement)
        return complement

    def _contains(self, other):
        or_args = [the_set.contains(other) for the_set in self.args]
        return Or(*or_args)

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
            return Union(set.evalf() for set in self.args)
        except:
            raise TypeError("Not all sets are evalf-able")

    def __iter__(self):
        import itertools
        if all(set.is_iterable for set in self.args):
            return itertools.chain(*(iter(arg) for arg in self.args))
        else:
            raise TypeError("Not all constituent sets are iterable")

    @property
    def is_real(self):
        return all(set.is_real for set in self.args)


class Intersection(Set):
    """
    Represents an intersection of sets as a Set.

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
    <http://en.wikipedia.org/wiki/Intersection_(set_theory)>
    """
    is_Intersection = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', True)

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

        # Intersection of no sets is everything
        if len(args) == 0:
            return S.UniversalSet

        args = sorted(args, key=default_sort_key)

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

    @property
    def _complement(self):
        raise NotImplementedError()

    def _eval_imageset(self, f):
        return Intersection(imageset(f, arg) for arg in self.args)

    def _contains(self, other):
        from sympy.logic.boolalg import And
        return And(*[set.contains(other) for set in self.args])

    def __iter__(self):
        for s in self.args:
            if s.is_iterable:
                other_sets = set(self.args) - set((s,))
                other = Intersection(other_sets, evaluate=False)
                return (x for x in s if x in other)

        raise ValueError("None of the constituent sets are iterable")

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

        # If any FiniteSets see which elements of that finite set occur within
        # all other sets in the intersection
        for s in args:
            if s.is_FiniteSet:
                return s.__class__(x for x in s
                        if all(x in other for other in args))

        # If any of the sets are unions, return a Union of Intersections
        for s in args:
            if s.is_Union:
                other_sets = set(args) - set((s,))
                other = Intersection(other_sets)
                return Union(Intersection(arg, other) for arg in s.args)

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
    http://en.wikipedia.org/wiki/Empty_set
    """
    is_EmptySet = True

    def _intersect(self, other):
        return S.EmptySet

    @property
    def _complement(self):
        return S.UniversalSet

    @property
    def _measure(self):
        return 0

    def _contains(self, other):
        return False

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

    @property
    def _boundary(self):
        return self

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
    http://en.wikipedia.org/wiki/Universal_set
    """

    is_UniversalSet = True

    def _intersect(self, other):
        return other

    @property
    def _complement(self):
        return S.EmptySet

    @property
    def _measure(self):
        return S.Infinity

    def _contains(self, other):
        return True

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

    References
    ==========
    http://en.wikipedia.org/wiki/Finite_set
    """
    is_FiniteSet = True
    is_iterable = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', True)
        if evaluate:
            if len(args) == 1 and iterable(args[0]):
                args = args[0]

            args = list(map(sympify, args))

            if len(args) == 0:
                return EmptySet()


        args = frozenset(args)  # remove duplicates
        obj = Basic.__new__(cls, *args)
        obj._elements = args
        return obj

    def __iter__(self):
        return iter(self.args)

    def _intersect(self, other):
        """
        This function should only be used internally

        See Set._intersect for docstring
        """
        if isinstance(other, self.__class__):
            return self.__class__(*(self._elements & other._elements))
        return self.__class__(el for el in self if el in other)

    def _union(self, other):
        """
        This function should only be used internally

        See Set._union for docstring
        """
        if other.is_FiniteSet:
            return FiniteSet(*(self._elements | other._elements))

        # If other set contains one of my elements, remove it from myself
        if any(other.contains(x) is True for x in self):
            return set((
                FiniteSet(x for x in self if other.contains(x) is not True),
                other))

        return None

    def _contains(self, other):
        """
        Tests whether an element, other, is in the set.

        Relies on Python's set class. This tests for object equality
        All inputs are sympified

        >>> from sympy import FiniteSet

        >>> 1 in FiniteSet(1, 2)
        True
        >>> 5 in FiniteSet(1, 2)
        False

        """
        return other in self._elements

    def _eval_imageset(self, f):
        return FiniteSet(*map(f, self))

    @property
    def _complement(self):
        """
        The complement of a real finite set is the Union of open Intervals
        between the elements of the set.

        >>> from sympy import FiniteSet
        >>> FiniteSet(1, 2, 3).complement
        (-oo, 1) U (1, 2) U (2, 3) U (3, oo)


        """
        if not all(elem.is_number for elem in self):
            raise ValueError("%s: Complement not defined for symbolic inputs"
                    % self)

        # as there are only numbers involved, a straight sort is sufficient;
        # default_sort_key is not needed
        args = sorted(self.args)

        intervals = []  # Build up a list of intervals between the elements
        intervals += [Interval(S.NegativeInfinity, args[0], True, True)]
        for a, b in zip(args[:-1], args[1:]):
            intervals.append(Interval(a, b, True, True))  # open intervals
        intervals.append(Interval(args[-1], S.Infinity, True, True))
        return Union(intervals, evaluate=False)

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

    def __sub__(self, other):
        return FiniteSet(el for el in self if el not in other)

    def as_relational(self, symbol):
        """Rewrite a FiniteSet in terms of equalities and logic operators. """
        from sympy.core.relational import Eq
        return Or(*[Eq(symbol, elem) for elem in self])

    @property
    def is_real(self):
        return all(el.is_real for el in self)

    def compare(self, other):
        return (hash(self) - hash(other))

    def _eval_evalf(self, prec):
        return FiniteSet(elem.evalf(prec) for elem in self)

    def _hashable_content(self):
        return (self._elements,)

    @property
    def _sorted_args(self):
        from sympy.utilities import default_sort_key
        return sorted(self.args, key=default_sort_key)

    def __ge__(self, other):
        return self.subset(other)

    def __gt__(self, other):
        return self != other and self >= other

    def __le__(self, other):
        return other.subset(self)

    def __lt__(self, other):
        return self != other and other >= self


def imageset(*args):
    """ Image of set under transformation ``f``

    .. math::
        { f(x) | x \in self }

    Examples
    ========

    >>> from sympy import Interval, Symbol, imageset
    >>> x = Symbol('x')

    >>> imageset(x, 2*x, Interval(0, 2))
    [0, 4]

    >>> imageset(lambda x: 2*x, Interval(0, 2))
    [0, 4]

    See Also:
        ImageSet
    """
    if len(args) == 3:
        from sympy import Lambda
        f = Lambda(*args[:2])
    else:
        f = args[0]
    set = args[-1]

    return set._eval_imageset(f)
