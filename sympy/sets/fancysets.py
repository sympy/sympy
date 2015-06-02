from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.compatibility import as_int, with_metaclass, range
from sympy.sets.sets import Set, Interval, Intersection, EmptySet
from sympy.core.singleton import Singleton, S
from sympy.core.sympify import _sympify
from sympy.core.decorators import deprecated
from sympy.core.function import Lambda


class Naturals(with_metaclass(Singleton, Set)):
    """
    Represents the natural numbers (or counting numbers) which are all
    positive integers starting from 1. This set is also available as
    the Singleton, S.Naturals.

    Examples
    ========

    >>> from sympy import S, Interval, pprint
    >>> 5 in S.Naturals
    True
    >>> iterable = iter(S.Naturals)
    >>> next(iterable)
    1
    >>> next(iterable)
    2
    >>> next(iterable)
    3
    >>> pprint(S.Naturals.intersect(Interval(0, 10)))
    {1, 2, ..., 10}

    See Also
    ========
    Naturals0 : non-negative integers (i.e. includes 0, too)
    Integers : also includes negative integers
    """

    is_iterable = True
    _inf = S.One
    _sup = S.Infinity

    def _intersect(self, other):
        if other.is_Interval:
            return Intersection(
                S.Integers, other, Interval(self._inf, S.Infinity))
        return None

    def _contains(self, other):
        if other.is_positive and other.is_integer:
            return S.true
        elif other.is_integer is False or other.is_positive is False:
            return S.false

    def __iter__(self):
        i = self._inf
        while True:
            yield i
            i = i + 1

    @property
    def _boundary(self):
        return self


class Naturals0(Naturals):
    """Represents the whole numbers which are all the non-negative integers,
    inclusive of zero.

    See Also
    ========
    Naturals : positive integers; does not include 0
    Integers : also includes the negative integers
    """
    _inf = S.Zero

    def _contains(self, other):
        if other.is_integer and other.is_nonnegative:
            return S.true
        elif other.is_integer is False or other.is_nonnegative is False:
            return S.false


class Integers(with_metaclass(Singleton, Set)):
    """
    Represents all integers: positive, negative and zero. This set is also
    available as the Singleton, S.Integers.

    Examples
    ========

    >>> from sympy import S, Interval, pprint
    >>> 5 in S.Naturals
    True
    >>> iterable = iter(S.Integers)
    >>> next(iterable)
    0
    >>> next(iterable)
    1
    >>> next(iterable)
    -1
    >>> next(iterable)
    2

    >>> pprint(S.Integers.intersect(Interval(-4, 4)))
    {-4, -3, ..., 4}

    See Also
    ========
    Naturals0 : non-negative integers
    Integers : positive and negative integers and zero
    """

    is_iterable = True

    def _intersect(self, other):
        from sympy.functions.elementary.integers import floor, ceiling
        if other is Interval(S.NegativeInfinity, S.Infinity) or other is S.Reals:
            return self
        elif other.is_Interval:
            s = Range(ceiling(other.left), floor(other.right) + 1)
            return s.intersect(other)  # take out endpoints if open interval
        return None

    def _contains(self, other):
        if other.is_integer:
            return S.true
        elif other.is_integer is False:
            return S.false

    def __iter__(self):
        yield S.Zero
        i = S(1)
        while True:
            yield i
            yield -i
            i = i + 1

    @property
    def _inf(self):
        return -S.Infinity

    @property
    def _sup(self):
        return S.Infinity

    @property
    def _boundary(self):
        return self

    def _eval_imageset(self, f):
        from sympy import Wild
        expr = f.expr
        if len(f.variables) > 1:
            return
        n = f.variables[0]

        a = Wild('a')
        b = Wild('b')

        match = expr.match(a*n + b)
        if match[a].is_negative:
            expr = -expr

        match = expr.match(a*n + b)
        if match[a] is S.One and match[b].is_integer:
            expr = expr - match[b]

        return ImageSet(Lambda(n, expr), S.Integers)


class Reals(with_metaclass(Singleton, Interval)):

    def __new__(cls):
        return Interval.__new__(cls, -S.Infinity, S.Infinity)


class ImageSet(Set):
    """
    Image of a set under a mathematical function

    Examples
    ========

    >>> from sympy import Symbol, S, ImageSet, FiniteSet, Lambda

    >>> x = Symbol('x')
    >>> N = S.Naturals
    >>> squares = ImageSet(Lambda(x, x**2), N) # {x**2 for x in N}
    >>> 4 in squares
    True
    >>> 5 in squares
    False

    >>> FiniteSet(0, 1, 2, 3, 4, 5, 6, 7, 9, 10).intersect(squares)
    {1, 4, 9}

    >>> square_iterable = iter(squares)
    >>> for i in range(4):
    ...     next(square_iterable)
    1
    4
    9
    16
    """
    def __new__(cls, lamda, base_set):
        return Basic.__new__(cls, lamda, base_set)

    lamda = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def __iter__(self):
        already_seen = set()
        for i in self.base_set:
            val = self.lamda(i)
            if val in already_seen:
                continue
            else:
                already_seen.add(val)
                yield val

    def _is_multivariate(self):
        return len(self.lamda.variables) > 1

    def _contains(self, other):
        from sympy.solvers import solve
        L = self.lamda
        if self._is_multivariate():
            solns = solve([expr - val for val, expr in zip(other, L.expr)],
                    L.variables)
        else:
            solns = solve(L.expr - other, L.variables[0])

        for soln in solns:
            try:
                if soln in self.base_set:
                    return S.true
            except TypeError:
                if soln.evalf() in self.base_set:
                    return S.true
        return S.false

    def _intersect(self, other):
        from sympy import solve, Dummy, imageset, Union

        if isinstance(other, ImageSet):
            return

        f = self.lamda
        if len(f.variables) > 1:
            return
        var = f.variables[0]
        expr = f.expr
        if not expr.is_polynomial(var):
            return
        d = Dummy()
        solns = solve(expr - d, var)
        return Union(imageset(f, self.base_set.intersect(
            imageset(d, f_inv, other))) for f_inv in solns)

    @property
    def is_iterable(self):
        return self.base_set.is_iterable

    def _intersect(self, other):
        from sympy import Dummy
        from sympy.solvers.diophantine import diophantine
        from sympy.sets.sets import imageset
        if self.base_set is S.Integers:
            if isinstance(other, ImageSet) and other.base_set is S.Integers:
                f, g = self.lamda.expr, other.lamda.expr
                n, m = self.lamda.variables[0], other.lamda.variables[0]

                # Diophantine sorts the solutions according to the alphabetic
                # order of the variable names, since the result should not depend
                # on the variable name, they are replaced by the dummy variables
                # below
                a, b = Dummy('a'), Dummy('b')
                f, g = f.subs(n, a), g.subs(m, b)
                solns_set = diophantine(f - g)
                if solns_set == set():
                    return EmptySet()
                solns = list(diophantine(f - g))
                if len(solns) == 1:
                    t = list(solns[0][0].free_symbols)[0]
                else:
                    return None

                # since 'a' < 'b'
                return imageset(Lambda(t, f.subs(a, solns[0][0])), S.Integers)

        if other == S.Reals:
            from sympy.solvers.solveset import solveset_real
            from sympy.core.function import expand_complex
            if len(self.lamda.variables) > 1:
                return None

            f = self.lamda.expr
            n = self.lamda.variables[0]

            n_ = Dummy(n.name, real=True)
            f_ = f.subs(n, n_)

            re, im = f_.as_real_imag()
            im = expand_complex(im)

            return imageset(Lambda(n_, re),
                            self.base_set.intersect(
                                solveset_real(im, n_)))


@deprecated(useinstead="ImageSet", issue=7057, deprecated_since_version="0.7.4")
def TransformationSet(*args, **kwargs):
    """Deprecated alias for the ImageSet constructor."""
    return ImageSet(*args, **kwargs)


class Range(Set):
    """
    Represents a range of integers.

    Examples
    ========

    >>> from sympy import Range
    >>> list(Range(5)) # 0 to 5
    [0, 1, 2, 3, 4]
    >>> list(Range(10, 15)) # 10 to 15
    [10, 11, 12, 13, 14]
    >>> list(Range(10, 20, 2)) # 10 to 20 in steps of 2
    [10, 12, 14, 16, 18]
    >>> list(Range(20, 10, -2)) # 20 to 10 backward in steps of 2
    [12, 14, 16, 18, 20]

    """

    is_iterable = True

    def __new__(cls, *args):
        from sympy.functions.elementary.integers import ceiling
        # expand range
        slc = slice(*args)
        start, stop, step = slc.start or 0, slc.stop, slc.step or 1
        try:
            start, stop, step = [w if w in [S.NegativeInfinity, S.Infinity] else S(as_int(w))
                                 for w in (start, stop, step)]
        except ValueError:
            raise ValueError("Inputs to Range must be Integer Valued\n" +
                    "Use ImageSets of Ranges for other cases")

        if not step.is_finite:
            raise ValueError("Infinite step is not allowed")
        if start == stop:
            return S.EmptySet

        n = ceiling((stop - start)/step)
        if n <= 0:
            return S.EmptySet

        # normalize args: regardless of how they are entered they will show
        # canonically as Range(inf, sup, step) with step > 0
        if n.is_finite:
            start, stop = sorted((start, start + (n - 1)*step))
        else:
            start, stop = sorted((start, stop - step))

        step = abs(step)
        if (start, stop) == (S.NegativeInfinity, S.Infinity):
            raise ValueError("Both the start and end value of "
                             "Range cannot be unbounded")
        else:
            return Basic.__new__(cls, start, stop + step, step)

    start = property(lambda self: self.args[0])
    stop = property(lambda self: self.args[1])
    step = property(lambda self: self.args[2])

    def _intersect(self, other):
        from sympy.functions.elementary.integers import floor, ceiling
        from sympy.functions.elementary.miscellaneous import Min, Max
        if other.is_Interval:
            osup = other.sup
            oinf = other.inf
            # if other is [0, 10) we can only go up to 9
            if osup.is_integer and other.right_open:
                osup -= 1
            if oinf.is_integer and other.left_open:
                oinf += 1

            # Take the most restrictive of the bounds set by the two sets
            # round inwards
            inf = ceiling(Max(self.inf, oinf))
            sup = floor(Min(self.sup, osup))
            # if we are off the sequence, get back on
            if inf.is_finite and self.inf.is_finite:
                off = (inf - self.inf) % self.step
            else:
                off = S.Zero
            if off:
                inf += self.step - off

            return Range(inf, sup + 1, self.step)

        if other == S.Naturals:
            return self._intersect(Interval(1, S.Infinity))

        if other == S.Integers:
            return self

        return None

    def _contains(self, other):
        if (((self.start - other)/self.step).is_integer or
            ((self.stop - other)/self.step).is_integer):
            return _sympify(other >= self.inf and other <= self.sup)
        elif (((self.start - other)/self.step).is_integer is False and
            ((self.stop - other)/self.step).is_integer is False):
            return S.false

    def __iter__(self):
        if self.start is S.NegativeInfinity:
            i = self.stop - self.step
            step = -self.step
        else:
            i = self.start
            step = self.step

        while(i < self.stop and i >= self.start):
            yield i
            i += step

    def __len__(self):
        return (self.stop - self.start)//self.step

    def __nonzero__(self):
        return True

    __bool__ = __nonzero__

    def _ith_element(self, i):
        return self.start + i*self.step

    @property
    def _last_element(self):
        if self.stop is S.Infinity:
            return S.Infinity
        elif self.start is S.NegativeInfinity:
            return self.stop - self.step
        else:
            return self._ith_element(len(self) - 1)

    @property
    def _inf(self):
        return self.start

    @property
    def _sup(self):
        return self.stop - self.step

    @property
    def _boundary(self):
        return self
