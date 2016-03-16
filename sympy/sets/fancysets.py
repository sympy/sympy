from __future__ import print_function, division

from sympy.logic.boolalg import And
from sympy.core.basic import Basic
from sympy.core.compatibility import as_int, with_metaclass, range, PY3
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)
from sympy.core.sympify import _sympify, sympify, converter
from sympy.core.singleton import Singleton, S
from sympy.core.function import Lambda
from sympy.utilities.misc import filldedent, func_name


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
        i = S.One
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

    def __eq__(self, other):
        return other == Interval(-S.Infinity, S.Infinity)

    def __hash__(self):
        return hash(Interval(-S.Infinity, S.Infinity))


class ImageSet(Set):
    """
    Image of a set under a mathematical function. The transformation
    must be given as a Lambda function which has as many arguments
    as the elements of the set upon which it operates, e.g. 1 argument
    when acting on the set of integers or 2 arguments when acting on
    a complex region.

    This function is not normally called directly, but is called
    from `imageset`.


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

    See Also
    ========
    sympy.sets.sets.imageset
    """
    def __new__(cls, lamda, base_set):
        if not isinstance(lamda, Lambda):
            raise ValueError('first argument must be a Lambda')
        if lamda is S.IdentityFunction:
            return base_set
        if not lamda.expr.free_symbols or not lamda.expr.args:
            return FiniteSet(lamda.expr)

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
        from sympy.matrices import Matrix
        from sympy.solvers.solveset import solveset, linsolve
        from sympy.utilities.iterables import iterable, cartes
        L = self.lamda
        if self._is_multivariate():
            if not iterable(L.expr):
                if iterable(other):
                    return S.false
                return other.as_numer_denom() in self.func(
                    Lambda(L.variables, L.expr.as_numer_denom()), self.base_set)
            if len(L.expr) != len(self.lamda.variables):
                raise NotImplementedError(filldedent('''
    Dimensions of input and output of Lambda are different.'''))
            eqs = [expr - val for val, expr in zip(other, L.expr)]
            variables = L.variables
            free = set(variables)
            if all(i.is_number for i in list(Matrix(eqs).jacobian(variables))):
                solns = list(linsolve([e - val for e, val in
                zip(L.expr, other)], variables))
            else:
                syms = [e.free_symbols & free for e in eqs]
                solns = {}
                for i, (e, s, v) in enumerate(zip(eqs, syms, other)):
                    if not s:
                        if e != v:
                            return S.false
                        solns[vars[i]] = [v]
                        continue
                    elif len(s) == 1:
                        sy = s.pop()
                        sol = solveset(e, sy)
                        if sol is S.EmptySet:
                            return S.false
                        elif isinstance(sol, FiniteSet):
                            solns[sy] = list(sol)
                        else:
                            raise NotImplementedError
                    else:
                        raise NotImplementedError
                solns = cartes(*[solns[s] for s in variables])
        else:
            # assume scalar -> scalar mapping
            solnsSet = solveset(L.expr - other, L.variables[0])
            if solnsSet.is_FiniteSet:
                solns = list(solnsSet)
            else:
                raise NotImplementedError(filldedent('''
                Determining whether an ImageSet contains %s has not
                been implemented.''' % func_name(other)))
        for soln in solns:
            try:
                if soln in self.base_set:
                    return S.true
            except TypeError:
                return self.base_set.contains(soln.evalf())
        return S.false

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
        if len(args) == 1:
            if isinstance(args[0], range if PY3 else xrange):
                args = args[0].__reduce__()[1]  # use pickle method

        # expand range
        slc = slice(*args)
        start, stop, step = slc.start or 0, slc.stop, slc.step or 1
        try:
            start, stop, step = [w if w in [S.NegativeInfinity, S.Infinity] else sympify(as_int(w))
                                 for w in (start, stop, step)]
        except ValueError:
            raise ValueError(filldedent('''
    Inputs to Range must be integers;
    use ImageSets of Ranges for other cases, e.g.
    `ImageSet(Lambda(i, i/10), Range(2))` to give [0, .1].'''))

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


if PY3:
    converter[range] = Range
else:
    converter[xrange] = Range

def normalize_theta_set(theta):
    """
    Normalize a Real Set `theta` in the Interval [0, 2*pi). It returns
    a normalized value of theta in the Set. For Interval, a maximum of
    one cycle [0, 2*pi], is returned i.e. for theta equal to [0, 10*pi],
    returned normalized value would be [0, 2*pi). As of now intervals
    with end points as non-multiples of `pi` is not supported.

    Raises
    ======

    NotImplementedError
        The algorithms for Normalizing theta Set are not yet
        implemented.
    ValueError
        The input is not valid, i.e. the input is not a real set.
    RuntimeError
        It is a bug, please report to the github issue tracker.

    Examples
    ========

    >>> from sympy.sets.fancysets import normalize_theta_set
    >>> from sympy import Interval, FiniteSet, pi
    >>> normalize_theta_set(Interval(9*pi/2, 5*pi))
    [pi/2, pi]
    >>> normalize_theta_set(Interval(-3*pi/2, pi/2))
    [0, 2*pi)
    >>> normalize_theta_set(Interval(-pi/2, pi/2))
    [0, pi/2] U [3*pi/2, 2*pi)
    >>> normalize_theta_set(Interval(-4*pi, 3*pi))
    [0, 2*pi)
    >>> normalize_theta_set(Interval(-3*pi/2, -pi/2))
    [pi/2, 3*pi/2]
    >>> normalize_theta_set(FiniteSet(0, pi, 3*pi))
    {0, pi}

    """
    from sympy.functions.elementary.trigonometric import _pi_coeff as coeff

    if theta.is_Interval:
        interval_len = theta.measure
        # one complete circle
        if interval_len >= 2*S.Pi:
            if interval_len == 2*S.Pi and theta.left_open and theta.right_open:
                k = coeff(theta.start)
                return Union(Interval(0, k*S.Pi, False, True),
                        Interval(k*S.Pi, 2*S.Pi, True, True))
            return Interval(0, 2*S.Pi, False, True)

        k_start, k_end = coeff(theta.start), coeff(theta.end)

        if k_start is None or k_end is None:
            raise NotImplementedError("Normalizing theta without pi as coefficient is "
                                    "not yet implemented")
        new_start = k_start*S.Pi
        new_end = k_end*S.Pi

        if new_start > new_end:
            return Union(Interval(S.Zero, new_end, False, theta.right_open),
                         Interval(new_start, 2*S.Pi, theta.left_open, True))
        else:
            return Interval(new_start, new_end, theta.left_open, theta.right_open)

    elif theta.is_FiniteSet:
        new_theta = []
        for element in theta:
            k = coeff(element)
            if k is None:
                raise NotImplementedError('Normalizing theta without pi as '
                                          'coefficient, is not Implemented.')
            else:
                new_theta.append(k*S.Pi)
        return FiniteSet(*new_theta)

    elif theta.is_Union:
        return Union(*[normalize_theta_set(interval) for interval in theta.args])

    elif theta.is_subset(S.Reals):
        raise NotImplementedError("Normalizing theta when, it is of type %s is not "
                                  "implemented" % type(theta))
    else:
        raise ValueError(" %s is not a real set" % (theta))


class ComplexRegion(Set):
    """
    Represents the Set of all Complex Numbers. It can represent a
    region of Complex Plane in both the standard forms Polar and
    Rectangular coordinates.

    * Polar Form
      Input is in the form of the ProductSet or Union of ProductSets
      of the intervals of r and theta, & use the flag polar=True.

    Z = {z in C | z = r*[cos(theta) + I*sin(theta)], r in [r], theta in [theta]}

    * Rectangular Form
      Input is in the form of the ProductSet or Union of ProductSets
      of interval of x and y the of the Complex numbers in a Plane.
      Default input type is in rectangular form.

    Z = {z in C | z = x + I*y, x in [Re(z)], y in [Im(z)]}

    Examples
    ========

    >>> from sympy.sets.fancysets import ComplexRegion
    >>> from sympy.sets import Interval
    >>> from sympy import S, I, Union
    >>> a = Interval(2, 3)
    >>> b = Interval(4, 6)
    >>> c = Interval(1, 8)
    >>> c1 = ComplexRegion(a*b)  # Rectangular Form
    >>> c1
    ComplexRegion([2, 3] x [4, 6], False)

    * c1 represents the rectangular region in complex plane
      surrounded by the coordinates (2, 4), (3, 4), (3, 6) and
      (2, 6), of the four vertices.

    >>> c2 = ComplexRegion(Union(a*b, b*c))
    >>> c2
    ComplexRegion([2, 3] x [4, 6] U [4, 6] x [1, 8], False)

    * c2 represents the Union of two rectangular regions in complex
      plane. One of them surrounded by the coordinates of c1 and
      other surrounded by the coordinates (4, 1), (6, 1), (6, 8) and
      (4, 8).

    >>> 2.5 + 4.5*I in c1
    True
    >>> 2.5 + 6.5*I in c1
    False

    >>> r = Interval(0, 1)
    >>> theta = Interval(0, 2*S.Pi)
    >>> c2 = ComplexRegion(r*theta, polar=True)  # Polar Form
    >>> c2  # unit Disk
    ComplexRegion([0, 1] x [0, 2*pi), True)

    * c2 represents the region in complex plane inside the
      Unit Disk centered at the origin.

    >>> 0.5 + 0.5*I in c2
    True
    >>> 1 + 2*I in c2
    False

    >>> unit_disk = ComplexRegion(Interval(0, 1)*Interval(0, 2*S.Pi), polar=True)
    >>> upper_half_unit_disk = ComplexRegion(Interval(0, 1)*Interval(0, S.Pi), polar=True)
    >>> intersection = unit_disk.intersect(upper_half_unit_disk)
    >>> intersection
    ComplexRegion([0, 1] x [0, pi], True)
    >>> intersection == upper_half_unit_disk
    True

    See Also
    ========

    Reals

    """
    is_ComplexRegion = True

    def __new__(cls, sets, polar=False):
        from sympy import symbols, Dummy, sympify, sin, cos

        x, y, r, theta = symbols('x, y, r, theta', cls=Dummy)
        I = S.ImaginaryUnit
        polar = sympify(polar)

        # Rectangular Form
        if polar == False:
            if all(_a.is_FiniteSet for _a in sets.args) and (len(sets.args) == 2):

                # ** ProductSet of FiniteSets in the Complex Plane. **
                # For Cases like ComplexRegion({2, 4}*{3}), It
                # would return {2 + 3*I, 4 + 3*I}
                complex_num = []
                for x in sets.args[0]:
                    for y in sets.args[1]:
                        complex_num.append(x + I*y)
                obj = FiniteSet(*complex_num)
            else:
                obj = ImageSet.__new__(cls, Lambda((x, y), x + I*y), sets)
            obj._variables = (x, y)
            obj._expr = x + I*y

        # Polar Form
        elif polar == True:
            new_sets = []
            # sets is Union of ProductSets
            if not sets.is_ProductSet:
                for k in sets.args:
                    new_sets.append(k)
            # sets is ProductSets
            else:
                new_sets.append(sets)
            # Normalize input theta
            for k, v in enumerate(new_sets):
                from sympy.sets import ProductSet
                new_sets[k] = ProductSet(v.args[0],
                                         normalize_theta_set(v.args[1]))
            sets = Union(*new_sets)
            obj = ImageSet.__new__(cls, Lambda((r, theta),
                                   r*(cos(theta) + I*sin(theta))),
                                   sets)
            obj._variables = (r, theta)
            obj._expr = r*(cos(theta) + I*sin(theta))

        else:
            raise ValueError("polar should be either True or False")

        obj._sets = sets
        obj._polar = polar
        return obj

    @property
    def sets(self):
        """
        Return raw input sets to the self.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, Union
        >>> a = Interval(2, 3)
        >>> b = Interval(4, 5)
        >>> c = Interval(1, 7)
        >>> C1 = ComplexRegion(a*b)
        >>> C1.sets
        [2, 3] x [4, 5]
        >>> C2 = ComplexRegion(Union(a*b, b*c))
        >>> C2.sets
        [2, 3] x [4, 5] U [4, 5] x [1, 7]

        """
        return self._sets

    @property
    def args(self):
        return (self._sets, self._polar)

    @property
    def variables(self):
        return self._variables

    @property
    def expr(self):
        return self._expr

    @property
    def psets(self):
        """
        Return a tuple of sets (ProductSets) input of the self.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, Union
        >>> a = Interval(2, 3)
        >>> b = Interval(4, 5)
        >>> c = Interval(1, 7)
        >>> C1 = ComplexRegion(a*b)
        >>> C1.psets
        ([2, 3] x [4, 5],)
        >>> C2 = ComplexRegion(Union(a*b, b*c))
        >>> C2.psets
        ([2, 3] x [4, 5], [4, 5] x [1, 7])

        """
        if self.sets.is_ProductSet:
            psets = ()
            psets = psets + (self.sets, )
        else:
            psets = self.sets.args
        return psets

    @property
    def a_interval(self):
        """
        Return the union of intervals of `x` when, self is in
        rectangular form, or the union of intervals of `r` when
        self is in polar form.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, Union
        >>> a = Interval(2, 3)
        >>> b = Interval(4, 5)
        >>> c = Interval(1, 7)
        >>> C1 = ComplexRegion(a*b)
        >>> C1.a_interval
        [2, 3]
        >>> C2 = ComplexRegion(Union(a*b, b*c))
        >>> C2.a_interval
        [2, 3] U [4, 5]

        """
        a_interval = []
        for element in self.psets:
            a_interval.append(element.args[0])

        a_interval = Union(*a_interval)
        return a_interval

    @property
    def b_interval(self):
        """
        Return the union of intervals of `y` when, self is in
        rectangular form, or the union of intervals of `theta`
        when self is in polar form.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, Union
        >>> a = Interval(2, 3)
        >>> b = Interval(4, 5)
        >>> c = Interval(1, 7)
        >>> C1 = ComplexRegion(a*b)
        >>> C1.b_interval
        [4, 5]
        >>> C2 = ComplexRegion(Union(a*b, b*c))
        >>> C2.b_interval
        [1, 7]

        """
        b_interval = []
        for element in self.psets:
            b_interval.append(element.args[1])

        b_interval = Union(*b_interval)
        return b_interval

    @property
    def polar(self):
        """
        Returns True if self is in polar form.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, Union, S
        >>> a = Interval(2, 3)
        >>> b = Interval(4, 5)
        >>> theta = Interval(0, 2*S.Pi)
        >>> C1 = ComplexRegion(a*b)
        >>> C1.polar
        False
        >>> C2 = ComplexRegion(a*theta, polar=True)
        >>> C2.polar
        True
        """
        return self._polar

    @property
    def _measure(self):
        """
        The measure of self.sets.

        Examples
        ========

        >>> from sympy import Interval, ComplexRegion, S
        >>> a, b = Interval(2, 5), Interval(4, 8)
        >>> c = Interval(0, 2*S.Pi)
        >>> c1 = ComplexRegion(a*b)
        >>> c1.measure
        12
        >>> c2 = ComplexRegion(a*c, polar=True)
        >>> c2.measure
        6*pi

        """
        return self.sets._measure

    def _contains(self, other):
        from sympy.functions import arg, Abs
        from sympy.core.containers import Tuple
        other = sympify(other)
        isTuple = isinstance(other, Tuple)
        if isTuple and len(other) != 2:
            raise ValueError('expecting Tuple of length 2')
        # self in rectangular form
        if not self.polar:
            re, im = other if isTuple else other.as_real_imag()
            for element in self.psets:
                if And(element.args[0]._contains(re),
                        element.args[1]._contains(im)):
                    return True
            return False

        # self in polar form
        elif self.polar:
            if isTuple:
                r, theta = other
            elif other.is_zero:
                r, theta = S.Zero, S.Zero
            else:
                r, theta = Abs(other), arg(other)
            for element in self.psets:
                if And(element.args[0]._contains(r),
                        element.args[1]._contains(theta)):
                    return True
                return False

    def _intersect(self, other):

        if other.is_ComplexRegion:
            # self in rectangular form
            if (not self.polar) and (not other.polar):
                return ComplexRegion(Intersection(self.sets, other.sets))

            # self in polar form
            elif self.polar and other.polar:
                r1, theta1 = self.a_interval, self.b_interval
                r2, theta2 = other.a_interval, other.b_interval
                new_r_interval = Intersection(r1, r2)
                new_theta_interval = Intersection(theta1, theta2)

                # 0 and 2*Pi means the same
                if ((2*S.Pi in theta1 and S.Zero in theta2) or
                   (2*S.Pi in theta2 and S.Zero in theta1)):
                    new_theta_interval = Union(new_theta_interval,
                                               FiniteSet(0))
                return ComplexRegion(new_r_interval*new_theta_interval,
                                    polar=True)

        if other is S.Reals:
            return other

        if other.is_subset(S.Reals):
            new_interval = []

            # self in rectangular form
            if not self.polar:
                for element in self.psets:
                    if S.Zero in element.args[0]:
                        new_interval.append(element.args[0])
                new_interval = Union(*new_interval)
                return Intersection(new_interval, other)

            # self in polar form
            elif self.polar:
                for element in self.psets:
                    if (0 in element.args[1]) or (S.Pi in element.args[1]):
                        new_interval.append(element.args[0])
                new_interval = Union(*new_interval)
                return Intersection(new_interval, other)

    def _union(self, other):

        if other.is_ComplexRegion:

            # self in rectangular form
            if (not self.polar) and (not other.polar):
                return ComplexRegion(Union(self.sets, other.sets))

            # self in polar form
            elif self.polar and other.polar:
                return ComplexRegion(Union(self.sets, other.sets), polar=True)

        if self == S.Complexes:
            return self

        return None


class Complexes(with_metaclass(Singleton, ComplexRegion)):

    def __new__(cls):
        return ComplexRegion.__new__(cls, S.Reals*S.Reals)

    def __eq__(self, other):
        return other == ComplexRegion(S.Reals*S.Reals)

    def __hash__(self):
        return hash(ComplexRegion(S.Reals*S.Reals))

    def __str__(self):
        return "S.Complexes"

    def __repr__(self):
        return "S.Complexes"
