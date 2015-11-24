from sympy import Order, S
from sympy.core.basic import Basic
from sympy.core import Add, Mul, Pow
from sympy.core.expr import AtomicExpr, Expr
from sympy.core.numbers import _sympifyit, oo
from sympy.core.sympify import _sympify
from sympy.sets.sets import (Interval, Intersection, FiniteSet, Union,
                             Complement, EmptySet)


def not_empty_in(finset_intersection, *syms):
    """ Finds the domain of the functions in `finite_set` in which the
    `finite_set` is not-empty

    Parameters
    ==========

    finset_intersection: The unevaluated intersection of FiniteSet containing
                        real-valued functions with Union of Sets
    syms: Tuple of symbols
            Symbol for which domain is to be found

    Raises
    ======

    NotImplementedError
        The algorithms to find the non-emptiness of the given FiniteSet are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report it to the github issue tracker
        (https://github.com/sympy/sympy/issues).

    Examples
    ========

    >>> from sympy import FiniteSet, Interval, not_empty_in, oo
    >>> from sympy.abc import x
    >>> not_empty_in(FiniteSet(x/2).intersect(Interval(0, 1)), x)
    [0, 2]
    >>> not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x)
    [-sqrt(2), -1] U [1, 2]
    >>> not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x)
    (-2, -1] U [2, oo)
    """

    # TODO: handle piecewise defined functions
    # TODO: handle transcendental functions
    # TODO: handle multivariate functions
    if len(syms) == 0:
        raise ValueError("A Symbol or a tuple of symbols must be given "
                         "as the third parameter")

    if finset_intersection.is_EmptySet:
        return EmptySet()

    if isinstance(finset_intersection, Union):
        elm_in_sets = finset_intersection.args[0]
        return Union(not_empty_in(finset_intersection.args[1], *syms),
                     elm_in_sets)

    if isinstance(finset_intersection, FiniteSet):
        finite_set = finset_intersection
        _sets = S.Reals
    else:
        finite_set = finset_intersection.args[1]
        _sets = finset_intersection.args[0]

    if not isinstance(finite_set, FiniteSet):
        raise ValueError('A FiniteSet must be given, not %s: %s' %
                         (type(finite_set), finite_set))

    if len(syms) == 1:
        symb = syms[0]
    else:
        raise NotImplementedError('more than one variables %s not handled' %
                                  (syms,))

    def elm_domain(expr, intrvl):
        """ Finds the domain of an expression in any given interval """
        from sympy.solvers.solveset import solveset

        _start = intrvl.start
        _end = intrvl.end
        _singularities = solveset(expr.as_numer_denom()[1], symb,
                                  domain=S.Reals)

        if intrvl.right_open:
            if _end is S.Infinity:
                _domain1 = S.Reals
            else:
                _domain1 = solveset(expr < _end, symb, domain=S.Reals)
        else:
            _domain1 = solveset(expr <= _end, symb, domain=S.Reals)

        if intrvl.left_open:
            if _start is S.NegativeInfinity:
                _domain2 = S.Reals
            else:
                _domain2 = solveset(expr > _start, symb, domain=S.Reals)
        else:
            _domain2 = solveset(expr >= _start, symb, domain=S.Reals)

        # domain in the interval
        expr_with_sing = Intersection(_domain1, _domain2)
        expr_domain = Complement(expr_with_sing, _singularities)
        return expr_domain

    if isinstance(_sets, Interval):
        return Union(*[elm_domain(element, _sets) for element in finite_set])

    if isinstance(_sets, Union):
        _domain = S.EmptySet
        for intrvl in _sets.args:
            _domain_element = Union(*[elm_domain(element, intrvl)
                                    for element in finite_set])
            _domain = Union(_domain, _domain_element)
        return _domain


class Limits(Interval, AtomicExpr):
    """
    Limits represent an interval `[a, b]`, which is always closed at the
    ends. Here `a` and `b` can be any value from extended real numbers.

    In many cases it would suffice to know that the limit set is bounded.
    However, in some other cases more exact information could be useful.
    For example, that all accumulation values of cos(x) + 1 are non-negative.
    (Limits(-1, 1) + 1 = Limits(0, 2))

    Let `a` and `b` be reals such that a <= b.

    `\langle a, b\rangle = \{x \in \mathbb{R} \mid a \le x \le b\}`

    `\langle -\infty, b\rangle = \{x \in \mathbb{R} \mid x \le b\}`

    `\langle a, \infty \rangle = \{x \in \mathbb{R} \mid, a \le x\}`

    `\langle -\infty, \infty \rangle = \mathbb{R} \cup \{-\infty, \infty\}`

    A Limits object is defined to be real Limits, if its end points are finite
    reals.

    Let `X`, `Y` be real Limits, then their sum, difference, product are
    defined to be the following sets:

    `X + Y = \{ x+y \mid x \in X \cap y \in Y\}`

    `X - Y = \{ x-y \mid x \in X \cap y \in Y\}`

    `X * Y = \{ x*y \mid x \in X \cap y \in Y\}`

    There is, however, no consensus on interval division.

    `X / Y = \{ z| \exists x \in X, y \in Y \mid y \neq 0, z = x/y\}`

    Note: According to this definition the quotient of two Limits may not
    be a Limits object but may be a union of Limits.

    `X^n = \{ x^n \mid x \in X\}`

    Note
    ====

    The main focus in the interval arithmetic is on the simplest way to calculate
    upper and lower endpoints for the range of values of a function in one or more
    variables. These barriers are not necessarily the supremum or infimum, since
    the precise calculation of those values can be difficult or impossible.

    Examples
    ========

    >>> from sympy import Limits, sin, exp, log, pi, E
    >>> from sympy.abc import x

    >>> Limits(0, 1) + Limits(1, 2)
    <1, 3>

    >>> Limits(0, 1) - Limits(0, 2)
    <-2, 1>

    >>> Limits(-2, 3)*Limits(-1, 1)
    <-3, 3>

    >>> Limits(1, 2)*Limits(3, 5)
    <3, 10>

    Note: `<a, b>^2` is not same as `<a, b>*<a, b>`

    >>> Limits(-1, 1)**2
    <0, 1>

    Some elementary functions can also take intervals as input.
    A function `f` evaluated for some interval `<a, b>` is defined as
    `f(<a, b>) = \{ f(x) | a \le x \le b \}`

    >>> sin(Limits(pi/6, pi/3))
    <1/2, sqrt(3)/2>

    >>> exp(Limits(0, 1))
    <1, E>

    >>> log(Limits(1, E))
    <0, 1>

    Some symbol in an experession can be substituted for an Interval.
    But it doesn't necessarily evaluate the Interval for that expression.

    >>> (x**2 + 2*x + 1).subs(x, Limits(-1, 1))
    <-1, 4>

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Interval_arithmetic

    .. [2] http://fab.cba.mit.edu/classes/S62.12/docs/Hickey_interval.pdf

    Notes
    =====

    Do not use ``Limits`` for floating point interval arithmetic calculations
    use ``mpmath.iv`` instead.
    """

    def __new__(cls, min, max):

        min = _sympify(min)
        max = _sympify(max)

        inftys = [S.Infinity, S.NegativeInfinity]
        # Only allow real intervals (use symbols with 'is_real=True').
        if not (min.is_real or min in inftys) \
           or not (max.is_real or max in inftys):
            raise ValueError("Only real intervals are supported")

        # Make sure that the created interval will be valid.
        if max.is_comparable and min.is_comparable:
            if max < min:
                raise ValueError("Lower limit should be smaller than upper limit")

        if max == min:
            return max

        return Basic.__new__(cls, min, max, S.false, S.false)

    # setting the operation priority
    _op_priority = 11.0

    def _intersect(self, other):
        res = Interval(self.min, self.max).intersect(
            Interval(other.min, other.max))
        return Limits(res.min, res.max)

    @property
    def min(self):
        return self.start

    @property
    def max(self):
        return self.end

    @property
    def delta(self):
        return self.max - self.min

    @property
    def mid(self):
        return (self.min + self.max)/2

    @_sympifyit('other', NotImplemented)
    def _eval_power(self, other):
        return self.__pow__(other)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, Limits):
                return Limits(Add(self.min, other.min), Add(self.max, other.max))
            if other is S.Infinity and self.min is S.NegativeInfinity or \
                    other is S.NegativeInfinity and self.max is S.Infinity:
                return Limits(-oo, oo)
            elif other.is_real:
                return Limits(Add(self.min, other), Add(self.max, other))
            return Add(self, other, evaluate=False)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return Limits(-self.max, -self.min)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, Limits):
                return Limits(Add(self.min, -other.max), Add(self.max, -other.min))
            if other is S.NegativeInfinity and self.min is S.NegativeInfinity or \
                    other is S.Infinity and self.max is S.Infinity:
                return Limits(-oo, oo)
            elif other.is_real:
                return Limits(Add(self.min, -other), Add(self.max, - other))
            return Add(self, -other, evaluate=False)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return self.__neg__() + other

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        from sympy.functions.elementary.miscellaneous import Min, Max
        if isinstance(other, Expr):
            if isinstance(other, Limits):
                return Limits(Min(Mul(self.min, other.min),
                            Mul(self.min, other.max),
                            Mul(self.max, other.min),
                            Mul(self.max, other.max)),
                        Max(Mul(self.min, other.min),
                            Mul(self.min, other.max),
                            Mul(self.max, other.min),
                            Mul(self.max, other.max)))
            if other is S.Infinity:
                if self.min.is_zero:
                    return Limits(0, oo)
                if self.max.is_zero:
                    return Limits(-oo, 0)
            if other is S.NegativeInfinity:
                if self.min.is_zero:
                    return Limits(-oo, 0)
                if self.max.is_zero:
                    return Limits(0, oo)
            if other.is_real:
                if other.is_zero:
                    if self == Limits(-oo, oo):
                        return Limits(-oo, oo)
                    if self.max is S.Infinity:
                        return Limits(0, oo)
                    if self.min is S.NegativeInfinity:
                        return Limits(-oo, 0)
                    return S.Zero
                if other.is_positive:
                    return Limits(Mul(self.min, other), Mul(self.max, other))
                elif other.is_negative:
                    return Limits(Mul(self.max, other), Mul(self.min, other))
            if isinstance(other, Order):
                return other
            return Mul(self, other, evaluate=False)
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        from sympy.functions.elementary.miscellaneous import Min, Max
        if isinstance(other, Expr):
            if isinstance(other, Limits):
                if not S.Zero in other:
                    return self*Limits(1/other.max, 1/other.min)

                if S.Zero in self and S.Zero in other:
                    if self.min.is_zero and other.min.is_zero:
                        return Limits(0, oo)
                    if self.max.is_zero and other.min.is_zero:
                        return Limits(-oo, 0)
                    return Limits(-oo, oo)

                if self.max.is_negative:
                    if other.min.is_negative:
                        if other.max.is_zero:
                            return Limits(self.max/other.min, oo)
                        if other.max.is_positive:
                            return Union(Limits(-oo, self.max/other.max),
                                    Limits(self.max/other.min, oo))

                    if other.min.is_zero and other.max.is_positive:
                        return Limits(-oo, self.max/other.max)

                if self.min.is_positive:
                    if other.min.is_negative:
                        if other.max.is_zero:
                            return Limits(-oo, self.min/other.min)
                        if other.max.is_positive:
                            return Union(Limits(-oo, self.min/other.min),
                                    Limits(self.min/other.max, oo))

                    if other.min.is_zero and other.max.is_positive:
                        return Limits(self.min/other.max, oo)

            elif other.is_real:
                if other is S.Infinity or other is S.NegativeInfinity:
                    if self == Limits(-oo, oo):
                        return Limits(-oo, oo)
                    if self.max is S.Infinity:
                        return Limits(Min(0, other), Max(0, other))
                    if self.min is S.NegativeInfinity:
                        return Limits(Min(0, -other), Max(0, -other))
                if other.is_positive:
                    return Limits(self.min/other, self.max/other)
                elif other.is_negative:
                    return Limits(self.max/other, self.min/other)
            return Mul(self, 1/other, evaluate=False)

        return NotImplemented

    __truediv__ = __div__

    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        from sympy.functions.elementary.miscellaneous import Min, Max
        if isinstance(other, Expr):
            if other.is_real:
                if other.is_zero:
                    return S.Zero
                if S.Zero in self:
                    if self.min == S.Zero:
                        if other.is_positive:
                            return Limits(Mul(other, 1/self.max), oo)
                        if other.is_negative:
                            return Limits(-oo, Mul(other, 1/self.max))
                    return Limits(-oo, oo)
                else:
                    return Limits(Min(other/self.min, other/self.max),
                            Max(other/self.min, other/self.max))
            return Mul(other, 1/self, evaluate=False)
        else:
            return NotImplemented

    __rtruediv__ = __rdiv__

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        from sympy.functions.elementary.miscellaneous import Max, Min
        from sympy.functions.elementary.miscellaneous import real_root
        if isinstance(other, Expr):
            if other.is_real:
                if other.is_zero:
                    return S.One
                if other.is_Rational:
                    if other.is_Integer:
                        if self.min.is_positive or self.max.is_negative:
                            return Limits(Min(self.min**other, self.max**other),
                                    Max(self.min**other, self.max**other))
                        else:
                            if self.min.is_zero:
                                return Limits(self.max**other, S.Infinity)
                            if self.max.is_zero:
                                if other % 2 == 0:
                                    return Limits(self.min**other, S.Infinity)
                                return Limits(-oo, self.min**other)
                            if other.is_negative:
                                if other % 2 == 0:
                                    return Limits(S.Zero, oo)
                                return Limits(-oo, oo)
                            if other % 2 == 0:
                                return Limits(S.Zero,
                                    Max(self.min**other, self.max**other))
                            return Limits(self.min**other, self.max**other)

                    num, den = other.as_numer_denom()
                    if num == S(1):
                        if den % 2 == 0:
                            if S.Zero in self:
                                if self.min.is_negative:
                                    return Limits(0, real_root(self.max, den))
                        return Limits(real_root(self.min, den), real_root(self.max, den))
                    num_pow = self**num
                    return num_pow**(1/den)
            return Pow(self, other, evaluate=False)

        return NotImplemented

    def __abs__(self):
        from sympy.functions.elementary.miscellaneous import Max
        if self.max.is_negative:
            return self.__neg__()
        elif self.min.is_negative:
            return Limits(S.Zero, Max(abs(self.min), self.max))
        else:
            return self
