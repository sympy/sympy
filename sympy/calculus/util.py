from sympy import Order, S, log, limit, lcm_list, pi, Abs
from sympy.core.basic import Basic
from sympy.core import Add, Mul, Pow
from sympy.logic.boolalg import And
from sympy.core.expr import AtomicExpr, Expr
from sympy.core.numbers import _sympifyit, oo
from sympy.core.sympify import _sympify
from sympy.sets.sets import (Interval, Intersection, FiniteSet, Union,
                             Complement, EmptySet)
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.utilities import filldedent


def continuous_domain(f, symbol, domain):
    """
    Returns the intervals in the given domain for which the function
    is continuous.
    This method is limited by the ability to determine the various
    singularities and discontinuities of the given function.

    Examples
    ========

    >>> from sympy import Symbol, S, tan, log, pi, sqrt
    >>> from sympy.sets import Interval
    >>> from sympy.calculus.util import continuous_domain
    >>> x = Symbol('x')
    >>> continuous_domain(1/x, x, S.Reals)
    Union(Interval.open(-oo, 0), Interval.open(0, oo))
    >>> continuous_domain(tan(x), x, Interval(0, pi))
    Union(Interval.Ropen(0, pi/2), Interval.Lopen(pi/2, pi))
    >>> continuous_domain(sqrt(x - 2), x, Interval(-5, 5))
    Interval(2, 5)
    >>> continuous_domain(log(2*x - 1), x, S.Reals)
    Interval.open(1/2, oo)

    """
    from sympy.solvers.inequalities import solve_univariate_inequality
    from sympy.solvers.solveset import solveset, _has_rational_power

    if domain.is_subset(S.Reals):
        constrained_interval = domain
        for atom in f.atoms(Pow):
            predicate, denom = _has_rational_power(atom, symbol)
            constraint = S.EmptySet
            if predicate and denom == 2:
                constraint = solve_univariate_inequality(atom.base >= 0,
                                                         symbol).as_set()
                constrained_interval = Intersection(constraint,
                                                    constrained_interval)

        for atom in f.atoms(log):
            constraint = solve_univariate_inequality(atom.args[0] > 0,
                                                     symbol).as_set()
            constrained_interval = Intersection(constraint,
                                                constrained_interval)

        domain = constrained_interval

    try:
        sings = S.EmptySet
        if f.has(Abs):
            sings = solveset(1/f, symbol, domain)
        else:
            for atom in f.atoms(Pow):
                predicate, denom = _has_rational_power(atom, symbol)
                if predicate and denom == 2:
                    sings = solveset(1/f, symbol, domain)
                    break
            else:
                sings = Intersection(solveset(1/f, symbol), domain)

    except BaseException:
        raise NotImplementedError(
            "Methods for determining the continuous domains"
            " of this function have not been developed.")

    return domain - sings


def function_range(f, symbol, domain):
    """
    Finds the range of a function in a given domain.
    This method is limited by the ability to determine the singularities and
    determine limits.

    Examples
    ========

    >>> from sympy import Symbol, S, exp, log, pi, sqrt, sin, tan
    >>> from sympy.sets import Interval
    >>> from sympy.calculus.util import function_range
    >>> x = Symbol('x')
    >>> function_range(sin(x), x, Interval(0, 2*pi))
    Interval(-1, 1)
    >>> function_range(tan(x), x, Interval(-pi/2, pi/2))
    Interval(-oo, oo)
    >>> function_range(1/x, x, S.Reals)
    Interval(-oo, oo)
    >>> function_range(exp(x), x, S.Reals)
    Interval.open(0, oo)
    >>> function_range(log(x), x, S.Reals)
    Interval(-oo, oo)
    >>> function_range(sqrt(x), x , Interval(-5, 9))
    Interval(0, 3)

    """
    from sympy.solvers.solveset import solveset

    vals = S.EmptySet
    period = periodicity(f, symbol)
    if not any(period is i for i in (None, S.Zero)):
        inf = domain.inf
        inf_period = S.Zero if inf.is_infinite else inf
        sup_period = inf_period + period
        periodic_interval = Interval(inf_period, sup_period)
        domain = domain.intersect(periodic_interval)

    intervals = continuous_domain(f, symbol, domain)
    range_int = S.EmptySet
    if isinstance(intervals, Interval):
        interval_iter = (intervals,)

    else:
        interval_iter = intervals.args

    for interval in interval_iter:
        critical_points = S.EmptySet
        critical_values = S.EmptySet
        bounds = ((interval.left_open, interval.inf, '+'),
                  (interval.right_open, interval.sup, '-'))

        for is_open, limit_point, direction in bounds:
            if is_open:
                critical_values += FiniteSet(limit(f, symbol, limit_point, direction))
                vals += critical_values

            else:
                vals += FiniteSet(f.subs(symbol, limit_point))

        critical_points += solveset(f.diff(symbol), symbol, domain)

        for critical_point in critical_points:
            vals += FiniteSet(f.subs(symbol, critical_point))

        left_open, right_open = False, False

        if critical_values is not S.EmptySet:
            if critical_values.inf == vals.inf:
                left_open = True

            if critical_values.sup == vals.sup:
                right_open = True

        range_int += Interval(vals.inf, vals.sup, left_open, right_open)

    return range_int


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
    Interval(0, 2)
    >>> not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x)
    Union(Interval(-sqrt(2), -1), Interval(1, 2))
    >>> not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x)
    Union(Interval.Lopen(-2, -1), Interval(2, oo))
    """

    # TODO: handle piecewise defined functions
    # TODO: handle transcendental functions
    # TODO: handle multivariate functions
    if len(syms) == 0:
        raise ValueError("One or more symbols must be given in syms.")

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


def periodicity(f, symbol, check=False):
    """
    Tests the given function for periodicity in the given symbol.

    Parameters
    ==========

    f : Expr.
        The concerned function.
    symbol : Symbol
        The variable for which the period is to be determined.
    check : Boolean
        The flag to verify whether the value being returned is a period or not.

    Returns
    =======

    period
        The period of the function is returned.
        `None` is returned when the function is aperiodic or has a complex period.
        The value of `0` is returned as the period of a constant function.

    Raises
    ======

    NotImplementedError
        The value of the period computed cannot be verified.


    Notes
    =====

    Currently, we do not support functions with a complex period.
    The period of functions having complex periodic values such
    as `exp`, `sinh` is evaluated to `None`.

    The value returned might not be the "fundamental" period of the given
    function i.e. it may not be the smallest periodic value of the function.

    The verification of the period through the `check` flag is not reliable
    due to internal simplification of the given expression. Hence, it is set
    to `False` by default.

    Examples
    ========
    >>> from sympy import Symbol, sin, cos, tan, exp
    >>> from sympy.calculus.util import periodicity
    >>> x = Symbol('x')
    >>> f = sin(x) + sin(2*x) + sin(3*x)
    >>> periodicity(f, x)
    2*pi
    >>> periodicity(sin(x)*cos(x), x)
    pi
    >>> periodicity(exp(tan(2*x) - 1), x)
    pi/2
    >>> periodicity(sin(4*x)**cos(2*x), x)
    pi
    >>> periodicity(exp(x), x)

    """
    from sympy import simplify, lcm_list
    from sympy.functions.elementary.trigonometric import TrigonometricFunction
    from sympy.solvers.decompogen import decompogen

    orig_f = f
    f = simplify(orig_f)
    period = None

    if not f.has(symbol):
        return S.Zero

    if isinstance(f, TrigonometricFunction):
        try:
            period = f.period(symbol)
        except NotImplementedError:
            pass

    if f.is_Pow:
        base, expo = f.args
        base_has_sym = base.has(symbol)
        expo_has_sym = expo.has(symbol)

        if base_has_sym and not expo_has_sym:
            period = periodicity(base, symbol)

        elif expo_has_sym and not base_has_sym:
            period = periodicity(expo, symbol)

        else:
            period = _periodicity(f.args, symbol)

    elif f.is_Mul:
        coeff, g = f.as_independent(symbol, as_Add=False)
        if isinstance(g, TrigonometricFunction) or coeff is not S.One:
            period = periodicity(g, symbol)

        else:
            period = _periodicity(g.args, symbol)

    elif f.is_Add:
        k, g = f.as_independent(symbol)
        if k is not S.Zero:
            return periodicity(g, symbol)

        period = _periodicity(g.args, symbol)

    elif period is None:
        from sympy.solvers.decompogen import compogen
        g_s = decompogen(f, symbol)
        num_of_gs = len(g_s)
        if num_of_gs > 1:
            for index, g in enumerate(reversed(g_s)):
                start_index = num_of_gs - 1 - index
                g = compogen(g_s[start_index:], symbol)
                if g != orig_f and g != f: # Fix for issue 12620
                    period = periodicity(g, symbol)
                    if period is not None:
                        break

    if period is not None:
        if check:
            if orig_f.subs(symbol, symbol + period) == orig_f:
                return period

            else:
                raise NotImplementedError(filldedent('''
                    The period of the given function cannot be verified.
                    Set check=False to obtain the value.'''))

        return period

    return None


def _periodicity(args, symbol):
    """Helper for periodicity to find the period of a list of simpler
    functions. It uses the `lcim` method to find the least common period of
    all the functions.
    """
    periods = []
    for f in args:
        period = periodicity(f, symbol)
        if period is None:
            return None

        if period is not S.Zero:
            periods.append(period)

    if len(periods) > 1:
        return lcim(periods)

    return periods[0]


def lcim(numbers):
    """Returns the least common integral multiple of a list of numbers.

    The numbers can be rational or irrational or a mixture of both.
    `None` is returned for incommensurable numbers.

    Examples
    ========
    >>> from sympy import S, pi
    >>> from sympy.calculus.util import lcim
    >>> lcim([S(1)/2, S(3)/4, S(5)/6])
    15/2
    >>> lcim([2*pi, 3*pi, pi, pi/2])
    6*pi
    >>> lcim([S(1), 2*pi])
    """
    result = None
    if all(num.is_irrational for num in numbers):
        factorized_nums = list(map(lambda num: num.factor(), numbers))
        factors_num = list(
            map(lambda num: num.as_coeff_Mul(),
                factorized_nums))
        term = factors_num[0][1]
        if all(factor == term for coeff, factor in factors_num):
            common_term = term
            coeffs = [coeff for coeff, factor in factors_num]
            result = lcm_list(coeffs) * common_term

    elif all(num.is_rational for num in numbers):
        result = lcm_list(numbers)

    else:
        pass

    return result


class AccumulationBounds(AtomicExpr):
    r"""
    # Note AccumulationBounds has an alias: AccumBounds

    AccumulationBounds represent an interval `[a, b]`, which is always closed
    at the ends. Here `a` and `b` can be any value from extended real numbers.

    The intended meaning of AccummulationBounds is to give an approximate
    location of the accumulation points of a real function at a limit point.

    Let `a` and `b` be reals such that a <= b.

    `\langle a, b\rangle = \{x \in \mathbb{R} \mid a \le x \le b\}`

    `\langle -\infty, b\rangle = \{x \in \mathbb{R} \mid x \le b\} \cup \{-\infty, \infty\}`

    `\langle a, \infty \rangle = \{x \in \mathbb{R} \mid a \le x\} \cup \{-\infty, \infty\}`

    `\langle -\infty, \infty \rangle = \mathbb{R} \cup \{-\infty, \infty\}`

    `oo` and `-oo` are added to the second and third definition respectively,
    since if either `-oo` or `oo` is an argument, then the other one should
    be included (though not as an end point). This is forced, since we have,
    for example, `1/AccumBounds(0, 1) = AccumBounds(1, oo)`, and the limit at
    `0` is not one-sided. As x tends to `0-`, then `1/x -> -oo`, so `-oo`
    should be interpreted as belonging to `AccumBounds(1, oo)` though it need
    not appear explicitly.

    In many cases it suffices to know that the limit set is bounded.
    However, in some other cases more exact information could be useful.
    For example, all accumulation values of cos(x) + 1 are non-negative.
    (AccumBounds(-1, 1) + 1 = AccumBounds(0, 2))

    A AccumulationBounds object is defined to be real AccumulationBounds,
    if its end points are finite reals.

    Let `X`, `Y` be real AccumulationBounds, then their sum, difference,
    product are defined to be the following sets:

    `X + Y = \{ x+y \mid x \in X \cap y \in Y\}`

    `X - Y = \{ x-y \mid x \in X \cap y \in Y\}`

    `X * Y = \{ x*y \mid x \in X \cap y \in Y\}`

    There is, however, no consensus on Interval division.

    `X / Y = \{ z \mid \exists x \in X, y \in Y \mid y \neq 0, z = x/y\}`

    Note: According to this definition the quotient of two AccumulationBounds
    may not be a AccumulationBounds object but rather a union of
    AccumulationBounds.

    Note
    ====

    The main focus in the interval arithmetic is on the simplest way to
    calculate upper and lower endpoints for the range of values of a
    function in one or more variables. These barriers are not necessarily
    the supremum or infimum, since the precise calculation of those values
    can be difficult or impossible.

    Examples
    ========

    >>> from sympy import AccumBounds, sin, exp, log, pi, E, S, oo
    >>> from sympy.abc import x

    >>> AccumBounds(0, 1) + AccumBounds(1, 2)
    <1, 3>

    >>> AccumBounds(0, 1) - AccumBounds(0, 2)
    <-2, 1>

    >>> AccumBounds(-2, 3)*AccumBounds(-1, 1)
    <-3, 3>

    >>> AccumBounds(1, 2)*AccumBounds(3, 5)
    <3, 10>

    The exponentiation of AccumulationBounds is defined
    as follows:

    If 0 does not belong to `X` or `n > 0` then

    `X^n = \{ x^n \mid x \in X\}`

    otherwise

    `X^n = \{ x^n \mid x \neq 0, x \in X\} \cup \{-\infty, \infty\}`

    Here for fractional `n`, the part of `X` resulting in a complex
    AccumulationBounds object is neglected.

    >>> AccumBounds(-1, 4)**(S(1)/2)
    <0, 2>

    >>> AccumBounds(1, 2)**2
    <1, 4>

    >>> AccumBounds(-1, oo)**(-1)
    <-oo, oo>

    Note: `<a, b>^2` is not same as `<a, b>*<a, b>`

    >>> AccumBounds(-1, 1)**2
    <0, 1>

    >>> AccumBounds(1, 3) < 4
    True

    >>> AccumBounds(1, 3) < -1
    False

    Some elementary functions can also take AccumulationBounds as input.
    A function `f` evaluated for some real AccumulationBounds `<a, b>`
    is defined as `f(\langle a, b\rangle) = \{ f(x) \mid a \le x \le b \}`

    >>> sin(AccumBounds(pi/6, pi/3))
    <1/2, sqrt(3)/2>

    >>> exp(AccumBounds(0, 1))
    <1, E>

    >>> log(AccumBounds(1, E))
    <0, 1>

    Some symbol in an expression can be substituted for a AccumulationBounds
    object. But it doesn't necessarily evaluate the AccumulationBounds for
    that expression.

    Same expression can be evaluated to different values depending upon
    the form it is used for substituion. For example:

    >>> (x**2 + 2*x + 1).subs(x, AccumBounds(-1, 1))
    <-1, 4>

    >>> ((x + 1)**2).subs(x, AccumBounds(-1, 1))
    <0, 4>

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Interval_arithmetic

    .. [2] http://fab.cba.mit.edu/classes/S62.12/docs/Hickey_interval.pdf

    Notes
    =====

    Do not use ``AccumulationBounds`` for floating point interval arithmetic
    calculations, use ``mpmath.iv`` instead.
    """

    is_real = True

    def __new__(cls, min, max):

        min = _sympify(min)
        max = _sympify(max)

        inftys = [S.Infinity, S.NegativeInfinity]
        # Only allow real intervals (use symbols with 'is_real=True').
        if not (min.is_real or min in inftys) \
           or not (max.is_real or max in inftys):
            raise ValueError("Only real AccumulationBounds are supported")

        # Make sure that the created AccumBounds object will be valid.
        if max.is_comparable and min.is_comparable:
            if max < min:
                raise ValueError(
                    "Lower limit should be smaller than upper limit")

        if max == min:
            return max

        return Basic.__new__(cls, min, max)

    # setting the operation priority
    _op_priority = 11.0

    @property
    def min(self):
        """
        Returns the minimum possible value attained by AccumulationBounds
        object.

        Examples
        ========

        >>> from sympy import AccumBounds
        >>> AccumBounds(1, 3).min
        1

        """
        return self.args[0]

    @property
    def max(self):
        """
        Returns the maximum possible value attained by AccumulationBounds
        object.

        Examples
        ========

        >>> from sympy import AccumBounds
        >>> AccumBounds(1, 3).max
        3

        """
        return self.args[1]

    @property
    def delta(self):
        """
        Returns the difference of maximum possible value attained by
        AccumulationBounds object and minimum possible value attained
        by AccumulationBounds object.

        Examples
        ========

        >>> from sympy import AccumBounds
        >>> AccumBounds(1, 3).delta
        2

        """
        return self.max - self.min

    @property
    def mid(self):
        """
        Returns the mean of maximum possible value attained by
        AccumulationBounds object and minimum possible value
        attained by AccumulationBounds object.

        Examples
        ========

        >>> from sympy import AccumBounds
        >>> AccumBounds(1, 3).mid
        2

        """
        return (self.min + self.max) / 2

    @_sympifyit('other', NotImplemented)
    def _eval_power(self, other):
        return self.__pow__(other)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                return AccumBounds(
                    Add(self.min, other.min),
                    Add(self.max, other.max))
            if other is S.Infinity and self.min is S.NegativeInfinity or \
                    other is S.NegativeInfinity and self.max is S.Infinity:
                return AccumBounds(-oo, oo)
            elif other.is_real:
                return AccumBounds(Add(self.min, other), Add(self.max, other))
            return Add(self, other, evaluate=False)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return AccumBounds(-self.max, -self.min)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                return AccumBounds(
                    Add(self.min, -other.max),
                    Add(self.max, -other.min))
            if other is S.NegativeInfinity and self.min is S.NegativeInfinity or \
                    other is S.Infinity and self.max is S.Infinity:
                return AccumBounds(-oo, oo)
            elif other.is_real:
                return AccumBounds(
                    Add(self.min, -other),
                    Add(self.max, -other))
            return Add(self, -other, evaluate=False)
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return self.__neg__() + other

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                return AccumBounds(Min(Mul(self.min, other.min),
                                       Mul(self.min, other.max),
                                       Mul(self.max, other.min),
                                       Mul(self.max, other.max)),
                                   Max(Mul(self.min, other.min),
                                       Mul(self.min, other.max),
                                       Mul(self.max, other.min),
                                       Mul(self.max, other.max)))
            if other is S.Infinity:
                if self.min.is_zero:
                    return AccumBounds(0, oo)
                if self.max.is_zero:
                    return AccumBounds(-oo, 0)
            if other is S.NegativeInfinity:
                if self.min.is_zero:
                    return AccumBounds(-oo, 0)
                if self.max.is_zero:
                    return AccumBounds(0, oo)
            if other.is_real:
                if other.is_zero:
                    if self == AccumBounds(-oo, oo):
                        return AccumBounds(-oo, oo)
                    if self.max is S.Infinity:
                        return AccumBounds(0, oo)
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(-oo, 0)
                    return S.Zero
                if other.is_positive:
                    return AccumBounds(
                        Mul(self.min, other),
                        Mul(self.max, other))
                elif other.is_negative:
                    return AccumBounds(
                        Mul(self.max, other),
                        Mul(self.min, other))
            if isinstance(other, Order):
                return other
            return Mul(self, other, evaluate=False)
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                if S.Zero not in other:
                    return self * AccumBounds(1/other.max, 1/other.min)

                if S.Zero in self and S.Zero in other:
                    if self.min.is_zero and other.min.is_zero:
                        return AccumBounds(0, oo)
                    if self.max.is_zero and other.min.is_zero:
                        return AccumBounds(-oo, 0)
                    return AccumBounds(-oo, oo)

                if self.max.is_negative:
                    if other.min.is_negative:
                        if other.max.is_zero:
                            return AccumBounds(self.max / other.min, oo)
                        if other.max.is_positive:
                            # the actual answer is a Union of AccumBounds,
                            # Union(AccumBounds(-oo, self.max/other.max),
                            #       AccumBounds(self.max/other.min, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_positive:
                        return AccumBounds(-oo, self.max / other.max)

                if self.min.is_positive:
                    if other.min.is_negative:
                        if other.max.is_zero:
                            return AccumBounds(-oo, self.min / other.min)
                        if other.max.is_positive:
                            # the actual answer is a Union of AccumBounds,
                            # Union(AccumBounds(-oo, self.min/other.min),
                            #       AccumBounds(self.min/other.max, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_positive:
                        return AccumBounds(self.min / other.max, oo)

            elif other.is_real:
                if other is S.Infinity or other is S.NegativeInfinity:
                    if self == AccumBounds(-oo, oo):
                        return AccumBounds(-oo, oo)
                    if self.max is S.Infinity:
                        return AccumBounds(Min(0, other), Max(0, other))
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(Min(0, -other), Max(0, -other))
                if other.is_positive:
                    return AccumBounds(self.min / other, self.max / other)
                elif other.is_negative:
                    return AccumBounds(self.max / other, self.min / other)
            return Mul(self, 1 / other, evaluate=False)

        return NotImplemented

    __truediv__ = __div__

    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        if isinstance(other, Expr):
            if other.is_real:
                if other.is_zero:
                    return S.Zero
                if S.Zero in self:
                    if self.min == S.Zero:
                        if other.is_positive:
                            return AccumBounds(Mul(other, 1 / self.max), oo)
                        if other.is_negative:
                            return AccumBounds(-oo, Mul(other, 1 / self.max))
                    if self.max == S.Zero:
                        if other.is_positive:
                            return AccumBounds(-oo, Mul(other, 1 / self.min))
                        if other.is_negative:
                            return AccumBounds(Mul(other, 1 / self.min), oo)
                    return AccumBounds(-oo, oo)
                else:
                    return AccumBounds(Min(other / self.min, other / self.max),
                                       Max(other / self.min, other / self.max))
            return Mul(other, 1 / self, evaluate=False)
        else:
            return NotImplemented

    __rtruediv__ = __rdiv__

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        from sympy.functions.elementary.miscellaneous import real_root
        if isinstance(other, Expr):
            if other is S.Infinity:
                if self.min.is_nonnegative:
                    if self.max < 1:
                        return S.Zero
                    if self.min > 1:
                        return S.Infinity
                    return AccumBounds(0, oo)
                elif self.max.is_negative:
                    if self.min > -1:
                        return S.Zero
                    if self.max < -1:
                        return FiniteSet(-oo, oo)
                    return AccumBounds(-oo, oo)
                else:
                    if self.min > -1:
                        if self.max < 1:
                            return S.Zero
                        return AccumBounds(0, oo)
                    return AccumBounds(-oo, oo)

            if other is S.NegativeInfinity:
                return (1 / self)**oo

            if other.is_real and other.is_number:
                if other.is_zero:
                    return S.One

                if other.is_Integer:
                    if self.min.is_positive:
                        return AccumBounds(
                            Min(self.min ** other, self.max ** other),
                            Max(self.min ** other, self.max ** other))
                    elif self.max.is_negative:
                        return AccumBounds(
                            Min(self.max ** other, self.min ** other),
                            Max(self.max ** other, self.min ** other))

                    if other % 2 == 0:
                        if other.is_negative:
                            if self.min.is_zero:
                                return AccumBounds(self.max**other, oo)
                            if self.max.is_zero:
                                return AccumBounds(self.min**other, oo)
                            return AccumBounds(0, oo)
                        return AccumBounds(
                            S.Zero, Max(self.min**other, self.max**other))
                    else:
                        if other.is_negative:
                            if self.min.is_zero:
                                return AccumBounds(self.max**other, oo)
                            if self.max.is_zero:
                                return AccumBounds(-oo, self.min**other)
                            return AccumBounds(-oo, oo)
                        return AccumBounds(self.min**other, self.max**other)

                num, den = other.as_numer_denom()
                if num == S(1):
                    if den % 2 == 0:
                        if S.Zero in self:
                            if self.min.is_negative:
                                return AccumBounds(0, real_root(self.max, den))
                    return AccumBounds(real_root(self.min, den),
                                       real_root(self.max, den))
                num_pow = self**num
                return num_pow**(1 / den)
            return Pow(self, other, evaluate=False)

        return NotImplemented

    def __abs__(self):
        if self.max.is_negative:
            return self.__neg__()
        elif self.min.is_negative:
            return AccumBounds(S.Zero, Max(abs(self.min), self.max))
        else:
            return self

    def __lt__(self, other):
        """
        Returns True if range of values attained by `self` AccumulationBounds
        object is less than the range of values attained by `other`, where
        other may be any value of type AccumulationBounds object or extended
        real number value, False is returned if `other` satisfies the same
        property,None if the values attained by AccumulationBounds object
        intersect.

        Examples
        ========

        >>> from sympy import AccumBounds, oo
        >>> AccumBounds(1, 3) < AccumBounds(4, oo)
        True
        >>> AccumBounds(1, 4) < AccumBounds(3, 4)

        >>> AccumBounds(1, oo) < -1
        False

        """
        other = _sympify(other)
        if isinstance(other, AccumBounds):
            if self.max < other.min:
                return True
            if self.min >= other.max:
                return False
            return None

        if not(other.is_real or other is S.Infinity or
               other is S.NegativeInfinity):
            raise TypeError(
                "Invalid comparison of %s %s" %
                (type(other), other))

        if other.is_comparable:
            if self.max < other:
                return True
            if self.min >= other:
                return False
        return None

    def __le__(self, other):
        """
        Returns True if range of values attained by `self` AccumulationBounds
        object is less than or equal to the range of values attained by
        `other`, where other may be any value of type AccumulationBounds
        object or extended real number value, AccumulationBounds object,
        False is returned if `other` satisfies the same property, None if
        the values attained by AccumulationBounds object intersect.

        Examples
        ========

        >>> from sympy import AccumBounds, oo
        >>> AccumBounds(1, 3) <= AccumBounds(4, oo)
        True
        >>> AccumBounds(1, 4) <= AccumBounds(3, 4)

        >>> AccumBounds(1, 3) <= 3
        True

        """
        other = _sympify(other)
        if isinstance(other, AccumBounds):
            if self.max <= other.min:
                return True
            if self.min > other.max:
                return False
            return None

        if not(other.is_real or other is S.Infinity or
               other is S.NegativeInfinity):
            raise TypeError(
                "Invalid comparison of %s %s" %
                (type(other), other))

        if other.is_comparable:
            if self.max <= other:
                return True
            if self.min > other:
                return False
        return None

    def __gt__(self, other):
        """
        Returns True if range of values attained by `self` AccumulationBounds
        object is greater than the range of values attained by `other`,
        where other may be any value of type AccumulationBounds object or
        extended real number value, False is returned if `other` satisfies
        the same property, None if the values attained by AccumulationBounds
        object intersect.

        Examples
        ========

        >>> from sympy import AccumBounds, oo
        >>> AccumBounds(1, 3) > AccumBounds(4, oo)
        False
        >>> AccumBounds(1, 4) > AccumBounds(3, 4)

        >>> AccumBounds(1, oo) > -1
        True

        """
        other = _sympify(other)
        if isinstance(other, AccumBounds):
            if self.min > other.max:
                return True
            if self.max <= other.min:
                return False
            return

        if not(other.is_real or other is S.Infinity or
               other is S.NegativeInfinity):
            raise TypeError(
                "Invalid comparison of %s %s" %
                (type(other), other))

        if other.is_comparable:
            if self.min > other:
                return True
            if self.max <= other:
                return False
        return None

    def __ge__(self, other):
        """
        Returns True if range of values attained by `self` AccumulationBounds
        object is less that the range of values attained by `other`, where
        other may be any value of type AccumulationBounds object or extended
        real number value, False is returned if `other` satisfies the same
        property, None if the values attained by AccumulationBounds object
        intersect.

        Examples
        ========

        >>> from sympy import AccumBounds, oo
        >>> AccumBounds(1, 3) >= AccumBounds(4, oo)
        False
        >>> AccumBounds(1, 4) >= AccumBounds(3, 4)

        >>> AccumBounds(1, oo) >= 1
        True

        """
        other = _sympify(other)
        if isinstance(other, AccumBounds):
            if self.min >= other.max:
                return True
            if self.max < other.min:
                return False
            return None

        if not(other.is_real or other is S.Infinity or
               other is S.NegativeInfinity):
            raise TypeError(
                "Invalid comparison of %s %s" %
                (type(other), other))

        if other.is_comparable:
            if self.min >= other:
                return True
            if self.max < other:
                return False
        return None

    def __contains__(self, other):
        """
        Returns True if other is contained in self, where other
        belongs to extended real numbers, False if not contained,
        otherwise TypeError is raised.

        Examples
        ========

        >>> from sympy import AccumBounds, oo
        >>> 1 in AccumBounds(-1, 3)
        True

        -oo and oo go together as limits (in AccumulationBounds).

        >>> -oo in AccumBounds(1, oo)
        True

        >>> oo in AccumBounds(-oo, 0)
        True

        """
        other = _sympify(other)
        if not (other.is_Symbol or other.is_number):
            raise TypeError("Input of type real symbol or Number expected")

        if other is S.Infinity or other is S.NegativeInfinity:
            if self.min is S.NegativeInfinity or self.max is S.Infinity:
                return True
            return False

        return And(self.min <= other and self.max >= other)

    def intersection(self, other):
        """
        Returns the intersection of 'self' and 'other'.
        Here other can be an instance of FiniteSet or AccumulationBounds.

        Examples
        ========

        >>> from sympy import AccumBounds, FiniteSet
        >>> AccumBounds(1, 3).intersection(AccumBounds(2, 4))
        <2, 3>

        >>> AccumBounds(1, 3).intersection(AccumBounds(4, 6))
        EmptySet()

        >>> AccumBounds(1, 4).intersection(FiniteSet(1, 2, 5))
        {1, 2}

        """
        if not isinstance(other, (AccumBounds, FiniteSet)):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")

        if isinstance(other, FiniteSet):
            fin_set = S.EmptySet
            for i in other:
                if i in self:
                    fin_set = fin_set + FiniteSet(i)
            return fin_set

        if self.max < other.min or self.min > other.max:
            return S.EmptySet

        if self.min <= other.min:
            if self.max <= other.max:
                return AccumBounds(other.min, self.max)
            if self.max > other.max:
                return other

        if other.min <= self.min:
            if other.max < self.max:
                return AccumBounds(self.min, other.max)
            if other.max > self.max:
                return self

    def union(self, other):
        # TODO : Devise a better method for Union of AccumBounds
        # this method is not actually correct and
        # can be made better
        if not isinstance(other, AccumBounds):
            raise TypeError(
                "Input must be AccumulationBounds or FiniteSet object")

        if self.min <= other.min and self.max >= other.min:
            return AccumBounds(self.min, Max(self.max, other.max))

        if other.min <= self.min and other.max >= self.min:
            return AccumBounds(other.min, Max(self.max, other.max))


# setting an alias for AccumulationBounds
AccumBounds = AccumulationBounds
