from sympy.functions.elementary.piecewise import Piecewise
from sympy.polys.polytools import lcm_list
from sympy.series.limits import limit
from sympy.series.order import Order
from sympy.core import Add, Mul, Pow, S
from sympy.core.basic import Basic
from sympy.core.expr import AtomicExpr, Expr
from sympy.core.function import diff, expand_mul
from sympy.core.kind import NumberKind
from sympy.core.mod import Mod
from sympy.core.numbers import _sympifyit, oo, zoo
from sympy.core.relational import is_le, is_lt, is_ge, is_gt, Relational
from sympy.core.symbol import Symbol, Dummy
from sympy.core.sympify import _sympify
from sympy.functions.elementary.complexes import Abs, im, re
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import (
    TrigonometricFunction, sin, cos, csc, sec)
from sympy.logic.boolalg import And
from sympy.polys.polytools import degree
from sympy.sets.sets import (Interval, Intersection, FiniteSet, Union,
                             Complement)
from sympy.sets.fancysets import ImageSet
from sympy.simplify.simplify import simplify
from sympy.solvers.decompogen import compogen, decompogen
from sympy.solvers.inequalities import solve_univariate_inequality
from sympy.utilities import filldedent
from sympy.utilities.iterables import iterable
from sympy.multipledispatch import dispatch


def continuous_domain(f, symbol, domain):
    """
    Returns the intervals in the given domain for which the function
    is continuous.
    This method is limited by the ability to determine the various
    singularities and discontinuities of the given function.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for which the intervals are to be determined.
    domain : :py:class:`~.Interval`
        The domain over which the continuity of the symbol has to be checked.

    Examples
    ========

    >>> from sympy import Interval, Symbol, S, tan, log, pi, sqrt
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

    Returns
    =======

    :py:class:`~.Interval`
        Union of all intervals where the function is continuous.

    Raises
    ======

    NotImplementedError
        If the method to determine continuity of such a function
        has not yet been developed.

    """
    from sympy.calculus.singularities import singularities

    if domain.is_subset(S.Reals):
        constrained_interval = domain
        for atom in f.atoms(Pow):
            den = atom.exp.as_numer_denom()[1]
            if den.is_even and den.is_nonzero:
                constraint = solve_univariate_inequality(atom.base >= 0,
                                                         symbol).as_set()
                constrained_interval = Intersection(constraint,
                                                    constrained_interval)

        for atom in f.atoms(log):
            constraint = solve_univariate_inequality(atom.args[0] > 0,
                                                     symbol).as_set()
            constrained_interval = Intersection(constraint,
                                                constrained_interval)


    return constrained_interval - singularities(f, symbol, domain)


def function_range(f, symbol, domain):
    """
    Finds the range of a function in a given domain.
    This method is limited by the ability to determine the singularities and
    determine limits.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for which the range of function is to be determined.
    domain : :py:class:`~.Interval`
        The domain under which the range of the function has to be found.

    Examples
    ========

    >>> from sympy import Interval, Symbol, S, exp, log, pi, sqrt, sin, tan
    >>> from sympy.calculus.util import function_range
    >>> x = Symbol('x')
    >>> function_range(sin(x), x, Interval(0, 2*pi))
    Interval(-1, 1)
    >>> function_range(tan(x), x, Interval(-pi/2, pi/2))
    Interval(-oo, oo)
    >>> function_range(1/x, x, S.Reals)
    Union(Interval.open(-oo, 0), Interval.open(0, oo))
    >>> function_range(exp(x), x, S.Reals)
    Interval.open(0, oo)
    >>> function_range(log(x), x, S.Reals)
    Interval(-oo, oo)
    >>> function_range(sqrt(x), x, Interval(-5, 9))
    Interval(0, 3)

    Returns
    =======

    :py:class:`~.Interval`
        Union of all ranges for all intervals under domain where function is
        continuous.

    Raises
    ======

    NotImplementedError
        If any of the intervals, in the given domain, for which function
        is continuous are not finite or real,
        OR if the critical points of the function on the domain cannot be found.
    """
    from sympy.solvers.solveset import solveset

    if domain is S.EmptySet:
        return S.EmptySet

    period = periodicity(f, symbol)
    if period == S.Zero:
        # the expression is constant wrt symbol
        return FiniteSet(f.expand())

    if period is not None:
        if isinstance(domain, Interval):
            if (domain.inf - domain.sup).is_infinite:
                domain = Interval(0, period)
        elif isinstance(domain, Union):
            for sub_dom in domain.args:
                if isinstance(sub_dom, Interval) and \
                ((sub_dom.inf - sub_dom.sup).is_infinite):
                    domain = Interval(0, period)

    intervals = continuous_domain(f, symbol, domain)
    range_int = S.EmptySet
    if isinstance(intervals,(Interval, FiniteSet)):
        interval_iter = (intervals,)

    elif isinstance(intervals, Union):
        interval_iter = intervals.args

    else:
            raise NotImplementedError(filldedent('''
                Unable to find range for the given domain.
                '''))

    for interval in interval_iter:
        if isinstance(interval, FiniteSet):
            for singleton in interval:
                if singleton in domain:
                    range_int += FiniteSet(f.subs(symbol, singleton))
        elif isinstance(interval, Interval):
            vals = S.EmptySet
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

            solution = solveset(f.diff(symbol), symbol, interval)

            if not iterable(solution):
                raise NotImplementedError(
                        'Unable to find critical points for {}'.format(f))
            if isinstance(solution, ImageSet):
                raise NotImplementedError(
                        'Infinite number of critical points for {}'.format(f))

            critical_points += solution

            for critical_point in critical_points:
                vals += FiniteSet(f.subs(symbol, critical_point))

            left_open, right_open = False, False

            if critical_values is not S.EmptySet:
                if critical_values.inf == vals.inf:
                    left_open = True

                if critical_values.sup == vals.sup:
                    right_open = True

            range_int += Interval(vals.inf, vals.sup, left_open, right_open)
        else:
            raise NotImplementedError(filldedent('''
                Unable to find range for the given domain.
                '''))

    return range_int


def not_empty_in(finset_intersection, *syms):
    """
    Finds the domain of the functions in ``finset_intersection`` in which the
    ``finite_set`` is not-empty

    Parameters
    ==========

    finset_intersection : Intersection of FiniteSet
                        The unevaluated intersection of FiniteSet containing
                        real-valued functions with Union of Sets
    syms : Tuple of symbols
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
    Union(Interval(1, 2), Interval(-sqrt(2), -1))
    >>> not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x)
    Union(Interval.Lopen(-2, -1), Interval(2, oo))
    """

    # TODO: handle piecewise defined functions
    # TODO: handle transcendental functions
    # TODO: handle multivariate functions
    if len(syms) == 0:
        raise ValueError("One or more symbols must be given in syms.")

    if finset_intersection is S.EmptySet:
        return S.EmptySet

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

    f : :py:class:`~.Expr`.
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for which the period is to be determined.
    check : bool, optional
        The flag to verify whether the value being returned is a period or not.

    Returns
    =======

    period
        The period of the function is returned.
        ``None`` is returned when the function is aperiodic or has a complex period.
        The value of $0$ is returned as the period of a constant function.

    Raises
    ======

    NotImplementedError
        The value of the period computed cannot be verified.


    Notes
    =====

    Currently, we do not support functions with a complex period.
    The period of functions having complex periodic values such
    as ``exp``, ``sinh`` is evaluated to ``None``.

    The value returned might not be the "fundamental" period of the given
    function i.e. it may not be the smallest periodic value of the function.

    The verification of the period through the ``check`` flag is not reliable
    due to internal simplification of the given expression. Hence, it is set
    to ``False`` by default.

    Examples
    ========
    >>> from sympy import periodicity, Symbol, sin, cos, tan, exp
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
    if symbol.kind is not NumberKind:
        raise NotImplementedError("Cannot use symbol of kind %s" % symbol.kind)
    temp = Dummy('x', real=True)
    f = f.subs(symbol, temp)
    symbol = temp

    def _check(orig_f, period):
        '''Return the checked period or raise an error.'''
        new_f = orig_f.subs(symbol, symbol + period)
        if new_f.equals(orig_f):
            return period
        else:
            raise NotImplementedError(filldedent('''
                The period of the given function cannot be verified.
                When `%s` was replaced with `%s + %s` in `%s`, the result
                was `%s` which was not recognized as being the same as
                the original function.
                So either the period was wrong or the two forms were
                not recognized as being equal.
                Set check=False to obtain the value.''' %
                (symbol, symbol, period, orig_f, new_f)))

    orig_f = f
    period = None

    if isinstance(f, Relational):
        f = f.lhs - f.rhs

    f = simplify(f)

    if symbol not in f.free_symbols:
        return S.Zero

    if isinstance(f, TrigonometricFunction):
        try:
            period = f.period(symbol)
        except NotImplementedError:
            pass

    if isinstance(f, Abs):
        arg = f.args[0]
        if isinstance(arg, (sec, csc, cos)):
            # all but tan and cot might have a
            # a period that is half as large
            # so recast as sin
            arg = sin(arg.args[0])
        period = periodicity(arg, symbol)
        if period is not None and isinstance(arg, sin):
            # the argument of Abs was a trigonometric other than
            # cot or tan; test to see if the half-period
            # is valid. Abs(arg) has behaviour equivalent to
            # orig_f, so use that for test:
            orig_f = Abs(arg)
            try:
                return _check(orig_f, period/2)
            except NotImplementedError as err:
                if check:
                    raise NotImplementedError(err)
            # else let new orig_f and period be
            # checked below

    if isinstance(f, exp) or (f.is_Pow and f.base == S.Exp1):
        f = Pow(S.Exp1, expand_mul(f.exp))
        if im(f) != 0:
            period_real = periodicity(re(f), symbol)
            period_imag = periodicity(im(f), symbol)
            if period_real is not None and period_imag is not None:
                period = lcim([period_real, period_imag])

    if f.is_Pow and f.base != S.Exp1:
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

    elif isinstance(f, Mod):
        a, n = f.args

        if a == symbol:
            period = n
        elif isinstance(a, TrigonometricFunction):
            period = periodicity(a, symbol)
        #check if 'f' is linear in 'symbol'
        elif (a.is_polynomial(symbol) and degree(a, symbol) == 1 and
            symbol not in n.free_symbols):
                period = Abs(n / a.diff(symbol))

    elif isinstance(f, Piecewise):
        pass  # not handling Piecewise yet as the return type is not favorable

    elif period is None:
        g_s = decompogen(f, symbol)
        num_of_gs = len(g_s)
        if num_of_gs > 1:
            for index, g in enumerate(reversed(g_s)):
                start_index = num_of_gs - 1 - index
                g = compogen(g_s[start_index:], symbol)
                if g not in (orig_f, f): # Fix for issue 12620
                    period = periodicity(g, symbol)
                    if period is not None:
                        break

    if period is not None:
        if check:
            return _check(orig_f, period)
        return period

    return None


def _periodicity(args, symbol):
    """
    Helper for `periodicity` to find the period of a list of simpler
    functions.
    It uses the `lcim` method to find the least common period of
    all the functions.

    Parameters
    ==========

    args : Tuple of :py:class:`~.Symbol`
        All the symbols present in a function.

    symbol : :py:class:`~.Symbol`
        The symbol over which the function is to be evaluated.

    Returns
    =======

    period
        The least common period of the function for all the symbols
        of the function.
        ``None`` if for at least one of the symbols the function is aperiodic.

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

    if periods:
        return periods[0]


def lcim(numbers):
    """Returns the least common integral multiple of a list of numbers.

    The numbers can be rational or irrational or a mixture of both.
    `None` is returned for incommensurable numbers.

    Parameters
    ==========

    numbers : list
        Numbers (rational and/or irrational) for which lcim is to be found.

    Returns
    =======

    number
        lcim if it exists, otherwise ``None`` for incommensurable numbers.

    Examples
    ========

    >>> from sympy.calculus.util import lcim
    >>> from sympy import S, pi
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

def is_convex(f, *syms, domain=S.Reals):
    r"""Determines the  convexity of the function passed in the argument.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    syms : Tuple of :py:class:`~.Symbol`
        The variables with respect to which the convexity is to be determined.
    domain : :py:class:`~.Interval`, optional
        The domain over which the convexity of the function has to be checked.
        If unspecified, S.Reals will be the default domain.

    Returns
    =======

    bool
        The method returns ``True`` if the function is convex otherwise it
        returns ``False``.

    Raises
    ======

    NotImplementedError
        The check for the convexity of multivariate functions is not implemented yet.

    Notes
    =====

    To determine concavity of a function pass `-f` as the concerned function.
    To determine logarithmic convexity of a function pass `\log(f)` as
    concerned function.
    To determine logartihmic concavity of a function pass `-\log(f)` as
    concerned function.

    Currently, convexity check of multivariate functions is not handled.

    Examples
    ========

    >>> from sympy import is_convex, symbols, exp, oo, Interval
    >>> x = symbols('x')
    >>> is_convex(exp(x), x)
    True
    >>> is_convex(x**3, x, domain = Interval(-1, oo))
    False

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Convex_function
    .. [2] http://www.ifp.illinois.edu/~angelia/L3_convfunc.pdf
    .. [3] https://en.wikipedia.org/wiki/Logarithmically_convex_function
    .. [4] https://en.wikipedia.org/wiki/Logarithmically_concave_function
    .. [5] https://en.wikipedia.org/wiki/Concave_function

    """

    if len(syms) > 1:
        raise NotImplementedError(
            "The check for the convexity of multivariate functions is not implemented yet.")

    f = _sympify(f)
    var = syms[0]
    condition = f.diff(var, 2) < 0
    if solve_univariate_inequality(condition, var, False, domain):
        return False
    return True


def stationary_points(f, symbol, domain=S.Reals):
    """
    Returns the stationary points of a function (where derivative of the
    function is 0) in the given domain.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for which the stationary points are to be determined.
    domain : :py:class:`~.Interval`
        The domain over which the stationary points have to be checked.
        If unspecified, ``S.Reals`` will be the default domain.

    Returns
    =======

    Set
        A set of stationary points for the function. If there are no
        stationary point, an :py:class:`~.EmptySet` is returned.

    Examples
    ========

    >>> from sympy import Interval, Symbol, S, sin, pi, pprint, stationary_points
    >>> x = Symbol('x')

    >>> stationary_points(1/x, x, S.Reals)
    EmptySet

    >>> pprint(stationary_points(sin(x), x), use_unicode=False)
              pi                              3*pi
    {2*n*pi + -- | n in Integers} U {2*n*pi + ---- | n in Integers}
              2                                2

    >>> stationary_points(sin(x),x, Interval(0, 4*pi))
    {pi/2, 3*pi/2, 5*pi/2, 7*pi/2}

    """
    from sympy.solvers.solveset import solveset

    if domain is S.EmptySet:
        return S.EmptySet

    domain = continuous_domain(f, symbol, domain)
    set = solveset(diff(f, symbol), symbol, domain)

    return set


def maximum(f, symbol, domain=S.Reals):
    """
    Returns the maximum value of a function in the given domain.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for maximum value needs to be determined.
    domain : :py:class:`~.Interval`
        The domain over which the maximum have to be checked.
        If unspecified, then the global maximum is returned.

    Returns
    =======

    number
        Maximum value of the function in given domain.

    Examples
    ========

    >>> from sympy import Interval, Symbol, S, sin, cos, pi, maximum
    >>> x = Symbol('x')

    >>> f = -x**2 + 2*x + 5
    >>> maximum(f, x, S.Reals)
    6

    >>> maximum(sin(x), x, Interval(-pi, pi/4))
    sqrt(2)/2

    >>> maximum(sin(x)*cos(x), x)
    1/2

    """
    if isinstance(symbol, Symbol):
        if domain is S.EmptySet:
            raise ValueError("Maximum value not defined for empty domain.")

        return function_range(f, symbol, domain).sup
    else:
        raise ValueError("%s is not a valid symbol." % symbol)


def minimum(f, symbol, domain=S.Reals):
    """
    Returns the minimum value of a function in the given domain.

    Parameters
    ==========

    f : :py:class:`~.Expr`
        The concerned function.
    symbol : :py:class:`~.Symbol`
        The variable for minimum value needs to be determined.
    domain : :py:class:`~.Interval`
        The domain over which the minimum have to be checked.
        If unspecified, then the global minimum is returned.

    Returns
    =======

    number
        Minimum value of the function in the given domain.

    Examples
    ========

    >>> from sympy import Interval, Symbol, S, sin, cos, minimum
    >>> x = Symbol('x')

    >>> f = x**2 + 2*x + 5
    >>> minimum(f, x, S.Reals)
    4

    >>> minimum(sin(x), x, Interval(2, 3))
    sin(3)

    >>> minimum(sin(x)*cos(x), x)
    -1/2

    """
    if isinstance(symbol, Symbol):
        if domain is S.EmptySet:
            raise ValueError("Minimum value not defined for empty domain.")

        return function_range(f, symbol, domain).inf
    else:
        raise ValueError("%s is not a valid symbol." % symbol)


class AccumulationBounds(AtomicExpr):
    r"""
    # Note AccumulationBounds has an alias: AccumBounds

    AccumulationBounds represent an interval `[a, b]`, which is always closed
    at the ends. Here `a` and `b` can be any value from extended real numbers.

    The intended meaning of AccummulationBounds is to give an approximate
    location of the accumulation points of a real function at a limit point.

    Let `a` and `b` be reals such that `a \le b`.

    `\left\langle a, b\right\rangle = \{x \in \mathbb{R} \mid a \le x \le b\}`

    `\left\langle -\infty, b\right\rangle = \{x \in \mathbb{R} \mid x \le b\} \cup \{-\infty, \infty\}`

    `\left\langle a, \infty \right\rangle = \{x \in \mathbb{R} \mid a \le x\} \cup \{-\infty, \infty\}`

    `\left\langle -\infty, \infty \right\rangle = \mathbb{R} \cup \{-\infty, \infty\}`

    ``oo`` and ``-oo`` are added to the second and third definition respectively,
    since if either ``-oo`` or ``oo`` is an argument, then the other one should
    be included (though not as an end point). This is forced, since we have,
    for example, ``1/AccumBounds(0, 1) = AccumBounds(1, oo)``, and the limit at
    `0` is not one-sided. As `x` tends to `0-`, then `1/x \rightarrow -\infty`, so `-\infty`
    should be interpreted as belonging to ``AccumBounds(1, oo)`` though it need
    not appear explicitly.

    In many cases it suffices to know that the limit set is bounded.
    However, in some other cases more exact information could be useful.
    For example, all accumulation values of `\cos(x) + 1` are non-negative.
    (``AccumBounds(-1, 1) + 1 = AccumBounds(0, 2)``)

    A AccumulationBounds object is defined to be real AccumulationBounds,
    if its end points are finite reals.

    Let `X`, `Y` be real AccumulationBounds, then their sum, difference,
    product are defined to be the following sets:

    `X + Y = \{ x+y \mid x \in X \cap y \in Y\}`

    `X - Y = \{ x-y \mid x \in X \cap y \in Y\}`

    `X \times Y = \{ x \times y \mid x \in X \cap y \in Y\}`

    When an AccumBounds is raised to a negative power, if 0 is contained
    between the bounds then an infinite range is returned, otherwise if an
    endpoint is 0 then a semi-infinite range with consistent sign will be returned.

    AccumBounds in expressions behave a lot like Intervals but the
    semantics are not necessarily the same. Division (or exponentiation
    to a negative integer power) could be handled with *intervals* by
    returning a union of the results obtained after splitting the
    bounds between negatives and positives, but that is not done with
    AccumBounds. In addition, bounds are assumed to be independent of
    each other; if the same bound is used in more than one place in an
    expression, the result may not be the supremum or infimum of the
    expression (see below). Finally, when a boundary is ``1``,
    exponentiation to the power of ``oo`` yields ``oo``, neither
    ``1`` nor ``nan``.

    Examples
    ========

    >>> from sympy import AccumBounds, sin, exp, log, pi, E, S, oo
    >>> from sympy.abc import x

    >>> AccumBounds(0, 1) + AccumBounds(1, 2)
    AccumBounds(1, 3)

    >>> AccumBounds(0, 1) - AccumBounds(0, 2)
    AccumBounds(-2, 1)

    >>> AccumBounds(-2, 3)*AccumBounds(-1, 1)
    AccumBounds(-3, 3)

    >>> AccumBounds(1, 2)*AccumBounds(3, 5)
    AccumBounds(3, 10)

    The exponentiation of AccumulationBounds is defined
    as follows:

    If 0 does not belong to `X` or `n > 0` then

    `X^n = \{ x^n \mid x \in X\}`

    >>> AccumBounds(1, 4)**(S(1)/2)
    AccumBounds(1, 2)

    otherwise, an infinite or semi-infinite result is obtained:

    >>> 1/AccumBounds(-1, 1)
    AccumBounds(-oo, oo)
    >>> 1/AccumBounds(0, 2)
    AccumBounds(1/2, oo)
    >>> 1/AccumBounds(-oo, 0)
    AccumBounds(-oo, 0)

    A boundary of 1 will always generate all nonnegatives:

    >>> AccumBounds(1, 2)**oo
    AccumBounds(0, oo)
    >>> AccumBounds(0, 1)**oo
    AccumBounds(0, oo)

    If the exponent is itself an AccumulationBounds or is not an
    integer then unevaluated results will be returned unless the base
    values are positive:

    >>> AccumBounds(2, 3)**AccumBounds(-1, 2)
    AccumBounds(1/3, 9)
    >>> AccumBounds(-2, 3)**AccumBounds(-1, 2)
    AccumBounds(-2, 3)**AccumBounds(-1, 2)

    >>> AccumBounds(-2, -1)**(S(1)/2)
    sqrt(AccumBounds(-2, -1))

    Note: `\left\langle a, b\right\rangle^2` is not same as `\left\langle a, b\right\rangle \times \left\langle a, b\right\rangle`

    >>> AccumBounds(-1, 1)**2
    AccumBounds(0, 1)

    >>> AccumBounds(1, 3) < 4
    True

    >>> AccumBounds(1, 3) < -1
    False

    Some elementary functions can also take AccumulationBounds as input.
    A function `f` evaluated for some real AccumulationBounds `\left\langle a, b \right\rangle`
    is defined as `f(\left\langle a, b\right\rangle) = \{ f(x) \mid a \le x \le b \}`

    >>> sin(AccumBounds(pi/6, pi/3))
    AccumBounds(1/2, sqrt(3)/2)

    >>> exp(AccumBounds(0, 1))
    AccumBounds(1, E)

    >>> log(AccumBounds(1, E))
    AccumBounds(0, 1)

    Some symbol in an expression can be substituted for a AccumulationBounds
    object. But it doesn't necessarily evaluate the AccumulationBounds for
    that expression.

    The same expression can be evaluated to different values depending upon
    the form it is used for substitution since each instance of an
    AccumulationBounds is considered independent. For example:

    >>> (x**2 + 2*x + 1).subs(x, AccumBounds(-1, 1))
    AccumBounds(-1, 4)

    >>> ((x + 1)**2).subs(x, AccumBounds(-1, 1))
    AccumBounds(0, 4)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Interval_arithmetic

    .. [2] http://fab.cba.mit.edu/classes/S62.12/docs/Hickey_interval.pdf

    Notes
    =====

    Do not use ``AccumulationBounds`` for floating point interval arithmetic
    calculations, use ``mpmath.iv`` instead.
    """

    is_extended_real = True

    def __new__(cls, min, max):

        min = _sympify(min)
        max = _sympify(max)

        # Only allow real intervals (use symbols with 'is_extended_real=True').
        if not min.is_extended_real or not max.is_extended_real:
            raise ValueError("Only real AccumulationBounds are supported")

        if max == min:
            return max

        # Make sure that the created AccumBounds object will be valid.
        if max.is_number and min.is_number:
            bad = max.is_comparable and min.is_comparable and max < min
        else:
            bad = (max - min).is_extended_negative
        if bad:
            raise ValueError(
                "Lower limit should be smaller than upper limit")

        return Basic.__new__(cls, min, max)

    # setting the operation priority
    _op_priority = 11.0

    def _eval_is_real(self):
        if self.min.is_real and self.max.is_real:
            return True

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
            elif other.is_extended_real:
                if self.min is S.NegativeInfinity and self.max is S.Infinity:
                    return AccumBounds(-oo, oo)
                elif self.min is S.NegativeInfinity:
                    return AccumBounds(-oo, self.max + other)
                elif self.max is S.Infinity:
                    return AccumBounds(self.min + other, oo)
                else:
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
            elif other.is_extended_real:
                if self.min is S.NegativeInfinity and self.max is S.Infinity:
                    return AccumBounds(-oo, oo)
                elif self.min is S.NegativeInfinity:
                    return AccumBounds(-oo, self.max - other)
                elif self.max is S.Infinity:
                    return AccumBounds(self.min - other, oo)
                else:
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
        if self.args == (-oo, oo):
            return self
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                if other.args == (-oo, oo):
                    return other
                v = set()
                for i in self.args:
                    vi = other*i
                    for i in vi.args or (vi,):
                        v.add(i)
                return AccumBounds(Min(*v), Max(*v))
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
            if other.is_extended_real:
                if other.is_zero:
                    if self.max is S.Infinity:
                        return AccumBounds(0, oo)
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(-oo, 0)
                    return S.Zero
                if other.is_extended_positive:
                    return AccumBounds(
                        Mul(self.min, other),
                        Mul(self.max, other))
                elif other.is_extended_negative:
                    return AccumBounds(
                        Mul(self.max, other),
                        Mul(self.min, other))
            if isinstance(other, Order):
                return other
            return Mul(self, other, evaluate=False)
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Expr):
            if isinstance(other, AccumBounds):
                if other.min.is_positive or other.max.is_negative:
                    return self * AccumBounds(1/other.max, 1/other.min)

                if (self.min.is_extended_nonpositive and self.max.is_extended_nonnegative and
                    other.min.is_extended_nonpositive and other.max.is_extended_nonnegative):
                    if self.min.is_zero and other.min.is_zero:
                        return AccumBounds(0, oo)
                    if self.max.is_zero and other.min.is_zero:
                        return AccumBounds(-oo, 0)
                    return AccumBounds(-oo, oo)

                if self.max.is_extended_negative:
                    if other.min.is_extended_negative:
                        if other.max.is_zero:
                            return AccumBounds(self.max / other.min, oo)
                        if other.max.is_extended_positive:
                            # if we were dealing with intervals we would return
                            # Union(Interval(-oo, self.max/other.max),
                            #       Interval(self.max/other.min, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_extended_positive:
                        return AccumBounds(-oo, self.max / other.max)

                if self.min.is_extended_positive:
                    if other.min.is_extended_negative:
                        if other.max.is_zero:
                            return AccumBounds(-oo, self.min / other.min)
                        if other.max.is_extended_positive:
                            # if we were dealing with intervals we would return
                            # Union(Interval(-oo, self.min/other.min),
                            #       Interval(self.min/other.max, oo))
                            return AccumBounds(-oo, oo)

                    if other.min.is_zero and other.max.is_extended_positive:
                        return AccumBounds(self.min / other.max, oo)

            elif other.is_extended_real:
                if other in (S.Infinity, S.NegativeInfinity):
                    if self == AccumBounds(-oo, oo):
                        return AccumBounds(-oo, oo)
                    if self.max is S.Infinity:
                        return AccumBounds(Min(0, other), Max(0, other))
                    if self.min is S.NegativeInfinity:
                        return AccumBounds(Min(0, -other), Max(0, -other))
                if other.is_extended_positive:
                    return AccumBounds(self.min / other, self.max / other)
                elif other.is_extended_negative:
                    return AccumBounds(self.max / other, self.min / other)
            if (1 / other) is S.ComplexInfinity:
                return Mul(self, 1 / other, evaluate=False)
            else:
                return Mul(self, 1 / other)

        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rtruediv__(self, other):
        if isinstance(other, Expr):
            if other.is_extended_real:
                if other.is_zero:
                    return S.Zero
                if (self.min.is_extended_nonpositive and self.max.is_extended_nonnegative):
                    if self.min.is_zero:
                        if other.is_extended_positive:
                            return AccumBounds(Mul(other, 1 / self.max), oo)
                        if other.is_extended_negative:
                            return AccumBounds(-oo, Mul(other, 1 / self.max))
                    if self.max.is_zero:
                        if other.is_extended_positive:
                            return AccumBounds(-oo, Mul(other, 1 / self.min))
                        if other.is_extended_negative:
                            return AccumBounds(Mul(other, 1 / self.min), oo)
                    return AccumBounds(-oo, oo)
                else:
                    return AccumBounds(Min(other / self.min, other / self.max),
                                       Max(other / self.min, other / self.max))
            return Mul(other, 1 / self, evaluate=False)
        else:
            return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        if isinstance(other, Expr):
            if other is S.Infinity:
                if self.min.is_extended_nonnegative:
                    if self.max < 1:
                        return S.Zero
                    if self.min > 1:
                        return S.Infinity
                    return AccumBounds(0, oo)
                elif self.max.is_extended_negative:
                    if self.min > -1:
                        return S.Zero
                    if self.max < -1:
                        return zoo
                    return S.NaN
                else:
                    if self.min > -1:
                        if self.max < 1:
                            return S.Zero
                        return AccumBounds(0, oo)
                    return AccumBounds(-oo, oo)

            if other is S.NegativeInfinity:
                return (1/self)**oo

            # generically true
            if (self.max - self.min).is_nonnegative:
                # well defined
                if self.min.is_nonnegative:
                    # no 0 to worry about
                    if other.is_nonnegative:
                        # no infinity to worry about
                        return self.func(self.min**other, self.max**other)

            if other.is_zero:
                return S.One  # x**0 = 1

            if other.is_Integer or other.is_integer:
                if self.min.is_extended_positive:
                    return AccumBounds(
                        Min(self.min**other, self.max**other),
                        Max(self.min**other, self.max**other))
                elif self.max.is_extended_negative:
                    return AccumBounds(
                        Min(self.max**other, self.min**other),
                        Max(self.max**other, self.min**other))

                if other % 2 == 0:
                    if other.is_extended_negative:
                        if self.min.is_zero:
                            return AccumBounds(self.max**other, oo)
                        if self.max.is_zero:
                            return AccumBounds(self.min**other, oo)
                        return AccumBounds(0, oo)
                    return AccumBounds(
                        S.Zero, Max(self.min**other, self.max**other))
                elif other % 2 == 1:
                    if other.is_extended_negative:
                        if self.min.is_zero:
                            return AccumBounds(self.max**other, oo)
                        if self.max.is_zero:
                            return AccumBounds(-oo, self.min**other)
                        return AccumBounds(-oo, oo)
                    return AccumBounds(self.min**other, self.max**other)

            # non-integer exponent
            # 0**neg or neg**frac yields complex
            if (other.is_number or other.is_rational) and (
                    self.min.is_extended_nonnegative or (
                    other.is_extended_nonnegative and
                    self.min.is_extended_nonnegative)):
                num, den = other.as_numer_denom()
                if num is S.One:
                    return AccumBounds(*[i**(1/den) for i in self.args])

                elif den is not S.One:  # e.g. if other is not Float
                    return (self**num)**(1/den)  # ok for non-negative base

            if isinstance(other, AccumBounds):
                if (self.min.is_extended_positive or
                        self.min.is_extended_nonnegative and
                        other.min.is_extended_nonnegative):
                    p = [self**i for i in other.args]
                    if not any(i.is_Pow for i in p):
                        a = [j for i in p for j in i.args or (i,)]
                        try:
                            return self.func(min(a), max(a))
                        except TypeError:  # can't sort
                            pass

            return Pow(self, other, evaluate=False)

        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        if other.is_real and other.is_extended_nonnegative and (
                self.max - self.min).is_extended_positive:
            if other is S.One:
                return S.One
            if other.is_extended_positive:
                a, b = [other**i for i in self.args]
                if min(a, b) != a:
                    a, b = b, a
                return self.func(a, b)
            if other.is_zero:
                if self.min.is_zero:
                    return self.func(0, 1)
                if self.min.is_extended_positive:
                    return S.Zero

        return Pow(other, self, evaluate=False)

    def __abs__(self):
        if self.max.is_extended_negative:
            return self.__neg__()
        elif self.min.is_extended_negative:
            return AccumBounds(S.Zero, Max(abs(self.min), self.max))
        else:
            return self


    def __contains__(self, other):
        """
        Returns ``True`` if other is contained in self, where other
        belongs to extended real numbers, ``False`` if not contained,
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

        if other in (S.Infinity, S.NegativeInfinity):
            if self.min is S.NegativeInfinity or self.max is S.Infinity:
                return True
            return False

        rv = And(self.min <= other, self.max >= other)
        if rv not in (True, False):
            raise TypeError("input failed to evaluate")
        return rv

    def intersection(self, other):
        """
        Returns the intersection of 'self' and 'other'.
        Here other can be an instance of :py:class:`~.FiniteSet` or AccumulationBounds.

        Parameters
        ==========

        other: AccumulationBounds
             Another AccumulationBounds object with which the intersection
             has to be computed.

        Returns
        =======

        AccumulationBounds
            Intersection of ``self`` and ``other``.

        Examples
        ========

        >>> from sympy import AccumBounds, FiniteSet
        >>> AccumBounds(1, 3).intersection(AccumBounds(2, 4))
        AccumBounds(2, 3)

        >>> AccumBounds(1, 3).intersection(AccumBounds(4, 6))
        EmptySet

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


@dispatch(AccumulationBounds, AccumulationBounds) # type: ignore # noqa:F811
def _eval_is_le(lhs, rhs): # noqa:F811
    if is_le(lhs.max, rhs.min):
        return True
    if is_gt(lhs.min, rhs.max):
        return False


@dispatch(AccumulationBounds, Basic) # type: ignore # noqa:F811
def _eval_is_le(lhs, rhs): # noqa: F811

    """
    Returns ``True `` if range of values attained by ``lhs`` AccumulationBounds
    object is greater than the range of values attained by ``rhs``,
    where ``rhs`` may be any value of type AccumulationBounds object or
    extended real number value, ``False`` if ``rhs`` satisfies
    the same property, else an unevaluated :py:class:`~.Relational`.

    Examples
    ========

    >>> from sympy import AccumBounds, oo
    >>> AccumBounds(1, 3) > AccumBounds(4, oo)
    False
    >>> AccumBounds(1, 4) > AccumBounds(3, 4)
    AccumBounds(1, 4) > AccumBounds(3, 4)
    >>> AccumBounds(1, oo) > -1
    True

    """
    if not rhs.is_extended_real:
            raise TypeError(
                "Invalid comparison of %s %s" %
                (type(rhs), rhs))
    elif rhs.is_comparable:
        if is_le(lhs.max, rhs):
            return True
        if is_gt(lhs.min, rhs):
            return False

@dispatch(AccumulationBounds, AccumulationBounds)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if is_ge(lhs.min, rhs.max):
        return True
    if is_lt(lhs.max, rhs.min):
        return False

@dispatch(AccumulationBounds, Expr)  # type:ignore
def _eval_is_ge(lhs, rhs): # noqa: F811
    """
    Returns ``True`` if range of values attained by ``lhs`` AccumulationBounds
    object is less that the range of values attained by ``rhs``, where
    other may be any value of type AccumulationBounds object or extended
    real number value, ``False`` if ``rhs`` satisfies the same
    property, else an unevaluated :py:class:`~.Relational`.

    Examples
    ========

    >>> from sympy import AccumBounds, oo
    >>> AccumBounds(1, 3) >= AccumBounds(4, oo)
    False
    >>> AccumBounds(1, 4) >= AccumBounds(3, 4)
    AccumBounds(1, 4) >= AccumBounds(3, 4)
    >>> AccumBounds(1, oo) >= 1
    True
    """

    if not rhs.is_extended_real:
        raise TypeError(
            "Invalid comparison of %s %s" %
            (type(rhs), rhs))
    elif rhs.is_comparable:
        if is_ge(lhs.min, rhs):
            return True
        if is_lt(lhs.max, rhs):
            return False


@dispatch(Expr, AccumulationBounds)  # type:ignore
def _eval_is_ge(lhs, rhs): # noqa:F811
    if not lhs.is_extended_real:
        raise TypeError(
            "Invalid comparison of %s %s" %
            (type(lhs), lhs))
    elif lhs.is_comparable:
        if is_le(rhs.max, lhs):
            return True
        if is_gt(rhs.min, lhs):
            return False


@dispatch(AccumulationBounds, AccumulationBounds)  # type:ignore
def _eval_is_ge(lhs, rhs): # noqa:F811
    if is_ge(lhs.min, rhs.max):
        return True
    if is_lt(lhs.max, rhs.min):
        return False

# setting an alias for AccumulationBounds
AccumBounds = AccumulationBounds
