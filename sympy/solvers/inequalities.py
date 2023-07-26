"""Tools for solving inequalities and systems of inequalities. """
import itertools

from sympy.calculus.util import (continuous_domain, periodicity,
    function_range)
from sympy.core import Symbol, Dummy, sympify
from sympy.core.containers import Dict
from sympy.core.exprtools import factor_terms
from sympy.core.relational import Relational, Lt, Ge, Eq, Le
from sympy.assumptions.ask import Q
from sympy.assumptions.relation.binrel import AppliedBinaryRelation
from sympy.sets.sets import Interval, FiniteSet, Union, Intersection
from sympy.core.singleton import S
from sympy.core.sorting import ordered
from sympy.core.function import expand_mul
from sympy.functions.elementary.complexes import im, Abs
from sympy.logic import And
from sympy.polys import Poly, PolynomialError, parallel_poly_from_expr
from sympy.polys.polyutils import _nsort
from sympy.solvers.solveset import solvify, solveset, linear_eq_to_matrix
from sympy.utilities.iterables import sift, iterable, numbered_symbols
from sympy.utilities.misc import filldedent


def solve_poly_inequality(poly, rel):
    """Solve a polynomial inequality with rational coefficients.

    Examples
    ========

    >>> from sympy import solve_poly_inequality, Poly
    >>> from sympy.abc import x

    >>> solve_poly_inequality(Poly(x, x, domain='ZZ'), '==')
    [{0}]

    >>> solve_poly_inequality(Poly(x**2 - 1, x, domain='ZZ'), '!=')
    [Interval.open(-oo, -1), Interval.open(-1, 1), Interval.open(1, oo)]

    >>> solve_poly_inequality(Poly(x**2 - 1, x, domain='ZZ'), '==')
    [{-1}, {1}]

    See Also
    ========
    solve_poly_inequalities
    """
    if not isinstance(poly, Poly):
        raise ValueError(
            'For efficiency reasons, `poly` should be a Poly instance')
    if poly.as_expr().is_number:
        t = Relational(poly.as_expr(), 0, rel)
        if t is S.true:
            return [S.Reals]
        elif t is S.false:
            return [S.EmptySet]
        else:
            raise NotImplementedError(
                "could not determine truth value of %s" % t)

    reals, intervals = poly.real_roots(multiple=False), []

    if rel == '==':
        for root, _ in reals:
            interval = Interval(root, root)
            intervals.append(interval)
    elif rel == '!=':
        left = S.NegativeInfinity

        for right, _ in reals + [(S.Infinity, 1)]:
            interval = Interval(left, right, True, True)
            intervals.append(interval)
            left = right
    else:
        if poly.LC() > 0:
            sign = +1
        else:
            sign = -1

        eq_sign, equal = None, False

        if rel == '>':
            eq_sign = +1
        elif rel == '<':
            eq_sign = -1
        elif rel == '>=':
            eq_sign, equal = +1, True
        elif rel == '<=':
            eq_sign, equal = -1, True
        else:
            raise ValueError("'%s' is not a valid relation" % rel)

        right, right_open = S.Infinity, True

        for left, multiplicity in reversed(reals):
            if multiplicity % 2:
                if sign == eq_sign:
                    intervals.insert(
                        0, Interval(left, right, not equal, right_open))

                sign, right, right_open = -sign, left, not equal
            else:
                if sign == eq_sign and not equal:
                    intervals.insert(
                        0, Interval(left, right, True, right_open))
                    right, right_open = left, True
                elif sign != eq_sign and equal:
                    intervals.insert(0, Interval(left, left))

        if sign == eq_sign:
            intervals.insert(
                0, Interval(S.NegativeInfinity, right, True, right_open))

    return intervals


def solve_poly_inequalities(polys):
    """Solve polynomial inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy import Poly
    >>> from sympy.solvers.inequalities import solve_poly_inequalities
    >>> from sympy.abc import x
    >>> solve_poly_inequalities(((
    ... Poly(x**2 - 3), ">"), (
    ... Poly(-x**2 + 1), ">")))
    Union(Interval.open(-oo, -sqrt(3)), Interval.open(-1, 1), Interval.open(sqrt(3), oo))
    """
    return Union(*[s for p in polys for s in solve_poly_inequality(*p)])


def solve_rational_inequalities(eqs):
    """Solve a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy import solve_rational_inequalities, Poly

    >>> solve_rational_inequalities([[
    ... ((Poly(-x + 1), Poly(1, x)), '>='),
    ... ((Poly(-x + 1), Poly(1, x)), '<=')]])
    {1}

    >>> solve_rational_inequalities([[
    ... ((Poly(x), Poly(1, x)), '!='),
    ... ((Poly(-x + 1), Poly(1, x)), '>=')]])
    Union(Interval.open(-oo, 0), Interval.Lopen(0, 1))

    See Also
    ========
    solve_poly_inequality
    """
    result = S.EmptySet

    for _eqs in eqs:
        if not _eqs:
            continue

        global_intervals = [Interval(S.NegativeInfinity, S.Infinity)]

        for (numer, denom), rel in _eqs:
            numer_intervals = solve_poly_inequality(numer*denom, rel)
            denom_intervals = solve_poly_inequality(denom, '==')

            intervals = []

            for numer_interval, global_interval in itertools.product(
                    numer_intervals, global_intervals):
                interval = numer_interval.intersect(global_interval)

                if interval is not S.EmptySet:
                    intervals.append(interval)

            global_intervals = intervals

            intervals = []

            for global_interval in global_intervals:
                for denom_interval in denom_intervals:
                    global_interval -= denom_interval

                if global_interval is not S.EmptySet:
                    intervals.append(global_interval)

            global_intervals = intervals

            if not global_intervals:
                break

        for interval in global_intervals:
            result = result.union(interval)

    return result


def reduce_rational_inequalities(exprs, gen, relational=True):
    """Reduce a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.solvers.inequalities import reduce_rational_inequalities

    >>> x = Symbol('x', real=True)

    >>> reduce_rational_inequalities([[x**2 <= 0]], x)
    Eq(x, 0)

    >>> reduce_rational_inequalities([[x + 2 > 0]], x)
    -2 < x
    >>> reduce_rational_inequalities([[(x + 2, ">")]], x)
    -2 < x
    >>> reduce_rational_inequalities([[x + 2]], x)
    Eq(x, -2)

    This function find the non-infinite solution set so if the unknown symbol
    is declared as extended real rather than real then the result may include
    finiteness conditions:

    >>> y = Symbol('y', extended_real=True)
    >>> reduce_rational_inequalities([[y + 2 > 0]], y)
    (-2 < y) & (y < oo)
    """
    exact = True
    eqs = []
    solution = S.Reals if exprs else S.EmptySet
    for _exprs in exprs:
        _eqs = []

        for expr in _exprs:
            if isinstance(expr, tuple):
                expr, rel = expr
            else:
                if expr.is_Relational:
                    expr, rel = expr.lhs - expr.rhs, expr.rel_op
                else:
                    expr, rel = expr, '=='

            if expr is S.true:
                numer, denom, rel = S.Zero, S.One, '=='
            elif expr is S.false:
                numer, denom, rel = S.One, S.One, '=='
            else:
                numer, denom = expr.together().as_numer_denom()

            try:
                (numer, denom), opt = parallel_poly_from_expr(
                    (numer, denom), gen)
            except PolynomialError:
                raise PolynomialError(filldedent('''
                    only polynomials and rational functions are
                    supported in this context.
                    '''))

            if not opt.domain.is_Exact:
                numer, denom, exact = numer.to_exact(), denom.to_exact(), False

            domain = opt.domain.get_exact()

            if not (domain.is_ZZ or domain.is_QQ):
                expr = numer/denom
                expr = Relational(expr, 0, rel)
                solution &= solve_univariate_inequality(expr, gen, relational=False)
            else:
                _eqs.append(((numer, denom), rel))

        if _eqs:
            eqs.append(_eqs)

    if eqs:
        solution &= solve_rational_inequalities(eqs)
        exclude = solve_rational_inequalities([[((d, d.one), '==')
            for i in eqs for ((n, d), _) in i if d.has(gen)]])
        solution -= exclude

    if not exact and solution:
        solution = solution.evalf()

    if relational:
        solution = solution.as_relational(gen)

    return solution


def reduce_abs_inequality(expr, rel, gen):
    """Reduce an inequality with nested absolute values.

    Examples
    ========

    >>> from sympy import reduce_abs_inequality, Abs, Symbol
    >>> x = Symbol('x', real=True)

    >>> reduce_abs_inequality(Abs(x - 5) - 3, '<', x)
    (2 < x) & (x < 8)

    >>> reduce_abs_inequality(Abs(x + 2)*3 - 13, '<', x)
    (-19/3 < x) & (x < 7/3)

    See Also
    ========

    reduce_abs_inequalities
    """
    if gen.is_extended_real is False:
         raise TypeError(filldedent('''
            Cannot solve inequalities with absolute values containing
            non-real variables.
            '''))

    def _bottom_up_scan(expr):
        exprs = []

        if expr.is_Add or expr.is_Mul:
            op = expr.func

            for arg in expr.args:
                _exprs = _bottom_up_scan(arg)

                if not exprs:
                    exprs = _exprs
                else:
                    exprs = [(op(expr, _expr), conds + _conds) for (expr, conds), (_expr, _conds) in
                            itertools.product(exprs, _exprs)]
        elif expr.is_Pow:
            n = expr.exp
            if not n.is_Integer:
                raise ValueError("Only Integer Powers are allowed on Abs.")

            exprs.extend((expr**n, conds) for expr, conds in _bottom_up_scan(expr.base))
        elif isinstance(expr, Abs):
            _exprs = _bottom_up_scan(expr.args[0])

            for expr, conds in _exprs:
                exprs.append(( expr, conds + [Ge(expr, 0)]))
                exprs.append((-expr, conds + [Lt(expr, 0)]))
        else:
            exprs = [(expr, [])]

        return exprs

    mapping = {'<': '>', '<=': '>='}
    inequalities = []

    for expr, conds in _bottom_up_scan(expr):
        if rel not in mapping.keys():
            expr = Relational( expr, 0, rel)
        else:
            expr = Relational(-expr, 0, mapping[rel])

        inequalities.append([expr] + conds)

    return reduce_rational_inequalities(inequalities, gen)


def reduce_abs_inequalities(exprs, gen):
    """Reduce a system of inequalities with nested absolute values.

    Examples
    ========

    >>> from sympy import reduce_abs_inequalities, Abs, Symbol
    >>> x = Symbol('x', extended_real=True)

    >>> reduce_abs_inequalities([(Abs(3*x - 5) - 7, '<'),
    ... (Abs(x + 25) - 13, '>')], x)
    (-2/3 < x) & (x < 4) & (((-oo < x) & (x < -38)) | ((-12 < x) & (x < oo)))

    >>> reduce_abs_inequalities([(Abs(x - 4) + Abs(3*x - 5) - 7, '<')], x)
    (1/2 < x) & (x < 4)

    See Also
    ========

    reduce_abs_inequality
    """
    return And(*[ reduce_abs_inequality(expr, rel, gen)
        for expr, rel in exprs ])


def solve_univariate_inequality(expr, gen, relational=True, domain=S.Reals, continuous=False):
    """Solves a real univariate inequality.

    Parameters
    ==========

    expr : Relational
        The target inequality
    gen : Symbol
        The variable for which the inequality is solved
    relational : bool
        A Relational type output is expected or not
    domain : Set
        The domain over which the equation is solved
    continuous: bool
        True if expr is known to be continuous over the given domain
        (and so continuous_domain() does not need to be called on it)

    Raises
    ======

    NotImplementedError
        The solution of the inequality cannot be determined due to limitation
        in :func:`sympy.solvers.solveset.solvify`.

    Notes
    =====

    Currently, we cannot solve all the inequalities due to limitations in
    :func:`sympy.solvers.solveset.solvify`. Also, the solution returned for trigonometric inequalities
    are restricted in its periodic interval.

    See Also
    ========

    sympy.solvers.solveset.solvify: solver returning solveset solutions with solve's output API

    Examples
    ========

    >>> from sympy import solve_univariate_inequality, Symbol, sin, Interval, S
    >>> x = Symbol('x')

    >>> solve_univariate_inequality(x**2 >= 4, x)
    ((2 <= x) & (x < oo)) | ((-oo < x) & (x <= -2))

    >>> solve_univariate_inequality(x**2 >= 4, x, relational=False)
    Union(Interval(-oo, -2), Interval(2, oo))

    >>> domain = Interval(0, S.Infinity)
    >>> solve_univariate_inequality(x**2 >= 4, x, False, domain)
    Interval(2, oo)

    >>> solve_univariate_inequality(sin(x) > 0, x, relational=False)
    Interval.open(0, pi)

    """
    from sympy.solvers.solvers import denoms

    if domain.is_subset(S.Reals) is False:
        raise NotImplementedError(filldedent('''
        Inequalities in the complex domain are
        not supported. Try the real domain by
        setting domain=S.Reals'''))
    elif domain is not S.Reals:
        rv = solve_univariate_inequality(
        expr, gen, relational=False, continuous=continuous).intersection(domain)
        if relational:
            rv = rv.as_relational(gen)
        return rv
    else:
        pass  # continue with attempt to solve in Real domain

    # This keeps the function independent of the assumptions about `gen`.
    # `solveset` makes sure this function is called only when the domain is
    # real.
    _gen = gen
    _domain = domain
    if gen.is_extended_real is False:
        rv = S.EmptySet
        return rv if not relational else rv.as_relational(_gen)
    elif gen.is_extended_real is None:
        gen = Dummy('gen', extended_real=True)
        try:
            expr = expr.xreplace({_gen: gen})
        except TypeError:
            raise TypeError(filldedent('''
                When gen is real, the relational has a complex part
                which leads to an invalid comparison like I < 0.
                '''))

    rv = None

    if expr is S.true:
        rv = domain

    elif expr is S.false:
        rv = S.EmptySet

    else:
        e = expr.lhs - expr.rhs
        period = periodicity(e, gen)
        if period == S.Zero:
            e = expand_mul(e)
            const = expr.func(e, 0)
            if const is S.true:
                rv = domain
            elif const is S.false:
                rv = S.EmptySet
        elif period is not None:
            frange = function_range(e, gen, domain)

            rel = expr.rel_op
            if rel in ('<', '<='):
                if expr.func(frange.sup, 0):
                    rv = domain
                elif not expr.func(frange.inf, 0):
                    rv = S.EmptySet

            elif rel in ('>', '>='):
                if expr.func(frange.inf, 0):
                    rv = domain
                elif not expr.func(frange.sup, 0):
                    rv = S.EmptySet

            inf, sup = domain.inf, domain.sup
            if sup - inf is S.Infinity:
                domain = Interval(0, period, False, True).intersect(_domain)
                _domain = domain

        if rv is None:
            n, d = e.as_numer_denom()
            try:
                if gen not in n.free_symbols and len(e.free_symbols) > 1:
                    raise ValueError
                # this might raise ValueError on its own
                # or it might give None...
                solns = solvify(e, gen, domain)
                if solns is None:
                    # in which case we raise ValueError
                    raise ValueError
            except (ValueError, NotImplementedError):
                # replace gen with generic x since it's
                # univariate anyway
                raise NotImplementedError(filldedent('''
                    The inequality, %s, cannot be solved using
                    solve_univariate_inequality.
                    ''' % expr.subs(gen, Symbol('x'))))

            expanded_e = expand_mul(e)
            def valid(x):
                # this is used to see if gen=x satisfies the
                # relational by substituting it into the
                # expanded form and testing against 0, e.g.
                # if expr = x*(x + 1) < 2 then e = x*(x + 1) - 2
                # and expanded_e = x**2 + x - 2; the test is
                # whether a given value of x satisfies
                # x**2 + x - 2 < 0
                #
                # expanded_e, expr and gen used from enclosing scope
                v = expanded_e.subs(gen, expand_mul(x))
                try:
                    r = expr.func(v, 0)
                except TypeError:
                    r = S.false
                if r in (S.true, S.false):
                    return r
                if v.is_extended_real is False:
                    return S.false
                else:
                    v = v.n(2)
                    if v.is_comparable:
                        return expr.func(v, 0)
                    # not comparable or couldn't be evaluated
                    raise NotImplementedError(
                        'relationship did not evaluate: %s' % r)

            singularities = []
            for d in denoms(expr, gen):
                singularities.extend(solvify(d, gen, domain))
            if not continuous:
                domain = continuous_domain(expanded_e, gen, domain)

            include_x = '=' in expr.rel_op and expr.rel_op != '!='

            try:
                discontinuities = set(domain.boundary -
                    FiniteSet(domain.inf, domain.sup))
                # remove points that are not between inf and sup of domain
                critical_points = FiniteSet(*(solns + singularities + list(
                    discontinuities))).intersection(
                    Interval(domain.inf, domain.sup,
                    domain.inf not in domain, domain.sup not in domain))
                if all(r.is_number for r in critical_points):
                    reals = _nsort(critical_points, separated=True)[0]
                else:
                    sifted = sift(critical_points, lambda x: x.is_extended_real)
                    if sifted[None]:
                        # there were some roots that weren't known
                        # to be real
                        raise NotImplementedError
                    try:
                        reals = sifted[True]
                        if len(reals) > 1:
                            reals = sorted(reals)
                    except TypeError:
                        raise NotImplementedError
            except NotImplementedError:
                raise NotImplementedError('sorting of these roots is not supported')

            # If expr contains imaginary coefficients, only take real
            # values of x for which the imaginary part is 0
            make_real = S.Reals
            if im(expanded_e) != S.Zero:
                check = True
                im_sol = FiniteSet()
                try:
                    a = solveset(im(expanded_e), gen, domain)
                    if not isinstance(a, Interval):
                        for z in a:
                            if z not in singularities and valid(z) and z.is_extended_real:
                                im_sol += FiniteSet(z)
                    else:
                        start, end = a.inf, a.sup
                        for z in _nsort(critical_points + FiniteSet(end)):
                            valid_start = valid(start)
                            if start != end:
                                valid_z = valid(z)
                                pt = _pt(start, z)
                                if pt not in singularities and pt.is_extended_real and valid(pt):
                                    if valid_start and valid_z:
                                        im_sol += Interval(start, z)
                                    elif valid_start:
                                        im_sol += Interval.Ropen(start, z)
                                    elif valid_z:
                                        im_sol += Interval.Lopen(start, z)
                                    else:
                                        im_sol += Interval.open(start, z)
                            start = z
                        for s in singularities:
                            im_sol -= FiniteSet(s)
                except (TypeError):
                    im_sol = S.Reals
                    check = False

                if im_sol is S.EmptySet:
                    raise ValueError(filldedent('''
                        %s contains imaginary parts which cannot be
                        made 0 for any value of %s satisfying the
                        inequality, leading to relations like I < 0.
                        '''  % (expr.subs(gen, _gen), _gen)))

                make_real = make_real.intersect(im_sol)

            sol_sets = [S.EmptySet]

            start = domain.inf
            if start in domain and valid(start) and start.is_finite:
                sol_sets.append(FiniteSet(start))

            for x in reals:
                end = x

                if valid(_pt(start, end)):
                    sol_sets.append(Interval(start, end, True, True))

                if x in singularities:
                    singularities.remove(x)
                else:
                    if x in discontinuities:
                        discontinuities.remove(x)
                        _valid = valid(x)
                    else:  # it's a solution
                        _valid = include_x
                    if _valid:
                        sol_sets.append(FiniteSet(x))

                start = end

            end = domain.sup
            if end in domain and valid(end) and end.is_finite:
                sol_sets.append(FiniteSet(end))

            if valid(_pt(start, end)):
                sol_sets.append(Interval.open(start, end))

            if im(expanded_e) != S.Zero and check:
                rv = (make_real).intersect(_domain)
            else:
                rv = Intersection(
                    (Union(*sol_sets)), make_real, _domain).subs(gen, _gen)

    return rv if not relational else rv.as_relational(_gen)


def _pt(start, end):
    """Return a point between start and end"""
    if not start.is_infinite and not end.is_infinite:
        pt = (start + end)/2
    elif start.is_infinite and end.is_infinite:
        pt = S.Zero
    else:
        if (start.is_infinite and start.is_extended_positive is None or
                end.is_infinite and end.is_extended_positive is None):
            raise ValueError('cannot proceed with unsigned infinite values')
        if (end.is_infinite and end.is_extended_negative or
                start.is_infinite and start.is_extended_positive):
            start, end = end, start
        # if possible, use a multiple of self which has
        # better behavior when checking assumptions than
        # an expression obtained by adding or subtracting 1
        if end.is_infinite:
            if start.is_extended_positive:
                pt = start*2
            elif start.is_extended_negative:
                pt = start*S.Half
            else:
                pt = start + 1
        elif start.is_infinite:
            if end.is_extended_positive:
                pt = end*S.Half
            elif end.is_extended_negative:
                pt = end*2
            else:
                pt = end - 1
    return pt


def _solve_inequality(ie, s, linear=False):
    """Return the inequality with s isolated on the left, if possible.
    If the relationship is non-linear, a solution involving And or Or
    may be returned. False or True are returned if the relationship
    is never True or always True, respectively.

    If `linear` is True (default is False) an `s`-dependent expression
    will be isolated on the left, if possible
    but it will not be solved for `s` unless the expression is linear
    in `s`. Furthermore, only "safe" operations which do not change the
    sense of the relationship are applied: no division by an unsigned
    value is attempted unless the relationship involves Eq or Ne and
    no division by a value not known to be nonzero is ever attempted.

    Examples
    ========

    >>> from sympy import Eq, Symbol
    >>> from sympy.solvers.inequalities import _solve_inequality as f
    >>> from sympy.abc import x, y

    For linear expressions, the symbol can be isolated:

    >>> f(x - 2 < 0, x)
    x < 2
    >>> f(-x - 6 < x, x)
    x > -3

    Sometimes nonlinear relationships will be False

    >>> f(x**2 + 4 < 0, x)
    False

    Or they may involve more than one region of values:

    >>> f(x**2 - 4 < 0, x)
    (-2 < x) & (x < 2)

    To restrict the solution to a relational, set linear=True
    and only the x-dependent portion will be isolated on the left:

    >>> f(x**2 - 4 < 0, x, linear=True)
    x**2 < 4

    Division of only nonzero quantities is allowed, so x cannot
    be isolated by dividing by y:

    >>> y.is_nonzero is None  # it is unknown whether it is 0 or not
    True
    >>> f(x*y < 1, x)
    x*y < 1

    And while an equality (or inequality) still holds after dividing by a
    non-zero quantity

    >>> nz = Symbol('nz', nonzero=True)
    >>> f(Eq(x*nz, 1), x)
    Eq(x, 1/nz)

    the sign must be known for other inequalities involving > or <:

    >>> f(x*nz <= 1, x)
    nz*x <= 1
    >>> p = Symbol('p', positive=True)
    >>> f(x*p <= 1, x)
    x <= 1/p

    When there are denominators in the original expression that
    are removed by expansion, conditions for them will be returned
    as part of the result:

    >>> f(x < x*(2/x - 1), x)
    (x < 1) & Ne(x, 0)
    """
    from sympy.solvers.solvers import denoms
    if s not in ie.free_symbols:
        return ie
    if ie.rhs == s:
        ie = ie.reversed
    if ie.lhs == s and s not in ie.rhs.free_symbols:
        return ie

    def classify(ie, s, i):
        # return True or False if ie evaluates when substituting s with
        # i else None (if unevaluated) or NaN (when there is an error
        # in evaluating)
        try:
            v = ie.subs(s, i)
            if v is S.NaN:
                return v
            elif v not in (True, False):
                return
            return v
        except TypeError:
            return S.NaN

    rv = None
    oo = S.Infinity
    expr = ie.lhs - ie.rhs
    try:
        p = Poly(expr, s)
        if p.degree() == 0:
            rv = ie.func(p.as_expr(), 0)
        elif not linear and p.degree() > 1:
            # handle in except clause
            raise NotImplementedError
    except (PolynomialError, NotImplementedError):
        if not linear:
            try:
                rv = reduce_rational_inequalities([[ie]], s)
            except PolynomialError:
                rv = solve_univariate_inequality(ie, s)
            # remove restrictions wrt +/-oo that may have been
            # applied when using sets to simplify the relationship
            okoo = classify(ie, s, oo)
            if okoo is S.true and classify(rv, s, oo) is S.false:
                rv = rv.subs(s < oo, True)
            oknoo = classify(ie, s, -oo)
            if (oknoo is S.true and
                    classify(rv, s, -oo) is S.false):
                rv = rv.subs(-oo < s, True)
                rv = rv.subs(s > -oo, True)
            if rv is S.true:
                rv = (s <= oo) if okoo is S.true else (s < oo)
                if oknoo is not S.true:
                    rv = And(-oo < s, rv)
        else:
            p = Poly(expr)

    conds = []
    if rv is None:
        e = p.as_expr()  # this is in expanded form
        # Do a safe inversion of e, moving non-s terms
        # to the rhs and dividing by a nonzero factor if
        # the relational is Eq/Ne; for other relationals
        # the sign must also be positive or negative
        rhs = 0
        b, ax = e.as_independent(s, as_Add=True)
        e -= b
        rhs -= b
        ef = factor_terms(e)
        a, e = ef.as_independent(s, as_Add=False)
        if (a.is_zero != False or  # don't divide by potential 0
                a.is_negative ==
                a.is_positive is None and  # if sign is not known then
                ie.rel_op not in ('!=', '==')): # reject if not Eq/Ne
            e = ef
            a = S.One
        rhs /= a
        if a.is_positive:
            rv = ie.func(e, rhs)
        else:
            rv = ie.reversed.func(e, rhs)

        # return conditions under which the value is
        # valid, too.
        beginning_denoms = denoms(ie.lhs) | denoms(ie.rhs)
        current_denoms = denoms(rv)
        for d in beginning_denoms - current_denoms:
            c = _solve_inequality(Eq(d, 0), s, linear=linear)
            if isinstance(c, Eq) and c.lhs == s:
                if classify(rv, s, c.rhs) is S.true:
                    # rv is permitting this value but it shouldn't
                    conds.append(~c)
        for i in (-oo, oo):
            if (classify(rv, s, i) is S.true and
                    classify(ie, s, i) is not S.true):
                conds.append(s < i if i is oo else i < s)

    conds.append(rv)
    return And(*conds)


def _reduce_inequalities(inequalities, symbols):
    # helper for reduce_inequalities

    poly_part, abs_part = {}, {}
    other = []

    for inequality in inequalities:

        expr, rel = inequality.lhs, inequality.rel_op  # rhs is 0

        # check for gens using atoms which is more strict than free_symbols to
        # guard against EX domain which won't be handled by
        # reduce_rational_inequalities
        gens = expr.atoms(Symbol)

        if len(gens) == 1:
            gen = gens.pop()
        else:
            common = expr.free_symbols & symbols
            if len(common) == 1:
                gen = common.pop()
                other.append(_solve_inequality(Relational(expr, 0, rel), gen))
                continue
            else:
                raise NotImplementedError(filldedent('''
                    inequality has more than one symbol of interest.
                    '''))

        if expr.is_polynomial(gen):
            poly_part.setdefault(gen, []).append((expr, rel))
        else:
            components = expr.find(lambda u:
                u.has(gen) and (
                u.is_Function or u.is_Pow and not u.exp.is_Integer))
            if components and all(isinstance(i, Abs) for i in components):
                abs_part.setdefault(gen, []).append((expr, rel))
            else:
                other.append(_solve_inequality(Relational(expr, 0, rel), gen))

    poly_reduced = [reduce_rational_inequalities([exprs], gen) for gen, exprs in poly_part.items()]
    abs_reduced = [reduce_abs_inequalities(exprs, gen) for gen, exprs in abs_part.items()]

    return And(*(poly_reduced + abs_reduced + other))


def reduce_inequalities(inequalities, symbols=[]):
    """Reduce a system of inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy import reduce_inequalities

    >>> reduce_inequalities(0 <= x + 3, [])
    (-3 <= x) & (x < oo)

    >>> reduce_inequalities(0 <= x + y*2 - 1, [x])
    (x < oo) & (x >= 1 - 2*y)
    """
    if not iterable(inequalities):
        inequalities = [inequalities]
    inequalities = [sympify(i) for i in inequalities]

    gens = set().union(*[i.free_symbols for i in inequalities])

    if not iterable(symbols):
        symbols = [symbols]
    symbols = (set(symbols) or gens) & gens
    if any(i.is_extended_real is False for i in symbols):
        raise TypeError(filldedent('''
            inequalities cannot contain symbols that are not real.
            '''))

    # make vanilla symbol real
    recast = {i: Dummy(i.name, extended_real=True)
        for i in gens if i.is_extended_real is None}
    inequalities = [i.xreplace(recast) for i in inequalities]
    symbols = {i.xreplace(recast) for i in symbols}

    # prefilter
    keep = []
    for i in inequalities:
        if isinstance(i, Relational):
            i = i.func(i.lhs.as_expr() - i.rhs.as_expr(), 0)
        elif i not in (True, False):
            i = Eq(i, 0)
        if i == True:
            continue
        elif i == False:
            return S.false
        if i.lhs.is_number:
            raise NotImplementedError(
                "could not determine truth value of %s" % i)
        keep.append(i)
    inequalities = keep
    del keep

    # solve system
    rv = _reduce_inequalities(inequalities, symbols)

    # restore original symbols and return
    return rv.xreplace({v: k for k, v in recast.items()})


class UnboundedLinearProgrammingError(Exception):
    """
    A linear programing problem is said to be unbounded if its objective
    function can assume arbitrarily large values.

    Example
    =======

    Suppose you want to maximize
        2x
    subject to
        x >= 0

    There's no upper limit that 2x can take.
    """
    pass


class InfeasibleLinearProgrammingError(Exception):
    """
    A linear programing problem is considered infeasible if its constraint set
    is empty. That is, if the set of all vectors satisfying the contraints is
    empty, then the problem is infeasible.

    Example
    =======

    Suppose you want to maximize
        x
    subject to
        x >= 10
        x <= 9

    It's not possible for x to satisfy the given contraints.
    """
    pass


def _pivot(M, i, j):
    """
    The pivot element `M[i, j]` is inverted and the rest of the matrix modified
    and returned as a new matrix; original is left unmodified.

    Example
    =======

    >>> from sympy.matrices.dense import Matrix
    >>> from sympy.solvers.inequalities import _pivot
    >>> from sympy import var
    >>> Matrix(3, 3, var('a:i'))
    Matrix([
    [a, b, c],
    [d, e, f],
    [g, h, i]])
    >>> _pivot(_, 1, 0)
    Matrix([
    [-a/d, -a*e/d + b, -a*f/d + c],
    [ 1/d,        e/d,        f/d],
    [-g/d,  h - e*g/d,  i - f*g/d]])
    """
    Mi, Mj, Mij = M[i,:], M[:,j], M[i,j]
    if Mij == 0:
        raise ZeroDivisionError("Tried to pivot about zero-valued entry.")
    A = M - Mj * (Mi / Mij)
    A[i, :] = Mi / Mij
    A[:, j] = -Mj / Mij
    A[i, j] = 1 / Mij
    return A


def _choose_pivot_row(A, B, candidate_rows, pivot_col, S):
    # Choose row with smallest ratio
    first_row = candidate_rows[0]
    min_ratio = B[first_row] / A[first_row, pivot_col]
    min_rows = [first_row]
    for i in candidate_rows[1:]:
        ratio = B[i] / A[i, pivot_col]
        if ratio < min_ratio:
            min_ratio = ratio
            min_rows = [i]
        elif ratio == min_ratio:
            min_rows.append(i)

    # If there are ties, pick using Bland's rule
    row = sorted(min_rows, key= lambda r: S[r])[0]
    return row


def _simplex(A, B, C, D=None, dual=False):
    """
    Return ``(o, x, y)`` obtained from the two-phase simplex method
    using Bland's rule where ``o`` is the minimum value of primal, ``Cx - D``,
    under constraints ``Ax <= B`` (with ``x >= 0``) and the maximum
    of the dual, ``y^{T}B - D``, under constraints
    ``A^{T}*y >= C^{T}`` (with ``y >= 0``). To compute the dual of
    the system, pass `dual=True` and ``(o, y, x)`` will be returned.

    Note: the nonnegative constraints for ``x`` and ``y`` supercede
    any values of ``A`` and ``B`` that are inconsistent with that
    assumption, so if a constraint of ``x >= -1`` is represented
    in ``A`` and ``B``, no value will be obtained that is negative; if
    a constraint of ``x <= -1`` is represented, an error will be
    raised since no solution is possible.

    Examples
    ========

    >>> from sympy.solvers.inequalities import _simplex
    >>> from sympy import Matrix

    Consider the simple minimization of ``f = x + y + 1`` under the constraint that
    ``y + 2*x >= 4``. In the nonnegative quadrant, this inequality describes a
    area above a triangle with vertices at (0, 4), (0, 0) and (2, 0). The minimum of
    f occurs at (2, 0). Defining A, B, C, D for the standard minimization gives

    >>> A, B, C, D = [Matrix(i) for i in [[[2, 1]], [4], [[1, 1]], [-1]]]

    Since `_simplex` will do a minimization for constraints given as ``A*x <= B``,
    the signs of each are negated (``-Ax <= -B``):

    >>> _simplex(-A, -B, C, D)
    (3, [2, 0], [1/2])

    The dual of minimizing ``f`` is maximizing ``F = c*y - d`` for ``a*y <= b`` where
    ``a``, ``b``, ``c``, ``d`` are derived from the transpose of the matrix
    representation of the standard minimization:

    >>> tr = lambda a, b, c, d: [i.T for i in (a, c, b, d)]
    >>> a, b, c, d = tr(A, B, C, D)

    This time ``a*x <= b`` is the expected inequality for the `_simplex` method, but
    to maximize ``F``, the sign of ``c`` and ``d`` must be inverted (so that minimizing
    the negative will give the negative of the maximum of ``F``):

    >>> _simplex(a, b, -c, -d)
    (-3, [1/2], [2, 0])

    The negated max shows that the max of ``F`` and the min of ``f`` are the same. The dual
    point `[1/2]` is the value of ``y`` that minimized ``F = c*y - d`` under
    constraints a*x <= b``:

    >>> y = Matrix(['y'])
    >>> (c*y - d)[0]
    4*y + 1
    >>> [i <= j for i, j in zip(a*y,b)]
    [2*y <= 1, y <= 1]

    In this 1-dimensional dual system, the more restrictive contraint is the first
    which limits ``y`` between 0 and 1/2 and the maximum of ``F`` is attained at the
    nonzero value, hence is ``4*(1/2) + 1 = 3``.

    In this case
    the values for ``x`` and ``y`` were the same when the dual representation was solved. This
    is not always the case (though the value of the function will be the same).

    >>> l = [[1, 1], [-1, 1], [0, 1], [-1, 0]], [5, 1, 2, -1], [[1, 1]], [-1]
    >>> A, B, C, D = [Matrix(i) for i in l]
    >>> _simplex(A, B, -C, -D)
    (-6, [3, 2], [1, 0, 0, 0])
    >>> _simplex(A, B, -C, -D, dual=True)  # dual point of [5, 0] != [3, 2]
    (-6, [1, 0, 0, 0], [5, 0])

    In both cases the function has the same value:

    >>> Matrix(C)*Matrix([3, 2]) == Matrix(C)*Matrix([5, 0])
    True

    See Also
    ========
    lp - poses min/max problem in form compatible with _simplex


    References
    ==========

    .. [1] Thomas S. Ferguson, LINEAR PROGRAMMING: A Concise Introduction
           web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
    """
    from sympy.matrices.dense import Matrix

    A, B, C, D = [Matrix(i) for i in (A, B, C, D or [0])]
    if dual:
        _o, d, p = _simplex(-A.T, C.T, B.T, -D)
        return -_o, d, p

    if A and B:
        M = Matrix([[A, B], [C, D]])
    else:
        assert not A and not B  # no constraints
        M = Matrix([[C, D]])
    n = M.cols - 1
    m = M.rows - 1

    if not all(i.is_Float or i.is_Rational for i in M):
        # with literal Float and Rational we are guaranteed the
        # ability of determining whether an expression is 0 or not
        raise TypeError("Only rationals and floats are allowed in the Simplex method.")

    # x variables are given priority over the y variables during Bland's rule
    # since False < True
    X = [(False, j) for j in range(n)]
    Y = [(True, i)  for i in range(m)]

    # Phase 1: find a feasible solution or determine none exist
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        if all(B[i] >= 0 for i in range(B.rows)):
            # We have found a feasible solution
            break

        # Find k: first row with a negative rightmost entry
        for k in range(B.rows):
            if B[k] < 0:
                break  # use current value of k below
        else:
            pass  # XXX is it an error if none was found?

        # Choose pivot column, c
        piv_cols = [_ for _ in range(A.cols) if A[k, _] < 0]
        if not piv_cols:
            raise InfeasibleLinearProgrammingError('The constraint set is empty!')
        c = sorted(piv_cols, key=lambda _: X[_])[0] # Bland's rule

        # Choose pivot row, r
        piv_rows = [_ for _ in range(A.rows) if A[_, c] > 0 and B[_] > 0]
        piv_rows.append(k)
        r = _choose_pivot_row(A, B, piv_rows, c, Y)

        M = _pivot(M, r, c)
        X[c], Y[r] = Y[r], X[c]

    # Phase 2: starting at a feasible solution, pivot until we reach optimal solution
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        C = M[-1, :-1]

        # Choose a pivot column, c
        piv_cols = []
        piv_cols = [_ for _ in range(n) if C[_] < 0]
        if not piv_cols:
            break
        c = sorted(piv_cols, key=lambda _: X[_])[0] # Bland's rule

        # Choose a pivot row, r
        piv_rows = [_ for _ in range(m) if A[_, c] > 0]
        if not piv_rows:
            raise UnboundedLinearProgrammingError('Objective function can assume arbitrarily large values!')
        r = _choose_pivot_row(A, B, piv_rows, c, Y)

        M = _pivot(M, r, c)
        X[c], Y[r] = Y[r], X[c]

    argmax = [None]*n
    argmin_dual = [None]*m

    for i, (v, n) in enumerate(X):
        if v == False:
            argmax[n] = 0
        else:
            argmin_dual[n] = M[-1, i]

    for i, (v, n) in enumerate(Y):
        if v == True:
            argmin_dual[n] = 0
        else:
            argmax[n] = M[i, -1]

    return -M[-1, -1], argmax, argmin_dual


def _np(constr, unbound=None):
    """return a (np, d, aux) where np is a list of nonpositive expressions
    that represent the given constraints (possibly rewritten in terms of
    auxilliary variables) expressible with nonnegative symbols, and d is
    a dictionary mapping a given symbols to an expression involving one
    or two auxilliary variables.

    If any constraint is False/empty, return None.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.solvers.inequalities import _np
    >>> from sympy.abc import x, y
    >>> _np([x >= y])
    ([-x + y], {}, [])
    >>> _np([x >= -oo])
    ([], {x: u1 - u2}, [u1, u2])
    >>> _np([x >= 3, x <= 5])
    ([u1 - 2], {x: u1 + 3}, [u1])
    >>> _np([x <= 5])
    ([], {x: 5 - u1}, [u1])
    >>> _np([x >= 1])
    ([], {x: u1 + 1}, [u1])
    """
    r = {}  # replacements to handle change of variables
    np = []  # nonpositive expressions
    aux = []  # will contain all symbols when done
    ui = numbered_symbols('u', start=1) # symbols to be introduced
    univariate = {}  # bound to be inferred from univariate constraints
    unbound = unbound or []  # symbols designated as unbound
    for i in constr:
        if isinstance(i, AppliedBinaryRelation):
            if i.function == Q.le:
                i = Le(i.lhs, i.rhs, evaluate=False)
            elif i.function == Q.ge:
                i = Le(i.rhs, i.lhs, evaluate=False)
            elif i.function == Q.eq:
                i = Eq(i.rhs, i.lhs, evaluate=False)
            else:
                assert None, i  # not allowed
        if i == True:
            continue  # ignore
        if i == False:
            return  # no solution
        if isinstance(i, Eq):  # could also collect k equalities and solve for k variables and eliminate them
            L = i.lhs - i.rhs
            np.extend([L, -L])
        elif isinstance(i, (Le, Ge)):
            i = i.canonical
            npi = i.lts - i.gts
            x = npi.free_symbols
            if len(x) == 1:
                x = x.pop()
                if x in unbound:
                    continue  # will handle later
                ivl = Le(npi, 0, evaluate=False).as_set()
                if x not in univariate:
                    univariate[x] = ivl
                else:
                    univariate[x] &= ivl
            else:
                np.append(npi)
        else:
            assert None, i  # not allowed

    # introduce auxilliary variables as needed for univariate
    # inequalities
    for x in univariate:
        i = univariate[x]
        if not i:
            return None  # no solution possible
        a, b = i.inf, i.sup
        if a.is_infinite and b.is_infinite:
            # this is unbound
            u = next(ui)
            v = next(ui)
            r[x] = u - v
            aux.extend([u, v])
        elif a.is_infinite:
            u = next(ui)
            r[x] = b - u
            aux.append(u)
        elif b.is_infinite:
            if a:
                u = next(ui)
                r[x] = a + u
                aux.append(u)
            else:
                # standard nonnegative relationship
                pass
        else:
            u = next(ui)
            aux.append(u)
            # shift so u = x - a => x = u + a
            r[x] = u + a
            # add constraint for u <= b - a
            # since when u = b-a then x = u + a = b - a + a = b:
            # the upper limit for x
            np.append(u - (b - a))

    # make change of variables for unbound variables
    for i in unbound:
        assert i not in univariate
        u, v = next(ui), next(ui)
        r[i] = u - v
        aux.extend([u, v])

    return np, r, aux


def lp(min_max, f, constr, unbound=None):
    """Return the optimization (min or max) of ``f`` with the given
    constraints; if any variables are unbound, pass them as a list for
    `unbound`.

    If `how=max` then the results corresponding to the maximization
    of ``f`` will be returned.

    Examples
    ========

    >>> from sympy.solvers.inequalities import lp
    >>> from sympy import symbols, Eq
    >>> from sympy.abc import x, y
    >>> x1, x2, x3, x4 = symbols('x1:5')
    >>> f = 5*x2 + x3 + 4*x4
    >>> c = [x1 + 5 >= 5*x2 + 2*x3 + 5*x4, Eq(3*x2 + x4, 2), Eq(-x1 + x3 + 2*x4, 1)]
    >>> lp(min, f, c)
    (9/2, {x1: 0, x2: 1/2, x3: 0, x4: 1/2})
    >>> lp(max, f, c, unbound=[x3])
    (6, {x1: 1, x2: 0, x3: -2, x4: 2})

    >>> lp(min, x - y, [x >= 3, x + y <= 7])
    (-1, {x: 3, y: 4})

    >>> lp(max, x - y, [x >= 3, x + y <= 7])
    (7, {x: 7, y: 0})

    >>> lp(max, x - y, [x <= 3, x + y <= 7])
    (3, {x: 3, y: 0})

    XXX how to define error so it prints better (without the __main__)

    >>> lp(min, x - y, [x <= 3, x + y <= 7], [])
    Traceback (most recent call last):
    ...
    __main__.UnboundedLinearProgrammingError: Objective function can assume arbitrarily large values!
    """
    how = min_max
    F = sympify(f)
    syms = F.free_symbols
    for i, j in enumerate(constr):
        constr[i] = sympify(j)
        syms |= constr[i].free_symbols
    if unbound is True:
        unbound = syms
    elif not unbound:
        unbound = []
    elif not all(i in syms for i in unbound):
        raise ValueError(filldedent('''
            all symbols in unbound should appear
            in the objective or constraints'''))

    # convert constraints to nonpositive expressions
    _ = _np(constr, unbound)
    if _ is None:
        raise InfeasibleLinearProgrammingError('inconsistent/False constraint')
    np, r, aux = _

    # do change of variables
    f = f.xreplace(r)
    np = [i.xreplace(r) for i in np]

    # convert to matrices and solve
    xx = list(ordered(syms)) + aux
    A, B = linear_eq_to_matrix(np, xx)
    C, D = linear_eq_to_matrix([f], xx)
    if how == max:
        _o, p, d = _simplex(A, B, -C, -D)
        o = -_o
    else:
        o, p, d = _simplex(A, B, C, D)

    # restore original variables and remove aux from p
    p = dict(zip(xx, p))
    if r:
        p.update(Dict(r).xreplace({k: p.pop(k) for k in aux}))
        # the o returned needs to be translated back to the
        # given function, not the one in terms of changes of
        # variables
        o = F.xreplace(p)

    # make p canonical
    p = {i: p[i] for i in ordered(p) if i}

    # XXX not sure that dual should be returned since new variables may
    # have been introduced; if the user wants the dual the should get it
    # from the standard min/max problem representation -- this function
    # does not enforce standard form, it gets the input into standard
    # form for simplex
    return o, p
