"""Tools for solving inequalities and systems of inequalities. """
import itertools

from sympy.calculus.util import (continuous_domain, periodicity,
    function_range)
from sympy.core import Symbol, Dummy, sympify
from sympy.core.exprtools import factor_terms
from sympy.core.relational import Relational, Ge, Lt, Gt, Eq, Ne, Le
from sympy.assumptions.ask import Q
from sympy.assumptions.relation.binrel import AppliedBinaryRelation
from sympy.matrices.immutable import ImmutableMatrix
from sympy.sets.sets import Interval, FiniteSet, Union, Intersection
from sympy.core.singleton import S
from sympy.core.function import expand_mul
from sympy.functions.elementary.complexes import im, Abs
from sympy.logic import And
from sympy.polys import Poly, PolynomialError, parallel_poly_from_expr
from sympy.polys.polyutils import _nsort
from sympy.solvers.solveset import solvify, solveset, linear_eq_to_matrix
from sympy.utilities.iterables import sift, iterable
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


def _simplex(A, B, C):
    """
    Two phase simplex method with Bland's rule

    References
    ==========

    .. [1] Thomas S. Ferguson, LINEAR PROGRAMMING: A Concise Introduction
           web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
    """
    from sympy.matrices.dense import Matrix

    D = ImmutableMatrix([0])
    M = Matrix([[A, B], [-C, D]])

    if not all(i.is_Float or i.is_Rational for i in M):
            raise TypeError("Only rationals and floats are allowed in the Simplex method.")

    # It's important that False < True so that x variables are given priority
    # over the y variables during Bland's rule
    r_orig = [(False, j) for j in range(M.cols - 1)] # what Ferguson's introduction calls 'x' variables
    s_orig = [(True, i)  for i in range(M.rows - 1)] # what Ferguson's introduction calls 'y' variables

    R = r_orig.copy()
    S = s_orig.copy()

    # Phase 1: find a feasible solution or determine none exist
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        if all(B[i] >= 0 for i in range(B.rows)):
            # We have found a feasible solution
            break

        # Find k, first row with a negative rightmost entry
        for k in range(B.rows):
            if B[k] < 0:
                break

        # Choose pivot column, j0
        piv_cols = [j for j in range(A.cols) if A[k, j] < 0]
        if not piv_cols:
            raise InfeasibleLinearProgrammingError('The constraint set is empty!')
        j0 = sorted(piv_cols, key=lambda c: R[c])[0] # Bland's rule

        # Choose pivot row, i0
        piv_rows = [i for i in range(A.rows) if A[i, j0] > 0 and B[i] > 0]
        piv_rows.append(k)
        i0 = _choose_pivot_row(A, B, piv_rows, j0, S)

        M = _pivot(M, i0, j0)
        R[j0], S[i0] = S[i0], R[j0]

    # Phase 2: starting at a feasible solution, pivot until we reach optimal solution
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        C = M[-1, :-1]

        # Choose a pivot column, j0
        piv_cols = []
        piv_cols = [j for j in range(C.cols) if C[j] < 0]
        if not piv_cols:
            break
        j0 = sorted(piv_cols, key=lambda c: R[c])[0] # Bland's rule

        # Choose a pivot row, i0
        piv_rows = [i for i in range(A.rows) if A[i, j0] > 0]
        if not piv_rows:
            raise UnboundedLinearProgrammingError('Objective function can assume arbitrarily large values!')
        i0 = _choose_pivot_row(A, B, piv_rows, j0, S)

        M = _pivot(M, i0, j0)
        R[j0], S[i0] = S[i0], R[j0]

    argmax = [None]*(M.cols-1)
    argmin_dual = [None]*(M.rows-1)

    for i, var in enumerate(R):
        v = var[0]
        n = var[1]
        if v == False:
            argmax[n] = 0
        else:
            argmin_dual[n] = M[-1, i]

    for i, var in enumerate(S):
        v = var[0]
        n = var[1]
        if v == True:
            argmin_dual[n] = 0
        else:
            argmax[n] = M[i, -1]

    return M[-1, -1], argmax, argmin_dual


def linprog_from_matrices(A, B, C):
    """
    Because of the duality of linear programming, there are two valid ways to
    interpret this function. The first is as a solver for the standard
    maximization problem:
        Maximizing Cx constrained to Ax <= B and x >= 0.
    The second as is a solver for the standard minimization problem:
        Minimizing y^{T}B constrained to y^{T}A >= C^{T} and y >= 0.

    x is a column vector of variables, and y is a column vector of dual
    variables.

    If we keep the maximization interpretation in mind, then matrix A and
    column vector B represent a set of constraints. For example,
        1*x + 0*y <= 2
        0*x + 0*y <= 4
    can be represented as:

    >>> from sympy.matrices.dense import Matrix
    >>> A = Matrix([[1, 0],
    ...             [0, 1]])
    >>> B = Matrix([[2],
    ...             [4]])

    In this interpretation, row vector C is an objective function to be
    maximized. For example,
        3x + y
    can be represented as:

    >>> C = Matrix([[3, 1]])

    Parameters
    ==========

    Under the maximization interpretation, m is the number of constraints and n
    is the number of variables.

    A : Matrix with shape (m, n)
        Must contain only Rational and float coefficients.

    B : Matrix with shape (m, 1)
        Must contain only Rational and float coefficients.

    C: Matrix with shape (1, n)
        Must contain only Rational and float coefficients.

    Returns
    =======

    Returns `(optimum, argmax, and argmax_dual)`:

    optimum : float or Rational
        maximum value of objective function under the constraints

    argmax : list of floats and/or Rationals
        x

    argmax_dual : list of floats and/or Rationals
        y

    Examples
    ========

    Suppose we want to maximize
        3x + y
    subject to
        x + y <= 2

    >>> from sympy.matrices.dense import Matrix
    >>> from sympy.solvers.inequalities import linprog_from_matrices
    >>> A = Matrix([[1, 1]])
    >>> B = Matrix([[2]])
    >>> C = Matrix([[3, 1]])
    >>> optimum, argmax, argmax_dual = linprog_from_matrices(A, B, C)
    >>> optimum
    6
    >>> argmax
    [2, 0]
    >>> argmax_dual
    [3]

    Suppose we want to maximize
        3x + y
    subject to
        x <= 2
        y <= 4

    >>> A = Matrix([[1 ,0], [0, 1]])
    >>> B = Matrix([[2], [4]])
    >>> C = Matrix([[3, 1]])
    >>> optimum, argmax, argmax_dual = linprog_from_matrices(A, B, C)
    >>> optimum
    10
    >>> argmax
    [2, 4]
    >>> argmax_dual
    [3, 1]

    See Also
    ========

    linprog_maximize_from_equations
    find_feasible
    """
    from sympy.matrices.dense import Matrix

    A, B, C = [Matrix(i) for i in (A, B, C)]

    m, n = A.shape
    if B.shape != (m, 1) or C.shape != (1, n):
        raise ValueError(f"The shape of matrix B ({B.shape}) or C ({C.shape})" \
                         f"does not match the shape of matrix A ({A.shape}).")
    return _simplex(A, B, C)


def _linear_programming_to_matrix(constraints, objective, variables):
    """
    Converts a list of constraints and an objective function into the standard
    form for linear programming.

    Examples
    ========

    >>> from sympy.solvers.inequalities import _linear_programming_to_matrix
    >>> from sympy.abc import x, y, z
    >>> from sympy import Eq
    >>> r1 = y + 2*z <= 3
    >>> r2 = -x - 3*z <= -2
    >>> r3 = 5 - y >= 2*x + 7*z
    >>> A, B, C, constraints = _linear_programming_to_matrix([r1, r2, r3], x + y + 5*z, [x, y, z])
    >>> A
    Matrix([
    [ 0, 1,  2],
    [-1, 0, -3],
    [ 2, 1,  7]])
    >>> B
    Matrix([
    [ 3],
    [-2],
    [ 5]])
    >>> C
    Matrix([[1, 1, 5]])
    >>> constraints
    [y + 2*z <= 3, -x - 3*z <= -2, 2*x + y + 7*z <= 5]
    >>> A, B, C, constraints = _linear_programming_to_matrix([Eq(x, 3)], x*10, [x])  # x = 3 become x >= 3 and x <= 3
    >>> A
    Matrix([
    [ 1],
    [-1]])
    >>> B
    Matrix([
    [ 3],
    [-3]])
    >>> C
    Matrix([[10]])
    >>> constraints
    [x <= 3, -x <= -3]
    """
    standard_constraints = []
    eqns = []

    for rel in constraints:
        if not isinstance(rel, (Relational, AppliedBinaryRelation)):
            raise TypeError(f"{rel} is not relational.")
        if type(rel) in [Lt, Gt] or isinstance(rel, AppliedBinaryRelation) and rel.function in [Q.lt, Q.gt]:
            raise TypeError("Strict inequalities are not allowed in linear programming.")
        if type(rel) == Ne or isinstance(rel, AppliedBinaryRelation) and rel.function == Q.ne:
            raise TypeError("'not equal to' is not allowed in linear programming.")

        if type(rel) == Le or isinstance(rel, AppliedBinaryRelation) and rel.function == Q.le:
            eqns.append(rel.lhs - rel.rhs)
            standard_constraints.append(rel.lhs <= rel.rhs)
        elif type(rel) == Ge or isinstance(rel, AppliedBinaryRelation) and rel.function == Q.ge:
            eqns.append(rel.rhs - rel.lhs)
            standard_constraints.append(-rel.lhs <= -rel.rhs)
        elif type(rel) == Eq or isinstance(rel, AppliedBinaryRelation) and rel.function == Q.eq:
            # x = 3 can be represented as x <= 3 and x >= 3
            eqns.append(rel.lhs - rel.rhs)
            standard_constraints.append(rel.lhs <= rel.rhs)
            eqns.append(rel.rhs - rel.lhs)
            standard_constraints.append(-rel.lhs <= -rel.rhs)
        else:
            raise TypeError(f"Unrecognized relation: {rel}")

    A, B = linear_eq_to_matrix(eqns, *variables)
    C = linear_eq_to_matrix(objective, *variables)[0] # constant terms can be safely ignored here

    return A, B, C, standard_constraints


def linprog_maximize_from_equations(constraints, objective, variables):
    """
    A function to maximize a linear objective function subject to linear
    constraints with the simplex method. All coefficients must be Rational or
    floats.

    Parameters
    ==========

    constraints : list of linear inequalities
        Stict inequalities (>, <) and not equals are not allowed.

    objective : a linear function to maximize
        If you want to minimize a function, f, simply make the objective
        function -f.

    variables : a list of variables included in constraints and objective

    Returns
    =======

    Returns `(optimum, argmax, argmax_dual)`:

    optimum : Rational or Float

    argmax : dictionary of variables to Rationals or Floats

    argmax_dual : dictionary of constraints to Rationals or Floats

    Examples
    ========

    >>> from sympy.abc import x, y, z
    >>> from sympy.solvers.inequalities import linprog_maximize_from_equations
    >>> r1 = y+2*z <= 3
    >>> r2 = -x-3*z <= -2
    >>> r3 = 2*x+y+7*z <= 5
    >>> optimum, argmax, argmax_dual  = linprog_maximize_from_equations([r1,r2,r3], x+y+5*z, [x, y, z])
    >>> optimum
    11/3
    >>> argmax
    {x: 0, y: 1/3, z: 2/3}
    >>> argmax_dual
    {-x - 3*z <= -2: 2/3, y + 2*z <= 3: 0, 2*x + y + 7*z <= 5: 1}

    See Also
    ========

    solve_univariate_inequality
    linprog_from_matrices
    find_feasible
    """
    A, B, C, standard_constraints = _linear_programming_to_matrix(constraints, objective, variables)
    optimum, argmax, argmax_dual = _simplex(A, B, C)
    argmax = {variables[i] : argmax[i] for i in range(len(variables))}
    argmax_dual = {standard_constraints[i] : argmax_dual[i] for i in range(len(standard_constraints))}
    return optimum, argmax, argmax_dual
