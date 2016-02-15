"""Tools for solving inequalities and systems of inequalities. """

from __future__ import print_function, division

from sympy.core import Symbol, Dummy, sympify
from sympy.core.compatibility import iterable, reduce
from sympy.sets import Interval
from sympy.core.relational import Relational, Eq, Ge, Lt
from sympy.sets.sets import FiniteSet, Union
from sympy.core.singleton import S

from sympy.functions import Abs
from sympy.logic import And
from sympy.polys import Poly, PolynomialError, parallel_poly_from_expr
from sympy.polys.polyutils import _nsort
from sympy.utilities.misc import filldedent

def solve_poly_inequality(poly, rel):
    """Solve a polynomial inequality with rational coefficients.

    Examples
    ========

    >>> from sympy import Poly
    >>> from sympy.abc import x
    >>> from sympy.solvers.inequalities import solve_poly_inequality

    >>> solve_poly_inequality(Poly(x, x, domain='ZZ'), '==')
    [{0}]

    >>> solve_poly_inequality(Poly(x**2 - 1, x, domain='ZZ'), '!=')
    [(-oo, -1), (-1, 1), (1, oo)]

    >>> solve_poly_inequality(Poly(x**2 - 1, x, domain='ZZ'), '==')
    [{-1}, {1}]

    See Also
    ========
    solve_poly_inequalities
    """
    if not isinstance(poly, Poly):
        raise ValueError(
            'For efficiency reasons, `poly` should be a Poly instance')
    if poly.is_number:
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

    >>> from sympy.solvers.inequalities import solve_poly_inequalities
    >>> from sympy.polys import Poly
    >>> from sympy.abc import x
    >>> solve_poly_inequalities(((
    ... Poly(x**2 - 3), ">"), (
    ... Poly(-x**2 + 1), ">")))
    (-oo, -sqrt(3)) U (-1, 1) U (sqrt(3), oo)
    """
    from sympy import Union
    return Union(*[solve_poly_inequality(*p) for p in polys])


def solve_rational_inequalities(eqs):
    """Solve a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy import Poly
    >>> from sympy.solvers.inequalities import solve_rational_inequalities

    >>> solve_rational_inequalities([[
    ... ((Poly(-x + 1), Poly(1, x)), '>='),
    ... ((Poly(-x + 1), Poly(1, x)), '<=')]])
    {1}

    >>> solve_rational_inequalities([[
    ... ((Poly(x), Poly(1, x)), '!='),
    ... ((Poly(-x + 1), Poly(1, x)), '>=')]])
    (-oo, 0) U (0, 1]

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

            for numer_interval in numer_intervals:
                for global_interval in global_intervals:
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

    >>> from sympy import Poly, Symbol
    >>> from sympy.solvers.inequalities import reduce_rational_inequalities

    >>> x = Symbol('x', real=True)

    >>> reduce_rational_inequalities([[x**2 <= 0]], x)
    Eq(x, 0)

    >>> reduce_rational_inequalities([[x + 2 > 0]], x)
    And(-2 < x, x < oo)
    >>> reduce_rational_inequalities([[(x + 2, ">")]], x)
    And(-2 < x, x < oo)
    >>> reduce_rational_inequalities([[x + 2]], x)
    Eq(x, -2)
    """
    exact = True
    eqs = []
    solution = S.EmptySet
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
                    only polynomials and
                    rational functions are supported in this context'''))

            if not opt.domain.is_Exact:
                numer, denom, exact = numer.to_exact(), denom.to_exact(), False

            domain = opt.domain.get_exact()

            if not (domain.is_ZZ or domain.is_QQ):
                expr = numer/denom
                expr = Relational(expr, 0, rel)
                solution = Union(solution, solve_univariate_inequality(expr, gen, relational=False))
            else:
                _eqs.append(((numer, denom), rel))

        eqs.append(_eqs)

    solution = Union(solution, solve_rational_inequalities(eqs))

    if not exact:
        solution = solution.evalf()

    if relational:
        solution = solution.as_relational(gen)

    return solution


def reduce_abs_inequality(expr, rel, gen):
    """Reduce an inequality with nested absolute values.

    Examples
    ========

    >>> from sympy import Abs, Symbol
    >>> from sympy.solvers.inequalities import reduce_abs_inequality
    >>> x = Symbol('x', real=True)

    >>> reduce_abs_inequality(Abs(x - 5) - 3, '<', x)
    And(2 < x, x < 8)

    >>> reduce_abs_inequality(Abs(x + 2)*3 - 13, '<', x)
    And(-19/3 < x, x < 7/3)

    See Also
    ========

    reduce_abs_inequalities
    """
    """if gen.is_real is False:
         raise TypeError(filldedent('''
            can't solve inequalities with absolute
            values containing non-real variables'''))"""

    def _bottom_up_scan(expr):
        exprs = []

        if expr.is_Add or expr.is_Mul:
            op = expr.func

            for arg in expr.args:
                _exprs = _bottom_up_scan(arg)

                if not exprs:
                    exprs = _exprs
                else:
                    args = []

                    for expr, conds in exprs:
                        for _expr, _conds in _exprs:
                            args.append((op(expr, _expr), conds + _conds))

                    exprs = args
        elif expr.is_Pow:
            n = expr.exp
            
            if not n.is_Integer:
                raise ValueError("Only Integer Powers are allowed on Abs.")

            _exprs = _bottom_up_scan(expr.base)

            for expr, conds in _exprs:
                exprs.append((expr**n, conds))
        elif isinstance(expr, Abs):
            _exprs = _bottom_up_scan(expr.args[0])

            for expr, conds in _exprs:
                exprs.append(( expr, conds + [Ge(expr, 0)]))
                exprs.append((-expr, conds + [Lt(expr, 0)]))
        else:
            exprs = [(expr, [])]

        return exprs

    exprs = _bottom_up_scan(expr)

    mapping = {'<': '>', '<=': '>='}
    inequalities = []

    for expr, conds in exprs:
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

    >>> from sympy import Abs, Symbol
    >>> from sympy.abc import x
    >>> from sympy.solvers.inequalities import reduce_abs_inequalities
    >>> x = Symbol('x', real=True)

    >>> reduce_abs_inequalities([(Abs(3*x - 5) - 7, '<'),
    ... (Abs(x + 25) - 13, '>')], x)
    And(-2/3 < x, Or(And(-12 < x, x < oo), And(-oo < x, x < -38)), x < 4)

    >>> reduce_abs_inequalities([(Abs(x - 4) + Abs(3*x - 5) - 7, '<')], x)
    And(1/2 < x, x < 4)

    See Also
    ========

    reduce_abs_inequality
    """
    return And(*[ reduce_abs_inequality(expr, rel, gen)
        for expr, rel in exprs ])


def solve_univariate_inequality(expr, gen, relational=True):
    """Solves a real univariate inequality.

    Examples
    ========

    >>> from sympy.solvers.inequalities import solve_univariate_inequality
    >>> from sympy.core.symbol import Symbol
    >>> x = Symbol('x')

    >>> solve_univariate_inequality(x**2 >= 4, x)
    Or(And(-oo < x, x <= -2), And(2 <= x, x < oo))

    >>> solve_univariate_inequality(x**2 >= 4, x, relational=False)
    (-oo, -2] U [2, oo)

    """

    from sympy.solvers.solvers import solve, denoms

    # This keeps the function independent of the assumptions about `gen`.
    # `solveset` makes sure this function is called only when the domain is
    # real.
    d = Dummy(real=True)
    expr = expr.subs(gen, d)
    _gen = gen
    gen = d

    if expr is S.true:
        rv = S.Reals
    elif expr is S.false:
        rv = S.EmptySet
    else:
        e = expr.lhs - expr.rhs
        parts = n, d = e.as_numer_denom()
        if all(i.is_polynomial(gen) for i in parts):
            solns = solve(n, gen, check=False)
            singularities = solve(d, gen, check=False)
        else:
            solns = solve(e, gen, check=False)
            singularities = []
            for d in denoms(e):
                singularities.extend(solve(d, gen))

        include_x = expr.func(0, 0)

        def valid(x):
            v = e.subs(gen, x)
            try:
                r = expr.func(v, 0)
            except TypeError:
                r = S.false
            if r in (S.true, S.false):
                return r
            if v.is_real is False:
                return S.false
            else:
                v = v.n(2)
                if v.is_comparable:
                    return expr.func(v, 0)
                return S.false

        start = S.NegativeInfinity
        sol_sets = [S.EmptySet]
        try:
            reals = _nsort(set(solns + singularities), separated=True)[0]
        except NotImplementedError:
            raise NotImplementedError('sorting of these roots is not supported')
        for x in reals:
            end = x

            if end in [S.NegativeInfinity, S.Infinity]:
                if valid(S(0)):
                    sol_sets.append(Interval(start, S.Infinity, True, True))
                    break

            if valid((start + end)/2 if start != S.NegativeInfinity else end - 1):
                sol_sets.append(Interval(start, end, True, True))

            if x in singularities:
                singularities.remove(x)
            elif include_x:
                sol_sets.append(FiniteSet(x))

            start = end

        end = S.Infinity

        # in case start == -oo then there were no solutions so we just
        # check a point between -oo and oo (e.g. 0) else pick a point
        # past the last solution (which is start after the end of the
        # for-loop above
        if valid(start + 1 if start is not S.NegativeInfinity else 0):
            sol_sets.append(Interval(start, end, True, True))

        rv = Union(*sol_sets).subs(gen, _gen)

    return rv if not relational else rv.as_relational(_gen)


def _solve_inequality(ie, s):
    expr = ie.lhs - ie.rhs
    try:
        p = Poly(expr, s)
        if p.degree() != 1:
            raise NotImplementedError
    except (PolynomialError, NotImplementedError):
        try:
            return reduce_rational_inequalities([[ie]], s)
        except PolynomialError:
            return solve_univariate_inequality(ie, s)
    a, b = p.all_coeffs()
    if a.is_positive or ie.rel_op in ('!=', '=='):
        return ie.func(s, -b/a)
    elif a.is_negative:
        return ie.reversed.func(s, -b/a)
    else:
        raise NotImplementedError


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
                    inequality has more than one
                    symbol of interest'''))

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

    poly_reduced = []
    abs_reduced = []

    for gen, exprs in poly_part.items():
        poly_reduced.append(reduce_rational_inequalities([exprs], gen))

    for gen, exprs in abs_part.items():
        abs_reduced.append(reduce_abs_inequalities(exprs, gen))

    return And(*(poly_reduced + abs_reduced + other))


def reduce_inequalities(inequalities, symbols=[]):
    """Reduce a system of inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy import sympify as S, Symbol
    >>> from sympy.abc import x, y
    >>> from sympy.solvers.inequalities import reduce_inequalities

    >>> reduce_inequalities(0 <= x + 3, [])
    And(-3 <= x, x < oo)

    >>> reduce_inequalities(0 <= x + y*2 - 1, [x])
    x >= -2*y + 1
    """
    if not iterable(inequalities):
        inequalities = [inequalities]
    inequalities = [sympify(i) for i in inequalities]

    gens = set().union(*[i.free_symbols for i in inequalities])

    if not iterable(symbols):
        symbols = [symbols]
    symbols = (set(symbols) or gens) & gens
    if any(i.is_real is False for i in symbols):
        raise TypeError(filldedent('''
            inequalities cannot contain symbols that are not real.'''))

    # make vanilla symbol real
    recast = dict([(i, Dummy(i.name, real=True))
        for i in gens if i.is_real is None])
    inequalities = [i.xreplace(recast) for i in inequalities]
    symbols = set([i.xreplace(recast) for i in symbols])

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
    return rv.xreplace(dict([(v, k) for k, v in recast.items()]))
