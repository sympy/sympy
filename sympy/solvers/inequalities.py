"""Tools for solving inequalities and systems of inequalities. """

from __future__ import print_function, division

from sympy.core import Symbol
from sympy.sets import Interval
from sympy.core.relational import Relational, Eq, Ge, Lt
from sympy.sets.sets import FiniteSet, Union
from sympy.core.singleton import S

from sympy.functions import re, im, Abs
from sympy.logic import And
from sympy.polys import Poly, PolynomialError, parallel_poly_from_expr
from sympy.simplify import simplify


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
        global_intervals = [Interval(S.NegativeInfinity, S.Infinity)]

        for (numer, denom), rel in _eqs:
            numer_intervals = solve_poly_inequality(numer*denom, rel)
            denom_intervals = solve_poly_inequality(denom, '==')

            if global_intervals is None:
                global_intervals = numer_intervals
            else:
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
    x == 0

    >>> reduce_rational_inequalities([[x + 2 > 0]], x)
    And(-2 < x, x < oo)
    >>> reduce_rational_inequalities([[(x + 2, ">")]], x)
    And(-2 < x, x < oo)
    >>> reduce_rational_inequalities([[x + 2]], x)
    x == -2
    """
    exact = True
    eqs = []

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
                continue
            elif expr is S.false:
                return S.EmptySet if not relational else expr

            try:
                (numer, denom), opt = parallel_poly_from_expr(
                    expr.together().as_numer_denom(), gen)
            except PolynomialError:
                raise PolynomialError("only polynomials and "
                    "rational functions are supported in this context")

            if not opt.domain.is_Exact:
                numer, denom, exact = numer.to_exact(), denom.to_exact(), False

            domain = opt.domain.get_exact()

            if not (domain.is_ZZ or domain.is_QQ):
                raise NotImplementedError(
                    "inequality solving is not supported over %s" % opt.domain)

            _eqs.append(((numer, denom), rel))

        eqs.append(_eqs)

    solution = solve_rational_inequalities(eqs)

    if not exact:
        solution = solution.evalf()

    if not relational:
        return solution

    real = gen.is_real

    if not real:
        result = And(solution.as_relational(re(gen)), Eq(im(gen), 0))
    else:
        result = solution.as_relational(gen)

    return result


def reduce_abs_inequality(expr, rel, gen):
    """Reduce an inequality with nested absolute values.

    Examples
    ========

    >>> from sympy import Q, Abs, Symbol
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
    if not gen.is_real:
        raise NotImplementedError("can't solve inequalities with absolute "
                                  "values of a complex variable")

    def _bottom_up_scan(expr):
        exprs = []

        if expr.is_Add or expr.is_Mul:
            op = expr.__class__

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

            if not n.is_Integer or n < 0:
                raise ValueError(
                    "only non-negative integer powers are allowed")

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

    >>> from sympy import Q, Abs, Symbol
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
    >>> x = Symbol('x', real=True)

    >>> solve_univariate_inequality(x**2 >= 4, x)
    Or(And(-oo < x, x <= -2), And(2 <= x, x < oo))

    >>> solve_univariate_inequality(x**2 >= 4, x, relational=False)
    (-oo, -2] U [2, oo)

    """

    # Implementation for continous functions

    from sympy.solvers.solvers import solve

    solns = solve(expr.lhs - expr.rhs, gen)
    oo = S.Infinity

    start = -oo

    sol_sets = [S.EmptySet]

    for x in sorted(s for s in solns if s.is_real):
        end = x
        if simplify(expr.subs(gen, (start + end)/2 if start != -oo else end - 1)):
            sol_sets.append(Interval(start, end, True, True))

        if simplify(expr.subs(gen, x)):
            sol_sets.append(FiniteSet(x))

        start = end

    end = oo

    if simplify(expr.subs(gen, start + 1)):
        sol_sets.append(Interval(start, end, True, True))

    rv = Union(*sol_sets)
    return rv if not relational else rv.as_relational(gen)


def _solve_inequality(ie, s):
    """ A hacky replacement for solve, since the latter only works for
        univariate inequalities. """
    if not ie.rel_op in ('>', '>=', '<', '<='):
        raise NotImplementedError
    expr = ie.lhs - ie.rhs
    try:
        p = Poly(expr, s)
        if p.degree() != 1:
            raise NotImplementedError
    except (PolynomialError, NotImplementedError):
        try:
            n, d = expr.as_numer_denom()
            return reduce_rational_inequalities([[ie]], s)
        except PolynomialError:
            return solve_univariate_inequality(ie, s)
    a, b = p.all_coeffs()
    if a.is_positive:
        return ie.func(s, -b/a)
    elif a.is_negative:
        return ie.func(-b/a, s)
    else:
        raise NotImplementedError


def reduce_inequalities(inequalities, symbols=[]):
    """Reduce a system of inequalities with rational coefficients.

    Examples
    ========

    >>> from sympy import Q, sympify as S, Symbol
    >>> from sympy.abc import x, y
    >>> from sympy.solvers.inequalities import reduce_inequalities

    >>> x = Symbol('x', real=True)
    >>> reduce_inequalities(S(0) <= x + 3, [])
    And(-3 <= x, x < oo)

    >>> x = Symbol('x')
    >>> reduce_inequalities(S(0) <= x + y*2 - 1, [x])
    -2*y + 1 <= x
    """
    if not hasattr(inequalities, '__iter__'):
        inequalities = [inequalities]

    if len(inequalities) == 1 and len(symbols) == 1 \
            and inequalities[0].is_Relational:
        try:
            return _solve_inequality(inequalities[0], symbols[0])
        except NotImplementedError:
            pass

    poly_part, abs_part = {}, {}

    for inequality in inequalities:
        if inequality == True:
            continue
        elif inequality == False:
            return False

        if inequality.is_Relational:
            expr, rel = inequality.lhs - inequality.rhs, inequality.rel_op
        else:
            expr, rel = inequality, '=='

        gens = expr.free_symbols

        if not gens:
            return False
        elif len(gens) == 1:
            gen = gens.pop()
        else:
            raise NotImplementedError(
                "only univariate inequalities are supported")

        components = expr.find(lambda u: u.is_Function)

        if not components:
            if gen in poly_part:
                poly_part[gen].append((expr, rel))
            else:
                poly_part[gen] = [(expr, rel)]
        else:
            if all(isinstance(comp, Abs) for comp in components):
                if gen in abs_part:
                    abs_part[gen].append((expr, rel))
                else:
                    abs_part[gen] = [(expr, rel)]
            else:
                raise NotImplementedError("can't reduce %s" % inequalities)

    poly_reduced = []
    abs_reduced = []

    for gen, exprs in poly_part.items():
        poly_reduced.append(reduce_rational_inequalities([exprs], gen))

    for gen, exprs in abs_part.items():
        abs_reduced.append(reduce_abs_inequalities(exprs, gen))

    return And(*(poly_reduced + abs_reduced))
