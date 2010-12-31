"""Tools for solving inequalities and systems of inequalities. """

from sympy.core import Symbol, Interval, Union
from sympy.core.relational import Relational, Eq, Ge, Lt
from sympy.core.singleton import S
from sympy.assumptions import ask, Assume
from sympy.functions import re, im, Abs
from sympy.logic import And, Or
from sympy.polys import Poly
from sympy.utilities import all

def interval_evalf(interval):
    """Proper implementation of evalf() on Interval. """
    return Interval(interval.left.evalf(), interval.right.evalf(),
        left_open=interval.left_open, right_open=interval.right_open)

def solve_poly_inequality(poly, rel):
    """Solve a polynomial inequality with rational coefficients.  """
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
                    intervals.insert(0, Interval(left, right, not equal, right_open))

                sign, right, right_open = -sign, left, not equal
            else:
                if sign == eq_sign and not equal:
                    intervals.insert(0, Interval(left, right, True, right_open))
                    right, right_open = left, True
                elif sign != eq_sign and equal:
                    intervals.insert(0, Interval(left, left))

        if sign == eq_sign:
            intervals.insert(0, Interval(S.NegativeInfinity, right, True, right_open))

    return intervals

def solve_poly_inequalities(polys):
    """Solve a system of polynomial inequalities with rational coefficients. """
    result = S.EmptySet

    for _polys in polys:
        global_intervals = None

        for poly, rel in _polys:
            local_intervals = solve_poly_inequality(poly, rel)

            if global_intervals is None:
                global_intervals = local_intervals
            else:
                intervals = []

                for local_interval in local_intervals:
                    for global_interval in global_intervals:
                        interval = local_interval.intersect(global_interval)

                        if interval is not S.EmptySet:
                            intervals.append(interval)

                global_intervals = intervals

            if not global_intervals:
                break

        for interval in global_intervals:
            result = result.union(interval)

    return result

def reduce_poly_inequalities(exprs, gen, assume=True, relational=True):
    """Reduce a system of polynomial inequalities with rational coefficients. """
    exact = True
    polys = []

    for _exprs in exprs:
        _polys = []

        for expr in _exprs:
            if isinstance(expr, tuple):
                expr, rel = expr
            else:
                if expr.is_Relational:
                    expr, rel = expr.lhs - expr.rhs, expr.rel_op
                else:
                    expr, rel = expr, '=='

            poly = Poly(expr, gen)

            if not poly.get_domain().is_Exact:
                poly, exact = poly.to_exact(), False

            domain = poly.get_domain()

            if not (domain.is_ZZ or domain.is_QQ):
                raise NotImplementedError("inequality solving is not supported over %s" % domain)

            _polys.append((poly, rel))

        polys.append(_polys)

    solution = solve_poly_inequalities(polys)

    if isinstance(solution, Union):
        intervals = list(solution.args)
    elif isinstance(solution, Interval):
        intervals = [solution]
    else:
        intervals = []

    if not exact:
        intervals = map(interval_evalf, intervals)

    if not relational:
        return intervals

    real = ask(gen, 'real', assume)

    def relationalize(gen):
        return Or(*[ i.as_relational(gen) for i in intervals ])

    if not real:
        result = And(relationalize(re(gen)), Eq(im(gen), 0))
    else:
        result = relationalize(gen)

    return result

def reduce_abs_inequality(expr, rel, gen, assume=True):
    """Reduce an inequality with nested absolute values. """
    if not ask(gen, 'real', assume):
        raise NotImplementedError("can't solve inequalities with absolute values")

    def bottom_up_scan(expr):
        exprs = []

        if expr.is_Add or expr.is_Mul:
            op = expr.__class__

            for arg in expr.args:
                _exprs = bottom_up_scan(arg)

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
                raise ValueError("only non-negative integer powers are allowed")

            _exprs = bottom_up_scan(expr.base)

            for expr, conds in _exprs:
                exprs.append((expr**n, conds))
        elif isinstance(expr, Abs):
            _exprs = bottom_up_scan(expr.args[0])

            for expr, conds in _exprs:
                exprs.append(( expr, conds + [Ge(expr, 0)]))
                exprs.append((-expr, conds + [Lt(expr, 0)]))
        else:
            exprs = [(expr, [])]

        return exprs

    exprs = bottom_up_scan(expr)

    mapping = {'<': '>', '<=': '>='}
    inequalities = []

    for expr, conds in exprs:
        if rel not in mapping.keys():
            expr = Relational( expr, 0, rel)
        else:
            expr = Relational(-expr, 0, mapping[rel])

        inequalities.append([expr] + conds)

    return reduce_poly_inequalities(inequalities, gen, assume)

def reduce_abs_inequalities(exprs, gen, assume=True):
    """Reduce a system of inequalities with nested absolute values. """
    return And(*[ reduce_abs_inequality(expr, rel, gen, assume) for expr, rel in exprs ])

def reduce_inequalities(inequalities, assume=True):
    """Reduce a system of inequalities with rational coefficients. """
    if not hasattr(inequalities, '__iter__'):
        inequalities = [inequalities]

    poly_part, abs_part, extra_assume = {}, {}, []

    for inequality in inequalities:
        if isinstance(inequality, bool):
            if inequality is False:
                return False
            else:
                continue

        if isinstance(inequality, Assume):
            extra_assume.append(inequality)
            continue

        if inequality.is_Relational:
            expr, rel = inequality.lhs - inequality.rhs, inequality.rel_op
        else:
            expr, rel = inequality, '=='

        gens = expr.atoms(Symbol)

        if not gens:
            return False
        elif len(gens) == 1:
            gen = gens.pop()
        else:
            raise NotImplementedError("only univariate inequalities are supported")

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

    extra_assume = And(*extra_assume)

    if assume is not None:
        assume = And(assume, extra_assume)
    else:
        assume = extra_assume

    poly_reduced = []
    abs_reduced = []

    for gen, exprs in poly_part.iteritems():
        poly_reduced.append(reduce_poly_inequalities([exprs], gen, assume))

    for gen, exprs in abs_part.iteritems():
        abs_reduced.append(reduce_abs_inequalities(exprs, gen, assume))

    return And(*(poly_reduced + abs_reduced))

