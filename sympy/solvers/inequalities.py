"""Tools for solving inequalities and systems of inequalities. """

from sympy import S, Poly, Interval, And, Or, Eq, DomainError, ask, re, im, Assume

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

        if rel == '<':
            eq_sign = -1

        if rel == '>=':
            eq_sign, equal = +1, True

        if rel == '<=':
            eq_sign, equal = -1, True

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

def solve_poly_inequalities(inequalities, relational=True):
    """Solve a system of polynomial inequalities with rational coefficients.  """
    if not hasattr(inequalities, '__iter__'):
        inequalities = [inequalities]

    polys, exact, assume = {}, {}, []

    for inequality in inequalities:
        if isinstance(inequality, bool):
            if inequality is False:
                return False
            else:
                continue

        if isinstance(inequality, Assume):
            assume.append(inequality)
            continue

        if inequality.is_Relational:
            expr, rel = inequality.lhs - inequality.rhs, inequality.rel_op
        else:
            expr, rel = inequality, '=='

        poly = Poly(expr, greedy=False)

        if not poly.gen.is_Symbol:
            raise NotImplementedError("only polynomial inequalities are supported")

        if not poly.is_univariate:
            raise NotImplementedError("only univariate inequalities are supported")

        _exact = True

        if not poly.get_domain().is_Exact:
            poly, _exact = poly.to_exact(), False

        domain = poly.get_domain()

        if not (domain.is_ZZ or domain.is_QQ):
            raise DomainError("inequality solving is not supported over %s" % domain)

        if poly.gen in polys:
            polys[poly.gen].append((poly, rel))
            exact[poly.gen] &= _exact
        else:
            polys[poly.gen] = [(poly, rel)]
            exact[poly.gen] = _exact

    results = {}

    for gen, polys_group in polys.iteritems():
        global_intervals = None

        for poly, rel in polys_group:
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

        intervals = global_intervals

        if not exact[gen]:
            intervals = [ Interval(i.left.evalf(), i.right.evalf(),
                left_open=i.left_open, right_open=i.right_open) for i in intervals ]

        real = ask(gen, 'real', And(*assume))

        if relational or not real:
            def relationalize(gen):
                return Or(*[ i.as_relational(gen) for i in intervals ])

            if not real:
                result = And(relationalize(re(gen)), Eq(im(gen), 0))
            else:
                result = relationalize(gen)
        else:
            result = intervals

        results[gen] = result

    if relational or not real:
        solution = And(*results.values())
    else:
        if len(results) == 1:
            solution = results.popitem()[1]
        else:
            solution = results

    return solution

