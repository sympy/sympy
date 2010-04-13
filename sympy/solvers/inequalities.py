"""Tools for solving inequalities and systems of inequalities. """

from sympy import S, Poly, Interval, Or, Eq

def solve_poly_inequality(inequality):
    """Solve a polynomial inequality with rational coefficients.  """
    poly = Poly(inequality.lhs - inequality.rhs, greedy=False)

    if not poly.is_univariate or not poly.gen.is_real:
        raise ValueError("only real univariate inequalities are supported")

    exact = True

    if not poly.get_domain().is_Exact:
        poly, exact = poly.to_exact(), False

    domain = poly.get_domain()

    if not (domain.is_ZZ or domain.is_QQ):
        raise ValueError("inequality solving is not supported over %s" % domain)

    reals, result = poly.real_roots(multiple=False), []

    if inequality.rel_op == '==':
        for root, _ in reals:
            interval = Interval(root, root)
            result.append(interval)
    elif inequality.rel_op == '!=':
        left = S.NegativeInfinity

        for right, _ in reals + [(S.Infinity, 1)]:
            interval = Interval(left, right, True, True)
            result.append(interval)
            left = right
    else:
        if poly.LC() > 0:
            sign = +1
        else:
            sign = -1

        eq_sign, equal = None, False

        if inequality.rel_op == '>':
            eq_sign = +1

        if inequality.rel_op == '<':
            eq_sign = -1

        if inequality.rel_op == '>=':
            eq_sign, equal = +1, True

        if inequality.rel_op == '<=':
            eq_sign, equal = -1, True

        right, right_open = S.Infinity, True

        for left, multiplicity in reversed(reals):
            if multiplicity % 2:
                if sign == eq_sign:
                    result.insert(0, Interval(left, right, not equal, right_open))

                sign, right, right_open = -sign, left, not equal
            else:
                if sign == eq_sign and not equal:
                    result.insert(0, Interval(left, right, True, right_open))
                    right, right_open = left, True
                elif sign != eq_sign and equal:
                    result.insert(0, Interval(left, left))

        if sign == eq_sign:
            result.insert(0, Interval(S.NegativeInfinity, right, True, right_open))

    return result, exact, poly.gen

def solve_poly_inequalities(inequalities, relational=False):
    """Solve a system of polynomial inequalities with rational coefficients.  """
    results, exact, gen = None, True, None

    if not hasattr(inequalities, '__iter__'):
        inequalities = [inequalities]

    for inequality in inequalities:
        if not inequality.is_Relational:
            inequality = Eq(inequality, S.Zero)

        intervals, _exact, _gen = solve_poly_inequality(inequality)

        if gen is None:
            gen = _gen
        elif gen != _gen:
            raise ValueError("only real univariate inequality systems are supported")

        if results is None:
            results = intervals
        else:
            _results = []

            for interval in intervals:
                for result in results:
                    _interval = interval.intersect(result)

                    if _interval is not S.EmptySet:
                        _results.append(_interval)

            results = _results

        exact &= _exact

        if not results:
            break

    if not exact:
        results = [ Interval(r.left.evalf(), r.right.evalf(),
            left_open=r.left_open, right_open=r.right_open) for r in results ]

    if relational:
        results = Or(*[ r.as_relational(gen) for r in results ])

    return results

