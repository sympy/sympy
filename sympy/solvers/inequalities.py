"""Tools for solving inequalities and systems of inequalities. """

from sympy import S, Poly, RootOf, intervals, Interval, Or

def solve_poly_inequality(eq, relational=False):
    """Solve a polynomial inequality with rational coefficients.  """
    poly = Poly(eq.lhs - eq.rhs, greedy=False)

    if not poly.is_univariate or not poly.gen.is_real:
        raise ValueError("only real univariate inequalities are supported")

    exact = True

    if not poly.get_domain().is_Exact:
        poly, exact = poly.to_exact(), False

    domain = poly.get_domain()

    if not (domain.is_ZZ or domain.is_QQ):
        raise ValueError("inequality solving is not supported over %s" % domain)

    result = []

    if eq.rel_op in ['==', '!=']:
        n, roots = poly.count_roots(), []

        for i in xrange(0, n):
            roots.append(RootOf(poly, i))

        if eq.rel_op == '==':
            for root in roots:
                interval = Interval(root, root)
                result.append(interval)
        else:
            left, roots = S.NegativeInfinity, roots + [S.Infinity]

            for right in roots:
                interval = Interval(left, right, True, True)
                result.append(interval)
                left = right
    else:
        _, factors = poly.factor_list()

        reals = intervals([ factor for factor, _ in factors])

        visited, roots = {}, []

        for interval, indices in reals:
            (ith,) = indices.keys()

            if ith in visited:
                index = visited[ith] + 1
            else:
                index = 0

            roots.append(RootOf(factors[ith][0], index))

        if poly.LC() > 0:
            sign = +1
        else:
            sign = -1

        eq_sign, equal = None, False

        if eq.rel_op == '>':
            eq_sign = +1

        if eq.rel_op == '<':
            eq_sign = -1

        if eq.rel_op == '>=':
            eq_sign, equal = +1, True

        if eq.rel_op == '<=':
            eq_sign, equal = -1, True

        right, right_open = S.Infinity, True

        for (interval, indices), left in reversed(zip(reals, roots)):
            (ith,) = indices.keys()

            if factors[ith][1] % 2:
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

    if not exact:
        result = [ r.evalf() for r in result ]

    if not relational:
        return result
    else:
        return Or(*[ interval.as_relational(poly.gen) for interval in result ])
