"""Tools for solving inequalities and systems of inequalities. """

from sympy import S, Poly, Interval, Or

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

    reals, result = poly.real_roots(multiple=False), []

    if eq.rel_op == '==':
        for root, _ in reals:
            interval = Interval(root, root)
            result.append(interval)
    elif eq.rel_op == '!=':
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

        if eq.rel_op == '>':
            eq_sign = +1

        if eq.rel_op == '<':
            eq_sign = -1

        if eq.rel_op == '>=':
            eq_sign, equal = +1, True

        if eq.rel_op == '<=':
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

    if not exact:
        result = [ r.evalf() for r in result ]

    if not relational:
        return result
    else:
        return Or(*[ interval.as_relational(poly.gen) for interval in result ])
