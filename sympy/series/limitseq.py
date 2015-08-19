"""Limits of sequences"""

from __future__ import print_function, division

from sympy.core.sympify import sympify
from sympy.core.singleton import S
from sympy.core.add import Add
from sympy.core.function import PoleError
from sympy.series.limits import limit


def difference_delta(expr, n=None, step=1):
    """Difference Operator.

    Discrete analogous to differential operator.

    Examples
    ========

    >>> from sympy import difference_delta as dd
    >>> from sympy.abc import n
    >>> dd(n*(n + 1), n)
    2*n + 2
    >>> dd(n*(n + 1), n, 2)
    4*n + 6

    References
    ==========

    .. [1] https://reference.wolfram.com/language/ref/DifferenceDelta.html
    """
    expr = sympify(expr)

    if n is None:
        f = expr.free_symbols
        if len(f) == 1:
            n = f.pop()
        elif len(f) == 0:
            return S.Zero
        else:
            raise ValueError("Since there is more than one variable in the"
                             " expression, a variable must be supplied to"
                             " take the difference of %s" % expr)
    step = sympify(step)
    if step.is_number is False:
        raise ValueError("Step should be a number.")
    elif step in [S.Infinity, -S.Infinity]:
        raise ValueError("Step should be bounded.")

    if hasattr(expr, '_eval_difference_delta'):
        result = expr._eval_difference_delta(n, step)
        if result:
            return result

    return expr.subs(n, n + step) - expr


def dominant(expr, n):
    """Finds the most dominating term in an expression.

    if limit(a/b, n, oo) is oo then a dominates b.
    if limit(a/b, n, oo) is 0 then b dominates a.
    else a and b are comparable.

    returns the most dominant term.
    If no unique domiant term, then returns ``None``.

    Examples
    ========

    >>> from sympy import Sum
    >>> from sympy.series.limitseq import dominant
    >>> from sympy.abc import n, k
    >>> dominant(5*n**3 + 4*n**2 + n + 1, n)
    5*n**3
    >>> dominant(2**n + Sum(k, (k, 0, n)), n)
    2**n

    See Also
    ========

    sympy.series.limitseq.dominant
    """
    terms = Add.make_args(expr.expand(func=True))
    term0 = terms[-1]
    comp = [term0]  # comparable terms
    for t in terms[:-1]:
        e = (term0 / t).combsimp()
        l = limit_seq(e, n)
        if l is S.Zero:
            term0 = t
            comp = [term0]
        elif l is None:
            return None
        elif l not in [S.Infinity, -S.Infinity]:
            comp.append(t)
    if len(comp) > 1:
        return None
    return term0


def _limit_inf(expr, n):
    try:
        return limit(expr, n, S.Infinity)
    except (NotImplementedError, PoleError):
        return None


def limit_seq(expr, n, trials=5):
    """Finds limits of terms having sequences at infinity.

    Admissible Terms
    ================

    The terms should be built from rational functions, indefinite sums,
    and indefinite products over an indeterminate n. A term is admissible
    if the scope of all product quantifiers are asymptotically positive.
    Every admissible term is asymptoticically monotonous.

    Examples
    ========

    >>> from sympy import limit_seq, Sum, binomial
    >>> from sympy.abc import n, k, m
    >>> limit_seq((5*n**3 + 3*n**2 + 4) / (3*n**3 + 4*n - 5), n)
    5/3
    >>> limit_seq(binomial(2*n, n) / Sum(binomial(2*k, k), (k, 1, n)), n)
    3/4
    >>> limit_seq(Sum(k**2 * Sum(2**m/m, (m, 1, k)), (k, 1, n)) / (2**n*n), n)
    4

    See Also
    ========

    sympy.series.limitseq.dominant

    References
    ==========

    .. [1] Computing Limits of Sequences - Manuel Kauers
    """
    from sympy.concrete.summations import Sum

    for i in range(trials):
        if not expr.has(Sum):
            result = _limit_inf(expr, n)
            if result is not None:
                return result

        num, den = expr.as_numer_denom()
        if not den.has(n):
            result = _limit_inf(expr.doit(), n)
            if result is not None:
                return result
            return None

        num, den = map(lambda t: difference_delta(t.expand(), n), [num, den])

        expr = (num / den).combsimp()

        if not expr.has(Sum):
            result = _limit_inf(expr, n)
            if result is not None:
                return result

        num, den = expr.as_numer_denom()

        num = dominant(num, n)
        if num is None:
            return None

        den = dominant(den, n)
        if den is None:
            return None

        expr = (num / den).combsimp()
