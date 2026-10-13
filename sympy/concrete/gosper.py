"""Gosper's algorithm for hypergeometric summation. """
from __future__ import annotations

from sympy.core import S, Dummy, Mul, symbols, sympify
from sympy.polys import Poly, parallel_poly_from_expr, factor, factor_list
from sympy.utilities.iterables import is_sequence


def gosper_normal(f, g, n, polys=True):
    r"""
    Compute the Gosper's normal form of ``f`` and ``g``.

    Explanation
    ===========

    Given relatively prime univariate polynomials ``f`` and ``g``,
    rewrite their quotient to a normal form defined as follows:

    .. math::
        \frac{f(n)}{g(n)} = Z \cdot \frac{A(n) C(n+1)}{B(n) C(n)}

    where ``Z`` is an arbitrary constant and ``A``, ``B``, ``C`` are
    monic polynomials in ``n`` with the following properties:

    1. `\gcd(A(n), B(n+h)) = 1 \forall h \in \mathbb{N}`
    2. `\gcd(B(n), C(n+1)) = 1`
    3. `\gcd(A(n), C(n)) = 1`

    This normal form, or rational factorization in other words, is a
    crucial step in Gosper's algorithm and in solving of difference
    equations. It can be also used to decide if two hypergeometric
    terms are similar or not.

    This procedure will return a tuple containing elements of this
    factorization in the form ``(Z*A, B, C)``.

    Examples
    ========

    >>> from sympy.concrete.gosper import gosper_normal
    >>> from sympy.abc import n

    >>> gosper_normal(4*n+5, 2*(4*n+1)*(2*n+3), n, polys=False)
    (1/4, n + 3/2, n + 1/4)

    """
    (p, q), opt = parallel_poly_from_expr(
        (f, g), n, field=True, extension=True)

    a, A = p.LC(), p.monic()
    b, B = q.LC(), q.monic()

    C, Z = A.one, a/b
    h = Dummy('h')

    D = Poly(n + h, n, h, domain=opt.domain)

    R = A.resultant(B.compose(D))
    roots = {r for r in R.ground_roots().keys() if r.is_Integer and r >= 0}
    for i in sorted(roots):
        d = A.gcd(B.shift(+i))

        A = A.quo(d)
        B = B.quo(d.shift(-i))

        for j in range(1, i + 1):
            C *= d.shift(-j)

    A = A.mul_ground(Z)

    if not polys:
        A = A.as_expr()
        B = B.as_expr()
        C = C.as_expr()

    return A, B, C


def gosper_term(f, n):
    r"""
    Compute Gosper's hypergeometric term for ``f``.

    Explanation
    ===========

    Suppose ``f`` is a hypergeometric term such that:

    .. math::
        s_n = \sum_{k=0}^{n-1} f_k

    and `f_k` does not depend on `n`. Returns a hypergeometric
    term `g_n` such that `g_{n+1} - g_n = f_n`.

    Examples
    ========

    >>> from sympy.concrete.gosper import gosper_term
    >>> from sympy import factorial
    >>> from sympy.abc import n

    >>> gosper_term((4*n + 1)*factorial(n)/factorial(2*n + 1), n)
    (-n - 1/2)/(n + 1/4)

    """
    from sympy.simplify import hypersimp
    r = hypersimp(f, n)

    if r is None:
        return None    # 'f' is *not* a hypergeometric term

    p, q = r.as_numer_denom()

    A, B, C = gosper_normal(p, q, n)
    B = B.shift(-1)

    N = A.degree()
    M = B.degree()
    K = C.degree()

    if (N != M) or (A.LC() != B.LC()):
        d = K - max(N, M)
    elif not N:
        d = K + 1
    else:
        d1 = K - N + 1
        d2 = (B.nth(N - 1) - A.nth(N - 1))/A.LC()
        if not d2.is_Integer:
            # could be symbolic or Rational
            d2 = -1 # invalid degree candidate
        d = max(d1, d2)

    if d < 0:
        return None    # 'f(n)' is *not* Gosper-summable

    coeffs = symbols('c:%s' % (d + 1), cls=Dummy)
    domain = A.get_domain().inject(*coeffs)

    x = Poly(coeffs, n, domain=domain)
    H = A*x.shift(1) - B*x - C

    from sympy.solvers.solvers import solve
    solution = solve(H.coeffs(), coeffs)

    if solution is None:
        return None    # 'f(n)' is *not* Gosper-summable

    x = x.as_expr().subs(solution)

    for coeff in coeffs:
        if coeff not in solution:
            x = x.subs(coeff, 0)

    if x.is_zero:
        return None    # 'f(n)' is *not* Gosper-summable
    else:
        return B.as_expr()*x/C.as_expr()


def _denominator_factors(expr):
    """Return the irreducible factors of the denominator of ``expr``."""
    from sympy.polys.rationaltools import together

    den = together(expr).as_numer_denom()[1]
    if den is S.One:
        return set()
    if not den.is_polynomial(*den.free_symbols):
        # A denominator that is not a polynomial in its symbols (it
        # can e.g. contain a Piecewise) is not factored; its
        # multiplicative parts are used as they are.
        return set(Mul.make_args(den))
    return {fac for fac, _exp in factor_list(den)[1]}


def _sum_at_degenerate_value(f, k, a, b):
    """Evaluate a definite sum whose parameters were specialised.

    The sum is evaluated directly when the specialised bounds are
    concrete integers, and with ``gosper_sum`` otherwise. ``None`` is
    returned when neither gives a value.
    """
    if a.is_Integer and b.is_Integer and 0 <= b - a <= 100:
        total = S.Zero
        for j in range(int(a), int(b) + 1):
            term = f.subs(k, j)
            if term.has(S.NaN, S.Infinity, -S.Infinity, S.ComplexInfinity):
                return None
            total += term
        return total
    try:
        return gosper_sum(f, (k, a, b))
    except Exception:
        return None


def _check_degenerate_parameters(f, g, result, k, a, b):
    """Check a definite Gosper sum at the poles of its certificate.

    The certificate ``g`` computed by Gosper's algorithm is derived
    under the assumption that the parameters of ``f`` avoid the zeros
    of the denominator of ``g``. If such a factor is still present in
    the denominator of ``result`` then the result makes no claim at
    the pole and is returned unchanged. If the factor has cancelled
    out of ``result`` then the result can silently be wrong at the
    pole (see issue 29731), so the sum is recomputed with the
    parameter replaced by each root of the factor that can be found.
    If any of those values disagrees with ``result`` then ``None``
    is returned instead of the incorrect closed form.
    """
    from sympy.simplify import simplify
    from sympy.solvers.solvers import solve

    params = sympify(f).free_symbols - {k}
    if not params:
        return result

    poles = {fac for fac in _denominator_factors(g)
             if fac.free_symbols & params}
    poles -= _denominator_factors(result)
    if not poles:
        return result

    for pole in poles:
        for param in pole.free_symbols & params:
            if not pole.is_polynomial(param):
                continue
            for solution in solve(pole, param, dict=True):
                root = solution[param]
                if param in root.free_symbols:
                    continue
                subs = {param: root}
                value = _sum_at_degenerate_value(
                    f.subs(subs), k, sympify(a).subs(subs),
                    sympify(b).subs(subs))
                if value is None:
                    continue
                claimed = result.subs(subs)
                if claimed.has(S.NaN, S.Infinity, -S.Infinity,
                               S.ComplexInfinity):
                    continue
                if simplify(value - claimed) != 0:
                    return None
    return result


def gosper_sum(f, k):
    r"""
    Gosper's hypergeometric summation algorithm.

    Explanation
    ===========

    Given a hypergeometric term ``f`` such that:

    .. math ::
        s_n = \sum_{k=0}^{n-1} f_k

    and `f(n)` does not depend on `n`, returns `g_{n} - g(0)` where
    `g_{n+1} - g_n = f_n`, or ``None`` if `s_n` cannot be expressed
    in closed form as a sum of hypergeometric terms.

    Examples
    ========

    >>> from sympy.concrete.gosper import gosper_sum
    >>> from sympy import factorial
    >>> from sympy.abc import n, k

    >>> f = (4*k + 1)*factorial(k)/factorial(2*k + 1)
    >>> gosper_sum(f, (k, 0, n))
    (-factorial(n) + 2*factorial(2*n + 1))/factorial(2*n + 1)
    >>> _.subs(n, 2) == sum(f.subs(k, i) for i in [0, 1, 2])
    True
    >>> gosper_sum(f, (k, 3, n))
    (-60*factorial(n) + factorial(2*n + 1))/(60*factorial(2*n + 1))
    >>> _.subs(n, 5) == sum(f.subs(k, i) for i in [3, 4, 5])
    True

    References
    ==========

    .. [1] Marko Petkovsek, Herbert S. Wilf, Doron Zeilberger, A = B,
           AK Peters, Ltd., Wellesley, MA, USA, 1997, pp. 73--100

    """
    indefinite = False

    if is_sequence(k):
        k, a, b = k
    else:
        indefinite = True

    g = gosper_term(f, k)

    if g is None:
        return None

    if indefinite:
        result = f*g
    else:
        result = (f*(g + 1)).subs(k, b) - (f*g).subs(k, a)

        if result is S.NaN:
            try:
                result = (f*(g + 1)).limit(k, b) - (f*g).limit(k, a)
            except NotImplementedError:
                result = None

        if result is not None:
            result = _check_degenerate_parameters(f, g, result, k, a, b)

    if result is None:
        return None

    return factor(result)
