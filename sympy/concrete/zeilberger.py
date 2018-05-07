from sympy.concrete.gosper import gosper_normal
from sympy import Poly, degree, factor, fraction, cancel, LC, linsolve, rsolve
from sympy import symbols, Dummy, Function, combsimp, hypersimp, product
"""
zb_recur provides an implementation of the algorithm described on pages 106 - 109
of 'A = B' https://www.math.upenn.edu/~wilf/AeqB.pdf, this aids in solving
sums of hypergeometric functions.

zb_sum uses this function to attempt to compute some definite hypergeometric sum.
"""

def _find_b_deg(p_2, p_3, p):
    """
    Given polynomials p_2, p_3, p in k, finds a bound on the degree of a
    polynomial b in k such that p_2(k)b(k+1)-p_3(k)b(k) = p(k).
    """
    d_2 = p_2.degree()
    d_3 = p_3.degree()
    d = p.degree()

    if d_2 == 0 and d_3 == 0:
        return d + 1

    if d_2 != d_3 or p_2.LC() != p_3.LC():
        return d - max(d_2, d_3)

    A = p_2.nth(d_2 - 1)
    B = p_3.nth(d_3 - 1)

    if (B - A) % p_2.LC() == 0:
        return max((B - A) // p_2.LC(), d - d_2 + 1)

    return d - d_2 + 1

def _zb_gosper(F, J, n, k):
    """
    Given a hypergeometric function F of n and k (and possibly other variables),
    we take t_k = a_0F(n, k) + a_1F(n + 1, k) + ... + a_JF(n + J, k) with a_i
    some variables and determine the gosper normal form of t_(k+1)/t_k.
    """
    i = symbols('i', cls = Dummy)
    a = symbols('a', cls = Function)

    r_1, r_2 = fraction(hypersimp(F, k))
    s_1, s_2 = fraction(hypersimp(F.subs(n, n-1), n))

    p_0 = Poly(sum(a(j) * (product(s_1.subs(n, n + j - i), (i, 0, j-1)) *
                    product(s_2.subs(n, n + i), (i, j + 1, J)))
                    for j in range(J + 1)), k)

    # We can factor out a p_0(k+1)/p_0(k) from the gosper
    # normal form calculation

    r = r_1 * product(s_2.subs(n, n + i), (i, 1, J))
    s = r_2 * product(s_2.subs([(n, n + i), (k, k + 1)]), (i, 1, J))
    r, s = fraction(cancel(r / s))

    p_2, p_3, p_1 = gosper_normal(r, s, k)

    return p_0 * p_1, p_2, p_3.shift(-1)

def zb_recur(F, J, n, k):
    """
    Given a hypergeometric function F of n and k (and possibly other variables),
    determines if there  exists a_0, a_1, ... , a_J rational functions in n and a
    function G of n and k such that G / F is rational and

    a_0F(n, k) + a_1F(n + 1, k) + ... + a_JF(n + j, k) = G(n, k + 1) - G(n, k).

    If such things exists this function returns a pair consisting of
    [a_0, a_1, ... , a_J] and the function G / F, otherwise it raised a value
    error 'try higher order', suggesting you try to find a higher order recurrence,
    i.e. increase J. If F is 'proper hypergeometric' (see page 64 of 'A = B'
    for a definition) then such a recurrence always exists for some J.

    The main application of this is to sum k over some suitable values so that
    the G telescopes and you obtain a recurrence for the sum of F(n, k).

    Examples
    ========
    Helping us discover that sum of binomial(n - k - 1, k) over k is fibonacci

    >>> from sympy.concrete.zeilberger import zb_recur
    >>> from sympy.abc import n, k, x
    >>> from sympy import binomial
    >>> F = binomial(n - k - 1, k)
    >>> zb_recur(F, 2, n, k)
    ([-1, -1, 1], k*(k - n)/((2*k - n)*(2*k - n - 1)))

    Discovering a recurrence for binomial(n, k)**3, the resulting recurrence for
    the sum can't be solved nicely, we don't care to look at G.

    >>> F = binomial(n, k)**3
    >>> zb_recur(F, 2, n, k)[0]
    [-2, -(7*n**2 + 21*n + 16)/(4*(n + 1)**2), (n + 2)**2/(4*(n + 1)**2)]

    F may have other symbols than n and k, here we discover the sum from k = 0 to n to be
    x / (x + n).

    >>> F = (-1)**k * binomial(n, k) / binomial(x + k, k)
    >>> zb_recur(F, 1, n, k)
    ([-n - x, n + x + 1], k*(k + x)/(k - n - 1))
    """
    a, w = symbols('a, w', cls = Function)

    p, p_2, p_3 = _zb_gosper(F, J, n, k)
    # This provides us the polynomials appearing in (6.3.11) of 'A = B', p_3
    # here is p_3(k-1)

    b_deg = _find_b_deg(p_2, p_3, p)

    if b_deg < 0:
        raise ValueError("try higher order")

    b = Poly([w(i) for i in range(b_deg, -1, -1)], k)

    P = p_2 * b.shift(1) - p_3 * b - p

    syms = [a(i) for i in range(J + 1)] + [w(i) for i in range(b_deg + 1)]

    solns = list(linsolve(P.coeffs(), syms))[0]

    if all(c == 0 for c in solns):
        raise ValueError("try higher order")

    save = a(0)

    for c in solns:
        for s in syms:
            if c.coeff(s):
                save = s
                break
        else:
            continue
        break

    sub_pairs = [(s, 1) if s == save else (s, 0) for s in syms]
    coeffs = list(zip(syms, factor(solns.subs(sub_pairs))))

    t = sum(a(i) * F.subs(n, n + i) for i in range(J + 1))

    G = combsimp((p_3 * b * t / p).subs(coeffs))

    return [c[1] for c in coeffs[0:J + 1]], combsimp(G / F)

def zb_sum(F, U, J, n, k, vanishes = True):
    """
    Attempts to compute the sum of a hypergometric function F(n, k) from
    k = 0 to k = U by first finding a recurrence of order J in n that
    F(n, k) satisfies. If it fails to find a recurrence it will raise a
    value error - "try higher order", suggesting you try again with a
    larger J.

    By default it is assumed that F(n, k) is zero for k < 0 and
    k > U. To tell zb_sum that F(n, k) does not vanish for
    k > U you should use the keyword vanishes = False

    Examples
    ========

    >>> from sympy.concrete.zeilberger import zb_sum
    >>> from sympy.abc import n, k, x
    >>> from sympy import binomial
    >>> F = (-1)**k * binomial(x - k + 1, k) * binomial(x - 2 * k, n - k)
    >>> zb_sum(F, n, 2, n, k)
    (-1)**n/2 + 1/2

    Here F doesn't vanish for k > n
    >>> F = binomial(n + k, k) / 2**k
    >>> zb_sum(F, n, 1, n, k, vanishes = False)
    2**n
    """
    f = symbols('f', cls = Function)

    F_rec, R = zb_recur(F, J, n, k)

    sum_rec = sum(F_rec[i] * f(n + i) for i in range(J + 1))

    initial = { f(i): sum(F.subs([(n, i), (k,  j)])
        for j in range(U.subs(n, i) + 1)) for i in range(1, J + 1) }

    if vanishes:
        return combsimp(rsolve(sum_rec, f(n), initial))

    G = F * R

    boundary = sum(F.subs([(n, n + i), (k, U + j)])
                    for i in range(J + 1) for j in range(i + 1, J + 1))

    G_0, G_U = G.subs(k, 0), G.subs(k, U + J)

    rec = sum_rec + boundary + G_0 - G_U
    # Your sum satisfies the reccurence rec = 0

    return combsimp(rsolve(combsimp(rec), f(n), initial))
