from sympy.concrete.gosper import gosper_normal
from sympy.core import Dummy, Function, symbols, sympify
from sympy.polys import Poly, factor, factor_list, cancel, degree, LC
from sympy.solvers import linsolve, rsolve
from sympy.simplify import hypersimp, combsimp, fraction, denom
from sympy.concrete import product
from sympy.functions import gamma

"""
zb_recur provides an implementation of the algorithm described on pages 106 - 109
of 'A = B' https://www.math.upenn.edu/~wilf/AeqB.pdf, this aids in solving
hypergeometric sums.

zb_sum uses this function to attempt to compute some definite hypergeometric sums.
"""

def _vanishes(F, b, n, k):
    """
    Determines if F(n, k) hypergeometric term vanishes for k > b,
    where k is integer assuming we already know F(n, b + 1) = 0
    by looking at F(n, k + 1) / F(n, k) = P / Q, and considering roots of Q.
    """
    L = []
    for f in factor_list(denom(hypersimp(F, k)))[1]:
        if degree(f[0], k) != 1 or not f[0].free_symbols.issubset({n, k}):
            continue
        L.append((LC(f[0], k) * k - f[0]) / LC(f[0], k))
    if L == [] or max(L) <= b:
        return True
    return False

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

    if (B - A) % p_2.LC() == 0 and (B - A).free_symbols == set():
        return max((B - A) // p_2.LC(), d - d_2 + 1)
    return d - d_2 + 1

def _zb_gosper(F, J, n, k):
    """
    Given a hypergeometric term F of n and k (and possibly other variables),
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

def zb_recur(F, n, k, J = 1):
    """
    Given a hypergeometric term F of n and k (and possibly other variables),
    determines if there  exists a_0, a_1, ... , a_J rational functions in n and a
    function G of n and k such that G / F is rational and

    a_0F(n, k) + a_1F(n + 1, k) + ... + a_JF(n + j, k) = G(n, k + 1) - G(n, k).

    If such things exists this function returns a pair consisting of
    [a_0, a_1, ... , a_J] and the function G / F, otherwise it will return
    None.

    Examples
    ========
    Helping us discover that sum of binomial(n - k - 1, k) over k is fibonacci

    >>> from sympy.concrete.zeilberger import zb_recur
    >>> from sympy.abc import n, k, x
    >>> from sympy import binomial
    >>> F = binomial(n - k - 1, k)
    >>> zb_recur(F, n, k, J = 2)
    ([-1, -1, 1], k*(k - n)/((2*k - n)*(2*k - n - 1)))

    Discovering a recurrence for binomial(n, k)**3, the resulting recurrence for
    the sum can't be solved nicely, we don't care to look at G.

    >>> F = binomial(n, k)**3
    >>> zb_recur(F, n, k, J = 2)[0]
    [-2, -(7*n**2 + 21*n + 16)/(4*(n + 1)**2), (n + 2)**2/(4*(n + 1)**2)]

    F may have other symbols than n and k, here we discover the sum from k = 0 to n to be
    x / (x + n).

    >>> F = (-1)**k * binomial(n, k) / binomial(x + k, k)
    >>> zb_recur(F, n, k)
    ([-n - x, n + x + 1], k*(k + x)/(k - n - 1))
    """
    a, w = symbols('a, w', cls = Function)

    p, p_2, p_3 = _zb_gosper(F, J, n, k)
    # This provides us the polynomials appearing in (6.3.11) of 'A = B', p_3
    # here is p_3(k-1)

    b_deg = _find_b_deg(p_2, p_3, p)

    if b_deg < 0:
        return None

    b = Poly([w(i) for i in range(b_deg, -1, -1)], k)

    P = p_2 * b.shift(1) - p_3 * b - p

    syms = [a(i) for i in range(J + 1)] + [w(i) for i in range(b_deg + 1)]
    solns = list(linsolve(P.coeffs(), syms))[0]
    if all(c == 0 for c in solns):
        return None

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

def zb_sum(F, k_a_b, J = 1):
    """
    Attempts to compute the sum of a hypergometric term F(n, k) from
    k = a (finite) to k = b (finite) by first finding a recurrence of order J in n that
    F(n, k) satisfies. If it fails to find a recurrence or solve the recurrence
    it will return None.

    Examples
    ========

    >>> from sympy.concrete.zeilberger import zb_sum
    >>> from sympy.core import symbols
    >>> n = symbols('n', positive = True, integer = True)
    >>> k, x = symbols('k, x')
    >>> from sympy import binomial
    >>> F = (-1)**k * binomial(x - k + 1, k) * binomial(x - 2 * k, n - k)
    >>> zb_sum(F, (k, 0, n), J = 2)[0]
    (-1)**n/2 + 1/2

    Rediscovering the binomial formula
    >>> F = binomial(n, k) * x**k
    >>> zb_sum(F, (k, 0, n))[0]
    (x + 1)**n

    Here F doesn't vanish for k > n
    >>> F = binomial(n + k, k) / 2**k
    >>> zb_sum(F, (k, 0, n))[0]
    2**n

    Here we need a large J
    >>> F = binomial(n, 3 * k)
    >>> zb_sum(F, (k, 0, n), J = 3)[0]
    2**n/3 + (1/2 - sqrt(3)*I/2)**n/3 + (1/2 + sqrt(3)*I/2)**n/3
    """
    k, a, b = k_a_b
    a, b = sympify(a), sympify(b)

    F = F.rewrite(gamma)

    for w in F.free_symbols.intersection((b - a).free_symbols):
        if  not (w.is_positive is False or w.is_integer is False) and not (
            hypersimp(F, k) is None or hypersimp(F, w) is None):
            break
    else:
        return None

    n = symbols('n', positive = True, integer = True)
    F = F.subs(w, n)
    a, b = a.subs(w, n), b.subs(w, n)

    pair = zb_recur(F, n, k, J)

    if pair is None:
        return None

    F_rec, G = pair[0], F * pair[1]

    G_a = G.subs(k, a)

    f = symbols('f', cls = Function)
    sum_rec = sum(F_rec[i] * f(n + i) for i in range(J + 1))

    i = symbols('i', cls = Dummy)
    if not (b.subs(n, i) - a.subs(n, n + i)).free_symbols.issubset({i}):
        return None

    initial = { f(i): sum(F.subs([(n, i), (k, a.subs(n, i) + j)])
            for j in range(b.subs(n, i) - a.subs(n, i) + 1)) for i in range(1, J + 2) }

    vanishes = (combsimp(F.subs(k, b + 1)) == 0) and _vanishes(F, b, n, k)
    if vanishes:
        try:
            res = rsolve(combsimp(sum_rec + G_a), f(n), initial)
            if res == None:
                return None
            return res.subs(n, w), w
        except ValueError:
            return None

    if not (b.subs(n, n + J) - b.subs(n, n + i)).free_symbols.issubset({i}):
        return None

    boundary = sum(sum(F.subs([(n, n + i), (k, b.subs(n, n + i) + 1 + j)])
            for j in range(b.subs(n, n + J) - b.subs(n, n + i))) for i in range(J + 1))

    G_b = G.subs(k, b.subs(n, n + J))
    try:
        res = rsolve(combsimp(sum_rec + boundary + G_a - G_b), f(n), initial)
        if res == None:
            return None
        return res.subs(n, w), w
    except ValueError:
        return None
