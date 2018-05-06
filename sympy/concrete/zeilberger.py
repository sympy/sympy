from sympy.concrete.gosper import gosper_normal
from sympy import Poly, degree, factor, fraction, cancel, LC, linsolve
from sympy import symbols, Dummy, Function, combsimp, hypersimp, product
"""
zb_recur provides an implementation of the algorithm described on pages 106 - 109 
of 'A = B' https://www.math.upenn.edu/~wilf/AeqB.pdf, this aids in solving 
sums of hypergeometric functions.
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
    determines if there  exists a_0, a_1, ... , a_J polynomials in n and a 
    function G of n and k such that G / F is rational and 

    a_0F(n, k) + a_1F(n + 1, k) + ... + a_JF(n + j, k) = G(n, k + 1) - G(n, k).

    If such things exists this function returns a pair consisting of 
    [a_0, a_1, ... , a_J] and the function G / F, otherwise it returns 
    'try higher order', suggesting you try to find a higher order recurrence, 
    i.e. increase J. If F is 'proper hypergeometric' (see page 64 of 'A = B' 
    for a definition) then such a recurrence always exists for some J.

    The main application of this is to sum k over some suitable values so that
    the G telescopes and you obtain a recurrence for the sum of F(n, k).

    Examples
    ========
    Helping us discover that sum of binomial(n - k - 1, k) over k is fibonacci

    >>> F = binomial(n - k - 1, k)
    >>> zb_recur(F, 1, n, k)
    try higher order
    >>> zb_recur(F, 2, n, k)
    ([-1, -1, 1], k*(k - n)/((2*k - n)*(2*k - n - 1)))

    Discovering a recurrence for binomial(n, k)**3, the resulting recurrence for
    the sum can't be solved nicely, we don't care to look at G.

    >>> F = binomial(n, k)**3
    >>> zb_recur(F, 2, n, k)[0]
    ([-2, -(7*n**2 + 21*n + 16)/(4*(n + 1)**2), (n + 2)**2/(4*(n + 1)**2)]

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
        return "try higher order"

    b = Poly(sum(w(i) * k**i for i in range(b_deg + 1)), k)

    P = p_2 * b.shift(1) - p_3 * b - p

    syms = [a(i) for i in range(J + 1)] + [w(i) for i in range(b_deg + 1)]

    solns = list(linsolve(P.coeffs(), syms))[0]

    if all(c == 0 for c in solns):
        return "try higher order"

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

"""
Tests
========
from sympy import factorial, binomial, simplify

n, k, x, y, i, j = symbols('n, k, x, y, i, j')

def test_zb_recur():

    F_1 = (-1)**k * binomial(x - k + 1, k) * binomial(x - 2*k, n - k)
    F_2 = binomial(n, k) / 2**n
    F_3 = (-1)**k * binomial(n, k) * x * binomial(x + n, n) / (k + x)
    F_4 = binomial(n - k - 1, k)
    F_5 = (factorial(n - i) * factorial(n - j) * factorial(i - 1) * 
        factorial(j - 1)) / (factorial(n - 1) * factorial(k - 1) * 
        factorial(n - i - j + k) * factorial(i - k) * factorial(j - k))

    R_1 = k*(-2*k + x + 1)*(-k + x + 2)/(2*(k - n - 2)*(k - n - 1))
    R_2 = k/(2*(k - n - 1))
    R_3 = k*(k + x)/(k - n - 1)
    R_4 = k*(k - n)/((2*k - n)*(2*k - n - 1))
    R_5 = -(k - 1)/n

    assert simplify(zb_recur(F_1, 2, n, k)[1]) == simplify(R_1)
    assert simplify(zb_recur(F_2, 1, n, k)[1]) == simplify(R_2)
    assert simplify(zb_recur(F_3, 1, n, k)[1]) == simplify(R_3)
    assert zb_recur(F_4, 1, n, k) == "try higher order"
    assert simplify(zb_recur(F_4, 2, n, k)[1]) == simplify(R_4)
    assert simplify(zb_recur(F_5, 1, n, k)[1]) == simplify(R_5)

test_zb_recur()
"""
