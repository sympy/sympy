
def npartitions(n):
    """
    Calculate the partition function P(n), i.e. the number of ways that
    n can be written as a sum of positive integers.

    P(n) is computed using a straightforward implementation of the
    Hardy-Ramanujan-Rademacher formula, described e.g. at
    http://mathworld.wolfram.com/PartitionFunctionP.html

    The speed is decent up to n = 10**5 or so. The function has
    been tested to give the correct result for n = 10**6.
    """

    n = int(n)
    if n < 0:
        return 0
    if n <= 5:
        return [1, 1, 2, 3, 5, 7][n]

    from sympy.core.numbers import gcd
    from sympy.thirdparty.mpmath import mp, mpf, pi, sqrt, exp, log, cos

    def frac(x):
        return x - int(x)

    def D(n, j):
        a = sqrt(mpf(2)/3) * pi / j
        b = mpf(n) - mpf(1)/24
        c = sqrt(b)
        expa = exp(a*c)
        iexpa = mpf(1)/expa
        ch = (expa + iexpa)*0.5
        sh = (expa - iexpa)*0.5
        return sqrt(j) / (2*sqrt(2)*b*pi) * (a*ch-sh/c)

    def A(n, j):
        if j == 1:
            return mpf(1)
        s = mpf(0)
        for h in xrange(1, j):
            if gcd(h,j) == 1:
                s += cos((g(h,j)-2*h*n)*pi/j)
        return s

    def g(h, j):
        if j < 3:
            return mpf(0)
        s = mpf(0)
        for k in xrange(1, j):
            s += k*(frac(h*mpf(k)/j)-0.5)
        return s

    # estimate number of digits in p(n)
    pdigits = int((pi*sqrt(2.0*n/3)-log(4*n))/log(10)+1)

    dps, mp.dps = mp.dps, pdigits*1.1 + 10

    s = mpf(0)

    M = max(6, int(0.24*sqrt(n)+4))
    for q in xrange(1, M):
        s += A(n,q) * D(n,q)
    p = int(s + 0.5)
    mp.dps = dps

    return p
