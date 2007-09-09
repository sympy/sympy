
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
    from sympy.numerics import Float
    from sympy.numerics.functions import pi_float, sqrt, exp, log, cos

    def frac(x):
        return x - int(x)

    def D(n, j):
        pi = pi_float()
        a = sqrt(Float(2)/3) * pi / j
        b = Float(n) - Float(1)/24
        c = sqrt(b)
        expa = exp(a*c)
        iexpa = Float(1)/expa
        ch = (expa + iexpa)*0.5
        sh = (expa - iexpa)*0.5
        return sqrt(j) / (2*sqrt(2)*b*pi) * (a*ch-sh/c)

    def A(n, j):
        if j == 1:
            return Float(1)
        s = Float(0)
        pi = pi_float()
        for h in xrange(1, j):
            if gcd(h,j) == 1:
                s += cos((g(h,j)-2*h*n)*pi/j)
        return s

    def g(h, j):
        if j < 3:
            return Float(0)
        s = Float(0)
        for k in xrange(1, j):
            s += k*(frac(h*Float(k)/j)-0.5)
        return s

    # estimate number of digits in p(n)
    pdigits = int((pi_float()*sqrt(2.0*n/3)-log(4*n))/log(10)+1)
    Float.store()
    Float.setdps(pdigits*1.1 + 10)
    s = Float(0)

    M = max(6, int(0.24*sqrt(n)+4))
    for q in xrange(1, M):
        s += A(n,q) * D(n,q)
    p = int(s + 0.5)
    Float.revert()

    return p
