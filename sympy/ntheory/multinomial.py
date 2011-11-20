from sympy.mpmath.libmp import bitcount
from sympy.core.compatibility import combinations_with_replacement
from collections import defaultdict

def binomial_coefficients(n):
    """Return a dictionary containing pairs {(k1,k2) : C_kn} where
    C_kn are binomial coefficients and n=k1+k2."""
    d = {(0, n):1, (n, 0):1}
    a = 1
    for k in xrange(1, n//2+1):
        a = (a * (n-k+1))//k
        d[k, n-k] = d[n-k, k] = a
    return d

def binomial_coefficients_list(n):
    """ Return a list of binomial coefficients as rows of the Pascal's
    triangle.
    """
    d = [1] * (n+1)
    a = 1
    for k in xrange(1, n//2+1):
        a = (a * (n-k+1))//k
        d[k] = d[n-k] = a
    return d

def _code_t(t, bits_exp):
    """encode a tuple in an integer to do fast addition
    """
    expv = 0
    for i, n in enumerate(t):
        expv += n<<(i*bits_exp)
    return expv

def _bmpnl_calc(m, n, _tuple):
    bits_exp = bitcount(n)
    mask_exp = (1 << bits_exp) - 1
    def t(i):
        rv = [0]*m
        if i==0: return rv
        rv[0]=-1
        rv[i]=1
        return rv
    p0 = [_code_t(_tuple(t(i)), bits_exp) for i in xrange(m)]
    n0 = [0]*m
    n0[0] = n
    l = [0]*(n*(m - 1) + 1)
    l[0] = [(_code_t(_tuple(n0), bits_exp), 1)]
    return bits_exp, mask_exp, p0, n0, l

def _dcalc(m, n, k, p0, l):
    d = defaultdict(int)
    for i in xrange(1, min(m, k+1)):
        nn = (n + 1)*i - k
        if not nn:
            continue
        t = p0[i]
        for t2, c2 in l[k - i]:
            if not c2:
                continue
            tt = t2 + t
            d[tt] += nn * c2
    return d

def multinomial_coefficients(m, n, _tuple=tuple):
    """Return a dictionary containing pairs ``{(k1,k2,..,km) : C_kn}``
    where ``C_kn`` are multinomial coefficients such that
    ``n=k1+k2+..+km``.

    For example:

    >>> from sympy.ntheory import multinomial_coefficients
    >>> multinomial_coefficients(2, 5)
    {(0, 5): 1, (1, 4): 5, (2, 3): 10, (3, 2): 10, (4, 1): 5, (5, 0): 1}

    The algorithm is based on the following result:

       Consider a polynomial and its ``n``-th exponent::

         P(x) = sum_{i=0}^m p_i x^i
         P(x)^n = sum_{k=0}^{m n} a(n,k) x^k

       The coefficients ``a(n,k)`` can be computed using the
       J.C.P. Miller Pure Recurrence [see D.E.Knuth, Seminumerical
       Algorithms, The art of Computer Programming v.2, Addison
       Wesley, Reading, 1981;]::

         a(n,k) = 1/(k p_0) sum_{i=1}^m p_i ((n+1)i-k) a(n,k-i),

       where ``a(n,0) = p_0^n``.

    The following optimizations have been made:
    i) monomial tuples are packed in integers to speed up their sums;
    ii) for `m` large with respect to `n` the monomial tuples `t` have
    many zeroes, and have the same coefficient as
    in `monomial_coefficients(n,n)` for `t` stripped of its zeroes;
    therefore by precomputing the latter coefficients memory and time
    are saved.
    """
    if m == 2:
        return binomial_coefficients(n)
    if m >= 2*n and n > 1:
        return dict(multinomial_coefficients_iterator(m, n))
    if m == 3 and n == 2:
        return {(0, 0, 2): 1, (0, 1, 1): 2, (0, 2, 0): 1, (1, 0, 1): 2,
                (1, 1, 0): 2, (2, 0, 0): 1}
    # The monomial tuples have entries between 0 and n;
    # in the algorithm new monomial tuples are obtained summing
    # them with tuples in `p0`, with the form (0,..,-1,0..,0,1,0...);
    # in these sums the sums and differences of the entries are
    # between 0 and n, so that one can pack a tuple in an int,
    # giving bitcount(n) bits for each entry; in this way sums of
    # tuples become sums of ints, which are fast.
    bits_exp, mask_exp, p0, n0, l = _bmpnl_calc(m, n, _tuple)
    r = {_tuple(n0):1}
    for k in xrange(1, n*(m-1)+1):
        d = _dcalc(m, n, k, p0, l)
        b = []
        for t, c in d.iteritems():
            if not c:
                continue
            c //= k
            if (t&mask_exp) > 0:
                b.append((t, c))
            # decode t
            a = [0]*m
            i = 0
            while t:
                exp = t&mask_exp
                t >>= bits_exp
                a[i] = exp
                i += 1
            r[_tuple(a)] = c
        l[k] = b
    return r


def multinomial_coefficients_iterator(m, n, _tuple=tuple):
    """multinomial coefficient iterator

    Examples:
    >>> from sympy.ntheory.multinomial import multinomial_coefficients_iterator
    >>> it = multinomial_coefficients_iterator(20,3)
    >>> it.next()
    ((3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 1)
    """
    if n == 0:
        d = multinomial_coefficients(m, 0)
        for x in d.items():
            yield x
    elif m < 2*n or n == 1:
        if m == 2:
            yield ((0, n), 1)
            yield ((n, 0), 1)
            a = 1
            for k in range(1, n//2):
                a = (a * (n-k+1))//k
                yield ((k, n-k), a)
                yield ((n-k, k), a)
            if n < 2:
                return
            k = n//2
            a = (a * (n-k+1))//k
            yield ((k, n-k), a)
            if n%2 == 1:
                yield ((n-k, k), a)
            return

        bits_exp, mask_exp, p0, n0, l = _bmpnl_calc(m, n, _tuple)
        yield (_tuple(n0), 1)
        for k in xrange(1, n*(m-1)+1):
            d = _dcalc(m, n, k, p0, l)
            b = []
            for t, c in d.iteritems():
                if not c:
                    continue
                c //= k
                if (t&mask_exp) > 0:
                    b.append((t, c))
                # decode t
                a = [0]*m
                i = 0
                while t:
                    exp = t&mask_exp
                    t >>= bits_exp
                    a[i] = exp
                    i += 1
                yield (_tuple(a), c)
            l[k] = b
    else:
        mc = multinomial_coefficients(n, n)
        mc1 = {}
        for k, v in mc.iteritems():
            mc1[filter(None, k)] = v
        mc = mc1
        comb_it = combinations_with_replacement(range(m), n)
        # one can iterate the monomial tuples by iterating
        # over the combinations with repetition of `n` objects
        # in `range(m)`
        # e.g. for `m=10, n=3 t=(5, 7, 7)` corresponds to
        # the monomial tuple `a = (0, 0, 0, 0, 0, 1, 0, 2, 0, 0) `
        # which in stripped form becomes `b = (1, 2)
        for t in comb_it:
            a = [0]*m
            for i in t:
                a[i] += 1
            a = _tuple(a)
            b = filter(None, a)
            yield (a, mc[b])
