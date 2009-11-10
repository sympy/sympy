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

def multinomial_coefficients(m, n, _tuple=tuple, _zip=zip):
    """Return a dictionary containing pairs ``{(k1,k2,..,km) : C_kn}``
    where ``C_kn`` are multinomial coefficients such that
    ``n=k1+k2+..+km``.

    For example:

    >>> from sympy.ntheory import multinomial_coefficients
    >>> print multinomial_coefficients(2,5)
    {(3, 2): 10, (1, 4): 5, (2, 3): 10, (5, 0): 1, (0, 5): 1, (4, 1): 5}

    The algorithm is based on the following result:

       Consider a polynomial and it's ``m``-th exponent::

         P(x) = sum_{i=0}^m p_i x^k
         P(x)^n = sum_{k=0}^{m n} a(n,k) x^k

       The coefficients ``a(n,k)`` can be computed using the
       J.C.P. Miller Pure Recurrence [see D.E.Knuth, Seminumerical
       Algorithms, The art of Computer Programming v.2, Addison
       Wesley, Reading, 1981;]::

         a(n,k) = 1/(k p_0) sum_{i=1}^m p_i ((n+1)i-k) a(n,k-i),

       where ``a(n,0) = p_0^n``.
    """

    if m==2:
        return binomial_coefficients(n)
    symbols = [(0,)*i + (1,) + (0,)*(m-i-1) for i in range(m)]
    s0 = symbols[0]
    p0 = [_tuple(aa-bb for aa,bb in _zip(s,s0)) for s in symbols]
    r = {_tuple(aa*n for aa in s0):1}
    r_get = r.get
    r_update = r.update
    l = [0] * (n*(m-1)+1)
    l[0] = r.items()
    for k in xrange(1, n*(m-1)+1):
        d = {}
        d_get = d.get
        for i in xrange(1, min(m,k+1)):
            nn = (n+1)*i-k
            if not nn:
                continue
            t = p0[i]
            for t2, c2 in l[k-i]:
                tt = _tuple([aa+bb for aa,bb in _zip(t2,t)])
                cc = nn * c2
                b = d_get(tt)
                if b is None:
                    d[tt] = cc
                else:
                    cc = b + cc
                    if cc:
                        d[tt] = cc
                    else:
                        del d[tt]
        r1 = [(t, c//k) for (t, c) in d.iteritems()]
        l[k] = r1
        r_update(r1)
    return r
