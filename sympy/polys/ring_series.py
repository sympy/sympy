from sympy.polys.domains import QQ
from sympy.polys.rings import PolyElement, ring
from sympy.polys.monomials import monomial_min, monomial_mul
from mpmath.libmp.libintmath import ifac
from sympy.core.numbers import Rational
from sympy.core.compatibility import as_int, range
from mpmath.libmp.libintmath import giant_steps
import math

def _invert_monoms(p1):
    """
    Compute ``x**n * p1(1/x)`` for ``p1`` univariate polynomial.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import _invert_monoms
    >>> R, x = ring('x', ZZ)
    >>> p = x**2 + 2*x + 3
    >>> _invert_monoms(p)
    3*x**2 + 2*x + 1

    See Also
    ========
    sympy.polys.densebasic.dup_reverse

    """
    terms = list(p1.items())
    terms.sort()
    deg = p1.degree()
    ring = p1.ring
    p = ring.zero
    cv = p1.listcoeffs()
    mv = p1.listmonoms()
    for i in range(len(mv)):
        p[(deg - mv[i][0],)] = cv[i]
    return p

def _giant_steps(target):
    """
    list of precision steps for the Newton's method

    """
    res = giant_steps(2, target)
    if res[0] != 2:
        res = [2] + res
    return res

def rs_trunc(p1, x, prec):
    """
    truncate the series in the ``x`` variable with precision ``prec``,
    that is modulo ``O(x**prec)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_trunc
    >>> R, x = ring('x', QQ)
    >>> p = x**10 + x**5 + x + 1
    >>> rs_trunc(p, x, 12)
    x**10 + x**5 + x + 1
    >>> rs_trunc(p, x, 10)
    x**5 + x + 1
    """

    ring = p1.ring
    p = ring.zero
    i = ring.gens.index(x)
    for exp1 in p1:
        if exp1[i] >= prec:
            continue
        p[exp1] = p1[exp1]
    return p

def rs_mul(p1, p2, x, prec):
    """
    product of series modulo ``O(x**prec)``

    ``x`` is the series variable or its position in the generators.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_mul
    >>> R, x = ring('x', QQ)
    >>> p1 = x**2 + 2*x + 1
    >>> p2 = x + 1
    >>> rs_mul(p1, p2, x, 3)
    3*x**2 + 3*x + 1
    """

    ring = p1.ring
    p = ring.zero
    if ring.__class__ != p2.ring.__class__ or ring != p2.ring:
        raise ValueError('p1 and p2 must have the same ring')
    iv = ring.gens.index(x)
    if not isinstance(p2, PolyElement):
        raise ValueError('p1 and p2 must have the same ring')
    if ring == p2.ring:
        get = p.get
        items2 = list(p2.items())
        items2.sort(key=lambda e: e[0][iv])
        if ring.ngens == 1:
            for exp1, v1 in p1.items():
                for exp2, v2 in items2:
                    exp = exp1[0] + exp2[0]
                    if exp < prec:
                        exp = (exp, )
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break
        else:
            monomial_mul = ring.monomial_mul
            for exp1, v1 in p1.items():
                for exp2, v2 in items2:
                    if exp1[iv] + exp2[iv] < prec:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                    else:
                        break

    p.strip_zero()
    return p

def rs_square(p1, x, prec):
    """
    square modulo ``O(x**prec)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_square
    >>> R, x = ring('x', QQ)
    >>> p = x**2 + 2*x + 1
    >>> rs_square(p, x, 3)
    6*x**2 + 4*x + 1
    """
    ring = p1.ring
    p = ring.zero
    iv = ring.gens.index(x)
    get = p.get
    items = list(p1.items())
    items.sort(key=lambda e: e[0][iv])
    monomial_mul = ring.monomial_mul
    for i in range(len(items)):
        exp1, v1 = items[i]
        for j in range(i):
            exp2, v2 = items[j]
            if exp1[iv] + exp2[iv] < prec:
                exp = monomial_mul(exp1, exp2)
                p[exp] = get(exp, 0) + v1*v2
            else:
                break
    p = p.imul_num(2)
    get = p.get
    for expv, v in p1.items():
        if 2*expv[iv] < prec:
            e2 = monomial_mul(expv, expv)
            p[e2] = get(e2, 0) + v**2
    p.strip_zero()
    return p

def rs_pow(p1, n, x, prec):
    """
    return ``p1**n`` modulo ``O(x**prec)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_pow
    >>> R, x = ring('x', QQ)
    >>> p = x + 1
    >>> rs_pow(p, 4, x, 3)
    6*x**2 + 4*x + 1
    """
    R = p1.ring
    p = R.zero
    if isinstance(n, Rational):
        raise NotImplementedError('to be implemented')

    n = as_int(n)
    if n == 0:
        if p1:
            return R(1)
        else:
            raise ValueError('0**0 is undefined')
    if n < 0:
        p1 = rs_pow(p1, -n, x, prec)
        return rs_series_inversion(p1, x, prec)
    if n == 1:
        return rs_trunc(p1, x, prec)
    if n == 2:
        return rs_square(p1, x, prec)
    if n == 3:
        p2 = rs_square(p1, x, prec)
        return rs_mul(p1, p2, x, prec)
    p = R(1)
    while 1:
        if n&1:
            p = rs_mul(p1, p, x, prec)
            n -= 1
            if not n:
                break
        p1 = rs_square(p1, x, prec)
        n = n // 2
    return p

def _has_constant_term(p, x):
    """
    Check if ``p`` has a constant term in ``x``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import _has_constant_term
    >>> R, x = ring('x', QQ)
    >>> p = x**2 + x + 1
    >>> _has_constant_term(p, x)
    True
    """
    ring = p.ring
    iv = ring.gens.index(x)
    zm = ring.zero_monom
    a = [0]*ring.ngens
    a[iv] = 1
    miv = tuple(a)
    for expv in p:
        if monomial_min(expv, miv) == zm:
            return True
    return False

def _series_inversion1(p, x, prec):
    """
    univariate series inversion ``1/p`` modulo ``O(x**prec)``

    The Newton method is used.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import _series_inversion1
    >>> R, x = ring('x', QQ)
    >>> p = x + 1
    >>> _series_inversion1(p, x, 4)
    -x**3 + x**2 - x + 1

    """
    ring = p.ring
    zm = ring.zero_monom
    if zm not in p:
        raise ValueError('no constant term in series')
    if _has_constant_term(p - p[zm], x):
        raise ValueError('p cannot contain a constant term depending on parameters')
    if p[zm] != ring(1):
        # TODO add check that it is a unit
        p1 = ring(1)/p[zm]
    else:
        p1 = ring(1)
    for precx in _giant_steps(prec):
        tmp = p1.square()
        tmp = rs_mul(tmp, p, x, precx)
        p1 = 2*p1 - tmp
    return p1

def rs_series_inversion(p, x, prec):
    """
    multivariate series inversion ``1/p`` modulo ``O(x**prec)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_series_inversion
    >>> R, x, y = ring('x, y', QQ)
    >>> rs_series_inversion(1 + x*y**2, x, 4)
    -x**3*y**6 + x**2*y**4 - x*y**2 + 1
    >>> rs_series_inversion(1 + x*y**2, y, 4)
    -x*y**2 + 1
    """
    ring = p.ring
    zm = ring.zero_monom
    ii = ring.gens.index(x)
    m = min(p, key=lambda k: k[ii])[ii]
    if m:
        raise NotImplementedError('no constant term in series')
    if zm not in p:
        raise NotImplementedError('no constant term in series')
    if _has_constant_term(p - p[zm], x):
        raise NotImplementedError('p - p[0] must not have a constant term in the series variables')
    return _series_inversion1(p, x, prec)

def rs_series_from_list(p, c, x, prec, concur=1):
    """
    series ``sum c[n]*p**n`` modulo ``O(x**prec)``

    reduce the number of multiplication summing concurrently
    ``ax = [1, p, p**2, .., p**(J - 1)]``
    ``s = sum(c[i]*ax[i] for i in range(r, (r + 1)*J))*p**((K - 1)*J)``
    with ``K >= (n + 1)/J``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_series_from_list, rs_trunc
    >>> R, x = ring('x', QQ)
    >>> p = x**2 + x + 1
    >>> c = [1, 2, 3]
    >>> rs_series_from_list(p, c, x, 4)
    6*x**3 + 11*x**2 + 8*x + 6
    >>> rs_trunc(1 + 2*p + 3*p**2, x, 4)
    6*x**3 + 11*x**2 + 8*x + 6
    >>> pc = R.from_list(list(reversed(c)))
    >>> rs_trunc(pc.compose(x, p), x, 4)
    6*x**3 + 11*x**2 + 8*x + 6

    See Also
    ========
    sympy.polys.ring.compose

    """

    ring = p.ring
    n = len(c)
    if not concur:
        q = ring(1)
        s = c[0]*q
        for i in range(1, n):
            q = rs_mul(q, p, x, prec)
            s += c[i]*q
        return s
    J = int(math.sqrt(n) + 1)
    K, r = divmod(n, J)
    if r:
        K += 1
    ax = [ring(1)]
    b = 1
    q = ring(1)
    if len(p) < 20:
        for i in range(1, J):
            q = rs_mul(q, p, x, prec)
            ax.append(q)
    else:
        for i in range(1, J):
            if i % 2 == 0:
                q = rs_square(ax[i//2], x, prec)
            else:
                q = rs_mul(q, p, x, prec)
            ax.append(q)
    # optimize using rs_square
    pj = rs_mul(ax[-1], p, x, prec)
    b = ring(1)
    s = ring(0)
    for k in range(K - 1):
        r = J*k
        s1 = c[r]
        for j in range(1, J):
            s1 += c[r + j]*ax[j]
        s1 = rs_mul(s1, b, x, prec)
        s += s1
        b = rs_mul(b, pj, x, prec)
        if not b:
            break
    k = K - 1
    r = J*k
    if r < n:
        s1 = c[r]*ring(1)
        for j in range(1, J):
            if r + j >= n:
                break
            s1 += c[r + j]*ax[j]
        s1 = rs_mul(s1, b, x, prec)
        s += s1
    return s

def rs_integrate(self, x):
    """
    integrate ``p`` with respect to ``x``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_integrate
    >>> R, x, y = ring('x, y', QQ)
    >>> p = x + x**2*y**3
    >>> rs_integrate(p, x)
    1/3*x**3*y**3 + 1/2*x**2
    """
    ring = self.ring
    p1 = ring.zero
    n = ring.gens.index(x)
    mn = [0]*ring.ngens
    mn[n] = 1
    mn = tuple(mn)

    for expv in self:
        e = monomial_mul(expv, mn)
        p1[e] = self[expv]/(expv[n] + 1)
    return p1

def fun(p, f, *args):
    """
    Function of a multivariate series computed by substitution

      p: multivariate series
      f: method name or function
      args[:-2] arguments of f, apart from the first one
      args[-2] = iv: names of the series variables
      args[-1] = prec: list of the precisions of the series variables

    The case with f method name is used to compute tan and nth_root
    of a multivariate series:

      fun(p, tan, iv, prec):
      tan series is first computed for a dummy variable _x, ie, tan(_x, iv, prec)
      Then we substitute _x with p to get the desired series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import fun, _tan1
    >>> R, x, y = ring('x, y', QQ)
    >>> p = x + x*y + x**2*y + x**3*y**2
    >>> fun(p, _tan1, x, 4)
    1/3*x**3*y**3 + 2*x**3*y**2 + x**3*y + 1/3*x**3 + x**2*y + x*y + x
    """
    _ring = p.ring
    ring1, _x = ring('_x', _ring)
    h = int(args[-1])
    args1 = args[:-2] + (_x, h)
    zm = _ring.zero_monom
    # separate the constant term of the series
    # compute the univariate series f(_x, .., 'x', sum(nv))
    # or _x.f(..., 'x', sum(nv)
    if zm in p:
        x1 = _x + p[zm]
        p1 = p - p[zm]
    else:
        x1 = _x
        p1 = p
    if isinstance(f, str):
        q = getattr(x1, f)(*args1)
    else:
        q = f(x1, *args1)
    a = sorted(q.items())
    c = [0]*h
    for x in a:
        c[x[0][0]] = x[1]
    p1 = rs_series_from_list(p1, c, args[-2], args[-1])
    return p1

def mul_xin(p, i, n):
    """
    Computes p*x_i**n

    x_i is the ith variable
    """

    n = as_int(n)
    ring = p.ring
    q = ring(0)
    for k, v in p.items():
        k1 = list(k)
        k1[i] += n
        q[tuple(k1)] = v
    return q

def nth_root(p, n, iv, prec):
    """
    Multivariate series of nth root of p
    Computes p**(1/n)

    n: Integer
    iv: List of variable names or variable indices
    prec: List of truncations for these variables

    In the case of one variable it can also be used as:

    iv: Variable name or variable index (0)
    prec: Truncation value for the variable

    p is a series with O(x_1**n_1*..x_m**n_m) in
    variables x_k with index or name iv[k - 1]

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x, y = lgens('x, y', QQ)
    >>> (1 + x + x*y).nth_root(-3, x, 3)
    2/9*x**2*y**2 + 4/9*x**2*y + 2/9*x**2 - 1/3*x*y - 1/3*x + 1
    """
    if n == 0:
        if p == 0:
            raise ValueError('0**0 expression')
        else:
            return p.ring(1)
    if n == 1:
        return rs_trunc(p, iv, prec)
    ring = p.ring
    ii = ring.gens.index(iv)
    m = min(p, key=lambda k: k[ii])[ii]
    try:
        mq, mr = divmod(m, n)
    except TypeError:
        if not gmpy_mode:
            raise ValueError
    if mr:
        raise TaylorEvalError('Not analytic in the series variable')
    p = mul_xin(p, ii, -m)
    prec -= mq
    if _has_constant_term(p-1, iv):
        if ring.zero_monom in p:
            c = p[ring.zero_monom]
    #####
            if isinstance(c, PythonRationalType):
                c1 = Rational(c.p, c.q)
                cn = Pow(c1, S.One/n)
            else:
                cn = Pow(c, S.One/n)
            if cn.is_Rational:
                if not lp.SR:
                    cn = lp.ring(cn.p, cn.q)
                res = cn*(p/c).nth_root(n, iv, prec)
                if mq:
                    res = res.mul_xin(iv, mq)
                return res
        if not lp.SR:
            raise TaylorEvalError('p - 1 must not have a constant term in the series variables')
        else:
            if lp.zero_mon in p:
                c = p[lp.zero_mon]
                if c.is_positive:
                    res = (p/c).nth_root(n, iv, prec)*c**Rational(1, n)
                    if mq:
                        res = res.mul_xin(iv, mq)
                    return res
                else:
                    raise NotImplementedError
    if lp.commuting and lp.ngens == 1:
        res = p._nth_root1(n, iv, prec)
    else:
        res = p.fun('_nth_root1', n, iv, prec)
    if mq:
        res = res.mul_xin(ii, mq)
    return res
#######

def rs_log(p, x, prec):
    """
    logarithm of ``p`` modulo ``O(x**prec)``

    Notes
    =====

    truncation of ``integral dx p**-1*d p/dx`` is used.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_log
    >>> R, x = ring('x', QQ)
    >>> rs_log(1 + x, x, 8)
    1/7*x**7 - 1/6*x**6 + 1/5*x**5 - 1/4*x**4 + 1/3*x**3 - 1/2*x**2 + x
    """
    ring = p.ring
    if p == 1:
        return 0
    if _has_constant_term(p - 1, x):
        raise NotImplementedError('p - 1 must not have a constant term in the series variables')
    dlog = p.diff(x)
    dlog = rs_mul(dlog, _series_inversion1(p, x, prec), x, prec - 1)
    return rs_integrate(dlog, x)

def rs_LambertW(p, iv, prec):
    """
    Calculates the series expansion of principal branch of the Lambert W function. For further
    imformation, refer to the doc of LambertW function.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.lambert(x, 8)
    16807/720*x**7 - 54/5*x**6 + 125/24*x**5 - 8/3*x**4 + 3/2*x**3 - x**2 + x
    """
    ring = p.ring
    p1 = ring(0)
    if _has_constant_term(p, iv):
        raise NotImplementedError('polynomial must not have constant term in the series variables')
    if iv in ring.gens:
        for precx in _giant_steps(prec):
            e = rs_exp(p1, iv, precx)
            p2 = rs_mul(e, p1, iv, precx) - p
            p3 = rs_mul(e, p1 + 1, iv, precx)
            p3 = rs_series_inversion(p3, iv, precx)
            tmp = rs_mul(p2, p3, iv, precx)
            p1 -= tmp
        return p1
    else:
        raise NotImplementedError

def _exp1(p, x, prec):
    """
    helper function for ``rs_exp``
    """
    ring = p.ring
    p1 = ring(1)
    for precx in _giant_steps(prec):
        pt = p - rs_log(p1, x, precx)
        tmp = rs_mul(pt, p1, x, precx)
        p1 += tmp
    return p1

def rs_exp(p, x, prec):
    """
    exponentiation of a series modulo ``O(x**prec)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_exp
    >>> R, x = ring('x', QQ)
    >>> rs_exp(x**2, x, 7)
    1/6*x**6 + 1/2*x**4 + x**2 + 1
    """
    ring = p.ring
    if _has_constant_term(p, x):
        raise NotImplementedError
    if len(p) > 20:
        return _exp1(p, x, prec)
    one = ring(1)
    n = 1
    k = 1
    c = []
    for k in range(prec):
        c.append(one/n)
        k += 1
        n *= k

    r = rs_series_from_list(p, c, x, prec)
    return r

def asin(p, iv, prec):
    """
    arcsine of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.asin(x, 8)
    5/112*x**7 + 3/40*x**5 + 1/6*x**3 + x
    """
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in the series variables')
    ring = p.ring
    if iv in ring.gens:
        # get a good value
        if len(p) > 20:
            dp = p.diff(iv)
            p1 = 1 - rs_square(p, iv, prec - 1)
           ##### p1 = p1.nth_root(-2, iv, prec - 1)
            p1 = rs_mul(dp, p1, iv, prec - 1)
            return rs_integrate(p1, iv)
        one = ring(1)
        c = [0, one, 0]
        for k in range(3, prec, 2):
            c.append((k - 2)**2*c[-2]/(k*(k - 1)))
            c.append(0)
        return rs_series_from_list(p, c, iv, prec)

    else:
        raise NotImplementedError

def _atan_series(p, iv, prec):
    ring = p.ring
    mo = ring(-1)
    c = [-mo]
    p2 = rs_square(p, iv, prec)
    for k in range(1, prec):
        c.append(mo**k/(2*k + 1))
    s = rs_series_from_list(p2, c, iv, prec)
    s = rs_mul(s, p, iv, prec)
    return s

def rs_atan(p, x, prec):
    """
    arctangent of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_atan
    >>> R, x = ring('x', QQ)
    >>> rs_atan(x, x, 8)
    -1/7*x**7 + 1/5*x**5 - 1/3*x**3 + x
    """
    if _has_constant_term(p, x):
        raise NotImplementedError('polynomial must not have constant term in the series variables')
    ring = p.ring
    if x in ring.gens:
        dp = p.diff(x)
        p1 = rs_square(p, x, prec) + ring(1)
        p1 = rs_series_inversion(p1, x, prec - 1)
        p1 = rs_mul(dp, p1, x, prec - 1)
        return rs_integrate(p1, x)
    else:
        return _atan_series(p, x, prec)

def _tan1(p, x, prec):
    """
    Helper function of ``rs_tan``
    """

    ring = p.ring
    p1 = ring(0)
    for precx in _giant_steps(prec):
        tmp = p - rs_atan(p1, x, precx)
        tmp = rs_mul(tmp, 1 + p1.square(), x, precx)
        p1 += tmp
    return p1

def rs_tan(p, x, prec):
    """
    Tangent of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_tan
    >>> R, x = ring('x', QQ)
    >>> rs_tan(x, x, 8)
    17/315*x**7 + 2/15*x**5 + 1/3*x**3 + x
    """
    ring = p.ring
    if _has_constant_term(p, x):
        raise NotImplementedError('p must not have constant part in series variables')
    if ring.ngens == 1:
        return _tan1(p, x, prec)
    return fun(p, rs_tan, x, prec)

def rs_sin(p, x, prec):
    """
    sine of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_sin
    >>> R, x = ring('x', QQ)
    >>> rs_sin(x, x, 6)
    1/120*x**5 - 1/6*x**3 + x
    """
    ring = x.ring
    if not p:
        return ring(0)
    if _has_constant_term(p, x):
        raise NotImplementedError
    # get a good value
    if len(p) > 20 and p.ngens == 1:
        t = rs_tan(p/2, x, prec)
        t2 = rs_square(t, x, prec)
        p1 = rs_series_inversion(1 + t2, x, prec)
        return rs_mul(p1, 2*t, x, prec)
    one = ring(1)
    n = 1
    c = [0]
    for k in range(2, prec + 2, 2):
        c.append(one/n)
        c.append(0)
        n *= -k*(k + 1)
    return rs_series_from_list(p, c, x, prec)

def rs_cos(p, iv, prec):
    """
    cosine of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.cos(x, 6)
    1/24*x**4 - 1/2*x**2 + 1
    """
    ring = p.ring
    if _has_constant_term(p, iv):
        zm = ring.zero_monom
        #if not lp.SR:    Needs to be checked
        #    raise TaylorEvalError
        c = p[zm]
        if not c.is_real:
            raise NotImplementedError
        p1 = p - c
        return sympy.functions.cos(c)*p1.cos(iv, prec) - \
                sympy.functions.sin(c)*p1.sin(iv, prec)
    # get a good value
    if len(p) > 20 and ring.ngens == 1:
        t = rs_tan(p/2, iv, prec)
        t2 = rs_square(t, iv, prec)
        p1 = rs_series_inversion(1+t2, iv, prec)
        return rs_mul_trunc([p1 ,1 - t2, iv, prec)
    one = ring(1)
    n = 1
    c = []
    for k in range(2, prec + 2, 2):
        c.append(one/n)
        c.append(0)
        n *= -k*(k - 1)
    return rs_series_from_list(p, c, iv, prec)

def rs_cos_sin(p, iv, prec):
    """
    Returns the tuple (rs_cos(p, iv, iv), rs_sin(p, iv, iv))
    Is faster than calling rs_cos and rs_sin separately
    """
    t = rs_tan(p/2, iv, prec)
    t2 = rs_square(t, iv, prec)
    p1 = rs_series_inversion(1+t2, iv, prec)
    return (rs_mul(p1, 1 - t2, iv, prec), rs_mul(p1, 2*t, iv, prec))

def check_series_var(p, iv):
    # TODO 
    #Double check whether an error needs to be raised
    ii = p.ring.gens.index(iv)
    m = min(p, key=lambda k: k[ii])[ii]
    return ii, m

def rs_cot(p, iv, prec):
    """
    cotangent of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.cot(x, 6)
    -2/945*x**5 - 1/45*x**3 - 1/3*x + x**-1
    """

    # i, m = p.check_series_var(iv, PoleError, 'cot')
    # Probably not needed
    i, m = check_series_var(p, iv)
    # see _taylor_term1 comment  sin(x**m) = x**pw*sin(x**m)/x**n0
    # prec1 = prec + pw + n0, with pw = n0 = m
    prec1 = prec + 2*m
    c, s = rs_cos_sin(p, iv, prec1)
    s = mul_xin(s, i, -1)
    s = rs_series_inversion(s, iv, prec1)
    res = rs_mul(c, s, iv, prec1)
    res = rs_mul(res, i, -1)
    res = rs_trunc(res, iv, prec)
    return res

def _atanh(p, iv, prec):
    ring = p.ring
    one = ring(1)
    c = [one]
    p2 = rs_square(p, iv, prec)
    for k in range(1, prec):
        c.append(one/(2*k + 1))
    s = rs_series_from_list(p2, c, iv, prec)
    s = rs_mul(s, p, iv, prec)
    return s

def atanh(p, iv, prec):
    """
    Hyperbolic arctangent of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.atanh(x, 8)
    1/7*x**7 + 1/5*x**5 + 1/3*x**3 + x
    """
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in the series variables')
    ring = p.ring
    if iv in ring.gens:
        dp = rs_diff(p, iv)
        p1 = - rs_square(p, iv, prec) + 1
        p1 = rs_series_inversion(p1, iv, prec - 1)
        p1 = rs_mul(dp, p1, iv, prec - 1)
        return rs_integrate(p1, iv)
    else:
        return _atanh(iv, prec)

def rs_sinh(p, iv, prec):
    """
    Hyperbolic sine of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> p = x.sinh(x, 8)
    >>> p
    1/5040*x**7 + 1/120*x**5 + 1/6*x**3 + x
    """
    #self.check_series_var(iv, NotImplementedError, 'sinh')
    t = rs_exp(p, iv, prec)
    t1 = rs_series_inversion(t, iv, prec)
    return (t - t1)/2

def rs_cosh(p, iv, prec):
    """
    Hyperbolic cosine of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.cosh(x, 8)
    1/720*x**6 + 1/24*x**4 + 1/2*x**2 + 1
    """
    #p.check_series_var(iv, NotImplementedError, 'cosh')
    t = rs_exp(p, iv, prec)
    t1 = rs_series_inversion(t, iv, prec)
    return (t + t1)/2

def _tanh(p, iv, prec):
    ring = p.ring
    p1 = ring(0)
    for precx in _giant_steps(prec):
        tmp = p - atanh(p1, iv, precx)
        tmp = rs_mul(tmp, 1 - p1.square(), iv, precx)
        p1 += tmp
    return p1

def tanh(p, iv, prec):
    """
    Hyperbolic tangent of a series

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> lp, x = lgens('x', QQ)
    >>> x.tanh(x, 8)
    -17/315*x**7 + 2/15*x**5 - 1/3*x**3 + x
    """
    ring = p.ring
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in the series variables')
    if and ring.ngens == 1:
        return _tanh(p, iv, prec)
    return fun(p, '_tanh', iv, prec)

def rs_newton(p, x, prec):
    """
    Compute the truncated Newton sum of the polynomial ``p``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_newton
    >>> R, x = ring('x', QQ)
    >>> p = x**2 - 2
    >>> rs_newton(p, x, 5)
    8*x**4 + 4*x**2 + 2
    """
    deg = p.degree()
    p1 = _invert_monoms(p)
    p2 = rs_series_inversion(p1, x, prec)
    p3 = rs_mul(p1.diff(x), p2, x, prec)
    res = deg - p3*x
    return res

def rs_hadamard_exp(p1, inverse=False):
    """
    return ``sum f_i/i!*x**i`` from ``sum f_i*x**i``,
    where ``x`` is the first variable.

    If ``invers=True`` return ``sum f_i*i!*x**i``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_hadamard_exp
    >>> R, x = ring('x', QQ)
    >>> p = 1 + x + x**2 + x**3
    >>> rs_hadamard_exp(p)
    1/6*x**3 + 1/2*x**2 + x + 1
    """
    R = p1.ring
    if R.domain != QQ:
        raise NotImplementedError
    p = R.zero
    if not inverse:
        for exp1, v1 in p1.items():
            p[exp1] = v1/int(ifac(exp1[0]))
    else:
        for exp1, v1 in p1.items():
            p[exp1] = v1*int(ifac(exp1[0]))
    return p

def rs_compose_add(p1, p2):
    """
    compute the composed sum ``prod(p2(x - beta) for beta root of p1)``

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_compose_add
    >>> R, x = ring('x', QQ)
    >>> f = x**2 - 2
    >>> g = x**2 - 3
    >>> rs_compose_add(f, g)
    x**4 - 10*x**2 + 1

    References
    ==========

    A. Bostan, P. Flajolet, B. Salvy and E. Schost
    "Fast Computation with Two Algebraic Numbers",
    (2002) Research Report 4579, Institut
    National de Recherche en Informatique et en Automatique
    """
    R = p1.ring
    x = R.gens[0]
    prec = p1.degree() * p2.degree() + 1
    np1 = rs_newton(p1, x, prec)
    np1e = rs_hadamard_exp(np1)
    np2 = rs_newton(p2, x, prec)
    np2e = rs_hadamard_exp(np2)
    np3e = rs_mul(np1e, np2e, x, prec)
    np3 = rs_hadamard_exp(np3e, True)
    np3a = (np3[(0,)] - np3)/x
    q = rs_integrate(np3a, x)
    q = rs_exp(q, x, prec)
    q = _invert_monoms(q)
    q = q.primitive()[1]
    dp = p1.degree() * p2.degree() - q.degree()
    # `dp` is the multiplicity of the zeroes of the resultant;
    # these zeroes are missed in this computation so they are put here.
    # if p1 and p2 are monic irreducible polynomials,
    # there are zeroes in the resultant
    # if and only if p1 = p2 ; in fact in that case p1 and p2 have a
    # root in common, so gcd(p1, p2) != 1; being p1 and p2 irreducible
    # this means p1 = p2
    if dp:
        q = q*x**dp
    return q
