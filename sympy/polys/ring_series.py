from sympy.polys.domains import QQ
from sympy.polys.rings import PolyElement, ring
from sympy.polys.monomials import monomial_min, monomial_mul, monomial_div
from mpmath.libmp.libintmath import ifac
from sympy.core.numbers import Rational
from sympy.core.compatibility import as_int, range
from sympy.core import S
from mpmath.libmp.libintmath import giant_steps
import math

def _invert_monoms(p1):
    """
    Compute ``x**n * p1(1/x)`` for a univariate polynomial ``p1`` in x.

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
    Return a list of precision steps for the Newton's method

    """
    res = giant_steps(2, target)
    if res[0] != 2:
        res = [2] + res
    return res

def rs_trunc(p1, x, prec):
    """
    Truncate the series in the ``x`` variable with precision ``prec``,
    that is, modulo ``O(x**prec)``

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
    Return the product of the given two series, modulo ``O(x**prec)``

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
    Square the series modulo ``O(x**prec)``

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
    Return ``p1**n`` modulo ``O(x**prec)``

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
    Univariate series inversion ``1/p`` modulo ``O(x**prec)``

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
    Multivariate series inversion ``1/p`` modulo ``O(x**prec)``

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
    Return a series ``sum c[n]*p**n`` modulo ``O(x**prec)``

    It reduces the number of multiplications by summing concurrently
    ``ax = [1, p, p**2, .., p**(J - 1)]``
    ``s = sum(c[i]*ax[i] for i in range(r, (r + 1)*J)) * p**((K - 1)*J)``
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

def rs_diff(p, x):
    """
    Computes partial derivative of p with respect to x

      `x`: variable with respect to which p is differentiated,

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_diff
    >>> R, x, y = ring('x, y', QQ)
    >>> p = x + x**2*y**3
    >>> rs_diff(p, x)
    2*x*y**3 + 1
    """
    ring = p.ring
    n = ring.gens.index(x)
    p1 = ring.zero
    mn = [0]*ring.ngens
    mn[n] = 1
    mn = tuple(mn)
    for expv in p:
        if expv[n]:
            e = monomial_div(expv, mn)
            p1[e] = p[expv]*expv[n]
    return p1

def rs_integrate(p, x):
    """
    Integrate ``p`` with respect to ``x``

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
    ring = p.ring
    p1 = ring.zero
    n = ring.gens.index(x)
    mn = [0]*ring.ngens
    mn[n] = 1
    mn = tuple(mn)

    for expv in p:
        e = monomial_mul(expv, mn)
        p1[e] = p[expv]/(expv[n] + 1)
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

      fun(p, tan, iv, prec)
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

    x_i is the ith variable in p
    """
    n = as_int(n)
    ring = p.ring
    q = ring(0)
    for k, v in p.items():
        k1 = list(k)
        k1[i] += n
        q[tuple(k1)] = v
    return q

def rs_log(p, x, prec):
    """
    The Logarithm of ``p`` modulo ``O(x**prec)``

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
    Calculates the series expansion of the principal branch of the Lambert W
    function.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_LambertW
    >>> R, x = ring('x', QQ)
    >>> rs_LambertW(x, x, 8)
    16807/720*x**7 - 54/5*x**6 + 125/24*x**5 - 8/3*x**4 + 3/2*x**3 - x**2 + x

    See Also
    ========

    LambertW
    """
    ring = p.ring
    p1 = ring(0)
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in \
              the series variables')
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
    Helper function for ``rs_exp``
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
    Exponentiation of a series modulo ``O(x**prec)``

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
    The arctangent of a series

    Returns the series expansion of the atan of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_atan
    >>> R, x = ring('x', QQ)
    >>> rs_atan(x, x, 8)
    -1/7*x**7 + 1/5*x**5 - 1/3*x**3 + x

    See Also
    ========

    atan
    """
    if _has_constant_term(p, x):
        raise NotImplementedError('Polynomial must not have constant term in \
              the series variables')
    ring = p.ring

    # Instead of using a closed form formula, we differentiate atan(p) to get
    # `1/(1+p**2) * dp`, whose series expansion is much easier to calculate.
    # Finally we integrate to get back atan
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

    Returns the series expansion of the tan of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_tan
    >>> R, x = ring('x', QQ)
    >>> rs_tan(x, x, 8)
    17/315*x**7 + 2/15*x**5 + 1/3*x**3 + x

   See Also
   ========

   tan
   """
    ring = p.ring
    if _has_constant_term(p, x):
        raise NotImplementedError('Polynomial must not have constant term in \
              series variables')
    if ring.ngens == 1:
        return _tan1(p, x, prec)
    return fun(p, rs_tan, x, prec)

def rs_sin(p, x, prec):
    """
    Sine of a series

    Returns the series expansion of the sin of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_sin
    >>> R, x = ring('x', QQ)
    >>> rs_sin(x, x, 6)
    1/120*x**5 - 1/6*x**3 + x

    See Also
    ========

    sin
    """
    ring = x.ring
    if not p:
        return ring(0)
    # Support for constant term can be extended on the lines of rs_cos
    if _has_constant_term(p, x):
        raise NotImplementedError
    # Series is calcualed in terms of tan as its evaluation is fast.
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
    Cosine of a series

    Returns the series expansion of the cos of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_cos
    >>> R, x = ring('x', QQ)
    >>> rs_cos(x, x, 6)
    1/24*x**4 - 1/2*x**2 + 1

    See Also
    ========

    cos
    """
    ring = p.ring
    if _has_constant_term(p, iv):
        zm = ring.zero_monom
        c = S(p[zm])
        if not c.is_real:
            raise NotImplementedError
        p1 = p - c

    # Makes use of sympy cos, sin fuctions to evaluate the values of the cos/sin
    # of the constant term. Should it be left unevaluated?
        from sympy.functions import cos, sin
        return cos(c)*rs_cos(p1, iv, prec) -  sin(c)*rs_sin(p1, iv, prec)

    # Series is calculated in terms of tan as its evaluation is fast.
    if len(p) > 20 and ring.ngens == 1:
        t = rs_tan(p/2, iv, prec)
        t2 = rs_square(t, iv, prec)
        p1 = rs_series_inversion(1+t2, iv, prec)
        return rs_mul(p1 ,1 - t2, iv, prec)
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
    ii = p.ring.gens.index(iv)
    m = min(p, key=lambda k: k[ii])[ii]
    return ii, m

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

def rs_atanh(p, iv, prec):
    """
    Hyperbolic arctangent of a series

    Returns the series expansion of the atanh of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_atanh
    >>> R, x = ring('x', QQ)
    >>> rs_atanh(x, x, 8)
    1/7*x**7 + 1/5*x**5 + 1/3*x**3 + x

    See Also
    ========

    atanh
    """
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in \
              the series variables')
    ring = p.ring

    # Instead of using a closed form formula, we differentiate atanh(p) to get
    # `1/(1-p**2) * dp`, whose series expansion is much easier to calculate.
    # Finally we integrate to get back atanh
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

    Returns the series expansion of the sinh of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_sinh
    >>> R, x = ring('x', QQ)
    >>> rs_sinh(x, x, 8)
    1/5040*x**7 + 1/120*x**5 + 1/6*x**3 + x

    See Also
    ========

    sinh
    """
    # Check for negative exponent
    t = rs_exp(p, iv, prec)
    t1 = rs_series_inversion(t, iv, prec)
    return (t - t1)/2

def rs_cosh(p, iv, prec):
    """
    Hyperbolic cosine of a series

    Returns the series expansion of the cosh of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_cosh
    >>> R, x = ring('x', QQ)
    >>> rs_cosh(x, x, 8)
    1/720*x**6 + 1/24*x**4 + 1/2*x**2 + 1

    See Also
    ========

    cosh
    """
    # Check for negative exponent
    t = rs_exp(p, iv, prec)
    t1 = rs_series_inversion(t, iv, prec)
    return (t + t1)/2

def _tanh(p, iv, prec):
    ring = p.ring
    p1 = ring(0)
    for precx in _giant_steps(prec):
        tmp = p - rs_atanh(p1, iv, precx)
        tmp = rs_mul(tmp, 1 - p1.square(), iv, precx)
        p1 += tmp
    return p1

def rs_tanh(p, iv, prec):
    """
    Hyperbolic tangent of a series

    Returns the series expansion of the tanh of p, about 0.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.ring_series import rs_tanh
    >>> R, x = ring('x', QQ)
    >>> rs_tanh(x, x, 8)
    -17/315*x**7 + 2/15*x**5 - 1/3*x**3 + x

    See Also
    ========

    tanh
    """
    ring = p.ring
    if _has_constant_term(p, iv):
        raise NotImplementedError('Polynomial must not have constant term in \
              the series variables')
    if ring.ngens == 1:
        return _tanh(p, iv, prec)
    return fun(p, _tanh, iv, prec)

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
    Return ``sum f_i/i!*x**i`` from ``sum f_i*x**i``,
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
