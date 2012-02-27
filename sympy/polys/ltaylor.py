"""
This module contains taylor() and different helper functions that it uses.

taylor() attempts taylor expansion using lpoly; if it fails it calls the
series method.
"""

from sympy.polys.lpoly import (LPoly, monomial_as_expr, TaylorEvalError)
from sympy.series.order import O
from sympy.core.singleton import S
from sympy.polys.domains import QQ, ZZ
from sympy.polys.monomialtools import lex
from sympy.functions.elementary.trigonometric import (cos, sin, tan, asin, atan, acos, acot)
from sympy.functions.elementary.exponential import (exp, log, LambertW)
from sympy.functions.elementary.hyperbolic import (sinh, cosh, tanh, atanh, asinh, acosh, acoth)
from sympy.core.numbers import (Number, Rational, Integer)
from sympy.core import pi, expand_multinomial, expand_mul
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy import I, collect
from sympy.core.sympify import sympify
from sympy.polys.domains import PythonRationalType
from sympy.polys.polyutils import basic_from_dict
from sympy.utilities.iterables import sift

# TODO cot(x), coth(x)

def _iceil(prec):
    i = int(prec)
    if i != prec:
        i += 1
    return i

def _is_monomial(p, var):
    """if p is a monomial c*var**n, where n is a real number,
    return the tuple (n, c)
    else return None

    Examples
    ========

    >>> from sympy.polys.ltaylor import _is_monomial
    >>> from sympy import Symbol, log, S
    >>> x = Symbol('x')
    >>> _is_monomial(-x, x)
    (1, -1)
    >>> _is_monomial(log(2)*3*x**5, x)
    (5, 3*log(2))
    >>> _is_monomial(S(2), x)
    (0, 2)
    >>> _is_monomial(x + 1, x) == None
    True
    """
    c, e = p.as_coeff_exponent(var)
    if not any(var in i.free_symbols for i in (c, e)):
        return e, c

def polynomial_degree(p, x):
    """return an upper bound to the degree of p if p is a polynomial in x;
    else return -1

    This function does not attempt any nontrivial simplifications that may
    result in an expression that does not appear to be a polynomial to
    become one; or in a polynomial that appears to have higher degree
    than its expanded form.

    Examples
    ========

    >>> from sympy import Symbol, S
    >>> from sympy.polys.ltaylor import polynomial_degree
    >>> x = Symbol('x')
    >>> polynomial_degree((x**2 + 1)*(1 + x**2) + x*2 + 1, x)
    4
    >>> polynomial_degree(S.One, x)
    0
    >>> polynomial_degree(x/(1 + x), x)
    -1
    """
    if not p.has(x):
        return 0
    if p == x:
        return 1
    if p.is_Pow:
        if p.exp < 0 or not p.exp.is_Integer:
            return -1
        deg_base = polynomial_degree(p.base, x)
        if deg_base < 0:
            return -1
        return deg_base*p.exp
    if p.is_Add:
        m = -1
        for q in p.args:
            n = polynomial_degree(q, x)
            if n < 0:
                return -1
            m = max(m, n)
        return m
    if p.is_Mul:
        a = [polynomial_degree(px, x) for px in p.args]
        if -1 in a:
            return -1
        else:
            s = sum([polynomial_degree(px, x) for px in p.args])
            return s
    else:
        return -1

def poly_truncate(p, x, prec):
    """truncate from a polynomial the monomials of x with power >= prec;
    if the degree of the polynomial is less than prec the output is as
    p.as_expr(x)

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.polys.ltaylor import poly_truncate
    >>> x = Symbol('x')
    >>> p = ((1 + x)**10).as_poly()
    >>> poly_truncate(p, x, 4)
    120*x**3 + 45*x**2 + 10*x + 1
    """
    d = p.rep.to_sympy_dict()
    d1 = {}
    for k, v in d.iteritems():
        if k[0] < prec:
            d1[k] = v
    return basic_from_dict(d1, x)

def _get_var_from_term(p, var):
    """factor the negative power of var from a term
    e.g.  x**-4*sin(x)*cos(x) -> (-4, sin(x)*cos(x))

    Examples
    ========

    >>> from sympy import Symbol, sin, cos
    >>> from sympy.polys.ltaylor import _get_var_from_term
    >>> x = Symbol('x')
    >>> _get_var_from_term(x**-4*sin(x)*cos(x), x)
    (-4, sin(x)*cos(x))
    >>> _get_var_from_term(x**4*sin(x)*cos(x), x)
    (0, x**4*sin(x)*cos(x))
    """
    if p.is_Mul:
        d=sift(Mul.make_args(p), lambda w: \
          w.is_Pow and w.base == var and w.exp.is_integer and w.exp.is_negative)
        rest = Mul(*d[False])
        pw = int(Add(*[q.exp for q in d[True]]))
        return (pw, rest)
    elif p.is_Pow:
        if p.base == var:
            n = p.exp
            if n.is_integer and n.is_negative:
                return (n, 1)
            else:
                return (0, p)
        else:
            return (0, p)
    else:
        return (0, p)

def _factor_var(q, var):
    """factor from a sum of terms the lowest negative power in var,
       without expanding

    return (n, q1), where q = var**n*q1

    Examples
    ========

    >>> from sympy import sin, Symbol
    >>> from sympy.polys.ltaylor import _factor_var
    >>> x = Symbol("x")
    >>> p = x + x**2
    >>> _factor_var(p, x)
    (0, x**2 + x)
    >>> p = 1 + sin(x)/x**2
    >>> _factor_var(p, x)
    (-2, x**2 + sin(x))
    >>> p = 1 + (x + x**2)/x**2
    >>> _factor_var(p, x)
    (-2, 2*x**2 + x)
    """
    if not q.is_Add:
        return (0, q)
    a = []
    num = S.One
    for q1 in q.args:
        a.append(_get_var_from_term(q1, var))
    pwn = min(a)[0]
    if pwn:
        p2 = S.Zero
        for n1, q1 in a:
            if n1 != pwn:
                p2 += var**(n1 - pwn)*q1
            else:
                p2 += q1
        num = num*p2
    else:
        num = num*q
    return (pwn, num)

def _as_expr(num, tev, typn):
    """convert lpoly polynomial num to SymPy expression
    tev list of TaylorEval objects
    typn index ot tev

    Examples
    ========

    >>> from sympy.polys.ltaylor import TaylorEval, _as_expr
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.lpoly import lgens
    >>> from sympy import Symbol, sympify, cos
    >>> x = Symbol('x')
    >>> lpq = LPoly('X', QQ)
    >>> lps = LPoly('X', sympify)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> X = lpq.gens[0]
    >>> _as_expr(X**2 + X + 1, tev, 0)
    x**2 + x + 1
    >>> X = lps.gens[0]
    >>> _as_expr(X**2 + cos(2)*X + 1, tev, 1)
    x**2 + x*cos(2) + 1
    """
    a = []
    gens = tev[typn].gens
    if typn == 0:
        for m1, c1 in num.iteritems():
            c1 = QQ.to_sympy(c1)
            m1 = monomial_as_expr(m1, *gens)
            a.append(c1*m1)
    else:
        for m1, c1 in num.iteritems():
            c1 = c1.expand()
            m1 = monomial_as_expr(m1, *gens)
            a.append(c1*m1)
    return Add(*a)


_PWMAX = [8, 4]
def _taylor_decompose(p, var, tev, typ, rdeco):
    """decompose p(x) in x**pw*pr(x)/x**n0
    where pr is a regular taylor expandible function with tev[typ1],
    and pr**n/x**n0 has constant limit different from 0 for x -> 0;
    return None if the decomposition fails.

    Output: pr, typ1, pw, n0

    Examples
    ========

    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_decompose
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy import Symbol, S, sympify, sin, cos
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ), LPoly('X0', sympify)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> _taylor_decompose(1 + sin(x)/x**2, x, tev, 0, 1)
    (x**2 + sin(x), 0, -1, 1)
    >>> _taylor_decompose(x + sin(x)/x, x, tev, 0, 1)
    (x + sin(x)/x, 0, 0, 0)
    >>> _taylor_decompose(x**2 + sin(x), x, tev, 0, 1)
    (x**2 + sin(x), 0, 1, 1)
    >>> _taylor_decompose(1 - cos(x), x, tev, 0, 1)
    (-cos(x) + 1, 0, 2, 2)
    >>> _taylor_decompose(1 - cos(x**3), x, tev, 0, 1)
    (-cos(x**3) + 1, 0, 6, 6)

    if n0 >= _PWMAX[typ] the decomposition fails
    >>> _taylor_decompose(1 - cos(x**4), x, tev, 0, 1)
    """
    # TODO
    # >>> _taylor_decompose(sqrt(x) + x**(S(3)/2), x, tev, 0)
    # (x + 1, 0, S.Half, 0)
    pw = 0
    s2, typ2 = _taylor_eval(p, _PWMAX[typ]*rdeco, tev, typ)
    n0 = 0
    # if this fails factor 1/x terms in p
    if not s2 or typ2 == 2:
        # p = 1 + sin(x)/x**2
        # pw, pr = -2, sin(x) + x**2
        pw, pr = _factor_var(p, var)
        s2, typ2 = _taylor_eval(pr, _PWMAX[typ]*rdeco, tev, typ)
        # pr = sin(x) + x**2 -> s2 = ... +X0**2 +X0, typ2=0
        if not s2 or typ2 == 2:
            return None
        n0 = min(s2.keys())[0]
        # 1 + sin(x)/x**2 = x**-2*(x**2 + sin(x)) = x**-1*(x**2 + sin(x))/x
        # pw = -1, n0 = 1  p = x**pw*pr/x**n0
        pw += n0
    else:
        # p = x**2 + sin(x)
        # p = x**pw*pr/x**n0, pr=p, n0 = 1, pw = 1
        n0 = min(s2.keys())[0]
        pw = n0
        pr = p
    return pr, typ2, pw, n0

def taylor(p, var=None, x0=0, n=6, dir="+", pol_pars=[], analytic=False, rdeco=1):
    """taylor series expansion of p

    uses the same arguments as series, with the addition
    of three default arguments

      var      series variable
      start    var=start point of expansion
      prec     precision of the series
      dir      direction from which the series is calculated
      pol_pars polynomial parameters
      analytic use the series method for the series expansion of a function
               for which there is no taylor expansion avalable in lpoly;
               only for regular functions.
      rdeco    ratio parameter to determine the precision used
               in _taylor_decompose

    ALGORITHM  separate p in c_1*p_1 + .. + c_n*p_n
    where p_1, .., p_n are product of terms depending on var,
    c1, .., cn are coefficients independent of var

    Try to compute p_i using polynomials; if if fails, compute it with
    the series method.

    The basic case which can be computed with polynomials is
    p_i = n_1*..*nh/(d_1*...*d_k), where n_i and d_j have a regular taylor
    expansion and d_j has constant limit for var -> 0;
    taylor(p_i,var,0,prec) = product taylor(n_i,var,0,prec)*
    product inversion series(d_j,var,0,prec); these can be done fast
    using polynomials; taylor has a very dumb interpreter transforming
    Mul, Pow, and the elementary functions in the corresponding operations
    in polynomials on a ring; it starts with the QQ ring; if the computation
    fails, it passes to the symbolic ring (SR) consisting of the SymPy
    expressions not having `var` and pol_pars in its atoms;
    if also that fails, it passes it to the series method.
    Finally convert the polynomial to a SymPy expression.

    This procedure can be extended to some other cases, e.g. when
    p_i = n_1*..*nh/(d_1*...*d_k) and
    d_i = var**pw_i*dr_i/var**n0_i, where dr_i/var**n0_i
    has regular taylor expansion with constant limit for var -> 0
    (computed by using taylor expansion at order PWMAX; if pw_k >= PWMAX
    fall back to the series method), and similarly but a bit simpler for n_i;
    then using lpoly do the taylor expansion of
    nr_1*..*nr_h/((dr_1/var**n0_1)*...*(dr_k/var**n0_k)) with precision
    depending on the factored power of `var`.

    Examples
    ========

    >>> from sympy.core.symbol import symbols
    >>> from sympy.functions.elementary.trigonometric import (sin, tan, atan)
    >>> from sympy.functions.elementary.exponential import (exp, log)
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.polys.ltaylor import taylor
    >>> from sympy import pi
    >>> from sympy.functions.special.error_functions import erf
    >>> x, y = symbols('x, y')
    >>> taylor(sin(x*tan(x)), x, 0, 10)
    x**2 + x**4/3 - x**6/30 - 71*x**8/630 + O(x**10)

    >>> p = x*exp(x)/(sin(x)*tan(2*x))
    >>> taylor(p, x, 0, 5)
    1/(2*x) + 1/2 - x/3 - x**2/2 - 11*x**3/20 - 67*x**4/180 + O(x**5)

    >>> taylor(sqrt(1 + x*sin(pi*x)), x, 0, 6)
    1 + pi*x**2/2 + x**4*(-pi**3/12 - pi**2/8) + O(x**6)

    >>> taylor(exp(x*log(x)), x, 0, 3)
    1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    In these examples y is first treated internally
    as a SymPy symbol, then as a polynomial parameter;
    the latter version is faster

    >>> taylor(atan(x*y + x**2), x, 0, 5)
    x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)
    >>> taylor(atan(x*y + x**2), x, 0, 5, pol_pars=[y])
    x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)

    example with `analytic=True`; the taylor expansion of erf is done
    using the series method
    >>> taylor(log(sin(2*x)*erf(x)*sqrt(pi)), x, 0, 10, analytic=True)
    2*log(2) + 2*log(x) - x**2 - 2*x**4/45 - 8*x**6/315 - 4*x**8/567 + O(x**10)
    """
    prec = n
    start = x0
    p = sympify(p)
    # TODO deal with some of these cases within taylor
    if var == None or not prec or start != 0 or dir != "+" or \
        prec in [S.Infinity, S.NegativeInfinity]:
        return p.series(var, start, prec, dir)

    prec = int(prec)
    if var not in p.free_symbols:
        return p

    gens = [var] + pol_pars
    ngens = len(gens)
    lpol_vars = ['X%d' % i for i in range(ngens)]
    # initialize the polynomial classes on QQ and on SR
    lpq = LPoly(lpol_vars, QQ, lex)
    lps = LPoly(lpol_vars, sympify, lex)
    tev = (TaylorEval(gens, lpq, analytic), TaylorEval(gens, lps, analytic))

    # taylor(p1 + p2, ...) = taylor(p1, ...) + taylor(p2, ...)
    # in the case in which p2=O(x**prec1), if prec1 < prec
    # series gives a value error
    if p.is_Add:
        orders = []
        # consider first the Order; in case the order is less than
        # prec series raises ValueError, so it is not necessary
        # to do the rest of the computation
        res = S.Zero
        pol_part = []
        non_pol_part = []
        for q in p.args:
            if q.is_Order:
                orders.append(q)
                continue
            m = _is_monomial(q, var)
            if m:
                n0, c0 = m
                if n0 < 0:
                    res += q
                    continue
                if not sympify(n0).is_Integer:
                    res += q
                    continue
            if q.is_polynomial(var):
                pol_part.append(q)
            else:
                non_pol_part.append(q)

        if not non_pol_part:
            res += Add(*(pol_part + orders)).series(var, 0, prec)
            return res

        if pol_part:
            pol_part = Add(*pol_part)
            res += pol_part.series(var, 0, prec)

        for p1 in non_pol_part:
            c, p1 = p1.as_independent(var, as_Add = False)
            if not p1:
                c = 0
            tres = _taylor_term(p1, var, tev, 0, start, prec, dir, pol_pars, rdeco)
            if not tres:
                ts = p1.series(var, 0, prec)
                p2, ord2 = ts.removeO(), ts.getO()
            else:
                p2, ord2 = tres
            if c == 1:
                res += p2
            else:
                res += c*p2
            orders.append(ord2)
        res = res.expand()
        ores = 0
        for o in orders:
            if o:
                ores += o
        return res + ores
    # end of Add case

    # log(p) = log(var**pw*c*p2/(var**n0*c))
    #        = pw*log(var) + log(c) + log(p2/(var**n0*c))
    # where p2/(var**n0*c) has taylor expansion starting with 1
    if p.__class__ == log:
        q = p.args[0]
        rq = _is_monomial(q, var)
        if rq:
            return log(rq[1]) + rq[0]*log(var)
        rx = _taylor_decompose(q, var, tev, 0, rdeco)
        if not rx:
            p1 = p.series(var, start, prec, dir)
            return p1

        p2, typ, pw, n0 = rx
        den, typd = _taylor_eval(p2, prec + n0, tev, typ)
        den = den/den.lp.gens[0]**n0
        c = den[den.lp.zero_mon]
        if c != 1:
            den = den/c
        num = den.log('X0', prec)
        p2 = _as_expr(num, tev, typd)
        p2 = pw*log(var) + p2
        if c != 1:
            if isinstance(c, PythonRationalType):
                c = S(c.numerator)/S(c.denominator)
            p2 += log(c)
        return p2 + O(var**prec)

    # factor out powers of log(x) and replace them with a polynomial
    # variable `tlog`
    # e.g. p = sin(x)**10*log(x)*log(sin(x) + x) ->
    # sin(x)**10*log((sin(x) + x)/x)*tlog**2
    classes = [q.__class__ for q in p.args]
    if p.is_Mul and log in classes:
        # initialise polynomials in tlog
        lpl = LPoly('tlog', sympify, lex)
        tlog = lpl.gens[0]
        p12 = lpl(1)
        for q in p.args:
            if q.__class__ is not log:
                p12 = p12*q
            else:
                qarg = q.args[0]
                rq = _is_monomial(qarg, var)
                if rq:
                    p12 = p12*(log(rq[1]) + rq[0]*tlog)
                else:
                    rx = _taylor_decompose(qarg, var, tev, 0, rdeco)
                    if not rx:
                        p1 = p.series(var, start, prec, dir)
                        return p1

                    p2, typ2, pw2, n2 = rx
                    # if pw2=0 there is no log(x) factor (e.g. for q = log(x + 1));
                    # qarg = sin(x) + x; pw2 = 1; replace qarg by qarg/x**pw2
                    # log(qarg) = log(qarg/x**pw2) + pw2*log(x) ->
                    # log(qarg/x**pw2) + pw2*tlog
                    if pw2 == 0:
                        p12 = p12*q
                    else:
                        p12 = p12*(log(qarg/var**pw2) + pw2*tlog)
        nlog = max(p12.keys())[0]
        p11 = [S.Zero]*(nlog + 1)
        for i in range(nlog + 1):
            if (i, ) in p12:
                p11[i] = p12[(i, )]
    else:
        p11 = [p]
        nlog = 0
    # for p = log(x)*log(sin(x))
    # p12 = (1)*tlog**2 + (log(sin(x)/x))*tlog
    # p11 = [0, log(sin(x)/x), 1]
    ps = 0
    logi = 1
    ords = S.Zero
    for i in range(nlog + 1):
        px = p11[i]
        c, px = px.as_independent(var, as_Add = False)
        if not px:
            c = 0
        if px == 0:
            pass
        if px == var:
            if prec > 1 or i > 0:
                ps = ps + c*var*logi
            else:
                ps = ps + O(var**prec*logi)
            continue
        elif px.is_Pow and \
            px.args[0] == var and var not in px.args[1].free_symbols:
            n = px.args[1]
            if n.is_number:
                if n < prec or (prec == n and i > 0):
                    ps = ps + c*var**n*logi
                else:
                    ps += O(var**prec)
            else:
                if n.is_positive:
                    ps += c*px*logi
                else:
                    return p.series(var, start, prec, dir)
            continue
        else:
            if px == 0:
                pass
            else:
                pxdeg = polynomial_degree(px, var)
                # if px is a polynomial in var it is faster to use the series
                # method or SymPy polynomials
                if pxdeg >= 0:
                    if px.is_Pow:
                        # for polynomials on rationals,
                        # e.g. px = (1 + 3*x + 2*x**2)**100, prec=10
                        # using SymPy polynomials is very fast
                        # px = (pi + 3*x + 2*x**2)**100 is faster
                        # using the series method; in this case
                        # pb.gens = (x, pi)
                        pb = px.base.as_poly()
                        if len(pb.gens) == 1:
                            px1 = pb**px.exp
                            if i == 0:
                                pxp = poly_truncate(px1, var, prec)
                            else:
                                pxp = poly_truncate(px1, var, prec + 1)
                            if pxdeg < prec or pxdeg == prec and i > 0:
                                ps += expand_mul(c*pxp*logi)
                            else:
                                ps += expand_mul(c*pxp*logi) + O(var**prec)
                    if not px.is_Pow or len(pb.gens) != 1:
                        if pxdeg < prec or pxdeg == prec and i > 0:
                            pxp = expand_multinomial(px)
                            ps += expand_mul(c*pxp*logi)
                        else:
                            ps += (c*px*logi).series(var, 0, prec)
                else:
                    prec1 = prec if i == 0 else prec + 1
                    tres =  _taylor_term(px, var, tev, 0, start, prec1, dir, pol_pars, rdeco)
                    if not tres:
                        ts = px.series(var, 0, prec1)
                        p1, ord1 = ts.removeO(), ts.getO()
                    else:
                        p1, ord1 = tres

                    if c*logi == 1:
                        ps += p1
                    else:
                        ps += expand_mul(c*p1*logi, deep=False)
                    if ord1:
                        if i > 0 and tres and not ord1.has(log):
                            ord1 = ord1/var
                        ords += ord1
            logi *= log(var)
    ps = collect(ps, var)
    ps += ords
    return ps

def _taylor_term(p, var, tev, typ=0, start=0, prec=6, dir="+", pol_pars=[], rdeco=1):
    """taylor expansion of a single term p

      tev = (TaylorEval(gens, lpq), TaylorEval(gens, lps))
        where lpq, lps have first variable name 'X0'
      typ index of tev selected

    p is assumed not to be a sum of terms
    if p.is_Mul, p.args does not have terms independent of var
    e.g. p = sin(x)*cos(x)/x
         p = (1 + x)**2
    return (p1, order)

    Examples
    ========

    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_term
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy import Symbol, sympify, cos, sin
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ), LPoly('X0', sympify)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> p = sin(x)*cos(x)/x
    >>> _taylor_term(p, x, tev, 0, 0, 6)
    (2*x**4/15 - 2*x**2/3 + 1, O(x**6))
    """

    assert p.has(var)
    res = _taylor_term1(p, var, tev, typ, prec, rdeco)
    if res == None:
        return None

    c, num, pw, typn = res
    p2 = _as_expr(num, tev, typn)
    if c != 1:
        p2 = (p2*c).expand()
    if pw:
        p2 = (p2*var**pw).expand() # multinomial
    return p2, O(var**prec)

def _taylor_term1(p, var, tev, typ, prec, rdeco):
    """taylor expansion of a single term p

      tev = (TaylorEval(gens, lpq), TaylorEval(gens, lps))
        where lpq, lps have first variable name 'X0'
      typ index of tev selected

    Output: c, p2, pw, typ2

      p2 LPolyElement
      typ2 computed with tev[typ2]

    The taylor expansion is _as_expr(p2, tev, typ2)*c*var**pw

    If the taylor expansion fails it returns None

    Examples
    ========

    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_term1, _as_expr
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy import Symbol, sympify, cos, sin, expand
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ), LPoly('X0', sympify)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> p = (1/x + x**2)**2
    >>> _taylor_term1(p, x, tev, 0, 6, 1)
    (1, X0**6 + 2*X0**3 + 1, -2, 0)
    >>> p = (1/x + cos(x))**2
    >>> c, p2, pw, typ2 = _taylor_term1(p, x, tev, 0, 4, 1); c, p2, pw, typ2
    (1, 1/12*X0**5 - X0**4 - X0**3 + X0**2 + 2*X0 + 1, -2, 0)
    >>> expand(p.series(x,0,4) - _as_expr(p2, tev, typ2)*c*x**pw)
    O(x**4)
    """

    if p.is_Pow:
        # x**a, where a is not real number, e.g. x**I, x**x, x**y
        n = p.args[1]
        p2 = p.args[0]
        if p2 == var:
            # var**number should not be passed to _taylor_term
            assert not(n.is_real and n.is_number)
            # case x**x
            return None
        if n.is_negative or n.is_Rational and not n.is_Integer:
            pw = 0
            rx = _taylor_decompose(p2, var, tev, 0, rdeco)
            if not rx:
                return None
            p2, typ2, pw0, n0 = rx
            pw -= pw0*n
            prec1 = _iceil(prec + pw + n0)
            den, typd = _taylor_eval(p2, prec1, tev, typ2)
            den = den/den.lp.gens[0]**n0
            prec1 = _iceil(prec + pw)
            if not n.is_Integer:
                c = den[den.lp.zero_mon]
                if c != 1:
                    den = den/c
            else:
                c = 1
            num = den.pow_trunc(n, 'X0', prec1)
            if c != 1:
                if isinstance(c, PythonRationalType):
                    c = S(c.p)/S(c.q)
                c = S(c)**n
            return c, num, -pw, typd

        if n.is_positive and n.is_integer:
            pw = 0
            p2 = p.args[0]
            n = int(n)
            pw2, r2 = _factor_var(p2, var)
            prec1 = int(prec - pw2*n)
            p2, typ2 = _taylor_eval(r2**n, prec1, tev, typ)
            if typ2 == 2:
                return None
            pw = pw2*n
            return 1, p2, pw, typ2

    if p.is_Mul:
        # product of terms depending on var
        # p.__class__ == Mul and each of p.args depends on var
        # collect positive and negative powers of the product
        # num = product of the positive powers
        # collect the negative powers in a list
        # collect the leading terms of the negative powers;
        # if the leading terms has power _PWMAX or greater,
        # fall back to the series method
        args = p.args
        num = S.One  # numerator terms of p
        pw = 0   # factor var**-pw of p extracted from p; then
                 # p must be taylor expanded with precision prec + pw
        denv = [] # list of denominator terms to be taylor expanded
        # (pr,typ,n0,n)
        c = S.One
        pwn = 0   #
        for q in args:
            if q.is_Pow:
                # case q = sin(x)**2  :  num = num*sin(x)**2
                # case q = sin(x)**-2 :  denv.append((...))
                qargs = q.args
                n = qargs[1]
                if n.is_Integer and n.is_positive:
                    num = num*q
                elif n.is_negative or n.is_Rational:
                    p2 = qargs[0]
                    rx = _taylor_decompose(p2, var, tev, 0, rdeco)
                    if not rx:
                        return None
                    p2x, typ2x, pw0x, n0x = rx
                    denv.append((p2x, typ2x, n0x, n))
                    pw -= pw0x*n
                else:
                    return None
            elif q.is_Add:
                # case (1 + 1/x**2)/cos(x) -> x**-2*(x**2 + 1)/cos(x)
                a = []
                for q1 in q.args:
                    a.append(_get_var_from_term(q1, var))
                pwn1 = min(a)[0]
                if pwn1:
                    p2 = S.Zero
                    for n1, q1 in a:
                        if n1 != pwn1:
                            p2 += var**(n1 - pwn1)*q1
                        else:
                            p2 += q1
                    num = num*p2
                    pwn += pwn1
                else:
                    num = num*q
            else:
                # case q = cos(x)  :  num = num*q
                num = num*q
        pw -= pwn
        den = 1
        # for each term (pr, typ, n0, n) in denv make the taylor expansion and
        # divide by var**n0, so that inversions can be done within lpoly
        # A nice feature of this approach is that each term can be
        # expanded in its proper ring, e.g. taking an extreme case
        # taylor(1/(exp(log(2)*x**100)*sin(sin(x))), x, 0, 101)
        # exp(log(2)*x**100) is expanded on SR
        if denv:
            prec1 = _iceil(prec + pw)
            px, typx, n0x, nx = denv[0]
            den, typd = _taylor_eval(px, prec1 + n0x, tev, typx)
            lp = den.lp
            ring, lvar = lp.ring, lp.gens[0]
            if n0x:
                den = den/lvar**n0x
            if not nx.is_Integer:
                cx = den[den.lp.zero_mon]
                if cx != 1:
                    den = den/cx
                    if isinstance(cx, PythonRationalType):
                        cx = S(cx.p)/S(cx.q)
                    c = c*cx**nx
            den = den.pow_trunc(nx, 'X0', prec1)
            for i in range(1, len(denv)):
                px, typx, n0x, nx = denv[i]
                den1, typd1 = _taylor_eval(px, prec1 + n0x, tev, typx)
                if n0x:
                    den1 = den1/den1.lp.gens[0]**n0x
                if not nx.is_Integer:
                    cx = den1[den1.lp.zero_mon]
                    if cx != 1:
                        den1 = den1/cx
                        if isinstance(cx, PythonRationalType):
                            cx = S(cx.p)/S(cx.q)
                        c = c*cx**nx
                den1 = den1.pow_trunc(nx, 'X0', prec1)
                if typd1 != typd:
                    if typd1 == 0:
                        den1 = den1.toSR(tev[1].lp)
                    else:
                        den = den.toSR(tev[1].lp)
                        typd = 1
                    lp = tev[1].lp
                den = den.mul_trunc(den1, 'X0', prec1)

            num, typn = _taylor_eval(num, prec1, tev, 0)
            if typn == 2:
                return None
            if den != 1:
                if typn < typd:
                    num = num.toSR(tev[1].lp)
                    typn = 1
                elif typn > typd:
                    den = den.toSR(tev[1].lp)
                num = num.mul_trunc(den, 'X0', prec1)
            return c, num, -pw, typn
        # end of denv clause
        elif pw:
            num, typn = _taylor_eval(num, prec + pw, tev, typ)
            if typn == 2:
                return None
            return 1, num, -pw, typn

    # end of p.is_Mul

    s1, typ1 = _taylor_eval(p, prec, tev, typ)
    if typ1 == 2:
        return None
    return 1, s1, 0, typ1


def _taylor_eval(p, prec, tev, typ):
    """taylor expansion attempted with lpoly

    p SymPy expression
    tev tuple of TaylorEval objects of the form
    tev = (TaylorEval(gens, lpq, prec), TaylorEval(gens, lps, prec))
      lpq LPoly on QQ, lps LPoly on SR
    typ = 0 try first the taylor expansion on QQ
          1 try the evaluation on SR

    output: (res, typ)
      res LPoly.LPolyElement object if successful, else None
      typ = 0 evaluation done in QQ
            1 evaluation done in SR
            2 evaluation failed

    try the taylor expansion first on [QQ, sympify][typ];
    if it fails try it in SR (if typ=0)
    if it succeeds, it returns (res, typ), where res
    is a lpoly polynomial on lpq it typ=0, on lps if typ=1
    if it fails with SR, return (None, 2)

    Examples
    ========

    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_eval
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy import Symbol, sympify, cos
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ), LPoly('X0', sympify)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> _taylor_eval(cos(x + 1), 3, tev, 0)
    ((-cos(1)/2)*X0**2 + (-sin(1))*X0 + (cos(1)), 1)
    """
    typ1 = typ
    for _ in range(typ, 2):
        te = tev[typ1]
        try:
            p1 = te(p, prec)
            return (p1, typ1)
        except TaylorEvalError:
            typ1 += 1
            if typ1 == 2:
                return (None, 2)
        except NotImplementedError:
            return (None, 2)


class TaylorEval:
    """evaluation of SymPy expressions as lpoly polynomials
    """
    def __init__(self, gens, lp, analytic=False):
        """
        prec  precision
        gens polynomial variables in SymPy (the series variable
             and the variables in pol_pars
        var = gens[0] is the series variable in SymPy
        ngens number of polynomial variables
        lp LPoly on which the computation is done
          lp can be on the ring QQ or on SR
          SR is the ring of the expressions in SymPy which
          do not contain the gens variables
        lvname = 'X0'
        lgens polynomial variables in lpoly
        lvar = X0 is the series variable in lpoly
        dgens dictionary to do from SymPy to lpoly variabled

        Examples
        ========

        >>> from sympy.polys.ltaylor import TaylorEval
        >>> from sympy.polys.lpoly import LPoly
        >>> from sympy.polys.domains import QQ
        >>> from sympy import Symbol, cos
        >>> lp = LPoly('X0', QQ)
        >>> x = Symbol('x')
        >>> te = TaylorEval([x], lp)
        >>> te(cos(x), 6)
        1/24*X0**4 - 1/2*X0**2 + 1
        """
        self.gens = gens
        self.var = gens[0]
        self.ngens = len(gens)
        self.ring = lp.ring
        self.lp = lp
        self.lvname = 'X0'
        self.lgens = lp.gens
        self.lvar = self.lgens[0]
        self.dgens = dict(zip(self.gens, self.lgens))
        self.analytic = analytic

    def coerce_number(self, elem):
        """conversion of terms independent from var and not in gens
        """
        if self.lp.SR:
            return elem
        if isinstance(elem, Rational):
            return QQ(elem.p, elem.q)
        else:
            raise TaylorEvalError

    def eval_polynomial(self, p, prec):
        gens = self.gens
        if p in gens:
            return self.dgens[p]
        if not any(p.has(w) for w in gens):
             return self.coerce_number(p)
        if p.is_Mul:
            s = self.lp(1)
            for q in p.args:
                q1 = self.eval_polynomial(q, prec)
                s = s*q1
            return s
        if p.is_Pow:
            s = self.eval_polynomial(p.base, prec)
            s1 = s.pow_trunc(p.exp, 'X0', prec)
            return s1
        if p.is_Add:
            s = self.lp(0)
            for q in p.args:
                s += self.eval_polynomial(q, prec)
            return s

    def __call__(self, f, prec):
        """evaluate the SymPy expression self as a lpoly polynomial
        representing a series of precision `prec`
        """
        prec = int(prec)

        if f in self.gens:
            return self.dgens[f]
        if isinstance(f, Number):
            return self.lp(self.coerce_number(f))

        head = f.__class__
        if head == Add:
            s = self.lp(0)
            for x in f.args:
                if x.is_polynomial(*self.gens):
                    x = self.eval_polynomial(x, prec)
                else:
                    x = self(x, prec)
                s += x
            return s
        if head == Mul:
            s = self.lp(1)
            # sin(x)*cos(x)/x
            pw = 0
            rest = []
            for q in f.args:
                if q.is_Pow:
                    qb, qp = q.args
                    if qb == self.var:
                        if qp.is_integer:
                            pw = -int(qp)
                        else:
                            rest.append(q)
                    else:
                        rest.append(q)
                else:
                    rest.append(q)
            prec1 = prec + pw
            for q in rest:
                if self.var in q.free_symbols:
                    q = self(q, prec1)
                    s = s.mul_trunc(q, self.lvname, prec1)
                elif q in self.gens:
                    q = self.dgens[q]
                    s = s.mul_trunc(q, self.lvname, prec1)
                else:
                    q = self.coerce_number(q)
                    s = s*q
            if pw > 0:
                s = s/self.lvar**pw
            elif pw < 0:
                s = s*self.lvar**-pw

            if min(s)[0] < 0:
                raise NotImplementedError
            return s
        if head == Pow:
            args = f.args
            pw = args[1]
            # f = q*pw(x) = exp(pw(x)*log(q))
            if self.var in pw.free_symbols:
                base = f.args[0]
                # log(x) is not dealt by lpoly
                if base == self.var:
                    raise NotImplementedError
                f = exp(log(f.args[0])*pw)
                x1 = self(f, prec)
                return x1
            x = args[0]
            if self.var in x.free_symbols:
                if x == self.var:
                    x = self.lvar
                else:
                    x = self(x, prec)
                if pw.is_Integer:
                    pw = int(pw)
                    if pw < 0:
                        raise NotImplementedError
                    x1 = x.pow_trunc(pw, self.lvname, prec)
                else:
                    if isinstance(pw, Rational):
                        num = int(pw.p)
                        den = int(pw.q)
                        x1 = x.pow_trunc(num, self.lvname, prec)
                        x1 = x1.nth_root(den, self.lvname, prec)
                    else:
                        raise NotImplementedError
                return x1
        if head == cos:
            q = self(f.args[0], prec)
            return q.cos(self.lvname, prec)
        if head == sin:
            q = self(f.args[0], prec)
            return q.sin(self.lvname, prec)
        if head == exp:
            q = self(f.args[0], prec)
            return q.exp(self.lvname, prec)
        if head == log:
            q = self(f.args[0], prec)
            return q.log(self.lvname, prec)
        if head == atan:
            q = self(f.args[0], prec)
            return q.atan(self.lvname, prec)
        if head == tan:
            q = self(f.args[0], prec)
            return q.tan(self.lvname, prec)
        if head == cosh:
            q = self(f.args[0], prec)
            return q.cosh(self.lvname, prec)
        if head == sinh:
            q = self(f.args[0], prec)
            return q.sinh(self.lvname, prec)
        if head == tanh:
            q = self(f.args[0], prec)
            return q.tanh(self.lvname, prec)
        if head == atanh:
            q = self(f.args[0], prec)
            return q.atanh(self.lvname, prec)
        if head == asin:
            q = self(f.args[0], prec)
            return q.asin(self.lvname, prec)
        if head == asinh:
            q = self(f.args[0], prec)
            return q.asinh(self.lvname, prec)
        if head == acos:
            if not self.lp.SR:
                raise TaylorEvalError
            q = self(f.args[0], prec)
            return pi/2 - q.asin(self.lvname, prec)
        if head == acosh:
            if not self.lp.SR:
                raise TaylorEvalError
            q = self(f.args[0], prec)
            return I*pi/2 - I*q.asin(self.lvname, prec)
        if head == acot:
            if not self.lp.SR:
                raise TaylorEvalError
            q = self(f.args[0], prec)
            return pi/2 + q.acot1(self.lvname, prec)
        if head == acoth:
            # see issue 564
            # sage gives taylor(acoth(x), x, 0, 9)
            # -1/2*I*pi + 1/9*x^9 + 1/7*x^7 + 1/5*x^5 + 1/3*x^3 + x
            # sage: acoth(0)
            # arccoth(0)
            if not self.lp.SR:
                raise TaylorEvalError
            q = self(f.args[0], prec)
            return I*pi/2 + q.re_acoth(self.lvname, prec)
        if head == LambertW:
            q = self(f.args[0], prec)
            return q.lambert(self.lvname, prec)

        if not self.analytic:
            raise NotImplementedError('case in __call__ not considered f=%s' % f)
        # only functions of one variable allowed
        if len(f.args) > 1:
            raise NotImplementedError

        # to compute taylor(f(p(var)), var, 0 ,prec) compute first the series
        # f_series(x) = a[0] + a[1]*x + .. + a[prec - 1]*x**(prec - 1)
        # then replace x with p_series(var); it must be p_series(0) = 0

        # e.g. f = erf(sin(var)); g = erf(var)
        g = head(self.var)
        # do the series expansion of the argument of the function,
        # in this example sin(var), and convert it to lpoly
        # require that the expansion is regular and that it vanishes for
        # var -> 0
        # q1a= (-1/6)*X0**3 + (1)*X0  for prec = 4
        q1a = self(f.args[0], prec)
        if min(q1a)[0] <= 0:
            raise NotImplementedError
        f1 = g.series(self.var, 0, prec)
        f1 = f1.removeO()
        # expansion of g = erf(var); `a` list of coefficients of its series
        # f1= 2*var/sqrt(pi) - 2*var**3/(3*sqrt(pi))
        #   = a[0] + a[1]*var + a[2]*var**2 + a[3]*var**3
        a = [QQ(0)]*prec
        var = self.var
        # case f = var or f = constant
        if not f1.args:
            if f1 == var:
                a[1] += 1
            else:
                a[0] += self.coerce_number(f1)
        for y in f1.args:
            if y == var:
                a[1] += 1
            else:
                if y.is_Pow:
                    if y.base == var and y.exp.is_Integer and y.exp.is_positive:
                        a[int(y.exp)] += 1
                    else:
                        raise NotImplementedError
                elif y.is_Mul:
                    cv = []
                    i = None
                    for yx in y.args:
                        if not yx.has(var):
                            cv.append(yx)
                        else:
                            if yx == var:
                                i = 1
                            elif yx.is_Pow and yx.base == var and yx.exp.is_Integer and yx.exp.is_positive:
                                i = int(yx.exp)
                            else:
                                raise NotImplementedError
                    c = Mul(*cv)
                    c = self.coerce_number(c)
                    if i != None:
                        a[i] += c
                else:
                    raise NotImplementedError
        # compute the taylor expansion of f = erf(sin(var)) using lpoly
        q = q1a.series_from_list(a, 'X0', prec)
        return q

def series_reversion(p, gens):
    """compute the series reversion r = q(w) of p(x) = w

    p(q(w)) = w + O(w**prec) where p is a series with order O(x**prec)

      p    taylor series in series variable x = gens[0] starting with c*x
           where c is an invertible constant
      gens x = gens[0] is the series variable;
           gens[1:-1] are parameters of the series
           w = gens[-1] is the series variable of the reversion series

    Examples
    ========

    >>> from sympy.abc import x, y, z
    >>> from sympy.series.order import O
    >>> from sympy.polys.ltaylor import series_reversion
    >>> gens = [x, y]
    >>> p = x - x**2/2 + x**3/3 - x**4/4 + O(x**5)
    >>> series_reversion(p, [x, y])
    y + y**2/2 + y**3/6 + y**4/24 + O(y**5)
    >>> gens = [x, y, z]
    >>> p = x + x**2 + x**2*y + x**3*y + O(x**4)
    >>> p1 = series_reversion(p, gens); p1
    z - z**2 + 2*z**3 - y*z**2 + 3*y*z**3 + 2*y**2*z**3 + O(z**4)
    >>> p.removeO().subs(x,p1).series(z, 0, 4)
    z + O(z**4)
    """
    order = p.getO()
    prec = order.args[0]
    if not prec.is_Pow or prec.base != gens[0] or not prec.exp.is_Integer or not prec.exp.is_positive:
        raise ValueError('the series must be O(var**n), n > 0')
    prec = prec.exp

    p = p.removeO()
    d = p.as_poly(*gens).as_dict()
    ngens = len(gens)
    lp = LPoly(['X%d' % i for i in range(ngens)], QQ)
    p1 = lp.from_dict(d)
    p2 = p1.series_reversion('X0', prec, 'X%d' % (ngens - 1))
    p3 = p2.as_expr(*gens)
    return p3 + O(gens[-1]**prec)
