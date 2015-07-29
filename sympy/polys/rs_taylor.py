
from sympy.polys.rings import ring, xring
from sympy.polys.ring_series import TaylorEvalError, rs_pow, rs_mul, \
rs_nth_root, rs_cos, rs_sin, rs_exp, rs_log, rs_tan, rs_atan, rs_cosh, \
rs_sinh, rs_tanh, rs_atanh, rs_asin, toEX, rs_series_from_list

from sympy.series.order import O
from sympy.core.singleton import S
from sympy.polys.domains import QQ, ZZ, EX
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
from sympy.polys.polyutils import basic_from_dict
from sympy.utilities.iterables import sift

_verbose = 0
def _print_message(p, a, _verbose):
    if _verbose:
        for x in a:
            if x == 's':
                print('series used in', p)

def rational2EX(c):
    if hasattr(c, 'p'):
        c = S(c.p)/S(c.q)
    elif hasattr(c, 'numerator'):
        c = S(c.numerator)/S(c.denominator)
    return c.as_expr()

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

def monomial_as_expr(monom, *gens):
    """Represent a monomial tuple as a SymPy expression in variables
    given in gens.

    Examples
    ========

    >>> from sympy.polys.lpoly import monomial_as_expr
    >>> from sympy import symbols
    >>> x, y = symbols('x, y')
    >>> monomial_as_expr((2, 1), x, y)
    x**2*y
    """
    assert len(monom) == len(gens)
    return Mul(*[Pow(g, m) for (g, m) in zip(gens, monom)])

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
        if not p.exp.is_Integer or p.exp < 0 :
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
    pwn = min([yy[0] for yy in a])
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
        for m1, c1 in num.items():
            c1 = QQ.to_sympy(c1)
            m1 = monomial_as_expr(m1, *gens)
            a.append(c1*m1)
    else:
        for m1, c1 in num.items():
            #c1 = c1.expand()
            m1 = monomial_as_expr(m1, *gens)
            a.append(c1*m1)
    return sympify(sum(a)).as_expr()

_PWMAX = [8, 4]
def _taylor_decompose(p, var, tev, typ, rdeco):
    pw = 0
    s2 = 0
    prec1 = _PWMAX[typ]*rdeco
    typ2 = 0
    while not s2 and not typ2:
        s2, typ2 = _taylor_eval(p, prec1, tev, typ)
        prec1 *= 2
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
    prec = n
    start = x0
    p = sympify(p)
    # TODO deal with some of these cases within taylor
    if var == None or not prec or start != 0 or dir != "+" or \
        prec in [S.Infinity, S.NegativeInfinity]:
        _print_message(p, ['s'], _verbose)
        return p.series(var, start, prec, dir)

    prec = int(prec)
    if var not in p.free_symbols:
        return p

    gens = [var] + pol_pars
    ngens = len(gens)
    lpol_vars = ['X%d' % i for i in range(ngens)]
    # initialize the polynomial classes on QQ and on SR
    lpq = xring(lpol_vars, QQ, lex)[0]
    lps = xring(lpol_vars, EX, lex)[0]
    tev = (TaylorEval(gens, lpq, analytic), TaylorEval(gens, lps, analytic))

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
            p1 = Add(*(pol_part + orders))
            _print_message(p1, ['s'], _verbose)
            res += Add(*(pol_part + orders)).series(var, 0, prec)
            return res

        if pol_part:
            pol_part = Add(*pol_part)
            _print_message(pol_part, ['s'], _verbose)
            res += pol_part.series(var, 0, prec)

        for p1 in non_pol_part:
            c, p1 = p1.as_independent(var, as_Add = False)
            if not p1:
                c = 0
            tres = _taylor_term(p1, var, tev, 0, start, prec, dir, pol_pars, rdeco)
            if not tres:
                _print_message(p1, ['s'], _verbose)
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
            _print_message(p, ['s'], _verbose)
            p1 = p.series(var, start, prec, dir)
            return p1

        p2, typ, pw, n0 = rx
        den, typd = _taylor_eval(p2, prec + n0, tev, typ)
        den = den/den.ring.gens[0]**n0
        c = den[den.ring.zero_monom]
        if c != 1:
            den = den/c
        num = rs_log(den, den.ring.gens[0], prec)
        p2 = _as_expr(num, tev, typd)
        p2 = pw*log(var) + p2
        if c != 1:
            #if isinstance(c, PythonRationalType):
            c = rational2EX(c)
            c = c.as_expr()
            p2 += log(c)
        return p2 + O(var**prec)

    # factor out powers of log(x) and replace them with a polynomial
    # variable `tlog`
    # e.g. p = sin(x)**10*log(x)*log(sin(x) + x) ->
    # sin(x)**10*log((sin(x) + x)/x)*tlog**2
    classes = [q.__class__ for q in p.args]
    if p.is_Mul and log in classes:
        # initialise polynomials in tlog
        lpl = xring('tlog', EX, lex)[0]
        tlog = lpl.gens[0]
        #p12 = S(1)
        p12 = lpl.domain.one
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
                        _print_message(p, ['s'], _verbose)
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
        # FIXME added as_expr because EX has not as_independent
        c, px = px.as_expr().as_independent(var, as_Add = False)
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
                    _print_message(p, ['s'], _verbose)
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
                            pm1 = c*px*logi
                            _print_message(pm1, ['s'], _verbose)
                            ps += pm1.series(var, 0, prec)
                else:
                    prec1 = prec if i == 0 else prec + 1
                    tres =  _taylor_term(px, var, tev, 0, start, prec1, dir, pol_pars, rdeco)
                    if not tres:
                        _print_message(px, ['s'], _verbose)
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
            den = den/den.ring.gens[0]**n0
            prec1 = _iceil(prec + pw)
            if not n.is_Integer:
                c = den[den.ring.zero_monom]
                if c != 1:
                    den = den/c
            else:
                c = 1
            num = rs_pow(den, n, den.ring.gens[0], prec1)
            if c != 1:
                #if isinstance(c, PythonRationalType):
                c = rational2EX(c)
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
                pwn1 = min([yy[0] for yy in a])
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
            lp = den.ring
            ring, lvar = lp.domain, lp.gens[0]
            if n0x:
                den = den/lvar**n0x
            if not nx.is_Integer:
                cx = den[den.ring.zero_monom]
                if cx != 1:
                    den = den/cx
                    #if isinstance(cx, PythonRationalType):
                    cx = rational2EX(cx)
                    c = c*cx**nx
            den = rs_pow(den, nx, lvar, prec1)
            for i in range(1, len(denv)):
                px, typx, n0x, nx = denv[i]
                den1, typd1 = _taylor_eval(px, prec1 + n0x, tev, typx)
                if n0x:
                    den1 = den1/den1.ring.gens[0]**n0x
                if not nx.is_Integer:
                    cx = den1[den1.ring.zero_monom]
                    if cx != 1:
                        den1 = den1/cx
                        #if isinstance(cx, PythonRationalType):
                        cx = rational2EX(cx)
                        c = c*cx**nx
                lvar1 = den1.ring.gens[0]
                den1 = rs_pow(den1, nx, lvar1, prec1)
                if typd1 != typd:
                    if typd1 == 0:
                        den1 = toEX(den1, tev[1].lp)
                    else:
                        den = toEX(den, tev[1].lp)
                        typd = 1
                    lp = tev[1].lp
                den = rs_mul(den, den1, den1.ring.gens[0], prec1)

            num, typn = _taylor_eval(num, prec1, tev, 0)
            if typn == 2:
                return None
            if den != 1:
                if typn < typd:
                    num = toEX(num, tev[1].lp)
                    typn = 1
                elif typn > typd:
                    den = toEX(den, tev[1].lp)
                lvar1 = den.ring.gens[0]
                num = rs_mul(num, den, lvar1, prec1)
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
    def __init__(self, gens, lp, analytic=False):
        self.gens = gens
        self.var = gens[0]
        self.ngens = len(gens)
        self.domain = lp.domain
        self.lp = lp
        self.lvname = 'X0'
        self.lgens = lp.gens
        self.lvar = self.lgens[0]
        self.dgens = dict(zip(self.gens, self.lgens))
        self.analytic = analytic

    def coerce_number(self, elem):
        """conversion of terms independent from var and not in gens
        """
        if self.lp.domain is EX:
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
            s1 = rs_pow(s, p.exp, self.lvar, prec)
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
                    s = rs_mul(s, q, self.lvar, prec1)
                elif q in self.gens:
                    q = self.dgens[q]
                    s = s.mul_trunc(q, self.lvar, prec1)
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
            # TODO explain when this part is used, test it
            if self.var in x.free_symbols:
                if x == self.var:
                    x = self.lvar
                else:
                    x = self(x, prec)
                if pw.is_Integer:
                    pw = int(pw)
                    if pw < 0:
                        raise NotImplementedError
                    x1 = rs_pow(x, pw, self.lvar, prec)
                else:
                    if isinstance(pw, Rational):
                        num = int(pw.p)
                        den = int(pw.q)
                        x1 = rs_pow(x, pw, self.lvar, prec)
                        #x1 = rs_pow(x, num, self.lvar, prec)
                        #x1 = rs_nth_root(x1, den, self.lvar, prec)
                    else:
                        raise NotImplementedError
                return x1
        if head == cos:
            q = self(f.args[0], prec)
            return rs_cos(q, self.lvar, prec)
        if head == sin:
            q = self(f.args[0], prec)
            return rs_sin(q, self.lvar, prec)
        if head == exp:
            q = self(f.args[0], prec)
            return rs_exp(q, self.lvar, prec)
        if head == log:
            q = self(f.args[0], prec)
            return rs_log(q, self.lvar, prec)
        if head == atan:
            q = self(f.args[0], prec)
            return rs_atan(q, self.lvar, prec)
        if head == tan:
            q = self(f.args[0], prec)
            return rs_tan(q, self.lvar, prec)
        if head == cosh:
            q = self(f.args[0], prec)
            return rs_cosh(q, self.lvar, prec)
        if head == sinh:
            q = self(f.args[0], prec)
            return rs_sinh(q, self.lvar, prec)
        if head == tanh:
            q = self(f.args[0], prec)
            return rs_tanh(q, self.lvar, prec)
        if head == atanh:
            q = self(f.args[0], prec)
            return rs_atanh(q, self.lvar, prec)
        if head == asin:
            q = self(f.args[0], prec)
            return rs_asin(q, self.lvar, prec)
        #if head == asinh:
        #    TODO
        #    q = self(f.args[0], prec)
        #    return rs_asinh(q, self.lvar, prec)
        if head == acos:
            #if not self.lp.ring is EX:
            if not self.domain is EX:
                raise TaylorEvalError
            q = self(f.args[0], prec)
            return pi/2 - rs_asin(q, self.lvar, prec)
        #if head == acoth:
            # see issue 564
            # sage gives taylor(acoth(x), x, 0, 9)
            # -1/2*I*pi + 1/9*x^9 + 1/7*x^7 + 1/5*x^5 + 1/3*x^3 + x
            # sage: acoth(0)
            # arccoth(0)
            #if not self.lp.ring is EX:
            #    raise TaylorEvalError
            #q = self(f.args[0], prec)
            #return I*pi/2 + rs_re_acoth(q, self.lvar, prec)
        #if head == LambertW:
        #    q = self(f.args[0], prec)
        #    return rs_lambert(q, self.lvar, prec)
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
        _print_message(g, ['s'], _verbose)
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
        q = rs_series_from_list(q1a, a, q1a.ring.gens[0], prec)
        return q





