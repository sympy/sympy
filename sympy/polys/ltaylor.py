"""ltaylor taylor expansion sing lpoly
"""

from sympy.polys.lpoly import (LPoly, monomial_tobasic, TaylorEvalError)
from sympy.series.order import O
from sympy.core.singleton import S
from sympy.polys.domains import QQ
from sympy.polys.monomialtools import lex
from sympy.functions.elementary.trigonometric import (cos, sin, tan, asin, atan, acos, acot)
from sympy.functions.elementary.exponential import (exp, log, LambertW)
from sympy.functions.elementary.hyperbolic import (sinh, cosh, tanh, atanh, asinh, acosh, acoth)
from sympy.core.numbers import (Number, Rational, Integer)
from sympy.core import pi
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy import I
from sympy.core.sympify import sympify
from sympy.polys.domains import PythonRationalType


# TODO cot(x), coth(x)

def _round(prec):
    if int(prec) != prec:
        prec = int(prec) + 1
    else:
        prec = int(prec)
    return prec

def _is_monomial(p, var):
    """if p is a monomial c*var**n, where n is a real number,
    return the tuple (n, c)
    else return None

    Examples
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
    if p == var:
        return (1, 1)
    # p = x
    if var not in p.atoms():
        return (0, p)
    # p = -x -> (-1, 1)
    # p = log(2)*3*x**5  -> (log(2)*3, 5)
    # p = 2              -> (2, 0)
    if p.__class__ == Mul:
        num_terms_with_var = 0
        n = 0
        c = 1
        for q in p.args:
            if var in q.atoms():
                num_terms_with_var += 1
                if num_terms_with_var > 1:
                    return None
                if q == var:
                    n = 1
                elif q.__class__ == Pow:
                    if q.args[0] != var:
                        return None
                    if q.args[1].is_real and q.args[1].is_number:
                        n = q.args[1]
                    else:
                        return None
                else:
                    return None
            else:
                c = c*q
        return (n, c)
    # x**3 -> (1, 3)
    if p.__class__ == Pow:
        if p.args[0] != var:
            return None
        if p.args[1].is_real:
            return (p.args[1], 1)
        else:
            return None
    return None

def _split_constant_part(p, var):
    """
    Split a product in the part dependent on var and the indendent one;
    else return (p, 1)

    Examples
    >>> from sympy import Symbol, log
    >>> from sympy.polys.ltaylor import _split_constant_part
    >>> x = Symbol('x')
    >>> _split_constant_part(2*log(2)*x**2*log(x+1), x)
    (x**2*log(x + 1), 2*log(2))
    """
    c = 1
    if p.__class__ == Mul:
        p1 = 1
        for q in p.args:
            if var not in q.atoms():
                c *= q
            else:
                p1 *= q
        if c != 1:
            p = p1
    return (p, c)

def _get_var_from_term(p, var):
    """factor the negative power of var from a term
    e.g.  x**-4*sin(x)*cos(x) -> (-4, sin(x)*cos(x))

    Examples
    >>> from sympy import Symbol, sin, cos
    >>> from sympy.polys.ltaylor import _get_var_from_term
    >>> x = Symbol('x')
    >>> _get_var_from_term(x**-4*sin(x)*cos(x), x)
    (-4, sin(x)*cos(x))
    >>> _get_var_from_term(x**4*sin(x)*cos(x), x)
    (0, x**4*sin(x)*cos(x))
    """
    if p.__class__ == Mul:
        pw = 0
        rest = 1
        for q in p.args:
            if var in q.atoms():
                if q.__class__ == Pow:
                    qb, qp = q.args
                    if qb == var:
                        if qp.is_integer and qp.is_negative:
                            pw += int(qp)
                        else:
                            rest *= q
                    else:
                        rest *= q
                else:
                    rest *= q
            else:
                rest *= q
        return (pw, rest)
    elif p.__class__ == Pow:
        if p.args[0] == var:
            n = p.args[1]
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
    return (n, q1), where q = var**n * q1

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
    if q.__class__.is_Pow:
        if q.base == var:
            return (q.exp, S.One)
    if q.__class__.is_Mul:
        cv = []
        pw = 0
        for x in q.args:
            if not x.has(var):
                cv.append(x)
            elif x.is_Pow:
                pw += x.exp
            else:
                break
        else:
            c = Mul(*cv)
            return (pw, c)
    if q.__class__ != Add:
        return (0, q)
    a = []
    num = sympify(1)
    for q1 in q.args:
        a.append(_get_var_from_term(q1, var))
    pwn = min(a)[0]
    if pwn:
        p2 = sympify(0)
        for n1, q1 in a:
            if n1 != pwn:
                p2 += var**(n1-pwn)*q1
            else:
                p2 += q1
        num = num*p2
    else:
        num = num*q
    return (pwn, num)

def tobasic(num, tev, typn):
    """
    convert lpoly polynomial num to Sympy expression
    tev list of TaylorEval objects
    typn index ot tev

    Examples
    >>> from sympy.polys.ltaylor import TaylorEval, tobasic
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy.polys.lpoly import lgens
    >>> from sympy import Symbol, sympify, cos
    >>> x = Symbol('x')
    >>> lpq = LPoly('X', QQ, lex)
    >>> lps = LPoly('X', sympify, lex)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> X = lpq.gens()[0]
    >>> tobasic(X**2 + X + 1, tev, 0)
    x**2 + x + 1
    >>> X = lps.gens()[0]
    >>> tobasic(X**2 + cos(2)*X + 1, tev, 1)
    x**2 + x*cos(2) + 1
    """
    a = []
    gens = tev[typn].gens
    if typn == 0:
        for m1, c1 in num.iteritems():
            c1 = QQ.to_sympy(c1)
            m1 = monomial_tobasic(m1, *gens)
            a.append(c1*m1)
    else:
        for m1, c1 in num.iteritems():
            c1 = c1.expand()
            m1 = monomial_tobasic(m1, *gens)
            a.append(c1*m1)
    return Add(*a)


_PWMAX = [8, 4]
def taylor_decompose(p, var, tev, typ):
    """
    decompose p(x) in pr(x)/x**n0 * x**pw
    where p2r is a regular taylor expandible function with tev[typ1],
    and p2r**n/x**n0 has constant limit different from 0;
    return None if the decomposition fails.

    Output: pr, typ1, pw, n0

    Examples
    >>> from sympy.polys.ltaylor import TaylorEval, taylor_decompose
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy import Symbol, S, sympify, sin
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ, lex), LPoly('X0', sympify, lex)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> taylor_decompose(1 + sin(x)/x**2, x, tev, 0)
    (x**2 + sin(x), 0, -1, 1)
    >>> taylor_decompose(x + sin(x)/x, x, tev, 0)
    (x + sin(x)/x, 0, 0, 0)
    >>> taylor_decompose(x**2 + sin(x), x, tev, 0)
    (x**2 + sin(x), 0, 1, 1)
    """
    # TODO
    # >>> taylor_decompose(sqrt(x) + x**(S(3)/2), x, tev, 0)
    # (x + 1, 0, S.Half, 0)
    pw = 0
    s2, typ2 = taylor_eval(p, _PWMAX[typ], tev, typ)
    n0 = 0
    # if this fails factor 1/x terms in p
    if not s2 or typ2 == 2:
        # p = 1 + sin(x)/x**2
        # pw, pr = -2, sin(x) + x**2
        pw, pr = _factor_var(p, var)
        s2, typ2 = taylor_eval(pr, _PWMAX[typ], tev, typ)
        # pr = sin(x) + x**2->s2 = ... +X0^2 +X0, typ2=0
        if not s2 or typ2 == 2:
            return None
        n0 = min(s2.keys())[0]
        # 1 + sin(x)/x**2 = x**-2*(x**2 + sin(x)) = x**-1*(x**2 + sin(x))/x
        # pw = -1, n0 = 1  p = x**pw * pr/x**n0
        pw += n0
    else:
        # p = x**2 + sin(x)
        # p = x**pw * pr/x**n0, pr=p, n0 = 1, pw = 1
        n0 = min(s2.keys())[0]
        pw = n0
        pr = p
    return pr, typ2, pw, n0


def taylor(p, var=None, x0=0, n=6, dir="+", pol_pars=[], analytic=False):
    """
    taylor series expansion of p
    kept the same arguments as series, with the addition
    of a few default argument for internal use
    var series variable
    start  var=start point of expansion
    prec precision of the series
    dir ...
    pol_pars polynomial parameters


    ALGORITHM separate p in c_1*p_1 + .. + c_n*p_n
    where p_1, .., p_n are product of terms depending on var,
    c1, .., cn are coefficients independent of var

    p_i = num/(d_1*...*d_k)
    den_v = [d_1, .., d_k]
    pwv =   [pw_1, .., pw_k], pw_k smallest power of var in d_i
    (computed by using taylor at order _PWMAX; if pw_k >= _PWMAX
    fall back to the series method)
    pw = sum(pwd)
    then do the taylor expansion of d_i = d_i/var**pw_i at order prec+pwd+pw_i
    finally do the division num/(d_1*...*d_k) within lpoly;
    finally convert to a sympy expression and restore the var**pw factor

    To compute the taylor expansions first compute the series in the
    QQ ring; if this fails compute it
    in the symbolic ring SR consisting
    of the sympy expressions which do not depend on
    var and pol_pars; if also this fails compute it
    using the series function

    EXAMPLES

    >>> from sympy.core.symbol import symbols
    >>> from sympy.functions.elementary.trigonometric import (sin, tan, atan)
    >>> from sympy.functions.elementary.exponential import (exp, log)
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.polys.ltaylor import taylor
    >>> from sympy import pi
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
    as a Sympy symbol, then as a polynomial parameter;
    the latter version is faster

    >>> taylor(atan(x*y + x**2), x, 0, 5)
    x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)
    >>> taylor(atan(x*y + x**2), x, 0, 5, pol_pars=[y])
    x*y + x**2 - x**3*y**3/3 - x**4*y**2 + O(x**5)
    >>> taylor(sum(sin(sin(n*x)) for n in range(1, 4)), x, 0, 10)
    6*x - 12*x**3 + 138*x**5/5 - 6176*x**7/105 + 7293*x**9/70 + O(x**10)
    """
    prec = n
    start = x0
    p = sympify(p)
    # TODO deal with some of these cases within taylor
    if var == None or not prec or start != 0 or dir != "+" or \
        prec in [S.Infinity, S.NegativeInfinity]:
        return p.series(var, start, prec, dir)

    prec = int(prec)
    if var not in p.atoms():
        return p

    gens = [var] + pol_pars
    ngens = len(gens)
    lpol_vars = ['X%d' % i for i in range(ngens)]
    # taylor(p1+p2, ...) = taylor(p1, ...) + taylor(p2, ...)
    # in the case in which p2=O(x**prec1), if prec1 < prec
    # series gives a value error
    if p.__class__ == Add:
        orders = []
        # consider first the Order; in case the order is less than
        # prec series raises ValueError, so it is not necessary
        # to do the rest of the computation
        pol_part = []
        non_pol_part = []
        for q in p.args:
            if q.is_Order:
                orders.append(q)
            q = sympify(q)
            m = _is_monomial(q, var)
            if m:
                pol_part.append(m)
            else:
                non_pol_part.append(q)

        if not non_pol_part:
            if max(pol_part)[0] < prec:
                return p
            else:
                p2 = 0
                for n, c in pol_part:
                    if n < prec:
                        p2 += c*var**n
                p2 += O(var**prec)
                return p2

        res = sympify(0)
        if pol_part:
            for n0, c0 in pol_part:
                res += c0*var**n0

        lpq = LPoly(lpol_vars, QQ, lex)
        lps = LPoly(lpol_vars, sympify, lex)
        tev = (TaylorEval(gens, lpq, analytic), TaylorEval(gens, lps, analytic))

        for p1 in non_pol_part:
            p1, c = _split_constant_part(p1, var)
            p2, ord2 = _taylor_term(p1, var, tev, 0, start, prec, dir, pol_pars)
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

    lpq = LPoly(lpol_vars, QQ, lex)
    lps = LPoly(lpol_vars, sympify, lex)
    tev = (TaylorEval(gens, lpq, analytic), TaylorEval(gens, lps, analytic))

    if p.__class__ == log:
        q = p.args[0]
        rx = taylor_decompose(q, var, tev, 0)
        if not rx:
            p1 = p.series(var, start, prec, dir)
            return p1

        p2, typ, pw, n0 = rx
        den, typd = taylor_eval(p2, prec+n0, tev, typ)
        if typd == 2:
            p1 = p.series(var, start, prec, dir)
            return p1
        den = den/den.lp.gens()[0]**n0
        c = den[den.lp.zero_mon]
        if c != 1:
            den = den/c
        num = den.log('X0', prec)
        p2 = tobasic(num, tev, typd)
        p2 = pw*log(var) + p2
        if c != 1:
              if isinstance(c, PythonRationalType):
                  c = S(c.numerator)/S(c.denominator)
              p2 += log(c)
        return p2 + O(var**prec)  # TOFIX ?

    # factor out powers of log(x) and replace them with a polynomial
    # variable `tlog`
    # e.g. p = sin(x)**10 * log(x)*log(sin(x)+x) ->
    # sin(x)**10 * log((sin(x)+x)/x) * tlog**2
    classes = [q.__class__ for q in p.args]
    if p.__class__ == Mul and log in classes:
        lpl = LPoly('tlog', sympify, lex)
        tlog = lpl.gens()[0]
        p12 = lpl(1)
        for q in p.args:
            if q.__class__ is not log:
                p12 = p12*q
            else:
                qarg = q.args[0]
                rx = taylor_decompose(qarg, var, tev, 0)
                if not rx:
                    p1 = p.series(var, start, prec, dir)
                    return p1

                p2, typ2, pw2, n2 = rx
                # if pw2=0 there is no log(x) factor
                # q = sin(x) + x; pw2 = 1; replace q by q/x**pw2
                # log(q) = log(q/x**pw2) + pw2*log(x) -> log(q/x**pw2) + pw2*tlog
                if pw2 == 0:
                    p12 = p12 * q
                else:
                    p12 = p12*(log(qarg/var**pw2) + pw2*tlog)
        nlog = max(p12.keys())[0]
        #p11 = [p12[(i,)] for i in range(nlog+1)]
        p11 = [0]*(nlog+1)
        for i in range(nlog+1):
            if (i, ) in p12:
                p11[i] = p12[(i, )]
    else:
        p11 = [p]
        nlog = 0
    ps = 0
    logi = 1

    for i in range(nlog+1):
        p = p11[i]
        p, c = _split_constant_part(p, var)
        if p == var:
            if prec > 1:
                ps = ps + c*var*logi
            else:
                ps = ps + O(var**prec*logi)
            continue
        elif p.__class__ == Pow and \
            p.args[0] == var and var not in p.args[1].atoms():
            n = p.args[1]
            if n.is_number:
                if n < prec:
                    ps = ps + c*var**n*logi
                else:
                    ps += c*p*logi + O(var**prec*logi)
            else:
                ps += c*p*logi
            continue
        else:
            p1, ord1 = _taylor_term(p, var, tev, 0, start, prec, dir, pol_pars)
            if ord1:
                if c == 1:
                    ps += p1*logi + ord1*logi
                else:
                    ps += c*p1*logi + ord1*logi
            else:
                if c == 1:
                    ps += p1*logi
                else:
                    ps += c*p1*logi
            logi *= log(var)
    return ps

def _taylor_term(p, var, tev, typ=0, start=0, prec=6, dir="+", pol_pars=[]):
    """taylor expansion of a single term p
    tev = (TaylorEval(gens, lpq), TaylorEval(gens, lps))
    where lpq, lps have first variable name 'X0'
    typ index of tev selected
    p.__class__ != Add
    if p.__class__ == Mul, p.args does not have terms independent of var
    e.g. p = sin(x)*cos(x)/x
         p = (1+x)**2
    return (p1, order)

    Examples
    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_term
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy import Symbol, sympify, cos, sin
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ, lex), LPoly('X0', sympify, lex)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> p = sin(x)*cos(x)/x
    >>> _taylor_term(p, x, tev, 0, 0, 6)
    (2*x**4/15 - 2*x**2/3 + 1, O(x**6))
    """

    res = _taylor_term1(p, var, tev, typ, prec)
    if res == None:
       p1 = p.series(var, start, prec, dir)
       return (p1.removeO(), p1.getO())

    c, num, pw, typn = res
    p2 = tobasic(num, tev, typn)
    if c != 1:
        if isinstance(c, PythonRationalType):
            c = S(c.p)/S(c.q)
        p2 = (p2*c).expand()
    if pw:
        p2 = (p2*var**pw).expand() # multinomial
    return p2, O(var**prec)

def _taylor_term1(p, var, tev, typ, prec):
    """taylor expansion of a single term p
    tev = (TaylorEval(gens, lpq), TaylorEval(gens, lps))
    where lpq, lps have first variable name 'X0'
    typ index of tev selected

    Output: c, p2, pw, typ2
    p2 LPolyElement
    typ2 computed with tev[typ2]
    The taylor expansion is tobasic(p2, tev, typ2) * c * var**pw

    Examples
    >>> from sympy.polys.ltaylor import TaylorEval, _taylor_term1
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy import Symbol, sympify, cos, sin
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ, lex), LPoly('X0', sympify, lex)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> p = (1/x + x**2)**2
    >>> _taylor_term1(p, x, tev, 0, 6)
    (1, X0**6 + 2*X0**3 + 1, -2, 0)
    """

    if p.__class__ == Pow:
        # x**a, where a is not real number, e.g. x**I, x**x, x**y
        n = p.args[1]
        p2 = p.args[0]
        if p2 == var:
            if n.is_real and n.is_number:
                raise ValueError('p should not be passed to _taylor_term')
            else:
                # case x**x
                return None
        if n.is_negative or n.is_positive and n.is_Rational and not n.is_Integer:
            pw = 0
            rx = taylor_decompose(p2, var, tev, 0)
            if not rx:
                return None
            p2, typ2, pw0, n0 = rx
            pw -= pw0*n
            prec1 = _round(prec+pw+n0)
            den, typd = taylor_eval(p2, prec1, tev, typ2)
            if typd == 2:
                return None
            den = den/den.lp.gens()[0]**n0
            # TODO check similar code elsewere
            prec1 = _round(prec + pw)
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
            n = p.args[1]
            n = int(n)
            pw2, r2 = _factor_var(p2, var)
            prec1 = int(prec-pw2*n)
            p2, typ2 = taylor_eval(r2**n, prec1, tev, typ)
            if typ2 == 2:
                return None
            pw = pw2*n
            return 1, p2, pw, typ2

    if p.__class__ == Mul:
        # product of terms depending on var
        # p.__class__ == Mul and each of p.args depends on var
        # collect positive and negative powers of the product
        # num = product of the positive powers
        # collect the negative powers in a list
        # collect the leading terms of the negative powers;
        # if the leading terms has power _PWMAX or greater,
        # fall back to the series method
        # if there is a power which is not integer, or a rational
        # number under certain conditions (TODO), use the series method
        args = p.args
        num = 1  # numerator terms of p
        pw = 0   # factor var**-pw of p extracted from p; then
                 # p must be taylor expanded with precision prec+pw
        denv = [] # list of denominator terms to be taylor expanded
        # (pr,typ,n0,n)
        c = S.One
        pwn = 0   #
        for q in args:
            if q.__class__ == Pow:
                # case q = sin(x)**2  :  num = num*sin(x)**2
                # case q = sin(x)**-2 :  denv.append((...))
                qargs = q.args
                n = qargs[1]
                if n.is_Integer and n.is_positive:
                    num = num*q
                elif n.is_negative or n.is_positive and n.is_Rational and not n.is_Integer:
                    p2 = qargs[0]
                    rx = taylor_decompose(p2, var, tev, 0)
                    if not rx:
                        return None
                    p2x, typ2x, pw0x, n0x = rx
                    denv.append((p2x, typ2x, n0x, n))
                    pw -= pw0x*n
                else:
                    return None
            elif q.__class__ == Add:
                # case (1 + 1/x**2)/cos(x) -> x**-2 * (x**2 + 1)/cos(x)
                a = []
                for q1 in q.args:
                    a.append(_get_var_from_term(q1, var))
                pwn1 = min(a)[0]
                if pwn1:
                    p2 = sympify(0)
                    for n1, q1 in a:
                        if n1 != pwn1:
                            p2 += var**(n1-pwn1)*q1
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
        # for each term (pr,typ,n0,n) in denv make the taylor expansion and
        # divide by var**n0, so that inversions can be done within lpoly
        # A nice feature of this approach is that each term can be
        # expanded in its proper ring, e.g. taking an extreme case
        # taylor(1/(exp(log(2)*x**100)*sin(sin(x))), x, 0, 101)
        # exp(log(2)*x**100) is expanded on SR
        if denv:
            prec1 = _round(prec + pw)
            px, typx, n0x, nx = denv[0]
            den, typd = taylor_eval(px, prec1+n0x, tev, typx)
            if typd == 2:
                return None
            lp = den.lp
            ring, lvar = lp.ring, lp.gens()[0]
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
                den1, typd1 = taylor_eval(px, prec1+n0x, tev, typx)
                if typd1 == 2:
                    return None
                if n0x:
                    den1 = den1/den1.lp.gens()[0]**n0x
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
                    lp = tev[1].lp  # TODO test this line
                den = den.mul_trunc(den1, 'X0', prec1)

            num, typn = taylor_eval(num, prec1, tev, 0)
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
            num, typn = taylor_eval(num, prec+pw, tev, typ)
            if typn == 2:
                return None
            return 1, num, -pw, typn

    # end of p.__class__ == Mul

    s1, typ1 = taylor_eval(p, prec, tev, typ)
    if typ1 == 2:
        return None
    return 1, s1, 0, typ1


def ev_args(te, a, prec):
    if len(a) == 1:
        a = a[0]
        if a == te.var:
            return te.lvar
        if isinstance(a, Number):
            return te.coerce_number(a)
        return te(a, prec)
    else:
        raise NotImplementedError

def taylor_eval(p, prec, tev, typ):
    """taylor expansion attempted with lpoly

    p Sympy expression
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
    >>> from sympy.polys.ltaylor import TaylorEval, taylor_eval
    >>> from sympy.polys.lpoly import LPoly
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy import Symbol, sympify, cos
    >>> x = Symbol('x')
    >>> lpq, lps = LPoly('X0', QQ, lex), LPoly('X0', sympify, lex)
    >>> tev = (TaylorEval([x], lpq), TaylorEval([x], lps))
    >>> taylor_eval(cos(x+1), 3, tev, 0)
    ((-cos(1)/2)*X0**2 + (-sin(1))*X0 + (cos(1)), 1)
    """
    typ1 = typ
    for ring in [QQ, sympify][typ:]:
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
        """
        prec  precision
        gens polynomial variables in Sympy (the series variable
             and the variables in pol_pars
        var = gens[0] is the series variable in Sympy
        ngens number of polynomial variables
        lp LPoly on which the computation is done
          lp can be lpq, with ring QQ, or lps, on SR
          SR is the ring of the expressions in Sympy which
          do not contain the gens variables
        lvname = 'X0'
        lgens polynomial variables in lpoly
        lvar = X0 is the series variable in lpoly
        dgens dictionary to do from Sympy to lpoly variabled
        """
        self.gens = gens
        self.var = gens[0]
        self.ngens = len(gens)
        # try first with ring QQ
        self.ring = lp.ring
        self.lp = lp
        self.lvname = 'X0'
        self.lgens = lp.gens()
        self.lvar = self.lgens[0]
        self.dgens = dict(zip(self.gens, self.lgens))
        self.analytic = analytic

    def coerce_number(self, a):
        ring = self.ring
        if self.lp.SR:
            return a
        if isinstance(a, Rational):
            if ring == QQ:
                return QQ(a.p, a.q)
            return a
        else:
            raise TaylorEvalError


    def __call__(self, f, prec):
        # sometimes prec is Integer, which is not accepted
        # my lpoly
        prec = int(prec)
        if isinstance(f, (int, Number)):
            return self.lp(self.coerce_number(f))
        if f in self.gens:
            return self.dgens[f]

        head = f.__class__
        if head == Add:
            s = self.lp(0)
            for x in f.args:
                if self.var in x.atoms():
                    x = self(x, prec)
                elif x in self.gens:
                    x = self.dgens[x]
                else:
                    x = self.coerce_number(x)
                s += x
            return s
        if head == Mul:
            s = self.lp(1)
            # sin(x)*cos(x)/x
            pw = 0
            rest = []
            for q in f.args:
                if q.__class__ == Pow:
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
                if self.var in q.atoms():
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

            return s
        if head == Pow:
            args = f.args
            pw = args[1]
            # f = q*pw(x) = exp(pw(x)*log(q))
            if self.var in pw.atoms():
                base = f.args[0]
                # log(x) is not dealt by lpoly
                if base == self.var:
                    raise NotImplementedError
                f = exp(log(f.args[0])*pw)
                x1 = self(f, prec)
                return x1
            x = args[0]
            if self.var in x.atoms():
                if x == self.var:
                    x = self.lvar
                else:
                    x = self(x, prec)
                if pw.__class__ == Integer:
                    pw = int(pw)
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
            raise NotImplementedError
        if head == cos:
            q = ev_args(self, f.args, prec)
            return q.cos(self.lvname, prec)
        if head == sin:
            q = ev_args(self, f.args, prec)
            return q.sin(self.lvname, prec)
        if head == exp:
            q = ev_args(self, f.args, prec)
            return q.exp(self.lvname, prec)
        if head == log:
            q = ev_args(self, f.args, prec)
            return q.log(self.lvname, prec)
        if head == atan:
            q = ev_args(self, f.args, prec)
            return q.atan(self.lvname, prec)
        if head == tan:
            q = ev_args(self, f.args, prec)
            return q.tan(self.lvname, prec)
        if head == cosh:
            q = ev_args(self, f.args, prec)
            return q.cosh(self.lvname, prec)
        if head == sinh:
            q = ev_args(self, f.args, prec)
            return q.sinh(self.lvname, prec)
        if head == tanh:
            q = ev_args(self, f.args, prec)
            return q.tanh(self.lvname, prec)
        if head == atanh:
            q = ev_args(self, f.args, prec)
            return q.atanh(self.lvname, prec)
        if head == asin:
            q = ev_args(self, f.args, prec)
            return q.asin(self.lvname, prec)
        if head == asinh:
            q = ev_args(self, f.args, prec)
            return q.asinh(self.lvname, prec)
        if head == acos:
            if not self.lp.SR:
                raise TaylorEvalError
            q = ev_args(self, f.args, prec)
            return pi/2 - q.asin(self.lvname, prec)
        if head == acosh:
            if not self.lp.SR:
                raise TaylorEvalError
            q = ev_args(self, f.args, prec)
            return I*pi/2 - I*q.asin(self.lvname, prec)
        if head == acot:
            if not self.lp.SR:
                raise TaylorEvalError
            q = ev_args(self, f.args, prec)
            return pi/2 + q.acot1(self.lvname, prec)
        if head == acoth:
            # see issue 564
            # sage gives taylor(acoth(x), x, 0, 9)
            # -1/2*I*pi + 1/9*x^9 + 1/7*x^7 + 1/5*x^5 + 1/3*x^3 + x
            # sage: acoth(0)
            # arccoth(0)
            if not self.lp.SR:
                raise TaylorEvalError
            q = ev_args(self, f.args, prec)
            if q == 0:
                return acoth(0)
            return I*pi/2 + q.re_acoth(self.lvname, prec)
        if head == LambertW:
            q = ev_args(self, f.args, prec)
            return q.lambert(self.lvname, prec)

        if not self.analytic:
            raise NotImplementedError('case in __call__ not considered f=%s' % f)
        if len(f.args) > 1:
            raise NotImplementedError

        g = head(self.var)
        q1a = self(f.args[0], prec)
        if min(q1a)[0] <= 0:
            raise NotImplementedError
        f1 = g.series(self.var, 0, prec)
        f1 = f1.removeO()
        a = [QQ(0)]*prec
        var = self.var
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
        q = q1a.series_from_list(a, 'X0', prec)
        return q
