from sympy import I, Add, Mul, Dummy, sympify, Symbol, QQ
from sympy.polys.polytools import resultant, degree, factor_list, Poly
from sympy.polys.polyerrors import NotAlgebraic

def _choose_factor(factors, Y, Y0, dom, prec=200, bound=5):
    if isinstance(factors[0], tuple):
        factors = [f[0] for f in factors]
    if len(factors) == 1:
        return factors[0]

    eps = 1.0/10**prec
    points = {}
    for n in range(bound**len(dom.gens)):
        prec1 = 10
        n_temp = n
        for x in dom.gens:
            points[x] = n_temp % bound
            n_temp = int(n_temp/bound)

        while 1:
            candidates = []
            for f in factors:
                if abs(f.as_expr().subs({Y: Y0}).evalf(prec1, subs=points)) < eps:
                    candidates.append(f)
            if candidates:
                factors = candidates
            if len(factors) == 1:
                return factors[0]
            if prec1 > prec:
                break
            prec1 *= 2
    return None

def _minpoly_op_algebraic_function(ex1, ex2, Y, dom, mp1=None, mp2=None, op=Add):
    from sympy.simplify.simplify import _mexpand

    Z = Dummy(str(Y))
    if mp1 is None:
        mp1 = _minpoly(ex1, Y, dom)
    if mp2 is None:
        mp2 = _minpoly(ex2, Z, dom)
    else:
        mp2 = mp2.subs({Y: Z})

    if op is Add:
        mp1a = mp1.subs({Y: Y - Z})
    elif op is Mul:
        mp1a = _mexpand(Z**degree(mp1, Y)*mp1.subs({Y: Y / Z}))
    else:
        raise NotImplementedError('option not available')

    r = Poly(resultant(mp1a, mp2, gens=[Z, Y]), Y, domain=dom)
    _, factors = r.factor_list()
    res = _choose_factor(factors, Y, dom, op(ex1, ex2))
    return res.as_expr()

def _minpoly_add(Y, dom, *a):
    mp = _minpoly_op_algebraic_function(a[0], a[1], Y, dom, op=Add)
    p = a[0] + a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_function(p, px, Y, dom, mp1=mp, op=Add)
        p = p + px
    return mp

def _minpoly_mul(Y, dom, *a):
    mp = _minpoly_op_algebraic_function(a[0], a[1], Y, dom, op=Mul)
    p = a[0] * a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_function(p, px, Y, dom, mp1=mp, op=Mul)
        p = p * px
    return mp

def _minpoly_pow(ex, pw, Y, dom, mp=None):
    from sympy.simplify.simplify import _mexpand

    pw = sympify(pw)
    if not mp:
        mp = _minpoly(ex, Y, dom)
    if not pw.is_rational:
        raise NotAlgebraic("%s doesn't seem to be an algebraic function" % ex)
    if pw < 0:
        if mp == Y:
            raise ZeroDivisionError('%s is zero' % ex)
        mp = _mexpand(Y**degree(mp, Y)*mp.subs(Y, 1/Y))
        if pw == -1:
            return mp
        pw = -pw
        ex = 1/ex
    Z = Dummy(str(Y))
    mp = mp.subs({Y: Z})
    n, d = pw.as_numer_denom()
    res = Poly(resultant(mp, Y**d - Z**n, gens=[Z]), Y, domain=dom)
    _, factors = factor_list(res, Y)
    res = _choose_factor(factors, Y, dom, ex**pw)
    return res.as_expr()

def _minpoly(ex, Y, dom):
    if ex.is_Rational:
        return ex.q*Y - ex.p
    if ex is I:
        return Y**2 + 1
    if ex in dom.gens:
        return Y - ex

    if ex.is_Add:
        res = _minpoly_add(Y, dom, *ex.args)
    elif ex.is_Mul:
        res = _minpoly_mul(Y, dom, *ex.args)
    elif ex.is_Pow:
        res = _minpoly_pow(ex.base, ex.exp, Y, dom)
    else:
        raise NotAlgebraic("%s doesn't seem to be an algebraic function" % ex)
    return res

def minpoly(ex, Y, dom=None, polys=False):
    from sympy.polys.domains import FractionField
    if dom == None:
        gens = ex.atoms(Symbol)
        dom = FractionField(QQ, *gens)
    if polys:
        return Poly(_minpoly(ex, Y, dom), Y, domain=dom)
    else:
        return _minpoly(ex, Y, dom).collect(Y)
