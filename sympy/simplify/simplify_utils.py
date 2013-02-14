from collections import defaultdict
from sympy.core.compatibility import ordered
from sympy.core import (Mul, Add, S)
from sympy.functions import cos, sin
from sympy.utilities.misc import debug

############### polynomial operations ################
# Monomials in a function ``f`` of one variable are represented by a
# dictionary with items ``(pos_i, exp_i)`` representing  ``f(pos_i)**exp_i``

def _fg_div(p1, p2):
    """
    division of two monomial terms

    Parameters
    ==========

    p1, p2 : dictionaries with items ``(pos_i, exp_i)``
    """
    p = {}
    for k, v in p1.iteritems():
        if k in p2:
            r = v - p2[k]
            if r:
                p[k] = r
        else:
            p[k] = v
    for k, v in p2.iteritems():
        if k not in p1:
            p[k] = -v
    return p

def _fg_factor(p1, p2):
    """
    factor the terms with minimum exponent between two monomial

    Parameters
    ==========

    p1, p2 : dictionaries with items ``(pos_i, exp_i)``

    Notes
    =====

    In ``replace_add_gen(expr, f, g, ...)`` a term
    ``f(x_1)**e_1*...f(x_k)**e_k``
    is represented internally by a dictionary
    ``{x_1:e_1, ..., x_k:e_k}``
    Comparing two terms one collects a common factor, and then
    sees if the two terms can be replaced with one term containing
    ``h1`` or ``h2``.
    If all exponents ``e_i`` are non-negative, the factored term
    is the monomial gcd.

    See Also
    ========

    sympy.polys.monomialtools.monomial_gcd

    """
    fact = {}
    for k1, e1 in p1.iteritems():
        if k1 in p2:
            e2 = p2[k1]
            if e1 <= e2:
                fact[k1] = e1
            else:
                fact[k1] = e2
        else:
            if e1 < 0:
                fact[k1] = e1
    for k2, e2 in p2.iteritems():
        if k2 not in p1:
            if e2 < 0:
                fact[k2] = e2

    p1 = _fg_div(p1, fact)
    p2 = _fg_div(p2, fact)
    return p1, p2, fact


################ helper functions ###############

def _splitfg(expr, f, g):
    """
    util function for ``_replace_add_get_a``

    Notes
    =====

    ``expr = c1*f(x11)**e11*f(x12)**e12*...*g(x21)**e21*g(x22)**e22``
    This function gives ``(c1, {x11:e11, x12:e12,...}, {x21:e21, x22:e22,...})
    """
    if not expr.is_Mul:
        args = [expr]
    else:
        args = expr.args
    base_exp = []
    for x in args:
        if x.is_Pow:
            base_exp.append((x.base, x.exp))
        else:
            base_exp.append((x, S.One))
    a0, a1, a2 = [], {}, {}
    for i, (base, exp) in enumerate(base_exp):
        if base.__class__ is f:
            a1[base.args[0]] = exp
        elif base.__class__ is g:
            a2[base.args[0]] = exp
        else:
            a0.append(args[i])
    a0 = Mul(*a0)
    return a0, a1, a2

def _replace_add_get_a(args, f, g):
    """
    util function for ``_get_args_fg``

    Notes
    =====

    ``args = (c1*f(x11)**e11*f(x12)**e12*...*g(x21)**e21*g(x22)**e22*..., ...)``
    This function gives
    ``([[c1, {x11:e11, x12:e12,...}, {x21:e21, x22:e22,...}],...], ``

    Examples
    ========
    >>> from sympy import sin, cos
    >>> from sympy.abc import x, y
    >>> from sympy.simplify.simplify_utils import _replace_add_get_a
    >>> args = (2*sin(x)**2 * cos(x)**2, cos(x))
    >>> _replace_add_get_a(args, sin, cos)
    [[2, {x: 2}, {x: 2}], [1, {}, {x: 1}]]
    """
    a = []
    # separate each addend into a part without products of f and g
    # at top level, in a part with powers of f and a part with powers of g
    for x in args:
        if not x.is_Mul:
            if x.__class__ is f:
                v = [S.One, {x.args[0]:1}, {}]
            elif x.__class__ is g:
                v = [S.One, {}, {x.args[0]:1}]
            elif x.is_Pow:
                x, expx = x.base, x.exp
                if x.__class__ is f:
                    v = [S.One, {x.args[0]: expx}, {}]
                elif x.__class__ is g:
                    v = [S.One, {}, {x.args[0]: expx}]
                else:
                    v = [x, {}, {}]
            else:
                v = [x, {}, {}]
            a.append(v)
            continue
        # with f=sin, g=cos
        # exp(x)*cos(x) -> [[exp(x), [], [(x, 1)]]]
        # sin(x)**2*cos(x)**2 -> [[1, [], [(x, 1)]], [1, [(x, 2)], [(x, 2)]]]
        # [[c1, [{x11:e11, x12:e12,...}, {x21:e21, x22:e22,...}]],...]
        a0, a1, a2 = _splitfg(x, f, g)
        a.append([a0, a1, a2])
    return a

def _get_args_fg(expr, f, g):
    """
    helper function for ``replace_add_gen``
    """
    a  = _replace_add_get_a(expr.args, f, g)
    # Matching can occur only for terms with the same coefficient modulo the sign.
    # Sort ``a`` according to the absolute value of the coefficients,
    # to group terms by them;
    # the rest of the sorting is to get a unique ordering when the keys
    # can be sorted.
    # Reorder args accordingly
    a = zip(a, expr.args)
    a.sort(key=lambda y: (abs(y[0][0]), sorted(list(y[0][1])), sorted(list(y[0][2]))))
    a, args = ([x[0] for x in a], [x[1] for x in a])
    # separate `a` in regions indexed by pos_a
    pos_a = []
    prev = S.Zero
    for i, (a0, a1, a2) in enumerate(a):
        if abs(a0) != abs(prev):
            prev = a0
            pos_a.append(i)
    if pos_a[-1] != len(a):
        pos_a.append(len(a))
    return a, args, pos_a


################ replacement function acting on product of functions ##########
def replace_mul_fpowxgpow(expr, f, g, rexp, h, rexph, mf=None, mg=None):
    """Helper for _match_div_rewrite.

    Parameters
    ==========

    f, g : SymPy functions
    rexp, h, rexph : functions
    mf, mg : functions, if not None

    Notes
    =====

    Replace ``F(f(b_))**c_*G(g(b_))**(rexp(c_))`` with
    ``h(b, *args)**rexph(c)`` if ``F(f(b_))`` and ``G(g(b_))``
    are both positive or if c_ is an integer.

    ``F, G`` are implicitely defined by
    ``mf(F(f(x))) = x``; if ``mf(b, f)`` is ``None`` no match is found.
    Similarly for ``G`` and ``mg``.
    """
    fargs = defaultdict(int)
    gargs = defaultdict(int)
    args = []
    if not mf:
        def mf(x, f):
            if x.__class__ is f:
                return x.args[0]
    if not mg:
        def mg(x, g):
            if x.__class__ is g:
                return x.args[0]
    for x in expr.args:
        if x.is_Pow:
            b, e = x.as_base_exp()
        else:
            b, e = x, S.One

        if e.is_integer or b.is_positive:
            m = mf(b, f)
            if m is not None:
                fargs[m] += e
                continue
            m = mg(b, g)
            if m is not None:
                gargs[m] += e
                continue

        args.append(x)

    common = set(fargs) & set(gargs)
    hit = False
    while common:
        key = common.pop()
        fe = fargs.pop(key)
        ge = gargs.pop(key)
        if fe == rexp(ge):
            args.append(h(key)**rexph(fe))
            hit = True
        else:
            fargs[key] = fe
            gargs[key] = ge
    if not hit:
        return expr
    while fargs:
        key, e = fargs.popitem()
        args.append(f(key)**e)
    while gargs:
        key, e = gargs.popitem()
        args.append(g(key)**e)
    return Mul(*args)

def replace_mul_fpowf2(expr, f, rexp, sgn, h, mn, md):
    """
    replace ``(f(x) + sign*f(y))**e1*(1 - sgn*sign*f(x)*f(y))**e2``
    with ``h(sign, x, y)**e1`` where ``e2 = rexp(e1)``
    """
    nargs = defaultdict(list)
    dargs = defaultdict(list)
    args = []
    for i, x in enumerate(expr.args):
        if x.is_Pow:
            b, e = x.as_base_exp()
        else:
            b, e = x, S.One

        if e.is_integer or b.is_positive:
            m = mn(b, f)
            if m is not None:
                nargs[m[1:]].append((i, m[0], e))
                continue
            m = md(b, f)
            if m is not None:
                coeff, sign, x1, x2 = m
                t = (sgn*sign, x1, x2)
                dargs[t].append((i, coeff, e))
                continue

        args.append(x)

    common = set(nargs) & set(dargs)
    hit = False
    while common:
        key = common.pop()
        a1 = nargs.pop(key)
        a2 = dargs.pop(key)
        e1 = Add(*[x[2] for x in a1])
        e2 = Add(*[x[2] for x in a2])
        if e1 == rexp(e2):
            c1 = Mul(*[x[1]**x[2] for x in a1])
            c2 = Mul(*[x[1]**x[2] for x in a2])
            args.append(c1*c2*h(*key)**e1)
            hit = True
        else:
            for i, c, e in a1:
                args.append(expr.args[i])
            for i, c, e in a2:
                args.append(expr.args[i])
    if not hit:
        return expr
    return Mul(*args)

def __match_f_plus_sign(expr, f, sign):
    """match (sign + f(x)); return ``x`` is found, else None
    """
    if not expr.is_Add:
        return None
    args = expr.args
    if not len(args) == 2:
        return None
    if args[0] == sign:
        b = args[1]
    elif args[1] == sign:
        b = args[0]
    else:
        return None
    if b.__class__ is not f:
        return None
    return b.args[0]

def _match_f_plus_1(expr, f):
    return __match_f_plus_sign(expr, f, 1)
def _match_f_minus_1(expr, f):
    return __match_f_plus_sign(expr, f, -1)


def _match_f1_flus_f2(expr, f):
    """
    match ``c1*f(x1) + c2*f(x2)``; return ``(c1, sign, x1, x2)``
    if ``c2 = sign*c1``, else return ``None``
    """
    if not expr.is_Add:
        return None
    args = expr.args
    if not len(args) == 2:
        return None
    a, b = args
    sa = _splitfg(a, f, None)
    sb = _splitfg(b, f, None)
    sign = sa[0]/sb[0]
    if sign not in (S.One, -S.One):
        return None
    d1 = sa[1]
    d2 = sb[1]
    if len(d1) != 1 or len(d2) != 1:
        return None
    x1, e1 = d1.items()[0]
    x2, e2 = d2.items()[0]
    if e1 != 1 or e2 != 1:
        return None
    a = [x1, x2]
    #a1 = ordered(a)
    #a1 = sorted(a)
    a1 = list(ordered(a))
    if a == a1:
        return (sa[0], sign, x1, x2)
    else:
        return (sb[0], sign, x2, x1)

def _match_ff_plus_1(expr, f):
    """
    match ``(a + sign*f(x1)*f(x2)``, return ``(a, sign, x1, x2)``
    or ``None``
    """
    if not expr.is_Add:
        return None
    args = expr.args
    if not len(args) == 2:
        return None
    a, b = args
    sa = _splitfg(a, f, None)
    sb = _splitfg(b, f, None)
    if not sb[1]:
        sa, sb = sb, sa
    if sa[1]:
        return None
    sign = sa[0]/sb[0]
    if sign not in (S.One, -S.One):
        return None
    d2 = sb[1]
    if len(d2) != 2:
        return None
    (x1, e1), (x2, e2) = d2.items()
    if e1 != 1 or e2 != 1:
        return None
    a = [x1, x2]
    #a1 = sorted(a)
    a1 = list(ordered(a))
    if a == a1:
        return (sa[0], sign, x1, x2)
    else:
        return (sa[0], sign, x2, x1)

def replace_mul_f1(expr, f, g):
    """
    replace a*f(x)**2 with a*g(x), a*f(x)**3 with a*g(x)*f(x)
    """

    if expr.is_Pow:
        args = [expr]
    elif expr.is_Mul:
        args = expr.args
    else:
        return expr
    for i, x in enumerate(args):
        if x.is_Pow:
            if not x.base.__class__ is f:
                continue
            basex, expx = x.base, x.exp
        elif x.__class__ is f:
            basex, expx = x, S.One
        else:
            continue
        basex = basex.args[0]
        if expx == 2:
            r = g(basex)
        elif expx == 3:
            r = g(basex)
        else:
            continue
        if not r.is_Add:
            args1 = args[:]
            args1[i] = r
            if expx == 3:
                args1.append(f(basex))
            return Mul(*args1)
        else:
            a = [x for j, x in enumerate(args) if j != i]
            if expx == 3:
                a.append(f(basex))
            r1 = Mul(*a)
            return Add(*[r1*x for x in r.args])

    return expr

def replace_mul_fapb(expr, f, g):
    """
    replace f(x+y)**z with f(x+y)**(z-1)*g(x, y) for z >=1 and integer

    Do only one replacement
    """
    if not isinstance(expr, Mul):
        args = [expr]
    else:
        args = expr.args
    found = False
    for i, x in enumerate(args):
        if x.is_Pow:
            if not x.base.__class__ is f:
                continue
            basex, expx = x.base, x.exp
        elif x.__class__ is f:
            basex, expx = x, S.One
        else:
            continue
        if not expx.is_integer and expx >= 1:
            continue
        if not isinstance(basex.args[0], Add):
            continue
        if not len(basex.args) == 1:
            continue
        if not len(basex.args[0].args) == 2:
            continue
        x, y = basex.args[0].args
        args1 = list(args)
        r = g(x, y)
        args1[i] = r
        if expx > 1:
            args1.append(f(x+y)**(expx - 1))
        return Mul(*args1)

    return expr

def replace_mul_f2(expr, f, g):
    """
    replace f(2*x)**c/f(x)**c with g(x)**c

    ``f`` is ``sin`` or ``sinh``
    """
    fargs = defaultdict(int)
    args = []
    for x in expr.args:
        if x.is_Pow:
            b, e = x.as_base_exp()
            if e.is_integer or b.is_positive:
                if b.__class__ is f:
                    fargs[b.args[0]] += e
                    continue
        if x.__class__ is f:
            fargs[x.args[0]] += 1
            continue

        args.append(x)
    for x, e1 in fargs.items():
        if 2*x in fargs:
            e2 = fargs[2*x]
            if e1 == -e2:
                args.append(g(x)**e2)
                break
            elif e1 < 0 and e2 > -e1:
                args.append(g(x)**(e2 + e1)*f(2*x)**-e1)
                break
    else:
        return expr
    del fargs[x]
    del fargs[2*x]
    while fargs:
        key, e = fargs.popitem()
        args.append(f(key)**e)
    return Mul(*args)

################ replacement function acting on sum of functions ##########
def replace_add_gen(expr, f, g, h1, h2, h3, sgn, rule, full=True, to_tan=True):
    """
    helper function for replace_add_fgfg, replace_add1, replace_add2

    Parameters
    ==========

    f, g : SymPy functions
    h1, h2, h3 : functions or None
    sgn : sign or None
    rule : replacement rule function
    full : flag for recursive call

    Notes
    =====

    ``expr`` is a sum of terms with ``f`` and ``g`` factors;
    the replacement rule has a pattern which is a sum or difference
    of monomials in ``f`` and ``g`` with coefficient ``1``.
    It is therefore possible to collect the terms of ``expr`` in groups
    which have the same coefficient (the part independent from ``f``
    and ``g``) apart the sign.
    """
    if not expr.is_Add:
        return expr
    if f == g:
        return expr
    a, args, pos_a = _get_args_fg(expr, f, g)
    # convert to dicts, see docstring of _fg_factor
    a_dicts = [x[1:] for x in a]
    found = []
    used = set()
    pos0 = 0
    for pos1 in pos_a:
        for i in range(pos0, pos1 - 1):
            if i in used:
                continue
            for j in range(i + 1, pos1):
                if j in used:
                    continue
                if a[i][0] == a[j][0]:
                    sign = 1
                elif a[i][0] == -a[j][0]:
                    sign = -1
                else:
                    continue
                b1, b2 = a_dicts[i]
                c1, c2 = a_dicts[j]
                # factor the common terms
                b1a, c1a, fact1 = _fg_factor(b1, c1)
                b2a, c2a, fact2 = _fg_factor(b2, c2)

                res = rule(f, g, h1, h2, h3, b1a, b2a, c1a, c2a, fact1,
                           fact2, sgn, sign, to_tan)
                if res is None:
                    continue
                fv, gv, hv = res
                r = Mul(a[i][0], fv, gv, hv)
                found.append(r)
                used.add(i)
                used.add(j)
                break # end j loop
        pos0 = pos1
    if found:
        a_new = [args[i] for i in range(len(args)) if i not in used]
        a_new += found
        res = Add(*a_new)
        if full == True:
            res = replace_add_gen(res, f, g, h1, h2, h3, sgn, rule, full, to_tan)
        return res
    else:
        return expr

def _rule_replace_add1(f, g, h1, h2, h3, b1a, b2a, c1a, c2a, fact1, fact2, sgn, sign, to_tan):
    """
    rule function for replace_add1

    Replace:
    ``1 - sgn*f(x)**2            -> g(x)**2``
    ``1 - g(x)**2           -> sgn*f(x)**2``

    If ``to_tan`` is ``True``,
    use the replacements ``f(x)/g(x) -> h1(x)``, ``g(x)/f(x) -> h2(x)``
    in cases like
    ``f(x)**e1*g(x)**e2*(f(x)**2 - 1) -> -sgn*h2(x)**-e1*g(x)**e3``
    when ``e1 < 0`` and ``e3 = e1 + e2 + 2 >= 0``

    """

    # 1 - sgn*f(x)**2            -> g(x)**2
    # b1a = {} b2a = {}; c1a = {x:2} c2a={}
    # n1 = 0; n2 = 1

    # 1 - g(x)**2           -> sgn*f(x)**2
    # b1a = {} b2a = {}; c1a={} c2a={x:2}
    # n1 = 0; n2 = 1

    # similarly with b, c exchanged
    assert h3 is None
    n1 = len(b1a) + len(b2a)
    n2 = len(c1a) + len(c2a)
    if n1 + n2 != 1:
        return None
    if n1:
        if b1a:
            # 1 - sgn*f(x)**2       -> g(x)**2
            # f(x)**2 - 1 = - sgn*g(x)**2
            # f(x)**e1*g(x)**e2*(f(x)**2 - 1) -> -sgn*f(x)**e1*g(x)**(e2 + 2)
            x, e = b1a.items()[0]
            if e != 2:
                return None
            if sign == -sgn:
                ex1 = fact1[x] if x in fact1 else 0
                ex2 = fact2[x] + 2 if x in fact2 else 2
                if (ex1 >= 0 and ex2 >= 0) or (ex1 <= 0 and ex2 <= 0):
                    op = 1
                    hv = g(x)**2
                elif to_tan and ex1 < 0 and ex2 > 0:
                    op = 2
                    hv = -sgn*h2(x)**-ex1 * g(x)**(ex1 + ex2)
                else:
                    #assert ex1 > 0 and ex2 < 0
                    if not to_tan:
                        return None
                    op = 2
                    hv = -sgn*h1(x)**-ex2 * f(x)**(ex1 + ex2)
            else:
                return None
        else:
            x, e = b2a.items()[0]
            if e != 2:
                return None
            # g(x)**2 - 1       -> -sgn*f(x)**2
            # f(x)**e1*g**e2*(g(x)**2 - 1) ->
            # -sgn*f(x)**(e1 + 2)*g**e2
            if sign == -1:
                ex1 = fact1[x] + 2 if x in fact1 else 2
                ex2 = fact2[x] if x in fact2 else 0
                if (ex1 >= 0 and ex2 >= 0) or (ex1 <= 0 and ex2 <= 0):
                    op = 1
                    hv = -sgn*f(x)**2
                elif to_tan and ex1 < 0 and ex2 > 0:
                    # f(x)**ex1*g**ex2 = h2(x)**-ex1 * g**(ex1 + ex2)
                    op = 2
                    hv = -sgn*h2(x)**-ex1 * g(x)**(ex1 + ex2)
                else:
                    # ex1 > 0 and ex2 < 0:
                    if not to_tan:
                        return None
                    op = 2
                    hv = -sgn*h1(x)**-ex2 * f(x)**(ex1 + ex2)
            else:
                return None
    else:
        assert n2 == 1
        if c1a:
            x, e = c1a.items()[0]
            if e != 2:
                return None
            # 1 - sgn*f(x)**2       -> g(x)**2
            # f(x)**e1*g**e2*(1 - sgn*f(x)**2) ->
            # f(x)**e1*g**(e2 + 2)
            if sign == -sgn:
                ex1 = fact1[x] if x in fact1 else 0
                ex2 = fact2[x] + 2 if x in fact2 else 2
                if (ex1 >= 0 and ex2 >= 0) or (ex1 <= 0 and ex2 <= 0):
                    op = 1
                    hv = g(x)**2
                elif to_tan and ex1 < 0 and ex2 > 0:
                    op = 2
                    hv = h2(x)**-ex1 * g(x)**(ex1 + ex2)
                else:
                    #assert ex1 > 0 and ex2 < 0
                    if not to_tan:
                        return None
                    op = 2
                    hv = h1(x)**-ex2 * f(x)**(ex1 + ex2)
            else:
                return None
        else:
            x, e = c2a.items()[0]
            if e != 2:
                return None
            # 1 - g(x)**2           -> sgn*f(x)**2
            # f(x)**e1*g**e2*(1 - g(x)**2) -> sgn*f(x)**(e1 + 2)*g**e2
            if sign == -1:
                ex1 = fact1[x] + 2if x in fact1 else 2
                ex2 = fact2[x] if x in fact2 else 0

                if (ex1 >= 0 and ex2 >= 0) or (ex1 <= 0 and ex2 <= 0):
                    op = 1
                    hv = sgn*f(x)**2
                elif to_tan and ex1 < 0 and ex2 > 0:
                    op = 2
                    hv = sgn*h2(x)**-ex1 * g(x)**(ex1 + ex2)
                else:
                    #assert ex1 > 0 and ex2 < 0
                    if not to_tan:
                        return None
                    op = 2
                    hv = sgn*h1(x)**-ex2 * f(x)**(ex1 + ex2)

            else:
                return None
    if op == 1:
        fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
        gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
    else:
        fv = Mul(*[f(xx)**yy for xx, yy in fact1.items() if xx != x])
        gv = Mul(*[g(xx)**yy for xx, yy in fact2.items() if xx != x])
    return fv, gv, hv

def _rule_replace_add2(f, g, h1, h2, h3, b1a, b2a, c1a, c2a, fact1, fact2, sgn, sign, to_tan):
    """
    helper function for replace_add2
    """
    # f(x)**2 + sgn*g(x)**2 -> 1
    # b1a={x:2} b2a={}; c1a={}; c2a={x:2}
    # n1 = 1; n2 = 1

    # similarly with b, c exchanged
    assert h3 is None
    n1 = len(b1a) + len(b2a)
    n2 = len(c1a) + len(c2a)
    if n1 != 1 or n2 != 1:
        return None
    if (b1a != c2a or b2a != c1a):
        return None
    if b1a:
        # f(x)**2 + sgn*g(x)**2 -> 1
        x, e = b1a.items()[0]
        if e != 2:
            return None
        if sign == sgn:
            fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
            gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
            hv = S.One
        else:
            return None
    elif b2a:
        # g(x)**2 + sgn*f(x)**2 -> sgn
        x, e = b2a.items()[0]
        if e != 2:
            return None
        if sign == sgn:
            fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
            gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
            hv = sign
        else:
            return None
    return fv, gv, hv


def _rule_replace_add_fgfg(f, g, h1, h2, h3, b1a, b2a, c1a, c2a, fact1, fact2, sgs, sign, to_tan):
    """
    helper function for replace_add_fgfg
    """
    # f(x)*g(y) + sign*f(y)*g(x) -> h2(x, y, sign)
    # there must be two common terms, with exponent 1
    # f(x)*f(y) + sign*g(x)*g(y) -> h1(x, y, sign)
    # terms with two ``f`` or two ``g``
    assert sgs is None
    n1 = len(b1a) + len(b2a)
    if n1 not in (1, 2):
        return None
    n2 = len(c1a) + len(c2a)
    if n2 not in (1, 2):
        return None
    if n1 == 1 and not h3:
        return None

    if b1a != c2a or b2a != c1a:
        return None
    # The replacement is done only if ``x`` and ``y``
    # are different, so that the exponents ``exx``
    # are equal to 1.
    if b1a and any(exx != 1 for exx in b1a.values()):
        return None
    if b2a and any(exx != 1 for exx in b2a.values()):
        return None
    # collect common terms
    fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
    gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
    if n1 == 2:
        if len(b1a) == 2:
            # f(x)*f(y) + sign*g(x)*g(y) -> h1(x, y, sign)
            x, y = b1a.keys()
            hv = h1(x, y, sign)
        elif len(b2a) == 2:
            # g(x)*g(y) + sign*f(x)*f(y) -> sign*h1(x, y, sign)
            x, y = b2a.keys()
            hv = sign*h1(x, y, sign)
        else:
            # f(x)*g(y) + sign*g(x)*f(h) -> h2(x, y, sign)
            x = b1a.keys()[0]
            y = b2a.keys()[0]
            hv = h2(x, y, sign)
    elif n1 == 1:
        if len(b1a) == 1:
            # f(x) + sign*g(x) -> h3(x, sign)
            x = b1a.keys()[0]
            hv = h3(x, sign)
        else:
            # g(x) + sign*f(x) -> sign*h3(x, sign)
            x = b2a.keys()[0]
            hv = h3(x, sign)
    return fv, gv, hv

def replace_add_fgfg(expr, f, g, h1, h2, h3, full=True, to_tan=True):
    """
    Inverse of sum or difference of arguments

    Parameters
    ==========

    f, g, h1, h2, h3 : SymPy functions (or None in the case of ``h3``)
    full : if True continue to replace till possible

    Notes
    =====

    replace ``f(x)*f(y) + sign*g(y)*g(x) with h1(x, y, sign)``
    replace ``f(x)*g(y) + sign*f(y)*g(x) with h2(x, y, sign)``

    if ``h3`` is not None:
    replace ``f(x) + sign*g(x)``           with ``h3(x, y, sign)``
    """
    rule = _rule_replace_add_fgfg
    expr = replace_add_gen(expr, f, g, h1, h2, h3, None, rule, full)
    return expr

def replace_add1(expr, f, g, sgn, h1, h2, full=True, to_tan=True):
    """
    replace 1 - sgn*f(x)**2       -> g(x)**2
    replace 1 - g(x)**2           -> sgn*f(x)**2
    h1(x) = f(x)/g(x), h2(x) = g(x)/f(x) used in cases like
            1/f(x)**2 - 1         -> sgn*g(x)**2


    If ``to_tan`` is ``True``,
    use the replacements ``f(x)/g(x) -> h1(x)``, ``g(x)/f(x) -> h2(x)``
    in cases like
    ``f(x)**e1*g(x)**e2*(f(x)**2 - 1) -> -sgn*h2(x)**-e1*g(x)**e3``
    when ``e1 < 0`` and ``e3 = e1 + e2 + 2 >= 0``

    ``f, g`` are ``sin, cos`` or ``sinh, cosh``
    """
    rule = _rule_replace_add1
    expr = replace_add_gen(expr, f, g, h1, h2, None, sgn, rule, full, to_tan)
    return expr

def replace_add2(expr, f, g, sgn, full=True):
    """
    Replace f(x)**2 + sgn*g(x)**2 -> 1
    """
    rule = _rule_replace_add2
    expr = replace_add_gen(expr, f, g, None, None, None, sgn, rule, full)
    return expr

def _replace_mult_morrie(expr):
    """
    ``cos(x)*cos(2*x)*...*cos(2**(k-1)*x) -> sin(2**k*x)/(2**k*sin(x))``

    see http://en.wikipedia.org/wiki/Morrie%27s_law
    """
    if not expr.is_Mul:
        return expr
    #print 'DB0 expr=', expr
    c, d, _ = _splitfg(expr, cos, None)
    keys = d.keys()
    items = sorted(d.items())
    ks = sorted(keys)
    #print 'DB1 ks=', ks, c
    n = len(ks)
    used = set()
    args = [c]
    for i in range(n):
        if i in used:
            continue
        a = ks[i]
        e = d[a]
        v = [i]
        a1 = 2*a
        for j in range(i + 1, n):
            if j in used:
                continue
            if a1 == ks[j] and d[a1] == e:
                v.append(j)
                a1 *= 2
        if len(v) > 2:
            nv = len(v)
            used = used | set(v)
            args.append((sin(2**nv*a)/(2**nv*sin(a))**e))
        else:
            args.append(cos(a)**e)
    return Mul(*args)
