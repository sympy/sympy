from sympy.core import (Mul, Add, S)
from sympy.utilities.misc import debug


def _splitfg(expr, f, g):
    """
    util function for ``_replace_add_get_a``

    Notes
    =====

    ``expr = c1*f(x11)**e11*f(x12)**e12*...*g(x21)**e21*g(x22)**e22``
    This function gives ``(c1, {x11:e11, x12:e12,...}, {x21:e21, x22:e22,...})
    """
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
    util function for ``replace_add_fgfg``

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

    In ``replace_add_fgfg(expr, f, g, h1, h2)`` a term
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

def replace_add_fgfg(expr, f, g, h1, h2, h3, full=True):
    """Inverse of sum or difference of arguments

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
    if not expr.is_Add:
        return expr
    if f == g:
        return expr
    a  = _replace_add_get_a(expr.args, f, g)
    # Matching can occur only for terms with the same coefficient modulo the sign.
    # Sort ``a`` according to the absolute value of the coefficients,
    # reorder args accordingly
    a = zip(a, expr.args)
    a.sort(key=lambda y: abs(y[0][0]))
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
                # f(x)*g(y) + sign*f(y)*g(x) -> h2(x, y, sign)
                # there must be two common terms, with exponent 1
                # f(x)*f(y) + sign*g(x)*g(y) -> h1(x, y, sign)
                # terms with two ``f`` or two ``g``
                n1 = len(b1a) + len(b2a)
                if n1 not in (1, 2):
                    continue
                n2 = len(c1a) + len(c2a)
                if n2 not in (1, 2):
                    continue
                if n1 == 1 and not h3:
                    continue

                if b1a != c2a or b2a != c1a:
                    continue
                # The replacement is done only if ``x`` and ``y``
                # are different, so that the exponents ``exx``
                # are equal to 1.
                if b1a and any(exx != 1 for exx in b1a.values()):
                    continue
                if b2a and any(exx != 1 for exx in b2a.values()):
                    continue
                # collect common terms
                fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
                gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
                if n1 == 2:
                    if len(b1a) == 2:
                        # f(x)*f(y) + sign*g(x)*g(y) -> h1(x, y, sign)
                        x, y = b1a.keys()
                        hv = h1(x, y, sign)
                        r = Mul(a[i][0], fv, gv, hv)
                    elif len(b2a) == 2:
                        # g(x)*g(y) + sign*f(x)*f(y) -> sign*h1(x, y, sign)
                        x, y = b2a.keys()
                        hv = h1(x, y, sign)
                        r = Mul(sign, a[i][0], fv, gv, hv)
                    else:
                        # f(x)*g(y) + sign*g(x)*f(h) -> h2(x, y, sign)
                        x = b1a.keys()[0]
                        y = b2a.keys()[0]
                        hv = h2(x, y, sign)
                        r = Mul(a[i][0], fv, gv, hv)
                elif n1 == 1:
                    if len(b1a) == 1:
                        # f(x) + sign*g(x) -> h3(x, sign)
                        x = b1a.keys()[0]
                        hv = h3(x, sign)
                        r = Mul(a[i][0], fv, gv, hv)
                    else:
                        # g(x) + sign*f(x) -> sign*h3(x, sign)
                        x = b2a.keys()[0]
                        hv = sign*h3(x, sign)
                        r = Mul(sign, a[i][0], fv, gv, hv)
                debug('%s + %s -> %s' %(args[i], args[j], r))
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
            res = replace_add_fgfg(res, f, g, h1, h2, h3, full)
        return res
    else:
        return expr
