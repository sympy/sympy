from sympy.core import (Mul, Add, S)

def _splitfg(expr, f, g):
    """
    util function for ``_replace_add_get_a``
    """
    args = expr.args
    base_exp = []
    for x in args:
        if x.is_Pow:
            base_exp.append((x.base, x.exp))
        else:
            base_exp.append((x, S.One))
    a0, a1, a2 = [], [], []
    for base, exp in base_exp:
        if base.__class__ is f:
            a1.append((base, exp))
        elif base.__class__ is g:
            a2.append((base, exp))
        else:
            a0.append((base, exp))
    a0 = Mul(*[base**exp for base, exp in a0])
    a1 = [(x.args[0], y) for x, y in a1]
    a2 = [(x.args[0], y) for x, y in a2]
    return a0, a1, a2

def _replace_add_get_a(expr, f, g):
    """
    util function for ``replace_add_fgfg``
    """
    args = expr.args
    a = []
    # separate each addend into a part without products of f and g
    # at top level, in a part with powers of f and a part with powers of g
    for x in args:
        if not x.is_Mul:
            if x.__class__ is f:
                v = [S.One, [(x.args[0], 1)], []]
            elif x.__class__ is g:
                v = [S.One, [], [(x.args[0], 1)]]
            elif x.is_Pow:
                x, expx = x.base, x.exp
                if x.__class__ is f:
                    v = [S.One, [(x.args[0], expx)], []]
                elif x.__class__ is g:
                    v = [S.One, [], [(x.args[0], expx)]]
                else:
                    v = [x, [], []]
            else:
                v = [x, [], []]
            a.append(v)
            continue
        a0, a1, a2 = _splitfg(x, f, g)
        a.append([a0, a1, a2])
    a.sort(key=lambda y: abs(y[0])) # assume coefficients differ at most by sign
    return a

def fg_sub(p1, p2):
    """
    subtract two polynomials
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

def fg_factor(p1, p2):
    """
    factor the terms with minimum exponent between two polynomials
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

    p1 = fg_sub(p1, fact)
    p2 = fg_sub(p2, fact)
    return p1, p2, fact

def replace_add_fgfg(expr, f, g, h1, h2, full=False):
    """
    replace f(x)*f(y) + sign*g(y)*g(x) with h1(x, y, sign)
    replace f(x)*g(y) + sign*f(y)*g(x) with h2(x, y, sign)
    """
    if not expr.is_Add:
        return expr
    if f == g:
        return expr
    a = _replace_add_get_a(expr, f, g)
    # separate `a` in regions indexed by pos_a
    pos_a = []
    prev = S.Zero
    i = 0
    # TODO filter also with len(a1), len(a2)
    for i, (a0, a1, a2) in enumerate(a):
        if abs(a0) != abs(prev):
            prev = a0
            pos_a.append(i)
    if pos_a[-1] != len(a):
        pos_a.append(len(a))
    # convert to sets
    a_dicts = [[dict(y), dict(z)] for x, y, z in a]
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
                b1a, c1a, fact1 = fg_factor(b1, c1)
                b2a, c2a, fact2 = fg_factor(b2, c2)
                # f(x)*g(y) + sign*f(y)*g(x) -> h2(x, y, sign)
                # there must be two common terms, with exponent 1
                # f(x)*f(y) + sign*g(x)*g(y) -> h1(x, y, sign)
                if len(b1a) + len(b2a) != 2:
                    continue
                if len(c1a) + len(c2a) != 2:
                    continue
                if b1a != c2a or b2a != c1a:
                    continue
                if b1a and any(exx != 1 for exx in b1a.values()):
                    continue
                if b2a and any(exx != 1 for exx in b2a.values()):
                    continue
                # collect common terms
                fv = Mul(*[f(xx)**yy for xx, yy in fact1.items()])
                gv = Mul(*[g(xx)**yy for xx, yy in fact2.items()])
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
                found.append(r)
                used.add(i)
                used.add(j)
                break # end j loop
        pos0 = pos1
    if found:
        a_new = [a[i] for i in range(len(a)) if i not in used]
        if a_new:
            x0, x1, x2 = a_new[0]
        a_new = [Mul(*[x0, Mul(*[f(xx)**yy for xx, yy in x1]),
                 Mul(*[g(xx)**yy for xx, yy in x2])]) for x0, x1, x2 in a_new]
        a_new += found
        res = Add(*a_new)
        if full == True:
            res = replace_add_fgfg(res, f, g, h1, h2, full)
        return res
    else:
        return expr

