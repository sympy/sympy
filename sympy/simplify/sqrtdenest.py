from sympy.functions import sqrt, sign
from sympy.core import S, Wild, Rational, sympify, Mul, Add, Expr
from sympy.core.mul import prod
from sympy.core.function import expand_multinomial

def sqrtdenest(expr):
    """
    Denests sqrts in an expression that contain other square roots
    if possible, otherwise return the expr unchanged.

    This algorithm is based on
    <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    Examples
    ========

    >>> from sympy.simplify.sqrtdenest import sqrtdenest
    >>> from sympy import sqrt
    >>> sqrtdenest(sqrt(5 + 2 * sqrt(6)))
    sqrt(2) + sqrt(3)

    See also: unrad in sympy.solvers.solvers
    """
    expr = sympify(expr)
    if expr.is_Pow and expr.exp is S.Half: #If expr is a square root
        n, d = expr.as_numer_denom()
        if d is S.One:
            return _sqrtdenest(expr)
        else:
            return _sqrtdenest(n)/_sqrtdenest(d)
    elif isinstance(expr, Expr):
        args = expr.args
        if args:
            return expr.func(*[sqrtdenest(a) for a in args])
    return expr

def _sqrtdenest(expr):
    from sympy.simplify.simplify import radsimp
    if not expr.is_Pow or expr.exp != S.Half:
        val = None
    else:
        a = expr.base
        if not a.args:
            val = (a, S.Zero, S.Zero)
        val = sqrt_match(a)

    if val:
        a, b, r = val
        # try a quick denesting
        d2 = expand_multinomial(a**2 - b**2*r)
        if d2.is_Number and d2.is_positive and \
           max([sqrt_depth(a), sqrt_depth(b), sqrt_depth(r)]) >= 1:
            d = sqrt(d2)
            vad = a + d
            # if a = m1*sqrt(m) and d = m2*sqrt(m) then a*d is Number
            # and vad = (m1 + m2)*sqrt(m)
            if sqrt_depth(vad) < sqrt_depth(a) or (a*d).is_Number:
                vad1 = radsimp(1/vad)
                return (sqrt(vad/2) + sign(b)*sqrt((b**2*r*vad1/2).expand())).expand()
        else:
            d = sqrt(d2)
            vad = a + d
            vp0, vp1 = vad.as_content_primitive()
            rp0, rp1 = r.as_content_primitive()
            q = rp1/vp1
            if q.is_Number:
                c = (b**2 * q * rp0)/(2 * vp0)
                if sqrt_depth(r) > sqrt_depth(c):
                    z = (sqrt(vp0/2)*sqrt(vp1) + sign(b)*sqrt(c)).expand()
                    return z

    else:
        return expr
    z = _denester([radsimp(expr**2)], (a, b, r, d2), 0)[0]
    if z == None:
        return expr
    if z is expr or not z.is_Add:
        return z
    else:
        a = z.args
        a = [sqrtdenest(x) for x in a]
        expr = Add(*a)
        return expr

def sqrt_depth(p):
    """
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import sqrt_depth
    >>> sqrt_depth(1 + sqrt(2)*(1 + sqrt(3)))
    1
    >>> sqrt_depth(1 + sqrt(2)*sqrt(1 + sqrt(3)))
    2
    """
    if p.is_Atom:
        return 0
    elif p.is_Add or p.is_Mul:
        return max([sqrt_depth(x) for x in p.args])
    elif p.is_Pow and p.exp == S.Half:
        return sqrt_depth(p.base) + 1
    else:
        return 0

def sqrt_match(p):
    """return (a, b, r) for match p = a + b*sqrt(r) where
    sqrt(r) has maximal nested sqrt among addends of p

    Examples:
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import sqrt_match
    >>> sqrt_match(1 + sqrt(2) + sqrt(2)*sqrt(3) +  2*sqrt(1+sqrt(5)))
    (1 + sqrt(2) + sqrt(6), 2, 1 + sqrt(5))
    """
    p = p.expand()
    if p.is_Number:
        return (p, S.Zero, S.Zero)
    if p.is_Add:
        pargs = list(p.args)
        v = [(sqrt_depth(x), i) for i, x in enumerate(pargs)]
        if not v:
            return None
        nmax = max(v)
        if nmax[0] == 0:
            b, r = p.as_coeff_Mul()
            if p.is_Pow and p.exp == S.Half:
                return (S.Zero, b, p.base)
            else:
                return None
        n = nmax[1]
        p1 = pargs[n]
        del pargs[n]
        a = Add(*pargs)
        b, r = p1.as_coeff_Mul()
        res = (a, b, r**2)
    else:
        b, r = p.as_coeff_Mul()
        if r.is_Pow and r.exp == S.Half:
            res = (S.Zero, b, r**2)
        else:
            return None
    return res

def _denester (nested, av0, h, max_depth_level=4):
    """
    Denests a list of expressions that contain nested square roots.
    This method should not be called directly - use 'sqrtdenest' instead.
    This algorithm is based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    It is assumed that all of the elements of 'nested' share the same
    bottom-level radicand. (This is stated in the paper, on page 177, in
    the paragraph immediately preceding the algorithm.)

    When evaluating all of the arguments in parallel, the bottom-level
    radicand only needs to be denested once. This means that calling
    _denester with x arguments results in a recursive invocation with x+1
    arguments; hence _denester has polynomial complexity.

    However, if the arguments were evaluated separately, each call would
    result in two recursive invocations, and the algorithm would have
    exponential complexity.

    This is discussed in the paper in the middle paragraph of page 179.
    """
    from sympy.simplify.simplify import radsimp
    if h > max_depth_level:
        return None, None
    #if not av0 and all(n.is_Number and n.is_positive for n in nested): #If none of the arguments are nested
    if not av0 and all(n.is_Number for n in nested): #If none of the arguments are nested
        for f in subsets(len(nested)): #Test subset 'f' of nested
            p = prod(nested[i] for i in range(len(f)) if f[i]).expand()
            if 1 in f and f.count(1) > 1 and f[-1]:
                p = -p
            sqp = sqrt(p)
            if sqp.is_Number:
                return sqp, f #If we got a perfect square, return its square root.
        return sqrt(nested[-1]), [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
    else:
        R = None
        if av0:
            values = [av0[:2]]
            R = av0[2]
            nested2 = [av0[3], R]
        else:
            av0 = False
            values = filter(None, [sqrt_match(expr) for expr in nested])
            for v in values:
                if v[2]: #Since if b=0, r is not defined
                    if R is not None:
                        if R != v[2]:
                            return None, False
                    else:
                        R = v[2]
            if R is None:
                return sqrt(nested[-1]), [0]*len(nested) # return the radicand from the previous invocation
            nested2 = [(v[0]**2).expand()-(R*v[1]**2).expand() for v in values] + [R]
        d, f = _denester(nested2, False, h+1)
        if not f:
            return None, None
        if all(fi == 0 for fi in f):
            v = values[-1]
            return sqrt(v[0] + v[1]*d), f
        else:
            p = prod(nested[i] for i in range(len(nested)) if f[i])
            if p == 1:
                p = S(p)
            v = sqrt_match(p.expand())
            v = list(v)
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested)-1]:
                v[0] = -1 * v[0]
                v[1] = -1 * v[1]
            if not f[len(nested)]: #Solution denests with square roots
                vad = (v[0] + d).expand()
                if not vad:
                    return sqrt(nested[-1]), [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
                vad1 = radsimp(1/vad)
                return (sqrt(vad/2) + sign(v[1])*sqrt((v[1]**2*R*vad1/2).expand())).expand(), f
            else: #Solution requires a fourth root
                FR, s = (R.expand()**Rational(1,4)), sqrt((v[1]*R).expand()+d)
                if s == 0:
                    return sqrt(nested[-1]), [0]*len(nested)
                return (s/(sqrt(2)*FR) + v[0]*FR/(sqrt(2)*s)).expand(), f

def subsets(n):
    """
    Returns all possible subsets of the set (0, 1, ..., n-1) except the
    empty set, listed in reversed lexicographical order according to binary
    representation, so that the case of the fourth root is treated last.

    Examples
    ========

    >>> from sympy.simplify.sqrtdenest import subsets
    >>> subsets(3)
    [[0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

    """
    if n == 1:
        a = [[1]]
    elif n == 2:
        a = [[1, 0], [0, 1], [1, 1]]
    elif n == 3:
        a = [[1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    else:
        b = subsets(n-1)
        a0 = [x+[0] for x in b]
        a1 = [x+[1] for x in b]
        a = a0 + [[0]*(n-1) + [1]] + a1
    return a
