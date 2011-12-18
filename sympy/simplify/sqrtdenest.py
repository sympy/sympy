from sympy.functions import sqrt, sign
from sympy.core import S, Wild, Rational, sympify, Mul, Add, Expr
from sympy.core.function import expand_multinomial, expand_mul
from sympy.core.symbol import Dummy

def _mexpand(expr):
    return expand_mul(expand_multinomial(expr))

def _sqrt_symbolic_denest(a, b, r, d2=None):
    """
    Given an expression, sqrt(A + B*sqrt(R)), return the denested
    expression or None. This is called by _sqrtdenest when A**2 - B**2*R is
    not a Rational.

    Algorithm:
    If r = ra + rb*sqrt(rr), try replacing sqrt(rr) in A with (y**2 - ra)/rb,
    and if the result is a quadratic expression in y, see if it is a square,
    in which case write it in square form and take the square root

    Examples
    >>> from sympy.simplify.sqrtdenest import _sqrt_symbolic_denest, sqrtdenest
    >>> from sympy import sqrt
    >>> from sympy.abc import x

    >>> _sqrt_symbolic_denest(16 - 2*sqrt(29), 2, -10*sqrt(29) + 55)
    sqrt(-2*sqrt(29) + 11) + sqrt(5)

    if sympy decides that sqrt and square can be simplified, it does,
    otherwise it leaves it in that form

    >>> w = 1 + sqrt(2) + sqrt(1 + sqrt(1+sqrt(3)))
    >>> sqrtdenest(sqrt((w**2).expand()))
    1 + sqrt(2) + sqrt(1 + sqrt(1 + sqrt(3)))
    >>> w = 1 + sqrt(2) + sqrt(1 + sqrt(1+sqrt(3+x)))
    >>> sqrtdenest(sqrt((w**2).expand()))
    sqrt((sqrt(sqrt(sqrt(x + 3) + 1) + 1) + 1 + sqrt(2))**2)
    """
    a, b, r = S(a), S(b), S(r)
    if d2 == None:
        d2 = _mexpand(a**2 - b**2*r)
    rval = sqrt_match(r)
    if rval is None:
        return None
    ra, rb, rr = rval
    if rb != 0:
        y = Dummy('y', positive=True)
        y2 = y**2
        a2 = a.subs(sqrt(rr), (y**2 - ra)/rb).expand()
        cav, cbv, ccv = [], [], []
        for xx in a2.args:
            cx, qx = xx.as_coeff_Mul()
            if qx.is_Mul:
                qxa = list(qx.args)
                if y in qxa:
                    qxa.remove(y)
                    cbv.append( Mul(*(qxa+[cx])) )
                elif y2 in qxa:
                    qxa.remove(y2)
                    cav.append(Mul(*(qxa+[cx])))
                else:
                    ccv.append(xx)
            elif qx == y2:
                cav.append(cx)
            else:
                if y == qx:
                    cbv.append( cx )
                elif y2 == qx:
                    cav.append(cx)
                else:
                    ccv.append(xx)
        ca = Add(*cav)
        cb = Add(*cbv)
        cc = Add(*ccv)
        if ca and not cc.has(y):
            cb += b
            discr = _mexpand(cb**2 - 4*ca*cc)
            if not discr:
                z = sqrt(ca*(sqrt(r) + cb/(2*ca))**2)
                c, q = z.as_content_primitive()
                z = c*q
                if z.is_number:
                    z = _mexpand(z)
                return z


def sqrt_numeric_denest(a, b, r, d2):
    """
    denest expr = a + b*sqrt(r), with d2 = a**2 - b**2*r > 0

    If not denested return None
    """
    from sympy.simplify.simplify import radsimp
    depthr = sqrt_depth(r)
    d = sqrt(d2)
    vad = a + d
    # sqrt_depth(res) <= sqrt_depth(vad) + 1
    # sqrt_depth(expr) = depthr + 2
    # there is denesting if sqrt_depth(vad)+1 < depthr + 2
    # if vad**2 is Number there is a fourth root
    if sqrt_depth(vad) < depthr + 1 or (vad**2).is_Rational:
        vad1 = radsimp(1/vad)
        return (sqrt(vad/2) + sign(b)*sqrt((b**2*r*vad1/2).expand())).expand()


def _sqrt_four_terms_denest(expr):
    """denest the square root of three or four square root of rationals

    See D.J.Jeffrey and A.D.Rich
    'Symplifying Square Roots of Square Roots by Denesting'

    Examples
    >>> from sympy import sqrt
    >>> from sympy.simplify.sqrtdenest import _sqrt_four_terms_denest
    >>> _sqrt_four_terms_denest(sqrt(-72*sqrt(2) + 158*sqrt(5) + 498))
    -sqrt(10) + sqrt(2) + 9 + 9*sqrt(5)
    >>> _sqrt_four_terms_denest(sqrt(12+2*sqrt(6)+2*sqrt(14)+2*sqrt(21)))
    sqrt(2) + sqrt(3) + sqrt(7)
    """
    from sympy.simplify.simplify import radsimp
    if not (expr.is_Pow and expr.exp == S.Half):
        return expr
    if expr.base < 0:
        return sqrt(-1)*_sqrt_four_terms_denest(sqrt(-expr.base))
    a = Add(*expr.base.args[:2])
    b = Add(*expr.base.args[2:])
    if a < 0:
        a, b = b, a
    d2 = _mexpand(a**2 - b**2)
    if d2 < 0:
        a, b = b, a
        d2 = -d2
    d = _sqrtdenest(sqrt(d2))
    if sqrt_depth(d) > 1:
        return expr
    vad = a + d
    c = _sqrtdenest(sqrt(vad))
    if sqrt_depth(c) > 1:
        return expr
    return radsimp((c/sqrt(2) + b/(sqrt(2)*c))).expand()


def sqrtdenest(expr, max_iter=3):
    """
    Denests sqrts in an expression that contain other square roots
    if possible, otherwise return the expr unchanged.

    This algorithm is based on
    <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    Examples
    >>> from sympy.simplify.sqrtdenest import sqrtdenest
    >>> from sympy import sqrt
    >>> sqrtdenest(sqrt(5 + 2 * sqrt(6)))
    sqrt(2) + sqrt(3)

    See also: unrad in sympy.solvers.solvers
    """
    expr = expand_mul(sympify(expr))
    for i in range(max_iter):
        z = sqrtdenest0(expr)
        if expr == z:
            return expr
        expr = z
    return expr


def sqrtdenest0(expr):
    if expr.is_Atom:
        return expr
    if expr.is_Pow and expr.exp is S.Half: #If expr is a square root
        n, d = expr.as_numer_denom()
        if d is S.One:
            nn = len(n.base.args)
            if 3 <= nn <= 4 and all([(x**2).is_Rational for x in n.base.args]):
                return _sqrt_four_terms_denest(n)
            if n.base.is_Add:
                expr = sqrt(_mexpand(Add(*[sqrtdenest0(x) for x in n.base.args])))
            return _sqrtdenest(expr)
        else:
            nn = len(n.base.args)
            if 3 <= nn <= 4 and all([(x**2).is_Rational for x in n.base.args]):
                n1 = _sqrt_four_terms_denest(n)
            else:
                n1 = _sqrtdenest(n)

            if d.is_Pow and d.exp == S.Half and \
                3 <= len(d.base.args) <= 4 and all([(x**2).is_Rational for x in d.base.args]):
                d1 = _sqrt_four_terms_denest(d)
            else:
                d1 = _sqrtdenest(d)
            return n1/d1

    elif expr.is_Pow and expr.exp == -S.Half:
        return 1/sqrtdenest0(sqrt(expr.base))
    elif expr.is_Mul:
        return Mul(*[sqrtdenest0(x) for x in expr.args])
    elif isinstance(expr, Expr):
        args = expr.args
        if args:
            return expr.func(*[sqrtdenest0(a) for a in args])
    return expr

def _sqrtdenest(expr):
    from sympy.simplify.simplify import radsimp
    if not expr.is_Pow or expr.exp != S.Half:
        val = None
    else:
        a = expr.base
        if a.is_Number:
            return expr
        if not a.args:
            val = (a, S.Zero, S.Zero)
        val = sqrt_match(a)

    if val:
        a, b, r = val
        # try a quick numeric denesting
        d2 = _mexpand(a**2 - b**2*r)
        if d2.is_Rational:
            if d2.is_positive:
                z = sqrt_numeric_denest(a, b, r, d2)
                if z != None:
                    return z
            else:
                # d2 negative
                # fourth root case
                # sqrtdenest(sqrt(3 + 2*sqrt(3))) =
                # sqrt(2)*3**(1/4)/2 + sqrt(2)*3**(3/4)/2
                dr2 = _mexpand(-d2*r)
                dr = sqrt(dr2)
                if dr.is_Rational:
                    z = sqrt_numeric_denest(_mexpand(b*r), a, r, dr2)
                    if z != None:
                        return z/r**Rational(1,4)

        else:
            z = _sqrt_symbolic_denest(a, b, r, d2)
            if z != None:
                return z

    else:
        return expr
    if not expr.is_number:
        return expr
    if not is_algebraic(expr):
        return expr
    av0 = [a, b, r, d2]
    z = _denester([radsimp(expr**2)], av0, 0, sqrt_depth(expr)-1)[0]
    if av0[1] == None:
        return expr
    if z != None:
        return z
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

def is_algebraic(p):
    """
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import is_algebraic
    >>> from sympy import cos
    >>> is_algebraic(sqrt(2)*(3/(sqrt(7) + sqrt(5)*sqrt(2))))
    True
    >>> is_algebraic(sqrt(2)*(3/(sqrt(7) + sqrt(5)*cos(2))))
    False
    """
    if p.is_Rational:
        return True
    elif p.is_Atom:
        return False
    elif p.is_Pow and (p.exp == S.Half or p.exp == -1):
        return is_algebraic(p.base)
    elif p.is_Add or p.is_Mul:
        return all([is_algebraic(x) for x in p.args])
    else:
        return False

def sqrt_match(p):
    """return (a, b, r) for match p = a + b*sqrt(r) where
    sqrt(r) has maximal nested sqrt among addends of p

    # FIXME should one count also fourth roots, or maybe any root?

    Examples:
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import sqrt_match
    >>> sqrt_match(1 + sqrt(2) + sqrt(2)*sqrt(3) +  2*sqrt(1+sqrt(5)))
    (1 + sqrt(2) + sqrt(6), 2, 1 + sqrt(5))
    """
    p = _mexpand(p)
    if p.is_Number:
        return (p, S.Zero, S.Zero)
    if p.is_Add:
        pargs = list(p.args)
        v = [(sqrt_depth(x), i) for i, x in enumerate(pargs)]
        nmax = max(v)
        if nmax[0] == 0:
            b, r = p.as_coeff_Mul()
            if p.is_Pow and p.exp == S.Half:
                return (S.Zero, b, p.base)
            else:
                return None
        depth = nmax[0]
        n = nmax[1]
        p1 = pargs[n]
        del pargs[n]
        a = Add(*pargs)
        bv = []
        rv = []
        if p1.is_Mul:
            for x in p1.args:
                if sqrt_depth(x) < depth:
                    bv.append(x)
                else:
                    rv.append(x)

            b = Mul(*bv)
            r = Mul(*rv)
        else:
            b = S.One
            r = p1
        res = (a, b, r**2)
    else:
        b, r = p.as_coeff_Mul()
        if r.is_Pow and r.exp == S.Half:
            res = (S.Zero, b, r**2)
        else:
            return None
    return res

def _denester (nested, av0, h, max_depth_level):
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
    if av0[1] == None:
        return None, None
    #If none of the arguments are nested
    if av0[0] == None and all(n.is_Number for n in nested): #If none of the arguments are nested
        for f in subsets(len(nested)): #Test subset 'f' of nested
            p = _mexpand(Mul(*[nested[i] for i in range(len(f)) if f[i]]))
            if 1 in f and f.count(1) > 1 and f[-1]:
                p = -p
            sqp = sqrt(p)
            if sqp.is_Number:
                return sqp, f #If we got a perfect square, return its square root.
        return sqrt(nested[-1]), [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
    else:
        R = None
        if av0[0] != None:
            values = [av0[:2]]
            R = av0[2]
            nested2 = [av0[3], R]
            av0[0] = None
        else:
            values = filter(None, [sqrt_match(expr) for expr in nested])
            for v in values:
                if v[2]: #Since if b=0, r is not defined
                    if R is not None:
                        if R != v[2]:
                            av0[1] = None
                            return None, None
                        #assert R == v[2] #All the 'r's should be the same.
                    else:
                        R = v[2]
            if R is None:
                return sqrt(nested[-1]), [0]*len(nested) # return the radicand from the previous invocation
            nested2 = [_mexpand(v[0]**2) - _mexpand(R*v[1]**2) for v in values] + [R]
        d, f = _denester(nested2, av0, h+1, max_depth_level)
        if not f:
            return None, None
        if not any(f[i] for i in range(len(nested))):
            v = values[-1]
            return sqrt(v[0] + v[1]*d), f
        else:
            p = Mul(*[nested[i] for i in range(len(nested)) if f[i]])
            v = sqrt_match(_mexpand(p))
            v = list(v)
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested)-1]:
                v[0] = -1 * v[0]
                v[1] = -1 * v[1]
            if not f[len(nested)]: #Solution denests with square roots
                vad = _mexpand(v[0] + d)
                if vad <= 0:
                    return sqrt(nested[-1]), [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
                if not(sqrt_depth(vad) < sqrt_depth(R) + 1 or (vad**2).is_Number):
                    av0[1] = None
                    return None, None

                vad1 = radsimp(1/vad)
                return _mexpand(sqrt(vad/2) + sign(v[1])*sqrt(_mexpand(v[1]**2*R*vad1/2))), f
            else: #Solution requires a fourth root
                s2 = _mexpand(v[1]*R) + d
                if s2 <= 0:
                    return sqrt(nested[-1]), [0]*len(nested)
                FR, s = (_mexpand(R)**Rational(1,4)), sqrt(s2)
                return _mexpand(s/(sqrt(2)*FR) + v[0]*FR/(sqrt(2)*s)), f

def subsets(n):
    """
    Returns all possible subsets of the set of n elements except the empty,
    listed in reversed lexicographical order according to binary
    representation, so that the case of the fourth root is treated last.
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
