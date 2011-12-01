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
    val = sqrt_match(expr)
    if val:
        a, b, r = val
        # try a quick denesting
        d2 = expand_multinomial(a**2 - b**2*r)
        if d2.is_Number and \
            max([sqrt_depth(a), sqrt_depth(b), sqrt_depth(r)]) >= 1:
            d = sqrt(d2)
            vad = a + d
            vad1 = radsimp(1/vad)
            return (sqrt(vad/2) + sign(b)*sqrt((b**2*r*vad1/2).expand())).expand()

    z = denester([expr], 0)[0]
    if z is expr or not z.is_Add:
        return z
    else:
        a = z.args
        a = [sqrtdenest(x) for x in a]
        expr = Add(*a)
        return expr

def sqrt_match(p):
    if not p.is_Pow or p.args[1] != S.Half:
        return None
    else:
        a = p.args[0]
        if not a.args:
            return (a, S.Zero, S.Zero)
        return sqrt_match0(a)

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

def sqrt_match0(p):
    """return (a, b, r) for match a + b*sqrt(r) where
    sqrt(r) has maximal nested sqrt among addends of p

    Examples:
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import sqrt_match0
    >>> sqrt_match0(1 + sqrt(2) + sqrt(2)*sqrt(3) +  2*sqrt(1+sqrt(5)))
    (1 + sqrt(2) + sqrt(6), 2, 1 + sqrt(5))
    """
    p = p.expand()
    if p.is_Number:
        return (p, S.Zero, S.Zero)
    pargs = list(p.args)
    v = [(sqrt_depth(x), i) for i, x in enumerate(pargs)]
    if not v:
        return None
    nmax = max(v)
    if nmax[0] == 0:
        return None
    n = nmax[1]
    p1 = pargs[n]
    del pargs[n]
    a = Add(*pargs)
    b, r = p1.as_coeff_Mul()
    return (a, b, r**2)

def denester (nested, h):
    """
    Denests a list of expressions that contain nested square roots.
    This method should not be called directly - use 'denest' instead.
    This algorithm is based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    It is assumed that all of the elements of 'nested' share the same
    bottom-level radicand. (This is stated in the paper, on page 177, in
    the paragraph immediately preceding the algorithm.)

    When evaluating all of the arguments in parallel, the bottom-level
    radicand only needs to be denested once. This means that calling
    denester with x arguments results in a recursive invocation with x+1
    arguments; hence denester has polynomial complexity.

    However, if the arguments were evaluated separately, each call would
    result in two recursive invocations, and the algorithm would have
    exponential complexity.

    This is discussed in the paper in the middle paragraph of page 179.
    """
    from sympy.simplify.simplify import radsimp
    if all((n**2).is_Number for n in nested): #If none of the arguments are nested
        for f in subsets(len(nested)): #Test subset 'f' of nested
            p = prod(nested[i]**2 for i in range(len(f)) if f[i]).expand()
            if 1 in f and f.count(1) > 1 and f[-1]:
                p = -p
            if sqrt(p).is_Number:
                return sqrt(p), f #If we got a perfect square, return its square root.
        return nested[-1], [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
    else:
        R = None
        values = filter(None, [sqrt_match(expr) for expr in nested])
        for v in values:
            if v[2]: #Since if b=0, r is not defined
                if R is not None:
                    assert R == v[2] #All the 'r's should be the same.
                else:
                    R = v[2]
        if R is None:
            return nested[-1], [0]*len(nested) # return the radicand from the pravious invocation
        d, f = denester([sqrt((v[0]**2).expand()-(R*v[1]**2).expand()) for v in values] + [sqrt(R)], h+1)
        if not any(f[i] for i in range(len(nested))):
        #if all(fi == 0 for fi in f):
            v = values[-1]
            return sqrt(v[0] + v[1]*d), f
        else:
            p = prod(nested[i]**2 for i in range(len(nested)) if f[i])
            if p == 1:
                p = S(p)
            v = sqrt_match0(p.expand())
            v = list(v)
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested)-1]:
                v[0] = -1 * v[0]
                v[1] = -1 * v[1]
            if not f[len(nested)]: #Solution denests with square roots
                vad = (v[0] + d).expand()
                if not vad:
                    return nested[-1], [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
                vad1 = radsimp(1/vad)
                return (sqrt(vad/2) + sign(v[1])*sqrt((v[1]**2*R*vad1/2).expand())).expand(), f
            else: #Solution requires a fourth root
                FR, s = (R.expand()**Rational(1,4)), sqrt((v[1]*R).expand()+d)
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
