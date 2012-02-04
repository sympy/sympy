from sympy.functions import sqrt, sign, root
from sympy.core import S, Wild, sympify, Mul, Add, Expr
from sympy.core.function import expand_multinomial, expand_mul
from sympy.core.symbol import Dummy
from sympy.polys import Poly, PolynomialError

def _mexpand(expr):
    return expand_mul(expand_multinomial(expr))


def is_sqrt(expr):
    """Return True if expr is a sqrt, otherwise False."""

    return expr.is_Pow and expr.exp.is_Rational and abs(expr.exp) is S.Half


def sqrt_depth(p):
    """Return the maximum depth of any square root argument of p.

    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import sqrt_depth

    Neither of these square roots contains any other square roots
    so the depth is 1:

    >>> sqrt_depth(1 + sqrt(2)*(1 + sqrt(3)))
    1

    The sqrt(3) is contained within a square root so the depth is
    2:

    >>> sqrt_depth(1 + sqrt(2)*sqrt(1 + sqrt(3)))
    2
    """

    if p.is_Atom:
        return 0
    elif p.is_Add or p.is_Mul:
        return max([sqrt_depth(x) for x in p.args])
    elif is_sqrt(p):
        return sqrt_depth(p.base) + 1
    else:
        return 0


def is_algebraic(p):
    """Return True if p is comprised of only Rationals or square roots
    of Rationals and algebraic operations.

    Examples
    ========
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
    elif is_sqrt(p) or p.is_Pow and p.exp.is_Integer:
        return is_algebraic(p.base)
    elif p.is_Add or p.is_Mul:
        return all(is_algebraic(x) for x in p.args)
    else:
        return False


def subsets(n):
    """
    Returns all possible subsets of the set (0, 1, ..., n-1) except the
    empty set, listed in reversed lexicographical order according to binary
    representation, so that the case of the fourth root is treated last.

    Examples
    ========

    >>> from sympy.simplify.sqrtdenest import subsets
    >>> subsets(2)
    [[1, 0], [0, 1], [1, 1]]

    """
    if n == 1:
        a = [[1]]
    elif n == 2:
        a = [[1, 0], [0, 1], [1, 1]]
    elif n == 3:
        a = [[1, 0, 0], [0, 1, 0], [1, 1, 0],
             [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    else:
        b = subsets(n-1)
        a0 = [x+[0] for x in b]
        a1 = [x+[1] for x in b]
        a = a0 + [[0]*(n-1) + [1]] + a1
    return a


def sqrtdenest(expr, max_iter=3):
    """Denests sqrts in an expression that contain other square roots
    if possible, otherwise returns the expr unchanged. This is based on the
    algorithms of [1].

    Examples
    ========

    >>> from sympy.simplify.sqrtdenest import sqrtdenest
    >>> from sympy import sqrt
    >>> sqrtdenest(sqrt(5 + 2 * sqrt(6)))
    sqrt(2) + sqrt(3)

    See Also
    ========
    sympy.solvers.solvers.unrad

    References
    ==========
    [1] http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf
    """
    expr = expand_mul(sympify(expr))
    for i in range(max_iter):
        z = _sqrtdenest0(expr)
        if expr == z:
            return expr
        expr = z
    return expr


def _sqrt_match(p):
    """Return [a, b, r] for p.match(a + b*sqrt(r)) where, in addition to
    matching, sqrt(r) also has then maximal sqrt_depth among addends of p.

    Examples
    ========
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.simplify.sqrtdenest import _sqrt_match
    >>> _sqrt_match(1 + sqrt(2) + sqrt(2)*sqrt(3) +  2*sqrt(1+sqrt(5)))
    [1 + sqrt(2) + sqrt(6), 2, 1 + sqrt(5)]
    """

    p = _mexpand(p)
    if p.is_Number:
        res = (p, S.Zero, S.Zero)
    elif p.is_Add:
        pargs = list(p.args)
        # to make the process canonical, the argument is included in the tuple
        # so when the max is selected, it will be the largest arg having a
        # given depth
        v = [(sqrt_depth(x), x, i) for i, x in enumerate(pargs)]
        nmax = max(v)
        if nmax[0] == 0:
            res = []
        else:
            depth, _, i = nmax
            r = pargs.pop(i)
            a = Add._from_args(pargs)
            b = S.One
            if r.is_Mul:
                bv = []
                rv = []
                for x in r.args:
                    if sqrt_depth(x) < depth:
                        bv.append(x)
                    else:
                        rv.append(x)

                b = Mul._from_args(bv)
                r = Mul._from_args(rv)
            res = (a, b, r**2)
    else:
        b, r = p.as_coeff_Mul()
        if is_sqrt(r):
            res = (S.Zero, b, r**2)
        else:
            res = []
    return list(res)

class SqrtdenestStopIteration(StopIteration):
    pass

def _sqrtdenest0(expr):
    """Returns expr after denesting its arguments."""

    if is_sqrt(expr):
        n, d = expr.as_numer_denom()
        if d is S.One: # n is a square root
            if n.base.is_Add:
                args = n.base.args
                if len(args) > 2 and all((x**2).is_Integer for x in args):
                    try:
                        return _sqrtdenest_rec(n)
                    except SqrtdenestStopIteration:
                        pass
                expr = sqrt(_mexpand(Add(*[_sqrtdenest0(x) for x in args])))
            return _sqrtdenest1(expr)
        else:
            n, d = [_sqrtdenest0(i) for i in (n, d)]
            return n/d
    if isinstance(expr, Expr):
        args = expr.args
        if args:
            return expr.func(*[_sqrtdenest0(a) for a in args])
    return expr

def _sqrtdenest_rec(expr):
    """Helper that denests the square root of three or more surds.

    It returns the denested expression; if it cannot be denested it
    throws SqrtdenestStopIteration

    Algorithm: expr.base is in the extension Q_m = Q(sqrt(r_1),..,sqrt(r_k));
    split expr.base = a + b*sqrt(r_k), where `a` and `b` are on
    Q_(m-1) = Q(sqrt(r_1),..,sqrt(r_(k-1))); then a**2 - b**2*r_k is
    on Q_(m-1); denest sqrt(a**2 - b**2*r_k) and so on.
    See [1], section 6.

    Examples
    ========
    >>> from sympy import sqrt
    >>> from sympy.simplify.sqrtdenest import _sqrtdenest_rec
    >>> _sqrtdenest_rec(sqrt(-72*sqrt(2) + 158*sqrt(5) + 498))
    -sqrt(10) + sqrt(2) + 9 + 9*sqrt(5)
    >>> w=-6*sqrt(55)-6*sqrt(35)-2*sqrt(22)-2*sqrt(14)+2*sqrt(77)+6*sqrt(10)+65
    >>> _sqrtdenest_rec(sqrt(w))
    -sqrt(11) - sqrt(7) + sqrt(2) + 3*sqrt(5)
    """
    from sympy.simplify.simplify import radsimp, split_surds, rad_rationalize
    if expr.base < 0:
        return sqrt(-1)*_sqrtdenest_rec(sqrt(-expr.base))
    a, b = split_surds(expr.base)
    if a < b:
        a, b = b, a
    c2 = _mexpand(a**2 - b**2)
    if len(c2.args) > 2:
        a1, b1 = split_surds(c2)
        if a1 < b1:
            a1, b1 = b1, a1
        c2_1 = _mexpand(a1**2 - b1**2)
        c_1 = _sqrtdenest_rec(sqrt(c2_1))
        d_1 = _sqrtdenest_rec(sqrt(a1 + c_1))
        num, den = rad_rationalize(b1, d_1)
        c = _mexpand(d_1/sqrt(2) + num/(den*sqrt(2)))
    else:
        c = _sqrtdenest1(sqrt(_mexpand(a**2 - b**2)))

    if sqrt_depth(c) > 1:
        raise SqrtdenestStopIteration
    d = sqrtdenest(sqrt(a + c))
    if sqrt_depth(d) > 1:
        raise SqrtdenestStopIteration
    num, den = rad_rationalize(b, d)
    r = d/sqrt(2) + num/(den*sqrt(2))
    r = radsimp(r)
    return _mexpand(r)

def _sqrtdenest1(expr):
    """Return denested expr after denesting with simpler methods or, that
    failing, using the denester."""

    from sympy.simplify.simplify import radsimp

    if not is_sqrt(expr):
        return expr

    a = expr.base
    if a.is_Atom:
        return expr
    val = _sqrt_match(a)
    if not val:
        return expr

    a, b, r = val
    # try a quick numeric denesting
    d2 = _mexpand(a**2 - b**2*r)
    if d2.is_Rational:
        if d2.is_positive:
            z = _sqrt_numeric_denest(a, b, r, d2)
            if z is not None:
                return z
        else:
            # fourth root case
            # sqrtdenest(sqrt(3 + 2*sqrt(3))) =
            # sqrt(2)*3**(1/4)/2 + sqrt(2)*3**(3/4)/2
            dr2 = _mexpand(-d2*r)
            dr = sqrt(dr2)
            if dr.is_Rational:
                z = _sqrt_numeric_denest(_mexpand(b*r), a, r, dr2)
                if z is not None:
                    return z/root(r, 4)

    else:
        z = _sqrt_symbolic_denest(a, b, r)
        if z is not None:
            return z

    if not is_algebraic(expr):
        return expr

    # now call to the denester
    av0 = [a, b, r, d2]
    z = _denester([radsimp(expr**2)], av0, 0, sqrt_depth(expr) - 1)[0]
    if av0[1] is None:
        return expr
    if z is not None:
        return z
    return expr


def _sqrt_symbolic_denest(a, b, r):
    """Given an expression, sqrt(a + b*sqrt(b)), return the denested
    expression or None.

    Algorithm:
    If r = ra + rb*sqrt(rr), try replacing sqrt(rr) in ``a`` with
    (y**2 - ra)/rb, and if the result is a quadratic, ca*y**2 + cb*y + cc, and
    (cb + b)**2 - 4*ca*cc is 0, then sqrt(a + b*sqrt(r)) can be rewritten as
    sqrt(ca*(sqrt(r) + (cb + b)/(2*ca))**2).

    Examples
    ========
    >>> from sympy.simplify.sqrtdenest import _sqrt_symbolic_denest, sqrtdenest
    >>> from sympy import sqrt, Symbol, Poly
    >>> from sympy.abc import x

    >>> a, b, r = 16 - 2*sqrt(29), 2, -10*sqrt(29) + 55
    >>> _sqrt_symbolic_denest(a, b, r)
    sqrt(-2*sqrt(29) + 11) + sqrt(5)

    If the expression is numeric, it will be simplified:

    >>> w = sqrt(sqrt(sqrt(3) + 1) + 1) + 1 + sqrt(2)
    >>> sqrtdenest(sqrt((w**2).expand()))
    1 + sqrt(2) + sqrt(1 + sqrt(1 + sqrt(3)))

    Otherwise, it will only be simplified if assumptions allow:

    >>> w = w.subs(sqrt(3), sqrt(x + 3))
    >>> sqrtdenest(sqrt((w**2).expand()))
    sqrt((sqrt(sqrt(sqrt(x + 3) + 1) + 1) + 1 + sqrt(2))**2)

    Notice that the argument of the sqrt is a square. If x is made positive
    then the sqrt of the square is resolved:

    >>> _.subs(x, Symbol('x', positive=True))
    sqrt(sqrt(sqrt(x + 3) + 1) + 1) + 1 + sqrt(2)
    """

    a, b, r = sympify([a, b, r])
    rval = _sqrt_match(r)
    if not rval:
        return None
    ra, rb, rr = rval
    if rb:
        y = Dummy('y', positive=True)
        try:
            newa = Poly(a.subs(sqrt(rr), (y**2 - ra)/rb), y)
        except PolynomialError:
            return None
        if newa.degree() == 2:
            ca, cb, cc = newa.all_coeffs()
            cb += b
            if _mexpand(cb**2 - 4*ca*cc).equals(0):
                z = sqrt(ca*(sqrt(r) + cb/(2*ca))**2)
                if z.is_number:
                    z = _mexpand(Mul._from_args(z.as_content_primitive()))
                return z


def _sqrt_numeric_denest(a, b, r, d2):
    """Helper that denest expr = a + b*sqrt(r), with d2 = a**2 - b**2*r > 0
    or returns None if not denested.
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


def _denester(nested, av0, h, max_depth_level):
    """Denests a list of expressions that contain nested square roots.

    Algorithm based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

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
    if av0[1] is None:
        return None, None
    if (av0[0] is None and
        all(n.is_Number for n in nested)): # no arguments are nested
        for f in subsets(len(nested)): # test subset 'f' of nested
            p = _mexpand(Mul(*[nested[i] for i in range(len(f)) if f[i]]))
            if f.count(1) > 1 and f[-1]:
                p = -p
            sqp = sqrt(p)
            if sqp.is_Rational:
                return sqp, f # got a perfect square so return its square root.
        # Otherwise, return the radicand from the previous invocation.
        return sqrt(nested[-1]), [0]*len(nested)
    else:
        R = None
        if av0[0] is not None:
            values = [av0[:2]]
            R = av0[2]
            nested2 = [av0[3], R]
            av0[0] = None
        else:
            values = filter(None, [_sqrt_match(expr) for expr in nested])
            for v in values:
                if v[2]: #Since if b=0, r is not defined
                    if R is not None:
                        if R != v[2]:
                            av0[1] = None
                            return None, None
                    else:
                        R = v[2]
            if R is None:
                # return the radicand from the previous invocation
                return sqrt(nested[-1]), [0]*len(nested)
            nested2 = [_mexpand(v[0]**2) -
                       _mexpand(R*v[1]**2) for v in values] + [R]
        d, f = _denester(nested2, av0, h + 1, max_depth_level)
        if not f:
            return None, None
        if not any(f[i] for i in range(len(nested))):
            v = values[-1]
            return sqrt(v[0] + v[1]*d), f
        else:
            p = Mul(*[nested[i] for i in range(len(nested)) if f[i]])
            v = _sqrt_match(p)
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested) - 1]:
                v[0] = -v[0]
                v[1] = -v[1]
            if not f[len(nested)]: #Solution denests with square roots
                vad = _mexpand(v[0] + d)
                if vad <= 0:
                    # return the radicand from the previous invocation.
                    return sqrt(nested[-1]), [0]*len(nested)
                if not(sqrt_depth(vad) < sqrt_depth(R) + 1 or
                       (vad**2).is_Number):
                    av0[1] = None
                    return None, None

                vad1 = radsimp(1/vad)
                return _mexpand(sqrt(vad/2) +
                                sign(v[1])*sqrt(_mexpand(v[1]**2*R*vad1/2))), f
            else: #Solution requires a fourth root
                s2 = _mexpand(v[1]*R) + d
                if s2 <= 0:
                    return sqrt(nested[-1]), [0]*len(nested)
                FR, s = root(_mexpand(R), 4), sqrt(s2)
                return _mexpand(s/(sqrt(2)*FR) + v[0]*FR/(sqrt(2)*s)), f
