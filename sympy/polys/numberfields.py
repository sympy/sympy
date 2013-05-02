"""Computational algebraic number field theory. """

from sympy import (
    S, C, Expr, Rational,
    Symbol, Add, Mul, sympify, Q, ask, Dummy, Tuple, expand_mul, I, pi
)

from sympy.polys.polytools import (
    Poly, PurePoly, sqf_norm, invert, factor_list, groebner, resultant, degree
)

from sympy.polys.polyclasses import (
    ANP, DMP,
)

from sympy.polys.polyerrors import (
    IsomorphismFailed,
    CoercionFailed,
    NotAlgebraic,
)

from sympy.polys.rootoftools import RootOf

from sympy.printing.lambdarepr import LambdaPrinter

from sympy.utilities import (
    numbered_symbols, variations, lambdify,
)

from sympy.simplify.simplify import _mexpand, _is_sum_surds
from sympy.ntheory import sieve
from sympy.mpmath import pslq, mp


def _choose_factor(factors, x, v, prec=200):
    """
    Return a factor having root ``v``
    It is assumed that one of the factors has root ``v``.
    """
    if isinstance(factors[0], tuple):
        factors = [xx[0] for xx in factors]
    if len(factors) == 1:
        return factors[0]
    prec1 = 10
    reps = {x: v}
    eps = 1./10**prec1
    while 1:
        candidates = []
        for f in factors:
            if abs(f.evalf(prec1, subs=reps)) < eps:
                candidates.append(f)
        if candidates:
            factors = candidates
        if len(factors) == 1:
            return factors[0]
        if prec1 > prec:
            return None
        prec1 *= 2

def _separate_sq(p):
    """
    helper function for ``_minimal_polynomial_sq``

    It selects a rational ``g`` such that the polynomial ``p``
    consists of a sum of terms whose surds squared have gcd equal to ``g``
    and a sum of terms with surds squared prime with ``g``;
    then it takes the field norm to eliminate ``sqrt(g)``

    See simplify.simplify.split_surds and polytools.sqf_norm.

    Examples
    ========

    >>> from sympy import sqrt
    >>> from sympy.abc import x
    >>> from sympy.polys.numberfields import _separate_sq
    >>> p= -x + sqrt(2) + sqrt(3) + sqrt(7)
    >>> p = _separate_sq(p); p
    -x**2 + 2*sqrt(3)*x + 2*sqrt(7)*x - 2*sqrt(21) - 8
    >>> p = _separate_sq(p); p
    -x**4 + 4*sqrt(7)*x**3 - 32*x**2 + 8*sqrt(7)*x + 20
    >>> p = _separate_sq(p); p
    -x**8 + 48*x**6 - 536*x**4 + 1728*x**2 - 400

    """
    from sympy.simplify.simplify import _split_gcd, _mexpand
    from sympy.utilities.iterables import sift
    def is_sqrt(expr):
        return expr.is_Pow and expr.exp is S.Half
    # p = c1*sqrt(q1) + ... + cn*sqrt(qn) -> a = [(c1, q1), .., (cn, qn)]
    a = []
    for y in p.args:
        if not y.is_Mul:
            if is_sqrt(y):
                a.append((S.One, y**2))
            elif y.is_Atom:
                a.append((y, S.One))
            else:
                raise NotImplementedError
            continue
        sifted = sift(y.args, is_sqrt)
        a.append((Mul(*sifted[False]), Mul(*sifted[True])**2))
    a.sort(key=lambda z: z[1])
    if a[-1][1] is S.One:
        # there are no surds
        return p
    surds = [z for y, z in a]
    for i in range(len(surds)):
        if surds[i] != 1:
            break
    g, b1, b2 = _split_gcd(*surds[i:])
    a1 = []
    a2 = []
    for y, z in a:
        if z in b1:
            a1.append(y*z**S.Half)
        else:
            a2.append(y*z**S.Half)
    p1 = Add(*a1)
    p2 = Add(*a2)
    p = _mexpand(p1**2) - _mexpand(p2**2)
    return p

def _minimal_polynomial_sq(p, n, x, prec):
    """
    Returns the minimal polynomial for the ``nth-root`` of a sum of surds
    or ``None`` if it fails.

    Parameters
    ==========

    p : sum of surds
    n : positive integer
    x : variable of the returned polynomial

    Examples
    ========

    >>> from sympy.polys.numberfields import _minimal_polynomial_sq
    >>> from sympy import sqrt
    >>> from sympy.abc import x
    >>> q = 1 + sqrt(2) + sqrt(3)
    >>> _minimal_polynomial_sq(q, 3, x, 15)
    x**12 - 4*x**9 - 4*x**6 + 16*x**3 - 8

    """
    from sympy.simplify.simplify import _is_sum_surds

    p = sympify(p)
    n = sympify(n)
    r = _is_sum_surds(p)
    if not n.is_Integer or not n > 0 or not _is_sum_surds(p):
        return None
    pn = p**Rational(1, n)
    # eliminate the square roots
    p -= x
    while 1:
        p1 = _separate_sq(p)
        if p1 is p:
            p = p1.subs({x:x**n})
            break
        else:
            p = p1

    # _separate_sq eliminates field extensions in a minimal way, so that
    # if n = 1 then `p = constant*(minimal_polynomial(p))`
    # if n > 1 it contains the minimal polynomial as a factor.
    if n == 1:
        p1 = Poly(p)
        if p.coeff(x**p1.degree()) < 0:
            p = -p
        p = p.primitive()[1]
        return p
    # by construction `p` has root `pn`
    # the minimal polynomial is the factor vanishing in x = pn
    factors = factor_list(p)[1]

    result = _choose_factor(factors, x, pn, prec)
    return result

def _minpoly_neg(ex, x, mp):
    """
    return the minimal polynomial of ``-ex``.

    Examples
    ========

    >>> from sympy import Rational
    >>> from sympy.polys.numberfields import minpoly, _minpoly_neg
    >>> from sympy.abc import x
    >>> p = 2**Rational(1, 3)
    >>> mp = minpoly(p, x)
    >>> _minpoly_neg(p, x, mp=mp)
    -x**3 - 2
    """
    if mp is None:
        mp = minpoly1(ex, x)
    return mp.subs({x: -x})

def minpoly_op_algebraic_number(ex1, ex2, x, mp1=None, mp2=None, prec=200,
        method='resultant', op=Add):
    """
    return the minimal polinomial for ``op(ex1, ex2)``

    Parameters
    ==========

    ex1, ex2 : expressions for the algebraic numbers
    x : indeterminate of the polynomials
    mp1, mp2 : minimal polynomials for ``ex1`` and ``ex2`` or None
    prec : max precision used in identifying factors
    method : use 'resultant' or 'groebner' method
    op : operation ``Add``, ``Mul`` or 'div'

    Examples
    ========

    >>> from sympy import sqrt, Mul
    >>> from sympy.polys.numberfields import minpoly_op_algebraic_number
    >>> from sympy.abc import x
    >>> p1 = sqrt(sqrt(2) + 1)
    >>> p2 = sqrt(sqrt(2) - 1)
    >>> minpoly_op_algebraic_number(p1, p2, x, method='groebner', op=Mul)
    x - 1

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Resultant
    [2] J.V. Brawley and L. Carlitz, Discrete Math., 65(2):115 (1987)
    "Irreducibles and the composed product over a finite field".
    """
    from sympy import gcd
    y = Dummy(str(x))
    if mp1 is None:
        mp1 = minpoly1(ex1, x)
    if mp2 is None:
        mp2 = minpoly1(ex2, y)
    else:
        mp2 = mp2.subs({x:y})

    if method == 'resultant':
        if op is Add:
            mp1a = mp1.subs({x:x - y})
        elif op is Mul:
            mp1a = y**degree(mp1)*mp1.subs({x:x / y})
        elif op is 'div':
            mp1a = mp1.subs({x:x * y})
        else:
            raise NotImplementedError('option not available')
        mp1a = _mexpand(mp1a)
        r = resultant(mp1a, mp2, gens=[y])
    elif method == 'groebner':
        z = Dummy('z')
        if op in [Add, Mul]:
            g = groebner([mp1, mp2, op(x, y) - z], gens=[x, y, z], order='lex')
        elif op is 'div':
            g = groebner([mp1, mp2, x - z*y], gens=[x, y, z], order='lex')
        else:
            raise NotImplementedError('option not available')
        r = g[-1].subs({z:x})
    else:
        raise NotImplementedError('option not available')

    gcd_deg = gcd(degree(mp1), degree(mp2))
    if gcd_deg == 1:
        # in this case `r` is irreducible, see [2]
        return r
    _, factors = factor_list(r)
    if op in [Add, Mul]:
        ex = op(ex1, ex2)
    elif op is 'div':
        ex = ex1 / ex2
    res = _choose_factor(factors, x, ex, prec)
    return res

def _minpoly_pow(ex, pw, x, mp=None, prec=200):
    """
    Returns ``minpoly(ex**pw, x)``

    Parameters
    ==========

    p  : algebraic number
    mp : minimal polynomial of ``p``
    pw : rational number
    x : indeterminate of the polynomial

    Examples
    ========

    >>> from sympy import sqrt
    >>> from sympy.polys.numberfields import _minpoly_pow, minpoly
    >>> from sympy.abc import x
    >>> p = sqrt(1 + sqrt(2))
    >>> _minpoly_pow(p, 2, x)
    x**2 - 2*x - 1
    >>> minpoly(p**2, x)
    x**2 - 2*x - 1
    """
    pw = sympify(pw)
    if not mp:
        mp = minpoly1(ex, x)
    if not pw.is_rational:
        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)
    if pw < 0:
        mp = expand_mul(x**degree(mp)*mp.subs(x, 1/x))
        pw = -pw
        ex = 1/ex
    y = Dummy(str(x))
    mp = mp.subs({x:y})
    n, d = pw.as_numer_denom()
    res = resultant(mp, x**d - y**n, gens=[y])
    _, factors = factor_list(res)
    res = _choose_factor(factors, x, ex**pw, prec)
    return res

def _minpoly_add(x, *a):
    if not a:
        return S.Zero
    if len(a) == 1:
        return minpoly1(a[0], x)
    mp = minpoly_op_algebraic_number(a[0], a[1], x, op=Add)
    p = a[0] + a[1]
    for px in a[2:]:
        p1 = p + px
        mp1 = minpoly_op_algebraic_number(p, px, x, mp1=mp, op=Add)
        p = p1
        mp = mp1
    return mp

def _minpoly_mul(x, *a):
    if not a:
        return S.One
    if len(a) == 1:
        return minpoly1(a[0], x)
    mp = minpoly_op_algebraic_number(a[0], a[1], x, op=Mul)
    p = a[0] + a[1]
    for px in a[2:]:
        p1 = p * px
        mp1 = minpoly_op_algebraic_number(p, px, x, mp1=mp, op=Mul)
        p = p1
        mp = mp1
    return mp

def _minpoly_sin(ex, x):
    """
    Returns the minimal polynomial of ``sin(ex)``
    see http://mathworld.wolfram.com/TrigonometryAngles.html
    """
    from sympy.functions.combinatorial.factorials import binomial
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        if c.is_rational:
            if c.p == 1:
                q = sympify(c.q)
                if q == 7:
                    return 64*x**6 - 112*x**4 + 56*x**2 - 7
                if q == 9:
                    return 64*x**6 - 96*x**4 + 36*x**2 - 3
                if q.is_prime:
                    s = 0
                    q2 = q // 2
                    for k in range(q2 + 1):
                        s += (-1)**k*binomial(q, 2*k + 1)*(1 - x**2)**(q2 - k)*x**(2*k)
                    return _mexpand(s)
                else:
                    raise NotImplementedError('case not covered')
                    # too slow
                    #ex1 = (C.exp(I*ex.args[0]) - C.exp(-I*ex.args[0]))/(2*I)
                    #return minpoly1(ex1, x)
            else:
                ex1 = (C.exp(I*ex.args[0]) - C.exp(-I*ex.args[0]))/(2*I)
                return minpoly1(ex1, x)
        else:
            raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)
    raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

def _minpoly_cos(ex, x):
    """
    Returns the minimal polynomial of ``cos(ex)``
    see http://mathworld.wolfram.com/TrigonometryAngles.html
    """
    from sympy import sqrt
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        if c.is_rational:
            if c.p == 1:
                q = sympify(c.q)
                if q.is_prime:
                    s = _minpoly_sin(ex, x)
                    factors = [_mexpand(s.subs({x:sqrt((1 + x)/2)})),
                            _mexpand(s.subs({x:sqrt((1 - x)/2)}))]

                    s = _choose_factor(factors, x, ex)
                    return s
                else:
                    raise NotImplementedError('case not covered')
                    # too slow
                    #ex1 = (C.exp(I*ex.args[0]) + C.exp(-I*ex.args[0]))/2
                    #return minpoly1(ex1, x)
            else:
                ex1 = (C.exp(I*ex.args[0]) + C.exp(-I*ex.args[0]))/2
                return minpoly1(ex1, x)
        else:
            raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)
    raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

def _minpoly_exp(ex, x):
    """
    Returns the minimal polynomial of ``exp(ex)``
    """
    c, a = ex.args[0].as_coeff_Mul()
    p = sympify(c.p)
    q = sympify(c.q)
    if a == I*pi:
        if c.is_rational:
            if c.p == 1 or c.p == -1:
                if q == 3:
                    return x**2 - x + 1
                if q == 4:
                    return x**4 + 1
                if q == 6:
                    return x**4 - x**2 + 1
                if q == 8:
                    return x**8 + 1
                if q == 9:
                    return x**6 - x**3 + 1
                if q == 10:
                    return x**8 - x**6 + x**4 - x**2 + 1
                if q.is_prime:
                    s = 0
                    for i in range(q):
                        s += (-x)**i
                    return s
                else:
                    raise NotImplementedError('case not covered')
            else:
                ex1 = C.exp(I*pi/q)
                mp = _minpoly_pow(ex1, p, x)
                return mp
        else:
            raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)
    raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

def _minpoly_rootof(ex, x):
    # TODO is `ex` already irreducible?
    p = ex.expr
    p = p.subs({ex.poly.gens[0]:x})
    _, factors = factor_list(p, x)
    result = _choose_factor(factors, x, ex)
    return result


def minpoly1(ex, x):
    """
    Computes the minimal polynomial of an algebraic number
    using operations on minimal polynomials

    Examples
    ========

    >>> from sympy import minimal_polynomial, sqrt, Rational
    >>> from sympy.abc import x
    >>> minimal_polynomial(sqrt(2) + 3*Rational(1, 3), x, compose=True)
    x**2 - 2*x - 1
    """
    if ex.is_Rational:
        return ex.q*x - ex.p
    if ex is I:
        return x**2 + 1

    if _is_sum_surds(ex):
        # eliminate the square roots
        ex -= x
        while 1:
            ex1 = _separate_sq(ex)
            if ex1 is ex:
                return ex
            else:
                ex = ex1

    if ex.is_Add:
        res = _minpoly_add(x, *ex.args)
    elif ex.is_Mul:
        res = _minpoly_mul(x, *ex.args)
    elif ex.is_Pow:
        res = _minpoly_pow(ex.base, ex.exp, x)
    elif ex.__class__ is C.sin:
        res = _minpoly_sin(ex, x)
    elif ex.__class__ is C.cos:
        res = _minpoly_cos(ex, x)
    elif ex.__class__ is C.exp:
        res = _minpoly_exp(ex, x)
    elif ex.__class__ is RootOf:
        res = _minpoly_rootof(ex, x)
    else:
        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)
    return res

def has_algebraic_number(p):
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
        return False
    elif p.is_Atom:
        return False
    elif p.is_Pow:
        return all(has_algebraic_number(x) for x in (p.base, p.exp))
    elif p.is_Add or p.is_Mul:
        return all(has_algebraic_number(x) for x in p.args)
    elif p.__class__ is RootOf:
        return False
    else:
        return False

def minimal_polynomial(ex, x=None, **args):
    """
    Computes the minimal polynomial of an algebraic number.

    Examples
    ========

    >>> from sympy import minimal_polynomial, sqrt
    >>> from sympy.abc import x

    >>> minimal_polynomial(sqrt(2), x)
    x**2 - 2
    >>> minimal_polynomial(sqrt(2) + sqrt(3), x)
    x**4 - 10*x**2 + 1

    """
    from sympy.polys.polytools import degree
    from sympy.core.function import expand_multinomial
    from sympy.core.basic import preorder_traversal

    compose = args.get('compose', True)
    polys = args.get('polys', False)
    ex = sympify(ex)
    for expr in preorder_traversal(ex):
        if expr.is_AlgebraicNumber:
            compose = False
            break

    if ex.is_AlgebraicNumber or polys:
        compose = False

    if x is not None:
        x, cls = sympify(x), Poly
    else:
        x, cls = Dummy('x'), PurePoly

    if compose:
        result = minpoly1(ex, x)
        c = result.coeff(x**degree(result))
        if c < 0:
            result = expand_mul(-result)
            c = -c
        if c != 1:
            result = result.primitive()[1]
        return result
    generator = numbered_symbols('a', cls=Dummy)
    mapping, symbols, replace = {}, {}, []

    def update_mapping(ex, exp, base=None):
        a = generator.next()
        symbols[ex] = a

        if base is not None:
            mapping[ex] = a**exp + base
        else:
            mapping[ex] = exp.as_expr(a)

        return a

    def bottom_up_scan(ex):
        if ex.is_Atom:
            if ex is S.ImaginaryUnit:
                if ex not in mapping:
                    return update_mapping(ex, 2, 1)
                else:
                    return symbols[ex]
            elif ex.is_Rational:
                return ex
        elif ex.is_Add:
            return Add(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Mul:
            return Mul(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Pow:
            if ex.exp.is_Rational:
                if ex.exp < 0 and ex.base.is_Add:
                    coeff, terms = ex.base.as_coeff_add()
                    elt, _ = primitive_element(terms, polys=True)

                    alg = ex.base - coeff

                    # XXX: turn this into eval()
                    inverse = invert(elt.gen + coeff, elt).as_expr()
                    base = inverse.subs(elt.gen, alg).expand()

                    if ex.exp == -1:
                        return bottom_up_scan(base)
                    else:
                        ex = base**(-ex.exp)
                if not ex.exp.is_Integer:
                    base, exp = (
                        ex.base**ex.exp.p).expand(), Rational(1, ex.exp.q)
                else:
                    base, exp = ex.base, ex.exp
                base = bottom_up_scan(base)
                expr = base**exp

                if expr not in mapping:
                    return update_mapping(expr, 1/exp, -base)
                else:
                    return symbols[expr]
        elif ex.is_AlgebraicNumber:
            if ex.root not in mapping:
                return update_mapping(ex.root, ex.minpoly)
            else:
                return symbols[ex.root]

        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

    def simpler_inverse(ex):
        """
        Returns True if it is more likely that the minimal polynomial
        algorithm works better with the inverse
        """
        if ex.is_Pow:
            if (1/ex.exp).is_integer and ex.exp < 0:
                if ex.base.is_Add:
                    return True
        if ex.is_Mul:
            hit = True
            a = []
            for p in ex.args:
                if p.is_Add:
                    return False
                if p.is_Pow:
                    if p.base.is_Add and p.exp > 0:
                        return False

            if hit:
                return True
        return False

    prec = args.pop('prec', 10)

    inverted = False
    ex = expand_multinomial(ex)
    if ex.is_AlgebraicNumber:
        if not polys:
            return ex.minpoly.as_expr(x)
        else:
            return ex.minpoly.replace(x)
    elif ex.is_Rational:
        result = ex.q*x - ex.p
    else:
        inverted = simpler_inverse(ex)
        if inverted:
            ex = ex**-1
        res = None
        if ex.is_Pow and (1/ex.exp).is_Integer:
            n = 1/ex.exp
            res = _minimal_polynomial_sq(ex.base, n, x, prec)

        elif _is_sum_surds(ex):
            res = _minimal_polynomial_sq(ex, S.One, x, prec)

        if res is not None:
            result = res

        if res is None:
            bus = bottom_up_scan(ex)
            F = [x - bus] + mapping.values()
            G = groebner(F, symbols.values() + [x], order='lex')

            _, factors = factor_list(G[-1])
            # by construction G[-1] has root `ex`
            result = _choose_factor(factors, x, ex, prec)
            if result is None:
                raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % ex)
    if inverted:
        result = expand_mul(x**degree(result)*result.subs(x, 1/x))
        if result.coeff(x**degree(result)) < 0:
            result = expand_mul(-result)
    if polys:
        return cls(result, x, field=True)
    else:
        return result

minpoly = minimal_polynomial


def _coeffs_generator(n):
    """Generate coefficients for `primitive_element()`. """
    for coeffs in variations([1, -1], n, repetition=True):
        yield list(coeffs)


def primitive_element(extension, x=None, **args):
    """Construct a common number field for all extensions. """
    if not extension:
        raise ValueError("can't compute primitive element for empty extension")

    if x is not None:
        x, cls = sympify(x), Poly
    else:
        x, cls = Dummy('x'), PurePoly
    if not args.get('ex', False):
        extension = [ AlgebraicNumber(ext, gen=x) for ext in extension ]

        g, coeffs = extension[0].minpoly.replace(x), [1]

        for ext in extension[1:]:
            s, _, g = sqf_norm(g, x, extension=ext)
            coeffs = [ s*c for c in coeffs ] + [1]

        if not args.get('polys', False):
            return g.as_expr(), coeffs
        else:
            return cls(g), coeffs

    generator = numbered_symbols('y', cls=Dummy)

    F, Y = [], []

    for ext in extension:
        y = generator.next()

        if ext.is_Poly:
            if ext.is_univariate:
                f = ext.as_expr(y)
            else:
                raise ValueError("expected minimal polynomial, got %s" % ext)
        else:
            f = minpoly(ext, y)

        F.append(f)
        Y.append(y)

    coeffs_generator = args.get('coeffs', _coeffs_generator)

    for coeffs in coeffs_generator(len(Y)):
        f = x - sum([ c*y for c, y in zip(coeffs, Y)])
        G = groebner(F + [f], Y + [x], order='lex', field=True)

        H, g = G[:-1], cls(G[-1], x, domain='QQ')

        for i, (h, y) in enumerate(zip(H, Y)):
            try:
                H[i] = Poly(y - h, x,
                            domain='QQ').all_coeffs()  # XXX: composite=False
            except CoercionFailed:  # pragma: no cover
                break  # G is not a triangular set
        else:
            break
    else:  # pragma: no cover
        raise RuntimeError("run out of coefficient configurations")

    _, g = g.clear_denoms()

    if not args.get('polys', False):
        return g.as_expr(), coeffs, H
    else:
        return g, coeffs, H


def is_isomorphism_possible(a, b):
    """Returns `True` if there is a chance for isomorphism. """
    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if m % n != 0:
        return False

    if n == m:
        return True

    da = a.minpoly.discriminant()
    db = b.minpoly.discriminant()

    i, k, half = 1, m//n, db//2

    while True:
        p = sieve[i]
        P = p**k

        if P > half:
            break

        if ((da % p) % 2) and not (db % P):
            return False

        i += 1

    return True


def field_isomorphism_pslq(a, b):
    """Construct field isomorphism using PSLQ algorithm. """
    if not a.root.is_real or not b.root.is_real:
        raise NotImplementedError("PSLQ doesn't support complex coefficients")

    f = a.minpoly
    g = b.minpoly.replace(f.gen)

    n, m, prev = 100, b.minpoly.degree(), None

    for i in xrange(1, 5):
        A = a.root.evalf(n)
        B = b.root.evalf(n)

        basis = [1, B] + [ B**i for i in xrange(2, m) ] + [A]

        dps, mp.dps = mp.dps, n
        coeffs = pslq(basis, maxcoeff=int(1e10), maxsteps=1000)
        mp.dps = dps

        if coeffs is None:
            break

        if coeffs != prev:
            prev = coeffs
        else:
            break

        coeffs = [S(c)/coeffs[-1] for c in coeffs[:-1]]

        while not coeffs[-1]:
            coeffs.pop()

        coeffs = list(reversed(coeffs))
        h = Poly(coeffs, f.gen, domain='QQ')

        if f.compose(h).rem(g).is_zero:
            d, approx = len(coeffs) - 1, 0

            for i, coeff in enumerate(coeffs):
                approx += coeff*B**(d - i)

            if A*approx < 0:
                return [ -c for c in coeffs ]
            else:
                return coeffs
        elif f.compose(-h).rem(g).is_zero:
            return [ -c for c in coeffs ]
        else:
            n *= 2

    return None


def field_isomorphism_factor(a, b):
    """Construct field isomorphism via factorization. """
    _, factors = factor_list(a.minpoly, extension=b)

    for f, _ in factors:
        if f.degree() == 1:
            coeffs = f.rep.TC().to_sympy_list()
            d, terms = len(coeffs) - 1, []

            for i, coeff in enumerate(coeffs):
                terms.append(coeff*b.root**(d - i))

            root = Add(*terms)

            if (a.root - root).evalf(chop=True) == 0:
                return coeffs

            if (a.root + root).evalf(chop=True) == 0:
                return [ -c for c in coeffs ]
    else:
        return None


def field_isomorphism(a, b, **args):
    """Construct an isomorphism between two number fields. """
    a, b = sympify(a), sympify(b)

    if not a.is_AlgebraicNumber:
        a = AlgebraicNumber(a)

    if not b.is_AlgebraicNumber:
        b = AlgebraicNumber(b)

    if a == b:
        return a.coeffs()

    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if n == 1:
        return [a.root]

    if m % n != 0:
        return None

    if args.get('fast', True):
        try:
            result = field_isomorphism_pslq(a, b)

            if result is not None:
                return result
        except NotImplementedError:
            pass

    return field_isomorphism_factor(a, b)


def to_number_field(extension, theta=None, **args):
    """Express `extension` in the field generated by `theta`. """
    gen = args.get('gen')

    if hasattr(extension, '__iter__'):
        extension = list(extension)
    else:
        extension = [extension]

    if len(extension) == 1 and type(extension[0]) is tuple:
        return AlgebraicNumber(extension[0])

    minpoly, coeffs = primitive_element(extension, gen, polys=True)
    root = sum([ coeff*ext for coeff, ext in zip(coeffs, extension) ])

    if theta is None:
        return AlgebraicNumber((minpoly, root))
    else:
        theta = sympify(theta)

        if not theta.is_AlgebraicNumber:
            theta = AlgebraicNumber(theta, gen=gen)

        coeffs = field_isomorphism(root, theta)

        if coeffs is not None:
            return AlgebraicNumber(theta, coeffs)
        else:
            raise IsomorphismFailed(
                "%s is not in a subfield of %s" % (root, theta.root))


class AlgebraicNumber(Expr):
    """Class for representing algebraic numbers in SymPy. """

    __slots__ = ['rep', 'root', 'alias', 'minpoly']

    is_AlgebraicNumber = True

    def __new__(cls, expr, coeffs=Tuple(), alias=None, **args):
        """Construct a new algebraic number. """
        expr = sympify(expr)

        if isinstance(expr, (tuple, Tuple)):
            minpoly, root = expr

            if not minpoly.is_Poly:
                minpoly = Poly(minpoly)
        elif expr.is_AlgebraicNumber:
            minpoly, root = expr.minpoly, expr.root
        else:
            minpoly, root = minimal_polynomial(
                expr, args.get('gen'), polys=True), expr

        dom = minpoly.get_domain()

        if coeffs != Tuple():
            if not isinstance(coeffs, ANP):
                rep = DMP.from_sympy_list(sympify(coeffs), 0, dom)
                scoeffs = Tuple(*coeffs)
            else:
                rep = DMP.from_list(coeffs.to_list(), 0, dom)
                scoeffs = Tuple(*coeffs.to_list())

            if rep.degree() >= minpoly.degree():
                rep = rep.rem(minpoly.rep)

            sargs = (root, scoeffs)

        else:
            rep = DMP.from_list([1, 0], 0, dom)

            if ask(Q.negative(root)):
                rep = -rep

            sargs = (root, coeffs)

        if alias is not None:
            if not isinstance(alias, Symbol):
                alias = Symbol(alias)
            sargs = sargs + (alias,)

        obj = Expr.__new__(cls, *sargs)

        obj.rep = rep
        obj.root = root
        obj.alias = alias
        obj.minpoly = minpoly

        return obj

    def __hash__(self):
        return super(AlgebraicNumber, self).__hash__()

    def _eval_evalf(self, prec):
        return self.as_expr()._evalf(prec)

    @property
    def is_aliased(self):
        """Returns ``True`` if ``alias`` was set. """
        return self.alias is not None

    def as_poly(self, x=None):
        """Create a Poly instance from ``self``. """
        if x is not None:
            return Poly.new(self.rep, x)
        else:
            if self.alias is not None:
                return Poly.new(self.rep, self.alias)
            else:
                return PurePoly.new(self.rep, Dummy('x'))

    def as_expr(self, x=None):
        """Create a Basic expression from ``self``. """
        return self.as_poly(x or self.root).as_expr().expand()

    def coeffs(self):
        """Returns all SymPy coefficients of an algebraic number. """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.all_coeffs() ]

    def native_coeffs(self):
        """Returns all native coefficients of an algebraic number. """
        return self.rep.all_coeffs()

    def to_algebraic_integer(self):
        """Convert ``self`` to an algebraic integer. """
        f = self.minpoly

        if f.LC() == 1:
            return self

        coeff = f.LC()**(f.degree() - 1)
        poly = f.compose(Poly(f.gen/f.LC()))

        minpoly = poly*coeff
        root = f.LC()*self.root

        return AlgebraicNumber((minpoly, root), self.coeffs())


class IntervalPrinter(LambdaPrinter):
    """Use ``lambda`` printer but print numbers as ``mpi`` intervals. """

    def _print_Integer(self, expr):
        return "mpi('%s')" % super(IntervalPrinter, self)._print_Integer(expr)

    def _print_Rational(self, expr):
        return "mpi('%s')" % super(IntervalPrinter, self)._print_Rational(expr)

    def _print_Pow(self, expr):
        return super(IntervalPrinter, self)._print_Pow(expr, rational=True)


def isolate(alg, eps=None, fast=False):
    """Give a rational isolating interval for an algebraic number. """
    alg = sympify(alg)

    if alg.is_Rational:
        return (alg, alg)
    elif not ask(Q.real(alg)):
        raise NotImplementedError(
            "complex algebraic numbers are not supported")

    func = lambdify((), alg, modules="mpmath", printer=IntervalPrinter())

    poly = minpoly(alg, polys=True)
    intervals = poly.intervals(sqf=True)

    dps, done = mp.dps, False

    try:
        while not done:
            alg = func()

            for a, b in intervals:
                if a <= alg.a and alg.b <= b:
                    done = True
                    break
            else:
                mp.dps *= 2
    finally:
        mp.dps = dps

    if eps is not None:
        a, b = poly.refine_root(a, b, eps=eps, fast=fast)

    return (a, b)
