"""Computational algebraic number field theory. """

from sympy import (
    S, Expr, I, Integer, Rational, Real,
    Symbol, Add, Mul, sympify, Q, ask,
)

from sympy.polys.polytools import (
    Poly, sqf_norm, invert, factor_list, groebner,
)

from sympy.polys.polyutils import (
    basic_from_dict,
)

from sympy.polys.polyclasses import (
    ANP, DMP,
)

from sympy.polys.polyerrors import (
    IsomorphismFailed,
    CoercionFailed,
    NotAlgebraic,
)

from sympy.utilities import (
    any, all, numbered_symbols, variations, lambdify,
)

from sympy.ntheory import sieve
from sympy.mpmath import pslq, mp

def minimal_polynomial(ex, x=None, **args):
    """Computes the minimal polynomial of an algebraic number. """
    generator = numbered_symbols('a', dummy=True)
    mapping, symbols, replace = {}, {}, []

    ex = sympify(ex)

    if x is not None:
        x = sympify(x)
    else:
        x = S.Pure

    def update_mapping(ex, exp, base=None):
        a = generator.next()
        symbols[ex] = a

        if base is not None:
            mapping[ex] = a**exp + base
        else:
            mapping[ex] = exp.as_basic(a)

        return a

    def bottom_up_scan(ex):
        if ex.is_Atom:
            if ex is S.ImaginaryUnit:
                if ex not in mapping:
                    return update_mapping(ex, 2, 1)
                else:
                    return symbols[ex]
            elif ex.is_Rational and ex.q != 0:
                return ex
        elif ex.is_Add:
            return Add(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Mul:
            return Mul(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Pow:
            if ex.exp.is_Rational:
                if ex.exp < 0 and ex.base.is_Add:
                    coeff, terms = ex.base.as_coeff_factors()
                    elt, _ = primitive_element(terms, polys=True)

                    alg = ex.base - coeff

                    # XXX: turn this into eval()
                    inverse = invert(elt.gen + coeff, elt).as_basic()
                    base = inverse.subs(elt.gen, alg).expand()

                    if ex.exp == -1:
                        return bottom_up_scan(base)
                    else:
                        ex = base**(-ex.exp)

                if not ex.exp.is_Integer:
                    base, exp = (ex.base**ex.exp.p).expand(), Rational(1, ex.exp.q)
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

    polys = args.get('polys', False)

    if ex.is_AlgebraicNumber:
        if not polys:
            return ex.minpoly.as_basic(x)
        else:
            return ex.minpoly.replace(x)
    elif ex.is_Rational and ex.q != 0:
        result = ex.q*x - ex.p
    else:
        F = [x - bottom_up_scan(ex)] + mapping.values()
        G = groebner(F, symbols.values() + [x], order='lex')

        _, factors = factor_list(G[-1])

        if len(factors) == 1:
            ((result, _),) = factors
        else:
            for result, _ in factors:
                if result.subs(x, ex).evalf(chop=True) == 0:
                    break
            else: # pragma: no cover
                raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % ex)

    if polys:
        return Poly(result, x, field=True)
    else:
        return result

minpoly = minimal_polynomial

def _coeffs_generator(n):
    """Generate coefficients for `primitive_element()`. """
    for coeffs in variations([1,-1], n, repetition=True):
        yield coeffs

def primitive_element(extension, x=None, **args):
    """Construct a common number field for all extensions. """
    if not extension:
        raise ValueError("can't compute primitive element for empty extension")

    if x is not None:
        x = sympify(x)
    else:
        x = S.Pure

    if not args.get('ex', False):
        extension = [ AlgebraicNumber(ext, gen=x) for ext in extension ]

        g, coeffs = extension[0].minpoly, [1]

        for ext in extension[1:]:
            s, _, g = sqf_norm(g, x, extension=ext)
            coeffs = [ s*c for c in coeffs ] + [1]

        if not args.get('polys', False):
            return g.as_basic(), coeffs
        else:
            return g, coeffs

    generator = numbered_symbols('y', dummy=True)

    F, Y = [], []

    for ext in extension:
        y = generator.next()

        if ext.is_Poly:
            if ext.is_univariate:
                f = ext.as_basic(y)
            else:
                raise ValueError("expected minimal polynomial, got %s" % ext)
        else:
            f = minpoly(ext, y)

        F.append(f)
        Y.append(y)

    coeffs_generator = args.get('coeffs', _coeffs_generator)

    for coeffs in coeffs_generator(len(Y)):
        f = x - sum([ c*y for c, y in zip(coeffs, Y)])
        G = groebner(F + [f], Y + [x], order='lex')

        H, g = G[:-1], Poly(G[-1], x, domain='QQ')

        for i, (h, y) in enumerate(zip(H, Y)):
            try:
                H[i] = Poly(y - h, x, domain='QQ').all_coeffs()
            except CoercionFailed: # pragma: no cover
                break # G is not a triangular set
        else:
            break
    else: # pragma: no cover
        raise RuntimeError("run out of coefficient configurations")

    _, g = g.clear_denoms()

    if not args.get('polys', False):
        return g.as_basic(), coeffs, H
    else:
        return g, coeffs, H

primelt = primitive_element

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
            d, approx = len(coeffs)-1, 0

            for i, coeff in enumerate(coeffs):
                approx += coeff*B**(d-i)

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
            d, terms = len(coeffs)-1, []

            for i, coeff in enumerate(coeffs):
                terms.append(coeff*b.root**(d-i))

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
            raise IsomorphismFailed("%s is not in a subfield of %s" % (root, theta.root))

class AlgebraicNumber(Expr):
    """Class for representing algebraic numbers in SymPy. """

    __slots__ = ['rep', 'root', 'alias', 'minpoly']

    is_AlgebraicNumber = True

    def __new__(cls, expr, coeffs=None, **args):
        """Construct a new algebraic number. """
        expr = sympify(expr)

        if type(expr) is tuple:
            minpoly, root = expr

            if not minpoly.is_Poly:
                minpoly = Poly(minpoly)
        elif expr.is_AlgebraicNumber:
            minpoly, root = expr.minpoly, expr.root
        else:
            minpoly, root = minimal_polynomial(expr, args.get('gen'), polys=True), expr

        dom = minpoly.get_domain()

        if coeffs is not None:
            if not isinstance(coeffs, ANP):
                rep = DMP.from_sympy_list(sympify(coeffs), 0, dom)
            else:
                rep = DMP.from_list(coeffs.to_list(), 0, dom)

            if rep.degree() >= minpoly.degree():
                rep = rep.rem(minpoly.rep)
        else:
            rep = DMP.from_list([1, 0], 0, dom)

            if ask(root, Q.negative):
                rep = -rep

        alias = args.get('alias')

        if alias is not None:
            if not isinstance(alias, Symbol):
                alias = Symbol(alias)

        obj = Expr.__new__(cls)

        obj.rep = rep
        obj.root = root
        obj.alias = alias
        obj.minpoly = minpoly

        return obj

    def __eq__(a, b):
        if not b.is_AlgebraicNumber:
            try:
                b = to_number_field(b, a)
            except (NotAlgebraic, IsomorphismFailed):
                return False

        return a.rep == b.rep and \
            a.minpoly.all_coeffs() == b.minpoly.all_coeffs()

    def __ne__(a, b):
        if not b.is_AlgebraicNumber:
            try:
                b = to_number_field(b, a)
            except (NotAlgebraic, IsomorphismFailed):
                return True

        return a.rep != b.rep or \
            a.minpoly.all_coeffs() != b.minpoly.all_coeffs()

    def __hash__(self):
        return super(AlgebraicNumber, self).__hash__()

    def _eval_evalf(self, prec):
        return self.as_basic()._evalf(prec)

    @property
    def is_aliased(self):
        """Returns `True` if `alias` was set. """
        return self.alias is not None

    def as_poly(self, x=None):
        """Create a Poly instance from `self`. """
        if x is not None:
            return Poly.new(self.rep, x)
        else:
            if self.alias is not None:
                return Poly.new(self.rep, self.alias)
            else:
                return Poly.new(self.rep, Symbol('x', dummy=True))

    def as_basic(self, x=None):
        """Create a Basic expression from `self`. """
        return self.as_poly(x or self.root).as_basic().expand()

    def coeffs(self):
        """Returns all SymPy coefficients of an algebraic number. """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.all_coeffs() ]

    def native_coeffs(self):
        """Returns all native coefficients of an algebraic number. """
        return self.rep.all_coeffs()

    def to_algebraic_integer(self):
        """Convert `self` to an algebraic integer. """
        f = self.minpoly

        if f.LC() == 1:
            return self

        coeff = f.LC()**(f.degree()-1)
        poly  = f.compose(Poly(f.gen/f.LC()))

        minpoly = poly*coeff
        root    = f.LC()*self.root

        return AlgebraicNumber((minpoly, root), self.coeffs())

def isolate(alg, eps=None, fast=False):
    """Give a rational isolating interval for an algebraic number. """
    alg = sympify(alg)

    if alg.is_Rational:
        return (alg, alg)
    elif not ask(alg, Q.real):
        raise NotImplementedError("complex algebraic numbers are not supported")

    from sympy.printing.lambdarepr import LambdaPrinter

    class IntervalPrinter(LambdaPrinter):
        """Use ``lambda`` printer but print numbers as ``mpi`` intervals. """

        def _print_Integer(self, expr):
            return "mpi('%s')" % super(IntervalPrinter, self)._print_Integer(expr)

        def _print_Rational(self, expr):
            return "mpi('%s')" % super(IntervalPrinter, self)._print_Rational(expr)

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

