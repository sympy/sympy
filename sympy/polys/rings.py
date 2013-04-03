"""Sparse polynomial rings. """

from operator import add, mul

from sympy.core.expr import Expr
from sympy.core.symbol import symbols as _symbols
from sympy.core.numbers import igcd
from sympy.core.sympify import CantSympify, sympify
from sympy.core.compatibility import is_sequence
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.monomialtools import (monomial_mul, monomial_div,
    monomial_ldiv, monomial_pow, monomial_min, monomial_gcd, lex)
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.compatibility import IPolys
from sympy.polys.polyutils import expr_from_dict, _dict_reorder
from sympy.polys.polyerrors import CoercionFailed, GeneratorsError, GeneratorsNeeded
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.polynomialring import PolynomialRing
from sympy.printing.defaults import DefaultPrinting

def ring(symbols, domain, order=lex):
    """Construct new polynomial ring returning (ring, x1, ..., xn). """
    _ring = PolyRing(symbols, domain, order)
    return (_ring,) + _ring.gens

def xring(symbols, domain, order=lex):
    """Construct new polynomial ring returning (ring, (x1, ..., xn)). """
    _ring = PolyRing(symbols, domain, order)
    return (_ring, _ring.gens)

def vring(symbols, domain, order=lex):
    """Construct new polynomial ring and inject generators into global namespace. """
    from inspect import currentframe
    frame = currentframe().f_back

    try:
        _ring = PolyRing(symbols, domain, order)

        for sym, gen in zip(_ring.symbols, _ring.gens):
            frame.f_globals[sym.name] = gen
    finally:
        del frame  # break cyclic dependencies as stated in inspect docs

    return _ring

def _parse_symbols(symbols):
    if not symbols:
        raise GeneratorsNeeded("generators weren't specified")

    if isinstance(symbols, basestring):
        return _symbols(symbols, seq=True)
    elif isinstance(symbols, Expr):
        return (symbols,)
    elif is_sequence(symbols):
        if all(isinstance(s, basestring) for s in symbols):
            return _symbols(symbols)
        elif all(isinstance(s, Expr) for s in symbols):
            return symbols

    raise GeneratorsError("expected a string, Symbol or expression or a non-empty sequence of strings, Symbols or expressions")

class PolyRing(DefaultPrinting, IPolys):

    def __init__(self, symbols, domain, order):
        self.dtype = PolyElement
        self.symbols = tuple(_parse_symbols(symbols))
        self.ngens = len(self.symbols)
        self.domain = domain
        self.order = order

        self.zero_monom = (0,)*self.ngens
        self.gens = self._gens()

    def _gens(self):
        """Return a list of polynomial generators. """
        one = self.domain.one
        _gens = []
        for i in xrange(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            _gens.append(poly)
        return tuple(_gens)

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.symbols, self.domain, self.order))
        return _hash

    def __eq__(self, other):
        return isinstance(other, PolyRing) and \
               self.symbols == other.symbols and \
               self.domain == other.domain and \
               self.order == other.order

    def __ne__(self, other):
        return not self.__eq__(other)

    def clone(self, symbols=None, domain=None, order=None):
        return self.__class__(symbols or self.symbols, domain or self.domain, order or self.order)

    def monomial_basis(self, i):
        """Return the ith-basis element. """
        basis = [0]*self.ngens
        basis[i] = 1
        return tuple(basis)

    def domain_new(self, element, orig_domain=None):
        return self.domain.convert(element, orig_domain)

    def ground_new(self, coeff):
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff):
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    def ring_new(self, element):
        if isinstance(element, PolyElement):
            if self == element.ring:
                return element
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, basestring):
            raise NotImplementedError("parsing")
        elif isinstance(element, dict):
            return self.from_dict(element)
        elif isinstance(element, list):
            return self.from_terms(element)
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            return self.ground_new(element)

    __call__ = ring_new

    def _rebuild_expr(self, expr, mapping):
        domain = self.domain

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return reduce(add, map(_rebuild, expr.args))
            elif expr.is_Mul:
                return reduce(mul, map(_rebuild, expr.args))
            elif expr.is_Pow and expr.exp.is_Integer and expr.exp >= 0:
                return _rebuild(expr.base)**int(expr.exp)
            else:
                return domain.convert(expr)

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        mapping = dict(zip(self.symbols, self.gens))

        try:
            poly = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a polynomial in %s, got %s" % (self, expr))
        else:
            return self.ring_new(poly)

    @property
    def zero(self):
        return self.dtype(self)

    @property
    def one(self):
        poly = self.zero
        poly[self.zero_monom] = self.domain.one
        return poly

    def from_dict(self, d):
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in d.iteritems():
            coeff = domain_new(coeff)
            if coeff:
                poly[monom] = coeff

        return poly

    def from_terms(self, terms):
        return self.from_dict(dict(terms))

    def _drop(self, gen):
        if isinstance(gen, int):
            i = gen
            if not (0 <= i and i < self.ngens):
                raise ValueError("invalid generator index")
        else:
            if gen not in self.gens:
                raise ValueError("invalid generator")
            else:
                i = list(self.gens).index(gen)

        if self.ngens == 1:
            raise ValueError("univariate polynomial") # TODO: return ground domain
        else:
            symbols = list(self.symbols)
            del symbols[i]
            return i, self.__class__(symbols, self.domain, self.order)

    def drop(self, gen):
        return self._drop(gen)[1]

    def __getitem__(self, key):
        return self.__class__(self.symbols[key], self.domain, self.order)

    def to_ground(self):
        ground = getattr(self.domain, "dom", None) # TODO: use CompositeDomain
        if ground is not None:
            return self.__class__(self.symbols, ground, self.order)
        else:
            raise ValueError("%s is not a composite domain" % self.domain)

    def to_domain(self):
        return PolynomialRing(self)

    def to_field(self):
        from sympy.polys.fields import FracField
        return FracField(self.symbols, self.domain, self.order)

class PolyElement(DomainElement, DefaultPrinting, CantSympify, dict):
    def __init__(self, ring, init=[]):
        self.ring = ring
        dict.__init__(self, init)

    def new(self, init):
        return self.__class__(self.ring, init)

    def parent(self):
        return self.ring.to_domain()

    _hash = None

    def __hash__(self):
        # XXX: This computes a hash of a dictionary, but currently we don't
        # protect dictionary from being changed so any use site modifications
        # will make hashing go wrong. Use this feature with caution until we
        # figure out how to make a safe API without compromising speed of this
        # low-level class.
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.ring, frozenset(self.items())))
        return _hash

    def copy(self):
        """Return a copy of polynomial self.

        Polynomials are mutable; if one is interested in preserving
        a polynomial, and one plans to use inplace operations, one
        can copy the polynomial. This method makes a shallow copy.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> R, x, y = ring('x, y', ZZ)
        >>> p = (x + y)**2
        >>> p1 = p.copy()
        >>> p2 = p
        >>> p[R.zero_monom] = 3
        >>> p
        x**2 + 2*x*y + y**2 + 3
        >>> p1
        x**2 + 2*x*y + y**2
        >>> p2
        x**2 + 2*x*y + y**2 + 3

        """
        return self.new(self)

    def set_ring(self, new_ring):
        if self.ring == new_ring:
            return self
        elif self.ring.symbols != new_ring.symbols:
            terms = zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols))
            return new_ring.from_terms(terms)
        else:
            return new_ring.from_dict(self)

    def as_expr(self, *symbols):
        if symbols and len(symbols) != self.ring.ngens:
            raise ValueError("not enough symbols, expected %s got %s" % (self.ring.ngens, len(symbols)))
        else:
            symbols = self.ring.symbols

        return expr_from_dict(self.as_expr_dict(), *symbols)

    def as_expr_dict(self):
        to_sympy = self.ring.domain.to_sympy
        return dict([ (monom, to_sympy(coeff)) for monom, coeff in self.terms() ])

    def clear_denoms(self):
        domain = self.ring.domain
        denom = domain.denom

        ground_ring = domain.get_ring()
        common = ground_ring.one
        lcm = ground_ring.lcm

        for coeff in self.values():
            common = lcm(common, denom(coeff))

        poly = self.new([ (k, v*common) for k, v in self.items() ])
        return common, poly

    def strip_zero(self):
        """Eliminate monomials with zero coefficient. """
        for k, v in list(self.items()):
            if not v:
                del self[k]

    def variables(self):
        """Return the tuple of the indices of the variables in self.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y, z = ring('x, y, z', ZZ)
        >>> p = x + y**2 + x*y
        >>> p.variables()
        (0, 1)

        """
        indices = []
        for expv in self:
            for i, e in enumerate(expv):
                if e and i not in indices:
                    indices.append(i)
        indices.sort()
        return tuple(indices)

    def __eq__(p1, p2):
        """Equality test for polynomials.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> p1 = (x + y)**2 + (x - y)**2
        >>> p1 == 4*x*y
        False
        >>> p1 == 2*(x**2 + y**2)
        True

        """
        if not p2:
            return not p1
        elif isinstance(p2, PolyElement) and p1.ring == p2.ring:
            return dict.__eq__(p1, p2)
        else:
            if len(p1) > 1:
                return False
            else:
                return p1.get(p1.ring.zero_monom) == p2

    def __ne__(p1, p2):
        return not p1.__eq__(p2)

    def drop(self, gen):
        i, ring = self.ring._drop(gen)
        poly = ring.zero
        for k, v in self.iteritems():
            if k[i] == 0:
                K = list(k)
                del K[i]
                poly[tuple(K)] = v
            else:
                raise ValueError("can't drop %s" % gen)
        return poly

    def to_dense(self):
        from sympy.polys.densebasic import dmp_from_dict
        return dmp_from_dict(self, self.ring.ngens-1, self.ring.domain)

    def to_dict(self):
        return dict(self)

    @property
    def is_generator(self):
        return self in self.ring.gens

    @property
    def is_ground(self):
        return not self or (len(self) == 1 and self.ring.zero_monom in self)

    @property
    def is_monomial(self):
        return not self or (len(self) == 1 and self.LC == 1)

    @property
    def is_term(self):
        return len(self) <= 1

    @property
    def is_negative(self):
        return self.ring.domain.is_negative(self.LC)

    @property
    def is_positive(self):
        return self.ring.domain.is_positive(self.LC)

    @property
    def is_nonnegative(self):
        return self.ring.domain.is_nonnegative(self.LC)

    @property
    def is_nonpositive(self):
        return self.ring.domain.is_nonpositive(self.LC)

    def __neg__(self):
        return self*(-1)

    def __pos__(self):
        return self

    def __add__(p1, p2):
        """Add two polynomials.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> (x + y)**2 + (x - y)**2
        2*x**2 + 2*y**2

        """
        if not p2:
            return p1.copy()
        ring = p1.ring
        if isinstance(p2, PolyElement):
            if ring == p2.ring:
                if len(p2) == 1:
                    p = p1.copy()
                    [(m, c)] = p2.items()
                    nc = p.get(m, ring.domain.zero) + c
                    if nc:
                        p[m] = nc
                    else:
                        del p[m]
                    return p
                else:
                    p = ring.zero
                    for k, v in p1.iteritems():
                        if k in p2:
                            r = v + p2[k]
                            if r:
                                p[k] = r
                        else:
                            p[k] = v
                    for k, v in p2.iteritems():
                        if k not in p1:
                            p[k] = v
                    return p
            elif isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__radd__(p1)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            p = p1.copy()
            if not cp2:
                return p
            zm = ring.zero_monom
            if zm not in list(p1.keys()):
                p[zm] = cp2
            else:
                if p2 == -p[zm]:
                    del p[zm]
                else:
                    p[zm] += cp2
            return p

    def __radd__(p1, n):
        p = p1.copy()
        if not n:
            return p
        ring = p1.ring
        zm = ring.zero_monom
        if zm not in list(p1.keys()):
            p[zm] = ring.domain_new(n)
        else:
            if n == -p[zm]:
                del p[zm]
            else:
                p[zm] += n
        return p

    def __sub__(p1, p2):
        """Subtract polynomial p2 from p1.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> p1 = x + y**2
        >>> p2 = x*y + y**2
        >>> p1 - p2
        -x*y + x

        """
        if not p2:
            return p1.copy()
        ring = p1.ring
        if isinstance(p2, PolyElement):
            if ring == p2.ring:
                if len(p2) == 1:
                    p = p1.copy()
                    [(m, c)] = p2.items()
                    nc = p.get(m, ring.domain.zero) - c
                    if nc:
                        p[m] = nc
                    else:
                        del p[m]
                    return p
                else:
                    p = ring.zero
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
            elif isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rsub__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            p = p1.copy()
            zm = ring.zero_monom
            if zm not in list(p1.keys()):
                p[zm] = -p2
            else:
                if p2 == p[zm]:
                    del p[zm]
                else:
                    p[zm] -= p2
            return p

    def __rsub__(p1, n):
        """n - p1 with n convertible to the coefficient domain.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 - p
        -x - y + 4

        """
        p = p1.ring.zero
        for expv in p1:
            p[expv] = -p1[expv]
        p += n
        return p

    def __mul__(p1, p2):
        """Multiply two polynomials.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', QQ)
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p1*p2
        x**2 - y**2

        """
        ring = p1.ring
        p = ring.zero
        if not p2:
            return p
        if isinstance(p2, PolyElement):
            if ring == p2.ring:
                get = p.get
                p2it = p2.items()
                for exp1, v1 in p1.iteritems():
                    for exp2, v2 in p2it:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                p.strip_zero()
                return p
            elif isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rmul__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            for exp1, v1 in p1.iteritems():
                v = v1*p2
                if v:
                    p[exp1] = v
            return p

    def __rmul__(p1, p2):
        """p2 * p1 with p2 in the coefficient domain of p1.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 * p
        4*x + 4*y

        """
        p = p1.ring.zero
        if not isinstance(p2, PolyElement):
            if not p2:
                return p
        p2 = p.ring.domain_new(p2)
        for exp1, v1 in p1.iteritems():
            v = p2*v1
            if v:
                p[exp1] = v
        return p

    def __pow__(self, n):
        """raise polynomial to power `n`

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p**3
        x**3 + 3*x**2*y**2 + 3*x*y**4 + y**6

        """
        ring = self.ring
        n = int(n)
        if n < 0:
            if (len(self) == 1):
                p = ring.zero
                k, v = list(self.items())[0]
                kn = monomial_pow(k, n)
                p[kn] = v**n
                return p
            raise ValueError('n >= 0 is required')
        if n == 0:
            if self:
                return ring(1)
            else:
                raise ValueError
        elif len(self) == 1:
            p = ring.zero
            k, v = list(self.items())[0]
            # treat case abs(v) = 1 separately to deal with the case
            # in which n is too large to be allowed in v**n
            kn = monomial_pow(k, n)
            if v == 1:
                p[kn] = v
            elif v == -1:
                if n % 2 == 0:
                    p[kn] = -v
                else:
                    p[kn] = v
            else:
                p[kn] = v**n
            return p
        elif n == 1:
            return self.copy()
        elif n == 2:
            return self.square()
        elif n == 3:
            return self*self.square()
        elif len(self) <= 5: # TODO: use an actuall density measure
            return self._pow_multinomial(n)
        else:
            return self._pow_generic(n)

    def _pow_generic(self, n):
        p = self.ring.one

        while True:
            if n & 1:
                p = p*self
                n -= 1
                if not n:
                    break

            self = self.square()
            n = n // 2

        return p

    def _pow_multinomial(self, n):
        multinomials = multinomial_coefficients(len(self), n).items()
        zero_monom = self.ring.zero_monom
        terms = list(self.terms())
        poly = self.ring.zero

        for multinomial, multinomial_coeff in multinomials:
            product_monom = zero_monom
            product_coeff = multinomial_coeff

            for exp, (monom, coeff) in zip(multinomial, terms):
                if exp:
                    product_monom = [ a + b*exp for a, b in zip(product_monom, monom) ]
                    product_coeff *= coeff**exp

            poly[tuple(product_monom)] = product_coeff

        return poly

    def square(self):
        """square of a polynomial

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p.square()
        x**2 + 2*x*y**2 + y**4
        """
        ring = self.ring
        p = ring.zero
        get = p.get
        keys = self.keys()
        for i in range(len(keys)):
            k1 = keys[i]
            pk = self[k1]
            for j in range(i):
                k2 = keys[j]
                exp = monomial_mul(k1, k2)
                p[exp] = get(exp, 0) + pk*self[k2]
        p = p.imul_num(2)
        get = p.get
        for k, v in self.iteritems():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, 0) + v**2
        p.strip_zero()
        return p

    def __truediv__(p1, p2):
        """division by a term in the coefficient domain or
        exact division by a polynomial

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.rings import ring
        >>> _, x, y = ring('x, y', QQ)
        >>> p1 = (x**2 + x + y)*(x**2 - y**2)
        >>> p2 = x + y
        >>> p3 = p1/p2
        >>> p4 = (x**2 + x + y)*(x - y)
        >>> p3 == p4
        True

        """
        ring = p1.ring
        if isinstance(p2, PolyElement) and ring == p2.ring:
            if len(p2) == 1:
                term = list(p2.terms())[0]
                return p1.quo_term(term)
            else:
                return p1.quo(p2)
        elif not p2:
            raise ZeroDivisionError
        else:
            try:
                p2 = ring.domain.convert(p2)
            except CoercionFailed:
                return NotImplemented
            else:
                return p1.quo_ground(p2)

    __floordiv__ = __div__ = __truediv__

    def __mod__(self, other):
        return self.rem(other)

    def __divmod__(self, other):
        return self.div(other)

    def _term_div(self):
        zm = self.ring.zero_monom
        domain = self.ring.domain
        domain_quo = domain.quo

        if domain.has_Field or not domain.is_Exact:
            def term_div((a_lm, a_lc), (b_lm, b_lc)):
                if b_lm == zm: # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if monom is not None:
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return None
        else:
            def term_div((a_lm, a_lc), (b_lm, b_lc)):
                if b_lm == zm: # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if not (monom is None or a_lc % b_lc):
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return None

        return term_div

    def div(self, fv):
        """Division algorithm, see [CLO] p64.

        fv array of polynomials
           return qv, r such that
           self = sum(fv[i]*qv[i]) + r

        All polynomials are required not to be Laurent polynomials.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> _, x, y = ring('x, y', ZZ)
        >>> f = x**3
        >>> f0 = x - y**2
        >>> f1 = x - y
        >>> qv, r = f.div((f0, f1))
        >>> qv[0]
        x**2 + x*y**2 + y**4
        >>> qv[1]
        0
        >>> r
        y**6

        TODO restrict to positive exponents
        """
        ring = self.ring
        domain = ring.domain
        ret_single = False
        if isinstance(fv, PolyElement):
            ret_single = True
            fv = [fv]
        if not self:
            if ret_single:
                return ring.zero, ring.zero
            else:
                return [], ring.zero
        for f in fv:
            if f.ring != ring:
                raise ValueError('self and f must have the same ring')
        gens = ring.gens
        s = len(fv)
        qv = [ring.zero for i in range(s)]
        p = self.copy()
        r = ring.zero
        term_div = self._term_div()
        expvs = [fx.leading_expv() for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p.leading_expv()
                term = term_div((expv, p[expv]), (expvs[i], fv[i][expvs[i]]))
                if term is not None:
                    expv1, c = term
                    qv[i] = qv[i]._iadd_monom((expv1, c))
                    p = p._iadd_poly_monom(fv[i], (expv1, -c))
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                expv =  p.leading_expv()
                r = r._iadd_monom((expv, p[expv]))
                del p[expv]
        if expv == ring.zero_monom:
            r += p
        if ret_single:
            if not qv:
                return ring.zero, r
            else:
                return qv[0], r
        else:
            return qv, r

    def quo(f, G):
        return f.div(G)[0]

    def rem(f, G):
        if isinstance(G, PolyElement):
            G = [G]
        domain = f.ring.domain
        order = f.ring.order
        r = f.ring.zero
        term_div = f._term_div()
        ltf = f.LT
        f = f.copy()
        get = f.get
        while f:
            for g in G:
                tq = term_div(ltf, g.LT)
                if tq is not None:
                    m, c = tq
                    for mg, cg in g.terms():
                        m1 = monomial_mul(mg, m)
                        c1 = get(m1, 0) - c*cg
                        if not c1:
                            del f[m1]
                        else:
                            f[m1] = c1
                    if f:
                        if order is lex:
                            ltm = max(f)
                        else:
                            ltm = max(f, key=order)
                        ltf = ltm, f[ltm]

                    break
            else:
                ltm, ltc = ltf
                if ltm in r:
                    r[ltm] += ltc
                else:
                    r[ltm] = ltc
                del f[ltm]
                if f:
                    if order is lex:
                        ltm = max(f)
                    else:
                        ltm = max(f, key=order)
                    ltf = ltm, f[ltm]
        return r

    def _iadd_monom(self, mc):
        """add to self the monomial coeff*x0**i0*x1**i1*...
        unless self is a generator -- then just return the sum of the two.

        mc is a tuple, (monom, coeff), where monomial is (i0, i1, ...)

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x**4 + 2*y
        >>> m = (1, 2)
        >>> p1 = p._iadd_monom((m, 5))
        >>> p1
        x**4 + 5*x*y**2 + 2*y
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p._iadd_monom((m, 5))
        >>> p1
        5*x*y**2 + x
        >>> p1 is p
        False

        """
        if self in self.ring.gens:
            self = self.copy()
        expv, coeff = mc
        c = self.get(expv)
        if c is None:
            self[expv] = coeff
        else:
            c += coeff
            if c:
                self[expv] = c
            else:
                del self[expv]
        return self

    def _iadd_poly_monom(p1, p2, mc):
        """add to self the product of (p)*(coeff*x0**i0*x1**i1*...)
        unless self is a generator -- then just return the sum of the two.

        mc is a tuple, (monom, coeff), where monomial is (i0, i1, ...)

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y, z = ring('x, y, z', ZZ)
        >>> p1 = x**4 + 2*y
        >>> p2 = y + z
        >>> m = (1, 2, 3)
        >>> p1 = p1._iadd_poly_monom(p2, (m, 3))
        >>> p1
        x**4 + 3*x*y**3*z**3 + 3*x*y**2*z**4 + 2*y

        """
        if p1 in p1.ring.gens:
            p1 = p1.copy()
        (m, c) = mc
        get = p1.get
        zero = p1.ring.domain.zero
        for k, v in p2.iteritems():
            ka = monomial_mul(k, m)
            coeff = get(ka, zero) + v*c
            if coeff:
                p1[ka] = coeff
            else:
                del p1[ka]
        return p1

    def leading_expv(self):
        """leading monomial tuple according to the monomial ordering

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import QQ
        >>> _, x, y, z = ring('x, y, z', QQ)
        >>> p = x**4 + x**3*y + x**2*z**2 + z**7
        >>> p.leading_expv()
        (4, 0, 0)

        """
        if self:
            order = self.ring.order
            if order is lex:
                return max(self)
            else:
                return max(self, key=order)
        else:
            return None

    @property
    def LM(self):
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom
        else:
            return expv

    @property
    def LT(self):
        expv = self.leading_expv()
        if expv is None:
            return (self.ring.zero_monom, self.ring.domain.zero)
        else:
            return (expv, self._get_coeff(expv))

    @property
    def LC(self):
        return self._get_coeff(self.leading_expv())

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def coeff(self, element):
        if element == 1:
            return self._get_coeff(self.ring.zero_monom)
        elif isinstance(element, PolyElement):
            terms = list(element.terms())
            if len(terms) == 1:
                monom, coeff = terms[0]
                if coeff == self.ring.domain.one:
                    return self._get_coeff(monom)

        raise ValueError("expected a monomial, got %s" % element)

    @property
    def leading_monom(self):
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self.ring.one
        return p

    @property
    def leading_term(self):
        """Leading term according to the monomial ordering.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import QQ

        >>> _, x, y = ring('x, y', QQ)
        >>> p = (x + y)**4
        >>> p.leading_term
        x**4

        """
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self[expv]
        return p

    def coeffs(self):
        return self.itervalues()

    def monoms(self):
        return self.iterkeys()

    def terms(self):
        return self.iteritems()

    def imul_num(p, c):
        """multiply inplace the polynomial p by an element in the
        coefficient ring, provided p is not one of the generators;
        else multiply not inplace

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p1 = p.imul_num(3)
        >>> p1
        3*x + 3*y**2
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p.imul_num(3)
        >>> p1
        3*x
        >>> p1 is p
        False

        """
        if p in p.ring.gens:
            return p*c
        if not c:
            p.clear()
            return
        for exp in p:
            p[exp] *= c
        return p

    def mul_term(f, term):
        monom, coeff = term

        if not f or not coeff:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.mul_ground(coeff)

        terms = [ (monomial_mul(f_monom, monom), f_coeff*coeff) for f_monom, f_coeff in f.iteritems() ]
        return f.new(terms)

    def mul_monom(f, monom):
        terms = [ (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in f.iteritems() ]
        return f.new(terms)

    def monic(f):
        """Divides all coefficients by the leading coefficient. """
        if not f:
            return f
        else:
            return f.quo_ground(f.LC)

    def primitive(f):
        """Returns content and a primitive polynomial. """
        cont = f.content()
        return cont, f.quo_ground(cont)

    def content(f):
        """Returns GCD of polynomial's coefficients. """
        domain = f.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in f.coeffs():
            cont = gcd(cont, coeff)

        return cont

    def mul_ground(f, x):
        if not x:
            return f.ring.zero

        terms = [ (monom, coeff*x) for monom, coeff in f.terms() ]
        return f.new(terms)

    def quo_ground(f, x):
        domain = f.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not f or x == domain.one:
            return f

        if domain.has_Field or not domain.is_Exact:
            quo = domain.quo
            terms = [ (monom, quo(coeff, x)) for monom, coeff in f.terms() ]
        else:
            terms = [ (monom, coeff // x) for monom, coeff in f.terms() if not (coeff % x) ]

        return f.new(terms)

    def quo_term(f, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError
        elif not f:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.quo_ground(coeff)

        term_div = f._term_div()

        terms = [ term_div(t, term) for t in f.terms() ]
        return f.new([ t for t in terms if t is not None ])

    def trunc_ground(f, p):
        if f.ring.domain.is_ZZ:
            terms = []

            for monom, coeff in f.terms():
                coeff = coeff % p

                if coeff > p // 2:
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [ (monom, coeff % p) for monom, coeff in f.terms() ]

        poly = f.new(terms)
        poly.strip_zero()
        return poly

    def extract_ground(f, g):
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def max_norm(f):
        if not f:
            return f.ring.domain.zero
        else:
            ground_abs = f.ring.domain.abs
            return max([ ground_abs(coeff) for coeff in f.coeffs() ])

    def deflate(f, *G):
        ring = f.ring
        polys = [f] + list(G)

        J = [0]*ring.ngens

        for p in polys:
            for monom in p.monoms():
                for i, m in enumerate(monom):
                    J[i] = igcd(J[i], m)

        for i, b in enumerate(J):
            if not b:
                J[i] = 1

        J = tuple(J)

        if all(b == 1 for b in J):
            return J, polys

        H = []

        for p in polys:
            h = ring.zero

            for I, coeff in p.terms():
                N = [ i // j for i, j in zip(I, J) ]
                h[tuple(N)] = coeff

            H.append(h)

        return J, H

    def inflate(f, J):
        poly = f.ring.zero

        for I, coeff in f.terms():
            N = [ i*j for i, j in zip(I, J) ]
            poly[tuple(N)] = coeff

        return poly

    def lcm(f, g):
        domain = f.ring.domain

        if not domain.has_Field:
            fc, f = f.primitive()
            gc, g = g.primitive()
            c = domain.lcm(fc, gc)

        h = (f*g).quo(f.gcd(g))

        if not domain.has_Field:
            return h.mul_ground(c)
        else:
            return h.monic()

    def gcd(f, g):
        return f.cofactors(g)[0]

    def cofactors(f, g):
        if not f and not g:
            zero = f.ring.zero
            return zero, zero, zero
        elif not f:
            h, cff, cfg = f._gcd_zero(g)
            return h, cff, cfg
        elif not g:
            h, cfg, cff = g._gcd_zero(f)
            return h, cff, cfg
        elif len(f) == 1:
            h, cff, cfg = f._gcd_monom(g)
            return h, cff, cfg
        elif len(g) == 1:
            h, cfg, cff = g._gcd_monom(f)
            return h, cff, cfg

        J, (f, g) = f.deflate(g)
        h, cff, cfg = f._gcd(g)

        return (h.inflate(J), cff.inflate(J), cfg.inflate(J))

    def _gcd_zero(f, g):
        one, zero = f.ring.one, f.ring.zero
        if g.is_nonnegative:
            return g, zero, one
        else:
            return -g, zero, -one

    def _gcd_monom(f, g):
        ring = f.ring
        ground_gcd = ring.domain.gcd
        ground_quo = ring.domain.quo
        mf, cf = list(f.terms())[0]
        _mgcd, _cgcd = mf, cf
        for mg, cg in g.terms():
            _mgcd = monomial_gcd(_mgcd, mg)
            _cgcd = ground_gcd(_cgcd, cg)
        h = f.new([(_mgcd, _cgcd)])
        cff = f.new([(monomial_ldiv(mf, _mgcd), ground_quo(cf, _cgcd))])
        cfg = f.new([(monomial_ldiv(mg, _mgcd), ground_quo(cg, _cgcd)) for mg, cg in g.terms()])
        return h, cff, cfg

    def _gcd(f, g):
        ring = f.ring

        if ring.domain.is_QQ:
            return f._gcd_QQ(g)
        elif ring.domain.is_ZZ:
            return f._gcd_ZZ(g)
        else: # TODO: don't use dense representation (port PRS algorithms)
            return ring.dmp_inner_gcd(f, g)

    def _gcd_ZZ(f, g):
        return heugcd(f, g)

    def _gcd_QQ(f, g):
        ring = f.ring
        new_ring = ring.clone(domain=ring.domain.get_ring())

        cf, f = f.clear_denoms()
        cg, g = g.clear_denoms()

        f = f.set_ring(new_ring)
        g = g.set_ring(new_ring)

        h, cff, cfg = f._gcd_ZZ(g)

        h = h.set_ring(ring)
        c, h = h.LC, h.monic()

        cff = cff.set_ring(ring).mul_ground(ring.domain.quo(c, cf))
        cfg = cfg.set_ring(ring).mul_ground(ring.domain.quo(c, cg))

        return h, cff, cfg

    def cancel(f, g):
        """
        Cancel common factors in a rational function ``f/g``.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x,y = ring("x,y", ZZ)

        >>> (2*x**2 - 2).cancel(x**2 - 2*x + 1)
        (2*x + 2, x - 1)

        """
        ring = f.ring
        new_ring = None
        domain = ring.domain

        if domain.has_Field and domain.has_assoc_Ring:
            new_ring = ring.clone(domain=domain.get_ring())

            cq, f = f.clear_denoms()
            cp, g = g.clear_denoms()

            f = f.set_ring(new_ring)
            g = g.set_ring(new_ring)
        else:
            cp = cq = domain.one

        _, p, q = f.cofactors(g)

        if new_ring is not None:
            p = p.set_ring(ring)
            q = q.set_ring(ring)

        p_neg = p.is_negative
        q_neg = q.is_negative

        if p_neg and q_neg:
            p, q = -p, -q
        elif p_neg:
            cp, p = -cp, -p
        elif q_neg:
            cp, q = -cp, -q

        p = p.mul_ground(cp)
        q = q.mul_ground(cq)

        return p, q

    def diff(f, x):
        """Computes partial derivative in ``x``.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y = ring("x,y", ZZ)
        >>> p = x + x**2*y**3
        >>> p.diff(x)
        2*x*y**3 + 1

        """
        ring = f.ring
        i = list(ring.gens).index(x)
        m = ring.monomial_basis(i)
        g = ring.zero
        for expv, coeff in f.terms():
            if expv[i]:
                e = monomial_ldiv(expv, m)
                g[e] = coeff*expv[i]
        return g

    def evaluate(f, x, a=None):
        if isinstance(x, list) and a is None:
            (X, a), x = x[0], x[1:]
            f = f.evaluate(X, a)
            if not x:
                return f
            else:
                x = [ (Y.drop(X), a) for (Y, a) in x ]
                return f.evaluate(x)

        ring = f.ring
        i = list(ring.gens).index(x)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.terms():
                result += coeff*a**n

            return result
        else:
            poly = ring[1:].zero

            for monom, coeff in f.terms():
                n, monom = monom[i], monom[:i] + monom[i+1:]
                coeff = coeff*a**n

                if monom in poly:
                    coeff = coeff + poly[monom]

                    if coeff:
                        poly[monom] = coeff
                    else:
                        del poly[monom]
                else:
                    if coeff:
                        poly[monom] = coeff

            return poly

    def subs(f, x, a=None):
        if isinstance(x, list) and a is None:
            for X, a in x:
                f = f.subs(X, a)
            return f

        ring = f.ring
        i = list(ring.gens).index(x)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.terms():
                result += coeff*a**n

            return ring.ground_new(result)
        else:
            poly = ring.zero

            for monom, coeff in f.terms():
                n, monom = monom[i], monom[:i] + (0,) + monom[i+1:]
                coeff = coeff*a**n

                if monom in poly:
                    coeff = coeff + poly[monom]

                    if coeff:
                        poly[monom] = coeff
                    else:
                        del poly[monom]
                else:
                    if coeff:
                        poly[monom] = coeff

            return poly

    def compose(f, x, a=None):
        ring = f.ring
        poly = ring.zero
        gens_map = dict(zip(ring.gens, range(ring.ngens)))

        if a is not None:
            replacements = [(x, a)]
        else:
            if isinstance(x, list):
                replacements = list(x)
            elif isinstance(x, dict):
                replacements = sorted(x.items(), key=lambda (k, _): gens_map[k])
            else:
                raise ValueError("expected a generator, value pair a sequence of such pairs")

        for k, (x, g) in enumerate(replacements):
            replacements[k] = (gens_map[x], g)

        for monom, coeff in f.terms():
            monom = list(monom)
            subpoly = ring.one

            for i, g in replacements:
                n, monom[i] = monom[i], 0
                if n:
                    subpoly *= g**n

            subpoly = subpoly.mul_term((tuple(monom), coeff))
            poly += subpoly

        return poly
