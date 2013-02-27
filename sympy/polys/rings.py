"""Sparse polynomial rings. """

from copy import copy

from sympy.core import Symbol, symbols
from sympy.core.numbers import igcd
from sympy.core.sympify import CantSympify
from sympy.core.compatibility import is_sequence
from sympy.polys.monomialtools import monomial_mul, monomial_ldiv, monomial_pow, monomial_min, lex, term_div
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.compatibility import IPolys

def ring(sgens, domain, order=lex):
    """Construct new polynomial ring returning (ring, x1, ..., xn). """
    _ring = PolyRing(sgens, domain, order)
    return (_ring,) + _ring.gens

def xring(sgens, domain, order=lex):
    """Construct new polynomial ring returning (ring, (x1, ..., xn)). """
    _ring = PolyRing(sgens, domain, order)
    return (_ring, _ring.gens)

class PolyRing(IPolys):
    def __init__(self, sgens, domain, order):
        if not is_sequence(sgens, include=(basestring, Symbol)) or not sgens:
            raise ValueError('expecting a string, Symbol or an ordered iterable')

        if isinstance(sgens, basestring):
            sgens = [s.name for s in symbols(sgens, seq=True)]
        elif isinstance(sgens, Symbol):
            sgens = sgens.name
        elif isinstance(sgens[0], Symbol):
            sgens = [s.name for s in sgens]

        self.sgens = tuple(sgens)
        self.ngens = len(sgens)
        self.domain = domain
        self.order = order

        self.gens_dict = dict(zip(self.sgens, xrange(self.ngens)))
        self.zero_monom = (0,)*self.ngens

        self.gens = self._gens()

    def _gens(self):
        """Return a list of polynomial generators. """
        one = self.domain.one
        _gens = []
        for i in xrange(self.ngens):
            expv = self.monomial_basis(i)
            poly = PolyElement(self)
            poly[expv] = one
            _gens.append(poly)
        return tuple(_gens)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "Polynomial ring in %s over %s with %s order" % (", ".join(self.sgens), self.domain, self.order)

    def __eq__(self, other):
        return self.ngens == other.ngens and \
               self.domain == other.domain and \
               self.order == other.order

    def __ne__(self, other):
        return not self.__eq__(other)

    def clone(self, sgens=None, domain=None, order=None):
        return self.__class__(sgens or self.sgens, domain or self.domain, order or self.order)

    def monomial_basis(self, i):
        """Return the ith-basis element. """
        basis = [0]*self.ngens
        basis[i] = 1
        return tuple(basis)

    def domain_new(self, element):
        return self.domain.convert(element)

    def ground_new(self, coeff):
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff):
        element = self.domain_new(coeff)
        poly = PolyElement(self)
        if element:
            poly[monom] = coeff
        poly.strip_zero()
        return poly

    def ring_new(self, element):
        if isinstance(element, PolyElement):
            if self == element.ring:
                return element
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, basestring):
            raise NotImplementedError("parsing")
        else:
            return self.ground_new(element)

    __call__ = ring_new

    @property
    def zero(self):
        return PolyElement(self)

    @property
    def one(self):
        poly = PolyElement(self)
        poly[self.zero_monom] = self.domain.one
        return poly

    def from_dict(self, d):
        poly = PolyElement(self)
        for k, v in d.iteritems():
            poly[k] = self.domain_new(v)
        poly.strip_zero()
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
            sgens = list(self.sgens)
            del sgens[i]
            return i, self.__class__(sgens, self.domain, self.order)

    def drop(self, gen):
        return self._drop(gen)[1]

    def __getitem__(self, key):
        return self.__class__(self.sgens[key], self.domain, self.order)

    def to_ground(self):
        ground = getattr(self.domain, "dom", None) # TODO: use CompositeDomain
        if ground is not None:
            return self.__class__(self.sgens, ground, self.order)
        else:
            raise ValueError("%s is not a composite domain" % self.domain)

class PolyElement(dict, CantSympify):
    def __init__(self, ring, init=[]):
        self.ring = ring
        dict.__init__(self, init)

    @property
    def freeze(self):
        return tuple(self.items())

    def copy(self):
        """Return a copy of polynomial self.

        Polynomials are mutable; if one is interested in preserving
        a polynomial, and one plans to use inplace operations, one
        can copy the polynomial

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
        return copy(self)

    def set_ring(self, new):
        if self.ring.ngens != new.ngens:
            raise NotImplementedError
        else:
            return PolyElement(new, [ (k, new.domain_new(v)) for k, v in self.items() ])

    def clear_denoms(self):
        domain = self.ring.domain
        denom = domain.denom

        ground_ring = domain.get_ring()
        common = ground_ring.one
        lcm = ground_ring.lcm

        for coeff in self.values():
            common = lcm(common, denom(coeff))

        poly = PolyElement(self.ring, [ (k, v*common) for k, v in self.items() ])
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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.str()

    def str(self, detailed=True):
        if detailed:
            from sympy.printing import sstr
        else:
            sstr = str
        if not self:
            return sstr(self.ring.domain.zero)
        ring = self.ring
        sgens = ring.sgens
        ngens = ring.ngens
        zm = ring.zero_monom
        sexpvs = []
        expvs = list(self.keys())
        expvs.sort(key=ring.order, reverse=True)
        for expv in expvs:
            coeff = self[expv]
            if ring.domain.is_positive(coeff):
                sexpvs.append(' + ')
            else:
                sexpvs.append(' - ')
            if ring.domain.is_negative(coeff):
                coeff = -coeff
            if coeff != 1 or expv == zm:
                scoeff = sstr(coeff)
            else:
                scoeff = ''
            sexpv = []
            for i in xrange(ngens):
                exp = expv[i]
                if not exp:
                    continue
                if exp != 1:
                    sexpv.append('%s**%d' % (sgens[i], exp))
                else:
                    sexpv.append('%s' % sgens[i])
            if scoeff:
                sexpv = [scoeff] + sexpv
            sexpvs.append('*'.join(sexpv))
        if sexpvs[0] in [" + ", " - "]:
            head = sexpvs.pop(0)
            if head == " - ":
                sexpvs.insert(0, "-")
        return ''.join(sexpvs)

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
        if not p2 or p2 == 0:
            return not p1
        ring1 = p1.ring
        if isinstance(p2, PolyElement):
            return dict.__eq__(p1, p2)
        else:
            zm = ring1.zero_monom
            if zm not in p1 or len(p1) > 1:
                return False
            return p1[zm] == p2 # ring1.domain_new(p2)

    def drop(self, gen):
        i, ring = self.ring._drop(gen)
        poly = PolyElement(ring)
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
        ring1 = p1.ring
        zm = ring1.zero_monom
        if isinstance(p2, PolyElement):
            if ring1 == p2.ring:
                p = PolyElement(ring1)
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
            elif p2.ring.__class__ == ring1.domain.__class__ and p2.ring == ring1.domain:
                p = p1.copy()
                if zm not in list(p1.keys()):
                    p[zm] = ring1.domain_new(p2)
                else:
                    if p2 == -p[zm]:
                        del p[zm]
                    else:
                        p[zm] += p2
                return p
            elif ring1.__class__ == p2.ring.domain.__class__ and ring1 == p2.ring.domain:
                return p2 + p1
            else:
                raise ValueError('cannot sum p1 and p2')
        # assume p2 in a number
        else:
            p = p1.copy()
            cp2 = ring1.domain_new(p2)
            if not cp2:
                return p
            if zm not in list(p1.keys()):
                p[zm] = cp2
            else:
                if p2 == -p[zm]:
                    del p[zm]
                else:
                    p[zm] += cp2
            return p

    def __radd__(p1, n):
        # assume n is in p1.ring.domain
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
        ring1 = p1.ring
        mz = ring1.zero_monom
        p = PolyElement(ring1)
        if isinstance(p2, PolyElement):
            if ring1 == p2.ring:
                for k in p1:
                    if k in p2:
                        r = p1[k] - p2[k]
                        if r:
                            p[k] = r
                    else:
                        p[k] = p1[k]
                for k in p2:
                    if k not in p1:
                        p[k] = -p2[k]
                return p
            elif p2.ring.__class__ == ring1.domain.__class__ and p2.ring == ring1.domain:
                p = p1.copy()
                if mz not in list(p1.keys()):
                    p[mz] = -ring1.domain_new(p2)
                else:
                    if p2 == p[mz]:
                        del p[mz]
                    else:
                        p[mz] -= p2
                return p
            else:
                raise ValueError('cannot coerce p2')
        # assume p2 in a number
        else:
            p2 = ring1.domain_new(p2)
            p = copy(p1)
            if mz not in list(p1.keys()):
                p[mz] = -p2
            else:
                if p2 == p[mz]:
                    del p[mz]
                else:
                    p[mz] -= p2
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
        p = PolyElement(p1.ring)
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
        ring1 = p1.ring
        p = PolyElement(ring1)
        if not p2:
            return p
        if isinstance(p2, PolyElement):
            if ring1 == p2.ring:
                get = p.get
                p2it = p2.items()
                for exp1, v1 in p1.iteritems():
                    for exp2, v2 in p2it:
                        exp = monomial_mul(exp1, exp2)
                        p[exp] = get(exp, 0) + v1*v2
                p.strip_zero()
                return p
            ring2 = p2.ring
            if ring2.__class__ != ring1.domain.__class__ or ring2 != ring1.domain:
                if ring1.__class__ == ring2.domain.__class__ and ring1 == ring2.domain:
                    p = PolyElement(p2.ring)
                    for exp2, v2 in p2.iteritems():
                        p[exp2] = p1*v2
                    return p
                else:
                    raise ValueError('p1 and p2 must have the same ring')
        # assume p2 in a number
        p2 = ring1.domain_new(p2)
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
        p = PolyElement(p1.ring)
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
                p = PolyElement(ring)
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
            p = PolyElement(ring)
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
            return copy(self)
        elif n == 2:
            return self.square()
        elif n == 3:
            return self*self.square()
        # TODO if ring.SR then use in some cases multinomial coefficients
        if ring.ngens == 1 and n >= 20 and ring.domain in (ZZ, QQ):
            return self.pow_miller(n)
        p = ring(1)
        while 1:
            if n&1:
                p = p*self
                n -= 1
                if not n:
                    break
            self = self.square()
            n = n // 2
        return p

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
        #if not ring.commuting:
        #   return self*self
        p = PolyElement(ring)
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
        ring1 = p1.ring
        ground_quo = ring1.domain.quo
        if isinstance(p2, PolyElement):
            if ring1 == p2.ring:
                if len(p2) == 1:
                    m = p2.keys()[0]
                    p = ring1(0)
                    c = p2.values()[0]
                    for k, v in p1.iteritems():
                        k1 = monomial_ldiv(k, m)
                        p[tuple(k1)] = ground_quo(v, c)
                    return p
                q, r = p1.div([p2])
                if r:
                    raise NotImplementedError('__div__ performs only division without remainder')
                return q[0]
            elif p2.ring.__class__ == ring1.domain.__class__ and p2.ring == ring1.domain:
                zm = p2.ring.zero_monom
                p = PolyElement(ring1)
                # if p is not a constant, not implemented
                if p2.keys() != [zm]:
                    raise NotImplementedError
                else:
                    p2 = p2[zm]
            else:
                raise NotImplementedError('cannot divide p1 by p2')
        # assume p2 in a number
        p = PolyElement(ring1)
        if not p2:
            raise ZeroDivisionError
        for exp1, v in p1.iteritems():
            coeff = ground_quo(v, p2)
            if coeff:
                p[exp1] = coeff
        return p

    __div__ = __truediv__

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
        #if self.is_laurent(*gens):
        #    raise NotImplementedError('self is a Laurent polynomial')
        # if any(p.is_laurent(*gens) for p in fv):
        #    raise NotImplementedError('there is a Laurent polynomial in fv')
        s = len(fv)
        qv = [PolyElement(ring) for i in range(s)]
        p = self.copy()
        r = PolyElement(ring)
        expvs = [fx.leading_expv() for fx in fv]
        # rn = range(ring.ngens)
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p.leading_expv()
                term = term_div((expv, p[expv]), (expvs[i], fv[i][expvs[i]]), domain)
                if term is not None:
                    expv1, c = term
                # expv1 = monomial_ldiv(expv, expvs[i])
                # if all(expv1[j] >= 0 for j in rn):
                #     c = p[expv]/fv[i][expvs[i]]        # XXX: ground_quo
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

    def quo(f, fv):
        return f.div(fv)[0]

    def rem(f, G):
        domain = f.ring.domain
        order = f.ring.order
        r = PolyElement(f.ring)
        ltf = f.LT
        f = f.copy()
        get = f.get
        while f:
            for g in G:
                tq = term_div(ltf, g.LT, domain)
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

    @property
    def leading_monom(self):
        p = PolyElement(self.ring)
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
        p = PolyElement(self.ring)
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
            return PolyElement(f.ring)

        terms = [ (monomial_mul(f_monom, monom), f_coeff*coeff) for f_monom, f_coeff in f.iteritems() ]
        return PolyElement(f.ring, terms)

    def mul_monom(f, monom):
        terms = [ (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in f.iteritems() ]
        return PolyElement(f.ring, terms)

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
        return PolyElement(f.ring, terms)

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
            terms = [ (monom, coeff // x) for monom, coeff in f.terms() ]

        return PolyElement(f.ring, terms)

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

        poly = PolyElement(f.ring, terms)
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
            h = PolyElement(ring)

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

    def cofactors(f, g):
        one, zero = f.ring.one, f.ring.zero

        if not f and not g:
            return zero, zero, zero
        elif not f:
            if g.LC >= 0:
                return g, zero, one
            else:
                return -g, zero, -one
        elif not g:
            if f.LC >= 0:
                return f, one, zero
            else:
                return -f, -one, zero

        J, (f, g) = f.deflate(g)
        h, cff, cfg = f._gcd(g)

        return (h.inflate(J), cff.inflate(J), cfg.inflate(J))

    def gcd(f, g):
        return f.cofactors(g)[0]

    def _gcd(f, g):
        ring = f.ring

        if ring.domain.is_QQ:
            return f._gcd_QQ(g)
        elif ring.domain.is_ZZ:
            return f._gcd_ZZ(g)
        else: # TODO: don't use dense representation (port PRS algorithms)
            return map(ring.from_dense, ring.dmp_inner_gcd(f.to_dense(), g.to_dense()))

    def _gcd_ZZ(f, g):
        return heugcd(f, g)

    def _gcd_QQ(f, g):
        ring = f.ring
        new_ring = ring.clone(domain=ring.domain.get_ring())

        cf, f = f.clear_denoms(f)
        cg, g = g.clear_denoms(g)

        f = f.set_ring(new_ring)
        g = g.set_ring(new_ring)

        h, cff, cfg = f._gcd_ZZ()

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
            new_ring = f.ring.clone(domain=domain.get_ring())

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

        p_neg = p.LC < 0
        q_neg = q.LC < 0

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
        i = ring.gens.index(x)
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
        i = ring.gens.index(x)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.terms():
                result += coeff*a**n

            return result
        else:
            poly = PolyElement(ring[1:])

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
        i = ring.gens.index(x)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.terms():
                result += coeff*a**n

            return ring.ground_new(result)
        else:
            poly = PolyElement(ring)

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
