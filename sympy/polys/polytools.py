"""User-friendly public interface to polynomial functions. """

from sympy.core import (
    S, Basic, Integer, sympify,
)

from sympy.core.decorators import (
    _sympifyit,
)

from sympy.polys.polyclasses import (
    GFP, DMP, SDP, DMF,
)

from sympy.polys.polyutils import (
    dict_from_basic,
    basic_from_dict,
    _sort_gens,
    _unify_gens,
    _dict_reorder,
    _update_args,
    _analyze_gens,
    _analyze_modulus,
    _analyze_extension,
    _dict_from_basic_no_gens,
)

from sympy.polys.groebnertools import (
    sdp_from_dict, sdp_groebner,
)

from sympy.polys.monomialtools import (
    monomial_cmp,
)

from sympy.polys.polyerrors import (
    OperationNotSupported, DomainError,
    CoercionFailed, UnificationFailed,
    GeneratorsNeeded, PolynomialError,
)

from sympy.utilities import any, all

import re

re_dom_poly = re.compile("^(Z|ZZ|Q|QQ)\[([^\]]+)\]$")
re_dom_frac = re.compile("^(Z|ZZ|Q|QQ)\(([^\]]+)\)$")
re_dom_GFp  = re.compile("^GF\((\d+)\)$")

from sympy.polys.algebratools import Algebra, ZZ, QQ, EX

def _analyze_domain(args):
    """Convert `domain` to an internal representation. """
    domain = args.get('domain')

    if domain is not None:
        domain = _parse_domain(domain, args)

    return domain

def _parse_domain(dom, args=None):
    """Make an instance out of a string representation of an algebra. """
    if isinstance(dom, Algebra):
        return dom

    if isinstance(dom, basestring):
        if dom in ['Z', 'ZZ']:
            return ZZ

        if dom in ['Q', 'QQ']:
            return QQ

        if dom == 'EX':
            return EX

        r = re.match(re_dom_poly, dom)

        if r is not None:
            ground, gens = r.groups()

            gens = map(sympify, gens.split(','))

            if ground in ['Z', 'ZZ']:
                return ZZ.poly_ring(*gens)
            else:
                return QQ.poly_ring(*gens)

        r = re.match(re_dom_frac, dom)

        if r is not None:
            ground, gens = r.groups()

            gens = map(sympify, gens.split(','))

            if ground in ['Z', 'ZZ']:
                return ZZ.frac_field(*gens)
            else:
                return QQ.frac_field(*gens)

        r = re.match(re_dom_GFp, dom)

        if r is not None:
            (modulus,) = r.groups()

            if args is not None:
                if not args.has_key('modulus'):
                    args['modulus'] = int(modulus)
                else:
                    raise PolynomialError('`modulus` already specified')

            return ZZ

    raise PolynomialError('expected a valid domain specification, got %s' % dom)

def _find_out_domain(rep, **args):
    """Finds the minimal domain that the coefficients of `rep` fit in. """
    field = args.get('field', False)
    has_rational = False

    for coeff in rep.itervalues():
        coeff = sympify(coeff)

        if coeff.is_Rational:
            if not coeff.is_Integer:
                has_rational = True
        elif coeff.is_Real:
            raise NotImplementedError('inexact coefficients')
        else:
            numers, denoms = [], []

            for coeff in rep.itervalues():
                num, den = coeff.as_numer_denom()

                try:
                    numers.append(_dict_from_basic_no_gens(num))
                except GeneratorsNeeded:
                    numers.append((num, None))

                try:
                    denoms.append(_dict_from_basic_no_gens(den))
                except GeneratorsNeeded:
                    denoms.append((den, None))

            gens = set([])

            for _, n_gens in numers:
                if n_gens is not None:
                    gens.update(n_gens)

            has_fractions = False

            for _, d_gens in denoms:
                if d_gens is not None:
                    gens.update(d_gens)
                    has_fractions = True

            gens = _sort_gens(gens, **args)

            if any(gen.is_Pow for gen in gens): # or I in gens
                # XXX: this should really go into algebraic function fields
                K, coeffs = EX, []

                for coeff in rep.itervalues():
                    coeffs.append(K.from_sympy(coeff))

                return EX, coeffs

            k, coeffs = len(gens), []

            if not field and not has_fractions:
                if all(den is S.One for den, _ in denoms):
                    K = ZZ.poly_ring(*gens)

                    for num, n_gens in numers:
                        if n_gens is not None:
                            n_monoms, n_coeffs = _dict_reorder(num, n_gens, gens)
                        else:
                            n_monoms, n_coeffs = [(0,)*k], [num]

                        n_coeffs = [ K.dom.from_sympy(c) for c in n_coeffs ]
                        coeffs.append(K(dict(zip(n_monoms, n_coeffs))))
                else:
                    K = QQ.poly_ring(*gens)

                    for (num, n_gens), (den, _) in zip(numers, denoms):
                        if n_gens is not None:
                            n_monoms, n_coeffs = _dict_reorder(num, n_gens, gens)
                            n_coeffs = [ coeff/den for coeff in n_coeffs ]
                        else:
                            n_monoms, n_coeffs = [(0,)*k], [num/den]

                        n_coeffs = [ K.dom.from_sympy(c) for c in n_coeffs ]
                        coeffs.append(K(dict(zip(n_monoms, n_coeffs))))

            else:
                K = ZZ.frac_field(*gens)

                for (num, n_gens), (den, d_gens) in zip(numers, denoms):
                    if n_gens is not None:
                        n_monoms, n_coeffs = _dict_reorder(num, n_gens, gens)
                    else:
                        n_monoms, n_coeffs = [(0,)*k], [num]

                    if d_gens is not None:
                        d_monoms, d_coeffs = _dict_reorder(den, d_gens, gens)
                    else:
                        d_monoms, d_coeffs = [(0,)*k], [den]

                    n_coeffs = [ K.dom.from_sympy(c) for c in n_coeffs ]
                    d_coeffs = [ K.dom.from_sympy(c) for c in d_coeffs ]

                    num = dict(zip(n_monoms, n_coeffs))
                    den = dict(zip(d_monoms, d_coeffs))

                    coeffs.append(K((num, den)))

            return K, coeffs

    if field or has_rational:
        K = QQ
    else:
        K = ZZ

    coeffs = []

    for coeff in rep.itervalues():
        coeffs.append(K.from_sympy(coeff))

    return K, coeffs

def _gen_to_level(f, gen):
    """Returns level associated with the given generator. """
    if isinstance(gen, int):
        return gen
    else:
        try:
            return list(f.gens).index(gen)
        except ValueError:
            raise PolynomialError("a valid generator expected, got %s" % gen)

def _init_poly_from_dict(rep, *gens, **args):
    """Initialize a Poly given a ... Poly instance. """
    domain = args.get('domain')
    modulus = args.get('modulus')

    if modulus is not None:
        if len(gens) != 1:
            raise GeneratorsNeeded("can't init from %s" % rep)

        if domain is not None:
            rep, _rep = {}, rep

            for k, v in _rep.iteritems():
                rep[k] = domain.convert(v)
        else:
            domain, coeffs = _find_out_domain(rep, **args)
            rep = dict(zip(rep.keys(), coeffs))
    else:
        if not gens:
            raise GeneratorsNeeded("can't init from %s" % rep)

        rep, _rep = {}, rep

        if domain is not None:
            for k, v in _rep.iteritems():
                if type(k) is not tuple:
                    rep[(k,)] = domain.convert(v)
                else:
                    rep[k] = domain.convert(v)
        else:
            domain, coeffs = _find_out_domain(_rep, **args)

            for k, v in zip(_rep.keys(), coeffs):
                if type(k) is not tuple:
                    rep[(k,)] = v
                else:
                    rep[k] = v

    if modulus is not None:
        rep = GFP(rep, modulus, domain)
    else:
        rep = DMP(rep, domain, len(gens)-1)

    return rep

def _init_poly_from_poly(rep, *gens, **args):
    """Initialize a Poly given a ... Poly instance. """
    domain = args.get('domain')
    modulus = args.get('modulus')

    rep_gens = rep.gens

    if isinstance(rep.rep, DMP):
        if not gens or rep.gens == gens:
            if domain is not None or modulus is not None:
                rep = rep.rep
            else:
                return rep
        else:
            if set(gens) != set(rep.gens):
                return Poly(rep.as_basic(), *gens, **args)
            else:
                monoms, coeffs = _dict_reorder(
                    rep.rep.to_dict(), rep.gens, gens)

                lev = len(gens) - 1

                _rep = dict(zip(monoms, coeffs))
                rep = DMP(_rep, rep.rep.dom, lev)

        if domain is not None:
            rep = rep.convert(domain)

        if modulus is not None:
            if not rep.lev and rep.dom.is_ZZ:
                rep = GFP(rep.rep, modulus, rep.dom)
            else:
                raise PolynomialError("can't make GFP out of %s" % rep)
    else:
        if not gens or rep.gens == gens:
            if modulus is not None:
                rep = GFP(rep.rep, modulus, rep.rep.dom)

            if domain is not None:
                rep = rep.convert(domain)
        else:
            raise PolynomialError("GFP multivariate polynomials are not supported")

    return rep, gens or rep_gens

def _init_poly_from_basic(ex, *gens, **args):
    """Initialize a Poly given a Basic expression. """
    if not gens:
        try:
            rep, gens = dict_from_basic(ex, **args)
        except GeneratorsNeeded:
            return ex
    else:
        rep = dict_from_basic(ex, gens, **args)

    domain = args.get('domain')
    modulus = args.get('modulus')

    if modulus is not None:
        if len(gens) > 1:
            raise PolynomialError("GFP multivariate polynomials not supported")

        if domain is None:
            domain = ZZ

        rep, _rep = {}, rep

        for (k,), v in _rep.iteritems():
            rep[k] = domain.from_sympy(v)

        rep = GFP(rep, modulus, domain)
    else:
        if domain is not None:
            for k, v in rep.iteritems():
                rep[k] = domain.from_sympy(v)
        else:
            domain, coeffs = _find_out_domain(rep, **args)
            rep = dict(zip(rep.keys(), coeffs))

        rep = DMP(rep, domain, len(gens)-1)

    return rep, gens

class Poly(Basic):
    """General class for representing polynomials in SymPy. """

    __slots__ = ['rep', 'gens']

    is_Poly = True

    def __new__(cls, rep, *gens, **args):
        """Create a new polynomial instance out of something useful. """
        if len(set(gens)) != len(gens):
            raise PolynomialError("duplicated generators: %s" % gens)

        if isinstance(rep, (DMP, GFP)):
            if rep.lev != len(gens)-1 or args:
                raise PolynomialError("invalid data for a polynomial")
        else:
            domain = _analyze_domain(args)
            modulus = _analyze_modulus(args)

            args = dict(args)

            if domain is not None:
                args['domain'] = domain
            if modulus is not None:
                args['modulus'] = modulus

            if type(rep) is dict:
                rep = _init_poly_from_dict(rep, *gens, **args)
            else:
                rep = sympify(rep)

                if rep.is_Poly:
                    result = _init_poly_from_poly(rep, *gens, **args)
                else:
                    result = _init_poly_from_basic(rep, *gens, **args)

                if type(result) is tuple:
                    rep, gens = result
                else:
                    return result

        obj = Basic.__new__(cls)

        obj.rep = rep
        obj.gens = gens

        return obj

    def __getnewargs__(self):
        """Data used by pickling protocol version 2. """
        return (self.rep, self.gens)

    def _hashable_content(self):
        """Allow SymPy to hash Poly instances. """
        return (self.rep, self.gens)

    @property
    def args(self):
        """Don't mess up with the core. """
        return [self.as_basic()]

    def unify(f, g):
        """Make `f` and `g` belong to the same domain. """
        if not f.is_Poly or not g.is_Poly:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        if isinstance(f.rep, DMP) and isinstance(g.rep, DMP):
            gens = _unify_gens(f.gens, g.gens)

            dom, lev = f.rep.dom.unify(g.rep.dom, gens), len(gens)-1

            if f.gens != gens:
                f_monoms, f_coeffs = _dict_reorder(f.rep.to_dict(), f.gens, gens)

                if f.rep.dom != dom:
                    f_coeffs = [ dom.convert(c, f.rep.dom) for c in f_coeffs ]

                F = DMP(dict(zip(f_monoms, f_coeffs)), dom, lev)
            else:
                F = f.rep.convert(dom)

            if g.gens != gens:
                g_monoms, g_coeffs = _dict_reorder(g.rep.to_dict(), g.gens, gens)

                if g.rep.dom != dom:
                    g_coeffs = [ dom.convert(c, g.rep.dom) for c in g_coeffs ]

                G = DMP(dict(zip(g_monoms, g_coeffs)), dom, lev)
            else:
                G = g.rep.convert(dom)
        elif isinstance(f.rep, GFP) and isinstance(g.rep, GFP):
            if f.gens != g.gens or f.mod != g.mod:
                raise UnificationFailed("can't unify %s with %s" % (f, g))
            else:
                F, G, gens = f.rep, g.rep, f.gens
        elif isinstance(f.rep, GFP) and isinstance(g.rep, DMP):
            if f.gens != g.gens:
                raise UnificationFailed("can't unify %s with %s" % (f, g))
            else:
                G, F, gens = GFP(g.rep.convert(f.rep.dom), f.mod, g.rep.dom), f.rep, f.gens
        elif isinstance(f.rep, DMP) and isinstance(g.rep, GFP):
            if f.gens != g.gens:
                raise UnificationFailed("can't unify %s with %s" % (f, g))
            else:
                F, G, gens = GFP(f.rep.convert(g.rep.dom), g.mod, g.rep.dom), g.rep, g.gens
        else:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        def per(rep, dom=dom, gens=gens, remove=None):
            if remove is not None:
                gens = gens[:remove]+gens[remove+1:]

                if not gens:
                    return dom.to_sympy(rep)

            return Poly(rep, *gens)

        return dom, per, F, G

    def per(f, rep, remove=None):
        """Create a Poly out of the given representation. """
        if remove is not None:
            gens = f.gens[:remove]+f.gens[remove+1:]

            if not gens:
                return f.rep.dom.to_sympy(rep)
        else:
            gens = f.gens

        return Poly(rep, *gens)

    def set_modulus(f, modulus):
        """Set modulus of `f`. """
        raise NotImplementedError

    def get_modulus(f):
        """Get modulus of `f`. """
        if isinstance(f.rep, GFP):
            return Integer(f.rep.mod)
        else:
            return None

    def set_domain(f, domain):
        """Set ground domain of `f`. """
        domain = _parse_domain(domain)
        return f.per(f.rep.convert(domain))

    def get_domain(f):
        """Get ground domain of `f`. """
        return f.rep.dom

    def to_ring(f):
        """Make the ground domain a field. """
        try:
            result = f.rep.to_ring()
        except AttributeError:
            raise OperationNotSupported(f, 'to_ring')

        return f.per(result)

    def to_field(f):
        """Make the ground domain a field. """
        try:
            result = f.rep.to_field()
        except AttributeError:
            raise OperationNotSupported(f, 'to_field')

        return f.per(result)

    def coeffs(f):
        """Returns all non-zero coefficients from `f` in lex order. """
        return [ f.rep.dom.to_sympy(c) for c in f.rep.coeffs() ]

    def monoms(f):
        """Returns all non-zero monomials from `f` in lex order. """
        return f.rep.monoms()

    def terms(f):
        """Returns all non-zero terms from `f` in lex order. """
        return [ (m, f.rep.dom.to_sympy(c)) for m, c in f.rep.terms() ]

    def all_coeffs(f):
        """Returns all coefficients from a univariate polynomial `f`. """
        return [ f.rep.dom.to_sympy(c) for c in f.rep.all_coeffs() ]

    def all_monoms(f):
        """Returns all monomials from a univariate polynomial `f`. """
        return f.rep.all_monoms()

    def all_terms(f):
        """Returns all terms from a univariate polynomial `f`. """
        return [ (m, f.rep.dom.to_sympy(c)) for m, c in f.rep.all_terms() ]

    def length(f):
        """Returns the number of non-zero terms in `f`. """
        return len(f.as_dict())

    def as_dict(f):
        """Switch to a dict representation with SymPy coefficients. """
        return f.rep.to_sympy_dict()

    def as_basic(f):
        """Convert a polynomial instance to a SymPy expression. """
        return basic_from_dict(f.rep.to_sympy_dict(), *f.gens)

    def deflate(f):
        """Reduce degree of `f` by mapping `x_i**m` to `y_i`. """
        try:
            J, result = f.rep.deflate()
        except AttributeError:
            raise OperationNotSupported(f, 'deflate')

        return J, f.per(result)

    def terms_gcd(f):
        """Remove GCD of terms from the polynomial `f`. """
        try:
            J, result = f.rep.terms_gcd()
        except AttributeError:
            raise OperationNotSupported(f, 'terms_gcd')

        return J, f.per(result)

    def abs(f):
        """Make all coefficients in `f` positive. """
        try:
            result = f.rep.abs()
        except AttributeError:
            raise OperationNotSupported(f, 'abs')

        return f.per(result)

    def neg(f):
        """Negate all cefficients in `f`. """
        try:
            result = f.rep.neg()
        except AttributeError:
            raise OperationNotSupported(f, 'neg')

        return f.per(result)

    def add(f, g):
        """Add two polynomials `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.add(G)
        except AttributeError:
            raise OperationNotSupported(f, 'add')

        return per(result)

    def sub(f, g):
        """Subtract two polynomials `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.sub(G)
        except AttributeError:
            raise OperationNotSupported(f, 'sub')

        return per(result)

    def mul(f, g):
        """Multiply two polynomials `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.mul(G)
        except AttributeError:
            raise OperationNotSupported(f, 'mul')

        return per(result)

    def sqr(f):
        """Square a polynomial `f`. """
        try:
            result = f.rep.sqr()
        except AttributeError:
            raise OperationNotSupported(f, 'sqr')

        return f.per(result)

    def pow(f, n):
        """Raise `f` to a non-negative power `n`. """
        n = int(n)

        try:
            result = f.rep.pow(n)
        except AttributeError:
            raise OperationNotSupported(f, 'pow')

        return f.per(result)

    def pdiv(f, g):
        """Polynomial pseudo-division of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            q, r = F.pdiv(G)
        except AttributeError:
            raise OperationNotSupported(f, 'pdiv')

        return per(q), per(r)

    def prem(f, g):
        """Polynomial pseudo-remainder of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.prem(G)
        except AttributeError:
            raise OperationNotSupported(f, 'prem')

        return per(result)

    def pquo(f, g):
        """Polynomial pseudo-quotient of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.pquo(G)
        except AttributeError:
            raise OperationNotSupported(f, 'pquo')

        return per(result)

    def pexquo(f, g):
        """Polynomial exact pseudo-quotient of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.pexquo(G)
        except AttributeError:
            raise OperationNotSupported(f, 'pexquo')

        return per(result)

    def div(f, g):
        """Polynomial division with remainder of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            q, r = F.div(G)
        except AttributeError:
            raise OperationNotSupported(f, 'div')

        return per(q), per(r)

    def rem(f, g):
        """Computes polynomial remainder of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.rem(G)
        except AttributeError:
            raise OperationNotSupported(f, 'rem')

        return per(result)

    def quo(f, g):
        """Computes polynomial quotient of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.quo(G)
        except AttributeError:
            raise OperationNotSupported(f, 'quo')

        return per(result)

    def exquo(f, g):
        """Computes polynomial exact quotient of `f` by `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.exquo(G)
        except AttributeError:
            raise OperationNotSupported(f, 'exquo')

        return per(result)

    def degree(f, gen=0):
        """Returns degree of `f` in `x_j`. """
        j = _gen_to_level(f, gen)

        try:
            return f.rep.degree(j)
        except AttributeError:
            raise OperationNotSupported(f, 'degree')

    def degree_list(f):
        """Returns a list of degrees of `f`. """
        try:
            return f.rep.degree_list()
        except AttributeError:
            raise OperationNotSupported(f, 'degree_list')

    def total_degree(f):
        """Returns the total degree of `f`. """
        try:
            return f.rep.total_degree()
        except AttributeError:
            raise OperationNotSupported(f, 'total_degree')

    def LC(f):
        """Returns the leading coefficent of `f`. """
        try:
            result = f.rep.LC()
        except AttributeError:
            raise OperationNotSupported(f, 'LC')

        return f.rep.dom.to_sympy(result)

    def TC(f):
        """Returns the trailing coefficent of `f`. """
        try:
            result = f.rep.TC()
        except AttributeError:
            raise OperationNotSupported(f, 'TC')

        return f.rep.dom.to_sympy(result)

    def EC(f):
        """Returns the last non-zero coefficent of `f`. """
        try:
            return f.coeffs()[-1]
        except AttributeError:
            raise OperationNotSupported(f, 'EC')

    def nth(f, *N):
        """Returns the `n`-th coefficient of `f`. """
        try:
            result = f.rep.nth(*map(int, N))
        except AttributeError:
            raise OperationNotSupported(f, 'nth')

        return f.rep.dom.to_sympy(result)

    def LM(f):
        """Returns the leading monomial of `f`. """
        return f.monoms()[0]

    def EM(f):
        """Returns the last non-zero monomial of `f`. """
        return f.monoms()[-1]

    def LT(f):
        """Returns the leading term of `f`. """
        return f.terms()[0]

    def ET(f):
        """Returns the last non-zero term of `f`. """
        return f.terms()[-1]

    def max_norm(f):
        """Returns maximum norm of `f`. """
        try:
            result = f.rep.max_norm()
        except AttributeError:
            raise OperationNotSupported(f, 'max_norm')

        return f.rep.dom.to_sympy(result)

    def l1_norm(f):
        """Returns l1 norm of `f`. """
        try:
            result = f.rep.l1_norm()
        except AttributeError:
            raise OperationNotSupported(f, 'l1_norm')

        return f.rep.dom.to_sympy(result)

    def ground_to_ring(f):
        """Clear denominators, but keep the ground domain. """
        try:
            coeff, result = f.rep.ground_to_ring()
        except AttributeError:
            raise OperationNotSupported(f, 'ground_to_ring')

        return f.rep.dom.to_sympy(coeff), f.per(result)

    def integrate(f, *specs, **args):
        """Computes indefinite integral of `f`. """
        if args.get('auto', True) and not f.rep.dom.has_Field:
            f = f.to_field()

        try:
            if not specs:
                return f.per(f.rep.integrate(m=1))

            rep = f.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.integrate(int(m), _gen_to_level(f, gen))

            return f.per(rep)
        except AttributeError:
            raise OperationNotSupported(f, 'integrate')

    def diff(f, *specs):
        """Computes partial derivative of `f`. """
        try:
            if not specs:
                return f.per(f.rep.diff(m=1))

            rep = f.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.diff(int(m), _gen_to_level(f, gen))

            return f.per(rep)
        except AttributeError:
            raise OperationNotSupported(f, 'diff')

    def eval(f, a, gen=0):
        """Evaluates `f` at `a` in the given variable. """
        j = _gen_to_level(f, gen)

        try:
            result = f.rep.eval(a, j)
        except AttributeError:
            raise OperationNotSupported(f, 'eval')

        return f.per(result, remove=j)

    def half_gcdex(f, g, **args):
        """Half extended Euclidean algorithm of `f` and `g`. """
        dom, per, F, G = f.unify(g)

        if args.get('auto', True) and not dom.has_Field:
            F, G = F.to_field(), G.to_field()

        try:
            s, h = F.half_gcdex(G)
        except AttributeError:
            raise OperationNotSupported(f, 'half_gcdex')

        return per(s), per(h)

    def gcdex(f, g, **args):
        """Extended Euclidean algorithm of `f` and `g`. """
        dom, per, F, G = f.unify(g)

        if args.get('auto', True) and not dom.has_Field:
            F, G = F.to_field(), G.to_field()

        try:
            s, t, h = F.gcdex(G)
        except AttributeError:
            raise OperationNotSupported(f, 'gcdex')

        return per(s), per(t), per(h)

    def invert(f, g, **args):
        """Invert `f` modulo `g`, if possible. """
        dom, per, F, G = f.unify(g)

        if args.get('auto', True) and not dom.has_Field:
            F, G = F.to_field(), G.to_field()

        try:
            result = F.invert(G)
        except AttributeError:
            raise OperationNotSupported(f, 'invert')

        return per(result)

    def subresultants(f, g):
        """Computes subresultant PRS sequence of `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.subresultants(G)
        except AttributeError:
            raise OperationNotSupported(f, 'subresultants')

        return map(per, result)

    def resultant(f, g):
        """Computes resultant of `f` and `g` via PRS. """
        _, per, F, G = f.unify(g)

        try:
            result = F.resultant(G)
        except AttributeError:
            raise OperationNotSupported(f, 'resultant')

        return per(result, remove=0)

    def discriminant(f):
        """Computes discriminant of `f`. """
        try:
            result = f.rep.discriminant()
        except AttributeError:
            raise OperationNotSupported(f, 'discriminant')

        return f.per(result, remove=0)

    def cofactors(f, g):
        """Returns GCD of `f` and `g` and their cofactors. """
        _, per, F, G = f.unify(g)

        try:
            h, cff, cfg = F.cofactors(G)
        except AttributeError:
            raise OperationNotSupported(f, 'cofactors')

        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """Returns polynomial GCD of `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.gcd(G)
        except AttributeError:
            raise OperationNotSupported(f, 'gcd')

        return per(result)

    def lcm(f, g):
        """Returns polynomial LCM of `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.lcm(G)
        except AttributeError:
            raise OperationNotSupported(f, 'lcm')

        return per(result)

    def monic(f):
        """Divides all coefficients by `LC(f)`. """
        try:
            result = f.rep.monic()
        except AttributeError:
            raise OperationNotSupported(f, 'monic')

        return f.per(result)

    def content(f):
        """Returns GCD of polynomial coefficients. """
        try:
            result = f.rep.monic()
        except AttributeError:
            raise OperationNotSupported(f, 'content')

        return f.rep.dom.to_sympy(result)

    def primitive(f):
        """Returns content and a primitive form of `f`. """
        try:
            cont, result = f.rep.primitive()
        except AttributeError:
            raise OperationNotSupported(f, 'primitive')

        return f.rep.dom.to_sympy(cont), f.per(result)

    def compose(f, g):
        """Computes functional composition of `f` and `g`. """
        _, per, F, G = f.unify(g)

        try:
            result = F.compose(G)
        except AttributeError:
            raise OperationNotSupported(f, 'compose')

        return per(result)

    def decompose(f):
        """Computes functional decomposition of `f`. """
        try:
            result = f.rep.decompose()
        except AttributeError:
            raise OperationNotSupported(f, 'decompose')

        return map(f.per, result)

    def sturm(f, **args):
        """Computes the Sturm sequence of `f`. """
        if args.get('auto', True):
            f = f.to_field()

        try:
            result = f.rep.sturm()
        except AttributeError:
            raise OperationNotSupported(f, 'sturm')

        return map(f.per, result)

    def sqf_part(f):
        """Computes square-free part of `f`. """
        try:
            result = f.rep.sqf_part()
        except AttributeError:
            raise OperationNotSupported(f, 'sqf_part')

        return f.per(result)

    def sqf_list(f, **args):
        """Returns a list of square-free factors of `f`. """
        try:
            result = f.rep.sqf_list(**args)
        except AttributeError:
            raise OperationNotSupported(f, 'sqf_list')

        if type(result) is not tuple:
            return [ (f.per(g), k) for g, k in result ]
        else:
            coeff = f.rep.dom.to_sympy(result[0])
            factors = [ (f.per(g), k) for g, k in result[1] ]

            return coeff, factors

    def factor_list(f, **args):
        """Returns a list of irreducible factors of `f`. """
        try:
            coeff, factors = f.rep.factor_list(**args)
        except DomainError:
            return S.One, [(f, 1)]
        except AttributeError:
            raise OperationNotSupported(f, 'factor_list')

        coeff = f.rep.dom.to_sympy(coeff)
        factors = [ (f.per(g), k) for g, k in factors ]

        return coeff, factors

    def cancel(f, g):
        """Cancel common factors in a rational function `f/g`.  """
        dom, per, F, G = f.unify(g)

        if F.is_zero or G.is_zero:
            return S.One, per(F), per(G)

        if dom.has_Field and not dom.is_EX:
            cF, F = F.ground_to_ring()
            cG, G = G.ground_to_ring()

            F = F.to_ring()
            G = G.to_ring()

        try:
            _, P, Q = F.cofactors(G)
        except AttributeError:
            raise OperationNotSupported(f, 'cofactors')

        if dom.has_Field:
            if not dom.is_EX:
                P, Q = P.to_field(), Q.to_field()

                cF = dom.to_sympy(cF)
                cG = dom.to_sympy(cG)

                return cG/cF, per(P), per(Q)
            else:
                P, Q = P.monic(), Q.monic()

        return S.One, per(P), per(Q)

    @property
    def is_zero(f):
        """Returns `True` if `f` is a zero polynomial. """
        return f.rep.is_zero

    @property
    def is_one(f):
        """Returns `True` if `f` is a unit polynomial. """
        return f.rep.is_one

    @property
    def is_sqf(f):
        """Returns `True` if `f` is a square-free polynomial. """
        return f.rep.is_sqf

    @property
    def is_monic(f):
        """Returns `True` if the leading coefficient of `f` is one. """
        return f.rep.is_monic

    @property
    def is_primitive(f):
        """Returns `True` if GCD of coefficients of `f` is one. """
        return f.rep.is_primitive

    @property
    def is_ground(f):
        """Returns `True` if `f` is an element of the ground domain. """
        return f.rep.is_ground

    @property
    def is_linear(f):
        """Returns `True` if `f` is linear in all its variables. """
        return f.rep.is_linear

    @property
    def is_monomial(f):
        """Returns `True` if `f` is zero or has only one term. """
        return f.length() <= 1

    @property
    def is_homogeneous(f):
        """Returns `True` if `f` has zero trailing coefficient. """
        return f.rep.is_homogeneous

    @property
    def is_univariate(f):
        """Returns `True` if `f` is an univariate polynomial. """
        return len(f.gens) == 1

    @property
    def is_multivariate(f):
        """Returns `True` if `f` is a multivariate polynomial. """
        return len(f.gens) != 1

    def __abs__(f):
        return f.abs()

    def __neg__(f):
        return f.neg()

    @_sympifyit('g', NotImplemented)
    def __add__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_basic() + g

        return f.add(g)

    @_sympifyit('g', NotImplemented)
    def __radd__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g + f.as_basic()

        return g.add(f)

    @_sympifyit('g', NotImplemented)
    def __sub__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_basic() - g

        return f.sub(g)

    @_sympifyit('g', NotImplemented)
    def __rsub__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g - f.as_basic()

        return g.sub(f)

    @_sympifyit('g', NotImplemented)
    def __mul__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_basic()*g

        return f.mul(g)

    @_sympifyit('g', NotImplemented)
    def __rmul__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g*f.as_basic()

        return g.mul(f)

    @_sympifyit('n', NotImplemented)
    def __pow__(f, n):
        if n.is_Integer and n >= 0:
            return f.pow(n)
        else:
            return f.as_basic()**n

    @_sympifyit('g', NotImplemented)
    def __divmod__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return divmod(f.as_basic(), g)

        return f.div(g)

    @_sympifyit('g', NotImplemented)
    def __rdivmod__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return divmod(g, f.as_basic())

        return g.div(f)

    @_sympifyit('g', NotImplemented)
    def __mod__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_basic() % g

        return f.rem(g)

    @_sympifyit('g', NotImplemented)
    def __mod__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g % f.as_basic()

        return g.rem(f)

    @_sympifyit('g', NotImplemented)
    def __floordiv__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_basic() // g

        return f.exquo(g)

    @_sympifyit('g', NotImplemented)
    def __rfloordiv__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g // f.as_basic()

        return g.exquo(f)

    @_sympifyit('g', NotImplemented)
    def __eq__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens) # XXX: args
            except PolynomialError:
                return False

        return f.rep == g.rep

    @_sympifyit('g', NotImplemented)
    def __ne__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens) # XXX: args
            except PolynomialError:
                return True

        return f.rep != g.rep

    def __nonzero__(f):
        return not f.is_zero

def poly_coerce(f, g, *gens, **args):
    """Cooperatively make polynomials out of `f` and `g`. """
    if gens:
        return (Poly(f, *gens, **args),
                Poly(g, *gens, **args))
    else:
        F = Poly(f, **args)
        G = Poly(g, **args)

        if F.is_Poly:
            if G.is_Poly:
                return F, G
            else:
                return F, Poly(g, *F.gens, **args)
        else:
            if G.is_Poly:
                return Poly(f, *G.gens, **args), G
            else:
                raise CoercionFailed(F, G)

def _all_is_poly(*polys):
    """Returns `True` if given at least one Poly instance. """
    return all(isinstance(poly, Poly) for poly in polys)

def _should_return_basic(*polys, **args):
    """Figure out if results should be returned as basic. """
    arg = args.get('polys')

    if arg is not None:
        return not arg
    else:
        return not _all_is_poly(*polys)

def pdiv(f, g, *gens, **args):
    """Polynomial pseudo-division of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute pseudo division of %s and %s without generators" % (f, g))

    q, r = F.pdiv(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic(), r.as_basic()
    else:
        return q, r

def prem(f, g, *gens, **args):
    """Polynomial pseudo-remainder of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute pseudo remainder of %s and %s without generators" % (f, g))

    r = F.prem(G)

    if _should_return_basic(f, g, **args):
        return r.as_basic()
    else:
        return r

def pquo(f, g, *gens, **args):
    """Polynomial pseudo-quotient of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute pseudo quotient of %s and %s without generators" % (f, g))

    q = F.pquo(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic()
    else:
        return q

def pexquo(f, g, *gens, **args):
    """Polynomial exact pseudo-quotient of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute pseudo quotient of %s and %s without generators" % (f, g))

    q = F.pexquo(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic()
    else:
        return q

def div(f, g, *gens, **args):
    """Polynomial division with remainder of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute division of %s and %s without generators" % (f, g))

    q, r = F.div(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic(), r.as_basic()
    else:
        return q, r

def rem(f, g, *gens, **args):
    """Computes polynomial remainder of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute remainder of %s and %s without generators" % (f, g))

    r = F.rem(G)

    if _should_return_basic(f, g, **args):
        return r.as_basic()
    else:
        return r

def quo(f, g, *gens, **args):
    """Computes polynomial quotient of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute quotient of %s and %s without generators" % (f, g))

    q = F.quo(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic()
    else:
        return q

def exquo(f, g, *gens, **args):
    """Computes polynomial exact quotient of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute quotient of %s and %s without generators" % (f, g))

    q = F.exquo(G)

    if _should_return_basic(f, g, **args):
        return q.as_basic()
    else:
        return q

def half_gcdex(f, g, *gens, **args):
    """Half extended Euclidean algorithm of `f` and `g`. """
    args = _update_args(args, 'field', True)
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.half_gcdex(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute half extended GCD of %s and %s without generators" % (f, g))

    s, t, h = F.half_gcdex(G, **args)

    if _should_return_basic(f, g, **args):
        return q.as_basic(), t.as_basic(), h.as_basic()
    else:
        return s, t, h

def gcdex(f, g, *gens, **args):
    """Extended Euclidean algorithm of `f` and `g`. """
    args = _update_args(args, 'field', True)
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.gcdex(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute extended GCD of %s and %s without generators" % (f, g))

    s, h = F.gcdex(G, **args)

    if _should_return_basic(f, g, **args):
        return s.as_basic(), h.as_basic()
    else:
        return s, h

def invert(f, g, *gens, **args):
    """Invert `f` modulo `g`, if possible. """
    args = _update_args(args, 'field', True)
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.invert(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute inversion of %s modulo %s without generators" % (f, g))

    h = F.invert(G, **args)

    if _should_return_basic(f, g, **args):
        return h.as_basic()
    else:
        return h

def subresultants(f, g, *gens, **args):
    """Computes subresultant PRS sequence of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute subresultants of %s and %s without generators" % (f, g))

    result = F.subresultants(G)

    if _should_return_basic(f, g, **args):
        return [ r.as_basic() for r in result ]
    else:
        return result

def resultant(f, g, *gens, **args):
    """Computes resultant of `f` and `g` via PRS. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute resultant of %s and %s without generators" % (f, g))

    result = F.resultant(G)

    if _should_return_basic(f, g, **args):
        return result.as_basic()
    else:
        return result

def discriminant(f, *gens, **args):
    """Computes discriminant of `f`. """
    gens = _analyze_gens(gens)
    F = Poly(f, *gens, **args)

    result = F.discriminant()

    if _should_return_basic(f, **args):
        return result.as_basic()
    else:
        return result

def cofactors(f, g, *gens, **args):
    """Returns GCD of `f` and `g` and their cofactors. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.cofactors(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute cofactors of %s and %s without generators" % (f, g))

    h, cff, cfg = F.cofactors(G)

    if _should_return_basic(f, g, **args):
        return h.as_basic(), cff.as_basic(), cfg.as_basic()
    else:
        return h, cff, cfg

def gcd(f, g, *gens, **args):
    """Returns polynomial GCD of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.gcd(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute GCD of %s and %s without generators" % (f, g))

    result = F.gcd(G)

    if _should_return_basic(f, g, **args):
        return result.as_basic()
    else:
        return result

def lcm(f, g, *gens, **args):
    """Returns polynomial LCM of `f` and `g`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        try:
            return f.lcm(g)
        except (AttributeError, TypeError):
            raise GeneratorsNeeded("can't compute LCM of %s and %s without generators" % (f, g))

    result = F.lcm(G)

    if _should_return_basic(f, g, **args):
        return result.as_basic()
    else:
        return result

def monic(f, *gens, **args):
    """Divides all coefficients by `LC(f)`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if _should_return_basic(f, **args):
        return F.monic().as_basic()
    else:
        return F.monic()

def content(f, *gens, **args):
    """Returns GCD of polynomial coefficients. """
    return Poly(f, *gens, **args).content()

def primitive(f, *gens, **args):
    """Returns content and a primitive form of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    cont, result = F.primitive()

    if _should_return_basic(f, **args):
        return cont, result.as_basic()
    else:
        return cont, result

def compose(f, g, *gens, **args):
    """Returns functional composition `f(g)`. """
    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(f, g, *gens, **args)
    except CoercionFailed, (f, g):
        raise GeneratorsNeeded("can't compute composition of %s and %s without generators" % (f, g))

    result = F.compose(G)

    if _should_return_basic(f, g, **args):
        return result.as_basic()
    else:
        return result

def decompose(f, *gens, **args):
    """Computes functional decomposition of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        raise GeneratorsNeeded("can't compute functional decomposition of %s without generators" % f)

    result = F.decompose()

    if _should_return_basic(f, **args):
        return [ r.as_basic() for r in result ]
    else:
        return result

def sturm(f, *gens, **args):
    """Computes the Sturm sequence of `f`. """
    args = _update_args(args, 'field', True)
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        raise GeneratorsNeeded("can't compute Sturm sequence of %s without generators" % f)

    result = F.sturm()

    if _should_return_basic(f, **args):
        return [ r.as_basic() for r in result ]
    else:
        return result

def sqf_part(f, *gens, **args):
    """Computes square-free part of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        raise GeneratorsNeeded("can't compute square-free part of %s without generators" % f)

    result = F.sqf_part()

    if _should_return_basic(f, **args):
        return result.as_basic()
    else:
        return result

def sqf_list(f, *gens, **args):
    """Returns a list of square-free factors of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        raise GeneratorsNeeded("can't compute square-free decomposition of %s without generators" % f)

    result = F.sqf_list(**args)

    if _should_return_basic(f, **args):
        if type(result) is tuple:
            return result[0], [ (g.as_basic(), k) for g, k in result[1] ]
        else:
            return [ (g.as_basic(), k) for g, k in result ]
    else:
        return result

def sqf(f, *gens, **args):
    """Returns square-free decomposition of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        raise GeneratorsNeeded("can't compute square-free decomposition of %s without generators" % f)

    (coeff, factors), result = F.sqf_list(), S.One

    for g, k in factors:
        result *= g.as_basic()**k

    return coeff*result

def factor_list(f, *gens, **args):
    """Returns a list of irreducible factors of `f`. """
    F = Poly(f, *_analyze_gens(gens), **args)

    if not F.is_Poly:
        return F, []

    coeff, factors = F.factor_list(**args)

    if _should_return_basic(f, **args):
        return coeff, [ (g.as_basic(), k) for g, k in factors ]
    else:
        return coeff, factors

def factor(f, *gens, **args):
    """Returns factorization into irreducibles of `f`. """
    frac = args.get('frac', False)
    gens = _analyze_gens(gens)

    f = sympify(f)

    if not frac:
        F = Poly(f, *gens, **args)

        if not F.is_Poly:
            return f

        coeff, factors = F.factor_list(**args)

        result = S.One

        for g, k in factors:
            result *= g.as_basic()**k

        return coeff*result
    else:
        p, q = f.as_numer_denom()

        try:
            p, q = poly_coerce(p, q, *gens, **args)
        except CoercionFailed:
            return f

        coeff_p, factors_p = p.factor_list(**args)
        coeff_q, factors_q = q.factor_list(**args)

        numer, denom = S.One, S.One

        for g, k in factors_p:
            numer *= g.as_basic()**k

        for g, k in factors_q:
            denom *= g.as_basic()**k

        return (coeff_p/coeff_q)*(numer/denom)

def cancel(f, *gens, **args):
    """Cancel common factors in a rational function `f`.  """
    f = sympify(f)

    if type(f) is not tuple:
        if f.is_Number:
            return f
        else:
            p, q = f.as_numer_denom()
    else:
        p, q = f

    gens = _analyze_gens(gens)

    try:
        F, G = poly_coerce(p, q, *gens, **args)
    except CoercionFailed:
        if type(f) is not tuple:
            return f
        else:
            return S.One, p, q

    c, P, Q = F.cancel(G)

    if type(f) is not tuple:
        return c*(P.as_basic()/Q.as_basic())
    else:
        if _should_return_basic(p, q, **args):
            return c, P.as_basic(), Q.as_basic()
        else:
            return c, P, Q

def groebner(F, *gens, **args):
    """Computes reduced Groebner basis for a set of polynomials. """
    if not F:
        return []

    order = args.pop('order', 'lex')
    basic = _should_return_basic(*F, **args)

    args = _update_args(args, 'field', True)
    gens = _analyze_gens(gens)

    if gens:
        F = [ Poly(f, *gens, **args) for f in F ]

        if any(f.rep.dom.is_EX for f in F):
            F = [ f.set_domain('EX') for f in F ]
    else:
        raise GeneratorsNeeded("can't compute Groebner basis without generators")

    gens, lev, dom = F[0].gens, F[0].rep.lev, F[0].rep.dom

    if isinstance(order, basestring):
        order = monomial_cmp(order)

    F = [ sdp_from_dict(f.rep.to_dict(), order) for f in F ]

    G = [ Poly(DMP(dict(g), dom, lev), *gens)
        for g in sdp_groebner(F, lev, order, dom) ]

    if basic:
        return [ g.as_basic() for g in G ]
    else:
        return G

