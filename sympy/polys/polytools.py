"""User-friendly public interface to polynomial functions. """

from sympy.core import (
    S, Basic, Expr, I, Integer, Add, Mul,
)

from sympy.core.sympify import (
    sympify, SympifyError,
)

from sympy.core.decorators import (
    _sympifyit,
)

from sympy.polys.polyclasses import (
    DMP, ANP, DMF,
)

from sympy.polys.polyutils import (
    basic_from_dict,
    _sort_gens,
    _unify_gens,
    _dict_reorder,
    _dict_from_expr,
    _parallel_dict_from_expr,
)

from sympy.polys.rationaltools import (
    together,
)

from sympy.polys.rootisolation import (
    dup_isolate_real_roots_list,
)

from sympy.polys.groebnertools import (
    sdp_from_dict, sdp_div, sdp_groebner,
)

from sympy.polys.monomialtools import (
    Monomial, monomial_key,
)

from sympy.polys.polyerrors import (
    OperationNotSupported, DomainError,
    CoercionFailed, UnificationFailed,
    GeneratorsNeeded, PolynomialError,
    PolificationFailed, FlagError,
    MultivariatePolynomialError,
    ExactQuotientFailed,
    ComputationFailed,
    GeneratorsError,
)

from sympy.polys.polycontext import (
    register_context,
)

from sympy.mpmath import (
    polyroots as npolyroots,
)

from sympy.utilities import (
    any, all, group,
)

from sympy.ntheory import isprime

import sympy.polys

from sympy.polys.domains import FF, QQ
from sympy.polys.constructor import construct_domain

from sympy.polys import polyoptions as options

class Poly(Expr):
    """Generic class for representing polynomials in SymPy. """

    __slots__ = ['rep', 'gens']

    is_Poly = True

    def __new__(cls, rep, *gens, **args):
        """Create a new polynomial instance out of something useful. """
        opt = options.build_options(gens, args)

        if 'order' in opt:
            raise NotImplementedError("'order' keyword is not implemented yet")

        if isinstance(rep, (dict, list)):
            if isinstance(rep, dict):
                return cls._from_dict(rep, opt)
            else:
                return cls._from_list(rep, opt)
        else:
            rep = sympify(rep)

            if rep.is_Poly:
                return cls._from_poly(rep, opt)
            else:
                return cls._from_expr(rep, opt)

    @classmethod
    def new(cls, rep, *gens):
        """Construct :class:`Poly` instance from raw representation. """
        if not isinstance(rep, DMP):
            raise PolynomialError("invalid polynomial representation: %s" % rep)
        elif rep.lev != len(gens)-1:
            raise PolynomialError("invalid arguments: %s, %s" % (rep, gens))

        obj = Basic.__new__(cls)

        obj.rep = rep
        obj.gens = gens

        return obj

    @classmethod
    def from_dict(cls, rep, *gens, **args):
        """Construct a polynomial from a ``dict``. """
        opt = options.build_options(gens, args)
        return cls._from_dict(rep, opt)

    @classmethod
    def from_list(cls, rep, *gens, **args):
        """Construct a polynomial from a ``list``. """
        opt = options.build_options(gens, args)
        return cls._from_list(rep, opt)

    @classmethod
    def from_poly(cls, rep, *gens, **args):
        """Construct a polynomial from a polynomial. """
        opt = options.build_options(gens, args)
        return cls._from_poly(rep, opt)

    @classmethod
    def from_expr(cls, rep, *gens, **args):
        """Construct a polynomial from an expression. """
        opt = options.build_options(gens, args)
        return cls._from_expr(rep, opt)

    @classmethod
    def _from_dict(cls, rep, opt):
        """Construct a polynomial from a ``dict``. """
        gens = opt.gens

        if not gens:
            raise GeneratorsNeeded("can't initialize from 'dict' without generators")

        level = len(gens)-1
        domain = opt.domain

        if domain is None:
            domain, rep = construct_domain(rep, opt=opt)
        else:
            for monom, coeff in rep.iteritems():
                rep[monom] = domain.convert(coeff)

        return cls.new(DMP.from_dict(rep, level, domain), *gens)

    @classmethod
    def _from_list(cls, rep, opt):
        """Construct a polynomial from a ``list``. """
        gens = opt.gens

        if not gens:
            raise GeneratorsNeeded("can't initialize from 'list' without generators")
        elif len(gens) != 1:
            raise MultivariatePolynomialError("'list' representation not supported")

        level = len(gens)-1
        domain = opt.domain

        if domain is None:
            domain, rep = construct_domain(rep, opt=opt)
        else:
            rep = map(domain.convert, rep)

        return cls.new(DMP.from_list(rep, level, domain), *gens)

    @classmethod
    def _from_poly(cls, rep, opt):
        """Construct a polynomial from a polynomial. """
        gens = opt.gens
        order = opt.order
        field = opt.field
        domain = opt.domain

        if gens and rep.gens != gens:
            if set(rep.gens) != set(gens):
                return cls._from_expr(rep.as_expr(), opt)
            else:
                rep = rep.reorder(*gens)

        if 'order' in opt:
            rep = rep.set_order(order)

        if 'domain' in opt and domain:
            rep = rep.set_domain(domain)
        elif field is True:
            rep = rep.to_field()

        return rep

    @classmethod
    def _from_expr(cls, rep, opt):
        """Construct a polynomial from an expression. """
        rep, opt = _dict_from_expr(rep, opt)
        return cls._from_dict(rep, opt)

    def __getnewargs__(self):
        """Data used by pickling protocol version 2. """
        return (self.rep, self.gens)

    def _hashable_content(self):
        """Allow SymPy to hash Poly instances. """
        return (self.rep, self.gens)

    @property
    def args(self):
        """
        Don't mess up with the core.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).args
        [1 + x**2]

        """
        return [self.as_expr()]

    @property
    def gen(self):
        """
        Return the principal generator.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).gen
        x

        """
        return self.gens[0]

    @property
    def domain(self):
        """Get the ground domain of ``self``. """
        return self.get_domain()

    @property
    def zero(self):
        """Return zero polynomial with ``self``'s properties. """
        return self.new(self.rep.zero(self.rep.lev, self.rep.dom), *self.gens)

    @property
    def one(self):
        """Return one polynomial with ``self``'s properties. """
        return self.new(self.rep.one(self.rep.lev, self.rep.dom), *self.gens)

    @property
    def unit(f):
        """Return unit polynomial with ``self``'s properties. """
        return self.new(self.rep.unit(self.rep.lev, self.rep.dom), *self.gens)

    def unify(f, g):
        """
        Make ``f`` and ``g`` belong to the same domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f, g = Poly(x/2 + 1), Poly(2*x + 1)

        >>> f
        Poly(1/2*x + 1, x, domain='QQ')
        >>> g
        Poly(2*x + 1, x, domain='ZZ')

        >>> F, G = f.unify(g)

        >>> F
        Poly(1/2*x + 1, x, domain='QQ')
        >>> G
        Poly(2*x + 1, x, domain='QQ')

        """
        _, per, F, G = f._unify(g)
        return per(F), per(G)

    def _unify(f, g):
        g = sympify(g)

        if not g.is_Poly:
            try:
                return f.rep.dom, f.per, f.rep, f.rep.per(f.rep.dom.from_sympy(g))
            except CoercionFailed:
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
        else:
            raise UnificationFailed("can't unify %s with %s" % (f, g))

        def per(rep, dom=dom, gens=gens, remove=None):
            if remove is not None:
                gens = gens[:remove]+gens[remove+1:]

                if not gens:
                    return dom.to_sympy(rep)

            return Poly.new(rep, *gens)

        return dom, per, F, G

    def per(f, rep, gens=None, remove=None):
        """
        Create a Poly out of the given representation.

        **Examples**

        >>> from sympy import Poly, ZZ
        >>> from sympy.abc import x, y

        >>> from sympy.polys.polyclasses import DMP

        >>> a = Poly(x**2 + 1)

        >>> a.per(DMP([ZZ(1), ZZ(1)], ZZ), gens=[y])
        Poly(y + 1, y, domain='ZZ')

        """
        if gens is None:
            gens = f.gens

        if remove is not None:
            gens = gens[:remove]+gens[remove+1:]

            if not gens:
                return f.rep.dom.to_sympy(rep)

        return Poly.new(rep, *gens)

    def set_domain(f, domain):
        """Set the ground domain of ``f``. """
        opt = options.build_options(f.gens, {'domain': domain})
        return f.per(f.rep.convert(opt.domain))

    def get_domain(f):
        """Get the ground domain of ``f``. """
        return f.rep.dom

    def set_modulus(f, modulus):
        """
        Set the modulus of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(5*x**2 + 2*x - 1, x).set_modulus(2)
        Poly(x**2 + 1, x, modulus=2)

        """
        modulus = options.Modulus.preprocess(modulus)
        return f.set_domain(FF(modulus))

    def get_modulus(f):
        """
        Get the modulus of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, modulus=2).get_modulus()
        2

        """
        domain = f.get_domain()

        if not domain.has_CharacteristicZero:
            return Integer(domain.characteristic())
        else:
            raise PolynomialError("not a polynomial over a Galois field")

    def _eval_subs(f, old, new):
        """Internal implementation of :func:`subs`. """
        if old in f.gens:
            if new.is_number:
                return f.eval(old, new)
            else:
                try:
                    return f.replace(old, new)
                except PolynomialError:
                    pass

        return f.as_expr().subs(old, new)

    def replace(f, x, y=None):
        """
        Replace ``x`` with ``y`` in generators list.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 1, x).replace(x, y)
        Poly(y**2 + 1, y, domain='ZZ')

        """
        if y is None:
            if f.is_univariate:
                x, y = f.gen, x
            else:
                raise PolynomialError("syntax supported only in univariate case")

        if x == y:
            return f

        if x in f.gens and y not in f.gens:
            dom = f.get_domain()

            if not dom.is_Composite or y not in dom.gens:
                gens = list(f.gens)
                gens[gens.index(x)] = y
                return f.per(f.rep, gens=gens)

        raise PolynomialError("can't replace %s with %s in %s" % (x, y, f))

    def reorder(f, *gens, **args):
        """
        Efficiently apply new order of generators.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + x*y**2, x, y).reorder(y, x)
        Poly(y**2*x + x**2, y, x, domain='ZZ')

        """
        opt = options.Options((), args)

        if not gens:
            gens = _sort_gens(f.gens, opt=opt)
        elif set(f.gens) != set(gens):
            raise PolynomialError("generators list can differ only up to order of elements")

        rep = dict(zip(*_dict_reorder(f.rep.to_dict(), f.gens, gens)))

        return f.per(DMP(rep, f.rep.dom, len(gens)-1), gens=gens)

    def ltrim(f, gen):
        """
        Remove dummy generators from the "left" of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y, z

        >>> Poly(y**2 + y*z**2, x, y, z).ltrim(y)
        Poly(y**2 + y*z**2, y, z, domain='ZZ')

        """
        rep = f.as_dict(native=True)
        j = f._gen_to_level(gen)
        terms = {}

        for monom, coeff in rep.iteritems():
            monom = monom[j:]

            if monom not in terms:
                terms[monom] = coeff
            else:
                raise PolynomialError("can't left trim %s" % f)

        gens = f.gens[j:]

        return f.new(DMP.from_dict(terms, len(gens)-1, f.rep.dom), *gens)

    def has_only_gens(f, *gens):
        """
        Return ``True`` if ``Poly(f, *gens)`` retains ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y, z

        >>> Poly(x*y + 1, x, y, z).has_only_gens(x, y)
        True
        >>> Poly(x*y + z, x, y, z).has_only_gens(x, y)
        False

        """
        f_gens = list(f.gens)
        indices = set([])

        for gen in gens:
            try:
                index = f_gens.index(gen)
            except ValueError:
                raise GeneratorsError("%s doesn't have %s as generator" % (f, gen))
            else:
                indices.add(index)

        for monom in f.monoms():
            for i, elt in enumerate(monom):
                if i not in indices and elt:
                    return False

        return True

    def to_ring(f):
        """
        Make the ground domain a ring.

        **Examples**

        >>> from sympy import Poly, QQ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, domain=QQ).to_ring()
        Poly(x**2 + 1, x, domain='ZZ')

        """
        if hasattr(f.rep, 'to_ring'):
            result = f.rep.to_ring()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'to_ring')

        return f.per(result)

    def to_field(f):
        """
        Make the ground domain a field.

        **Examples**

        >>> from sympy import Poly, ZZ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x, domain=ZZ).to_field()
        Poly(x**2 + 1, x, domain='QQ')

        """
        if hasattr(f.rep, 'to_field'):
            result = f.rep.to_field()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'to_field')

        return f.per(result)

    def to_exact(f):
        """
        Make the ground domain exact.

        **Examples**

        >>> from sympy import Poly, RR
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1.0, x, domain=RR).to_exact()
        Poly(x**2 + 1, x, domain='QQ')

        """
        if hasattr(f.rep, 'to_exact'):
            result = f.rep.to_exact()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'to_exact')

        return f.per(result)

    def retract(f, field=None):
        """
        Recalculate the ground domain of a polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2 + 1, x, domain='QQ[y]')
        >>> f
        Poly(x**2 + 1, x, domain='QQ[y]')

        >>> f.retract()
        Poly(x**2 + 1, x, domain='ZZ')
        >>> f.retract(field=True)
        Poly(x**2 + 1, x, domain='QQ')

        """
        dom, rep = construct_domain(f.as_dict(), field=field)
        return f.from_dict(rep, f.gens, domain=dom)

    def slice(f, x, m, n=None):
        """Take a continuous subsequence of terms of ``f``. """
        if n is None:
            j, m, n = 0, x, m
        else:
            j = f._gen_to_level(x)

        m, n = int(m), int(n)

        if hasattr(f.rep, 'slice'):
            result = f.rep.slice(m, n, j)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'slice')

        return f.per(result)

    def coeffs(f, order=None):
        """
        Returns all non-zero coefficients from ``f`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x + 3, x).coeffs()
        [1, 2, 3]

        """
        return [ f.rep.dom.to_sympy(c) for c in f.rep.coeffs(order=order) ]

    def monoms(f, order=None):
        """
        Returns all non-zero monomials from ``f`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).monoms()
        [(2, 0), (1, 2), (1, 1), (0, 1)]

        """
        return f.rep.monoms(order=order)

    def terms(f, order=None):
        """
        Returns all non-zero terms from ``f`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).terms()
        [((2, 0), 1), ((1, 2), 2), ((1, 1), 1), ((0, 1), 3)]

        """
        return [ (m, f.rep.dom.to_sympy(c)) for m, c in f.rep.terms(order=order) ]

    def all_coeffs(f):
        """
        Returns all coefficients from a univariate polynomial ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_coeffs()
        [1, 0, 2, -1]

        """
        return [ f.rep.dom.to_sympy(c) for c in f.rep.all_coeffs() ]

    def all_monoms(f):
        """
        Returns all monomials from a univariate polynomial ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_monoms()
        [(3,), (2,), (1,), (0,)]

        """
        return f.rep.all_monoms()

    def all_terms(f):
        """
        Returns all terms from a univariate polynomial ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_terms()
        [((3,), 1), ((2,), 0), ((1,), 2), ((0,), -1)]

        """
        return [ (m, f.rep.dom.to_sympy(c)) for m, c in f.rep.all_terms() ]

    def termwise(f, func, *gens, **args):
        """
        Apply a function to all terms of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> def func((k,), coeff):
        ...     return coeff//10**(2-k)

        >>> Poly(x**2 + 20*x + 400).termwise(func)
        Poly(x**2 + 2*x + 4, x, domain='ZZ')

        """
        terms = {}

        for monom, coeff in f.terms():
            result = func(monom, coeff)

            if isinstance(result, tuple):
                monom, coeff = result
            else:
                coeff = result

            if coeff:
                if monom not in terms:
                    terms[monom] = coeff
                else:
                    raise PolynomialError("%s monomial was generated twice" % monom)

        return f.from_dict(terms, *(gens or f.gens), **args)

    def length(f):
        """
        Returns the number of non-zero terms in ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 2*x - 1).length()
        3

        """
        return len(f.as_dict())

    def as_dict(f, native=False):
        """
        Switch to a ``dict`` representation.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 - y, x, y).as_dict()
        {(0, 1): -1, (1, 2): 2, (2, 0): 1}

        """
        if native:
            return f.rep.to_dict()
        else:
            return f.rep.to_sympy_dict()

    def as_list(f, native=False):
        """Switch to a ``list`` representation. """
        if native:
            return f.rep.to_list()
        else:
            return f.rep.to_sympy_list()

    def as_expr(f, *gens):
        """
        Convert a polynomial an expression.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2 + 2*x*y**2 - y, x, y)

        >>> f.as_expr()
        -y + x**2 + 2*x*y**2
        >>> f.as_expr({x: 5})
        25 - y + 10*y**2
        >>> f.as_expr(5, 6)
        379

        """
        if not gens:
            gens = f.gens
        elif len(gens) == 1 and isinstance(gens[0], dict):
            mapping = gens[0]
            gens = list(f.gens)

            for gen, value in mapping.iteritems():
                try:
                    index = gens.index(gen)
                except ValueError:
                    raise GeneratorsError("%s doesn't have %s as generator" % (f, gen))
                else:
                    gens[index] = value

        return basic_from_dict(f.rep.to_sympy_dict(), *gens)

    def lift(f):
        """
        Convert algebraic coefficients to rationals.

        **Examples**

        >>> from sympy import Poly, I
        >>> from sympy.abc import x

        >>> Poly(x**2 + I*x + 1, x, extension=I).lift()
        Poly(x**4 + 3*x**2 + 1, x, domain='QQ')

        """
        if hasattr(f.rep, 'lift'):
            result = f.rep.lift()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'lift')

        return f.per(result)

    def deflate(f):
        """
        Reduce degree of ``f`` by mapping ``x_i**m`` to ``y_i``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**6*y**2 + x**3 + 1, x, y).deflate()
        ((3, 2), Poly(x**2*y + x + 1, x, y, domain='ZZ'))

        """
        if hasattr(f.rep, 'deflate'):
            J, result = f.rep.deflate()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'deflate')

        return J, f.per(result)

    def inject(f, front=False):
        """
        Inject ground domain generators into ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x)

        >>> f.inject()
        Poly(x**2*y + x*y**3 + x*y + 1, x, y, domain='ZZ')
        >>> f.inject(front=True)
        Poly(y**3*x + y*x**2 + y*x + 1, y, x, domain='ZZ')

        """
        dom = f.rep.dom

        if dom.is_Numerical:
            return f
        elif not dom.is_Poly:
            raise DomainError("can't inject generators over %s" % dom)

        if hasattr(f.rep, 'inject'):
            result = f.rep.inject(front=front)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'inject')

        if front:
            gens = dom.gens + f.gens
        else:
            gens = f.gens + dom.gens

        return f.new(result, *gens)

    def eject(f, *gens):
        """
        Eject selected generators into the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x, y)

        >>> f.eject(x)
        Poly(x*y**3 + (x + x**2)*y + 1, y, domain='ZZ[x]')
        >>> f.eject(y)
        Poly(y*x**2 + (y + y**3)*x + 1, x, domain='ZZ[y]')

        """
        dom = f.rep.dom

        if not dom.is_Numerical:
            raise DomainError("can't eject generators over %s" % dom)

        n, k = len(f.gens), len(gens)

        if f.gens[:k] == gens:
            _gens, front = f.gens[n-k:], True
        elif f.gens[-k:] == gens:
            _gens, front = f.gens[:n-k], False
        else:
            raise NotImplementedError("can only eject front or back generators")

        dom = dom.inject(*gens)

        if hasattr(f.rep, 'eject'):
            result = f.rep.eject(dom, front=front)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'eject')

        return f.new(result, *_gens)

    def terms_gcd(f):
        """
        Remove GCD of terms from the polynomial ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**6*y**2 + x**3*y, x, y).terms_gcd()
        ((3, 1), Poly(x**3*y + 1, x, y, domain='ZZ'))

        """
        if hasattr(f.rep, 'terms_gcd'):
            J, result = f.rep.terms_gcd()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'terms_gcd')

        return J, f.per(result)

    def add_ground(f, coeff):
        """
        Add an element of the ground domain to ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).add_ground(2)
        Poly(x + 3, x, domain='ZZ')

        """
        if hasattr(f.rep, 'add_ground'):
            result = f.rep.add_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'add_ground')

        return f.per(result)

    def sub_ground(f, coeff):
        """
        Subtract an element of the ground domain from ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).sub_ground(2)
        Poly(x - 1, x, domain='ZZ')

        """
        if hasattr(f.rep, 'sub_ground'):
            result = f.rep.sub_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sub_ground')

        return f.per(result)

    def mul_ground(f, coeff):
        """
        Multiply ``f`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).mul_ground(2)
        Poly(2*x + 2, x, domain='ZZ')

        """
        if hasattr(f.rep, 'mul_ground'):
            result = f.rep.mul_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'mul_ground')

        return f.per(result)

    def quo_ground(f, coeff):
        """
        Quotient of ``f`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x + 4).quo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).quo_ground(2)
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2 does not divide 3 in ZZ

        """
        if hasattr(f.rep, 'quo_ground'):
            result = f.rep.quo_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'quo_ground')

        return f.per(result)

    def exquo_ground(f, coeff):
        """
        Exact quotient of ``f`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x + 4).exquo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).exquo_ground(2)
        Poly(x + 1, x, domain='ZZ')

        """
        if hasattr(f.rep, 'exquo_ground'):
            result = f.rep.exquo_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'exquo_ground')

        return f.per(result)

    def abs(f):
        """
        Make all coefficients in ``f`` positive.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).abs()
        Poly(x**2 + 1, x, domain='ZZ')

        """
        if hasattr(f.rep, 'abs'):
            result = f.rep.abs()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'abs')

        return f.per(result)

    def neg(f):
        """
        Negate all coefficients in ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).neg()
        Poly(-x**2 + 1, x, domain='ZZ')

        >>> -Poly(x**2 - 1, x)
        Poly(-x**2 + 1, x, domain='ZZ')

        """
        if hasattr(f.rep, 'neg'):
            result = f.rep.neg()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'neg')

        return f.per(result)

    def add(f, g):
        """
        Add two polynomials ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).add(Poly(x - 2, x))
        Poly(x**2 + x - 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x) + Poly(x - 2, x)
        Poly(x**2 + x - 1, x, domain='ZZ')

        """
        g = sympify(g)

        if not g.is_Poly:
            return f.add_ground(g)

        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'add'):
            result = F.add(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'add')

        return per(result)

    def sub(f, g):
        """
        Subtract two polynomials ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).sub(Poly(x - 2, x))
        Poly(x**2 - x + 3, x, domain='ZZ')

        >>> Poly(x**2 + 1, x) - Poly(x - 2, x)
        Poly(x**2 - x + 3, x, domain='ZZ')

        """
        g = sympify(g)

        if not g.is_Poly:
            return f.sub_ground(g)

        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'sub'):
            result = F.sub(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sub')

        return per(result)

    def mul(f, g):
        """
        Multiply two polynomials ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).mul(Poly(x - 2, x))
        Poly(x**3 - 2*x**2 + x - 2, x, domain='ZZ')

        >>> Poly(x**2 + 1, x)*Poly(x - 2, x)
        Poly(x**3 - 2*x**2 + x - 2, x, domain='ZZ')

        """
        g = sympify(g)

        if not g.is_Poly:
            return f.mul_ground(g)

        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'mul'):
            result = F.mul(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'mul')

        return per(result)

    def sqr(f):
        """
        Square a polynomial ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x - 2, x).sqr()
        Poly(x**2 - 4*x + 4, x, domain='ZZ')

        >>> Poly(x - 2, x)**2
        Poly(x**2 - 4*x + 4, x, domain='ZZ')

        """
        if hasattr(f.rep, 'sqr'):
            result = f.rep.sqr()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sqr')

        return f.per(result)

    def pow(f, n):
        """
        Raise ``f`` to a non-negative power ``n``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x - 2, x).pow(3)
        Poly(x**3 - 6*x**2 + 12*x - 8, x, domain='ZZ')

        >>> Poly(x - 2, x)**3
        Poly(x**3 - 6*x**2 + 12*x - 8, x, domain='ZZ')

        """
        n = int(n)

        if hasattr(f.rep, 'pow'):
            result = f.rep.pow(n)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'pow')

        return f.per(result)

    def pdiv(f, g):
        """
        Polynomial pseudo-division of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).pdiv(Poly(2*x - 4, x))
        (Poly(2*x + 4, x, domain='ZZ'), Poly(20, x, domain='ZZ'))

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'pdiv'):
            q, r = F.pdiv(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'pdiv')

        return per(q), per(r)

    def prem(f, g):
        """
        Polynomial pseudo-remainder of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).prem(Poly(2*x - 4, x))
        Poly(20, x, domain='ZZ')

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'prem'):
            result = F.prem(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'prem')

        return per(result)

    def pquo(f, g):
        """
        Polynomial pseudo-quotient of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).pquo(Poly(2*x - 2, x))
        Poly(2*x + 2, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).pquo(Poly(2*x - 4, x))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: -4 + 2*x does not divide 1 + x**2

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'pquo'):
            try:
                result = F.pquo(G)
            except ExactQuotientFailed, exc:
                raise exc.new(f.as_expr(), g.as_expr())
        else: # pragma: no cover
            raise OperationNotSupported(f, 'pquo')

        return per(result)

    def pexquo(f, g):
        """
        Polynomial exact pseudo-quotient of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).pexquo(Poly(2*x - 4, x))
        Poly(2*x + 4, x, domain='ZZ')

        >>> Poly(x**2 - 1, x).pexquo(Poly(2*x - 2, x))
        Poly(2*x + 2, x, domain='ZZ')

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'pexquo'):
            result = F.pexquo(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'pexquo')

        return per(result)

    def div(f, g, auto=False):
        """
        Polynomial division with remainder of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly, QQ, ZZ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x, domain=ZZ).div(Poly(2*x - 4, x, domain=ZZ))
        (Poly(0, x, domain='ZZ'), Poly(x**2 + 1, x, domain='ZZ'))

        >>> Poly(x**2 + 1, x, domain=QQ).div(Poly(2*x - 4, x, domain=QQ))
        (Poly(1/2*x + 1, x, domain='QQ'), Poly(5, x, domain='QQ'))

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'div'):
            q, r = F.div(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'div')

        return per(q), per(r)

    def rem(f, g, auto=False):
        """
        Computes the polynomial remainder of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly, ZZ, QQ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x, domain=ZZ).rem(Poly(2*x - 4, x, domain=ZZ))
        Poly(x**2 + 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x, domain=QQ).rem(Poly(2*x - 4, x, domain=QQ))
        Poly(5, x, domain='QQ')

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'rem'):
            result = F.rem(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'rem')

        return per(result)

    def quo(f, g, auto=False):
        """
        Computes polynomial quotient of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).quo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).quo(Poly(2*x - 4, x))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: -4 + 2*x does not divide 1 + x**2

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'quo'):
            try:
                result = F.quo(G)
            except ExactQuotientFailed, exc:
                raise exc.new(f.as_expr(), g.as_expr())
        else: # pragma: no cover
            raise OperationNotSupported(f, 'quo')

        return per(result)

    def exquo(f, g, auto=False):
        """
        Computes polynomial exact quotient of ``f`` by ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).exquo(Poly(2*x - 4, x))
        Poly(0, x, domain='ZZ')

        >>> Poly(x**2 - 1, x).exquo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'exquo'):
            result = F.exquo(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'exquo')

        return per(result)

    def _gen_to_level(f, gen):
        """Returns level associated with the given generator. """
        if isinstance(gen, int):
            length = len(f.gens)

            if -length <= gen < length:
                if gen < 0:
                    return length + gen
                else:
                    return gen
            else:
                raise PolynomialError("-%s <= gen < %s expected, got %s" % (length, length, gen))
        else:
            try:
                return list(f.gens).index(sympify(gen))
            except ValueError:
                raise PolynomialError("a valid generator expected, got %s" % gen)

    def degree(f, gen=0):
        """
        Returns degree of ``f`` in ``x_j``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).degree()
        2
        >>> Poly(x**2 + y*x + y, x, y).degree(y)
        1

        """
        j = f._gen_to_level(gen)

        if hasattr(f.rep, 'degree'):
            return f.rep.degree(j)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'degree')

    def degree_list(f):
        """
        Returns a list of degrees of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).degree_list()
        (2, 1)

        """
        if hasattr(f.rep, 'degree_list'):
            return f.rep.degree_list()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'degree_list')

    def total_degree(f):
        """
        Returns the total degree of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).total_degree()
        3

        """
        if hasattr(f.rep, 'total_degree'):
            return f.rep.total_degree()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'total_degree')

    def LC(f, order=None):
        """
        Returns the leading coefficient of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(4*x**3 + 2*x**2 + 3*x, x).LC()
        4

        """
        if order is not None:
            return f.coeffs(order)[0]

        if hasattr(f.rep, 'LC'):
            result = f.rep.LC()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'LC')

        return f.rep.dom.to_sympy(result)

    def TC(f):
        """
        Returns the trailing coefficent of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x**2 + 3*x, x).TC()
        0

        """
        if hasattr(f.rep, 'TC'):
            result = f.rep.TC()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'TC')

        return f.rep.dom.to_sympy(result)

    def EC(f, order=None):
        """
        Returns the last non-zero coefficient of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x**2 + 3*x, x).EC()
        3

        """
        if hasattr(f.rep, 'coeffs'):
            return f.coeffs(order)[-1]
        else: # pragma: no cover
            raise OperationNotSupported(f, 'EC')

    def nth(f, *N):
        """
        Returns the ``n``-th coefficient of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**3 + 2*x**2 + 3*x, x).nth(2)
        2
        >>> Poly(x**3 + 2*x*y**2 + y**2, x, y).nth(1, 2)
        2

        """
        if hasattr(f.rep, 'nth'):
            result = f.rep.nth(*map(int, N))
        else: # pragma: no cover
            raise OperationNotSupported(f, 'nth')

        return f.rep.dom.to_sympy(result)

    def LM(f, order=None):
        """
        Returns the leading monomial of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LM()
        (2, 0)

        """
        return f.monoms(order)[0]

    def EM(f, order=None):
        """
        Returns the last non-zero monomial of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).EM()
        (0, 1)

        """
        return f.monoms(order)[-1]

    def LT(f, order=None):
        """
        Returns the leading term of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LT()
        ((2, 0), 4)

        """
        return f.terms(order)[0]

    def ET(f, order=None):
        """
        Returns the last non-zero term of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).ET()
        ((0, 1), 3)

        """
        return f.terms(order)[-1]

    def max_norm(f):
        """
        Returns maximum norm of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(-x**2 + 2*x - 3, x).max_norm()
        3

        """
        if hasattr(f.rep, 'max_norm'):
            result = f.rep.max_norm()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'max_norm')

        return f.rep.dom.to_sympy(result)

    def l1_norm(f):
        """
        Returns l1 norm of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(-x**2 + 2*x - 3, x).l1_norm()
        6

        """
        if hasattr(f.rep, 'l1_norm'):
            result = f.rep.l1_norm()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'l1_norm')

        return f.rep.dom.to_sympy(result)

    def clear_denoms(f, convert=False):
        """
        Clear denominators, but keep the ground domain.

        **Examples**

        >>> from sympy import Poly, S, QQ
        >>> from sympy.abc import x

        >>> f = Poly(x/2 + S(1)/3, x, domain=QQ)

        >>> f.clear_denoms()
        (6, Poly(3*x + 2, x, domain='QQ'))
        >>> f.clear_denoms(convert=True)
        (6, Poly(3*x + 2, x, domain='ZZ'))

        """
        if not f.rep.dom.has_Field:
            return S.One, f

        dom = f.rep.dom.get_ring()

        if hasattr(f.rep, 'clear_denoms'):
            coeff, result = f.rep.clear_denoms()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'clear_denoms')

        coeff, f = dom.to_sympy(coeff), f.per(result)

        if not convert:
            return coeff, f
        else:
            return coeff, f.to_ring()

    def rat_clear_denoms(f, g):
        """
        Clear denominators in a rational function ``f/g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2/y + 1, x)
        >>> g = Poly(x**3 + y, x)

        >>> p, q = f.rat_clear_denoms(g)

        >>> p
        Poly(x**2 + y, x, domain='ZZ[y]')
        >>> q
        Poly(y*x**3 + y**2, x, domain='ZZ[y]')

        """
        dom, per, f, g = f._unify(g)

        f = per(f)
        g = per(g)

        if not (dom.has_Field and dom.has_assoc_Ring):
            return f, g

        a, f = f.clear_denoms(convert=True)
        b, g = g.clear_denoms(convert=True)

        f = f.mul_ground(b)
        g = g.mul_ground(a)

        return f, g

    def integrate(f, *specs, **args):
        """
        Computes indefinite integral of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x + 1, x).integrate()
        Poly(1/3*x**3 + x**2 + x, x, domain='QQ')

        >>> Poly(x*y**2 + x, x, y).integrate((0, 1), (1, 0))
        Poly(1/2*x**2*y**2 + 1/2*x**2, x, y, domain='QQ')

        """
        if args.get('auto', True) and f.rep.dom.has_Ring:
            f = f.to_field()

        if hasattr(f.rep, 'integrate'):
            if not specs:
                return f.per(f.rep.integrate(m=1))

            rep = f.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.integrate(int(m), f._gen_to_level(gen))

            return f.per(rep)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'integrate')

    def diff(f, *specs):
        """
        Computes partial derivative of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x + 1, x).diff()
        Poly(2*x + 2, x, domain='ZZ')

        >>> Poly(x*y**2 + x, x, y).diff((0, 0), (1, 1))
        Poly(2*x*y, x, y, domain='ZZ')

        """
        if hasattr(f.rep, 'diff'):
            if not specs:
                return f.per(f.rep.diff(m=1))

            rep = f.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.diff(int(m), f._gen_to_level(gen))

            return f.per(rep)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'diff')

    def eval(f, x, a=None, auto=True):
        """
        Evaluate ``f`` at ``a`` in the given variable.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y, z

        >>> Poly(x**2 + 2*x + 3, x).eval(2)
        11

        >>> Poly(2*x*y + 3*x + y + 2, x, y).eval(x, 2)
        Poly(5*y + 8, y, domain='ZZ')

        >>> f = Poly(2*x*y + 3*x + y + 2*z, x, y, z)

        >>> f.eval({x: 2})
        Poly(5*y + 2*z + 6, y, z, domain='ZZ')
        >>> f.eval({x: 2, y: 5})
        Poly(2*z + 31, z, domain='ZZ')
        >>> f.eval({x: 2, y: 5, z: 7})
        45

        """
        if a is None:
            if isinstance(x, (list, dict)):
                try:
                    mapping = x.items()
                except AttributeError:
                    mapping = x

                for gen, value in mapping:
                    f = f.eval(gen, value)

                return f
            else:
                j, a = 0, x
        else:
            j = f._gen_to_level(x)

        if not hasattr(f.rep, 'eval'): # pragma: no cover
            raise OperationNotSupported(f, 'eval')

        try:
            result = f.rep.eval(a, j)
        except CoercionFailed:
            if not auto:
                raise DomainError("can't evaluate at %s in %s" % (a, f.rep.dom))
            else:
                domain, [a] = construct_domain([a])
                f = f.set_domain(domain)
                result = f.rep.eval(a, j)

        return f.per(result, remove=j)

    def half_gcdex(f, g, auto=True):
        """
        Half extended Euclidean algorithm of ``f`` and ``g``.

        Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> Poly(f).half_gcdex(Poly(g))
        (Poly(-1/5*x + 3/5, x, domain='QQ'), Poly(x + 1, x, domain='QQ'))

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'half_gcdex'):
            s, h = F.half_gcdex(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'half_gcdex')

        return per(s), per(h)

    def gcdex(f, g, auto=True):
        """
        Extended Euclidean algorithm of ``f`` and ``g``.

        Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> Poly(f).gcdex(Poly(g))
        (Poly(-1/5*x + 3/5, x, domain='QQ'),
         Poly(1/5*x**2 - 6/5*x + 2, x, domain='QQ'),
         Poly(x + 1, x, domain='QQ'))

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'gcdex'):
            s, t, h = F.gcdex(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'gcdex')

        return per(s), per(t), per(h)

    def invert(f, g, auto=True):
        """
        Invert ``f`` modulo ``g`` when possible.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).invert(Poly(2*x - 1, x))
        Poly(-4/3, x, domain='QQ')

        >>> Poly(x**2 - 1, x).invert(Poly(x - 1, x))
        Traceback (most recent call last):
        ...
        NotInvertible: zero divisor

        """
        dom, per, F, G = f._unify(g)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(f.rep, 'invert'):
            result = F.invert(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'invert')

        return per(result)

    def revert(f, n):
        """Compute ``f**(-1)`` mod ``x**n``. """
        if hasattr(f.rep, 'revert'):
            result = f.rep.revert(int(n))
        else: # pragma: no cover
            raise OperationNotSupported(f, 'revert')

        return f.per(result)

    def subresultants(f, g):
        """
        Computes the subresultant PRS sequence of ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).subresultants(Poly(x**2 - 1, x))
        [Poly(x**2 + 1, x, domain='ZZ'),
         Poly(x**2 - 1, x, domain='ZZ'),
         Poly(-2, x, domain='ZZ')]

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'subresultants'):
            result = F.subresultants(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'subresultants')

        return map(per, result)

    def resultant(f, g):
        """
        Computes the resultant of ``f`` and ``g`` via PRS.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).resultant(Poly(x**2 - 1, x))
        4

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'resultant'):
            result = F.resultant(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'resultant')

        return per(result, remove=0)

    def discriminant(f):
        """
        Computes the discriminant of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 2*x + 3, x).discriminant()
        -8

        """
        if hasattr(f.rep, 'discriminant'):
            result = f.rep.discriminant()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'discriminant')

        return f.per(result, remove=0)

    def cofactors(f, g):
        """
        Returns the GCD of ``f`` and ``g`` and their cofactors.

        Returns polynomials ``(h, cff, cfg)`` such that ``h = gcd(f, g)``, and
        ``cff = quo(f, h)`` and ``cfg = quo(g, h)`` are, so called, cofactors
        of ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).cofactors(Poly(x**2 - 3*x + 2, x))
        (Poly(x - 1, x, domain='ZZ'),
         Poly(x + 1, x, domain='ZZ'),
         Poly(x - 2, x, domain='ZZ'))

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'cofactors'):
            h, cff, cfg = F.cofactors(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'cofactors')

        return per(h), per(cff), per(cfg)

    def gcd(f, g):
        """
        Returns the polynomial GCD of ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).gcd(Poly(x**2 - 3*x + 2, x))
        Poly(x - 1, x, domain='ZZ')

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'gcd'):
            result = F.gcd(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'gcd')

        return per(result)

    def lcm(f, g):
        """
        Returns polynomial LCM of ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).lcm(Poly(x**2 - 3*x + 2, x))
        Poly(x**3 - 2*x**2 - x + 2, x, domain='ZZ')

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'lcm'):
            result = F.lcm(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'lcm')

        return per(result)

    def trunc(f, p):
        """
        Reduce ``f`` modulo a constant ``p``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**3 + 3*x**2 + 5*x + 7, x).trunc(3)
        Poly(-x**3 - x + 1, x, domain='ZZ')

        """
        p = f.rep.dom.convert(p)

        if hasattr(f.rep, 'trunc'):
            result = f.rep.trunc(p)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'trunc')

        return f.per(result)

    def monic(f, auto=True):
        """
        Divides all coefficients by ``LC(f)``.

        **Examples**

        >>> from sympy import Poly, ZZ
        >>> from sympy.abc import x

        >>> Poly(3*x**2 + 6*x + 9, x, domain=ZZ).monic()
        Poly(x**2 + 2*x + 3, x, domain='QQ')

        >>> Poly(3*x**2 + 4*x + 2, x, domain=ZZ).monic()
        Poly(x**2 + 4/3*x + 2/3, x, domain='QQ')

        """
        if auto and f.rep.dom.has_Ring:
            f = f.to_field()

        if hasattr(f.rep, 'monic'):
            result = f.rep.monic()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'monic')

        return f.per(result)

    def content(f):
        """
        Returns the GCD of polynomial coefficients.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(6*x**2 + 8*x + 12, x).content()
        2

        """
        if hasattr(f.rep, 'content'):
            result = f.rep.content()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'content')

        return f.rep.dom.to_sympy(result)

    def primitive(f):
        """
        Returns the content and a primitive form of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 + 8*x + 12, x).primitive()
        (2, Poly(x**2 + 4*x + 6, x, domain='ZZ'))

        """
        if hasattr(f.rep, 'primitive'):
            cont, result = f.rep.primitive()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'primitive')

        return f.rep.dom.to_sympy(cont), f.per(result)

    def compose(f, g):
        """
        Computes the functional composition of ``f`` and ``g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + x, x).compose(Poly(x - 1, x))
        Poly(x**2 - x, x, domain='ZZ')

        """
        _, per, F, G = f._unify(g)

        if hasattr(f.rep, 'compose'):
            result = F.compose(G)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'compose')

        return per(result)

    def decompose(f):
        """
        Computes a functional decomposition of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**4 + 2*x**3 - x - 1, x, domain='ZZ').decompose()
        [Poly(x**2 - x - 1, x, domain='ZZ'), Poly(x**2 + x, x, domain='ZZ')]

        """
        if hasattr(f.rep, 'decompose'):
            result = f.rep.decompose()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'decompose')

        return map(f.per, result)

    def sturm(f, auto=True):
        """
        Computes the Sturm sequence of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 - 2*x**2 + x - 3, x).sturm()
        [Poly(x**3 - 2*x**2 + x - 3, x, domain='QQ'),
         Poly(3*x**2 - 4*x + 1, x, domain='QQ'),
         Poly(2/9*x + 25/9, x, domain='QQ'),
         Poly(-2079/4, x, domain='QQ')]

        """
        if auto and f.rep.dom.has_Ring:
            f = f.to_field()

        if hasattr(f.rep, 'sturm'):
            result = f.rep.sturm()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sturm')

        return map(f.per, result)

    def gff_list(f):
        """
        Computes greatest factorial factorization of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**5 + 2*x**4 - x**3 - 2*x**2

        >>> Poly(f).gff_list()
        [(Poly(x, x, domain='ZZ'), 1), (Poly(x + 2, x, domain='ZZ'), 4)]

        """
        if hasattr(f.rep, 'gff_list'):
            result = f.rep.gff_list()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'gff_list')

        return [ (f.per(g), k) for g, k in result ]

    def sqf_norm(f):
        """
        Computes square-free norm of ``f``.

        Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and
        ``r(x) = Norm(g(x))`` is a square-free polynomial over ``K``,
        where ``a`` is the algebraic extension of the ground domain.

        **Examples**

        >>> from sympy import Poly, sqrt
        >>> from sympy.abc import x

        >>> s, f, r = Poly(x**2 + 1, x, extension=[sqrt(3)]).sqf_norm()

        >>> s
        1
        >>> f
        Poly(x**2 - 2*3**(1/2)*x + 4, x, domain='QQ<3**(1/2)>')
        >>> r
        Poly(x**4 - 4*x**2 + 16, x, domain='QQ')

        """
        if hasattr(f.rep, 'sqf_norm'):
            s, g, r = f.rep.sqf_norm()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sqf_norm')

        return s, f.per(g), f.per(r)

    def sqf_part(f):
        """
        Computes square-free part of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 - 3*x - 2, x).sqf_part()
        Poly(x**2 - x - 2, x, domain='ZZ')

        """
        if hasattr(f.rep, 'sqf_part'):
            result = f.rep.sqf_part()
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sqf_part')

        return f.per(result)

    def sqf_list(f, all=False):
        """
        Returns a list of square-free factors of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

        >>> Poly(f).sqf_list()
        (2, [(Poly(x + 1, x, domain='ZZ'), 2),
             (Poly(x + 2, x, domain='ZZ'), 3)])

        >>> Poly(f).sqf_list(all=True)
        (2, [(Poly(1, x, domain='ZZ'), 1),
             (Poly(x + 1, x, domain='ZZ'), 2),
             (Poly(x + 2, x, domain='ZZ'), 3)])

        """
        if hasattr(f.rep, 'sqf_list'):
            coeff, factors = f.rep.sqf_list(all)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sqf_list')

        return f.rep.dom.to_sympy(coeff), [ (f.per(g), k) for g, k in factors ]

    def sqf_list_include(f, all=False):
        """
        Returns a list of square-free factors of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

        >>> Poly(f).sqf_list_include()
        [(Poly(2*x + 2, x, domain='ZZ'), 2),
         (Poly(x + 2, x, domain='ZZ'), 3)]

        >>> Poly(f).sqf_list_include(all=True)
        [(Poly(2, x, domain='ZZ'), 1),
         (Poly(x + 1, x, domain='ZZ'), 2),
         (Poly(x + 2, x, domain='ZZ'), 3)]

        """
        if hasattr(f.rep, 'sqf_list_include'):
            factors = f.rep.sqf_list_include(all)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'sqf_list_include')

        return [ (f.per(g), k) for g, k in factors ]

    def factor_list(f):
        """
        Returns a list of irreducible factors of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = 2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y

        >>> Poly(f).factor_list()
        (2, [(Poly(x + y, x, y, domain='ZZ'), 1),
             (Poly(x**2 + 1, x, y, domain='ZZ'), 2)])

        """
        if hasattr(f.rep, 'factor_list'):
            try:
                coeff, factors = f.rep.factor_list()
            except DomainError:
                return S.One, [(f, 1)]
        else: # pragma: no cover
            raise OperationNotSupported(f, 'factor_list')

        return f.rep.dom.to_sympy(coeff), [ (f.per(g), k) for g, k in factors ]

    def factor_list_include(f):
        """
        Returns a list of irreducible factors of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = 2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y

        >>> Poly(f).factor_list_include()
        [(Poly(2*x + 2*y, x, y, domain='ZZ'), 1),
         (Poly(x**2 + 1, x, y, domain='ZZ'), 2)]

        """
        if hasattr(f.rep, 'factor_list_include'):
            try:
                factors = f.rep.factor_list_include()
            except DomainError:
                return [(f, 1)]
        else: # pragma: no cover
            raise OperationNotSupported(f, 'factor_list_include')

        return [ (f.per(g), k) for g, k in factors ]

    def intervals(f, all=False, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """
        Compute isolating intervals for roots of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 3, x).intervals()
        [((-2, -1), 1), ((1, 2), 1)]
        >>> Poly(x**2 - 3, x).intervals(eps=1e-2)
        [((-26/15, -19/11), 1), ((19/11, 26/15), 1)]

        """
        if eps is not None:
            eps = QQ.convert(eps)

        if inf is not None:
            inf = QQ.convert(inf)
        if sup is not None:
            sup = QQ.convert(sup)

        if hasattr(f.rep, 'intervals'):
            result = f.rep.intervals(all=all, eps=eps, inf=inf, sup=sup, fast=fast, sqf=sqf)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'intervals')

        if sqf:
            def _real(interval):
                s, t = interval
                return (QQ.to_sympy(s), QQ.to_sympy(t))

            if not all:
                return map(_real, result)

            def _complex(rectangle):
                (u, v), (s, t) = rectangle
                return (QQ.to_sympy(u) + I*QQ.to_sympy(v),
                        QQ.to_sympy(s) + I*QQ.to_sympy(t))

            real_part, complex_part = result

            return map(_real, real_part), map(_complex, complex_part)
        else:
            def _real(interval):
                (s, t), k = interval
                return ((QQ.to_sympy(s), QQ.to_sympy(t)), k)

            if not all:
                return map(_real, result)

            def _complex(rectangle):
                ((u, v), (s, t)), k = rectangle
                return ((QQ.to_sympy(u) + I*QQ.to_sympy(v),
                         QQ.to_sympy(s) + I*QQ.to_sympy(t)), k)

            real_part, complex_part = result

            return map(_real, real_part), map(_complex, complex_part)

    def refine_root(f, s, t, eps=None, steps=None, fast=False, check_sqf=False):
        """
        Refine an isolating interval of a root to the given precision.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 3, x).refine_root(1, 2, eps=1e-2)
        (19/11, 26/15)

        """
        if check_sqf and not f.is_sqf:
            raise PolynomialError("only square-free polynomials supported")

        s, t = QQ.convert(s), QQ.convert(t)

        if eps is not None:
            eps = QQ.convert(eps)

        if steps is not None:
            steps = int(steps)
        elif eps is None:
            steps = 1

        if hasattr(f.rep, 'refine_root'):
            S, T = f.rep.refine_root(s, t, eps=eps, steps=steps, fast=fast)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'refine_root')

        return QQ.to_sympy(S), QQ.to_sympy(T)

    def count_roots(f, inf=None, sup=None):
        """
        Return the number of roots of ``f`` in ``[inf, sup]`` interval.

        **Examples**

        >>> from sympy import Poly, I
        >>> from sympy.abc import x

        >>> Poly(x**4 - 4, x).count_roots(-3, 3)
        2
        >>> Poly(x**4 - 4, x).count_roots(0, 1 + 3*I)
        1

        """
        inf_real, sup_real = True, True

        if inf is not None:
            inf = sympify(inf)

            if inf is S.NegativeInfinity:
                inf = None
            else:
                re, im = inf.as_real_imag()

                if not im:
                    inf = QQ.convert(inf)
                else:
                    inf, inf_real = map(QQ.convert, (re, im)), False

        if sup is not None:
            sup = sympify(sup)

            if sup is S.Infinity:
                sup = None
            else:
                re, im = sup.as_real_imag()

                if not im:
                    sup = QQ.convert(sup)
                else:
                    sup, sup_real = map(QQ.convert, (re, im)), False

        if inf_real and sup_real:
            if hasattr(f.rep, 'count_real_roots'):
                count = f.rep.count_real_roots(inf=inf, sup=sup)
            else: # pragma: no cover
                raise OperationNotSupported(f, 'count_real_roots')
        else:
            if inf_real and inf is not None:
                inf = (inf, QQ.zero)

            if sup_real and sup is not None:
                sup = (sup, QQ.zero)

            if hasattr(f.rep, 'count_complex_roots'):
                count = f.rep.count_complex_roots(inf=inf, sup=sup)
            else: # pragma: no cover
                raise OperationNotSupported(f, 'count_complex_roots')

        return Integer(count)

    def real_roots(f, multiple=True):
        """
        Return a list of real roots with multiplicities.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**3 - 7*x**2 + 4*x + 4, x).real_roots()
        [-1/2, 2, 2]

        """
        reals = sympy.polys.rootoftools.RootOf(f)

        if multiple:
            return reals
        else:
            return group(reals, multiple=False)

    def nroots(f, maxsteps=50, cleanup=True, error=False):
        """
        Compute numerical approximations of roots of ``f``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 3).nroots()
        [-1.73205080756888, 1.73205080756888]

        """
        if f.is_multivariate:
            raise MultivariatePolynomialError("can't compute numerical roots of %s" % f)

        if f.degree() <= 0:
            return []

        try:
            coeffs = [ complex(c) for c in f.all_coeffs() ]
        except ValueError:
            raise DomainError("numerical domain expected, got %s" % f.rep.dom)

        return sympify(npolyroots(coeffs, maxsteps=maxsteps, cleanup=cleanup, error=error))

    def ground_roots(f):
        """
        Compute roots of ``f`` by factorization in the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**6 - 4*x**4 + 4*x**3 - x**2).ground_roots()
        {1: 2, 0: 2}

        """
        if f.is_multivariate:
            raise MultivariatePolynomialError("can't compute ground roots of %s" % f)

        roots = {}

        for factor, k in f.factor_list()[1]:
            if factor.is_linear:
                a, b = factor.all_coeffs()
                roots[-b/a] = k

        return roots

    def cancel(f, g):
        """
        Cancel common factors in a rational function ``f/g``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 - 2, x).cancel(Poly(x**2 - 2*x + 1, x))
        (1, Poly(2*x + 2, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'))

        """
        dom, per, F, G = f._unify(g)

        if hasattr(F, 'cancel'):
            cp, cq, p, q = F.cancel(G, multout=False)
        else: # pragma: no cover
            raise OperationNotSupported(f, 'cancel')

        if dom.has_assoc_Ring:
            dom = dom.get_ring()

        cp = dom.to_sympy(cp)
        cq = dom.to_sympy(cq)

        return cp/cq, per(p), per(q)

    @property
    def is_zero(f):
        """
        Returns ``True`` if ``f`` is a zero polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(0, x).is_zero
        True
        >>> Poly(1, x).is_zero
        False

        """
        return f.rep.is_zero

    @property
    def is_one(f):
        """
        Returns ``True`` if ``f`` is a unit polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(0, x).is_one
        False
        >>> Poly(1, x).is_one
        True

        """
        return f.rep.is_one

    @property
    def is_sqf(f):
        """
        Returns ``True`` if ``f`` is a square-free polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 2*x + 1, x).is_sqf
        False
        >>> Poly(x**2 - 1, x).is_sqf
        True

        """
        return f.rep.is_sqf

    @property
    def is_monic(f):
        """
        Returns ``True`` if the leading coefficient of ``f`` is one.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 2, x).is_monic
        True
        >>> Poly(2*x + 2, x).is_monic
        False

        """
        return f.rep.is_monic

    @property
    def is_primitive(f):
        """
        Returns ``True`` if GCD of the coefficients of ``f`` is one.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 + 6*x + 12, x).is_primitive
        False
        >>> Poly(x**2 + 3*x + 6, x).is_primitive
        True

        """
        return f.rep.is_primitive

    @property
    def is_ground(f):
        """
        Returns ``True`` if ``f`` is an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x, x).is_ground
        False
        >>> Poly(2, x).is_ground
        True
        >>> Poly(y, x).is_ground
        True

        """
        return f.rep.is_ground

    @property
    def is_linear(f):
        """
        Returns ``True`` if ``f`` is linear in all its variables.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x + y + 2, x, y).is_linear
        True
        >>> Poly(x*y + 2, x, y).is_linear
        False

        """
        return f.rep.is_linear

    @property
    def is_quadratic(f):
        """
        Returns ``True`` if ``f`` is quadratic in all its variables.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x*y + 2, x, y).is_quadratic
        True
        >>> Poly(x*y**2 + 2, x, y).is_quadratic
        False

        """
        return f.rep.is_quadratic

    @property
    def is_monomial(f):
        """
        Returns ``True`` if ``f`` is zero or has only one term.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(3*x**2, x).is_monomial
        True
        >>> Poly(3*x**2 + 1, x).is_monomial
        False

        """
        return f.length() <= 1

    @property
    def is_homogeneous(f):
        """
        Returns ``True`` if ``f`` has zero trailing coefficient.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x*y + x + y, x, y).is_homogeneous
        True
        >>> Poly(x*y + x + y + 1, x, y).is_homogeneous
        False

        """
        return f.rep.is_homogeneous

    @property
    def is_irreducible(f):
        """
        Returns ``True`` if ``f`` has no factors over its domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >> Poly(x**2 + x + 1, x, modulus=2).is_irreducible
        True
        >> Poly(x**2 + 1, x, modulus=2).is_irreducible
        False

        """
        return f.rep.is_irreducible

    @property
    def is_univariate(f):
        """
        Returns ``True`` if ``f`` is an univariate polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + x + 1, x).is_univariate
        True
        >>> Poly(x*y**2 + x*y + 1, x, y).is_univariate
        False
        >>> Poly(x*y**2 + x*y + 1, x).is_univariate
        True
        >>> Poly(x**2 + x + 1, x, y).is_univariate
        False

        """
        return len(f.gens) == 1

    @property
    def is_multivariate(f):
        """
        Returns ``True`` if ``f`` is a multivariate polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + x + 1, x).is_multivariate
        False
        >>> Poly(x*y**2 + x*y + 1, x, y).is_multivariate
        True
        >>> Poly(x*y**2 + x*y + 1, x).is_multivariate
        False
        >>> Poly(x**2 + x + 1, x, y).is_multivariate
        True

        """
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
                return f.as_expr() + g

        return f.add(g)

    @_sympifyit('g', NotImplemented)
    def __radd__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g + f.as_expr()

        return g.add(f)

    @_sympifyit('g', NotImplemented)
    def __sub__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_expr() - g

        return f.sub(g)

    @_sympifyit('g', NotImplemented)
    def __rsub__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g - f.as_expr()

        return g.sub(f)

    @_sympifyit('g', NotImplemented)
    def __mul__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return f.as_expr()*g

        return f.mul(g)

    @_sympifyit('g', NotImplemented)
    def __rmul__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens)
            except PolynomialError:
                return g*f.as_expr()

        return g.mul(f)

    @_sympifyit('n', NotImplemented)
    def __pow__(f, n):
        if n.is_Integer and n >= 0:
            return f.pow(n)
        else:
            return f.as_expr()**n

    @_sympifyit('g', NotImplemented)
    def __divmod__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return f.div(g)

    @_sympifyit('g', NotImplemented)
    def __rdivmod__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return g.div(f)

    @_sympifyit('g', NotImplemented)
    def __mod__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return f.rem(g)

    @_sympifyit('g', NotImplemented)
    def __rmod__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return g.rem(f)

    @_sympifyit('g', NotImplemented)
    def __floordiv__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return f.exquo(g)

    @_sympifyit('g', NotImplemented)
    def __rfloordiv__(f, g):
        if not g.is_Poly:
            g = Poly(g, *f.gens)

        return g.exquo(f)

    @_sympifyit('g', NotImplemented)
    def __div__(f, g):
        return f.as_expr()/g.as_expr()

    @_sympifyit('g', NotImplemented)
    def __rdiv__(f, g):
        return g.as_expr()/f.as_expr()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    @_sympifyit('g', NotImplemented)
    def __eq__(f, g):
        if not g.is_Poly:
            try:
                g = Poly(g, *f.gens, **{'domain': f.get_domain()})
            except (PolynomialError, DomainError, CoercionFailed):
                return False

        if f.gens != g.gens:
            return False

        if f.rep.dom != g.rep.dom:
            try:
                dom = f.rep.dom.unify(g.rep.dom, f.gens)
            except UnificationFailed:
                return False

            f = f.set_domain(dom)
            g = g.set_domain(dom)

        return f.rep == g.rep

    @_sympifyit('g', NotImplemented)
    def __ne__(f, g):
        return not f.__eq__(g)

    def __hash__(self):
        return super(Poly, self).__hash__()

    def __nonzero__(f):
        return not f.is_zero

def poly_from_expr(expr, *gens, **args):
    """Construct a polynomial from an expression. """
    opt = options.build_options(gens, args)
    return _poly_from_expr(expr, opt)

def _poly_from_expr(expr, opt):
    """Construct a polynomial from an expression. """
    orig, expr = expr, sympify(expr)

    if not isinstance(expr, Basic):
        raise PolificationFailed(opt, orig, expr)
    elif expr.is_Poly:
        poly = Poly(expr, opt=opt)

        opt['gens'] = poly.gens
        opt['domain'] = poly.domain

        if opt.polys is None:
            opt['polys'] = True

        return poly, opt
    elif opt.expand:
        expr = expr.expand()

    try:
        rep, opt = _dict_from_expr(expr, opt)
    except GeneratorsNeeded:
        raise PolificationFailed(opt, orig, expr)

    monoms, coeffs = zip(*rep.items())
    domain = opt.domain

    if domain is None:
        domain, coeffs = construct_domain(coeffs, opt=opt)
    else:
        coeffs = map(domain.from_sympy, coeffs)

    level = len(opt.gens)-1

    poly = Poly.new(DMP.from_monoms_coeffs(monoms, coeffs, level, domain), *opt.gens)

    opt['domain'] = domain

    if opt.polys is None:
        opt['polys'] = False

    return poly, opt

def parallel_poly_from_expr(exprs, *gens, **args):
    """Construct polynomials from expressions. """
    opt = options.build_options(gens, args)
    return _parallel_poly_from_expr(exprs, opt)

def _parallel_poly_from_expr(exprs, opt):
    """Construct polynomials from expressions. """
    if len(exprs) == 2:
        f, g = exprs

        if isinstance(f, Poly) and isinstance(g, Poly):
            f = Poly._from_poly(f, opt)
            g = Poly._from_poly(g, opt)

            f, g = f.unify(g)

            opt['gens'] = f.gens
            opt['domain'] = f.domain

            if opt.polys is None:
                opt['polys'] = True

            return [f, g], opt

    origs, exprs = list(exprs), []
    _exprs, _polys = [], []

    failed = False

    for i, expr in enumerate(origs):
        expr = sympify(expr)

        if isinstance(expr, Basic):
            if expr.is_Poly:
                _polys.append(i)
            else:
                _exprs.append(i)

                if opt.expand:
                    expr = expr.expand()
        else:
            failed = True

        exprs.append(expr)

    if failed:
        raise PolificationFailed(opt, origs, exprs, True)

    if _polys:
        # XXX: this is a temporary solution
        for i in _polys:
            exprs[i] = exprs[i].as_expr()

    try:
        reps, opt = _parallel_dict_from_expr(exprs, opt)
    except GeneratorsNeeded:
        raise PolificationFailed(opt, origs, exprs, True)

    coeffs_list, lengths = [], []

    all_monoms = []
    all_coeffs = []

    for rep in reps:
        monoms, coeffs = zip(*rep.items())

        coeffs_list.extend(coeffs)
        all_monoms.append(monoms)

        lengths.append(len(coeffs))

    domain = opt.domain

    if domain is None:
        domain, coeffs_list = construct_domain(coeffs_list, opt=opt)
    else:
        coeffs_list = map(domain.from_sympy, coeffs_list)

    for k in lengths:
        all_coeffs.append(coeffs_list[:k])
        coeffs_list = coeffs_list[k:]

    polys, level = [], len(opt.gens)-1

    for monoms, coeffs in zip(all_monoms, all_coeffs):
        rep = DMP.from_monoms_coeffs(monoms, coeffs, level, domain)
        polys.append(Poly.new(rep, *opt.gens))

    opt['domain'] = domain

    if opt.polys is None:
        opt['polys'] = bool(_polys)

    return polys, opt

def _update_args(args, key, value):
    """Add a new ``(key, value)`` pair to arguments ``dict``. """
    args = dict(args)

    if key not in args:
        args[key] = value

    return args

def _keep_coeff(coeff, factors):
    """Return ``coeff*factors`` unevaluated if necessary. """
    if coeff == 1:
        return factors
    elif coeff == -1:
        return -factors
    elif not factors.is_Add:
        return coeff*factors
    else:
        return Mul(coeff, factors, evaluate=False)

def degree(f, *gens, **args):
    """
    Return the degree of ``f`` in the given variable.

    **Examples**

    >>> from sympy import degree
    >>> from sympy.abc import x, y

    >>> degree(x**2 + y*x + 1, gen=x)
    2
    >>> degree(x**2 + y*x + 1, gen=y)
    1

    """
    options.allowed_flags(args, ['gen', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('degree', 1, exc)

    return Integer(F.degree(opt.gen))

def degree_list(f, *gens, **args):
    """
    Return a list of degrees of ``f`` in all variables.

    **Examples**

    >>> from sympy import degree_list
    >>> from sympy.abc import x, y

    >>> degree_list(x**2 + y*x + 1)
    (2, 1)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('degree_list', 1, exc)

    degrees = F.degree_list()

    return tuple(map(Integer, degrees))

def LC(f, *gens, **args):
    """
    Return the leading coefficient of ``f``.

    **Examples**

    >>> from sympy import LC
    >>> from sympy.abc import x, y

    >>> LC(4*x**2 + 2*x*y**2 + x*y + 3*y)
    4

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('LC', 1, exc)

    return F.LC(order=opt.order)

def LM(f, *gens, **args):
    """
    Return the leading monomial of ``f``.

    **Examples**

    >>> from sympy import LM
    >>> from sympy.abc import x, y

    >>> LM(4*x**2 + 2*x*y**2 + x*y + 3*y)
    x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('LM', 1, exc)

    monom = Monomial(*F.LM(order=opt.order))

    return monom.as_expr(*opt.gens)

def LT(f, *gens, **args):
    """
    Return the leading term of ``f``.

    **Examples**

    >>> from sympy import LT
    >>> from sympy.abc import x, y

    >>> LT(4*x**2 + 2*x*y**2 + x*y + 3*y)
    4*x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('LT', 1, exc)

    monom, coeff = F.LT(order=opt.order)

    return coeff*Monomial(*monom).as_expr(*opt.gens)

def pdiv(f, g, *gens, **args):
    """
    Compute polynomial pseudo-division of ``f`` and ``g``.

    **Examples**

    >>> from sympy import pdiv
    >>> from sympy.abc import x

    >>> pdiv(x**2 + 1, 2*x - 4)
    (4 + 2*x, 20)

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('pdiv', 2, exc)

    q, r = F.pdiv(G)

    if not opt.polys:
        return q.as_expr(), r.as_expr()
    else:
        return q, r

def prem(f, g, *gens, **args):
    """
    Compute polynomial pseudo-remainder of ``f`` and ``g``.

    **Examples**

    >>> from sympy import prem
    >>> from sympy.abc import x

    >>> prem(x**2 + 1, 2*x - 4)
    20

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('prem', 2, exc)

    r = F.prem(G)

    if not opt.polys:
        return r.as_expr()
    else:
        return r

def pquo(f, g, *gens, **args):
    """
    Compute polynomial pseudo-quotient of ``f`` and ``g``.

    **Examples**

    >>> from sympy import pquo
    >>> from sympy.abc import x

    >>> pquo(x**2 - 1, 2*x - 2)
    2 + 2*x

    >>> pquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: -4 + 2*x does not divide 1 + x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('pquo', 2, exc)

    q = F.pquo(G)

    if not opt.polys:
        return q.as_expr()
    else:
        return q

def pexquo(f, g, *gens, **args):
    """
    Compute polynomial exact pseudo-quotient of ``f`` and ``g``.

    **Examples**

    >>> from sympy import pexquo
    >>> from sympy.abc import x

    >>> pexquo(x**2 + 1, 2*x - 4)
    4 + 2*x
    >>> pexquo(x**2 - 1, 2*x - 1)
    1 + 2*x

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('pexquo', 2, exc)

    q = F.pexquo(G)

    if not opt.polys:
        return q.as_expr()
    else:
        return q

def div(f, g, *gens, **args):
    """
    Compute polynomial division of ``f`` and ``g``.

    **Examples**

    >>> from sympy import div, ZZ, QQ
    >>> from sympy.abc import x

    >>> div(x**2 + 1, 2*x - 4, domain=ZZ)
    (0, 1 + x**2)
    >>> div(x**2 + 1, 2*x - 4, domain=QQ)
    (1 + x/2, 5)

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('div', 2, exc)

    q, r = F.div(G)

    if not opt.polys:
        return q.as_expr(), r.as_expr()
    else:
        return q, r

def rem(f, g, *gens, **args):
    """
    Compute polynomial remainder of ``f`` and ``g``.

    **Examples**

    >>> from sympy import rem, ZZ, QQ
    >>> from sympy.abc import x

    >>> rem(x**2 + 1, 2*x - 4, domain=ZZ)
    1 + x**2
    >>> rem(x**2 + 1, 2*x - 4, domain=QQ)
    5

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('rem', 2, exc)

    r = F.rem(G)

    if not opt.polys:
        return r.as_expr()
    else:
        return r

def quo(f, g, *gens, **args):
    """
    Compute polynomial quotient of ``f`` and ``g``.

    **Examples**

    >>> from sympy import quo
    >>> from sympy.abc import x

    >>> quo(x**2 - 1, x - 1)
    1 + x

    >>> quo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: -4 + 2*x does not divide 1 + x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('quo', 2, exc)

    q = F.quo(G)

    if not opt.polys:
        return q.as_expr()
    else:
        return q

def exquo(f, g, *gens, **args):
    """
    Compute polynomial exact quotient of ``f`` and ``g``.

    **Examples**

    >>> from sympy import exquo
    >>> from sympy.abc import x

    >>> exquo(x**2 + 1, 2*x - 4)
    0
    >>> exquo(x**2 - 1, x - 1)
    1 + x

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('exquo', 2, exc)

    q = F.exquo(G)

    if not opt.polys:
        return q.as_expr()
    else:
        return q

def half_gcdex(f, g, *gens, **args):
    """
    Half extended Euclidean algorithm of ``f`` and ``g``.

    Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

    **Examples**

    >>> from sympy import half_gcdex
    >>> from sympy.abc import x

    >>> half_gcdex(x**4 - 2*x**3 - 6*x**2 + 12*x + 15, x**3 + x**2 - 4*x - 4)
    (3/5 - x/5, 1 + x)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'half_gcdex'):
            try:
                return f.half_gcdex(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('half_gcdex', 2, exc)

    s, h = F.half_gcdex(G, auto=opt.auto)

    if not opt.polys:
        return s.as_expr(), h.as_expr()
    else:
        return s, h

def gcdex(f, g, *gens, **args):
    """
    Extended Euclidean algorithm of ``f`` and ``g``.

    Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

    **Examples**

    >>> from sympy import gcdex
    >>> from sympy.abc import x

    >>> gcdex(x**4 - 2*x**3 - 6*x**2 + 12*x + 15, x**3 + x**2 - 4*x - 4)
    (3/5 - x/5, 2 - 6*x/5 + x**2/5, 1 + x)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'gcdex'):
            try:
                return f.gcdex(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('gcdex', 2, exc)

    s, t, h = F.gcdex(G, auto=opt.auto)

    if not opt.polys:
        return s.as_expr(), t.as_expr(), h.as_expr()
    else:
        return s, t, h

def invert(f, g, *gens, **args):
    """
    Invert ``f`` modulo ``g`` when possible.

    **Examples**

    >>> from sympy import invert
    >>> from sympy.abc import x

    >>> invert(x**2 - 1, 2*x - 1)
    -4/3

    >>> invert(x**2 - 1, x - 1)
    Traceback (most recent call last):
    ...
    NotInvertible: zero divisor

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'invert'):
            try:
                return f.invert(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('invert', 2, exc)

    h = F.invert(G, auto=opt.auto)

    if not opt.polys:
        return h.as_expr()
    else:
        return h

def subresultants(f, g, *gens, **args):
    """
    Compute subresultant PRS of ``f`` and ``g``.

    **Examples**

    >>> from sympy import subresultants
    >>> from sympy.abc import x

    >>> subresultants(x**2 + 1, x**2 - 1)
    [1 + x**2, -1 + x**2, -2]

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('subresultants', 2, exc)

    result = F.subresultants(G)

    if not opt.polys:
        return [ r.as_expr() for r in result ]
    else:
        return result

def resultant(f, g, *gens, **args):
    """
    Compute resultant of ``f`` and ``g``.

    **Examples**

    >>> from sympy import resultant
    >>> from sympy.abc import x

    >>> resultant(x**2 + 1, x**2 - 1)
    4

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('resultant', 2, exc)

    result = F.resultant(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def discriminant(f, *gens, **args):
    """
    Compute discriminant of ``f``.

    **Examples**

    >>> from sympy import discriminant
    >>> from sympy.abc import x

    >>> discriminant(x**2 + 2*x + 3)
    -8

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('discriminant', 1, exc)

    result = F.discriminant()

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def cofactors(f, g, *gens, **args):
    """
    Compute GCD and cofactors of ``f`` and ``g``.

    Returns polynomials ``(h, cff, cfg)`` such that ``h = gcd(f, g)``, and
    ``cff = quo(f, h)`` and ``cfg = quo(g, h)`` are, so called, cofactors
    of ``f`` and ``g``.

    **Examples**

    >>> from sympy import cofactors
    >>> from sympy.abc import x

    >>> cofactors(x**2 - 1, x**2 - 3*x + 2)
    (-1 + x, 1 + x, -2 + x)

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'cofactors'):
            try:
                return f.cofactors(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('cofactors', 2, exc)

    h, cff, cfg = F.cofactors(G)

    if not opt.polys:
        return h.as_expr(), cff.as_expr(), cfg.as_expr()
    else:
        return h, cff, cfg

def gcd_list(seq, *gens, **args):
    """
    Compute GCD of a list of polynomials.

    **Examples**

    >>> from sympy import gcd_list
    >>> from sympy.abc import x

    >>> gcd_list([x**3 - 1, x**2 - 1, x**2 - 3*x + 2])
    -1 + x

    """
    if not gens and not args:
        if not seq:
            return S.Zero

        seq = sympify(seq)

        if all(s.is_Number for s in seq):
            result, numbers = seq[0], seq[1:]

            for number in numbers:
                result = result.gcd(number)

                if result is S.One:
                    break

            return result

    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('gcd_list', len(seq), exc)

    if not polys:
        if not opt.polys:
            return S.Zero
        else:
            return Poly(0, opt=opt)

    result, polys = polys[0], polys[1:]

    for poly in polys:
        result = result.gcd(poly)

        if result.is_one:
            break

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def gcd(f, g=None, *gens, **args):
    """
    Compute GCD of ``f`` and ``g``.

    **Examples**

    >>> from sympy import gcd
    >>> from sympy.abc import x

    >>> gcd(x**2 - 1, x**2 - 3*x + 2)
    -1 + x

    """
    if hasattr(f, '__iter__'):
        if g is not None:
            gens = (g,) + gens

        return gcd_list(f, *gens, **args)

    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'gcd'):
            try:
                return f.gcd(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('gcd', 2, exc)

    result = F.gcd(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def lcm_list(seq, *gens, **args):
    """
    Compute LCM of a list of polynomials.

    **Examples**

    >>> from sympy import lcm_list
    >>> from sympy.abc import x

    >>> lcm_list([x**3 - 1, x**2 - 1, x**2 - 3*x + 2])
    2 + x - x**2 - 2*x**3 - x**4 + x**5

    """
    if not gens and not args:
        if not seq:
            return S.One

        seq = sympify(seq)

        if all(s.is_Number for s in seq):
            result, numbers = seq[0], seq[1:]

            for number in numbers:
                result = result.lcm(number)

            return result

    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('lcm_list', len(seq), exc)

    if not polys:
        if not opt.polys:
            return S.One
        else:
            return Poly(1, opt=opt)

    result, polys = polys[0], polys[1:]

    for poly in polys:
        result = result.lcm(poly)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def lcm(f, g=None, *gens, **args):
    """
    Compute LCM of ``f`` and ``g``.

    **Examples**

    >>> from sympy import lcm
    >>> from sympy.abc import x

    >>> lcm(x**2 - 1, x**2 - 3*x + 2)
    2 - x - 2*x**2 + x**3

    """
    if hasattr(f, '__iter__'):
        if g is not None:
            gens = (g,) + gens

        return lcm_list(f, *gens, **args)

    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        f, g = exc.exprs

        if hasattr(f, 'lcm'):
            try:
                return f.lcm(g)
            except (SympifyError, ValueError):
                pass

        raise ComputationFailed('lcm', 2, exc)

    result = F.lcm(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def terms_gcd(f, *gens, **args):
    """
    Remove GCD of terms from ``f``.

    **Examples**

    >>> from sympy import terms_gcd
    >>> from sympy.abc import x, y

    >>> terms_gcd(x**6*y**2 + x**3*y, x, y)
    y*x**3*(1 + y*x**3)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        return exc.expr

    J, f = F.terms_gcd()

    if opt.domain.has_Ring:
        if opt.domain.has_Field:
            denom, f = f.clear_denoms(convert=True)

        coeff, f = f.primitive()

        if opt.domain.has_Field:
            coeff /= denom
    else:
        coeff = S.One

    term = Mul(*[ x**j for x, j in zip(f.gens, J) ])

    return _keep_coeff(coeff, term*f.as_expr())

def trunc(f, p, *gens, **args):
    """
    Reduce ``f`` modulo a constant ``p``.

    **Examples**

    >>> from sympy import trunc
    >>> from sympy.abc import x

    >>> trunc(2*x**3 + 3*x**2 + 5*x + 7, 3)
    1 - x - x**3

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('trunc', 1, exc)

    result = F.trunc(sympify(p))

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def monic(f, *gens, **args):
    """
    Divide all coefficients of ``f`` by ``LC(f)``.

    **Examples**

    >>> from sympy import monic
    >>> from sympy.abc import x

    >>> monic(3*x**2 + 4*x + 2)
    2/3 + 4*x/3 + x**2

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('monic', 1, exc)

    result = F.monic(auto=opt.auto)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def content(f, *gens, **args):
    """
    Compute GCD of coefficients of ``f``.

    **Examples**

    >>> from sympy import content
    >>> from sympy.abc import x

    >>> content(6*x**2 + 8*x + 12)
    2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('content', 1, exc)

    return F.content()

def primitive(f, *gens, **args):
    """
    Compute content and the primitive form of ``f``.

    **Examples**

    >>> from sympy import primitive
    >>> from sympy.abc import x

    >>> primitive(6*x**2 + 8*x + 12)
    (2, 6 + 4*x + 3*x**2)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('primitive', 1, exc)

    cont, result = F.primitive()

    if not opt.polys:
        return cont, result.as_expr()
    else:
        return cont, result

def compose(f, g, *gens, **args):
    """
    Compute functional composition ``f(g)``.

    **Examples**

    >>> from sympy import compose
    >>> from sympy.abc import x

    >>> compose(x**2 + x, x - 1)
    -x + x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('compose', 2, exc)

    result = F.compose(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def decompose(f, *gens, **args):
    """
    Compute functional decomposition of ``f``.

    **Examples**

    >>> from sympy import decompose
    >>> from sympy.abc import x

    >>> decompose(x**4 + 2*x**3 - x - 1)
    [-1 - x + x**2, x + x**2]

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('decompose', 1, exc)

    result = F.decompose()

    if not opt.polys:
        return [ r.as_expr() for r in result ]
    else:
        return result

def sturm(f, *gens, **args):
    """
    Compute Sturm sequence of ``f``.

    **Examples**

    >>> from sympy import sturm
    >>> from sympy.abc import x

    >>> sturm(x**3 - 2*x**2 + x - 3)
    [-3 + x - 2*x**2 + x**3, 1 - 4*x + 3*x**2, 25/9 + 2*x/9, -2079/4]

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('sturm', 1, exc)

    result = F.sturm(auto=opt.auto)

    if not opt.polys:
        return [ r.as_expr() for r in result ]
    else:
        return result

def gff_list(f, *gens, **args):
    """
    Compute a list of greatest factorial factors of ``f``.

    **Examples**

    >>> from sympy import gff_list, ff
    >>> from sympy.abc import x

    >>> f = x**5 + 2*x**4 - x**3 - 2*x**2

    >>> gff_list(f)
    [(x, 1), (2 + x, 4)]

    >>> (ff(x, 1)*ff(x + 2, 4)).expand() == f
    True

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('gff_list', 1, exc)

    factors = F.gff_list()

    if not opt.polys:
        return [ (g.as_expr(), k) for g, k in factors ]
    else:
        return factors

def gff(f, *gens, **args):
    """Compute greatest factorial factorization of ``f``. """
    raise NotImplementedError('symbolic falling factorial')

def sqf_norm(f, *gens, **args):
    """
    Compute square-free norm of ``f``.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and
    ``r(x) = Norm(g(x))`` is a square-free polynomial over ``K``,
    where ``a`` is the algebraic extension of the ground domain.

    **Examples**

    >>> from sympy import sqf_norm, sqrt
    >>> from sympy.abc import x

    >>> sqf_norm(x**2 + 1, extension=[sqrt(3)])
    (1, 4 - 2*x*3**(1/2) + x**2, 16 - 4*x**2 + x**4)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('sqf_norm', 1, exc)

    s, g, r = F.sqf_norm()

    if not opt.polys:
        return Integer(s), g.as_expr(), r.as_expr()
    else:
        return Integer(s), g, r

def sqf_part(f, *gens, **args):
    """
    Compute square-free part of ``f``.

    **Examples**

    >>> from sympy import sqf_part
    >>> from sympy.abc import x

    >>> sqf_part(x**3 - 3*x - 2)
    -2 - x + x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('sqf_part', 1, exc)

    result = F.sqf_part()

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def _sorted_factors(factors, method):
    """Sort a list of ``(expr, exp)`` pairs. """
    if method == 'sqf':
        def key(obj):
            poly, exp = obj
            rep = poly.rep.rep
            return (exp, len(rep), rep)
    else:
        def key(obj):
            poly, exp = obj
            rep = poly.rep.rep
            return (len(rep), exp, rep)

    return sorted(factors, key=key)

def _factors_product(factors):
    """Multiply a list of ``(expr, exp)`` pairs. """
    return Mul(*[ f.as_expr()**k for f, k in factors ])

def _symbolic_factor_list(expr, opt, method):
    """Helper function for :func:`_symbolic_factor`. """
    coeff, factors = S.One, []

    for arg in Mul.make_args(expr):
        if arg.is_Pow:
            base, exp = arg.args
        else:
            base, exp = arg, S.One

        if base.is_Number:
            coeff *= arg
        else:
            try:
                poly, _ = _poly_from_expr(base, opt)
            except PolificationFailed, exc:
                coeff *= exc.expr**exp
            else:
                func = getattr(poly, method + '_list')

                _coeff, _factors = func()
                coeff *= _coeff**exp

                if exp is S.One:
                    factors.extend(_factors)
                else:
                    for factor, k in _factors:
                        factors.append((factor, k*exp))

    return coeff, factors

def _symbolic_factor(expr, opt, method):
    """Helper function for :func:`_factor`. """
    if isinstance(expr, Expr) and not expr.is_Relational:
        numer, denom = together(expr).as_numer_denom()

        cp, fp = _symbolic_factor_list(numer, opt, method)
        cq, fq = _symbolic_factor_list(denom, opt, method)

        fp = _factors_product(fp)
        fq = _factors_product(fq)

        coeff, expr = cp/cq, fp/fq

        return _keep_coeff(coeff, expr)
    elif hasattr(expr, 'args'):
        return expr.new(*[ _symbolic_factor(arg, opt, method) for arg in expr.args ])
    elif hasattr(expr, '__iter__'):
        return expr.__class__([ _symbolic_factor(arg, opt, method) for arg in expr ])
    else:
        return expr

def _generic_factor_list(expr, gens, args, method):
    """Helper function for :func:`sqf_list` and :func:`factor_list`. """
    options.allowed_flags(args, ['frac', 'polys'])
    opt = options.build_options(gens, args)

    expr = sympify(expr)

    if isinstance(expr, Expr) and not expr.is_Relational:
        numer, denom = together(expr).as_numer_denom()

        cp, fp = _symbolic_factor_list(numer, opt, method)
        cq, fq = _symbolic_factor_list(denom, opt, method)

        if fq and not opt.frac:
            raise PolynomialError("a polynomial expected, got %s" % expr)

        fp = _sorted_factors(fp, method)
        fq = _sorted_factors(fq, method)

        if not opt.polys:
            fp = [ (f.as_expr(), k) for f, k in fp ]
            fq = [ (f.as_expr(), k) for f, k in fq ]

        coeff = cp/cq

        if not opt.frac:
            return coeff, fp
        else:
            return coeff, fp, fq
    else:
        raise PolynomialError("a polynomial expected, got %s" % expr)

def _generic_factor(expr, gens, args, method):
    """Helper function for :func:`sqf` and :func:`factor`. """
    options.allowed_flags(args, [])
    opt = options.build_options(gens, args)
    return _symbolic_factor(sympify(expr), opt, method)

def sqf_list(f, *gens, **args):
    """
    Compute a list of square-free factors of ``f``.

    **Examples**

    >>> from sympy import sqf_list
    >>> from sympy.abc import x

    >>> sqf_list(2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16)
    (2, [(1 + x, 2), (2 + x, 3)])

    """
    return _generic_factor_list(f, gens, args, method='sqf')

def sqf(f, *gens, **args):
    """
    Compute square-free factorization of ``f``.

    **Examples**

    >>> from sympy import sqf
    >>> from sympy.abc import x

    >>> sqf(2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16)
    2*(1 + x)**2*(2 + x)**3

    """
    return _generic_factor(f, gens, args, method='sqf')

def factor_list(f, *gens, **args):
    """
    Compute a list of irreducible factors of ``f``.

    **Examples**

    >>> from sympy import factor_list
    >>> from sympy.abc import x, y

    >>> factor_list(2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y)
    (2, [(x + y, 1), (1 + x**2, 2)])

    """
    return _generic_factor_list(f, gens, args, method='factor')

def factor(f, *gens, **args):
    """
    Compute the factorization of ``f`` into irreducibles.

    There two modes implemented: symbolic and formal. If ``f`` is not an
    instance of :class:`Poly` and generators are not specified, then the
    former mode is used. Otherwise, the formal mode is used.

    In symbolic mode, :func:`factor` will traverse the expression tree and
    factor its components without any prior expansion, unless an instance
    of :class:`Add` is encountered (in this case formal factorization is
    used). This way :func:`factor` can handle large or symbolic exponents.

    By default, the factorization is computed over the rationals. To factor
    over other domain, e.g. an algebraic or finite field, use appropriate
    options: ``extension``, ``modulus`` or ``domain``.

    **Examples**

    >>> from sympy import factor, sqrt
    >>> from sympy.abc import x, y

    >>> factor(2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y)
    2*(1 + x**2)**2*(x + y)

    >>> factor(x**2 + 1)
    1 + x**2
    >>> factor(x**2 + 1, modulus=2)
    (1 + x)**2
    >>> factor(x**2 + 1, gaussian=True)
    (x + I)*(x - I)

    >>> factor(x**2 - 2, extension=sqrt(2))
    (x + 2**(1/2))*(x - 2**(1/2))

    >>> factor((x**2 - 1)/(x**2 + 4*x + 4))
    -(1 + x)*(1 - x)/(2 + x)**2
    >>> factor((x**2 + 4*x + 4)**10000000*(x**2 + 1))
    (2 + x)**20000000*(1 + x**2)

    """
    return _generic_factor(f, gens, args, method='factor')

def intervals(F, all=False, eps=None, inf=None, sup=None, strict=False, fast=False, sqf=False):
    """
    Compute isolating intervals for roots of ``f``.

    **Examples**

    >>> from sympy import intervals
    >>> from sympy.abc import x

    >>> intervals(x**2 - 3)
    [((-2, -1), 1), ((1, 2), 1)]
    >>> intervals(x**2 - 3, eps=1e-2)
    [((-26/15, -19/11), 1), ((19/11, 26/15), 1)]

    """
    if not hasattr(F, '__iter__'):
        try:
            F = Poly(F)
        except GeneratorsNeeded:
            return []

        return F.intervals(all=all, eps=eps, inf=inf, sup=sup, fast=fast, sqf=sqf)
    else:
        polys, opt = parallel_poly_from_expr(F, domain='QQ')

        if len(opt.gens) > 1:
            raise MultivariatePolynomialError

        for i, poly in enumerate(polys):
            polys[i] = poly.rep.rep

        if eps is not None:
            eps = opt.domain.convert(eps)
        if inf is not None:
            inf = opt.domain.convert(inf)
        if sup is not None:
            sup = opt.domain.convert(sup)

        intervals = dup_isolate_real_roots_list(polys, opt.domain,
            eps=eps, inf=inf, sup=sup, strict=strict, fast=fast)

        result = []

        for (s, t), indices in intervals:
            s, t = opt.domain.to_sympy(s), opt.domain.to_sympy(t)
            result.append(((s, t), indices))

        return result

def refine_root(f, s, t, eps=None, steps=None, fast=False, check_sqf=False):
    """
    Refine an isolating interval of a root to the given precision.

    **Examples**

    >>> from sympy import refine_root
    >>> from sympy.abc import x

    >>> refine_root(x**2 - 3, 1, 2, eps=1e-2)
    (19/11, 26/15)

    """
    try:
        F = Poly(f)
    except GeneratorsNeeded:
        raise PolynomialError("can't refine a root of %s, not a polynomial" % f)

    return F.refine_root(s, t, eps=eps, steps=steps, fast=fast, check_sqf=check_sqf)

def count_roots(f, inf=None, sup=None):
    """
    Return the number of roots of ``f`` in ``[inf, sup]`` interval.

    If one of ``inf`` or ``sup`` is complex, it will return the number of roots
    in the complex rectangle with corners at ``inf`` and ``sup``.

    **Examples**

    >>> from sympy import count_roots, I
    >>> from sympy.abc import x

    >>> count_roots(x**4 - 4, -3, 3)
    2
    >>> count_roots(x**4 - 4, 0, 1 + 3*I)
    1

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError("can't count roots of %s, not a polynomial" % f)

    return F.count_roots(inf=inf, sup=sup)

def real_roots(f, multiple=True):
    """
    Return a list of real roots with multiplicities of ``f``.

    **Examples**

    >>> from sympy import real_roots
    >>> from sympy.abc import x

    >>> real_roots(2*x**3 - 7*x**2 + 4*x + 4)
    [-1/2, 2, 2]

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError("can't compute real roots of %s, not a polynomial" % f)

    return F.real_roots(multiple=multiple)

def nroots(f, maxsteps=50, cleanup=True, error=False):
    """
    Compute numerical approximations of roots of ``f``.

    **Examples**

    >>> from sympy import nroots
    >>> from sympy.abc import x

    >>> nroots(x**2 - 3)
    [-1.73205080756888, 1.73205080756888]

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError("can't compute numerical roots of %s, not a polynomial" % f)

    return F.nroots(maxsteps=maxsteps, cleanup=cleanup, error=error)

def ground_roots(f, *gens, **args):
    """
    Compute roots of ``f`` by factorization in the ground domain.

    **Examples**

    >>> from sympy import ground_roots
    >>> from sympy.abc import x

    >>> ground_roots(x**6 - 4*x**4 + 4*x**3 - x**2)
    {1: 2, 0: 2}

    """
    options.allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('ground_roots', 1, exc)

    return F.ground_roots()

def cancel(f, *gens, **args):
    """
    Cancel common factors in a rational function ``f``.

    **Examples**

    >>> from sympy import cancel
    >>> from sympy.abc import x

    >>> cancel((2*x**2 - 2)/(x**2 - 2*x + 1))
    -(2 + 2*x)/(1 - x)

    """
    options.allowed_flags(args, ['polys'])

    f = sympify(f)

    if type(f) is not tuple:
        if f.is_Number:
            return f
        else:
            p, q = f.as_numer_denom()
    else:
        p, q = f

    try:
        (F, G), opt = parallel_poly_from_expr((p, q), *gens, **args)
    except PolificationFailed, exc:
        if type(f) is not tuple:
            return f
        else:
            return S.One, p, q

    c, P, Q = F.cancel(G)

    if type(f) is not tuple:
        return c*(P.as_expr()/Q.as_expr())
    else:
        if not opt.polys:
            return c, P.as_expr(), Q.as_expr()
        else:
            return c, P, Q

def reduced(f, G, *gens, **args):
    """
    Reduces a polynomial ``f`` modulo a set of polynomials ``G``.

    Given a polynomial ``f`` and a set of polynomials ``G = (g_1, ..., g_n)``,
    computes a set of quotients ``q = (q_1, ..., q_n)`` and the remainder ``r``
    such that ``f = q_1*f_1 + ... + q_n*f_n + r``, where ``r`` vanishes or ``r``
    is a completely reduced polynomial with respect to ``G``.

    **Examples**

    >>> from sympy import reduced
    >>> from sympy.abc import x, y

    >>> reduced(2*x**4 + y**2 - x**2 + y**3, [x**3 - x, y**3 - y])
    ([2*x, 1], y + x**2 + y**2)

    """
    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr([f] + list(G), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('reduced', 0, exc)

    for i, poly in enumerate(polys):
        polys[i] = sdp_from_dict(poly.rep.to_dict(), opt.order)

    level = len(opt.gens)-1

    Q, r = sdp_div(polys[0], polys[1:], level, opt.order, opt.domain)

    Q = [ Poly.new(DMP(dict(q), opt.domain, level), *opt.gens) for q in Q ]
    r =   Poly.new(DMP(dict(r), opt.domain, level), *opt.gens)

    if not opt.polys:
        return [ q.as_expr() for q in Q ], r.as_expr()
    else:
        return Q, r

def groebner(F, *gens, **args):
    """
    Computes the reduced Groebner basis for a set of polynomials.

    Use the ``order`` argument to set the monomial ordering that will be
    used to compute the basis. Allowed orders are ``lex``, ``grlex`` and
    ``grevlex``. If no order is specified, it defaults to ``lex``.

    **Examples**

    >>> from sympy import groebner
    >>> from sympy.abc import x, y

    >>> groebner([x*y - 2*y, 2*y**2 - x**2], order='lex')
    [x**2 - 2*y**2, -2*y + x*y, -2*y + y**3]
    >>> groebner([x*y - 2*y, 2*y**2 - x**2], order='grlex')
    [-2*y + y**3, x**2 - 2*y**2, -2*y + x*y]
    >>> groebner([x*y - 2*y, 2*y**2 - x**2], order='grevlex')
    [-2*x**2 + x**3, y**2 - x**2/2, -2*y + x*y]

    **References**

    1. [Buchberger01]_
    2. [Cox97]_

    """
    options.allowed_flags(args, ['monic', 'polys'])
    args = _update_args(args, 'field', True)

    try:
        polys, opt = parallel_poly_from_expr(F, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('groebner', len(F), exc)

    for i, poly in enumerate(polys):
        polys[i] = sdp_from_dict(poly.rep.to_dict(), opt.order)

    level = len(opt.gens)-1

    G = sdp_groebner(polys, level, opt.order, opt.domain, monic=opt.monic)
    G = [ Poly.new(DMP(dict(g), opt.domain, level), *opt.gens) for g in G ]

    if not opt.polys:
        return [ g.as_expr() for g in G ]
    else:
        return G

def poly(expr, **args):
    """
    Efficiently transform an expression into a polynomial.

    **Examples**

    >>> from sympy import poly
    >>> from sympy.abc import x

    >>> poly(x*(x**2 + x - 1)**2)
    Poly(x**5 + 2*x**4 - x**3 - 2*x**2 + x, x, domain='ZZ')

    """
    expr = sympify(expr)

    if expr.is_Poly:
        return expr.reorder(**args)

    terms, poly_terms = [], []

    for term in Add.make_args(expr):
        factors, poly_factors = [], []

        for factor in Mul.make_args(term):
            if factor.is_Add:
                poly_factors.append(poly(factor))
            elif factor.is_Pow and factor.base.is_Add and factor.exp.is_Integer:
                poly_factors.append(poly(factor.base).pow(factor.exp))
            else:
                factors.append(factor)

        if not poly_factors:
            terms.append(term)
        else:
            product = poly_factors[0]

            for factor in poly_factors[1:]:
                product = product.mul(factor)

            if factors:
                factor = Mul(*factors)

                if factor.is_Number:
                    product = product.mul(factor)
                else:
                    product = product.mul(Poly(factor, expand=False))

            poly_terms.append(product)

    if not poly_terms:
        result = Poly(expr, expand=False)
    else:
        result = poly_terms[0]

        for term in poly_terms[1:]:
            result = result.add(term)

        if terms:
            term = Add(*terms)

            if term.is_Number:
                result = result.add(term)
            else:
                result = result.add(Poly(term, expand=False))

    return result.reorder(**args)

