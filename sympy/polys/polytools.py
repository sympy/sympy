"""User-friendly public interface to polynomial functions. """

from sympy.core import (
    S, Basic, Expr, I, Integer, Add, Mul, Dummy,
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

        if hasattr(rep, '__iter__'):
            if isinstance(rep, dict):
                return cls._from_dict(rep, opt)
            else:
                return cls._from_list(list(rep), opt)
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
        [x**2 + 1]

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
    def unit(self):
        """Return unit polynomial with ``self``'s properties. """
        return self.new(self.rep.unit(self.rep.lev, self.rep.dom), *self.gens)

    def unify(self, other):
        """
        Make ``self`` and ``other`` belong to the same domain.

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
        _, per, F, G = self._unify(other)
        return per(F), per(G)

    def _unify(self, other):
        other = sympify(other)

        if not other.is_Poly:
            try:
                return self.rep.dom, self.per, self.rep, self.rep.per(self.rep.dom.from_sympy(other))
            except CoercionFailed:
                raise UnificationFailed("can't unify %s with %s" % (self, other))

        if isinstance(self.rep, DMP) and isinstance(other.rep, DMP):
            gens = _unify_gens(self.gens, other.gens)

            dom, lev = self.rep.dom.unify(other.rep.dom, gens), len(gens)-1

            if self.gens != gens:
                f_monoms, f_coeffs = _dict_reorder(self.rep.to_dict(), self.gens, gens)

                if self.rep.dom != dom:
                    f_coeffs = [ dom.convert(c, self.rep.dom) for c in f_coeffs ]

                F = DMP(dict(zip(f_monoms, f_coeffs)), dom, lev)
            else:
                F = self.rep.convert(dom)

            if other.gens != gens:
                g_monoms, g_coeffs = _dict_reorder(other.rep.to_dict(), other.gens, gens)

                if other.rep.dom != dom:
                    g_coeffs = [ dom.convert(c, other.rep.dom) for c in g_coeffs ]

                G = DMP(dict(zip(g_monoms, g_coeffs)), dom, lev)
            else:
                G = other.rep.convert(dom)
        else:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        def per(rep, dom=dom, gens=gens, remove=None):
            if remove is not None:
                gens = gens[:remove]+gens[remove+1:]

                if not gens:
                    return dom.to_sympy(rep)

            return Poly.new(rep, *gens)

        return dom, per, F, G

    def per(self, rep, gens=None, remove=None):
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
            gens = self.gens

        if remove is not None:
            gens = gens[:remove]+gens[remove+1:]

            if not gens:
                return self.rep.dom.to_sympy(rep)

        return Poly.new(rep, *gens)

    def set_domain(self, domain):
        """Set the ground domain of ``self``. """
        opt = options.build_options(self.gens, {'domain': domain})
        return self.per(self.rep.convert(opt.domain))

    def get_domain(self):
        """Get the ground domain of ``self``. """
        return self.rep.dom

    def set_modulus(self, modulus):
        """
        Set the modulus of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(5*x**2 + 2*x - 1, x).set_modulus(2)
        Poly(x**2 + 1, x, modulus=2)

        """
        modulus = options.Modulus.preprocess(modulus)
        return self.set_domain(FF(modulus))

    def get_modulus(self):
        """
        Get the modulus of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, modulus=2).get_modulus()
        2

        """
        domain = self.get_domain()

        if not domain.has_CharacteristicZero:
            return Integer(domain.characteristic())
        else:
            raise PolynomialError("not a polynomial over a Galois field")

    def _eval_subs(self, old, new):
        """Internal implementation of :func:`subs`. """
        if old in self.gens:
            if new.is_number:
                return self.eval(old, new)
            else:
                try:
                    return self.replace(old, new)
                except PolynomialError:
                    pass

        return self.as_expr().subs(old, new)

    def exclude(self):
        """
        Remove unnecessary generators from ``self``.

        **Example**

        >>> from sympy import Poly
        >>> from sympy.abc import a, b, c, d, x

        >>> Poly(a + x, a, b, c, d, x).exclude()
        Poly(a + x, a, x, domain='ZZ')

        """
        J, new = self.rep.exclude()
        gens = []

        for j in range(len(self.gens)):
            if j not in J:
                gens.append(self.gens[j])

        return self.per(new, gens=gens)

    def replace(self, x, y=None):
        """
        Replace ``x`` with ``y`` in generators list.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 1, x).replace(x, y)
        Poly(y**2 + 1, y, domain='ZZ')

        """
        if y is None:
            if self.is_univariate:
                x, y = self.gen, x
            else:
                raise PolynomialError("syntax supported only in univariate case")

        if x == y:
            return self

        if x in self.gens and y not in self.gens:
            dom = self.get_domain()

            if not dom.is_Composite or y not in dom.gens:
                gens = list(self.gens)
                gens[gens.index(x)] = y
                return self.per(self.rep, gens=gens)

        raise PolynomialError("can't replace %s with %s in %s" % (x, y, self))

    def reorder(self, *gens, **args):
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
            gens = _sort_gens(self.gens, opt=opt)
        elif set(self.gens) != set(gens):
            raise PolynomialError("generators list can differ only up to order of elements")

        rep = dict(zip(*_dict_reorder(self.rep.to_dict(), self.gens, gens)))

        return self.per(DMP(rep, self.rep.dom, len(gens)-1), gens=gens)

    def ltrim(self, gen):
        """
        Remove dummy generators from the "left" of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y, z

        >>> Poly(y**2 + y*z**2, x, y, z).ltrim(y)
        Poly(y**2 + y*z**2, y, z, domain='ZZ')

        """
        rep = self.as_dict(native=True)
        j = self._gen_to_level(gen)
        terms = {}

        for monom, coeff in rep.iteritems():
            monom = monom[j:]

            if monom not in terms:
                terms[monom] = coeff
            else:
                raise PolynomialError("can't left trim %s" % self)

        gens = self.gens[j:]

        return self.new(DMP.from_dict(terms, len(gens)-1, self.rep.dom), *gens)

    def has_only_gens(self, *gens):
        """
        Return ``True`` if ``Poly(self, *gens)`` retains ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y, z

        >>> Poly(x*y + 1, x, y, z).has_only_gens(x, y)
        True
        >>> Poly(x*y + z, x, y, z).has_only_gens(x, y)
        False

        """
        f_gens = list(self.gens)
        indices = set([])

        for gen in gens:
            try:
                index = f_gens.index(gen)
            except ValueError:
                raise GeneratorsError("%s doesn't have %s as generator" % (self, gen))
            else:
                indices.add(index)

        for monom in self.monoms():
            for i, elt in enumerate(monom):
                if i not in indices and elt:
                    return False

        return True

    def to_ring(self):
        """
        Make the ground domain a ring.

        **Examples**

        >>> from sympy import Poly, QQ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, domain=QQ).to_ring()
        Poly(x**2 + 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'to_ring'):
            result = self.rep.to_ring()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'to_ring')

        return self.per(result)

    def to_field(self):
        """
        Make the ground domain a field.

        **Examples**

        >>> from sympy import Poly, ZZ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x, domain=ZZ).to_field()
        Poly(x**2 + 1, x, domain='QQ')

        """
        if hasattr(self.rep, 'to_field'):
            result = self.rep.to_field()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'to_field')

        return self.per(result)

    def to_exact(self):
        """
        Make the ground domain exact.

        **Examples**

        >>> from sympy import Poly, RR
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1.0, x, domain=RR).to_exact()
        Poly(x**2 + 1, x, domain='QQ')

        """
        if hasattr(self.rep, 'to_exact'):
            result = self.rep.to_exact()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'to_exact')

        return self.per(result)

    def retract(self, field=None):
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
        dom, rep = construct_domain(self.as_dict(zero=True), field=field)
        return self.from_dict(rep, self.gens, domain=dom)

    def slice(self, x, m, n=None):
        """Take a continuous subsequence of terms of ``self``. """
        if n is None:
            j, m, n = 0, x, m
        else:
            j = self._gen_to_level(x)

        m, n = int(m), int(n)

        if hasattr(self.rep, 'slice'):
            result = self.rep.slice(m, n, j)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'slice')

        return self.per(result)

    def coeffs(self, order=None):
        """
        Returns all non-zero coefficients from ``self`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x + 3, x).coeffs()
        [1, 2, 3]

        """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.coeffs(order=order) ]

    def monoms(self, order=None):
        """
        Returns all non-zero monomials from ``self`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).monoms()
        [(2, 0), (1, 2), (1, 1), (0, 1)]

        """
        return self.rep.monoms(order=order)

    def terms(self, order=None):
        """
        Returns all non-zero terms from ``self`` in lex order.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).terms()
        [((2, 0), 1), ((1, 2), 2), ((1, 1), 1), ((0, 1), 3)]

        """
        return [ (m, self.rep.dom.to_sympy(c)) for m, c in self.rep.terms(order=order) ]

    def all_coeffs(self):
        """
        Returns all coefficients from a univariate polynomial ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_coeffs()
        [1, 0, 2, -1]

        """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.all_coeffs() ]

    def all_monoms(self):
        """
        Returns all monomials from a univariate polynomial ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_monoms()
        [(3,), (2,), (1,), (0,)]

        """
        return self.rep.all_monoms()

    def all_terms(self):
        """
        Returns all terms from a univariate polynomial ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x - 1, x).all_terms()
        [((3,), 1), ((2,), 0), ((1,), 2), ((0,), -1)]

        """
        return [ (m, self.rep.dom.to_sympy(c)) for m, c in self.rep.all_terms() ]

    def termwise(self, func, *gens, **args):
        """
        Apply a function to all terms of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> def func((k,), coeff):
        ...     return coeff//10**(2-k)

        >>> Poly(x**2 + 20*x + 400).termwise(func)
        Poly(x**2 + 2*x + 4, x, domain='ZZ')

        """
        terms = {}

        for monom, coeff in self.terms():
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

        return self.from_dict(terms, *(gens or self.gens), **args)

    def length(self):
        """
        Returns the number of non-zero terms in ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 2*x - 1).length()
        3

        """
        return len(self.as_dict())

    def as_dict(self, native=False, zero=False):
        """
        Switch to a ``dict`` representation.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x*y**2 - y, x, y).as_dict()
        {(0, 1): -1, (1, 2): 2, (2, 0): 1}

        """
        if native:
            return self.rep.to_dict(zero=zero)
        else:
            return self.rep.to_sympy_dict(zero=zero)

    def as_list(self, native=False):
        """Switch to a ``list`` representation. """
        if native:
            return self.rep.to_list()
        else:
            return self.rep.to_sympy_list()

    def as_expr(self, *gens):
        """
        Convert a polynomial an expression.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2 + 2*x*y**2 - y, x, y)

        >>> f.as_expr()
        x**2 + 2*x*y**2 - y
        >>> f.as_expr({x: 5})
        10*y**2 - y + 25
        >>> f.as_expr(5, 6)
        379

        """
        if not gens:
            gens = self.gens
        elif len(gens) == 1 and isinstance(gens[0], dict):
            mapping = gens[0]
            gens = list(self.gens)

            for gen, value in mapping.iteritems():
                try:
                    index = gens.index(gen)
                except ValueError:
                    raise GeneratorsError("%s doesn't have %s as generator" % (self, gen))
                else:
                    gens[index] = value

        return basic_from_dict(self.rep.to_sympy_dict(), *gens)

    def lift(self):
        """
        Convert algebraic coefficients to rationals.

        **Examples**

        >>> from sympy import Poly, I
        >>> from sympy.abc import x

        >>> Poly(x**2 + I*x + 1, x, extension=I).lift()
        Poly(x**4 + 3*x**2 + 1, x, domain='QQ')

        """
        if hasattr(self.rep, 'lift'):
            result = self.rep.lift()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'lift')

        return self.per(result)

    def deflate(self):
        """
        Reduce degree of ``self`` by mapping ``x_i**m`` to ``y_i``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**6*y**2 + x**3 + 1, x, y).deflate()
        ((3, 2), Poly(x**2*y + x + 1, x, y, domain='ZZ'))

        """
        if hasattr(self.rep, 'deflate'):
            J, result = self.rep.deflate()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'deflate')

        return J, self.per(result)

    def inject(self, front=False):
        """
        Inject ground domain generators into ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x)

        >>> f.inject()
        Poly(x**2*y + x*y**3 + x*y + 1, x, y, domain='ZZ')
        >>> f.inject(front=True)
        Poly(y**3*x + y*x**2 + y*x + 1, y, x, domain='ZZ')

        """
        dom = self.rep.dom

        if dom.is_Numerical:
            return self
        elif not dom.is_Poly:
            raise DomainError("can't inject generators over %s" % dom)

        if hasattr(self.rep, 'inject'):
            result = self.rep.inject(front=front)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'inject')

        if front:
            gens = dom.gens + self.gens
        else:
            gens = self.gens + dom.gens

        return self.new(result, *gens)

    def eject(self, *gens):
        """
        Eject selected generators into the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x, y)

        >>> f.eject(x)
        Poly(x*y**3 + (x**2 + x)*y + 1, y, domain='ZZ[x]')
        >>> f.eject(y)
        Poly(y*x**2 + (y**3 + y)*x + 1, x, domain='ZZ[y]')

        """
        dom = self.rep.dom

        if not dom.is_Numerical:
            raise DomainError("can't eject generators over %s" % dom)

        n, k = len(self.gens), len(gens)

        if self.gens[:k] == gens:
            _gens, front = self.gens[n-k:], True
        elif self.gens[-k:] == gens:
            _gens, front = self.gens[:n-k], False
        else:
            raise NotImplementedError("can only eject front or back generators")

        dom = dom.inject(*gens)

        if hasattr(self.rep, 'eject'):
            result = self.rep.eject(dom, front=front)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'eject')

        return self.new(result, *_gens)

    def terms_gcd(self):
        """
        Remove GCD of terms from the polynomial ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**6*y**2 + x**3*y, x, y).terms_gcd()
        ((3, 1), Poly(x**3*y + 1, x, y, domain='ZZ'))

        """
        if hasattr(self.rep, 'terms_gcd'):
            J, result = self.rep.terms_gcd()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'terms_gcd')

        return J, self.per(result)

    def add_ground(self, coeff):
        """
        Add an element of the ground domain to ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).add_ground(2)
        Poly(x + 3, x, domain='ZZ')

        """
        if hasattr(self.rep, 'add_ground'):
            result = self.rep.add_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'add_ground')

        return self.per(result)

    def sub_ground(self, coeff):
        """
        Subtract an element of the ground domain from ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).sub_ground(2)
        Poly(x - 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'sub_ground'):
            result = self.rep.sub_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sub_ground')

        return self.per(result)

    def mul_ground(self, coeff):
        """
        Multiply ``self`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 1).mul_ground(2)
        Poly(2*x + 2, x, domain='ZZ')

        """
        if hasattr(self.rep, 'mul_ground'):
            result = self.rep.mul_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'mul_ground')

        return self.per(result)

    def quo_ground(self, coeff):
        """
        Quotient of ``self`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x + 4).quo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).quo_ground(2)
        Poly(x + 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'quo_ground'):
            result = self.rep.quo_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'quo_ground')

        return self.per(result)

    def exquo_ground(self, coeff):
        """
        Exact quotient of ``self`` by a an element of the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x + 4).exquo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).exquo_ground(2)
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2 does not divide 3 in ZZ

        """
        if hasattr(self.rep, 'exquo_ground'):
            result = self.rep.exquo_ground(coeff)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'exquo_ground')

        return self.per(result)

    def abs(self):
        """
        Make all coefficients in ``self`` positive.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).abs()
        Poly(x**2 + 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'abs'):
            result = self.rep.abs()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'abs')

        return self.per(result)

    def neg(self):
        """
        Negate all coefficients in ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).neg()
        Poly(-x**2 + 1, x, domain='ZZ')

        >>> -Poly(x**2 - 1, x)
        Poly(-x**2 + 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'neg'):
            result = self.rep.neg()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'neg')

        return self.per(result)

    def add(self, other):
        """
        Add two polynomials ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).add(Poly(x - 2, x))
        Poly(x**2 + x - 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x) + Poly(x - 2, x)
        Poly(x**2 + x - 1, x, domain='ZZ')

        """
        other = sympify(other)

        if not other.is_Poly:
            return self.add_ground(other)

        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'add'):
            result = F.add(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'add')

        return per(result)

    def sub(self, other):
        """
        Subtract two polynomials ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).sub(Poly(x - 2, x))
        Poly(x**2 - x + 3, x, domain='ZZ')

        >>> Poly(x**2 + 1, x) - Poly(x - 2, x)
        Poly(x**2 - x + 3, x, domain='ZZ')

        """
        other = sympify(other)

        if not other.is_Poly:
            return self.sub_ground(other)

        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'sub'):
            result = F.sub(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sub')

        return per(result)

    def mul(self, other):
        """
        Multiply two polynomials ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).mul(Poly(x - 2, x))
        Poly(x**3 - 2*x**2 + x - 2, x, domain='ZZ')

        >>> Poly(x**2 + 1, x)*Poly(x - 2, x)
        Poly(x**3 - 2*x**2 + x - 2, x, domain='ZZ')

        """
        other = sympify(other)

        if not other.is_Poly:
            return self.mul_ground(other)

        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'mul'):
            result = F.mul(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'mul')

        return per(result)

    def sqr(self):
        """
        Square a polynomial ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x - 2, x).sqr()
        Poly(x**2 - 4*x + 4, x, domain='ZZ')

        >>> Poly(x - 2, x)**2
        Poly(x**2 - 4*x + 4, x, domain='ZZ')

        """
        if hasattr(self.rep, 'sqr'):
            result = self.rep.sqr()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sqr')

        return self.per(result)

    def pow(self, n):
        """
        Raise ``self`` to a non-negative power ``n``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x - 2, x).pow(3)
        Poly(x**3 - 6*x**2 + 12*x - 8, x, domain='ZZ')

        >>> Poly(x - 2, x)**3
        Poly(x**3 - 6*x**2 + 12*x - 8, x, domain='ZZ')

        """
        n = int(n)

        if hasattr(self.rep, 'pow'):
            result = self.rep.pow(n)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'pow')

        return self.per(result)

    def pdiv(self, other):
        """
        Polynomial pseudo-division of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).pdiv(Poly(2*x - 4, x))
        (Poly(2*x + 4, x, domain='ZZ'), Poly(20, x, domain='ZZ'))

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'pdiv'):
            q, r = F.pdiv(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'pdiv')

        return per(q), per(r)

    def prem(self, other):
        """
        Polynomial pseudo-remainder of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).prem(Poly(2*x - 4, x))
        Poly(20, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'prem'):
            result = F.prem(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'prem')

        return per(result)

    def pquo(self, other):
        """
        Polynomial pseudo-quotient of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).pquo(Poly(2*x - 4, x))
        Poly(2*x + 4, x, domain='ZZ')

        >>> Poly(x**2 - 1, x).pquo(Poly(2*x - 2, x))
        Poly(2*x + 2, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'pquo'):
            result = F.pquo(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'pquo')

        return per(result)

    def pexquo(self, other):
        """
        Polynomial exact pseudo-quotient of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).pexquo(Poly(2*x - 2, x))
        Poly(2*x + 2, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).pexquo(Poly(2*x - 4, x))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'pexquo'):
            try:
                result = F.pexquo(G)
            except ExactQuotientFailed, exc:
                raise exc.new(self.as_expr(), other.as_expr())
        else: # pragma: no cover
            raise OperationNotSupported(self, 'pexquo')

        return per(result)

    def div(self, other, auto=True):
        """
        Polynomial division with remainder of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).div(Poly(2*x - 4, x))
        (Poly(1/2*x + 1, x, domain='QQ'), Poly(5, x, domain='QQ'))

        >>> Poly(x**2 + 1, x).div(Poly(2*x - 4, x), auto=False)
        (Poly(0, x, domain='ZZ'), Poly(x**2 + 1, x, domain='ZZ'))

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.has_Ring and not dom.has_Field:
            F, G = F.to_field(), G.to_field()
            retract = True

        if hasattr(self.rep, 'div'):
            q, r = F.div(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'div')

        if retract:
            try:
                Q, R = q.to_ring(), r.to_ring()
            except CoercionFailed:
                pass
            else:
                q, r = Q, R

        return per(q), per(r)

    def rem(self, other, auto=True):
        """
        Computes the polynomial remainder of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly, ZZ, QQ
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).rem(Poly(2*x - 4, x))
        Poly(5, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).rem(Poly(2*x - 4, x), auto=False)
        Poly(x**2 + 1, x, domain='ZZ')

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.has_Ring and not dom.has_Field:
            F, G = F.to_field(), G.to_field()
            retract = True

        if hasattr(self.rep, 'rem'):
            r = F.rem(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'rem')

        if retract:
            try:
                r = r.to_ring()
            except CoercionFailed:
                pass

        return per(r)

    def quo(self, other, auto=True):
        """
        Computes polynomial quotient of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).quo(Poly(2*x - 4, x))
        Poly(1/2*x + 1, x, domain='QQ')

        >>> Poly(x**2 - 1, x).quo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.has_Ring and not dom.has_Field:
            F, G = F.to_field(), G.to_field()
            retract = True

        if hasattr(self.rep, 'quo'):
            q = F.quo(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'quo')

        if retract:
            try:
                q = q.to_ring()
            except CoercionFailed:
                pass

        return per(q)

    def exquo(self, other, auto=True):
        """
        Computes polynomial exact quotient of ``self`` by ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).exquo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).exquo(Poly(2*x - 4, x))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.has_Ring and not dom.has_Field:
            F, G = F.to_field(), G.to_field()
            retract = True

        if hasattr(self.rep, 'exquo'):
            try:
                q = F.exquo(G)
            except ExactQuotientFailed, exc:
                raise exc.new(self.as_expr(), other.as_expr())
        else: # pragma: no cover
            raise OperationNotSupported(self, 'exquo')

        if retract:
            try:
                q = q.to_ring()
            except CoercionFailed:
                pass

        return per(q)

    def _gen_to_level(self, gen):
        """Returns level associated with the given generator. """
        if isinstance(gen, int):
            length = len(self.gens)

            if -length <= gen < length:
                if gen < 0:
                    return length + gen
                else:
                    return gen
            else:
                raise PolynomialError("-%s <= gen < %s expected, got %s" % (length, length, gen))
        else:
            try:
                return list(self.gens).index(sympify(gen))
            except ValueError:
                raise PolynomialError("a valid generator expected, got %s" % gen)

    def degree(self, gen=0):
        """
        Returns degree of ``self`` in ``x_j``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).degree()
        2
        >>> Poly(x**2 + y*x + y, x, y).degree(y)
        1

        """
        j = self._gen_to_level(gen)

        if hasattr(self.rep, 'degree'):
            return self.rep.degree(j)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'degree')

    def degree_list(self):
        """
        Returns a list of degrees of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).degree_list()
        (2, 1)

        """
        if hasattr(self.rep, 'degree_list'):
            return self.rep.degree_list()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'degree_list')

    def total_degree(self):
        """
        Returns the total degree of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + y*x + 1, x, y).total_degree()
        3

        """
        if hasattr(self.rep, 'total_degree'):
            return self.rep.total_degree()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'total_degree')

    def LC(self, order=None):
        """
        Returns the leading coefficient of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(4*x**3 + 2*x**2 + 3*x, x).LC()
        4

        """
        if order is not None:
            return self.coeffs(order)[0]

        if hasattr(self.rep, 'LC'):
            result = self.rep.LC()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'LC')

        return self.rep.dom.to_sympy(result)

    def TC(self):
        """
        Returns the trailing coefficent of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x**2 + 3*x, x).TC()
        0

        """
        if hasattr(self.rep, 'TC'):
            result = self.rep.TC()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'TC')

        return self.rep.dom.to_sympy(result)

    def EC(self, order=None):
        """
        Returns the last non-zero coefficient of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 + 2*x**2 + 3*x, x).EC()
        3

        """
        if hasattr(self.rep, 'coeffs'):
            return self.coeffs(order)[-1]
        else: # pragma: no cover
            raise OperationNotSupported(self, 'EC')

    def nth(self, *N):
        """
        Returns the ``n``-th coefficient of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**3 + 2*x**2 + 3*x, x).nth(2)
        2
        >>> Poly(x**3 + 2*x*y**2 + y**2, x, y).nth(1, 2)
        2

        """
        if hasattr(self.rep, 'nth'):
            result = self.rep.nth(*map(int, N))
        else: # pragma: no cover
            raise OperationNotSupported(self, 'nth')

        return self.rep.dom.to_sympy(result)

    def LM(self, order=None):
        """
        Returns the leading monomial of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LM()
        (2, 0)

        """
        return self.monoms(order)[0]

    def EM(self, order=None):
        """
        Returns the last non-zero monomial of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).EM()
        (0, 1)

        """
        return self.monoms(order)[-1]

    def LT(self, order=None):
        """
        Returns the leading term of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LT()
        ((2, 0), 4)

        """
        return self.terms(order)[0]

    def ET(self, order=None):
        """
        Returns the last non-zero term of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).ET()
        ((0, 1), 3)

        """
        return self.terms(order)[-1]

    def max_norm(self):
        """
        Returns maximum norm of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(-x**2 + 2*x - 3, x).max_norm()
        3

        """
        if hasattr(self.rep, 'max_norm'):
            result = self.rep.max_norm()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'max_norm')

        return self.rep.dom.to_sympy(result)

    def l1_norm(self):
        """
        Returns l1 norm of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(-x**2 + 2*x - 3, x).l1_norm()
        6

        """
        if hasattr(self.rep, 'l1_norm'):
            result = self.rep.l1_norm()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'l1_norm')

        return self.rep.dom.to_sympy(result)

    def clear_denoms(self, convert=False):
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
        if not self.rep.dom.has_Field:
            return S.One, self

        dom = self.rep.dom.get_ring()

        if hasattr(self.rep, 'clear_denoms'):
            coeff, result = self.rep.clear_denoms()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'clear_denoms')

        coeff, self = dom.to_sympy(coeff), self.per(result)

        if not convert:
            return coeff, self
        else:
            return coeff, self.to_ring()

    def rat_clear_denoms(self, other):
        """
        Clear denominators in a rational function ``self/other``.

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
        dom, per, self, other = self._unify(other)

        f = per(self)
        g = per(other)

        if not (dom.has_Field and dom.has_assoc_Ring):
            return f, g

        a, f = f.clear_denoms(convert=True)
        b, g = g.clear_denoms(convert=True)

        f = f.mul_ground(b)
        g = g.mul_ground(a)

        return f, g

    def integrate(self, *specs, **args):
        """
        Computes indefinite integral of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x + 1, x).integrate()
        Poly(1/3*x**3 + x**2 + x, x, domain='QQ')

        >>> Poly(x*y**2 + x, x, y).integrate((0, 1), (1, 0))
        Poly(1/2*x**2*y**2 + 1/2*x**2, x, y, domain='QQ')

        """
        if args.get('auto', True) and self.rep.dom.has_Ring:
            self = self.to_field()

        if hasattr(self.rep, 'integrate'):
            if not specs:
                return self.per(self.rep.integrate(m=1))

            rep = self.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.integrate(int(m), self._gen_to_level(gen))

            return self.per(rep)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'integrate')

    def diff(self, *specs):
        """
        Computes partial derivative of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x**2 + 2*x + 1, x).diff()
        Poly(2*x + 2, x, domain='ZZ')

        >>> Poly(x*y**2 + x, x, y).diff((0, 0), (1, 1))
        Poly(2*x*y, x, y, domain='ZZ')

        """
        if hasattr(self.rep, 'diff'):
            if not specs:
                return self.per(self.rep.diff(m=1))

            rep = self.rep

            for spec in specs:
                if type(spec) is tuple:
                    gen, m = spec
                else:
                    gen, m = spec, 1

                rep = rep.diff(int(m), self._gen_to_level(gen))

            return self.per(rep)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'diff')

    def eval(self, x, a=None, auto=True):
        """
        Evaluate ``self`` at ``a`` in the given variable.

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
                    self = self.eval(gen, value)

                return self
            else:
                j, a = 0, x
        else:
            j = self._gen_to_level(x)

        if not hasattr(self.rep, 'eval'): # pragma: no cover
            raise OperationNotSupported(self, 'eval')

        try:
            result = self.rep.eval(a, j)
        except CoercionFailed:
            if not auto:
                raise DomainError("can't evaluate at %s in %s" % (a, self.rep.dom))
            else:
                domain, [a] = construct_domain([a])
                self = self.set_domain(domain)
                result = self.rep.eval(a, j)

        return self.per(result, remove=j)

    def half_gcdex(self, other, auto=True):
        """
        Half extended Euclidean algorithm of ``self`` and ``other``.

        Returns ``(s, h)`` such that ``h = gcd(self, other)`` and ``s*self = h (mod other)``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> Poly(f).half_gcdex(Poly(g))
        (Poly(-1/5*x + 3/5, x, domain='QQ'), Poly(x + 1, x, domain='QQ'))

        """
        dom, per, F, G = self._unify(other)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(self.rep, 'half_gcdex'):
            s, h = F.half_gcdex(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'half_gcdex')

        return per(s), per(h)

    def gcdex(self, other, auto=True):
        """
        Extended Euclidean algorithm of ``self`` and ``other``.

        Returns ``(s, t, h)`` such that ``h = gcd(self, other)`` and ``s*self + t*other = h``.

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
        dom, per, F, G = self._unify(other)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(self.rep, 'gcdex'):
            s, t, h = F.gcdex(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'gcdex')

        return per(s), per(t), per(h)

    def invert(self, other, auto=True):
        """
        Invert ``self`` modulo ``other`` when possible.

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
        dom, per, F, G = self._unify(other)

        if auto and dom.has_Ring:
            F, G = F.to_field(), G.to_field()

        if hasattr(self.rep, 'invert'):
            result = F.invert(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'invert')

        return per(result)

    def revert(self, n):
        """Compute ``self**(-1)`` mod ``x**n``. """
        if hasattr(self.rep, 'revert'):
            result = self.rep.revert(int(n))
        else: # pragma: no cover
            raise OperationNotSupported(self, 'revert')

        return self.per(result)

    def subresultants(self, other):
        """
        Computes the subresultant PRS sequence of ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).subresultants(Poly(x**2 - 1, x))
        [Poly(x**2 + 1, x, domain='ZZ'),
         Poly(x**2 - 1, x, domain='ZZ'),
         Poly(-2, x, domain='ZZ')]

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'subresultants'):
            result = F.subresultants(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'subresultants')

        return map(per, result)

    def resultant(self, other):
        """
        Computes the resultant of ``self`` and ``other`` via PRS.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 1, x).resultant(Poly(x**2 - 1, x))
        4

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'resultant'):
            result = F.resultant(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'resultant')

        return per(result, remove=0)

    def discriminant(self):
        """
        Computes the discriminant of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + 2*x + 3, x).discriminant()
        -8

        """
        if hasattr(self.rep, 'discriminant'):
            result = self.rep.discriminant()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'discriminant')

        return self.per(result, remove=0)

    def cofactors(self, other):
        """
        Returns the GCD of ``self`` and ``other`` and their cofactors.

        Returns polynomials ``(h, cff, cfg)`` such that ``h = gcd(self, other)``, and
        ``cff = quo(self, h)`` and ``cfg = quo(other, h)`` are, so called, cofactors
        of ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).cofactors(Poly(x**2 - 3*x + 2, x))
        (Poly(x - 1, x, domain='ZZ'),
         Poly(x + 1, x, domain='ZZ'),
         Poly(x - 2, x, domain='ZZ'))

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'cofactors'):
            h, cff, cfg = F.cofactors(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'cofactors')

        return per(h), per(cff), per(cfg)

    def gcd(self, other):
        """
        Returns the polynomial GCD of ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).gcd(Poly(x**2 - 3*x + 2, x))
        Poly(x - 1, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'gcd'):
            result = F.gcd(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'gcd')

        return per(result)

    def lcm(self, other):
        """
        Returns polynomial LCM of ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 1, x).lcm(Poly(x**2 - 3*x + 2, x))
        Poly(x**3 - 2*x**2 - x + 2, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'lcm'):
            result = F.lcm(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'lcm')

        return per(result)

    def trunc(self, p):
        """
        Reduce ``self`` modulo a constant ``p``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**3 + 3*x**2 + 5*x + 7, x).trunc(3)
        Poly(-x**3 - x + 1, x, domain='ZZ')

        """
        p = self.rep.dom.convert(p)

        if hasattr(self.rep, 'trunc'):
            result = self.rep.trunc(p)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'trunc')

        return self.per(result)

    def monic(self, auto=True):
        """
        Divides all coefficients by ``LC(self)``.

        **Examples**

        >>> from sympy import Poly, ZZ
        >>> from sympy.abc import x

        >>> Poly(3*x**2 + 6*x + 9, x, domain=ZZ).monic()
        Poly(x**2 + 2*x + 3, x, domain='QQ')

        >>> Poly(3*x**2 + 4*x + 2, x, domain=ZZ).monic()
        Poly(x**2 + 4/3*x + 2/3, x, domain='QQ')

        """
        if auto and self.rep.dom.has_Ring:
            self = self.to_field()

        if hasattr(self.rep, 'monic'):
            result = self.rep.monic()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'monic')

        return self.per(result)

    def content(self):
        """
        Returns the GCD of polynomial coefficients.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(6*x**2 + 8*x + 12, x).content()
        2

        """
        if hasattr(self.rep, 'content'):
            result = self.rep.content()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'content')

        return self.rep.dom.to_sympy(result)

    def primitive(self):
        """
        Returns the content and a primitive form of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 + 8*x + 12, x).primitive()
        (2, Poly(x**2 + 4*x + 6, x, domain='ZZ'))

        """
        if hasattr(self.rep, 'primitive'):
            cont, result = self.rep.primitive()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'primitive')

        return self.rep.dom.to_sympy(cont), self.per(result)

    def compose(self, other):
        """
        Computes the functional composition of ``self`` and ``other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 + x, x).compose(Poly(x - 1, x))
        Poly(x**2 - x, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        if hasattr(self.rep, 'compose'):
            result = F.compose(G)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'compose')

        return per(result)

    def decompose(self):
        """
        Computes a functional decomposition of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**4 + 2*x**3 - x - 1, x, domain='ZZ').decompose()
        [Poly(x**2 - x - 1, x, domain='ZZ'), Poly(x**2 + x, x, domain='ZZ')]

        """
        if hasattr(self.rep, 'decompose'):
            result = self.rep.decompose()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'decompose')

        return map(self.per, result)

    def shift(self, a):
        """
        Efficiently compute Taylor shift ``self(x + a)``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 2*x + 1, x).shift(2)
        Poly(x**2 + 2*x + 1, x, domain='ZZ')

        """
        if hasattr(self.rep, 'shift'):
            result = self.rep.shift(a)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'shift')

        return self.per(result)

    def sturm(self, auto=True):
        """
        Computes the Sturm sequence of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 - 2*x**2 + x - 3, x).sturm()
        [Poly(x**3 - 2*x**2 + x - 3, x, domain='QQ'),
         Poly(3*x**2 - 4*x + 1, x, domain='QQ'),
         Poly(2/9*x + 25/9, x, domain='QQ'),
         Poly(-2079/4, x, domain='QQ')]

        """
        if auto and self.rep.dom.has_Ring:
            self = self.to_field()

        if hasattr(self.rep, 'sturm'):
            result = self.rep.sturm()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sturm')

        return map(self.per, result)

    def gff_list(self):
        """
        Computes greatest factorial factorization of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**5 + 2*x**4 - x**3 - 2*x**2

        >>> Poly(f).gff_list()
        [(Poly(x, x, domain='ZZ'), 1), (Poly(x + 2, x, domain='ZZ'), 4)]

        """
        if hasattr(self.rep, 'gff_list'):
            result = self.rep.gff_list()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'gff_list')

        return [ (self.per(other), k) for other, k in result ]

    def sqf_norm(self):
        """
        Computes square-free norm of ``self``.

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
        if hasattr(self.rep, 'sqf_norm'):
            s, g, r = self.rep.sqf_norm()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sqf_norm')

        return s, self.per(g), self.per(r)

    def sqf_part(self):
        """
        Computes square-free part of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**3 - 3*x - 2, x).sqf_part()
        Poly(x**2 - x - 2, x, domain='ZZ')

        """
        if hasattr(self.rep, 'sqf_part'):
            result = self.rep.sqf_part()
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sqf_part')

        return self.per(result)

    def sqf_list(self, all=False):
        """
        Returns a list of square-free factors of ``self``.

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
        if hasattr(self.rep, 'sqf_list'):
            coeff, factors = self.rep.sqf_list(all)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sqf_list')

        return self.rep.dom.to_sympy(coeff), [ (self.per(g), k) for g, k in factors ]

    def sqf_list_include(self, all=False):
        """
        Returns a list of square-free factors of ``self``.

        **Examples**

        >>> from sympy import Poly, expand
        >>> from sympy.abc import x

        >>> f = expand(2*(x + 1)**3*x**4)
        >>> f
        2*x**7 + 6*x**6 + 6*x**5 + 2*x**4

        >>> Poly(f).sqf_list_include()
        [(Poly(2, x, domain='ZZ'), 1),
         (Poly(x + 1, x, domain='ZZ'), 3),
         (Poly(x, x, domain='ZZ'), 4)]

        >>> Poly(f).sqf_list_include(all=True)
        [(Poly(2, x, domain='ZZ'), 1),
         (Poly(1, x, domain='ZZ'), 2),
         (Poly(x + 1, x, domain='ZZ'), 3),
         (Poly(x, x, domain='ZZ'), 4)]

        """
        if hasattr(self.rep, 'sqf_list_include'):
            factors = self.rep.sqf_list_include(all)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'sqf_list_include')

        return [ (self.per(g), k) for g, k in factors ]

    def factor_list(self):
        """
        Returns a list of irreducible factors of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = 2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y

        >>> Poly(f).factor_list()
        (2, [(Poly(x + y, x, y, domain='ZZ'), 1),
             (Poly(x**2 + 1, x, y, domain='ZZ'), 2)])

        """
        if hasattr(self.rep, 'factor_list'):
            try:
                coeff, factors = self.rep.factor_list()
            except DomainError:
                return S.One, [(self, 1)]
        else: # pragma: no cover
            raise OperationNotSupported(self, 'factor_list')

        return self.rep.dom.to_sympy(coeff), [ (self.per(g), k) for g, k in factors ]

    def factor_list_include(self):
        """
        Returns a list of irreducible factors of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> f = 2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y

        >>> Poly(f).factor_list_include()
        [(Poly(2*x + 2*y, x, y, domain='ZZ'), 1),
         (Poly(x**2 + 1, x, y, domain='ZZ'), 2)]

        """
        if hasattr(self.rep, 'factor_list_include'):
            try:
                factors = self.rep.factor_list_include()
            except DomainError:
                return [(self, 1)]
        else: # pragma: no cover
            raise OperationNotSupported(self, 'factor_list_include')

        return [ (self.per(g), k) for g, k in factors ]

    def intervals(self, all=False, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """
        Compute isolating intervals for roots of ``self``.

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

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

        if inf is not None:
            inf = QQ.convert(inf)
        if sup is not None:
            sup = QQ.convert(sup)

        if hasattr(self.rep, 'intervals'):
            result = self.rep.intervals(all=all, eps=eps, inf=inf, sup=sup, fast=fast, sqf=sqf)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'intervals')

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

    def refine_root(self, s, t, eps=None, steps=None, fast=False, check_sqf=False):
        """
        Refine an isolating interval of a root to the given precision.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 3, x).refine_root(1, 2, eps=1e-2)
        (19/11, 26/15)

        """
        if check_sqf and not self.is_sqf:
            raise PolynomialError("only square-free polynomials supported")

        s, t = QQ.convert(s), QQ.convert(t)

        if eps is not None:
            eps = QQ.convert(eps)

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

        if steps is not None:
            steps = int(steps)
        elif eps is None:
            steps = 1

        if hasattr(self.rep, 'refine_root'):
            S, T = self.rep.refine_root(s, t, eps=eps, steps=steps, fast=fast)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'refine_root')

        return QQ.to_sympy(S), QQ.to_sympy(T)

    def count_roots(self, inf=None, sup=None):
        """
        Return the number of roots of ``self`` in ``[inf, sup]`` interval.

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
            if hasattr(self.rep, 'count_real_roots'):
                count = self.rep.count_real_roots(inf=inf, sup=sup)
            else: # pragma: no cover
                raise OperationNotSupported(self, 'count_real_roots')
        else:
            if inf_real and inf is not None:
                inf = (inf, QQ.zero)

            if sup_real and sup is not None:
                sup = (sup, QQ.zero)

            if hasattr(self.rep, 'count_complex_roots'):
                count = self.rep.count_complex_roots(inf=inf, sup=sup)
            else: # pragma: no cover
                raise OperationNotSupported(self, 'count_complex_roots')

        return Integer(count)

    def real_roots(self, multiple=True):
        """
        Return a list of real roots with multiplicities.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**3 - 7*x**2 + 4*x + 4, x).real_roots()
        [-1/2, 2, 2]

        """
        reals = sympy.polys.rootoftools.RootOf(self)

        if multiple:
            return reals
        else:
            return group(reals, multiple=False)

    def nroots(self, maxsteps=50, cleanup=True, error=False):
        """
        Compute numerical approximations of roots of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 3).nroots()
        [-1.73205080756888, 1.73205080756888]

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError("can't compute numerical roots of %s" % self)

        if self.degree() <= 0:
            return []

        try:
            coeffs = [ complex(c) for c in self.all_coeffs() ]
        except ValueError:
            raise DomainError("numerical domain expected, got %s" % self.rep.dom)

        return sympify(npolyroots(coeffs, maxsteps=maxsteps, cleanup=cleanup, error=error))

    def ground_roots(self):
        """
        Compute roots of ``self`` by factorization in the ground domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**6 - 4*x**4 + 4*x**3 - x**2).ground_roots()
        {0: 2, 1: 2}

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError("can't compute ground roots of %s" % self)

        roots = {}

        for factor, k in self.factor_list()[1]:
            if factor.is_linear:
                a, b = factor.all_coeffs()
                roots[-b/a] = k

        return roots

    def nth_power_roots_poly(self, n):
        """
        Construct a polynomial with n-th powers of roots of ``self``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = Poly(x**4 - x**2 + 1)

        >>> f.nth_power_roots_poly(2)
        Poly(x**4 - 2*x**3 + 3*x**2 - 2*x + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(3)
        Poly(x**4 + 2*x**2 + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(4)
        Poly(x**4 + 2*x**3 + 3*x**2 + 2*x + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(12)
        Poly(x**4 - 4*x**3 + 6*x**2 - 4*x + 1, x, domain='ZZ')

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError("must be a univariate polynomial")

        N = sympify(n)

        if N.is_Integer and N >= 1:
            n = int(N)
        else:
            raise ValueError("'n' must an integer and n >= 1, got %s" % n)

        x = self.gen
        t = Dummy('t')

        r = self.resultant(self.__class__.from_expr(x**n - t, x, t))

        return r.replace(t, x)

    def cancel(self, other, include=False):
        """
        Cancel common factors in a rational function ``self/other``.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 - 2, x).cancel(Poly(x**2 - 2*x + 1, x))
        (1, Poly(2*x + 2, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'))

        >>> Poly(2*x**2 - 2, x).cancel(Poly(x**2 - 2*x + 1, x), include=True)
        (Poly(2*x + 2, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'))

        """
        dom, per, F, G = self._unify(other)

        if hasattr(F, 'cancel'):
            result = F.cancel(G, include=include)
        else: # pragma: no cover
            raise OperationNotSupported(self, 'cancel')

        if not include:
            if dom.has_assoc_Ring:
                dom = dom.get_ring()

            cp, cq, p, q = result

            cp = dom.to_sympy(cp)
            cq = dom.to_sympy(cq)

            return cp/cq, per(p), per(q)
        else:
            return tuple(map(per, result))

    @property
    def is_zero(self):
        """
        Returns ``True`` if ``self`` is a zero polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(0, x).is_zero
        True
        >>> Poly(1, x).is_zero
        False

        """
        return self.rep.is_zero

    @property
    def is_one(self):
        """
        Returns ``True`` if ``self`` is a unit polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(0, x).is_one
        False
        >>> Poly(1, x).is_one
        True

        """
        return self.rep.is_one

    @property
    def is_sqf(self):
        """
        Returns ``True`` if ``self`` is a square-free polynomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x**2 - 2*x + 1, x).is_sqf
        False
        >>> Poly(x**2 - 1, x).is_sqf
        True

        """
        return self.rep.is_sqf

    @property
    def is_monic(self):
        """
        Returns ``True`` if the leading coefficient of ``self`` is one.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(x + 2, x).is_monic
        True
        >>> Poly(2*x + 2, x).is_monic
        False

        """
        return self.rep.is_monic

    @property
    def is_primitive(self):
        """
        Returns ``True`` if GCD of the coefficients of ``self`` is one.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(2*x**2 + 6*x + 12, x).is_primitive
        False
        >>> Poly(x**2 + 3*x + 6, x).is_primitive
        True

        """
        return self.rep.is_primitive

    @property
    def is_ground(self):
        """
        Returns ``True`` if ``self`` is an element of the ground domain.

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
        return self.rep.is_ground

    @property
    def is_linear(self):
        """
        Returns ``True`` if ``self`` is linear in all its variables.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x + y + 2, x, y).is_linear
        True
        >>> Poly(x*y + 2, x, y).is_linear
        False

        """
        return self.rep.is_linear

    @property
    def is_quadratic(self):
        """
        Returns ``True`` if ``self`` is quadratic in all its variables.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x*y + 2, x, y).is_quadratic
        True
        >>> Poly(x*y**2 + 2, x, y).is_quadratic
        False

        """
        return self.rep.is_quadratic

    @property
    def is_monomial(self):
        """
        Returns ``True`` if ``self`` is zero or has only one term.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> Poly(3*x**2, x).is_monomial
        True
        >>> Poly(3*x**2 + 1, x).is_monomial
        False

        """
        return self.rep.is_monomial

    @property
    def is_homogeneous(self):
        """
        Returns ``True`` if ``self`` has zero trailing coefficient.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x, y

        >>> Poly(x*y + x + y, x, y).is_homogeneous
        True
        >>> Poly(x*y + x + y + 1, x, y).is_homogeneous
        False

        """
        return self.rep.is_homogeneous

    @property
    def is_irreducible(self):
        """
        Returns ``True`` if ``self`` has no factors over its domain.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >> Poly(x**2 + x + 1, x, modulus=2).is_irreducible
        True
        >> Poly(x**2 + 1, x, modulus=2).is_irreducible
        False

        """
        return self.rep.is_irreducible

    @property
    def is_univariate(self):
        """
        Returns ``True`` if ``self`` is a univariate polynomial.

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
        return len(self.gens) == 1

    @property
    def is_multivariate(self):
        """
        Returns ``True`` if ``self`` is a multivariate polynomial.

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
        return len(self.gens) != 1

    @property
    def is_cyclotomic(self):
        """
        Returns ``True`` if ``self`` is a cyclotomic polnomial.

        **Examples**

        >>> from sympy import Poly
        >>> from sympy.abc import x

        >>> f = x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1

        >>> Poly(f).is_cyclotomic
        False

        >>> g = x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1

        >>> Poly(g).is_cyclotomic
        True

        """
        return self.rep.is_cyclotomic

    def __abs__(self):
        return self.abs()

    def __neg__(self):
        return self.neg()

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return self.as_expr() + other

        return self.add(other)

    @_sympifyit('other', NotImplemented)
    def __radd__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return other + self.as_expr()

        return other.add(self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return self.as_expr() - other

        return self.sub(other)

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return other - self.as_expr()

        return other.sub(self)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return self.as_expr()*other

        return self.mul(other)

    @_sympifyit('other', NotImplemented)
    def __rmul__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens)
            except PolynomialError:
                return other*self.as_expr()

        return other.mul(self)

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        if other.is_Integer and other >= 0:
            return self.pow(other)
        else:
            return self.as_expr() ** other

    @_sympifyit('other', NotImplemented)
    def __divmod__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return self.div(other)

    @_sympifyit('other', NotImplemented)
    def __rdivmod__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return other.div(self)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return self.rem(other)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return other.rem(self)

    @_sympifyit('other', NotImplemented)
    def __floordiv__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return self.quo(other)

    @_sympifyit('other', NotImplemented)
    def __rfloordiv__(self, other):
        if not other.is_Poly:
            other = Poly(other, *self.gens)

        return other.quo(self)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        return self.as_expr()/other.as_expr()

    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        return other.as_expr()/self.as_expr()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        if not other.is_Poly:
            try:
                other = Poly(other, *self.gens, **{'domain': self.get_domain()})
            except (PolynomialError, DomainError, CoercionFailed):
                return False

        if self.gens != other.gens:
            return False

        if self.rep.dom != other.rep.dom:
            try:
                dom = self.rep.dom.unify(other.rep.dom, self.gens)
            except UnificationFailed:
                return False

            self = self.set_domain(dom)
            other = other.set_domain(dom)

        return self.rep == other.rep

    @_sympifyit('other', NotImplemented)
    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return super(Poly, self).__hash__()

    def __nonzero__(self):
        return not self.is_zero

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
    (2*x + 4, 20)

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

    >>> pquo(x**2 + 1, 2*x - 4)
    2*x + 4
    >>> pquo(x**2 - 1, 2*x - 1)
    2*x + 1

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

    >>> pexquo(x**2 - 1, 2*x - 2)
    2*x + 2

    >>> pexquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

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
    (0, x**2 + 1)
    >>> div(x**2 + 1, 2*x - 4, domain=QQ)
    (x/2 + 1, 5)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('div', 2, exc)

    q, r = F.div(G, auto=opt.auto)

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
    x**2 + 1
    >>> rem(x**2 + 1, 2*x - 4, domain=QQ)
    5

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('rem', 2, exc)

    r = F.rem(G, auto=opt.auto)

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

    >>> quo(x**2 + 1, 2*x - 4)
    x/2 + 1
    >>> quo(x**2 - 1, x - 1)
    x + 1

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('quo', 2, exc)

    q = F.quo(G, auto=opt.auto)

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

    >>> exquo(x**2 - 1, x - 1)
    x + 1

    >>> exquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('exquo', 2, exc)

    q = F.exquo(G, auto=opt.auto)

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
    (-x/5 + 3/5, x + 1)

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
    (-x/5 + 3/5, x**2/5 - 6*x/5 + 2, x + 1)

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
    [x**2 + 1, x**2 - 1, -2]

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
    (x - 1, x + 1, x - 2)

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
    x - 1

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
    x - 1

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
    x**5 - x**4 - 2*x**3 - x**2 + x + 2

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
    x**3 - 2*x**2 - x + 2

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
    x**3*y*(x**3*y + 1)

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
    -x**3 - x + 1

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
    x**2 + 4*x/3 + 2/3

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
    (2, 3*x**2 + 4*x + 6)

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
    x**2 - x

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
    [x**2 - x - 1, x**2 + x]

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
    [x**3 - 2*x**2 + x - 3, 3*x**2 - 4*x + 1, 2*x/9 + 25/9, -2079/4]

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
    [(x, 1), (x + 2, 4)]

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
    (1, x**2 - 2*3**(1/2)*x + 4, x**4 - 4*x**2 + 16)

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
    x**2 - x - 2

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
        coeff, factors = _symbolic_factor_list(together(expr), opt, method)
        return _keep_coeff(coeff, _factors_product(factors))
    elif hasattr(expr, 'args'):
        return expr.func(*[ _symbolic_factor(arg, opt, method) for arg in expr.args ])
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
    (2, [(x + 1, 2), (x + 2, 3)])

    """
    return _generic_factor_list(f, gens, args, method='sqf')

def sqf(f, *gens, **args):
    """
    Compute square-free factorization of ``f``.

    **Examples**

    >>> from sympy import sqf
    >>> from sympy.abc import x

    >>> sqf(2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16)
    2*(x + 1)**2*(x + 2)**3

    """
    return _generic_factor(f, gens, args, method='sqf')

def factor_list(f, *gens, **args):
    """
    Compute a list of irreducible factors of ``f``.

    **Examples**

    >>> from sympy import factor_list
    >>> from sympy.abc import x, y

    >>> factor_list(2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y)
    (2, [(x + y, 1), (x**2 + 1, 2)])

    """
    return _generic_factor_list(f, gens, args, method='factor')

def factor(f, *gens, **args):
    """
    Compute the factorization of ``f`` into irreducibles. (Use factorint to
    factor an integer.)

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
    2*(x + y)*(x**2 + 1)**2

    >>> factor(x**2 + 1)
    x**2 + 1
    >>> factor(x**2 + 1, modulus=2)
    (x + 1)**2
    >>> factor(x**2 + 1, gaussian=True)
    (x - I)*(x + I)

    >>> factor(x**2 - 2, extension=sqrt(2))
    (x - 2**(1/2))*(x + 2**(1/2))

    >>> factor((x**2 - 1)/(x**2 + 4*x + 4))
    (x - 1)*(x + 1)/(x + 2)**2
    >>> factor((x**2 + 4*x + 4)**10000000*(x**2 + 1))
    (x + 2)**20000000*(x**2 + 1)

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

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

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
    {0: 2, 1: 2}

    """
    options.allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('ground_roots', 1, exc)

    return F.ground_roots()

def nth_power_roots_poly(f, n, *gens, **args):
    """
    Construct a polynomial with n-th powers of roots of ``f``.

    **Examples**

    >>> from sympy import nth_power_roots_poly, factor, roots
    >>> from sympy.abc import x

    >>> f = x**4 - x**2 + 1
    >>> g = factor(nth_power_roots_poly(f, 2))

    >>> g
    (x**2 - x + 1)**2

    >>> R_f = [ (r**2).expand() for r in roots(f) ]
    >>> R_g = roots(g).keys()

    >>> set(R_f) == set(R_g)
    True

    """
    options.allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('nth_power_roots_poly', 1, exc)

    result = F.nth_power_roots_poly(n)

    if not opt.polys:
        return result.as_expr()
    else:
        return result

def cancel(f, *gens, **args):
    """
    Cancel common factors in a rational function ``f``.

    **Examples**

    >>> from sympy import cancel
    >>> from sympy.abc import x

    >>> cancel((2*x**2 - 2)/(x**2 - 2*x + 1))
    (2*x + 2)/(x - 1)

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
    ([2*x, 1], x**2 + y**2 + y)

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
    [x**2 - 2*y**2, x*y - 2*y, y**3 - 2*y]
    >>> groebner([x*y - 2*y, 2*y**2 - x**2], order='grlex')
    [y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y]
    >>> groebner([x*y - 2*y, 2*y**2 - x**2], order='grevlex')
    [x**3 - 2*x**2, -x**2 + 2*y**2, x*y - 2*y]

    **References**

    1. [Buchberger01]_
    2. [Cox97]_

    """
    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr(F, *gens, **args)
    except PolificationFailed, exc:
        raise ComputationFailed('groebner', len(F), exc)

    domain = opt.domain

    if domain.has_assoc_Field:
        opt.domain = domain.get_field()
    else:
        raise DomainError("can't compute a Groebner basis over %s" % domain)

    for i, poly in enumerate(polys):
        poly = poly.set_domain(opt.domain).rep.to_dict()
        polys[i] = sdp_from_dict(poly, opt.order)

    level = len(gens)-1

    G = sdp_groebner(polys, level, opt.order, opt.domain)
    G = [ Poly._from_dict(dict(g), opt) for g in G ]

    if not domain.has_Field:
        G = [ g.clear_denoms(convert=True)[1] for g in G ]

    if not opt.polys:
        return [ g.as_expr() for g in G ]
    else:
        return G

def poly(expr, *gens, **args):
    """
    Efficiently transform an expression into a polynomial.

    **Examples**

    >>> from sympy import poly
    >>> from sympy.abc import x

    >>> poly(x*(x**2 + x - 1)**2)
    Poly(x**5 + 2*x**4 - x**3 - 2*x**2 + x, x, domain='ZZ')

    """
    options.allowed_flags(args, [])

    def _poly(expr, opt):
        terms, poly_terms = [], []

        for term in Add.make_args(expr):
            factors, poly_factors = [], []

            for factor in Mul.make_args(term):
                if factor.is_Add:
                    poly_factors.append(_poly(factor, opt))
                elif factor.is_Pow and factor.base.is_Add and factor.exp.is_Integer:
                    poly_factors.append(_poly(factor.base, opt).pow(factor.exp))
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
                        product = product.mul(Poly._from_expr(factor, opt))

                poly_terms.append(product)

        if not poly_terms:
            result = Poly._from_expr(expr, opt)
        else:
            result = poly_terms[0]

            for term in poly_terms[1:]:
                result = result.add(term)

            if terms:
                term = Add(*terms)

                if term.is_Number:
                    result = result.add(term)
                else:
                    result = result.add(Poly._from_expr(term, opt))

        return result.reorder(**args)

    expr = sympify(expr)

    if expr.is_Poly:
        return Poly(expr, *gens, **args)

    if 'expand' not in args:
        args['expand'] = False

    opt = options.build_options(gens, args)

    return _poly(expr, opt)
