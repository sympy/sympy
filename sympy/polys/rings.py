"""Sparse polynomial rings. """

from __future__ import annotations

import sympy
from sympy.external.gmpy import GROUND_TYPES

from operator import add, mul, lt, le, gt, ge
from functools import reduce
from types import GeneratorType

from sympy.core.cache import cacheit
from sympy.core.expr import Expr
from sympy.core.intfunc import igcd
from sympy.core.symbol import Symbol, symbols as _symbols
from sympy.core.sympify import CantSympify, sympify
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.compatibility import IPolys
from sympy.polys.constructor import construct_domain
from sympy.polys.densebasic import ninf, dmp_to_dict, dmp_from_dict
from sympy.polys.domains import ZZ, QQ
from sympy.polys.domains.domain import Domain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.polynomialring import PolynomialRing
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.monomials import MonomialOps
from sympy.polys.orderings import lex, MonomialOrder, grlex, grevlex
from sympy.polys.polyerrors import (
    CoercionFailed, GeneratorsError,
    ExactQuotientFailed, MultivariatePolynomialError, DomainError)
from sympy.polys.polyoptions import (Domain as DomainOpt,
                                     Order as OrderOpt, build_options)
from sympy.polys.polyutils import (expr_from_dict, _dict_reorder,
                                   _parallel_dict_from_expr)
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public, subsets
from sympy.utilities.iterables import is_sequence
from sympy.utilities.magic import pollute


if GROUND_TYPES == 'flint':
    import flint
    def _supported_flint_domain(domain):
        return domain == ZZ or domain == QQ
else:
    flint = None
    def _supported_flint_domain(domain):
        return False


@public
def ring(symbols, domain, order: MonomialOrder|str = lex):
    """Construct a polynomial ring returning ``(ring, x_1, ..., x_n)``.

    Parameters
    ==========

    symbols : str
        Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~.Domain` or coercible
    order : :class:`~.MonomialOrder` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.orderings import lex

    >>> R, x, y, z = ring("x,y,z", ZZ, lex)
    >>> R
    Polynomial ring in x, y, z over ZZ with lex order
    >>> x + y + z
    x + y + z
    >>> type(_)
    <class 'sympy.polys.rings.PolyElement'>

    """
    _ring = PolyRing(symbols, domain, order)
    return (_ring,) + _ring.gens

@public
def xring(symbols, domain, order=lex):
    """Construct a polynomial ring returning ``(ring, (x_1, ..., x_n))``.

    Parameters
    ==========

    symbols : str
        Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~.Domain` or coercible
    order : :class:`~.MonomialOrder` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> from sympy.polys.rings import xring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.orderings import lex

    >>> R, (x, y, z) = xring("x,y,z", ZZ, lex)
    >>> R
    Polynomial ring in x, y, z over ZZ with lex order
    >>> x + y + z
    x + y + z
    >>> type(_)
    <class 'sympy.polys.rings.PolyElement'>

    """
    _ring = PolyRing(symbols, domain, order)
    return (_ring, _ring.gens)

@public
def vring(symbols, domain, order=lex):
    """Construct a polynomial ring and inject ``x_1, ..., x_n`` into the global namespace.

    Parameters
    ==========

    symbols : str
        Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~.Domain` or coercible
    order : :class:`~.MonomialOrder` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> from sympy.polys.rings import vring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.orderings import lex

    >>> vring("x,y,z", ZZ, lex)
    Polynomial ring in x, y, z over ZZ with lex order
    >>> x + y + z # noqa:
    x + y + z
    >>> type(_)
    <class 'sympy.polys.rings.PolyElement'>

    """
    _ring = PolyRing(symbols, domain, order)
    pollute([ sym.name for sym in _ring.symbols ], _ring.gens)
    return _ring

@public
def sring(exprs, *symbols, **options):
    """Construct a ring deriving generators and domain from options and input expressions.

    Parameters
    ==========

    exprs : :class:`~.Expr` or sequence of :class:`~.Expr` (sympifiable)
    symbols : sequence of :class:`~.Symbol`/:class:`~.Expr`
    options : keyword arguments understood by :class:`~.Options`

    Examples
    ========

    >>> from sympy import sring, symbols

    >>> x, y, z = symbols("x,y,z")
    >>> R, f = sring(x + 2*y + 3*z)
    >>> R
    Polynomial ring in x, y, z over ZZ with lex order
    >>> f
    x + 2*y + 3*z
    >>> type(_)
    <class 'sympy.polys.rings.PolyElement'>

    """
    single = False

    if not is_sequence(exprs):
        exprs, single = [exprs], True

    exprs = list(map(sympify, exprs))
    opt = build_options(symbols, options)

    # TODO: rewrite this so that it doesn't use expand() (see poly()).
    reps, opt = _parallel_dict_from_expr(exprs, opt)

    if opt.domain is None:
        coeffs = sum([ list(rep.values()) for rep in reps ], [])

        opt.domain, coeffs_dom = construct_domain(coeffs, opt=opt)

        coeff_map = dict(zip(coeffs, coeffs_dom))
        reps = [{m: coeff_map[c] for m, c in rep.items()} for rep in reps]

    _ring = PolyRing(opt.gens, opt.domain, opt.order)
    polys = list(map(_ring.from_dict, reps))

    if single:
        return (_ring, polys[0])
    else:
        return (_ring, polys)

def _parse_symbols(symbols):
    if isinstance(symbols, str):
        return _symbols(symbols, seq=True) if symbols else ()
    elif isinstance(symbols, Expr):
        return (symbols,)
    elif is_sequence(symbols):
        if all(isinstance(s, str) for s in symbols):
            return _symbols(symbols)
        elif all(isinstance(s, Expr) for s in symbols):
            return symbols

    raise GeneratorsError("expected a string, Symbol or expression or a non-empty sequence of strings, Symbols or expressions")

def _parse_symbols_flint(symbols):
    from sympy.core.symbol import Symbol

    if isinstance(symbols, str):
        if not symbols:
            return ()
        return tuple(symbols.replace(' ', '').split(','))

    if isinstance(symbols, Symbol):
        return (symbols.name,)

    if hasattr(symbols, '__iter__'):
        result = []
        for s in symbols:
            if isinstance(s, Symbol):
                result.append(s.name)
            elif isinstance(s, str):
                result.extend(s.replace(' ', '').split(','))
            else:
                raise ValueError(f"Invalid symbol: {s}")
        return tuple(result)

    raise ValueError(f"Cannot parse symbols: {symbols}")


class PolyRing(DefaultPrinting, IPolys):
    """Base class for multivariate distributed polynomial rings."""
    gens: tuple  # Will be set by subclasses
    symbols: tuple[Expr, ...]
    ngens: int
    domain: Domain
    order: MonomialOrder

    element_class = None
    dtype = None # Will be set by subclasses as a callable

    def __new__(cls, symbols, domain, order=lex):
        symbol_strings = tuple(_parse_symbols(symbols))
        #ngens = len(symbol_strings)
        preprocessed_domain = DomainOpt.preprocess(domain)
        preprocessed_order = OrderOpt.preprocess(order)

        if preprocessed_domain.is_Composite and set(symbol_strings) & set(preprocessed_domain.symbols):
            raise GeneratorsError("polynomial ring and its ground domain share generators")

        if flint is not None and _supported_flint_domain(preprocessed_domain):
            subclass = FlintPolyRing
        else:
            subclass = PythonPolyRing

        obj = object.__new__(subclass)
        obj._initialized = False
        obj.__init__(symbols, preprocessed_domain, preprocessed_order)
        obj._initialized = True
        return obj

    def __init__(self, symbols, domain, order=lex):
        if getattr(self, '_initialized', False):
            return

        self.symbols = tuple(_parse_symbols(symbols))
        self.ngens = len(self.symbols)
        self.domain = domain
        self.order = order

        self._hash_tuple = (self.__class__.__name__, self.symbols, self.ngens, self.domain, self.order)
        self._hash = hash(self._hash_tuple)

        self.zero_monom = (0,) * self.ngens

    def __hash__(self):
        return self._hash

    def _gens(self):
        raise NotImplementedError("Subclasses must implement _gens")

    def __getnewargs__(self):
        return (self.symbols, self.domain, self.order)

    def __getstate__(self):
        raise NotImplementedError("Subclasses must implement __getstate__")

    def __eq__(self, other):
        raise NotImplementedError("Subclasses must implement __eq__")

    def __ne__(self, other):
        return not self == other

    @cacheit
    def _clone(self, symbols, domain, order):
        return self.__class__(symbols or self.symbols, domain or self.domain, order or self.order)

    def clone(self, symbols=None, domain=None, order=None):
        if symbols is not None and isinstance(symbols, list):
            symbols = tuple(symbols)
        return self._clone(symbols, domain, order)

    def monomial_basis(self, i):
        """Return the ith-basis element."""
        basis = [0] * self.ngens
        basis[i] = 1
        return tuple(basis)

    @property
    def zero(self):
        raise NotImplementedError("Subclasses must implement zero")

    @property
    def one(self):
        raise NotImplementedError("Subclasses must implement one")

    def is_element(self, element):
        raise NotImplementedError("Subclasses must implement is_element")

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
        raise NotImplementedError("Subclasses must implement ring_new")

    __call__ = ring_new

    def from_dict(self, element, orig_domain=None):
        raise NotImplementedError("subclasses must implement from_dict")

    def from_terms(self, element, orig_domain=None):
        return self.from_dict(dict(element), orig_domain)

    def from_list(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1, self.domain))

    def from_expr(self, expr):
        raise NotImplementedError("subclasses must implement from_expr")

    def index(self, gen):
        raise NotImplementedError("subclasses must implement index")

    def drop(self, *gens):
        raise NotImplementedError("subclasses must implement drop")

    def __getitem__(self, key):
        symbols = self.symbols[key]
        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def to_ground(self):
        if self.domain.is_Composite or hasattr(self.domain, 'domain'):
            return self.clone(domain=self.domain.domain)
        else:
            raise ValueError("%s is not a composite domain" % self.domain)

    def to_domain(self):
        return PolynomialRing(self)

    def to_field(self):
        from sympy.polys.fields import FracField
        return FracField(self.symbols, self.domain, self.order)

    @property
    def is_univariate(self):
        raise NotImplementedError("subclasses must implement is_univariate")

    @property
    def is_multivariate(self):
        raise NotImplementedError("subclasses must implement is_multivariate")

    def compose(self, other):
        raise NotImplementedError("subclasses must implement compose")

    def add_gens(self, symbols):
        raise NotImplementedError("subclasses must implement add_gens")

    def symmetric_poly(self, n):
        raise NotImplementedError("subclasses must implement symmetric_poly")


class PolyElement(DomainElement, DefaultPrinting, CantSympify):
    """Base class for elements of multivariate distributed polynomial rings."""

    def __init__(self, ring):
        super().__init__()
        self.ring = ring

    def _check(self):
        raise NotImplementedError

    def new(self, init):
        """Create a new polynomial element of the same type."""
        raise NotImplementedError("Subclasses must implement new()")

    def parent(self):
        return self.ring.to_domain()

    def __getnewargs__(self):
        return (self.ring, list(self.iterterms()))

    def __hash__(self):
        raise NotImplementedError("Subclasses must implement __hash__")

    def copy(self):
        raise NotImplementedError("Subclasses must implement copy")

    def set_ring(self, new_ring):
        raise NotImplementedError("Subclasses must implement set_ring")

    def as_expr(self, *symbols):
        if not symbols:
            symbols = self.ring.symbols
        elif len(symbols) != self.ring.ngens:
            raise ValueError(
                "Wrong number of symbols, expected %s got %s" %
                (self.ring.ngens, len(symbols))
            )

        return expr_from_dict(self.as_expr_dict(), *symbols)

    def as_expr_dict(self):
        to_sympy = self.ring.domain.to_sympy
        return {monom: to_sympy(coeff) for monom, coeff in self.iterterms()}

    def clear_denoms(self):
        raise NotImplementedError("Subclasses must implement clear_denoms")

    def strip_zero(self):
        raise NotImplementedError("Subclasses must implement strip_zero")

    def __eq__(self, other):
        raise NotImplementedError("Subclasses must implement __eq__")

    def __ne__(self, other):
        return not self == other

    def almosteq(self, other, tolerance=None):
        raise NotImplementedError("Subclasses must implement almosteq")

    def sort_key(self):
        raise NotImplementedError("Subclasses must implement sort_key")

    def _cmp(self, other, op):
        if self.ring.is_element(other):
            return op(self.sort_key(), other.sort_key())
        return NotImplemented

    def __lt__(self, other):
        return self._cmp(other, lt)

    def __le__(self, other):
        return self._cmp(other, le)

    def __gt__(self, other):
        return self._cmp(other, gt)

    def __ge__(self, other):
        return self._cmp(other, ge)

    def drop(self, gen):
        raise NotImplementedError("Subclasses must implement drop")

    def drop_to_ground(self, gen):
        raise NotImplementedError("Subclasses must implement drop_to_ground")

    def to_dense(self):
        raise NotImplementedError("Subclasses must implement to_dense")

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        if not self:
            return printer._print(self.ring.domain.zero)
        prec_mul = precedence["Mul"]
        prec_atom = precedence["Atom"]
        ring = self.ring
        symbols = ring.symbols
        ngens = ring.ngens
        zm = ring.zero_monom
        sexpvs = []
        for expv, coeff in self.terms():
            negative = ring.domain.is_negative(coeff)
            sign = " - " if negative else " + "
            sexpvs.append(sign)
            if expv == zm:
                scoeff = printer._print(coeff)
                if negative and scoeff.startswith("-"):
                    scoeff = scoeff[1:]
            else:
                if negative:
                    coeff = -coeff
                if coeff != self.ring.domain.one:
                    scoeff = printer.parenthesize(coeff, prec_mul, strict=True)
                else:
                    scoeff = ''
            sexpv = []
            for i in range(ngens):
                exp = expv[i]
                if not exp:
                    continue
                symbol = printer.parenthesize(symbols[i], prec_atom, strict=True)
                if exp != 1:
                    if exp != int(exp) or exp < 0:
                        sexp = printer.parenthesize(exp, prec_atom, strict=False)
                    else:
                        sexp = exp
                    sexpv.append(exp_pattern % (symbol, sexp))
                else:
                    sexpv.append('%s' % symbol)
            if scoeff:
                sexpv = [scoeff] + sexpv
            sexpvs.append(mul_symbol.join(sexpv))
        if sexpvs[0] in [" + ", " - "]:
            head = sexpvs.pop(0)
            if head == " - ":
                sexpvs.insert(0, "-")
        return "".join(sexpvs)

    @property
    def is_generator(self):
        return self in self.ring._gens_set

    @property
    def is_ground(self):
        raise NotImplementedError("Subclasses must implement is_ground")

    @property
    def is_monomial(self):
        raise NotImplementedError("Subclasses must implement is_monomial")

    @property
    def is_term(self):
        return NotImplementedError("Subclasses must implement is_term")

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

    @property
    def is_zero(self):
        raise NotImplementedError("Subclasses must implement is_zero")

    @property
    def is_one(self):
        return self == self.ring.one

    @property
    def is_monic(self):
        return self.ring.domain.is_one(self.LC)

    @property
    def is_primitive(self):
        return self.ring.domain.is_one(self.content())

    @property
    def is_linear(self):
        return all(sum(monom) <= 1 for monom in self.itermonoms())

    @property
    def is_quadratic(self):
        return all(sum(monom) <= 2 for monom in self.itermonoms())

    @property
    def is_squarefree(self):
        raise NotImplementedError("Subclasses must implement is_squarefree")

    @property
    def is_irreducible(self):
        raise NotImplementedError("Subclasses must implement is_irreducible")

    @property
    def is_cyclotomic(self):
        raise NotImplementedError("Subclasses must implement is_cyclotomic")

    def __neg__(self):
        return (-1) * self

    def __pos__(self):
        return self

    def __add__(self, other):
        raise NotImplementedError("Subclasses must implement __add__")

    def __radd__(self, other):
        raise NotImplementedError("Subclasses must implement __radd__")

    def __sub__(self, other):
        raise NotImplementedError("Subclasses must implement __sub__")

    def __rsub__(self, other):
        raise NotImplementedError("Subclasses must implement __rsub__")

    def __mul__(self, other):
        raise NotImplementedError("Subclasses must implement __mul__")

    def __rmul__(self, other):
        raise NotImplementedError("Subclasses must implement __rmul__")

    def __pow__(self, n):
        raise NotImplementedError("Subclasses must implement __pow__")

    def square(self):
        raise NotImplementedError("Subclasses must implement square")

    def __divmod__(self, other):
        raise NotImplementedError("Subclasses must implement __divmod__")

    def __rdivmod__(self, other):
        raise NotImplementedError("Subclasses must implement __rdivmod__")

    def __mod__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1.rem(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rmod__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1.rem_ground(p2)

    def __rmod__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2.rem(p1)

    def __floordiv__(self, other):
        raise NotImplementedError("Subclasses must implement __floordiv__")

    def __truediv__(self, other):
        raise NotImplementedError("Subclasses must implement __truediv__")

    def __rtruediv__(self, other):
        raise NotImplementedError("Subclasses must implement __rtruediv__")

    def div(self, fv):
        raise NotImplementedError("Subclasses must implement div")

    def rem(self, G):
        raise NotImplementedError("Subclasses must implement rem")

    def quo(self, G):
        return self.div(G)[0]

    def exquo(self, G):
        q, r = self.div(G)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(self, G)

    def _iadd_monom(self, mc):
        raise NotImplementedError("Subclasses must implement _iadd_monom")

    def _iadd_poly_monom(self, p2, mc):
        raise NotImplementedError("Subclasses must implement _iadd_poly_monom")

    def degree(self):
        raise NotImplementedError("Subclasses must implement degree")

    def degrees(self):
        raise NotImplementedError("Subclasses must implement degrees")

    def tail_degree(self, x=None):
        raise NotImplementedError("Subclasses must implement tail_degree")

    def tail_degrees(self):
        raise NotImplementedError("Subclasses must implement tail_degrees")

    def coeff(self, element):
        raise NotImplementedError("Subclasses must implement coeff")

    def const(self):
        raise NotImplementedError("Subclasses must implement const")

    @property
    def LC(self):
        raise NotImplementedError("Subclasses must implement LC")

    @property
    def LM(self):
        raise NotImplementedError("Subclasses must implement LM")

    def leading_monom(self):
        raise NotImplementedError("Subclasses must implement leading_monom")

    @property
    def LT(self):
        raise NotImplementedError("Subclasses must implement LT")

    def leading_term(self):
        raise NotImplementedError("Subclasses must implement leading_term")

    def _sorted(self, seq, order):
        if order is None:
            order = self.ring.order
        else:
            order = OrderOpt.preprocess(order)

        if order is lex:
            return sorted(seq, key=lambda monom: monom[0], reverse=True)
        else:
            return sorted(seq, key=lambda monom: order(monom[0]), reverse=True)

    def coeffs(self, order=None):
        raise NotImplementedError("Subclasses must implement coeffs")

    def monoms(self, order=None):
        raise NotImplementedError("Subclasses must implement monoms")

    def terms(self, order=None):
        raise NotImplementedError("Subclasses must implement terms")

    def itercoeffs(self):
        """Iterator over coefficients of a polynomial. """
        raise NotImplementedError("Subclasses must implement itercoeffs")

    def itermonoms(self):
        """Iterator over monomials of a polynomial. """
        raise NotImplementedError("Subclasses must implement itermonoms")

    def iterterms(self):
        """Iterator over terms of a polynomial. """
        raise NotImplementedError("Subclasses must implement iterterms")

    def listcoeffs(self):
        """Unordered list of polynomial coefficients. """
        raise NotImplementedError("Subclasses must implement listcoeffs")

    def listmonoms(self):
        """Unordered list of polynomial monomials. """
        raise NotImplementedError("Subclasses must implement listmonoms")

    def listterms(self):
        """Unordered list of polynomial terms. """
        raise NotImplementedError("Subclasses must implement listterms")

    def imul_num(selfself, c):
        raise NotImplementedError("Subclasses must implement imul_num")

    def content(self):
        raise NotImplementedError("Subclasses must implement content")

    def primitive(self):
        raise NotImplementedError("Subclasses must implement primitive")

    def monic(self):
        raise NotImplementedError("Subclasses must implement monic")

    def mul_ground(self, x):
        raise NotImplementedError("Subclasses must implement mul_ground")

    def mul_monom(self, monom):
        raise NotImplementedError("Subclasses must implement mul_monom")

    def mul_term(self, term):
        raise NotImplementedError("Subclasses must implement mul_term")

    def quo_ground(self, x):
        raise NotImplementedError("Subclasses must implement quo_ground")

    def quo_term(self, term):
        raise NotImplementedError("Subclasses must implement quo_term")

    def trunc_ground(self, p):
        raise NotImplementedError("Subclasses must implement trunc_ground")

    rem_ground = trunc_ground

    def extract_ground(self, g):
        raise NotImplementedError("Subclasses must implement extract_ground")

    def max_norm(self):
        raise NotImplementedError("Subclasses must implement max_norm")

    def li_norm(self):
        raise NotImplementedError("Subclasses must implement li_norm")

    def deflate(self, *G):
        raise NotImplementedError("Subclasses must implement deflate")

    def inflate(self, J):
        return NotImplementedError("Subclasses must implement inflate")

    def lcm(self, g):
        raise NotImplementedError("Subclasses must implement lcm")

    def gcd(self, other):
        raise NotImplementedError("Subclasses must implement gcd")

    def cofactors(self, other):
        """Return the cofactors of self and other."""
        raise NotImplementedError("Subclasses must implement cofactors")

    def cancel(self, g):
        raise NotImplementedError("Subclasses must implement cancel")

    def canonical_unit(self):
        domain = self.ring.domain
        return domain.canonical_unit(self.LC)

    def diff(self, x):
        raise NotImplementedError("Subclasses must implement diff")

    def __call__(self, *values):
        if 0 < len(values) <= self.ring.ngens:
            return self.evaluate(list(zip(self.ring.gens, values)))
        else:
            raise ValueError("expected at least 1 and at most %s values, got %s" % (self.ring.ngens, len(values)))

    def evaluate(self, x, a=None):
        raise NotImplementedError("Subclasses must implement evaluate")

    def subs(self, x, a=None):
        raise NotImplementedError("Subclasses must implement subs")

    def symmetrize(self):
        raise NotImplementedError("Subclasses must implement symmetrize")

    def compose(self, x, a=None):
        raise NotImplementedError("Subclasses must implement compose")

    def coeff_wrt(self, x, deg):
        raise NotImplementedError("Subclasses must implement coeff_wrt")

    def prem(self, g, x=None):
        """
        Pseudo-remainder of the polynomial ``self`` with respect to ``g``.

        The pseudo-quotient ``q`` and pseudo-remainder ``r`` with respect to
        ``z`` when dividing ``f`` by ``g`` satisfy ``m*f = g*q + r``,
        where ``deg(r,z) < deg(g,z)`` and
        ``m = LC(g,z)**(deg(f,z) - deg(g,z)+1)``.

        See :meth:`pdiv` for explanation of pseudo-division.


        Parameters
        ==========

        g : :py:class:`~.PolyElement`
            The polynomial to divide ``self`` by.
        x : generator or generator index, optional
            The main variable of the polynomials and default is first generator.

        Returns
        =======

        :py:class:`~.PolyElement`
            The pseudo-remainder polynomial.

        Raises
        ======

        ZeroDivisionError : If ``g`` is the zero polynomial.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x, y = ring("x, y", ZZ)

        >>> f = x**2 + x*y
        >>> g = 2*x + 2
        >>> f.prem(g) # first generator is chosen by default if it is not given
        -4*y + 4
        >>> f.rem(g) # shows the difference between prem and rem
        x**2 + x*y
        >>> f.prem(g, y) # generator is given
        0
        >>> f.prem(g, 1) # generator index is given
        0

        See Also
        ========

        pdiv, pquo, pexquo, sympy.polys.domains.ring.Ring.rem

        """
        f = self
        x = f.ring.index(x)
        df = f.degree(x)
        dg = g.degree(x)

        if dg < 0:
            raise ZeroDivisionError('polynomial division')

        r, dr = f, df

        if df < dg:
            return r

        N = df - dg + 1

        lc_g = g.coeff_wrt(x, dg)

        xp = f.ring.gens[x]

        while True:

            lc_r = r.coeff_wrt(x, dr)
            j, N = dr - dg, N - 1

            R = r * lc_g
            G = g * lc_r * xp ** j
            r = R - G

            dr = r.degree(x)

            if dr < dg:
                break

        c = lc_g ** N

        return r * c

    def pdiv(self, g, x=None):
        """
        Computes the pseudo-division of the polynomial ``self`` with respect to ``g``.

        The pseudo-division algorithm is used to find the pseudo-quotient ``q``
        and pseudo-remainder ``r`` such that ``m*f = g*q + r``, where ``m``
        represents the multiplier and ``f`` is the dividend polynomial.

        The pseudo-quotient ``q`` and pseudo-remainder ``r`` are polynomials in
        the variable ``x``, with the degree of ``r`` with respect to ``x``
        being strictly less than the degree of ``g`` with respect to ``x``.

        The multiplier ``m`` is defined as
        ``LC(g, x) ^ (deg(f, x) - deg(g, x) + 1)``,
        where ``LC(g, x)`` represents the leading coefficient of ``g``.

        It is important to note that in the context of the ``prem`` method,
        multivariate polynomials in a ring, such as ``R[x,y,z]``, are treated
        as univariate polynomials with coefficients that are polynomials,
        such as ``R[x,y][z]``. When dividing ``f`` by ``g`` with respect to the
        variable ``z``, the pseudo-quotient ``q`` and pseudo-remainder ``r``
        satisfy ``m*f = g*q + r``, where ``deg(r, z) < deg(g, z)``
        and ``m = LC(g, z)^(deg(f, z) - deg(g, z) + 1)``.

        In this function, the pseudo-remainder ``r`` can be obtained using the
        ``prem`` method, the pseudo-quotient ``q`` can
        be obtained using the ``pquo`` method, and
        the function ``pdiv`` itself returns a tuple ``(q, r)``.


        Parameters
        ==========

        g : :py:class:`~.PolyElement`
            The polynomial to divide ``self`` by.
        x : generator or generator index, optional
            The main variable of the polynomials and default is first generator.

        Returns
        =======

        :py:class:`~.PolyElement`
            The pseudo-division polynomial (tuple of ``q`` and ``r``).

        Raises
        ======

        ZeroDivisionError : If ``g`` is the zero polynomial.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x, y = ring("x, y", ZZ)

        >>> f = x**2 + x*y
        >>> g = 2*x + 2
        >>> f.pdiv(g) # first generator is chosen by default if it is not given
        (2*x + 2*y - 2, -4*y + 4)
        >>> f.div(g) # shows the difference between pdiv and div
        (0, x**2 + x*y)
        >>> f.pdiv(g, y) # generator is given
        (2*x**3 + 2*x**2*y + 6*x**2 + 2*x*y + 8*x + 4, 0)
        >>> f.pdiv(g, 1) # generator index is given
        (2*x**3 + 2*x**2*y + 6*x**2 + 2*x*y + 8*x + 4, 0)

        See Also
        ========

        prem
            Computes only the pseudo-remainder more efficiently than
            `f.pdiv(g)[1]`.
        pquo
            Returns only the pseudo-quotient.
        pexquo
            Returns only an exact pseudo-quotient having no remainder.
        div
            Returns quotient and remainder of f and g polynomials.

        """
        f = self
        x = f.ring.index(x)

        df = f.degree(x)
        dg = g.degree(x)

        if dg < 0:
            raise ZeroDivisionError("polynomial division")

        q, r, dr = x, f, df

        if df < dg:
            return q, r

        N = df - dg + 1
        lc_g = g.coeff_wrt(x, dg)

        xp = f.ring.gens[x]

        while True:

            lc_r = r.coeff_wrt(x, dr)
            j, N = dr - dg, N - 1

            Q = q * lc_g

            q = Q + (lc_r) * xp ** j

            R = r * lc_g

            G = g * lc_r * xp ** j

            r = R - G

            dr = r.degree(x)

            if dr < dg:
                break

        c = lc_g ** N

        q = q * c
        r = r * c

        return q, r

    def pquo(self, g, x=None):
        """
        Polynomial pseudo-quotient in multivariate polynomial ring.

        Examples
        ========
        >>> from sympy.polys import ring, ZZ
        >>> R, x,y = ring("x,y", ZZ)

        >>> f = x**2 + x*y
        >>> g = 2*x + 2*y
        >>> h = 2*x + 2
        >>> f.pquo(g)
        2*x
        >>> f.quo(g) # shows the difference between pquo and quo
        0
        >>> f.pquo(h)
        2*x + 2*y - 2
        >>> f.quo(h) # shows the difference between pquo and quo
        0

        See Also
        ========

        prem, pdiv, pexquo, sympy.polys.domains.ring.Ring.quo

        """
        f = self
        return f.pdiv(g, x)[0]

    def pexquo(self, g, x=None):
        """
        Polynomial exact pseudo-quotient in multivariate polynomial ring.

        Examples
        ========
        >>> from sympy.polys import ring, ZZ
        >>> R, x,y = ring("x,y", ZZ)

        >>> f = x**2 + x*y
        >>> g = 2*x + 2*y
        >>> h = 2*x + 2
        >>> f.pexquo(g)
        2*x
        >>> f.exquo(g) # shows the difference between pexquo and exquo
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2*x + 2*y does not divide x**2 + x*y
        >>> f.pexquo(h)
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2*x + 2 does not divide x**2 + x*y

        See Also
        ========

        prem, pdiv, pquo, sympy.polys.domains.ring.Ring.exquo

        """
        f = self
        q, r = f.pdiv(g, x)

        if r.is_zero:
            return q
        else:
            raise ExactQuotientFailed(f, g)

    def subresultants(self, g, x=None):
        """
        Computes the subresultant PRS of two polynomials ``self`` and ``g``.

        Parameters
        ==========

        g : :py:class:`~.PolyElement`
            The second polynomial.
        x : generator or generator index
            The variable with respect to which the subresultant sequence is computed.

        Returns
        =======

        R : list
            Returns a list polynomials representing the subresultant PRS.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x, y = ring("x, y", ZZ)

        >>> f = x**2*y + x*y
        >>> g = x + y
        >>> f.subresultants(g) # first generator is chosen by default if not given
        [x**2*y + x*y, x + y, y**3 - y**2]
        >>> f.subresultants(g, 0) # generator index is given
        [x**2*y + x*y, x + y, y**3 - y**2]
        >>> f.subresultants(g, y) # generator is given
        [x**2*y + x*y, x + y, x**3 + x**2]

        """
        f = self
        x = f.ring.index(x)
        n = f.degree(x)
        m = g.degree(x)

        if n < m:
            f, g = g, f
            n, m = m, n

        if f == 0:
            return [0, 0]

        if g == 0:
            return [f, 1]

        R = [f, g]

        d = n - m
        b = (-1) ** (d + 1)

        # Compute the pseudo-remainder for f and g
        h = f.prem(g, x)
        h = h * b

        # Compute the coefficient of g with respect to x**m
        lc = g.coeff_wrt(x, m)

        c = lc ** d

        S = [1, c]

        c = -c

        while h:
            k = h.degree(x)

            R.append(h)
            f, g, m, d = g, h, k, m - k

            b = -lc * c ** d
            h = f.prem(g, x)
            h = h.exquo(b)

            lc = g.coeff_wrt(x, k)

            if d > 1:
                p = (-lc) ** d
                q = c ** (d - 1)
                c = p.exquo(q)
            else:
                c = -lc

            S.append(-c)

        return R

    # TODO: following methods should point to polynomial
    # representation independent algorithm implementations.

    def half_gcdex(f, g):
        raise NotImplementedError

    def gcdex(f, g):
        raise NotImplementedError

    def resultant(f, g):
        raise NotImplementedError("Subclasses must implement resultant")

    def discriminant(f):
        raise NotImplementedError("Subclasses must implement discriminant")

    def decompose(f):
        raise NotImplementedError("Subclasses must implement decompose")

    def shift(f, a):
        raise NotImplementedError("Subclasses must implement shift")

    def shift_list(f, a):
        raise NotImplementedError("Subclasses must implement shift_list")

    def sturm(f):
        raise NotImplementedError("Subclasses must implement sturm")

    def gff_list(f):
        raise NotImplementedError("Subclasses must implement gff_list")

    def norm(f):
        raise NotImplementedError("Subclasses must implement norm")

    def sqf_norm(f):
        raise NotImplementedError("Subclasses must implement sqf_norm")

    def sqf_part(f):
        raise NotImplementedError("Subclasses must implement sqf_part")

    def sqf_list(f, all=False):
        raise NotImplementedError("Subclasses must implement sqf_list")

    def factor_list(f):
        raise NotImplementedError("Subclasses must implement factor_list")


class PythonPolyElement(PolyElement, dict):

    def __init__(self, ring, init):
        super().__init__(ring)
        dict.__init__(self, init)

    def _check(self):
        assert isinstance(self, PolyElement)
        assert isinstance(self.ring, PolyRing)
        dom = self.ring.domain
        assert isinstance(dom, Domain)
        for monom, coeff in self.items():
            assert dom.of_type(coeff)
            assert len(monom) == self.ring.ngens
            assert all(isinstance(exp, int) and exp >= 0 for exp in monom)

    def new(self, init):
        return self.__class__(self.ring, init)

    def __getnewargs__(self):
        return (self.ring, list(self.iterterms()))

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
            terms = list(zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols)))
            return new_ring.from_terms(terms, self.ring.domain)
        else:
            return new_ring.from_dict(self, self.ring.domain)

    def clear_denoms(self):
        domain = self.ring.domain

        if not domain.is_Field or not domain.has_assoc_Ring:
            return domain.one, self

        ground_ring = domain.get_ring()
        common = ground_ring.one
        lcm = ground_ring.lcm
        denom = domain.denom

        for coeff in self.values():
            common = lcm(common, denom(coeff))

        poly = self.new([ (k, v*common) for k, v in self.items() ])
        return common, poly

    def strip_zero(self):
        """Eliminate monomials with zero coefficient. """
        for k, v in list(self.items()):
            if not v:
                del self[k]

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
        elif p1.ring.is_element(p2):
            return dict.__eq__(p1, p2)
        elif len(p1) > 1:
            return False
        else:
            return p1.get(p1.ring.zero_monom) == p2

    def almosteq(p1, p2, tolerance=None):
        """Approximate equality test for polynomials. """
        ring = p1.ring

        if ring.is_element(p2):
            if set(p1.keys()) != set(p2.keys()):
                return False

            almosteq = ring.domain.almosteq

            for k in p1.keys():
                if not almosteq(p1[k], p2[k], tolerance):
                    return False
            return True
        elif len(p1) > 1:
            return False
        else:
            try:
                p2 = ring.domain.convert(p2)
            except CoercionFailed:
                return False
            else:
                return ring.domain.almosteq(p1.const(), p2, tolerance)

    def sort_key(self):
        return (len(self), self.terms())

    def _drop(self, gen):
        ring = self.ring
        i = ring.index(gen)

        if ring.ngens == 1:
            return i, ring.domain
        else:
            symbols = list(ring.symbols)
            del symbols[i]
            return i, ring.clone(symbols=symbols)

    def drop(self, gen):
        i, ring = self._drop(gen)

        if self.ring.ngens == 1:
            if self.is_ground:
                return self.coeff(1)
            else:
                raise ValueError("Cannot drop %s" % gen)
        else:
            poly = ring.zero

            for k, v in self.items():
                if k[i] == 0:
                    K = list(k)
                    del K[i]
                    poly[tuple(K)] = v
                else:
                    raise ValueError("Cannot drop %s" % gen)

            return poly

    def _drop_to_ground(self, gen):
        ring = self.ring
        i = ring.index(gen)

        symbols = list(ring.symbols)
        del symbols[i]
        return i, ring.clone(symbols=symbols, domain=ring[i])

    def drop_to_ground(self, gen):
        if self.ring.ngens == 1:
            raise ValueError("Cannot drop only generator to ground")

        i, ring = self._drop_to_ground(gen)
        poly = ring.zero
        gen = ring.domain.gens[0]

        for monom, coeff in self.iterterms():
            mon = monom[:i] + monom[i+1:]
            if mon not in poly:
                poly[mon] = (gen**monom[i]).mul_ground(coeff)
            else:
                poly[mon] += (gen**monom[i]).mul_ground(coeff)

        return poly

    def to_dense(self):
        return dmp_from_dict(self, self.ring.ngens-1, self.ring.domain)

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
    def is_zero(f):
        return not f

    @property
    def is_squarefree(f):
        if not f.ring.ngens:
            return True
        return f.ring.dmp_sqf_p(f)

    @property
    def is_irreducible(f):
        if not f.ring.ngens:
            return True
        return f.ring.dmp_irreducible_p(f)

    @property
    def is_cyclotomic(f):
        if f.ring.is_univariate:
            return f.ring.dup_cyclotomic_p(f)
        else:
            raise MultivariatePolynomialError("cyclotomic polynomial")

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
        if ring.is_element(p2):
            p = p1.copy()
            get = p.get
            zero = ring.domain.zero
            for k, v in p2.items():
                v = get(k, zero) + v
                if v:
                    p[k] = v
                else:
                    del p[k]
            return p
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
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
            if zm not in p1.keys():
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
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            zm = ring.zero_monom
            if zm not in p1.keys():
                p[zm] = n
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
        if ring.is_element(p2):
            p = p1.copy()
            get = p.get
            zero = ring.domain.zero
            for k, v in p2.items():
                v = get(k, zero) - v
                if v:
                    p[k] = v
                else:
                    del p[k]
            return p
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
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
            if zm not in p1.keys():
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
        ring = p1.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            p = ring.zero
            for expv in p1:
                p[expv] = -p1[expv]
            p += n
            # p._check()
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
        if not p1 or not p2:
            return p
        elif ring.is_element(p2):
            get = p.get
            zero = ring.domain.zero
            monomial_mul = ring.monomial_mul
            p2it = list(p2.items())
            for exp1, v1 in p1.items():
                for exp2, v2 in p2it:
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, zero) + v1*v2
            p.strip_zero()
            # p._check()
            return p
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
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
            for exp1, v1 in p1.items():
                v = v1*p2
                if v:
                    p[exp1] = v
            # p._check()
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
        if not p2:
            return p
        try:
            p2 = p.ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            for exp1, v1 in p1.items():
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
        if not isinstance(n, int):
            raise TypeError("exponent must be an integer, got %s" % n)
        elif n < 0:
            raise ValueError("exponent must be a non-negative integer, got %s" % n)

        ring = self.ring

        if not n:
            if self:
                return ring.one
            else:
                raise ValueError("0**0")
        elif len(self) == 1:
            monom, coeff = list(self.items())[0]
            p = ring.zero
            if coeff == ring.domain.one:
                p[ring.monomial_pow(monom, n)] = coeff
            else:
                p[ring.monomial_pow(monom, n)] = coeff**n
            # p._check()
            return p

        # For ring series, we need negative and rational exponent support only
        # with monomials.
        n = int(n)
        if n < 0:
            raise ValueError("Negative exponent")

        elif n == 1:
            return self.copy()
        elif n == 2:
            return self.square()
        elif n == 3:
            return self*self.square()
        elif len(self) <= 5: # TODO: use an actual density measure
            return self._pow_multinomial(n)
        else:
            return self._pow_generic(n)

    def _pow_generic(self, n):
        p = self.ring.one
        c = self

        while True:
            if n & 1:
                p = p*c
                n -= 1
                if not n:
                    break

            c = c.square()
            n = n // 2

        return p

    def _pow_multinomial(self, n):
        multinomials = multinomial_coefficients(len(self), n).items()
        monomial_mulpow = self.ring.monomial_mulpow
        zero_monom = self.ring.zero_monom
        terms = self.items()
        zero = self.ring.domain.zero
        poly = self.ring.zero

        for multinomial, multinomial_coeff in multinomials:
            product_monom = zero_monom
            product_coeff = multinomial_coeff

            for exp, (monom, coeff) in zip(multinomial, terms):
                if exp:
                    product_monom = monomial_mulpow(product_monom, monom, exp)
                    product_coeff *= coeff**exp

            monom = tuple(product_monom)
            coeff = product_coeff

            coeff = poly.get(monom, zero) + coeff

            if coeff:
                poly[monom] = coeff
            elif monom in poly:
                del poly[monom]

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
        keys = list(self.keys())
        zero = ring.domain.zero
        monomial_mul = ring.monomial_mul
        for i in range(len(keys)):
            k1 = keys[i]
            pk = self[k1]
            for j in range(i):
                k2 = keys[j]
                exp = monomial_mul(k1, k2)
                p[exp] = get(exp, zero) + pk*self[k2]
        p = p.imul_num(2)
        get = p.get
        for k, v in self.items():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, zero) + v**2
        p.strip_zero()
        # p._check()
        return p

    def __divmod__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1.div(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rdivmod__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return (p1.quo_ground(p2), p1.rem_ground(p2))

    def __rdivmod__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2.div(p1)

    def __floordiv__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1.quo(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rtruediv__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1.quo_ground(p2)

    def __rfloordiv__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2.quo(p1)

    def __truediv__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1.exquo(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif isinstance(p2.ring.domain, PolynomialRing) and p2.ring.domain.ring == ring:
                return p2.__rtruediv__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1.quo_ground(p2)

    def __rtruediv__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2.exquo(p1)

    def _term_div(self):
        zm = self.ring.zero_monom
        domain = self.ring.domain
        domain_quo = domain.quo
        monomial_div = self.ring.monomial_div

        if domain.is_Field:
            def term_div(a_lm_a_lc, b_lm_b_lc):
                a_lm, a_lc = a_lm_a_lc
                b_lm, b_lc = b_lm_b_lc
                if b_lm == zm: # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if monom is not None:
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return None
        else:
            def term_div(a_lm_a_lc, b_lm_b_lc):
                a_lm, a_lc = a_lm_a_lc
                b_lm, b_lc = b_lm_b_lc
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

        """
        ring = self.ring
        ret_single = False
        if isinstance(fv, PolyElement):
            ret_single = True
            fv = [fv]
        if not all(fv):
            raise ZeroDivisionError("polynomial division")
        if not self:
            if ret_single:
                return ring.zero, ring.zero
            else:
                return [], ring.zero
        for f in fv:
            if f.ring != ring:
                raise ValueError('self and f must have the same ring')
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

    def rem(self, G):
        f = self
        if isinstance(G, PolyElement):
            G = [G]
        if not all(G):
            raise ZeroDivisionError("polynomial division")
        ring = f.ring
        domain = ring.domain
        zero = domain.zero
        monomial_mul = ring.monomial_mul
        r = ring.zero
        term_div = f._term_div()
        ltf = f.LT
        f = f.copy()
        get = f.get
        while f:
            for g in G:
                tq = term_div(ltf, g.LT)
                if tq is not None:
                    m, c = tq
                    for mg, cg in g.iterterms():
                        m1 = monomial_mul(mg, m)
                        c1 = get(m1, zero) - c*cg
                        if not c1:
                            del f[m1]
                        else:
                            f[m1] = c1
                    ltm = f.leading_expv()
                    if ltm is not None:
                        ltf = ltm, f[ltm]

                    break
            else:
                ltm, ltc = ltf
                if ltm in r:
                    r[ltm] += ltc
                else:
                    r[ltm] = ltc
                del f[ltm]
                ltm = f.leading_expv()
                if ltm is not None:
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
        if self in self.ring._gens_set:
            cpself = self.copy()
        else:
            cpself = self
        expv, coeff = mc
        c = cpself.get(expv)
        if c is None:
            cpself[expv] = coeff
        else:
            c += coeff
            if c:
                cpself[expv] = c
            else:
                del cpself[expv]
        return cpself

    def _iadd_poly_monom(self, p2, mc):
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
        p1 = self
        if p1 in p1.ring._gens_set:
            p1 = p1.copy()
        (m, c) = mc
        get = p1.get
        zero = p1.ring.domain.zero
        monomial_mul = p1.ring.monomial_mul
        for k, v in p2.items():
            ka = monomial_mul(k, m)
            coeff = get(ka, zero) + v*c
            if coeff:
                p1[ka] = coeff
            else:
                del p1[ka]
        return p1

    def degree(f, x=None):
        """
        The leading degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        i = f.ring.index(x)

        if not f:
            return ninf
        elif i < 0:
            return 0
        else:
            return max(monom[i] for monom in f.itermonoms())

    def degrees(f):
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not f:
            return (ninf,)*f.ring.ngens
        else:
            return tuple(map(max, list(zip(*f.itermonoms()))))

    def tail_degree(f, x=None):
        """
        The tail degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        i = f.ring.index(x)

        if not f:
            return ninf
        elif i < 0:
            return 0
        else:
            return min(monom[i] for monom in f.itermonoms())

    def tail_degrees(f):
        """
        A tuple containing tail degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not f:
            return (ninf,)*f.ring.ngens
        else:
            return tuple(map(min, list(zip(*f.itermonoms()))))

    def leading_expv(self):
        """Leading monomial tuple according to the monomial ordering.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y, z = ring('x, y, z', ZZ)
        >>> p = x**4 + x**3*y + x**2*z**2 + z**7
        >>> p.leading_expv()
        (4, 0, 0)

        """
        if self:
            return self.ring.leading_expv(self)
        else:
            return None

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def coeff(self, element):
        """
        Returns the coefficient that stands next to the given monomial.

        Parameters
        ==========

        element : PolyElement (with ``is_monomial = True``) or 1

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y, z = ring("x,y,z", ZZ)
        >>> f = 3*x**2*y - x*y*z + 7*z**3 + 23

        >>> f.coeff(x**2*y)
        3
        >>> f.coeff(x*y)
        0
        >>> f.coeff(1)
        23

        """
        if element == 1:
            return self._get_coeff(self.ring.zero_monom)
        elif self.ring.is_element(element):
            terms = list(element.iterterms())
            if len(terms) == 1:
                monom, coeff = terms[0]
                if coeff == self.ring.domain.one:
                    return self._get_coeff(monom)

        raise ValueError("expected a monomial, got %s" % element)

    def const(self):
        """Returns the constant coefficient. """
        return self._get_coeff(self.ring.zero_monom)

    @property
    def LC(self):
        return self._get_coeff(self.leading_expv())

    @property
    def LM(self):
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom
        else:
            return expv

    def leading_monom(self):
        """
        Leading monomial as a polynomial element.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y = ring('x, y', ZZ)
        >>> (3*x*y + y**2).leading_monom()
        x*y

        """
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self.ring.domain.one
        return p

    @property
    def LT(self):
        expv = self.leading_expv()
        if expv is None:
            return (self.ring.zero_monom, self.ring.domain.zero)
        else:
            return (expv, self._get_coeff(expv))

    def leading_term(self):
        """Leading term as a polynomial element.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> _, x, y = ring('x, y', ZZ)
        >>> (3*x*y + y**2).leading_term()
        3*x*y

        """
        p = self.ring.zero
        expv = self.leading_expv()
        if expv is not None:
            p[expv] = self[expv]
        return p

    def coeffs(self, order=None):
        """Ordered list of polynomial coefficients.

        Parameters
        ==========

        order : :class:`~.MonomialOrder` or coercible, optional

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.orderings import lex, grlex

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.coeffs()
        [2, 1]
        >>> f.coeffs(grlex)
        [1, 2]

        """
        return [coeff for _, coeff in self.terms(order)]

    def monoms(self, order=None):
        """Ordered list of polynomial monomials.

        Parameters
        ==========

        order : :class:`~.MonomialOrder` or coercible, optional

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.orderings import lex, grlex

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.monoms()
        [(2, 3), (1, 7)]
        >>> f.monoms(grlex)
        [(1, 7), (2, 3)]

        """
        return [monom for monom, _ in self.terms(order)]

    def itercoeffs(self):
        """Iterator over coefficients of a polynomial. """
        return iter(self.values())

    def itermonoms(self):
        """Iterator over monomials of a polynomial. """
        return iter(self.keys())

    def iterterms(self):
        """Iterator over terms of a polynomial. """
        return iter(self.items())

    def listcoeffs(self):
        """Unordered list of polynomial coefficients. """
        return list(self.values())

    def listmonoms(self):
        """Unordered list of polynomial monomials. """
        return list(self.keys())

    def listterms(self):
        """Unordered list of polynomial terms. """
        return list(self.items())

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
        if p in p.ring._gens_set:
            return p*c
        if not c:
            p.clear()
            return
        for exp in p:
            p[exp] *= c
        return p

    def content(f):
        """Returns GCD of polynomial's coefficients. """
        domain = f.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in f.itercoeffs():
            cont = gcd(cont, coeff)

        return cont

    def primitive(f):
        """Returns content and a primitive polynomial. """
        cont = f.content()
        if cont == f.ring.domain.zero:
            return (cont, f)
        return cont, f.quo_ground(cont)

    def monic(f):
        """Divides all coefficients by the leading coefficient. """
        if not f:
            return f
        else:
            return f.quo_ground(f.LC)

    def mul_ground(f, x):
        if not x:
            return f.ring.zero

        terms = [ (monom, coeff*x) for monom, coeff in f.iterterms() ]
        return f.new(terms)

    def mul_monom(f, monom):
        monomial_mul = f.ring.monomial_mul
        terms = [ (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in f.items() ]
        return f.new(terms)

    def mul_term(f, term):
        monom, coeff = term

        if not f or not coeff:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.mul_ground(coeff)

        monomial_mul = f.ring.monomial_mul
        terms = [ (monomial_mul(f_monom, monom), f_coeff*coeff) for f_monom, f_coeff in f.items() ]
        return f.new(terms)

    def quo_ground(f, x):
        domain = f.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not f or x == domain.one:
            return f

        if domain.is_Field:
            quo = domain.quo
            terms = [ (monom, quo(coeff, x)) for monom, coeff in f.iterterms() ]
        else:
            terms = [ (monom, coeff // x) for monom, coeff in f.iterterms() if not (coeff % x) ]

        return f.new(terms)

    def quo_term(f, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError("polynomial division")
        elif not f:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.quo_ground(coeff)

        term_div = f._term_div()

        terms = [ term_div(t, term) for t in f.iterterms() ]
        return f.new([ t for t in terms if t is not None ])

    def trunc_ground(f, p):
        if f.ring.domain.is_ZZ:
            terms = []

            for monom, coeff in f.iterterms():
                coeff = coeff % p

                if coeff > p // 2:
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [ (monom, coeff % p) for monom, coeff in f.iterterms() ]

        poly = f.new(terms)
        poly.strip_zero()
        return poly

    rem_ground = trunc_ground

    def extract_ground(self, g):
        f = self
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def _norm(f, norm_func):
        if not f:
            return f.ring.domain.zero
        else:
            ground_abs = f.ring.domain.abs
            return norm_func([ ground_abs(coeff) for coeff in f.itercoeffs() ])

    def max_norm(f):
        return f._norm(max)

    def l1_norm(f):
        return f._norm(sum)

    def deflate(f, *G):
        ring = f.ring
        polys = [f] + list(G)

        J = [0]*ring.ngens

        for p in polys:
            for monom in p.itermonoms():
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

            for I, coeff in p.iterterms():
                N = [ i // j for i, j in zip(I, J) ]
                h[tuple(N)] = coeff

            H.append(h)

        return J, H

    def inflate(f, J):
        poly = f.ring.zero

        for I, coeff in f.iterterms():
            N = [ i*j for i, j in zip(I, J) ]
            poly[tuple(N)] = coeff

        return poly

    def lcm(self, g):
        f = self
        domain = f.ring.domain

        if not domain.is_Field:
            fc, f = f.primitive()
            gc, g = g.primitive()
            c = domain.lcm(fc, gc)

        h = (f*g).quo(f.gcd(g))

        if not domain.is_Field:
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
        monomial_gcd = ring.monomial_gcd
        monomial_ldiv = ring.monomial_ldiv
        mf, cf = list(f.iterterms())[0]
        _mgcd, _cgcd = mf, cf
        for mg, cg in g.iterterms():
            _mgcd = monomial_gcd(_mgcd, mg)
            _cgcd = ground_gcd(_cgcd, cg)
        h = f.new([(_mgcd, _cgcd)])
        cff = f.new([(monomial_ldiv(mf, _mgcd), ground_quo(cf, _cgcd))])
        cfg = f.new([(monomial_ldiv(mg, _mgcd), ground_quo(cg, _cgcd)) for mg, cg in g.iterterms()])
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

    def _gcd_QQ(self, g):
        f = self
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

    def cancel(self, g):
        """
        Cancel common factors in a rational function ``f/g``.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x,y = ring("x,y", ZZ)

        >>> (2*x**2 - 2).cancel(x**2 - 2*x + 1)
        (2*x + 2, x - 1)

        """
        f = self
        ring = f.ring

        if not f:
            return f, ring.one

        domain = ring.domain

        if not (domain.is_Field and domain.has_assoc_Ring):
            _, p, q = f.cofactors(g)
        else:
            new_ring = ring.clone(domain=domain.get_ring())

            cq, f = f.clear_denoms()
            cp, g = g.clear_denoms()

            f = f.set_ring(new_ring)
            g = g.set_ring(new_ring)

            _, p, q = f.cofactors(g)
            _, cp, cq = new_ring.domain.cofactors(cp, cq)

            p = p.set_ring(ring)
            q = q.set_ring(ring)

            p = p.mul_ground(cp)
            q = q.mul_ground(cq)

        # Make canonical with respect to sign or quadrant in the case of ZZ_I
        # or QQ_I. This ensures that the LC of the denominator is canonical by
        # multiplying top and bottom by a unit of the ring.
        u = q.canonical_unit()
        if u == domain.one:
            pass
        elif u == -domain.one:
            p, q = -p, -q
        else:
            p = p.mul_ground(u)
            q = q.mul_ground(u)

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
        i = ring.index(x)
        m = ring.monomial_basis(i)
        g = ring.zero
        for expv, coeff in f.iterterms():
            if expv[i]:
                e = ring.monomial_ldiv(expv, m)
                g[e] = ring.domain_new(coeff*expv[i])
        return g

    def evaluate(self, x, a=None):
        f = self

        if isinstance(x, list) and a is None:
            (X, a), x = x[0], x[1:]
            f = f.evaluate(X, a)

            if not x:
                return f
            else:
                x = [ (Y.drop(X), a) for (Y, a) in x ]
                return f.evaluate(x)

        ring = f.ring
        i = ring.index(x)
        a = ring.domain.convert(a)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.iterterms():
                result += coeff*a**n

            return result
        else:
            poly = ring.drop(x).zero

            for monom, coeff in f.iterterms():
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

    def subs(self, x, a=None):
        f = self

        if isinstance(x, list) and a is None:
            for X, a in x:
                f = f.subs(X, a)
            return f

        ring = f.ring
        i = ring.index(x)
        a = ring.domain.convert(a)

        if ring.ngens == 1:
            result = ring.domain.zero

            for (n,), coeff in f.iterterms():
                result += coeff*a**n

            return ring.ground_new(result)
        else:
            poly = ring.zero

            for monom, coeff in f.iterterms():
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

    def symmetrize(self):
        r"""
        Rewrite *self* in terms of elementary symmetric polynomials.

        Explanation
        ===========

        If this :py:class:`~.PolyElement` belongs to a ring of $n$ variables,
        we can try to write it as a function of the elementary symmetric
        polynomials on $n$ variables. We compute a symmetric part, and a
        remainder for any part we were not able to symmetrize.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> R, x, y = ring("x,y", ZZ)

        >>> f = x**2 + y**2
        >>> f.symmetrize()
        (x**2 - 2*y, 0, [(x, x + y), (y, x*y)])

        >>> f = x**2 - y**2
        >>> f.symmetrize()
        (x**2 - 2*y, -2*y**2, [(x, x + y), (y, x*y)])

        Returns
        =======

        Triple ``(p, r, m)``
            ``p`` is a :py:class:`~.PolyElement` that represents our attempt
            to express *self* as a function of elementary symmetric
            polynomials. Each variable in ``p`` stands for one of the
            elementary symmetric polynomials. The correspondence is given
            by ``m``.

            ``r`` is the remainder.

            ``m`` is a list of pairs, giving the mapping from variables in
            ``p`` to elementary symmetric polynomials.

            The triple satisfies the equation ``p.compose(m) + r == self``.
            If the remainder ``r`` is zero, *self* is symmetric. If it is
            nonzero, we were not able to represent *self* as symmetric.

        See Also
        ========

        sympy.polys.polyfuncs.symmetrize

        References
        ==========

        .. [1] Lauer, E. Algorithms for symmetrical polynomials, Proc. 1976
            ACM Symp. on Symbolic and Algebraic Computing, NY 242-247.
            https://dl.acm.org/doi/pdf/10.1145/800205.806342

        """
        f = self.copy()
        ring = f.ring
        n = ring.ngens

        if not n:
            return f, ring.zero, []

        polys = [ring.symmetric_poly(i+1) for i in range(n)]

        poly_powers = {}
        def get_poly_power(i, n):
            if (i, n) not in poly_powers:
                poly_powers[(i, n)] = polys[i]**n
            return poly_powers[(i, n)]

        indices = list(range(n - 1))
        weights = list(range(n, 0, -1))

        symmetric = ring.zero

        while f:
            _height, _monom, _coeff = -1, None, None

            for i, (monom, coeff) in enumerate(f.terms()):
                if all(monom[i] >= monom[i + 1] for i in indices):
                    height = max(n*m for n, m in zip(weights, monom))

                    if height > _height:
                        _height, _monom, _coeff = height, monom, coeff

            if _height != -1:
                monom, coeff = _monom, _coeff
            else:
                break

            exponents = []
            for m1, m2 in zip(monom, monom[1:] + (0,)):
                exponents.append(m1 - m2)

            symmetric += ring.term_new(tuple(exponents), coeff)

            product = coeff
            for i, n in enumerate(exponents):
                product *= get_poly_power(i, n)
            f -= product

        mapping = list(zip(ring.gens, polys))

        return symmetric, f, mapping

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
                replacements = sorted(x.items(), key=lambda k: gens_map[k[0]])
            else:
                raise ValueError("expected a generator, value pair a sequence of such pairs")

        for k, (x, g) in enumerate(replacements):
            replacements[k] = (gens_map[x], ring.ring_new(g))

        for monom, coeff in f.iterterms():
            monom = list(monom)
            subpoly = ring.one

            for i, g in replacements:
                n, monom[i] = monom[i], 0
                if n:
                    subpoly *= g**n

            subpoly = subpoly.mul_term((tuple(monom), coeff))
            poly += subpoly

        return poly

    def coeff_wrt(self, x, deg):
        """
        Coefficient of ``self`` with respect to ``x**deg``.

        Treating ``self`` as a univariate polynomial in ``x`` this finds the
        coefficient of ``x**deg`` as a polynomial in the other generators.

        Parameters
        ==========

        x : generator or generator index
            The generator or generator index to compute the expression for.
        deg : int
            The degree of the monomial to compute the expression for.

        Returns
        =======

        :py:class:`~.PolyElement`
            The coefficient of ``x**deg`` as a polynomial in the same ring.

        Examples
        ========

        >>> from sympy.polys import ring, ZZ
        >>> R, x, y, z = ring("x, y, z", ZZ)

        >>> p = 2*x**4 + 3*y**4 + 10*z**2 + 10*x*z**2
        >>> deg = 2
        >>> p.coeff_wrt(2, deg) # Using the generator index
        10*x + 10
        >>> p.coeff_wrt(z, deg) # Using the generator
        10*x + 10
        >>> p.coeff(z**2) # shows the difference between coeff and coeff_wrt
        10

        See Also
        ========

        coeff, coeffs

        """
        p = self
        i = p.ring.index(x)
        terms = [(m, c) for m, c in p.iterterms() if m[i] == deg]

        if not terms:
            return p.ring.zero

        monoms, coeffs = zip(*terms)
        monoms = [m[:i] + (0,) + m[i + 1:] for m in monoms]
        return p.ring.from_dict(dict(zip(monoms, coeffs)))

    def half_gcdex(f, g):
        return f.ring.dmp_half_gcdex(f, g)

    def gcdex(f, g):
        return f.ring.dmp_gcdex(f, g)

    def resultant(f, g):
        return f.ring.dmp_resultant(f, g)

    def discriminant(f):
        return f.ring.dmp_discriminant(f)

    def decompose(f):
        if f.ring.is_univariate:
            return f.ring.dup_decompose(f)
        else:
            raise MultivariatePolynomialError("polynomial decomposition")

    def shift(f, a):
        if f.ring.is_univariate:
            return f.ring.dup_shift(f, a)
        else:
            raise MultivariatePolynomialError("shift: use shift_list instead")

    def shift_list(f, a):
        return f.ring.dmp_shift(f, a)

    def sturm(f):
        if f.ring.is_univariate:
            return f.ring.dup_sturm(f)
        else:
            raise MultivariatePolynomialError("sturm sequence")

    def gff_list(f):
        return f.ring.dmp_gff_list(f)

    def norm(f):
        return f.ring.dmp_norm(f)

    def sqf_norm(f):
        return f.ring.dmp_sqf_norm(f)

    def sqf_part(f):
        return f.ring.dmp_sqf_part(f)

    def sqf_list(f, all=False):
        return f.ring.dmp_sqf_list(f, all=all)

    def factor_list(f):
        return f.ring.dmp_factor_list(f)

    def terms(self, order=None):
        """Ordered list of polynomial terms.

        Parameters
        ==========

        order : :class:`~.MonomialOrder` or coercible, optional

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.orderings import lex, grlex

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.terms()
        [((2, 3), 2), ((1, 7), 1)]
        >>> f.terms(grlex)
        [((1, 7), 1), ((2, 3), 2)]

        """
        return self._sorted(list(self.items()), order)


class PythonPolyRing(PolyRing):
    """Multivariate distributed polynomial ring using pure Python."""

    element_class = PythonPolyElement

    def __init__(self, symbols, domain, order=lex):
        super().__init__(symbols, domain, order)

        self.dtype = lambda init: self.element_class(self, init)

        self._one = self.dtype([(self.zero_monom, self.domain.one)])

        if self.ngens:
            codegen = MonomialOps(self.ngens)
            self.monomial_mul = codegen.mul()
            self.monomial_pow = codegen.pow()
            self.monomial_mulpow = codegen.mulpow()
            self.monomial_ldiv = codegen.ldiv()
            self.monomial_div = codegen.div()
            self.monomial_lcm = codegen.lcm()
            self.monomial_gcd = codegen.gcd()
        else:
            monunit = lambda a, b: ()
            self.monomial_mul = monunit
            self.monomial_pow = monunit
            self.monomial_mulpow = lambda a, b, c: ()
            self.monomial_ldiv = monunit
            self.monomial_div = monunit
            self.monomial_lcm = monunit
            self.monomial_gcd = monunit

        if self.order is lex:
            self.leading_expv = max
        else:
            self.leading_expv = lambda f: max(f, key=self.order)

        self.gens = self._gens()
        self._gens_set = set(self.gens)

        for symbol, generator in zip(self.symbols, self.gens):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if not hasattr(self, name):
                    setattr(self, name, generator)

    __hash__ = PolyRing.__hash__

    def __eq__(self, other):
        return isinstance(other, PolyRing) and \
            (self.symbols, self.domain, self.ngens, self.order) == \
            (other.symbols, other.domain, other.ngens, other.order)

    def _gens(self):
        """Return a list of polynomial generators."""
        one = self.domain.one
        _gens = []
        for i in range(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            _gens.append(poly)
        return tuple(_gens)

    def __getstate__(self):
        state = self.__dict__.copy()
        if "leading_expv" in state:
            del state["leading_expv"]
        for key in list(state.keys()):
            if key.startswith("monomial_"):
                del state[key]
        return state

    @property
    def zero(self):
        return self.dtype([])

    @property
    def one(self):
        return self.dtype(self._one)

    def is_element(self, element):
        return isinstance(element, PolyElement) and element.ring == self

    def ring_new(self, element):
        if isinstance(element, (PolyElement, FlintPolyElement)):
            if self == element.ring:
                return element
            elif isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, str):
            raise NotImplementedError("parsing")
        elif isinstance(element, dict):
            return self.from_dict(element)
        elif isinstance(element, list):
            try:
                return self.from_terms(element)
            except ValueError:
                return self.from_list(element)
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            return self.ground_new(element)

    __call__ = ring_new

    def from_dict(self, element, orig_domain=None):
        domain_new = self.domain_new
        poly = self.zero
        for monom, coeff in element.items():
            coeff = domain_new(coeff, orig_domain)
            if coeff:
                poly[monom] = coeff
        return poly

    def _rebuild_expr(self, expr, mapping):
        domain = self.domain

        def _rebuild(expr):
            generator = mapping.get(expr)
            if generator is not None:
                return generator
            elif expr.is_Add:
                return reduce(add, list(map(_rebuild, expr.args)))
            elif expr.is_Mul:
                return reduce(mul, list(map(_rebuild, expr.args)))
            else:
                base, exp = expr.as_base_exp()
                if exp.is_Integer and exp > 1:
                    return _rebuild(base)**int(exp)
                else:
                    return self.ground_new(domain.convert(expr))

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        mapping = dict(list(zip(self.symbols, self.gens)))
        try:
            poly = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a polynomial in %s, got %s" % (self, expr))
        else:
            return self.ring_new(poly)

    def index(self, gen):
        if gen is None:
            if self.ngens:
                i = 0
            else:
                i = -1
        elif isinstance(gen, int):
            i = gen
            if 0 <= i and i < self.ngens:
                pass
            elif -self.ngens <= i and i <= -1:
                i = -i - 1
            else:
                raise ValueError("invalid generator index: %s" % gen)
        elif self.is_element(gen):
            try:
                i = self.gens.index(gen)
            except ValueError:
                raise ValueError("invalid generator: %s" % gen)
        elif isinstance(gen, str):
            try:
                i = self.symbols.index(gen)
            except ValueError:
                raise ValueError("invalid generator: %s" % gen)
        else:
            raise ValueError("expected a polynomial generator, an integer, a string or None, got %s" % gen)
        return i

    def drop(self, *gens):
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]
        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    @property
    def is_univariate(self):
        return len(self.gens) == 1

    @property
    def is_multivariate(self):
        return len(self.gens) > 1

    def compose(self, other):
        if self != other:
            syms = set(self.symbols).union(set(other.symbols))
            return self.clone(symbols=list(syms))
        else:
            return self

    def add_gens(self, symbols):
        syms = set(self.symbols).union(set(symbols))
        return self.clone(symbols=list(syms))

    def symmetric_poly(self, n):
        if n < 0 or n > self.ngens:
            raise ValueError("Cannot generate symmetric polynomial of order %s for %s" % (n, self.gens))
        elif not n:
            return self.one
        else:
            poly = self.zero
            for s in subsets(range(self.ngens), int(n)):
                monom = tuple(int(i in s) for i in range(self.ngens))
                poly += self.term_new(monom, self.domain.one)
            return poly

    def term_new(self, monom, coeff):
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    def add(self, *objs):
        """
        Add a sequence of polynomials or containers of polynomials.
        """
        p = self.zero
        for obj in objs:
            if is_sequence(obj, include=GeneratorType):
                p += self.add(*obj)
            else:
                p += obj
        return p

    def mul(self, *objs):
        """
        Multiply a sequence of polynomials or containers of polynomials.
        """
        p = self.one
        for obj in objs:
            if is_sequence(obj, include=GeneratorType):
                p *= self.mul(*obj)
            else:
                p *= obj
        return p

    def drop_to_ground(self, *gens):
        """
        Remove specified generators from the ring and inject them into its domain.
        """
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]
        gens = [gen for i, gen in enumerate(self.gens) if i not in indices]
        if not symbols:
            return self
        else:
            return self.clone(symbols=symbols, domain=self.drop(*gens))


class FlintPolyElement(PolyElement):

    def __init__(self, ring, flint_poly):
        super().__init__(ring)
        self.flint_poly = flint_poly

    def _check(self):
        assert isinstance(self, PolyElement)
        assert isinstance(self.ring, PolyRing)
        dom = self.ring.domain
        assert isinstance(dom, Domain)
        for monom, coeff in self.flint_poly.to_dict().items():
            assert dom.of_type(coeff)
            assert len(monom) == self.ring.ngens
            assert all(isinstance(exp, int) and exp >= 0 for exp in monom)

    def new(self, flint_poly):
        return self.__class__(self.ring, flint_poly)

    def __hash__(self):
        poly_dict = self.flint_poly.to_dict()
        poly_tuple = tuple(sorted(poly_dict.items()))
        return hash((self.ring, poly_tuple))

    def copy(self):
        new_poly = self.ring.flint_ctx.from_dict(self.flint_poly.to_dict())
        return self.new(new_poly)

    def set_ring(self, new_ring):
        if self.ring == new_ring:
            return self

        terms_dict = self.flint_poly.to_dict()

        if self.ring.symbols != new_ring.symbols:
            new_keys, new_vals = _dict_reorder(terms_dict, self.ring.symbols, new_ring.symbols)
            terms = list(zip(new_keys, new_vals))
            return new_ring.from_terms(terms, self.ring.domain)
        else:
            return new_ring.from_dict(terms_dict, self.ring.domain)

    def clear_denoms(self):
        domain = self.ring.domain

        if not domain.is_Field or not domain.has_assoc_Ring:
            return domain.one, self

        ground_ring = domain.get_ring()
        common = ground_ring.one
        lcm = ground_ring.lcm
        denom = domain.denom

        for coeff in self.flint_poly.to_dict().values():
            common = lcm(common, denom(coeff))

        scaled_poly = self.flint_poly * common
        result_poly = self.ring.ring_new(scaled_poly)

        return common, result_poly

    def strip_zero(self):
        raise NotImplementedError

    def __eq__(self, other):
        if isinstance(other, FlintPolyElement):
            return self.flint_poly == other.flint_poly
        elif not other:
            return self.flint_poly.is_zero()
        elif self.ring.is_element(other):
            return self == self.ring(other)
        elif self.flint_poly.is_zero():
            return False
        elif len(self.flint_poly.to_dict()) > 1:
            return False
        else:
            constant_term = self.flint_poly.to_dict().get(tuple([0] * self.ring.ngens), 0)
            return constant_term == other

    def almosteq(self, other, tolerance=None):
        raise NotImplementedError

    def sort_key(self):
        terms_dict = self.flint_poly.to_dict()
        return len(terms_dict), self.terms()

    def _drop(self, gen):
        ring = self.ring
        i = ring.index(gen)

        if ring.ngens == 1:
            return i, ring.domain
        else:
            symbols = list(ring.symbols)
            del symbols[i]
            return i, ring.clone(symbols=symbols)

    def drop(self, gen):
        i, ring = self._drop(gen)

        if self.ring.ngens == 1:
            terms = self.flint_poly.to_dict()
            if len(terms) <= 1 and (0,) in terms:
                return self.ring.domain(terms.get((0,), 0))
            else:
                raise ValueError("Cannot drop %s" % gen)
        else:
            poly = ring.zero
            terms = self.flint_poly.to_dict()

            for k, v in terms.items():
                if k[i] == 0:
                    K = [int(exp) for exp in k]
                    del K[i]
                    poly[tuple(K)] = v
                else:
                    raise ValueError("Cannot drop %s" % gen)

            return poly

    def drop_to_ground(self, gen):
        raise NotImplementedError

    def to_dense(self):
        raise NotImplementedError

    @property
    def is_generator(self):
        return self in self.ring._gens_set

    @property
    def is_ground(self):
        return self.flint_poly.is_constant()

    @property
    def is_monomial(self):
        return self.flint_poly.is_zero() or (len(self.flint_poly.to_dict()) == 1 and self.LC == 1)

    @property
    def is_term(self):
        return len(self.flint_poly.to_dict()) <= 1

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

    @property
    def is_zero(self):
        return self.flint_poly.is_zero()

    @property
    def is_one(self):
        return self == self.ring.one

    @property
    def is_monic(self):
        return self.ring.domain.is_one(self.LC)

    @property
    def is_primitive(self):
        content, _ = self.flint_poly.factor()
        return self.ring.domain.is_one(content)

    @property
    def is_linear(self):
        return all(sum(monom) <= 1 for monom in self.itermonoms())

    @property
    def is_quadratic(self):
        return all(sum(monom) <= 2 for monom in self.itermonoms())

    @property
    def is_squarefree(self):
        if not self.ring.ngens:
            return True
        _, factors = self.flint_poly.factor_squarefree()
        return all(exp == 1 for _, exp in factors)

    @property
    def is_irreducible(self):
        raise NotImplementedError

    @property
    def is_cyclotomic(self):
        if not self.ring.is_univariate:
            raise MultivariatePolynomialError("cyclotomic polynomial")
        else:
            raise NotImplementedError

    def __add__(self, other):
        if self.ring.is_element(other):
            return self.new(self.flint_poly + other.flint_poly)
        elif other in self.ring.domain:
            if self.ring.domain.is_ZZ:
                return self.new(self.flint_poly + flint.fmpz(other))
            elif self.ring.domain.is_QQ:
                return self.new(self.flint_poly + flint.fmpq(other))
        return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if self.ring.is_element(other):
            return self.new(self.flint_poly - other.flint_poly)
        elif other in self.ring.domain:
            if self.ring.domain.is_ZZ:
                return self.new(self.flint_poly - flint.fmpz(other))
            elif self.ring.domain.is_QQ:
                return self.new(self.flint_poly - flint.fmpq(other))
        return NotImplemented

    def __rsub__(self, other):
        if other in self.ring.domain:
            if self.ring.domain.is_ZZ:
                const = flint.fmpz(other)
            elif self.ring.domain.is_QQ:
                const = flint.fmpq(other)
            else:
                raise DomainError(f"Unsupported domain for FlintPolyElement: {self.ring.domain}")

            result = const - self.flint_poly
            return self.new(result)
        else:
            return NotImplemented

    def __mul__(self, other):
        if self.ring.is_element(other):
            return self.new(self.flint_poly * other.flint_poly)
        elif other in self.ring.domain:
            if self.ring.domain.is_ZZ:
                return self.new(self.flint_poly * flint.fmpz(other))
            elif self.ring.domain.is_QQ:
                return self.new(self.flint_poly * flint.fmpq(other))
        return NotImplemented

    def __rmul__(self, other):

        if other in self.ring.domain:
            if self.ring.domain.is_ZZ:
                return self.new(flint.fmpz(other) * self.flint_poly)
            elif self.ring.domain.is_QQ:
                return self.new(flint.fmpq(other) * self.flint_poly)
        return NotImplemented

    def __pow__(self, n):
        if not isinstance(n, int):
            raise TypeError("exponent must be an integer, got %s" % n)
        if n < 0:
            raise ValueError("exponent must be a non-negative integer, got %s" % n)
        if n == 0:
            if self.flint_poly != 0:
                return self.ring.one
            else:
                raise ValueError("0**0")
        return self.new(self.flint_poly ** n)

    def square(self):
        return NotImplementedError

    def __divmod__(self, other):
        raise NotImplementedError

    def __rdivmod__(self, other):
        raise NotImplementedError

    def __floordiv__(self, other):
        """Compute the quotient of polynomial division."""
        if isinstance(other, FlintPolyElement) and other.ring == self.ring:
            if other.flint_poly.is_zero():
                raise ZeroDivisionError("polynomial division by zero")
            result = self.flint_poly // other.flint_poly
            new_poly = object.__new__(FlintPolyElement)
            new_poly.ring = self.ring
            new_poly.flint_poly = result
            return new_poly
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, FlintPolyElement) and other.ring == self.ring:

            if other.flint_poly.is_zero():
                raise ZeroDivisionError("polynomial division by zero")
            result = self.flint_poly / other.flint_poly
            new_poly = object.__new__(FlintPolyElement)
            new_poly.ring = self.ring
            new_poly.flint_poly = result
            return new_poly
        elif isinstance(other, (int, float, flint.fmpq)):

            if other == 0:
                raise ZeroDivisionError("polynomial division by zero")
            if self.ring.domain.is_ZZ:
                raise DomainError("Division by scalar in ZZ ring requires QQ domain")
            elif self.ring.domain.is_QQ:
                if not isinstance(other, flint.fmpq):
                    other = flint.fmpq(other)
                result = self.flint_poly / other
                new_poly = object.__new__(FlintPolyElement)
                new_poly.ring = self.ring
                new_poly.flint_poly = result
                return new_poly
        return NotImplemented

    def __rtruediv__(self, other):
        return NotImplemented

    def div(self, fv):
        raise NotImplementedError

    def rem(self, G):
        raise NotImplementedError

    def _iadd_monom(self, mc):
        raise NotImplementedError

    def _iadd_poly_monom(self, p2, mc):
        raise NotImplementedError

    def degree(self, x=None):
        if self.flint_poly.is_zero():
            return ninf

        if self.ring.ngens == 0:
            return 0

        degs = self.flint_poly.degrees()
        if x is None:
            return degs[0]
        else:
            i = self.ring.index(x)
            return degs[i]

    def degrees(self):
        if self.flint_poly == 0:
            return tuple([-float('inf')] * self.ring.ngens)

        flint_degrees = self.flint_poly.degrees()
        return flint_degrees

    def tail_degree(self, x=None):

        ring = self.ring
        if x is None:
            i = 0
        else:
            i = ring.index(x)

        if self.flint_poly.is_zero():
            return ninf

        if ring.ngens == 0:
            return 0

        if i < 0:
            return 0

        return min(monom[i] for monom in self.itermonoms())

    def tail_degrees(self):
        if self.flint_poly.is_zero():
            return (ninf,) * self.ring.ngens
        else:
            return tuple(map(min, list(zip(*self.itermonoms()))))

    def coeff(self, element):
        if element == 1:
            expv = (0,) * self.ring.ngens
        elif self.ring.is_element(element):
            terms = list(element.flint_poly.terms())
            if len(terms) != 1:
                raise ValueError("expected a monomial, got %s" % element)
            monom, coeff = terms[0]
            if not self.ring.domain.is_one(coeff):
                raise ValueError("expected a monomial with coefficient 1, got %s" % element)
            expv = monom
        else:
            raise ValueError("expected a monomial, got %s" % element)

        for monom, coeff in self.flint_poly.terms():
            if monom == expv:
                return self.ring.domain(coeff)
        return self.ring.domain.zero

    def const(self):
        return self.flint_poly.to_dict().get(tuple([0] * self.ring.ngens), 0)

    @property
    def LC(self):
        return self.flint_poly.leading_coefficient()

    @property
    def LM(self):
        terms_list = list(self.flint_poly.terms())
        if not terms_list:
            return self.ring.zero_monom
        return terms_list[0][0]

    def leading_monom(self):
        raise NotImplementedError

    @property
    def LT(self):
        terms_list = list(self.flint_poly.terms())
        if not terms_list:
            return (self.ring.zero_monom, self.ring.domain.zero)
        expv, coeff = terms_list[0]
        return (expv, self.ring.domain(coeff))

    def leading_term(self):
        raise NotImplementedError

    def coeffs(self, order=None):
        if order is None:
            return self.flint_poly.coeffs()
        else:
            return [coeff for _, coeff in self.terms(order)]

    def monoms(self, order=None):
        if order is None:
            return self.flint_poly.monoms()
        else:
            return [monom for monom, _ in self.terms(order)]

    def terms(self, order=None):
        return self._sorted(list(self.flint_poly.to_dict().items()), order)

    def itercoeffs(self):
        return iter(self.flint_poly.coeffs())

    def itermonoms(self):
        return iter(self.flint_poly.monoms())

    def iterterms(self):
        return iter(self.flint_poly.to_dict().items())

    def listcoeffs(self):
        """Unordered list of polynomial coefficients. """
        return self.flint_poly.coeffs()

    def listmonoms(self):
        """Unordered list of polynomial monomials. """
        return self.flint_poly.monoms()

    def listterms(self):
        raise NotImplementedError

    def imul_num(self, other):
        raise NotImplementedError

    def content(self):
        """Return the GCD of the coefficients of self."""
        return self.flint_poly.content()

    def primitive(self):
        cont, prim = self.flint_poly.primitive()
        return cont, FlintPolyElement(self.ring, prim)

    def monic(self):
        raise NotImplementedError

    def mul_ground(self, other):
        raise NotImplementedError

    def mul_monom(self, monom):
        raise NotImplementedError

    def mul_term(self, term):
        raise NotImplementedError

    def quo_ground(self, x):
        raise NotImplementedError

    def quo_term(self, term):
        raise NotImplementedError

    def trunc_ground(self, other):
        raise NotImplementedError

    def extract_ground(self, g):
        raise NotImplementedError

    def max_norm(self):
        raise NotImplementedError

    def li_norm(self):
        raise NotImplementedError

    def deflate(f, *G):
        ring = f.ring
        polys = [f] + list(G)

        J = [0] * ring.ngens
        for p in polys:
            for monom in p.itermonoms():
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
            deflated_flint_poly = p.flint_poly.deflate(list(J))
            h = ring.ring_new(deflated_flint_poly)
            H.append(h)

        return J, H

    def inflate(f, J):
        J_list = list(J)
        inflated_flint_poly = f.flint_poly.inflate(J_list)
        return f.ring.ring_new(inflated_flint_poly)

    def lcm(self, other):
        raise NotImplementedError

    def gcd(self, g):
        """Return the GCD of self and g."""
        return self.new(self.flint_poly.gcd(g.flint_poly))

    def cofactors(self, other):
        raise NotImplementedError

    def cancel(self, g):
        c = self.gcd(g)
        if c.flint_poly.is_one():
            return self, g
        f_new = self.new(self.flint_poly // c.flint_poly)
        g_new = g.new(g.flint_poly // c.flint_poly)
        return f_new, g_new

    def diff(self, x):
        ring = self.ring
        i = ring.index(x)
        derivative_poly = self.flint_poly.derivative(i)
        return ring.ring_new(derivative_poly)

    def evaluate(selfself, x, a=None):
        raise NotImplementedError

    def subs(self, x, a=None):
        """Substitute a generator with a value in the polynomial."""
        if isinstance(x, list) and a is None:
            result = self
            for X, val in x:
                result = result.subs(X, val)
            return result

        ring = self.ring
        if ring.ngens == 0:
            raise ValueError("Cannot substitute in a ring with no generators")

        i = ring.index(x)
        var_name = ring.symbol_strings[i]

        a = ring.domain.convert(a)

        if ring.domain.is_ZZ:
            flint_val = flint.fmpz(int(a))
        elif ring.domain.is_QQ:
            flint_val = flint.fmpq(a.numerator, a.denominator)
        else:
            raise DomainError(f"Unsupported domain for substitution: {ring.domain}")

        subst_dict = {var_name: flint_val}
        result_poly = self.flint_poly.subs(subst_dict)

        return FlintPolyElement(ring, result_poly)

    def symmetrize(self):
        raise NotImplementedError

    def compose(self, x, a=None):
        """Compose this polynomial with other polynomials."""
        ring = self.ring
        gens_map = dict(zip(ring.gens, range(ring.ngens)))

        if a is not None:
            replacements = [(x, a)]
        else:
            if isinstance(x, list):
                replacements = list(x)
            elif isinstance(x, dict):
                replacements = sorted(x.items(), key=lambda k: gens_map[k[0]])
            else:
                raise ValueError("expected a generator, value pair or a sequence of such pairs")

        if not replacements:
            return self

        subst_polys = list(ring.gens)
        target_ring = ring

        for gen, poly in replacements:
            if gen not in gens_map:
                raise ValueError(f"Generator {gen} not in ring {ring}")

            idx = gens_map[gen]

            if not isinstance(poly, FlintPolyElement):
                poly = ring(poly)
            else:
                if poly.ring.ngens > 0:
                    target_ring = poly.ring
                    subst_polys = [target_ring(p) if isinstance(p, FlintPolyElement) else target_ring.gens[i]
                                   for i, p in enumerate(subst_polys)]

            subst_polys[idx] = poly

        subst_polys = [target_ring(p) if p.ring != target_ring else p for p in subst_polys]
        result_poly = self.flint_poly.compose(*[p.flint_poly for p in subst_polys], ctx=target_ring.flint_ctx)
        return FlintPolyElement(target_ring, result_poly)

    def coeff_wrt(self, x, deg):
        i = self.ring.index(x)
        terms = [(monom[:i] + (0,) + monom[i + 1:], coeff)
                 for monom, coeff in self.flint_poly.terms() if monom[i] == deg]
        if not terms:
            return self.ring.zero
        return self.new(self.ring.flint_ctx.from_dict(dict(terms)))

    def half_gcdex(self, g):
        raise NotImplementedError

    def gcdex(self, g):
        raise NotImplementedError

    def resultant(self, g):
        raise NotImplementedError

    def discriminant(self):
        """Compute the discriminant of the polynomial with respect to the main variable."""
        if self.ring.ngens == 0:
            raise ValueError("Cannot compute discriminant in a ring with no generators")

        var = self.ring.symbol_strings[0]
        return self.flint_poly.discriminant(var)

    def decompose(self):
        raise NotImplementedError

    def shift(self, a):
        raise NotImplementedError

    def shift_list(self, a):
        raise NotImplementedError

    def sturm(self):
        raise NotImplementedError

    def gff_list(self):
        raise NotImplementedError

    def norm(self):
        raise NotImplementedError

    def sqf_norm(self):
        raise NotImplementedError

    def sqf_list(self, all=False):
        raise NotImplementedError

    def factor_list(self):
        c, factors = self.flint_poly.factor()
        factor_list = [(FlintPolyElement(self.ring, poly), exp) for poly, exp in factors]
        return c, factor_list

    def _validate_monom(self, monom):
        if not isinstance(monom, tuple) or len(monom) != self.ring.ngens:
            raise ValueError(f"Monomial exponent vector must be a tuple of length {self.ring.ngens}")
        if not all(isinstance(exp, int) and exp >= 0 for exp in monom):
            raise ValueError("Monomial exponents must be non-negative integers")

    def __getitem__(self, monom):
        self._validate_monom(monom)
        terms = self.flint_poly.to_dict()
        return terms.get(monom, self.ring.domain.zero)

    def __setitem__(self, monom, coeff):
        self._validate_monom(monom)
        try:
            coeff_converted = self.ring.domain.convert(coeff)
        except (TypeError, ValueError, CoercionFailed):
            raise ValueError(f"Cannot convert coefficient {coeff} to {self.ring.domain}")
        terms = self.flint_poly.to_dict()
        if coeff_converted == self.ring.domain.zero:
            terms.pop(monom, None)
        else:
            terms[monom] = coeff_converted
        self.flint_poly = self.ring.flint_ctx.from_dict(terms)


class FlintPolyRing(PolyRing):
    """Multivariate distributed polynomial ring using FLINT via python-flint."""
    element_class = FlintPolyElement

    def __init__(self, symbols, domain, order=lex):
        super().__init__(symbols, domain, order)

        assert isinstance(self.order, MonomialOrder), f"Order not preprocessed: {self.order}"

        self.symbol_strings = tuple(str(s) for s in self.symbols)

        if order == lex:
            flint_order = 'lex'
        elif order == grlex:
            flint_order = 'deglex'
        elif order == grevlex:
            flint_order = 'degrevlex'
        else:
            raise ValueError(f"Unsupported monomial order for FlintPolyRing: {order}")

        if self.domain.is_ZZ:
            self.flint_ctx = flint.fmpz_mpoly_ctx.get(self.symbol_strings, flint_order)
            self._convert_coeff = flint.fmpz
        elif self.domain.is_QQ:
            self.flint_ctx = flint.fmpq_mpoly_ctx.get(self.symbol_strings, flint_order)
            self._convert_coeff = flint.fmpq
        else:
            raise DomainError(f"Unsupported domain for FlintPolyRing: {domain}")

        def create_element(init):
            if isinstance(init, list):
                terms = dict(init)
            elif isinstance(init, dict):
                terms = init
            else:
                raise ValueError(f"Unsupported init type for FlintPolyRing: {type(init)}")

            converted_terms = {monom: self._convert_coeff(coeff) for monom, coeff in terms.items()}

            flint_poly = self.flint_ctx.from_dict(converted_terms)
            return self.element_class(self, flint_poly)

        self.dtype = create_element

        self._one = self.dtype({self.zero_monom: self.domain.one})

        self.gens = self._gens()
        self._gens_set = set(self.gens)

        for symbol, generator in zip(self.symbols, self.gens):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if not hasattr(self, name):
                    setattr(self, name, generator)

    __hash__ = PolyRing.__hash__

    def _gens(self):
        """Return a list of polynomial generators."""
        gens = []
        for i in range(self.ngens):
            monom = [0] * self.ngens
            monom[i] = 1
            monom = tuple(monom)
            gen = self.dtype({monom: self.domain.one})
            gens.append(gen)
        return tuple(gens)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["flint_ctx"]
        del state["_gens_set"]
        for symbol in self.symbols:
            symbol_str = str(symbol)
            if symbol_str in state:
                del state[symbol_str]
        return state

    def __eq__(self, other):
        return isinstance(other, PolyRing) and \
            (self.symbols, self.domain, self.ngens, self.order) == \
            (other.symbols, other.domain, other.ngens, other.order)

    @property
    def zero(self):
        return self.dtype([])

    @property
    def one(self):
        return self._one

    def is_element(self, element):
        return isinstance(element, FlintPolyElement) and element.ring == self

    def ring_new(self, element):
        if isinstance(element, FlintPolyElement):
            if self == element.ring:
                return element
            elif isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                raise NotImplementedError("conversion from FlintPolyElement in a different ring")
        elif isinstance(element, (flint.fmpz_mpoly, flint.fmpq_mpoly)):
            return FlintPolyElement(self, element)
        elif isinstance(element, PolyElement):
            if isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                terms = element.terms()
                element_dict = dict(terms)
                return self.from_dict(element_dict)
        elif isinstance(element, str):
            raise NotImplementedError("parsing strings into polynomials")
        elif isinstance(element, dict):
            return self.from_dict(element)
        elif isinstance(element, list):
            try:
                return self.from_terms(element)
            except ValueError:
                return self.from_list(element)
        elif isinstance(element, Expr):
            raise NotImplementedError
        else:
            return self.ground_new(element)

    __call__ = ring_new

    def from_dict(self, element, orig_domain=None):
        domain_new = self.domain_new
        poly_dict = {}

        for monom, coeff in element.items():
            coeff = domain_new(coeff, orig_domain)
            if coeff:
                poly_dict[monom] = coeff

        if not poly_dict:
            return self.zero

        return FlintPolyElement(self, self.flint_ctx.from_dict(poly_dict))

    def from_expr(self, expr):
        raise NotImplementedError("from_expr is not yet implemented for FlintPolyRing")

    def index(self, gen):
        if gen is None:
            if self.ngens:
                return 0
            else:
                return -1
        elif isinstance(gen, int):
            i = gen
            if 0 <= i < self.ngens:
                pass
            elif -self.ngens <= i <= -1:
                i = self.ngens + i
            else:
                raise ValueError(f"invalid generator index: {gen}")
        elif self.is_element(gen):
            try:
                i = self.gens.index(gen)
            except ValueError:
                raise ValueError(f"invalid generator: {gen}")
        elif isinstance(gen, str):
            try:
                i = self.symbol_strings.index(gen)
            except ValueError:
                raise ValueError(f"invalid generator: {gen}")
        else:
            raise ValueError(f"expected a polynomial generator, an integer, a string or None, got {gen}")
        return i

    def drop(self, *gens):
        indices = {self.index(g) for g in gens}
        new_symbols = tuple(s for i, s in enumerate(self.symbols) if i not in indices)
        if not new_symbols:
            return self.domain
        return self.clone(symbols=new_symbols)

    @property
    def is_univariate(self):
        return self.ngens == 1

    @property
    def is_multivariate(self):
        return self.ngens > 1

    def compose(self, other):
        if not isinstance(other, PolyRing):
            raise ValueError("Can only compose with a PolyRing subclass")
        new_symbols = tuple(set(self.symbols).union(other.symbols))
        return self.clone(symbols=list(new_symbols))

    def add_gens(self, symbols):
        parsed_symbols = _parse_symbols_flint(symbols)
        parsed_sympy_symbols = sympy.symbols(",".join(parsed_symbols)) if parsed_symbols else ()
        if isinstance(parsed_sympy_symbols, Symbol):
            parsed_sympy_symbols = (parsed_sympy_symbols,)
        new_symbols = tuple(set(self.symbols).union(parsed_sympy_symbols))
        return self.clone(symbols=new_symbols)

    def symmetric_poly(self, n):
        from itertools import combinations
        if n < 0 or n > self.ngens:
            raise ValueError(f"Invalid degree {n} for {self.ngens} generators")
        result = self.zero
        for comb in combinations(self.gens, n):
            term = self.one
            for gen in comb:
                term = term * gen
            result = result + term
        return result

    def term_new(self, monom, coeff):
        coeff = self.domain_new(coeff)
        if not coeff:
            return self.zero
        poly_dict = {monom: coeff}
        return FlintPolyElement(self, self.flint_ctx.from_dict(poly_dict))

    def ring_new(self, element):
        if isinstance(element, FlintPolyElement):
            if self == element.ring:
                return element
            elif isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                raise NotImplementedError("conversion from FlintPolyElement in a different ring")
        elif isinstance(element, (flint.fmpz_mpoly, flint.fmpq_mpoly)):
            return FlintPolyElement(self, element)
        elif isinstance(element, PolyElement):
            if isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                terms = element.terms()
                element_dict = dict(terms)
                return self.from_dict(element_dict)
        elif isinstance(element, str):
            raise NotImplementedError("parsing strings into polynomials")
        elif isinstance(element, dict):
            return self.from_dict(element)
        elif isinstance(element, list):
            try:
                return self.from_terms(element)
            except ValueError:
                return self.from_list(element)
        elif isinstance(element, Expr):
            return self.from_expr(element)
        else:
            return self.ground_new(element)

    def add(self, *objs):
        result = self.zero
        for obj in objs:
            if isinstance(obj, (list, tuple)):
                result = result + self.add(*obj)
            elif self.is_element(obj):
                result = result + obj
            elif obj in self.domain:
                result = result + self.from_domain(obj)
        return result

    def mul(self, *objs):
        result = self.one
        for obj in objs:
            if isinstance(obj, (list, tuple)):
                result = result * self.mul(*obj)
            elif self.is_element(obj):
                result = result * obj
            elif obj in self.domain:
                result = result * self.from_domain(obj)
        return result

    def drop_to_ground(self, *gens):
        if not gens:
            return self
        indices = {self.index(g) for g in gens}
        new_symbols = tuple(s for i, s in enumerate(self.symbols) if i not in indices)

        return self.clone(symbols=new_symbols)

    def from_domain(self, elem):
        if self.domain.is_ZZ:
            return FlintPolyElement(self, self.flint_ctx.from_dict({self.zero_monom: flint.fmpz(elem)}))
        elif self.domain.is_QQ:
            return FlintPolyElement(self, self.flint_ctx.from_dict({self.zero_monom: flint.fmpq(elem)}))
        raise DomainError(f"Unsupported domain for from_domain: {self.domain}")
