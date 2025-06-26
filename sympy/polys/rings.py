"""Sparse polynomial rings."""

from __future__ import annotations

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
from sympy.polys.domains.domain import Domain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.polynomialring import PolynomialRing
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.monomials import MonomialOps
from sympy.polys.orderings import lex, MonomialOrder
from sympy.polys.polyerrors import (
    CoercionFailed,
    GeneratorsError,
    ExactQuotientFailed,
    MultivariatePolynomialError,
)
from sympy.polys.polyoptions import (
    Domain as DomainOpt,
    Order as OrderOpt,
    build_options,
)
from sympy.polys.polyutils import (
    expr_from_dict,
    _dict_reorder,
    _parallel_dict_from_expr,
)
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public, subsets
from sympy.utilities.iterables import is_sequence
from sympy.utilities.magic import pollute


@public
def ring(symbols, domain, order: MonomialOrder | str = lex):
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
    pollute([sym.name for sym in _ring.symbols], _ring.gens)
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
        coeffs = sum([list(rep.values()) for rep in reps], [])

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

    raise GeneratorsError(
        "expected a string, Symbol or expression or a non-empty sequence of strings, Symbols or expressions"
    )


class PolyRing(DefaultPrinting, IPolys):
    """Multivariate distributed polynomial ring."""

    gens: tuple[PolyElement, ...]
    symbols: tuple[Expr, ...]
    ngens: int
    domain: Domain
    order: MonomialOrder

    def __new__(cls, symbols, domain, order=lex):
        # Create a new ring instance.
        symbols = tuple(_parse_symbols(symbols))
        ngens = len(symbols)
        domain = DomainOpt.preprocess(domain)
        order = OrderOpt.preprocess(order)

        # Validate that symbols do not overlap with domain symbols
        if domain.is_Composite and set(symbols) & set(domain.symbols):
            raise GeneratorsError(
                "polynomial ring and its ground domain share generators"
            )

        # Create and initialize instance
        obj = object.__new__(cls)
        obj._hash_tuple = (cls.__name__, symbols, ngens, domain, order)
        obj._hash = hash(obj._hash_tuple)
        obj.symbols = symbols
        obj.ngens = ngens
        obj.domain = domain
        obj.order = order

        # Set up polynomial creation and basic elements
        obj.dtype = PolyElement(obj, ()).new
        obj.zero_monom = (0,) * ngens
        obj.gens = obj._gens()
        obj._gens_set = set(obj.gens)
        obj._one = [(obj.zero_monom, domain.one)]

        # Initialize monomial operations
        obj._init_monomial_operations()

        # Set up leading exponent function
        obj._init_leading_expv_function(order)

        # Add generator attributes for Symbol names
        obj._add_generator_attributes()

        return obj

    def _init_monomial_operations(self):
        # Initialize monomial operations based on number of generators.
        if self.ngens:
            # Operations for rings with at least one variable
            codegen = MonomialOps(self.ngens)
            self.monomial_mul = codegen.mul()
            self.monomial_pow = codegen.pow()
            self.monomial_mulpow = codegen.mulpow()
            self.monomial_ldiv = codegen.ldiv()
            self.monomial_div = codegen.div()
            self.monomial_lcm = codegen.lcm()
            self.monomial_gcd = codegen.gcd()
        else:
            # No variables, all operations return empty tuple
            monunit = lambda a, b: ()
            self.monomial_mul = monunit
            self.monomial_pow = monunit
            self.monomial_mulpow = lambda a, b, c: ()
            self.monomial_ldiv = monunit
            self.monomial_div = monunit
            self.monomial_lcm = monunit
            self.monomial_gcd = monunit

    def _init_leading_expv_function(self, order):
        # Initialize the leading exponent vector function.
        if order is lex:
            self.leading_expv = max
        else:
            self.leading_expv = lambda f: max(f, key=order)

    def _add_generator_attributes(self):
        """Add generator attributes for Symbol names."""
        for symbol, generator in zip(self.symbols, self.gens):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if not hasattr(self, name):
                    setattr(self, name, generator)

    # Pickle support
    def __getnewargs__(self):
        return self.symbols, self.domain, self.order

    # Hash and equality
    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return isinstance(other, PolyRing) and self._ring_equality(other)

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, key):
        # Get a subring with subset of symbols.
        symbols = self.symbols[key]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    # Properties
    @property
    def zero(self):
        """The zero polynomial."""
        return self.dtype([])

    @property
    def one(self):
        """The unit polynomial."""
        return self.dtype(self._one)

    @property
    def is_univariate(self):
        """True if this is a univariate ring."""
        return self.ngens == 1

    @property
    def is_multivariate(self):
        """True if this is a multivariate ring."""
        return self.ngens > 1

    # Ring operations and cloning
    def clone(self, symbols=None, domain=None, order=None):
        """Create a clone with modified parameters."""
        # Convert list to tuple for hashability
        if symbols is not None and isinstance(symbols, list):
            symbols = tuple(symbols)
        return self._clone(symbols, domain, order)

    @cacheit
    def _clone(self, symbols, domain, order):
        # Cached clone implementation.
        return self.__class__(
            symbols or self.symbols, domain or self.domain, order or self.order
        )

    def compose(self, other):
        """Add the generators of other ring to this ring."""
        if self != other:
            syms = set(self.symbols).union(set(other.symbols))
            return self.clone(symbols=list(syms))
        else:
            return self

    # Domain conversions
    def to_domain(self):
        """Convert to a domain."""
        return PolynomialRing(self)

    def to_field(self):
        """Convert to a field of fractions."""
        from sympy.polys.fields import FracField

        return FracField(self.symbols, self.domain, self.order)

    def to_ground(self):
        """Convert to ground domain."""
        if self.domain.is_Composite or hasattr(self.domain, "domain"):
            return self.clone(domain=self.domain.domain)
        else:
            raise ValueError(f"{self.domain} is not a composite domain")

    # Element creation and testing
    def is_element(self, element):
        """Check if element belongs to this ring."""
        return isinstance(element, PolyElement) and element.ring == self

    def domain_new(self, element, orig_domain=None):
        """Create a new element of the ground domain."""
        return self.domain.convert(element, orig_domain)

    def ground_new(self, coeff):
        """Create a constant polynomial with given coefficient."""
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff):
        """Create a polynomial with a single term."""
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    # Polynomial creation from various formats
    def from_dict(self, element, orig_domain=None):
        """Create polynomial from dictionary of monomials to coefficients."""
        if not isinstance(element, dict):
            raise TypeError(
                "Input must be a dictionary mapping monomials to coefficients"
            )
        return self._from_dict_ground(element, orig_domain)

    def from_terms(self, element, orig_domain=None):
        """Create polynomial from sequence of (monomial, coefficient) pairs."""
        return self.from_dict(dict(element), orig_domain)

    def from_list(self, element):
        """Create polynomial from list(dense) representation."""
        return self.from_dict(dmp_to_dict(element, self.ngens - 1, self.domain))

    def from_expr(self, expr):
        """Create polynomial from SymPy expression."""
        mapping = dict(zip(self.symbols, self.gens))

        try:
            poly = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError(
                f"expected an expression convertible to a polynomial in {self}, "
                f"got {expr}"
            )
        else:
            return self.ring_new(poly)

    def _rebuild_expr(self, expr, mapping):
        # Rebuild expression as polynomial.
        domain = self.domain

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return reduce(add, map(_rebuild, expr.args))
            elif expr.is_Mul:
                return reduce(mul, map(_rebuild, expr.args))
            else:
                # Handle powers and other expressions
                base, exp = expr.as_base_exp()
                if exp.is_Integer and exp > 1:
                    return _rebuild(base) ** int(exp)
                else:
                    return self.ground_new(domain.convert(expr))

        return _rebuild(sympify(expr))

    # Generator operations
    def monomial_basis(self, i):
        """Return the i-th basis element."""
        basis = [0] * self.ngens
        basis[i] = 1
        return tuple(basis)

    def index(self, gen):
        """Get index of generator in the ring."""
        if gen is None:
            return 0 if self.ngens else -1  # Impossible choice indicator
        elif isinstance(gen, (int, str)):
            return self._gen_index(gen)
        elif self.is_element(gen):
            try:
                return self.gens.index(gen)
            except ValueError:
                raise ValueError(f"invalid generator: {gen}")
        else:
            raise ValueError(
                f"expected a polynomial generator, an integer, a string or None, "
                f"got {gen}"
            )

    def _gen_index(self, gen):
        # Get generator index from int or string.
        if isinstance(gen, int):
            if 0 <= gen < self.ngens:
                return gen
            elif -self.ngens <= gen <= -1:
                return gen + self.ngens
            else:
                raise ValueError(f"invalid generator index: {gen}")
        else:  # gen is a string
            try:
                return self.symbols.index(gen)
            except ValueError:
                raise ValueError(f"invalid generator: {gen}")

    def add_gens(self, symbols):
        """Add new generators to the ring."""
        syms = set(self.symbols).union(set(symbols))
        return self.clone(symbols=list(syms))

    def drop(self, *gens):
        """Remove specified generators from the ring."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def drop_to_ground(self, *gens):
        """Remove generators and inject them into the ground domain."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]
        gens_to_drop = [gen for i, gen in enumerate(self.gens) if i not in indices]

        if not symbols:
            return self
        else:
            return self.clone(symbols=symbols, domain=self.drop(*gens_to_drop))

    # Polynomial operations
    def add(self, *objs):
        """
        Add a sequence of polynomials or containers of polynomials.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> R, x = ring("x", ZZ)
        >>> R.add([ x**2 + 2*i + 3 for i in range(4) ])
        4*x**2 + 24
        >>> _.factor_list()
        (4, [(x**2 + 6, 1)])

        """
        result = self.zero

        for obj in objs:
            if is_sequence(obj, include=GeneratorType):
                result += self.add(*obj)
            else:
                result += obj

        return result

    def mul(self, *objs):
        """
        Multiply a sequence of polynomials or containers of polynomials.

        Examples
        ========

        >>> from sympy.polys.rings import ring
        >>> from sympy.polys.domains import ZZ

        >>> R, x = ring("x", ZZ)
        >>> R.mul([ x**2 + 2*i + 3 for i in range(4) ])
        x**8 + 24*x**6 + 206*x**4 + 744*x**2 + 945
        >>> _.factor_list()
        (1, [(x**2 + 3, 1), (x**2 + 5, 1), (x**2 + 7, 1), (x**2 + 9, 1)])

        """
        result = self.one

        for obj in objs:
            if is_sequence(obj, include=GeneratorType):
                result *= self.mul(*obj)
            else:
                result *= obj

        return result

    def symmetric_poly(self, n):
        """Return the elementary symmetric polynomial of degree n."""
        if n < 0 or n > self.ngens:
            raise ValueError(
                f"Cannot generate symmetric polynomial of order {n} for {self.gens}"
            )
        elif not n:
            return self.one
        else:
            poly = self.zero
            for s in subsets(range(self.ngens), int(n)):
                monom = tuple(int(i in s) for i in range(self.ngens))
                poly += self.term_new(monom, self.domain.one)
            return poly

    # Main element creation method
    def ring_new(self, element):
        """Create a ring element from various input types."""
        if isinstance(element, PolyElement):
            if self == element.ring:
                return element
            elif (
                isinstance(self.domain, PolynomialRing)
                and self.domain.ring == element.ring
            ):
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

    # Serialization support
    def __getstate__(self):
        state = self.__dict__.copy()

        # Remove function objects that can't be pickled
        state.pop("leading_expv", None)

        # Remove monomial operation functions
        for key in list(state.keys()):
            if key.startswith("monomial_"):
                state.pop(key, None)

        return state

    # Internal helper methods
    def _gens(self):
        # Generate the polynomial generators.
        one = self.domain.one
        generators = []

        for i in range(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            generators.append(poly)

        return tuple(generators)

    def _ring_equality(self, other):
        # Check equality of two polynomial rings.
        return (self.symbols, self.domain, self.ngens, self.order) == (
            other.symbols,
            other.domain,
            other.ngens,
            other.order,
        )

    def _from_dict_ground(self, element, orig_domain=None):
        # Create polynomial from dictionary with ground domain conversion.
        poly = self.zero
        domain_new = self.domain_new

        for monom, coeff in element.items():
            if coeff:  # Skip zero coefficients
                coeff = domain_new(coeff, orig_domain)
                poly[monom] = coeff

        return poly


class PolyElement(DomainElement, DefaultPrinting, CantSympify, dict):
    """Element of multivariate distributed polynomial ring."""

    def __init__(self, ring, init):
        super().__init__(init)
        self.ring = ring
        # This check would be too slow to run every time:
        # self._check()

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

    def __ne__(self, other):
        return not self == other

    def __pos__(self):
        return self

    def __lt__(self, other):
        return self._cmp(other, lt)

    def __le__(self, other):
        return self._cmp(other, le)

    def __gt__(self, other):
        return self._cmp(other, gt)

    def __ge__(self, other):
        return self._cmp(other, ge)

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
            if p2.ring == ring:
                return p1._add(p2)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring
            ):
                return p1._add_ground(p2)
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2._add_ground(p1)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._add_ground(cp2)

    def __radd__(p1, n):
        ring = p1.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._add_ground(n)

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
            if p2.ring == ring:
                return p1._sub(p2)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring
            ):
                return p1._sub_ground(p2)
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rsub__(p1)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._sub_ground(cp2)

    def __rsub__(self, n):
        ring = self.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            p = self.__neg__()
            return p._add_ground(n)

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
        if not p1 or not p2:
            return ring.zero

        if isinstance(p2, PolyElement):
            if p2.ring == ring:
                return p1._mul(p2)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring
            ):
                return p1.mul_ground(p2)
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rmul__(p1)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1.mul_ground(cp2)

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
        ring = p1.ring
        if isinstance(p2, PolyElement):
            try:
                cp2 = ring.ring_new(p2)
            except CoercionFailed:
                return NotImplemented
            else:
                return p1._mul(cp2)

        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1.mul_ground(cp2)

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

        if not n:
            if self:
                return self.ring.one
            else:
                raise ValueError("0**0")

        return self._pow_int(n)

    def __divmod__(p1, p2):
        ring = p1.ring
        if not p2:
            raise ZeroDivisionError("polynomial division")
        if isinstance(p2, PolyElement):
            if p2.ring == ring:
                return p1._divmod(p2)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring
            ):
                pass
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rdivmod__(p1)
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._divmod_ground(cp2)

    def __rdivmod__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2._divmod(p1)

    def __mod__(p1, p2):
        ring = p1.ring
        if not p2:
            raise ZeroDivisionError("polynomial division")
        if isinstance(p2, PolyElement):
            if p2.ring == ring:
                return p1._mod(p2)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring
            ):
                pass
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rmod__(p1)
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._mod_ground(cp2)

    def __rmod__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2._mod(p1)

    def __floordiv__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1._floordiv(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rtruediv__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._floordiv_ground(p2)

    def __rfloordiv__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2._floordiv(p1)

    def __truediv__(p1, p2):
        ring = p1.ring

        if not p2:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(p2):
            return p1._truediv(p2)
        elif isinstance(p2, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == p2.ring:
                pass
            elif (
                isinstance(p2.ring.domain, PolynomialRing)
                and p2.ring.domain.ring == ring
            ):
                return p2.__rtruediv__(p1)
            else:
                return NotImplemented

        try:
            p2 = ring.domain_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p1._floordiv_ground(p2)

    def __rtruediv__(p1, p2):
        ring = p1.ring
        try:
            p2 = ring.ground_new(p2)
        except CoercionFailed:
            return NotImplemented
        else:
            return p2._truediv(p1)

    @property
    def is_generator(self):
        return self in self.ring._gens_set

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

    @property
    def is_monic(f):
        return f.ring.domain.is_one(f.LC)

    @property
    def is_primitive(f):
        return f.ring.domain.is_one(f.content())

    @property
    def is_linear(f):
        return all(sum(monom) <= 1 for monom in f.itermonoms())

    @property
    def is_quadratic(f):
        return all(sum(monom) <= 2 for monom in f.itermonoms())

    def _check(self):
        """Validate polynomial structure (expensive - use sparingly)."""
        assert isinstance(self, PolyElement)
        assert isinstance(self.ring, PolyRing)

        dom = self.ring.domain
        assert isinstance(dom, Domain)

        for monom, coeff in self.iterterms():
            assert dom.of_type(coeff)
            assert len(monom) == self.ring.ngens
            assert all(isinstance(exp, int) and exp >= 0 for exp in monom)

    def new(self, init):
        """Create a new polynomial element in the same ring."""
        return self.__class__(self.ring, init)

    def parent(self):
        """Return the parent domain of this polynomial."""
        return self.ring.to_domain()

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

    def as_expr(self, *symbols):
        """Convert polynomial to SymPy expression."""
        if not symbols:
            symbols = self.ring.symbols
        elif len(symbols) != self.ring.ngens:
            raise ValueError(
                f"Wrong number of symbols, expected {self.ring.ngens} "
                f"got {len(symbols)}"
            )

        return expr_from_dict(self._as_expr_dict(), *symbols)

    def set_ring(self, new_ring):
        """Change the ring of this polynomial."""
        if self.ring == new_ring:
            return self
        return self._change_ring(new_ring)

    def strip_zero(self):
        """Eliminate monomials with zero coefficient."""
        for monom, coeff in self.listterms():
            if not coeff:
                del self[monom]

    def almosteq(self, other, tolerance=None):
        """Approximate equality test for polynomials."""
        ring = self.ring

        if ring.is_element(other):
            if set(self.itermonoms()) != set(other.itermonoms()):
                return False

            almosteq = ring.domain.almosteq
            for monom in self.itermonoms():
                if not almosteq(self[monom], other[monom], tolerance):
                    return False
            return True
        elif len(self) > 1:
            return False
        else:
            try:
                other = ring.domain.convert(other)
            except CoercionFailed:
                return False
            else:
                return ring.domain.almosteq(self.const(), other, tolerance)

    def sort_key(self):
        """Return a key for sorting polynomials."""
        return len(self), self.terms()

    def _drop(self, gen):
        ring = self.ring
        i = ring.index(gen)

        if ring.ngens == 1:
            return i, ring.domain
        else:
            new_ring = ring.drop(gen)
            return i, new_ring

    def drop(self, gen):
        i, ring = self._drop(gen)

        if self.ring.ngens == 1:
            if self.is_ground:
                return self.coeff(1)
            else:
                raise ValueError("Cannot drop %s" % gen)
        else:
            poly = ring.zero

            for k, v in self.iterterms():
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
            mon = monom[:i] + monom[i + 1 :]
            term = (gen ** monom[i]).mul_ground(coeff)
            if mon not in poly:
                poly[mon] = term
            else:
                poly[mon] = poly[mon] + term

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
        return self._square()

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
            return f._degree(i)

    def degrees(f):
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not f:
            return (ninf,) * f.ring.ngens
        else:
            return f._degrees()

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
            return (ninf,) * f.ring.ngens
        else:
            return tuple(map(min, list(zip(*f.itermonoms()))))

    def monic(f):
        """Divides all coefficients by the leading coefficient."""
        if not f:
            return f
        else:
            return f.quo_ground(f.LC)

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
        return self._div(fv)

    def quo_ground(f, x):
        domain = f.ring.domain

        if not x:
            raise ZeroDivisionError("polynomial division")
        if not f or x == domain.one:
            return f
        return f._quo_ground(x)

    def extract_ground(self, g):
        f = self
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def quo_term(f, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError("polynomial division")
        elif not f:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.quo_ground(coeff)
        return f._quo_term(term)

    def _norm(f, norm_func):
        if not f:
            return f.ring.domain.zero
        else:
            ground_abs = f.ring.domain.abs
            return norm_func([ground_abs(coeff) for coeff in f.itercoeffs()])

    def max_norm(f):
        return f._norm(max)

    def l1_norm(f):
        return f._norm(sum)

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

        return J, f._deflate(J, polys)

    def canonical_unit(f):
        domain = f.ring.domain
        return domain.canonical_unit(f.LC)

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
        return f._diff(i)

    def trunc_ground(f, p):
        if f.ring.domain.is_ZZ:
            terms = []

            for monom, coeff in f.iterterms():
                coeff = coeff % p

                if coeff > p // 2:
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [(monom, coeff % p) for monom, coeff in f.iterterms()]

        poly = f.new(terms)
        poly.strip_zero()
        return poly

    rem_ground = trunc_ground

    def lcm(self, g):
        f = self
        domain = f.ring.domain

        if not domain.is_Field:
            fc, f = f.primitive()
            gc, g = g.primitive()
            c = domain.lcm(fc, gc)

        h = (f * g).quo(f.gcd(g))

        if not domain.is_Field:
            return h.mul_ground(c)
        else:
            return h.monic()

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
        monoms = [m[:i] + (0,) + m[i + 1 :] for m in monoms]
        return p.ring.from_dict(dict(zip(monoms, coeffs)))

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
                raise ValueError(
                    "expected a generator, value pair a sequence of such pairs"
                )

        replacements = [(gens_map[x], ring.ring_new(g)) for x, g in replacements]

        return f._compose(replacements, initial_poly=poly)

    def __call__(f, *values):
        if 0 < len(values) <= f.ring.ngens:
            return f.evaluate(list(zip(f.ring.gens, values)))
        else:
            raise ValueError(
                "expected at least 1 and at most %s values, got %s"
                % (f.ring.ngens, len(values))
            )

    def evaluate(self, *args, **kwargs):
        eval_dict = {}
        ring = self.ring

        if len(args) == 1 and isinstance(args[0], list) and not kwargs:
            for gen, val in args[0]:
                idx = ring.index(gen)
                eval_dict[idx] = ring.domain.convert(val)

        elif len(args) == 2 and not kwargs:
            x, a = args
            idx = ring.index(x)
            eval_dict[idx] = ring.domain.convert(a)
        else:
            raise ValueError("Invalid arguments for evaluate()")

        if not eval_dict:
            return self
        elif len(eval_dict) == ring.ngens:
            return self._evaluate(eval_dict)
        else:
            temp_result = self._subs(eval_dict)
            new_ring = ring.drop(*[ring.gens[i] for i in eval_dict.keys()])
            return temp_result.set_ring(new_ring)

    def subs(self, *args, **kwargs):
        subs_dict = {}
        ring = self.ring

        if len(args) == 1 and isinstance(args[0], list) and not kwargs:
            for gen, val in args[0]:
                idx = ring.index(gen)
                subs_dict[idx] = ring.domain.convert(val)

        elif len(args) == 2 and not kwargs:
            x, a = args
            idx = ring.index(x)
            subs_dict[idx] = ring.domain.convert(a)
        else:
            raise ValueError("Invalid arguments for subs()")

        if not subs_dict:
            return self
        elif len(subs_dict) == ring.ngens:
            result = self._evaluate(subs_dict)
            return ring.ground_new(result)
        else:
            return self._subs(subs_dict)

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
        return self._prem(g, x)

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
        return self._pdiv(g, x)

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
        return self._pquo(g, x)

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
        return self._pexquo(g, x)

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
        return self._subresultants(g, x)

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
        return self._symmetrize()

    # <--------------------------INTERFACE / IMPLEMENTATION DIVIDER-------------------------->#

    def __eq__(self, other):
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
        if not other:
            return not self
        elif self.ring.is_element(other):
            return dict.__eq__(self, other)
        elif len(self) > 1:
            return False
        else:
            return self.get(self.ring.zero_monom) == other

    def __neg__(self):
        # Return (-1) * self in case of python-flint
        return self.new([(monom, -coeff) for monom, coeff in self.iterterms()])

    def _add(p1, p2):
        p = p1.copy()
        get = p.get
        zero = p1.ring.domain.zero
        for k, v in p2.items():
            v = get(k, zero) + v
            if v:
                p[k] = v
            else:
                del p[k]
        return p

    def _add_ground(p1, cp2):
        p = p1.copy()
        if not cp2:
            return p
        ring = p1.ring
        zm = ring.zero_monom
        v = p1.get(zm, ring.domain.zero) + cp2
        if v:
            p[zm] = v
        else:
            del p[zm]
        return p

    def _sub(p1, p2):
        p = p1.copy()
        get = p.get
        zero = p1.ring.domain.zero
        for k, v in p2.items():
            v = get(k, zero) - v
            if v:
                p[k] = v
            else:
                del p[k]
        return p

    def _sub_ground(p1, cp2):
        p = p1.copy()
        #if not cp2:
            #return p
        ring = p1.ring
        zm = ring.zero_monom
        v = p1.get(zm, ring.domain.zero) - cp2
        if v:
            p[zm] = v
        else:
            del p[zm]
        return p

    def _mul(p1, p2):
        ring = p1.ring
        p = ring.zero
        for exp1, v1 in p1.iterterms():
            for exp2, v2 in p2.iterterms():
                exp = ring.monomial_mul(exp1, exp2)
                v = v1 * v2
                p[exp] = p.get(exp, ring.domain.zero) + v
        p.strip_zero()
        return p

    def mul_ground(f, x):
        if not x:
            return f.ring.zero

        terms = [(monom, coeff * x) for monom, coeff in f.iterterms()]
        return f.new(terms)

    def _pow_int(self, n):
        if n == 1:
            return self.copy()
        elif n == 2:
            return self.square()
        elif n == 3:
            return self * self.square()
        elif len(self) <= 5:  # TODO: use an actual density measure
            return self._pow_multinomial(n)
        else:
            return self._pow_generic(n)

    def _pow_generic(self, n):
        p = self.ring.one
        c = self

        while True:
            if n & 1:
                p = p * c
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

    def _square(self):
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
                p[exp] = get(exp, zero) + pk * self[k2]
        p = p.imul_num(2)
        get = p.get
        for k, v in self.items():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, zero) + v**2
        p.strip_zero()
        # p._check()
        return p

    def _divmod(self, other):
        return self.div(other)

    def _divmod_ground(self, x):
        return self.quo_ground(x), self.rem_ground(x)

    def _floordiv(self, p2):
        return self.quo(p2)

    def _floordiv_ground(self, p2):
        return self.quo_ground(p2)

    def _truediv(self, p2):
        return self.exquo(p2)

    def _term_div(self):
        zm = self.ring.zero_monom
        domain = self.ring.domain
        domain_quo = domain.quo
        monomial_div = self.ring.monomial_div

        if domain.is_Field:

            def term_div(a_lm_a_lc, b_lm_b_lc):
                a_lm, a_lc = a_lm_a_lc
                b_lm, b_lc = b_lm_b_lc
                if b_lm == zm:  # apparently this is a very common case
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
                if b_lm == zm:  # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if not (monom is None or a_lc % b_lc):
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return None

        return term_div

    def rem(self, G):
        f = self
        if isinstance(G, PolyElement):
            G = [G]
        if not all(G):
            raise ZeroDivisionError("polynomial division")
        return f._rem(G)

    def quo(f, G):
        return f.div(G)[0]

    def exquo(f, G):
        q, r = f.div(G)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(f, G)

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
            coeff = get(ka, zero) + v * c
            if coeff:
                p1[ka] = coeff
            else:
                del p1[ka]
        return p1

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
            return p * c
        if not c:
            p.clear()
            return
        for exp in p:
            p[exp] *= c
        return p

    def _rem(f, G):
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
                        c1 = get(m1, zero) - c * cg
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

    def _mod(self, other):
        return self.rem(other)

    def _mod_ground(self, x):
        return self.rem_ground(x)

    @property
    def is_ground(self):
        # Return self.flint_poly.is_constant() in case of python-flint
        return not self or (len(self) == 1 and self.ring.zero_monom in self)

    @property
    def is_zero(f):
        # Return self.flint_poly.is_zero() in case of python-flint
        return not f

    @property
    def is_one(f):
        # Return self.flint_poly.is_one() in case of python-flint
        return f == f.ring.one

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

    @property
    def LC(self):
        # Just use leafing_coefficient() in case of python-flint
        return self._get_coeff(self.leading_expv())

    @property
    def LM(self):
        # Use monomial(0) in case of python-flint
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom
        else:
            return expv

    @property
    def LT(self):
        # Use monomial(0) and leafing_coefficient() in case of python-flint
        expv = self.leading_expv()
        if expv is None:
            return (self.ring.zero_monom, self.ring.domain.zero)
        else:
            return (expv, self._get_coeff(expv))

    def clear_denoms(self):
        """Clear denominators from polynomial coefficients."""
        domain = self.ring.domain

        if not domain.is_Field or not domain.has_assoc_Ring:
            return domain.one, self

        ground_ring = domain.get_ring()
        common = ground_ring.one
        lcm = ground_ring.lcm
        denom = domain.denom

        for coeff in self.values():
            common = lcm(common, denom(coeff))

        poly = self.new([(monom, coeff * common) for monom, coeff in self.items()])
        return common, poly

    def _change_ring(self, new_ring):
        # Use fmpz_mpoly.compose() or fmpz_mpoly.compose() in case of python-flint
        if self.ring.symbols != new_ring.symbols:
            terms = list(zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols)))
            return new_ring.from_terms(terms, self.ring.domain)
        else:
            return new_ring.from_dict(self, self.ring.domain)

    def _as_expr_dict(self):
        # Can just use self.flint_poly.to_dict() in case of python-flint
        # Or this can just directly go into the baseclass as is since iterterms
        # will be implemented separately for pure python and flint versions anyways
        to_sympy = self.ring.domain.to_sympy
        return {monom: to_sympy(coeff) for monom, coeff in self.iterterms()}

    def _cmp(self, other, op):
        # We can override this for python-flint version
        # to use the native lt, le, gt, ge methods
        if self.ring.is_element(other):
            return op(self.sort_key(), other.sort_key())
        else:
            return NotImplemented

    def to_dense(self):
        return dmp_from_dict(self, self.ring.ngens - 1, self.ring.domain)

    def to_dict(self):
        # Return a self.flint_poly.to_dict() in case of python-flint
        return dict(self)

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        # Use str(self.flint_poly).replace("^", "**") in case of python-flint
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
                    scoeff = ""
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
                    sexpv.append("%s" % symbol)
            if scoeff:
                sexpv = [scoeff] + sexpv
            sexpvs.append(mul_symbol.join(sexpv))
        if sexpvs[0] in [" + ", " - "]:
            head = sexpvs.pop(0)
            if head == " - ":
                sexpvs.insert(0, "-")
        return "".join(sexpvs)

    def _degree(self, i):
        return max(monom[i] for monom in self.itermonoms())

    def _degrees(self):
        return tuple(map(max, list(zip(*self.itermonoms()))))

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
        # Use fmpz_mpoly.monomial(1) or fmpq_mpoly.monomial(1) in case of python-flint
        if self:
            return self.ring.leading_expv(self)
        else:
            return None

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def const(self):
        # Use
        """Returns the constant coefficient."""
        return self._get_coeff(self.ring.zero_monom)

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

    def _sorted(self, seq, order):
        if order is None:
            order = self.ring.order
        else:
            order = OrderOpt.preprocess(order)

        if order is lex:
            return sorted(seq, key=lambda monom: monom[0], reverse=True)
        else:
            return sorted(seq, key=lambda monom: order(monom[0]), reverse=True)

    def itercoeffs(self):
        """Iterator over coefficients of a polynomial."""
        return iter(self.values())

    def itermonoms(self):
        """Iterator over monomials of a polynomial."""
        return iter(self.keys())

    def iterterms(self):
        """Iterator over terms of a polynomial."""
        return iter(self.items())

    def listcoeffs(self):
        """Unordered list of polynomial coefficients."""
        return list(self.values())

    def listmonoms(self):
        """Unordered list of polynomial monomials."""
        return list(self.keys())

    def listterms(self):
        """Unordered list of polynomial terms."""
        return list(self.items())

    def content(f):
        """Returns GCD of polynomial's coefficients."""
        # In the flint version, we will have to override
        # this to use the native content() method for ZZ
        # and use the pure python technique for other domains
        domain = f.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in f.itercoeffs():
            cont = gcd(cont, coeff)

        return cont

    def primitive(f):
        """Returns content and a primitive polynomial."""
        cont = f.content()
        if cont == f.ring.domain.zero:
            return (cont, f)
        return cont, f.quo_ground(cont)

    def mul_monom(f, monom):
        monomial_mul = f.ring.monomial_mul
        terms = [
            (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in f.items()
        ]
        return f.new(terms)

    def mul_term(f, term):
        monom, coeff = term

        if not f or not coeff:
            return f.ring.zero
        elif monom == f.ring.zero_monom:
            return f.mul_ground(coeff)

        monomial_mul = f.ring.monomial_mul
        terms = [
            (monomial_mul(f_monom, monom), f_coeff * coeff)
            for f_monom, f_coeff in f.items()
        ]
        return f.new(terms)

    def _quo_ground(self, x):
        domain = self.ring.domain
        if domain.is_Field:
            quo = domain.quo
            terms = [(monom, quo(coeff, x)) for monom, coeff in self.iterterms()]
        else:
            terms = [
                (monom, coeff // x)
                for monom, coeff in self.iterterms()
                if not (coeff % x)
            ]
        return self.new(terms)

    def _quo_term(self, term):
        term_div = self._term_div()
        terms = [term_div(t, term) for t in self.iterterms()]
        return self.new([t for t in terms if t is not None])

    def _deflate(self, J, polys):
        ring = self.ring
        H = []
        for p in polys:
            h = ring.zero
            for I, coeff in p.iterterms():
                N = [i // j for i, j in zip(I, J)]
                h[tuple(N)] = coeff
            H.append(h)
        return H

    def inflate(f, J):
        poly = f.ring.zero

        for I, coeff in f.iterterms():
            N = [i * j for i, j in zip(I, J)]
            poly[tuple(N)] = coeff

        return poly

    def gcd(f, g):
        return f.cofactors(g)[0]

    def _diff(self, i):
        # Use the native derivative() method in case of python-flint
        ring = self.ring
        m = ring.monomial_basis(i)
        g = ring.zero
        for expv, coeff in self.iterterms():
            if expv[i]:
                e = ring.monomial_ldiv(expv, m)
                g[e] = ring.domain_new(coeff * expv[i])
        return g

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
        mf, cf = f.listterms()[0]
        _mgcd, _cgcd = mf, cf
        for mg, cg in g.iterterms():
            _mgcd = monomial_gcd(_mgcd, mg)
            _cgcd = ground_gcd(_cgcd, cg)
        h = f.new([(_mgcd, _cgcd)])
        cff = f.new([(monomial_ldiv(mf, _mgcd), ground_quo(cf, _cgcd))])
        cfg = f.new(
            [
                (monomial_ldiv(mg, _mgcd), ground_quo(cg, _cgcd))
                for mg, cg in g.iterterms()
            ]
        )
        return h, cff, cfg

    def _gcd(f, g):
        ring = f.ring

        if ring.domain.is_QQ:
            return f._gcd_QQ(g)
        elif ring.domain.is_ZZ:
            return f._gcd_ZZ(g)
        else:  # TODO: don't use dense representation (port PRS algorithms)
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

    def _compose(self, replacements, initial_poly):
        ring = self.ring
        poly = initial_poly

        for monom, coeff in self.iterterms():
            monom = list(monom)
            subpoly = ring.one

            for i, g in replacements:
                n, monom[i] = monom[i], 0
                if n:
                    subpoly *= g**n

            subpoly = subpoly.mul_term((tuple(monom), coeff))
            poly += subpoly

        return poly

    def _div(self, fv):
        # Implement the same algorithm from [CLO] p64. in python-flint
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
                raise ValueError("self and f must have the same ring")
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
                expv = p.leading_expv()
                r = r._iadd_monom((expv, p[expv]))
                del p[expv]
        if expv == ring.zero_monom:
            r += p
        if ret_single:
            if not qv:
                return ring.zero, r # This never happens since this branch can never be True
            else:
                return qv[0], r
        else:
            return qv, r

    # The following _p* and _subresultants methods can just be converted to pure python
    # methods in case of python-flint since their speeds don't exactly matter wrt the
    # flint version.
    def _prem(self, g, x):
        f = self
        x = f.ring.index(x)
        df = f.degree(x)
        dg = g.degree(x)

        if dg < 0:
            raise ZeroDivisionError("polynomial division")

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
            G = g * lc_r * xp**j
            r = R - G

            dr = r.degree(x)

            if dr < dg:
                break

        c = lc_g**N

        return r * c

    def _pdiv(self, g, x):
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

            q = Q + (lc_r) * xp**j

            R = r * lc_g

            G = g * lc_r * xp**j

            r = R - G

            dr = r.degree(x)

            if dr < dg:
                break

        c = lc_g**N

        q = q * c
        r = r * c

        return q, r

    def _pquo(self, g, x):
        f = self
        return f.pdiv(g, x)[0]

    def _pexquo(self, g, x):
        f = self
        q, r = f.pdiv(g, x)

        if r.is_zero:
            return q
        else:
            raise ExactQuotientFailed(f, g)

    def _subresultants(self, g, x):
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

        c = lc**d

        S = [1, c]

        c = -c

        while h:
            k = h.degree(x)

            R.append(h)
            f, g, m, d = g, h, k, m - k

            b = -lc * c**d
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

    def _subs(self, subs_dict):
        ring = self.ring
        result_poly = ring.zero

        for monom, coeff in self.iterterms():
            new_coeff = coeff
            new_monom = list(monom)

            for i, val in subs_dict.items():
                exp = monom[i]
                if exp > 0:
                    new_coeff *= val**exp
                new_monom[i] = 0

            if new_coeff:
                new_monom = tuple(new_monom)
                if new_monom in result_poly:
                    result_poly[new_monom] += new_coeff
                    if not result_poly[new_monom]:
                        del result_poly[new_monom]
                else:
                    result_poly[new_monom] = new_coeff

        return result_poly

    def _evaluate(self, eval_dict):
        result = self.ring.domain.zero

        for monom, coeff in self.iterterms():
            monom_value = self.ring.domain.one
            for i, exp in enumerate(monom):
                if exp > 0:
                    monom_value *= eval_dict[i] ** exp

            result += coeff * monom_value

        return result

    def _symmetrize(self):
        # For the python-flint override this can just be converted back to
        # the pure python version until python-flint provides some
        # equivalent functionality.
        f = self.copy()
        ring = f.ring
        n = ring.ngens

        if not n:
            return f, ring.zero, []

        polys = [ring.symmetric_poly(i + 1) for i in range(n)]

        poly_powers = {}

        def get_poly_power(i, n):
            if (i, n) not in poly_powers:
                poly_powers[(i, n)] = polys[i] ** n
            return poly_powers[(i, n)]

        indices = list(range(n - 1))
        weights = list(range(n, 0, -1))

        symmetric = ring.zero

        while f:
            _height, _monom, _coeff = -1, None, None

            for i, (monom, coeff) in enumerate(f.terms()):
                if all(monom[i] >= monom[i + 1] for i in indices):
                    height = max(n * m for n, m in zip(weights, monom))

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

    # TODO: following methods should point to polynomial
    # representation independent algorithm implementations.

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

