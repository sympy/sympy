"""Sparse polynomial rings."""

from __future__ import annotations

from typing import Generic, Callable, Iterable, TYPE_CHECKING, Any

from operator import add, mul, lt, le, gt, ge
from functools import reduce
from types import GeneratorType

from sympy.external.gmpy import GROUND_TYPES

from sympy.core.cache import cacheit
from sympy.core.expr import Expr
from sympy.core.intfunc import igcd
from sympy.core.symbol import Symbol, symbols as _symbols
from sympy.core.sympify import CantSympify, sympify
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.compatibility import IPolys
from sympy.polys.constructor import construct_domain
from sympy.polys.densebasic import ninf, dmp_to_dict, dmp_from_dict
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domain import Domain, Er, Es
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.polynomialring import PolynomialRing
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.monomials import MonomialOps
from sympy.polys.orderings import lex, MonomialOrder, LexOrder, GradedLexOrder, ReversedGradedLexOrder
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


if GROUND_TYPES == 'flint':
    import flint
    def _supports_flint(domain, order):
        supported_domains = domain.is_ZZ or domain.is_QQ
        supported_orders = (LexOrder, GradedLexOrder, ReversedGradedLexOrder)
        supported_order = isinstance(order, supported_orders)
        return supported_domains and supported_order

else:
    flint = None
    def _supports_flint(domain, order):
        return False

#flint = None

if TYPE_CHECKING:
    from typing import TypeIs
    from sympy.polys.fields import FracField


Mon = tuple[int, ...]


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
    >>> type(_) # doctest: +SKIP
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
    >>> type(_) # doctest: +SKIP
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
    >>> type(_)  # doctest: +SKIP
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
    >>> type(_) # doctest: +SKIP
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


class PolyRing(DefaultPrinting, IPolys, Generic[Er]):
    """Multivariate distributed polynomial rings."""

    symbols: tuple[Expr, ...]
    gens: tuple[PolyElement[Er], ...]
    ngens: int
    _gens_set: set[PolyElement]
    domain: Domain[Er]
    order: MonomialOrder
    _hash: int
    _hash_tuple: tuple
    _one: list[tuple[Mon, Er]]

    dtype: Callable[[Iterable[tuple[Mon, Er]] | dict[Mon, Er]], PolyElement[Er]]

    monomial_mul: Callable[[Mon, Mon], Mon]
    monomial_pow: Callable[[Mon, int], Mon]
    monomial_mulpow: Callable[[Mon, int, int], Mon]
    monomial_ldiv: Callable[[Mon, Mon], Mon]
    monomial_div: Callable[[Mon, Mon], Mon]
    monomial_lcm: Callable[[Mon, Mon], Mon]
    monomial_gcd: Callable[[Mon, Mon], Mon]
    leading_expv: Callable[[PolyElement[Er]], Mon]
    zero_monom: Mon

    def __new__(cls,
                symbols,
                domain: Domain[Er],
                order: str | MonomialOrder | None = lex
            ) -> PolyRing[Er]:

        # validate inputs
        symbols = tuple(_parse_symbols(symbols))
        preprocessed_domain = DomainOpt.preprocess(domain)
        morder = OrderOpt.preprocess(order)

        # Validate that symbols do not overlap with domain symbols
        if isinstance(preprocessed_domain, CompositeDomain) and set(symbols) & set(preprocessed_domain.symbols):
            raise GeneratorsError("polynomial ring and it's ground domain share generators")

        # delegation to appropriate subclass
        if flint is not None and _supports_flint(preprocessed_domain, morder):
            return FPolyRing._new(symbols, preprocessed_domain, morder)
        else:
            return PyPolyRing._new(symbols, preprocessed_domain, morder)


    @classmethod
    def _new(cls, symbols, domain, order):
        raise NotImplementedError("Must be implemented by subclass")

    def _init_instance(self, symbols, domain, order):
        # initialize the instance with common attributes.
        self.symbols = symbols
        self.ngens = len(symbols)
        self.domain = domain
        self.order = order
        self._hash_tuple = ("PolyRing", symbols, self.ngens, domain, order)
        self._hash = hash(self._hash_tuple)

        # set up polynomial creation and basic elements
        self.dtype = PolyElement(self, ()).new
        self.zero_monom = (0,) * self.ngens
        self.gens = self._gens()
        self._gens_set = set(self.gens)
        self._one = [(self.zero_monom, domain.one)]

        # initialize monomial operations
        self._init_monomial_operations()

        # set up leading exponent vector function
        self._init_leading_expv_function(order)

        # add generator attributes for Symbol names
        self._add_generator_attributes()

    def _init_monomial_operations(self):
        # initialize monomial operations based on number of generators
        if self.ngens:
            # operations for rings with at least one generator
            codegen = MonomialOps(self.ngens)
            self.monomial_mul = codegen.mul()
            self.monomial_pow = codegen.pow()
            self.monomial_mulpow = codegen.mulpow()
            self.monomial_ldiv = codegen.ldiv()
            self.monomial_div = codegen.div()
            self.monomial_lcm = codegen.lcm()
            self.monomial_gcd = codegen.gcd()
        else:
            # rings with no generator
            monunit = lambda a, b: ()
            self.monomial_mul = monunit
            self.monomial_pow = monunit
            self.monomial_mulpow = lambda a, b, c: ()
            self.monomial_ldiv = monunit
            self.monomial_div = monunit
            self.monomial_lcm = monunit
            self.monomial_gcd = monunit

    def _init_leading_expv_function(self, order):
        if order is lex:
            self.leading_expv = max
        else:
            self.leading_expv = lambda f: max(f, key=order)

    def _add_generator_attributes(self):
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
        return isinstance(other, PolyRing) and (self.symbols, self.domain, self.ngens, self.order) == (
            other.symbols,
            other.domain,
            other.ngens,
            other.order,
        )

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, key):
        """Get a subring with subset of symbols."""
        symbols = self.symbols[key]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def from_python(self, py_poly): # to satisfy mypy
        # makes sense in FPolyRing
        raise NotImplementedError

    # Properties
    @property
    def zero(self) -> PolyElement[Er]:
        """The zero polynomial."""
        return self.dtype([])

    @property
    def one(self) -> PolyElement[Er]:
        """The unit polynomial."""
        return self.dtype(self._one)

    @property
    def is_univariate(self) -> bool:
        """True if this is a univariate ring."""
        return self.ngens == 1

    @property
    def is_multivariate(self) -> bool:
        """True if this is a multivariate ring."""
        return self.ngens > 1

    # Ring operations and cloning
    def clone(self, symbols=None, domain=None, order=None) -> PolyRing:
        """Create a clone with modified parameters."""
        # convert list to tuple for hashability
        if symbols is not None and isinstance(symbols, list):
            symbols = tuple(symbols)
        return self._clone(symbols, domain, order)

    @cacheit
    def _clone(self, symbols, domain, order):
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
    def to_domain(self) -> PolynomialRing[Er]:
        """Convert to a domain."""
        return PolynomialRing(self)

    def to_field(self) -> FracField[Er]:
        """Convert to a field of fractions."""
        from sympy.polys.fields import FracField
        return FracField(self.symbols, self.domain, self.order)

    def to_ground(self) -> PolyRing:
        """Convert to ground domain."""
        if isinstance(self.domain, CompositeDomain) or hasattr(self.domain, 'domain'):
            return self.clone(domain=self.domain.domain)
        else:
            raise ValueError(f"{self.domain} is not a composite domain")

    # Element creation and testing
    def is_element(self, element) -> TypeIs[PolyElement[Er]]:
        """Check if element belongs to this ring."""
        return isinstance(element, PolyElement) and element.ring == self

    def domain_new(self, element, orig_domain=None) -> Er:
        """Create a new element of the ground domain."""
        return self.domain.convert(element, orig_domain)

    def ground_new(self, coeff) -> PolyElement[Er]:
        """Create a constant polynomial with given coefficient."""
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff) -> PolyElement[Er]:
        """Create a polynomial with a single term."""
        raise NotImplementedError

    # Polynomial creation from various formats
    def from_dict(self, element, orig_domain=None) -> PolyElement[Er]:
        """Create polynomial from dictionary of monomials to coefficients."""
        if not isinstance(element, dict):
            raise TypeError(
                "Input must be a dictionary mapping monomials to coefficients"
            )
        return self._from_dict_ground(element, orig_domain)

    def from_terms(self, element, orig_domain=None) -> PolyElement[Er]:
        """Create polynomial from sequence of (monomial, coefficient) pairs."""
        return self.from_dict(dict(element), orig_domain)

    def from_list(self, element) -> PolyElement[Er]:
        """Create polynomial from list(dense) representation."""
        return self.from_dict(dmp_to_dict(element, self.ngens - 1, self.domain))

    def from_expr(self, expr) -> PolyElement[Er]:
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

    def _rebuild_expr(self, expr, mapping) -> PolyElement[Er]:
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
    def monomial_basis(self, i) -> tuple[int, ...]:
        """Return the i-th basis element."""
        basis = [0] * self.ngens
        basis[i] = 1
        return tuple(basis)

    def index(self, gen) -> int:
        """Get index of generator in the ring."""
        if gen is None:
            return 0 if self.ngens else -1  # Impossible choice indicator
        if isinstance(gen, int):
            if 0 <= gen < self.ngens:
                return self._gen_index(gen)
            else:
                raise ValueError(f"invalid generator index: {gen}")

        elif isinstance(gen, str):
            if gen in self.symbols:
                return self._gen_index(gen)
            else:
                raise ValueError(f"invalid generator name: {gen}")

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

    def drop(self, *gens):
        """Remove specified generator(s) from the ring."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def add_gens(self, symbols):
        """Add new generator(s) to the ring."""
        syms = set(self.symbols).union(set(symbols))
        return self.clone(symbols=list(syms))

    # Polynomial operations
    def add(self, *objs):
        """Add a sequence of polynomials or containers of polynomials."""
        result = self.zero

        for obj in objs:
            if is_sequence(obj, include=GeneratorType) and not isinstance(obj, PolyElement):
                result += self.add(*obj)
            else:
                result += obj

        return result

    def mul(self, *objs):
        """Multiply a sequence of polynomials or containers of polynomials."""
        result = self.one

        for obj in objs:
            if is_sequence(obj, include=GeneratorType)and not isinstance(obj, PolyElement):
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
    def ring_new(self, element) -> PolyElement[Er]:
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
        del state["leading_expv"]

        for key in state:
            if key.startswith("monomial_"):
                del state[key]

        return state

    # internal helpers
    def _gens(self):
        """Generate the polynomial generators."""
        one = self.domain.one
        generators = []

        for i in range(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            generators.append(poly)

        return tuple(generators)

    # to be overridden
    def _gen_index(self, gen):
        # get generator index from int or str
        raise NotImplementedError("Must be implemented by subclass")

    def _from_dict_ground(self, element, orig_domain=None):
        # create polynomial from dict
        raise NotImplementedError("Must be implemented by subclass")

    def drop_to_ground(self, *gens):
        """Remove generators and inject them into the ground domain."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]
        gens_to_drop = [gen for i, gen in enumerate(self.gens) if i not in indices]

        if not symbols:
            return self
        else:
            return self.clone(symbols=symbols, domain=self.drop(*gens_to_drop))


class PyPolyRing(PolyRing):
    """Multivariate distributed polynomial rings."""

    @classmethod
    def _new(cls, symbols, domain, order):
        obj = object.__new__(cls)
        obj._init_instance(symbols, domain, order)
        return obj

    def _gen_index(self, gen):
        # get generator index from int or str
        return gen if isinstance(gen, int) else self.symbols.index(gen)

    def _from_dict_ground(self, element, orig_domain=None):
        # create polynomial from dict
        poly = self.zero
        domain_new = self.domain_new

        for monom, coeff in element.items():
            if coeff:  # Skip zero coefficients
                coeff = domain_new(coeff, orig_domain)
                poly[monom] = coeff

        return poly

    def term_new(self, monom, coeff) -> PolyElement[Er]:
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    def from_flint(self, flint_element):
        """Convert a FLINT polynomial(FPolyElement) to a PyPolyElement."""
        if not isinstance(flint_element, FPolyElement):
            raise TypeError("Expected FPolyElement, got %s" % type(flint_element))

        if (flint_element.ring.symbols != self.symbols or
            flint_element.ring.domain != self.domain or
            flint_element.ring.order != self.order):
            raise TypeError("Ring mismatch - cannot convert between incompatible rings")

        return self.from_dict(dict(flint_element.items()))


class FPolyRing(PolyRing):
    """Multivariate distributed polynomial rings backed by FLINT."""

    _python_ring: PyPolyRing  # type annotation to satisfy mypy
    flint_ctx: Any  # putting placeholder for flint context

    @classmethod
    def _new(cls, symbols, domain, order):
        obj = object.__new__(cls)

        obj.flint_ctx = None
        flint_names = cls._flint_mapping(symbols)
        str_order = cls._get_flint_order(order)

        if domain.is_ZZ:
            obj.flint_ctx = flint.fmpz_mpoly_ctx.get(flint_names, str_order)
        elif domain.is_QQ:
            obj.flint_ctx = flint.fmpq_mpoly_ctx.get(flint_names, str_order)
        else:
            raise NotImplementedError(f"Unsupported domain for FPolyRing: {domain}")

        obj._init_instance(symbols, domain, order)
        obj._python_ring = PyPolyRing._new(symbols, domain, order)

        return obj

    @staticmethod
    def _flint_mapping(symbols):
        flint_names = []
        used_names = set()

        for sym in symbols:
            base_name = str(sym)
            flint_name = base_name

            counter = 1
            while flint_name in used_names:
                flint_name = f"{base_name}_{counter}"
                counter += 1

            used_names.add(flint_name)
            flint_names.append(flint_name)

        return flint_names

    @staticmethod
    def _get_flint_order(order):
        """Convert SymPy monomial order to FLINT order string."""
        if isinstance(order, LexOrder):
            return "lex"
        elif isinstance(order, GradedLexOrder):
            return "deglex"
        else:
            return "degrevlex"

    def _gen_index(self, gen):
        # get generator index from int or str
        return self.flint_ctx.variable_to_index(gen)

    def _from_dict_ground(self, element, orig_domain=None):
        """Create polynomial from dict."""
        # For now, delegate to Python implementation
        # TODO: Implement direct polynomial creation from dict once FPolyElement is ready
        return self._python_ring._from_dict_ground(element, orig_domain)

    def term_new(self, monom, coeff) -> PolyElement[Er]:
        # TODO: Implement direct polynomial creation from term once FPolyElement is ready
        return self._python_ring.term_new(monom, coeff)

    def from_python(self, python_element):
        if not isinstance(python_element, PyPolyElement):
            raise TypeError("Expected PyPolyElement, got %s" % type(python_element))
        elif not _supports_flint(self.domain, self.order):
            raise NotImplementedError("Expected supported domain and order for FLINT conversion, got %s and %s" % self.domain, self.order)

        if (python_element.ring.symbols != self.symbols or
            python_element.ring.domain != self.domain or
            python_element.ring.order != self.order):
            raise TypeError("Ring mismatch, cannot convert between incompatible rings")

        return FPolyElement._new(self, python_element)


class PolyElement(DomainElement, DefaultPrinting, CantSympify, Generic[Er]):
    """Element of multivariate distributed polynomial ring."""

    def __new__(cls, ring, init):
        # Delegation to appropriate subclass based on ring type
        if isinstance(ring, FPolyRing):
            return FPolyElement._new(ring, init)
        else:
            return PyPolyElement._new(ring, init)

    @classmethod
    def _new(cls, ring, init):
        """Create a new polynomial element - to be implemented by subclasses."""
        raise NotImplementedError("Must be implemented by subclass")

    def __init__(self, ring: PolyRing[Er], init: dict[Mon, Er] | Iterable[tuple[Mon, Er]]):
        super().__init__()
        self.ring = ring
        # This check would be too slow to run every time:
        # self._check()

    def __getnewargs__(self):
        return (self.ring, list(self.iterterms()))

    def __setitem__(self, monom, coeff):
        raise NotImplementedError

    def __getitem__(self, monom):
        raise NotImplementedError

    def __len__(self):
        raise NotImplementedError

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

    def __pos__(self) -> PolyElement[Er]:
        return self

    def __lt__(self, other):
        return self._cmp(other, lt)

    def __le__(self, other):
        return self._cmp(other, le)

    def __gt__(self, other):
        return self._cmp(other, gt)

    def __ge__(self, other):
        return self._cmp(other, ge)

    def _get_python_ring(self) -> PyPolyRing: # to satisfy mypy
        # makes sense in FPolyElement
        raise NotImplementedError

    def as_expr(self, *symbols):
        if not symbols:
            symbols = self.ring.symbols
        elif len(symbols) != self.ring.ngens:
            raise ValueError(
                "Wrong number of symbols, expected %s got %s" %
                (self.ring.ngens, len(symbols))
            )
        return expr_from_dict(self.as_expr_dict(), *symbols)

    def __add__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        """Add two polynomials.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring

        >>> _, x, y = ring('x, y', ZZ)
        >>> (x + y)**2 + (x - y)**2
        2*x**2 + 2*y**2

        """
        if not other:
            return self.copy()

        ring = self.ring

        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._add(other)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring
            ):
                return self._add_ground(other)
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other._add_ground(self)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._add_ground(cp2)

    def __radd__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring
        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._add_ground(other)

    def __sub__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
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
        if not other:
            return self.copy()

        ring = self.ring

        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._sub(other)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring
            ):
                return self._sub_ground(other)
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rsub__(self)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._sub_ground(cp2)

    def __rsub__(self, n: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            p = self.__neg__()
            return p._add_ground(n)

    def __mul__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
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
        ring = self.ring
        if not self or not other:
            return ring.zero

        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._mul(other)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring
            ):
                return self.mul_ground(other)
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rmul__(self)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.mul_ground(cp2)

    def __rmul__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
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
        ring = self.ring
        if isinstance(other, PolyElement):
            try:
                p2 = ring.ring_new(other)
            except CoercionFailed:
                return NotImplemented # unreachable
            else:
                return self._mul(p2)

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.mul_ground(cp2)

    def __pow__(self, n: int) -> PolyElement[Er]:
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
        acceptable_exps: tuple[type, ...]  # to satisfy mypy
        if GROUND_TYPES == "flint":
            acceptable_exps = (int, flint.fmpz)
        else:
            acceptable_exps = (int,)
        if not isinstance(n, acceptable_exps):
            raise TypeError("exponent must be an integer, got %s" % n)
        elif n < 0:
            raise ValueError("exponent must be a non-negative integer, got %s" % n)

        if not n:
            if self:
                return self.ring.one
            else:
                raise ValueError("0**0")

        return self._pow_int(n)

    def __divmod__(self, other: PolyElement[Er] | Er | int) -> tuple[PolyElement[Er], PolyElement[Er]]:
        ring = self.ring
        if not other:
            raise ZeroDivisionError("polynomial division")
        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._divmod(other)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rdivmod__(self)
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._divmod_ground(cp2)

    def __rdivmod__(self, other: PolyElement[Er] | Er | int) -> tuple[PolyElement[Er], PolyElement[Er]]:
        ring = self.ring
        try:
            other = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other._divmod(self)

    def __mod__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring
        if not other:
            raise ZeroDivisionError("polynomial division")
        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._mod(other)
            elif (
                isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rmod__(self)
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._mod_ground(cp2)

    def __rmod__(self, other):
        ring = self.ring
        try:
            other = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other._mod(self)

    def __floordiv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(other):
            return self._floordiv(other)
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rtruediv__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._floordiv_ground(other)

    def __rfloordiv__(self, other):
        ring = self.ring
        try:
            other = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other._floordiv(self)

    def __truediv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(other):
            return self._truediv(other)
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rtruediv__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._floordiv_ground(other)

    def __rtruediv__(self, other):
        ring = self.ring
        try:
            other = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other._truediv(self)

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

    def _check(self):
        """Validate polynomial structure."""
        assert isinstance(self, PolyElement)
        assert isinstance(self.ring, PolyRing)

        dom = self.ring.domain
        assert isinstance(dom, Domain)

        for monom, coeff in self.iterterms():
            assert dom.of_type(coeff)
            assert len(monom) == self.ring.ngens
            assert all(isinstance(exp, int) and exp >= 0 for exp in monom)

    def new(self, init) -> PolyElement[Er]:
        """Create a new polynomial element in the same ring."""
        return self.__class__(self.ring, init)

    def parent(self) -> PolynomialRing[Er]:
        """Return the parent domain of this polynomial."""
        return self.ring.to_domain()

    def copy(self) -> PolyElement[Er]:
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

    def set_ring(self, new_ring: PolyRing[Es]) -> PolyElement[Es]:
        """Change the ring of this polynomial."""
        if self.ring == new_ring:
            return self # type: ignore
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

    def square(self) -> PolyElement[Er]:
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

    def degree(self, x=None):
        """
        The leading degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        i = self.ring.index(x)

        if not self:
            return ninf
        elif i < 0:
            return 0
        else:
            return self._degree(i)

    def degrees(self):
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not self:
            return (ninf,) * self.ring.ngens
        else:
            return self._degrees()

    def tail_degree(self, x=None):
        """
        The tail degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        i = self.ring.index(x)

        if not self:
            return ninf
        elif i < 0:
            return 0
        else:
            return min(monom[i] for monom in self.itermonoms())

    def tail_degrees(self):
        """
        A tuple containing tail degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not self:
            return (ninf,) * self.ring.ngens
        else:
            return tuple(map(min, list(zip(*self.itermonoms()))))

    def monic(self):
        """Divides all coefficients by the leading coefficient."""
        if not self:
            return self
        else:
            return self.quo_ground(self.LC)

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

    def quo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError("polynomial division")
        if not self or x == domain.one:
            return self
        return self._quo_ground(x)

    def extract_ground(self, g):
        f = self
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def quo_term(self, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError("polynomial division")
        elif not self:
            return self.ring.zero
        elif monom == self.ring.zero_monom:
            return self.quo_ground(coeff)
        return self._quo_term(term)

    def _norm(self, norm_func):
        if not self:
            return self.ring.domain.zero
        else:
            ground_abs = self.ring.domain.abs
            return norm_func([ground_abs(coeff) for coeff in self.itercoeffs()])

    def max_norm(self):
        return self._norm(max)

    def l1_norm(self):
        return self._norm(sum)

    def deflate(self, *G):
        ring = self.ring
        polys = [self] + list(G)

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

        return J, self._deflate(J, polys)

    def canonical_unit(self):
        domain = self.ring.domain
        return domain.canonical_unit(self.LC)

    def diff(self, x):
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
        ring = self.ring
        i = ring.index(x)
        return self._diff(i)

    def trunc_ground(self, p):
        if self.ring.domain.is_ZZ:
            terms = []

            for monom, coeff in self.iterterms():
                coeff = coeff % p

                if coeff > p // 2:
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [(monom, coeff % p) for monom, coeff in self.iterterms()]

        poly = self.new(terms)
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

    def compose(self, x, a=None):
        ring = self.ring
        poly = ring.zero

        if a is not None:
            replacements = [(x, a)]
        else:
            if isinstance(x, list):
                replacements = list(x)
            elif isinstance(x, dict):
                replacements = sorted(x.items(), key=lambda k: ring.index(k[0]))
            else:
                raise ValueError(
                    "expected a generator, value pair a sequence of such pairs"
                )

        replacements = [(ring.index(x), ring.ring_new(g)) for x, g in replacements]

        return self._compose(replacements, initial_poly=poly)

    def __call__(self, *values):
        if 0 < len(values) <= self.ring.ngens:
            return self.evaluate(list(zip(self.ring.gens, values)))
        else:
            raise ValueError(
                "expected at least 1 and at most %s values, got %s"
                % (self.ring.ngens, len(values))
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
        return self._equality(other)

    def _equality(self, other):
        raise NotImplementedError("_equality will be implemented in subclass")

    def __neg__(self) -> PolyElement[Er]:
        # Return (-1) * self in case of python-flint
        return self._negate()

    def _negate(self):
        raise NotImplementedError("_negate will be implemented in subclass")

    def _add(self, p2):
        raise NotImplementedError

    def _add_ground(self, cp2):
        raise NotImplementedError

    def _sub(self, p2):
        raise NotImplementedError

    def _sub_ground(self, cp2):
        raise NotImplementedError

    def _mul(self, other):
        raise NotImplementedError

    def mul_ground(self, x):
        raise NotImplementedError

    def _pow_int(self, n):
        raise NotImplementedError

    def _pow_generic(self, n):
        raise NotImplementedError

    def _pow_multinomial(self, n):
        raise NotImplementedError

    def _square(self):
        raise NotImplementedError

    def _divmod(self, other):
        raise NotImplementedError

    def _divmod_ground(self, x):
        raise NotImplementedError

    def _floordiv(self, p2):
        raise NotImplementedError

    def _floordiv_ground(self, p2):
        raise NotImplementedError

    def _truediv(self, p2):
        raise NotImplementedError

    def _term_div(self):
        raise NotImplementedError

    def rem(self, G):
        raise NotImplementedError

    def quo(self, G):
        raise NotImplementedError

    def exquo(self, G):
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

    def imul_num(self, c):
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
        raise NotImplementedError

    def _rem(self, G):
        raise NotImplementedError

    def _mod(self, other):
        raise NotImplementedError

    def _mod_ground(self, x):
        raise NotImplementedError

    @property
    def is_ground(self):
        # Return self.flint_poly.is_constant() in case of python-flint
        raise NotImplementedError

    @property
    def is_zero(self):
        # Return self.flint_poly.is_zero() in case of python-flint
        raise NotImplementedError

    @property
    def is_one(self):
        # Return self.flint_poly.is_one() in case of python-flint
        raise NotImplementedError

    @property
    def is_squarefree(self):
        raise NotImplementedError

    @property
    def is_irreducible(self):
        raise NotImplementedError

    @property
    def is_cyclotomic(self):
        raise NotImplementedError

    @property
    def LC(self):
        # Just use leafing_coefficient() in case of python-flint
        raise NotImplementedError

    @property
    def LM(self):
        # Use monomial(0) in case of python-flint
        raise NotImplementedError

    @property
    def LT(self):
        # Use monomial(0) and leafing_coefficient() in case of python-flint
        raise NotImplementedError

    def clear_denoms(self):
        """Clear denominators from polynomial coefficients."""
        raise NotImplementedError

    def _change_ring(self, new_ring):
        # Use fmpz_mpoly.compose() or fmpz_mpoly.compose() in case of python-flint
        raise NotImplementedError

    def as_expr_dict(self):
        # Can just use self.flint_poly.to_dict() in case of python-flint
        # Or this can just directly go into the baseclass as is since iterterms
        # will be implemented separately for pure python and flint versions anyways
        raise NotImplementedError

    def _cmp(self, other, op):
        # We can override this for python-flint version
        # to use the native lt, le, gt, ge methods
        raise NotImplementedError

    def to_dense(self):
        raise NotImplementedError

    def to_dict(self):
        # Return a self.flint_poly.to_dict() in case of python-flint
        raise NotImplementedError

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        # Use str(self.flint_poly).replace("^", "**") in case of python-flint
        raise NotImplementedError

    def _degree(self, i):
        raise NotImplementedError

    def _degrees(self):
        raise NotImplementedError

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
        raise NotImplementedError

    def _get_coeff(self, expv):
        raise NotImplementedError

    def const(self):
        # Use
        """Returns the constant coefficient."""
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

    def _sorted(self, seq, order):
        raise NotImplementedError

    def itercoeffs(self):
        """Iterator over coefficients of a polynomial."""
        raise NotImplementedError

    def itermonoms(self):
        """Iterator over monomials of a polynomial."""
        raise NotImplementedError

    def iterterms(self):
        """Iterator over terms of a polynomial."""
        raise NotImplementedError

    def listcoeffs(self):
        """Unordered list of polynomial coefficients."""
        raise NotImplementedError

    def listmonoms(self):
        """Unordered list of polynomial monomials."""
        raise NotImplementedError

    def listterms(self):
        """Unordered list of polynomial terms."""
        raise NotImplementedError

    def content(self):
        """Returns GCD of polynomial's coefficients."""
        # In the flint version, we will have to override
        # this to use the native content() method for ZZ
        # and use the pure python technique for other domains
        raise NotImplementedError

    def primitive(self):
        """Returns content and a primitive polynomial."""
        raise NotImplementedError

    def mul_monom(self, monom):
        raise NotImplementedError

    def mul_term(self, term):
        raise NotImplementedError

    def _quo_ground(self, x):
        raise NotImplementedError

    def _quo_term(self, term):
        raise NotImplementedError

    def _deflate(self, J, polys):
        raise NotImplementedError

    def inflate(self, J):
        raise NotImplementedError

    def gcd(self, other):
        return self.cofactors(other)[0]

    def _diff(self, i):
        # Use the native derivative() method in case of python-flint
        raise NotImplementedError

    def cofactors(self, other):
        raise NotImplementedError

    def _gcd_zero(self, other):
        raise NotImplementedError

    def _gcd_monom(self, other):
        raise NotImplementedError

    def _gcd(self, other):
        raise NotImplementedError

    def _gcd_ZZ(self, other):
        raise NotImplementedError

    def _gcd_QQ(self, g):
        raise NotImplementedError

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
        raise NotImplementedError

    def _compose(self, replacements, initial_poly):
        raise NotImplementedError

    def _div(self, fv):
        # Implement the same algorithm from [CLO] p64. in python-flint
        raise NotImplementedError

    # The following _p* and _subresultants methods can just be converted to pure python
    # methods in case of python-flint since their speeds don't exactly matter wrt the
    # flint version.
    def _prem(self, g, x):
        raise NotImplementedError

    def _pdiv(self, g, x):
        raise NotImplementedError

    def _pquo(self, g, x):
        raise NotImplementedError

    def _pexquo(self, g, x):
        raise NotImplementedError

    def _subresultants(self, g, x):
        raise NotImplementedError

    def _subs(self, subs_dict):
        raise NotImplementedError

    def _evaluate(self, eval_dict):
        raise NotImplementedError

    def _symmetrize(self):
        # For the python-flint override this can just be converted back to
        # the pure python version until python-flint provides some
        # equivalent functionality.
        raise NotImplementedError

    # TODO: following methods should point to polynomial
    # representation independent algorithm implementations.

    def half_gcdex(self, other):
        raise NotImplementedError

    def gcdex(self, other):
        raise NotImplementedError

    def resultant(self, other):
        raise NotImplementedError

    def discriminant(self):
        raise NotImplementedError

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

    def sqf_part(self):
        raise NotImplementedError

    def sqf_list(self, all=False):
        raise NotImplementedError

    def factor_list(self):
        raise NotImplementedError


class PyPolyElement(PolyElement, dict): # type: ignore
    """Python-based sparse multivariate polynomial element."""

    @classmethod
    def _new(cls, ring, init):
        obj = dict.__new__(cls)
        dict.__init__(obj, init)
        obj.ring = ring
        obj._hash = None
        return obj

    def __init__(self, ring, init):
        if not hasattr(self, 'ring'):  # Only initialize if not already done by _new
            super().__init__(init)
            self.ring = ring
            self._hash = None

    def __setitem__(self, monom, coeff):
        dict.__setitem__(self, monom, coeff)

    def __getitem__(self, monom):
        return dict.__getitem__(self, monom)

    def __len__(self):
        return dict.__len__(self)

    def _equality(self, other):
        if not other:
            return not self
        elif self.ring.is_element(other):
            return dict.__eq__(self, other)
        elif len(self) > 1:
            return False
        else:
            return self.get(self.ring.zero_monom) == other

    def _negate(self):
        return self.new([(monom, -coeff) for monom, coeff in self.iterterms()])

    def _add(self, p2):
        p = self.copy()
        get = p.get
        zero = self.ring.domain.zero
        for k, v in p2.items():
            v = get(k, zero) + v
            if v:
                p[k] = v
            else:
                del p[k]
        return p

    def _add_ground(self, cp2):
        p = self.copy()
        if not cp2:
            return p
        ring = self.ring
        zm = ring.zero_monom
        v = self.get(zm, ring.domain.zero) + cp2
        if v:
            p[zm] = v
        else:
            del p[zm]
        return p

    def _sub(self, p2):
        p = self.copy()
        get = p.get
        zero = self.ring.domain.zero
        for k, v in p2.items():
            v = get(k, zero) - v
            if v:
                p[k] = v
            else:
                del p[k]
        return p

    def _sub_ground(self, cp2):
        p = self.copy()
        #if not cp2:
            #return p
        ring = self.ring
        zm = ring.zero_monom
        v = self.get(zm, ring.domain.zero) - cp2
        if v:
            p[zm] = v
        else:
            del p[zm]
        return p

    def _mul(self, other):
        ring = self.ring
        p = ring.zero
        for exp1, v1 in self.iterterms():
            for exp2, v2 in other.iterterms():
                exp = ring.monomial_mul(exp1, exp2)
                v = v1 * v2
                p[exp] = p.get(exp, ring.domain.zero) + v
        p.strip_zero()
        return p

    def mul_ground(self, x):
        if not x:
            return self.ring.zero

        terms = [(monom, coeff * x) for monom, coeff in self.iterterms()]
        return self.new(terms)

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

    def quo(self, G):
        return self.div(G)[0]

    def exquo(self, G):
        q, r = self.div(G)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(self, G)

    def _iadd_monom(self, mc):
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

    def imul_num(self, c):
        if self in self.ring._gens_set:
            return self * c
        if not c:
            self.clear()
            return
        for exp in self:
            self[exp] *= c
        return self

    def _rem(self, G):
        ring = self.ring
        domain = ring.domain
        zero = domain.zero
        monomial_mul = ring.monomial_mul
        r = ring.zero
        term_div = self._term_div()
        ltf = self.LT
        f = self.copy()
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
        return not self or (len(self) == 1 and self.ring.zero_monom in self)

    @property
    def is_zero(self):
        return not self

    @property
    def is_one(self):
        return self == self.ring.one

    @property
    def is_squarefree(self):
        if not self.ring.ngens:
            return True
        return self.ring.dmp_sqf_p(self)

    @property
    def is_irreducible(self):
        if not self.ring.ngens:
            return True
        return self.ring.dmp_irreducible_p(self)

    @property
    def is_cyclotomic(self):
        if self.ring.is_univariate:
            return self.ring.dup_cyclotomic_p(self)
        else:
            raise MultivariatePolynomialError("cyclotomic polynomial")

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

    @property
    def LT(self):
        expv = self.leading_expv()
        if expv is None:
            return (self.ring.zero_monom, self.ring.domain.zero)
        else:
            return (expv, self._get_coeff(expv))

    def clear_denoms(self) -> tuple[Er, PolyElement[Er]]:
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
        if self.ring.symbols != new_ring.symbols:
            terms = list(zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols)))
            return new_ring.from_terms(terms, self.ring.domain)
        else:
            return new_ring.from_dict(self, self.ring.domain)

    def as_expr_dict(self) -> dict[tuple[int, ...], Expr]:
        # Can just use self.flint_poly.to_dict() in case of python-flint
        # Or this can just directly go into the baseclass as is since iterterms
        # will be implemented separately for pure python and flint versions anyways
        to_sympy = self.ring.domain.to_sympy
        return {monom: to_sympy(coeff) for monom, coeff in self.iterterms()}

    def _cmp(self, other, op):
        if self.ring.is_element(other):
            return op(self.sort_key(), other.sort_key())
        else:
            return NotImplemented

    def to_dense(self):
        return dmp_from_dict(self, self.ring.ngens - 1, self.ring.domain)

    def to_dict(self):
        return dict(self)

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
        return int(max(monom[i] for monom in self.itermonoms()))

    def _degrees(self):
        return tuple(map(max, list(zip(*self.itermonoms()))))

    def leading_expv(self):
        if self:
            return self.ring.leading_expv(self)
        else:
            return None

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def const(self):
        return self._get_coeff(self.ring.zero_monom)

    def coeff(self, element):
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
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self.ring.domain.one
        return p

    def leading_term(self):
        p = self.ring.zero
        expv = self.leading_expv()
        if expv is not None:
            p[expv] = self[expv]
        return p

    def coeffs(self, order=None):
        return [coeff for _, coeff in self.terms(order)]

    def monoms(self, order=None):
        return [monom for monom, _ in self.terms(order)]

    def terms(self, order=None):
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
        return iter(self.values())

    def itermonoms(self):
        return iter(self.keys())

    def iterterms(self):
        return iter(self.items())

    def listcoeffs(self):
        return list(self.values())

    def listmonoms(self):
        return list(self.keys())

    def listterms(self):
        return list(self.items())

    def content(self):
        domain = self.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in self.itercoeffs():
            cont = gcd(cont, coeff)

        return cont

    def primitive(self):
        cont = self.content()
        if cont == self.ring.domain.zero:
            return (cont, self)
        return cont, self.quo_ground(cont)

    def mul_monom(self, monom):
        monomial_mul = self.ring.monomial_mul
        terms = [
            (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in self.items()
        ]
        return self.new(terms)

    def mul_term(self, term):
        monom, coeff = term

        if not self or not coeff:
            return self.ring.zero
        elif monom == self.ring.zero_monom:
            return self.mul_ground(coeff)

        monomial_mul = self.ring.monomial_mul
        terms = [
            (monomial_mul(f_monom, monom), f_coeff * coeff)
            for f_monom, f_coeff in self.items()
        ]
        return self.new(terms)

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

    def inflate(self, J):
        poly = self.ring.zero

        for I, coeff in self.iterterms():
            N = [i * j for i, j in zip(I, J)]
            poly[tuple(N)] = coeff

        return poly

    def gcd(self, other):
        return self.cofactors(other)[0]

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

    def cofactors(self: PolyElement[Er], other: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        if not self and not other:
            zero = self.ring.zero
            return zero, zero, zero
        elif not self:
            h, cff, cfg = self._gcd_zero(other)
            return h, cff, cfg
        elif not other:
            h, cfg, cff = other._gcd_zero(self)
            return h, cff, cfg
        elif len(self) == 1:
            h, cff, cfg = self._gcd_monom(other)
            return h, cff, cfg
        elif len(other) == 1:
            h, cfg, cff = other._gcd_monom(self)
            return h, cff, cfg

        J, (self, other) = self.deflate(other)
        h, cff, cfg = self._gcd(other)

        return (h.inflate(J), cff.inflate(J), cfg.inflate(J))

    def _gcd_zero(self, other):
        one, zero = self.ring.one, self.ring.zero
        if other.is_nonnegative:
            return other, zero, one
        else:
            return -other, zero, -one

    def _gcd_monom(self, other: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        ring = self.ring
        ground_gcd = ring.domain.gcd
        ground_quo = ring.domain.quo
        monomial_gcd = ring.monomial_gcd
        monomial_ldiv = ring.monomial_ldiv
        mf, cf = self.listterms()[0]
        _mgcd, _cgcd = mf, cf
        for mg, cg in other.iterterms():
            _mgcd = monomial_gcd(_mgcd, mg)
            _cgcd = ground_gcd(_cgcd, cg)
        h = self.new([(_mgcd, _cgcd)])
        cff = self.new([(monomial_ldiv(mf, _mgcd), ground_quo(cf, _cgcd))])
        cfg = self.new(
            [
                (monomial_ldiv(mg, _mgcd), ground_quo(cg, _cgcd))
                for mg, cg in other.iterterms()
            ]
        )
        return h, cff, cfg

    def _gcd(self, other):
        ring = self.ring

        if ring.domain.is_QQ:
            return self._gcd_QQ(other)
        elif ring.domain.is_ZZ:
            return self._gcd_ZZ(other)
        else:  # TODO: don't use dense representation (port PRS algorithms)
            return ring.dmp_inner_gcd(self, other)

    def _gcd_ZZ(self, other):
        return heugcd(self, other)

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
            return qv[0], r
        else:
            return qv, r

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

    def half_gcdex(self, other):
        return self.ring.dmp_half_gcdex(self, other)

    def gcdex(self, other):
        return self.ring.dmp_gcdex(self, other)

    def resultant(self, other):
        return self.ring.dmp_resultant(self, other)

    def discriminant(self):
        return self.ring.dmp_discriminant(self)

    def decompose(self):
        if self.ring.is_univariate:
            return self.ring.dup_decompose(self)
        else:
            raise MultivariatePolynomialError("polynomial decomposition")

    def shift(self, a):
        if self.ring.is_univariate:
            return self.ring.dup_shift(self, a)
        else:
            raise MultivariatePolynomialError("shift: use shift_list instead")

    def shift_list(self, a):
        return self.ring.dmp_shift(self, a)

    def sturm(self):
        if self.ring.is_univariate:
            return self.ring.dup_sturm(self)
        else:
            raise MultivariatePolynomialError("sturm sequence")

    def gff_list(self):
        return self.ring.dmp_gff_list(self)

    def norm(self):
        return self.ring.dmp_norm(self)

    def sqf_norm(self):
        return self.ring.dmp_sqf_norm(self)

    def sqf_part(self):
        return self.ring.dmp_sqf_part(self)

    def sqf_list(self, all=False):
        return self.ring.dmp_sqf_list(self, all=all)

    def factor_list(self):
        return self.ring.dmp_factor_list(self)


class FPolyElement(PolyElement):
    """FLINT-backed sparse multivariate polynomial element."""

    @classmethod
    def _new(cls, ring, init):
        """Create a new FPolyElement efficiently."""
        obj = object.__new__(cls)
        obj.__init__(ring, init)
        return obj

    def __init__(self, ring, init):
        """Initialize FPolyElement with FLINT polynomial."""
        self.ring = ring
        self._hash = None

        # Determine FLINT polynomial class based on domain
        if ring.domain.is_ZZ:
            flint_poly_cls = flint.fmpz_mpoly
        elif ring.domain.is_QQ:
            flint_poly_cls = flint.fmpq_mpoly
        else:
            raise NotImplementedError(f"Unsupported domain for FPolyElement: {ring.domain}")

        # Create FLINT polynomial directly from init
        if isinstance(init, FPolyElement):
            init_dict = dict(init.items())  # Convert FPolyElement to monomial:coefficient dict
            self._flint_poly = flint_poly_cls(init_dict, ring.flint_ctx)
        elif isinstance(init, (PyPolyElement, dict)):
            self._flint_poly = flint_poly_cls(init, ring.flint_ctx)
        else:
            try:
                init_dict = dict(init)
                self._flint_poly = flint_poly_cls(init_dict, ring.flint_ctx)
            except (ValueError, TypeError) as e:
                raise TypeError(f"Cannot initialize FPolyElement from {type(init)}: {e}")

    @classmethod
    def from_flint_poly(cls, ring, flint_poly):
        """Create FPolyElement directly from an existing FLINT polynomial."""
        obj = object.__new__(cls)
        obj.ring = ring
        obj._hash = None
        obj._flint_poly = flint_poly
        return obj

    def __len__(self):
        return len(self._flint_poly)

    def __contains__(self, monom):
        return monom in self._flint_poly

    def __getitem__(self, monom):
        return self._flint_poly[monom]

    def __setitem__(self, monom, coeff):
        self._flint_poly[monom] = coeff
        self._hash = None

    def __delitem__(self, key):
        del self._flint_poly.to_dict()[key]
        self._hash = None

    def __iter__(self):
        return iter(self._flint_poly)

    def keys(self):
        return self._flint_poly.to_dict().keys()

    def items(self):
        return self._flint_poly.to_dict().items()

    def values(self):
        return self._flint_poly.to_dict().values()

    def get(self, key, default=None):
        try:
            return self._flint_poly[key]
        except KeyError:
            return default if default is not None else self.ring.domain.zero

    def __repr__(self):
        return str(self._flint_poly).replace("^", "**")

    def __str__(self):
        return str(self._flint_poly).replace("^", "**")

    def _get_python_ring(self):
        return PyPolyRing._new(self.ring.symbols, self.ring.domain, self.ring.order)

    def _equality(self, other):
        py_ring = self._get_python_ring()

        self_py = py_ring.from_flint(self)
        if isinstance(other, FPolyElement):
            other_python = py_ring.from_flint(other)
        else:
            other_python = other

        return self_py._equality(other_python)

    def _negate(self):
        py_ring = self._get_python_ring()
        self_py = py_ring.from_flint(self)

        result = self_py._negate()
        return self.ring.from_python(result)

    def _add(self, p2):
        python_ring = self._get_python_ring()

        self_python = python_ring.from_flint(self)

        if isinstance(p2, FPolyElement):
            p2_python = python_ring.from_flint(p2)
        else:
            p2_python = p2

        result_python = self_python._add(p2_python)

        return self.ring.from_python(result_python)

    def _add_ground(self, cp2):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._add_ground(cp2)

        return self.ring.from_python(result_py)

    def _sub(self, other):
        python_ring = self._get_python_ring()

        self_python = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other_python = python_ring.from_flint(other)
        else:
            other_python = other

        result_py = self_python._sub(other_python)

        return self.ring.from_python(result_py)

    def _sub_ground(self, cp2):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._sub_ground(cp2)

        return self.ring.from_python(result_py)

    def _mul(self, other):
        python_ring = self._get_python_ring()

        self_python = python_ring.from_flint(self)
        if isinstance(other, FPolyElement):
            other_python = python_ring.from_flint(other)
        else:
            other_python = other

        result_py = self_python._mul(other_python)
        return self.ring.from_python(result_py)

    def mul_ground(self, cp2):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.mul_ground(cp2)

        return self.ring.from_python(result)

    def _pow_int(self, n):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._pow_int(n)

        return self.ring.from_python(result_py)

    # not sure if the following exponentiation methods are needed
    # for the flint version or not

    def _pow_generic(self, n):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._pow_generic(n)

        return self.ring.from_python(result_py)

    def _pow_multinomial(self, n):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._pow_multinomial(n)

        return self.ring.from_python(result_py)

    def _square(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._square()

        return self.ring.from_python(result_py)

    def _divmod(self, other):
        """
        Bootstrapped implementation of polynomial division with remainder.
        Returns (quotient, remainder) tuple.
        """
        # Create corresponding Python ring for conversion
        python_ring = PyPolyRing._new(self.ring.symbols, self.ring.domain, self.ring.order)

        # Convert self to PyPolyElement
        self_python = python_ring.from_flint(self)

        # Convert other to PyPolyElement if it's also FPolyElement
        if isinstance(other, FPolyElement):
            other_python = python_ring.from_flint(other)
        else:
            # If other is already PyPolyElement, use it directly
            other_python = other

        # Perform the division using PyPolyElement's _divmod method
        quotient_python, remainder_python = self_python._divmod(other_python)

        # Convert results back to FPolyElement
        quotient_flint = self.ring.from_python(quotient_python)
        remainder_flint = self.ring.from_python(remainder_python)

        return quotient_flint, remainder_flint

    def _divmod_ground(self, x):
        """
        Bootstrapped implementation of division by ground domain element.
        Returns (quotient, remainder) tuple where coefficients are divided by x.
        """
        # Create corresponding Python ring for conversion
        python_ring = PyPolyRing._new(self.ring.symbols, self.ring.domain, self.ring.order)

        # Convert self to PyPolyElement
        self_python = python_ring.from_flint(self)

        # Perform the ground division using PyPolyElement's _divmod_ground method
        quotient_python, remainder_python = self_python._divmod_ground(x)

        # Convert results back to FPolyElement
        quotient_flint = self.ring.from_python(quotient_python)
        remainder_flint = self.ring.from_python(remainder_python)

        return quotient_flint, remainder_flint

    def _floordiv(self, p2):
        python_ring = self._get_python_ring()

        self_py = python_ring.from_flint(self)
        if isinstance(p2, FPolyElement):
            p2_py = python_ring.from_flint(p2)
        else:
            p2_py = p2

        result_py = self_py._floordiv(p2_py)

        return self.ring.from_python(result_py)

    def _floordiv_ground(self, cp2):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._floordiv_ground(cp2)

        return self.ring.from_python(result_py)

    def _truediv(self, p2):
        python_ring = self._get_python_ring()

        self_py = python_ring.from_flint(self)
        if isinstance(p2, FPolyElement):
            p2_py = python_ring.from_flint(p2)
        else:
            p2_py = p2

        result_py = self_py._truediv(p2_py)

        return self.ring.from_python(result_py)

    def rem(self, G):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        #print(type(G), " hahahaha")
        G_py = []

        if isinstance(G, FPolyElement):
            G_py = python_ring.from_flint(G)
        elif isinstance(G, PyPolyElement):
            G_py = G
        else:
            for poly in G:
                if isinstance(poly, FPolyElement):
                    G_py.append(python_ring.from_flint(poly))
                elif isinstance(poly, PyPolyElement):
                    G_py.append(poly)
                else:
                    raise TypeError(f"Expected PolyElement, got {type(poly)}")

        # Perform remainder computation using PyPolyElement's _rem method
        rem_py = self_py.rem(G_py)

        # Convert result back to FPolyElement
        if isinstance(rem_py, list):
            rem_flint = [self.ring.from_python(q) for q in rem_py]
        else:
            rem_flint = self.ring.from_python(rem_py)

        return rem_flint

    def quo(self, G):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        G_py = []

        if isinstance(G, FPolyElement):
            G_py = python_ring.from_flint(G)
        elif isinstance(G, PyPolyElement):
            G_py = G
        else:
            for poly in G:
                if isinstance(poly, FPolyElement):
                    G_py.append(python_ring.from_flint(poly))
                elif isinstance(poly, PyPolyElement):
                    G_py.append(poly)
                else:
                    raise TypeError(f"Expected PolyElement, got {type(poly)}")

        # Perform remainder computation using PyPolyElement's _rem method
        quo_py = self_py.quo(G_py)
        #print("QUO PYTHON", quo_py)
        # Convert result back to FPolyElement
        if isinstance(quo_py, list):
            quo_flint = [self.ring.from_python(q) for q in quo_py]
        else:
            quo_flint = self.ring.from_python(quo_py)

        return quo_flint

    def exquo(self, G):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        # print(type(G), " hahahaha")
        G_py = []

        if isinstance(G, FPolyElement):
            G_py = python_ring.from_flint(G)
        elif isinstance(G, PyPolyElement):
            G_py = G
        else:
            for poly in G:
                if isinstance(poly, FPolyElement):
                    G_py.append(python_ring.from_flint(poly))
                elif isinstance(poly, PyPolyElement):
                    G_py.append(poly)
                else:
                    raise TypeError(f"Expected PolyElement, got {type(poly)}")

        # Perform remainder computation using PyPolyElement's _rem method
        exquo_py = self_py.exquo(G_py)

        # Convert result back to FPolyElement
        if isinstance(exquo_py, list):
            exquo_flint = [self.ring.from_python(q) for q in exquo_py]
        else:
            exquo_flint = self.ring.from_python(exquo_py)

        return exquo_flint

    def _iadd_monom(self, mc):
        python_ring = self._get_python_ring()

        if self in self.ring._gens_set:
            cpself = self.copy()
        else:
            cpself = self

        cpself_py = python_ring.from_flint(cpself)

        expv, coeff = mc
        c = cpself_py.get(expv)
        if c is None:
            cpself_py[expv] = coeff
        else:
            c += coeff
            if c:
                cpself_py[expv] = c
            else:
                del cpself_py[expv]

        cpself._flint_poly = self.ring.from_python(cpself_py)._flint_poly
        return cpself

    def imul_num(self, c):
        if self in self.ring._gens_set:
            return self * c

        c = self.ring.domain.convert(c)

        if not c:
            python_ring = self._get_python_ring()
            self_python = python_ring.from_flint(self)
            self_python.clear()
            result = self.ring.from_python(self_python)
            self._flint_poly = result._flint_poly
            self._hash = None
            return

        python_ring = PyPolyRing._new(self.ring.symbols, self.ring.domain, self.ring.order)
        self_python = python_ring.from_flint(self)

        for exp in list(self_python.keys()):
            self_python[exp] *= c

        result = self.ring.from_python(self_python)

        self._flint_poly = result._flint_poly
        self._hash = None

        return self

    def _mod(self, other):
        python_ring = self._get_python_ring()

        self_py = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other = python_ring.from_flint(other)
        else:
            other = other

        result_py = self_py._mod(other)

        return self.ring.from_python(result_py)

    def _mod_ground(self, x):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._mod_ground(x)

        return self.ring.from_python(result_py)

    @property
    def is_ground(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_ground

    @property
    def is_zero(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_zero

    @property
    def is_one(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_one

    @property
    def is_squarefree(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_squarefree

    @property
    def is_irreducible(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_irreducible

    @property
    def is_cyclotomic(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.is_cyclotomic

    @property
    def LC(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.LC

    @property
    def LM(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.LM

    @property
    def LT(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.LT

    def clear_denoms(self):
        # XXXXXXX
        python_ring = self._get_python_ring()

        self_py = python_ring.from_flint(self)

        content, poly_py = self_py.clear_denoms()

        if isinstance(poly_py, PyPolyElement):
            if _supports_flint(poly_py.ring.domain, poly_py.ring.order):
                poly_flint = self.ring.from_python(poly_py)
                return content, poly_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % poly_py.ring.domain, poly_py.ring.order)
        elif isinstance(poly_py, FPolyElement):
            return content, poly_py

    def _change_ring(self, new_ring):
        if self.ring.symbols != new_ring.symbols:
            terms = list(zip(*_dict_reorder(self._flint_poly.to_dict(), self.ring.symbols, new_ring.symbols)))
            return new_ring.from_terms(terms, self.ring.domain)
        else:
            return new_ring.from_dict(self._flint_poly.to_dict(), self.ring.domain)

    def as_expr_dict(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.as_expr_dict()

    def _cmp(self, other, op):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other = python_ring.from_flint(other)
        else:
            other = other

        return self_py._cmp(other, op)

    def to_dense(self):
        return dmp_from_dict(self._flint_poly.to_dict(), self.ring.ngens - 1, self.ring.domain)

    def to_dict(self):
        return self._flint_poly.to_dict()

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.str(printer, precedence, exp_pattern, mul_symbol)

    def _degree(self, i):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return int(self_py._degree(i))

    def _degrees(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        degrees = self_py._degrees()
        for deg in degrees:
            deg = int(deg)

        return degrees

    def leading_expv(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        return self_py.leading_expv()

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def const(self):
        return self._get_coeff(self.ring.zero_monom)

    def coeff(self, element):
        python_ring = self._get_python_ring()

        self_py = python_ring.from_flint(self)

        if isinstance(element, FPolyElement):
            element_py = python_ring.from_flint(element)
        else:
            element_py = element

        return self_py.coeff(element_py)

    def leading_monom(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        monom_py = self_py.leading_monom()

        if isinstance(monom_py, FPolyElement):
            return monom_py
        else:
            return self.ring.from_python(monom_py)

    def leading_term(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        term_py = self_py.leading_term()

        if isinstance(term_py, FPolyElement):
            return term_py
        else:
            return self.ring.from_python(term_py)

    def coeffs(self, order=None):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        if order:
            assert _supports_flint(self_py.ring.domain, self_py.ring.order)
        return self_py.coeffs(order)

    def monoms(self, order=None):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        if order:
            assert _supports_flint(self_py.ring.domain, self_py.ring.order)
        return self_py.monoms(order)

    def terms(self, order=None):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)
        if order:
            assert _supports_flint(self_py.ring.domain, self_py.ring.order)
        return self_py.terms(order)

    def itercoeffs(self):
        return iter(self.values())

    def itermonoms(self):
        return iter(self.keys())

    def iterterms(self):
        return iter(self.items())

    def listcoeffs(self):
        return list(self.values())

    def listmonoms(self):
        return list(self.keys())

    def listterms(self):
        return list(self.items())

    def content(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        return self_py.content()

    def primitive(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        content, primitive_poly = self_py.primitive()

        if isinstance(primitive_poly, PyPolyElement):
            if _supports_flint(primitive_poly.ring.domain, primitive_poly.ring.order):
                primitive_flint = self.ring.from_python(primitive_poly)
                return content, primitive_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % primitive_poly.ring.domain, primitive_poly.ring.order)
        elif isinstance(primitive_poly, FPolyElement):
            return content, primitive_poly

    def mul_monom(self, monom):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.mul_monom(monom)

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain, result_py.ring.order)
        else:
            return result_py


    def mul_term(self, term):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.mul_term(term)

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain,
                                          result_py.ring.order)
        else:
            return result_py

    def _quo_ground(self, x):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._quo_ground(x)

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain, result_py.ring.order)
        else:
            return result_py

    def _quo_term(self, term):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._quo_term(term)

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain, result_py.ring.order)
        else:
            return result_py

    def _deflate(self, J, polys):
        """Deflate polynomials by dividing exponents by J, using PyPolyElement implementation."""
        # Get the Python ring for conversions
        python_ring = self._get_python_ring()

        # Convert self and all polynomials in polys to PyPolyElement
        self_py = python_ring.from_flint(self)
        polys_py = []
        for poly in polys:
            if isinstance(poly, FPolyElement):
                polys_py.append(python_ring.from_flint(poly))
            elif isinstance(poly, PyPolyElement):
                polys_py.append(poly)
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly)}")

        # Call PyPolyElement's _deflate method
        result_py = self_py._deflate(J, polys_py)

        # Convert the resulting deflated polynomials back to FPolyElement
        result_flint = []

        for poly in result_py:
            if isinstance(poly, PyPolyElement):
                assert _supports_flint(poly.ring.domain, poly.ring.order)
                result_flint.append(self.ring.from_python(poly))

            elif isinstance(poly, FPolyElement):
                result_flint.append(poly)
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly)}")

        return result_flint

    def inflate(self, J):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.inflate(J)

        return self.ring.from_python(result_py)

    def gcd(self, other):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other_py = python_ring.from_flint(other)
        else:
            other_py = other

        result_py = self_py.gcd(other_py)

        return self.ring.from_python(result_py)

    def _diff(self, i):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py._diff(i)

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain, result_py.ring.order)
        else:
            return result_py


    def cofactors(self: PolyElement[Er], other: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:

        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other_py = python_ring.from_flint(other)
        elif isinstance(other, PyPolyElement):
            assert _supports_flint(other.ring.domain, other.ring.order)
            other_py = other

        result = self_py.cofactors(other_py)

        for poly in result:
            if isinstance(poly, PyPolyElement):
                if _supports_flint(poly.ring.domain, poly.ring.order):
                    poly = self.ring.from_python(poly)
                else:
                    raise NotImplementedError(f"Unsupported domain or order by flint got {poly.ring.domain}, {poly.ring.order}")
            elif isinstance(poly, FPolyElement):
                pass
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly)}")

        return result

    def _gcd_monom(self, other: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:

        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        # Convert other to python representation
        other_py = python_ring.from_flint(other) if isinstance(other, FPolyElement) else other

        h_py, cff_py, cfg_py = self_py._gcd_monom(other_py)

        # Convert results back to flint with _supports_flint check
        if isinstance(h_py, PyPolyElement):
            if _supports_flint(h_py.ring.domain, h_py.ring.order):
                h_flint = self.ring.from_python(h_py)
            else:
                raise NotImplementedError(
                    f"Unsupported domain or order by flint: {h_py.ring.domain}, {h_py.ring.order}")
        else:
            h_flint = h_py

        if isinstance(cff_py, PyPolyElement):
            if _supports_flint(cff_py.ring.domain, cff_py.ring.order):
                cff_flint = self.ring.from_python(cff_py)
            else:
                raise NotImplementedError(
                    f"Unsupported domain or order by flint: {cff_py.ring.domain}, {cff_py.ring.order}")
        else:
            cff_flint = cff_py

        if isinstance(cfg_py, PyPolyElement):
            if _supports_flint(cfg_py.ring.domain, cfg_py.ring.order):
                cfg_flint = self.ring.from_python(cfg_py)
            else:
                raise NotImplementedError(
                    f"Unsupported domain or order by flint: {cfg_py.ring.domain}, {cfg_py.ring.order}")
        else:
            cfg_flint = cfg_py

        return (h_flint, cff_flint, cfg_flint)

    def cancel(self, g):
        # Bootstrap remaining
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
        """Compose polynomial by substituting generators with given polynomials."""
        python_ring = self._get_python_ring()

        # Convert self to PyPolyElement
        self_py = python_ring.from_flint(self)

        # Convert initial_poly to PyPolyElement
        if isinstance(initial_poly, FPolyElement):
            initial_poly_py = python_ring.from_flint(initial_poly)
        elif isinstance(initial_poly, PyPolyElement):
            initial_poly_py = initial_poly
        else:
            raise TypeError(f"Expected PolyElement for initial_poly, got {type(initial_poly)}")

        # Convert replacement polynomials to PyPolyElement
        replacements_py = []
        for idx, poly in replacements:
            if isinstance(poly, FPolyElement):
                poly_py = python_ring.from_flint(poly)
            elif isinstance(poly, PyPolyElement):
                poly_py = poly
            else:
                raise TypeError(f"Expected PolyElement for replacement, got {type(poly)}")
            replacements_py.append((idx, poly_py))

        # Perform composition using PyPolyElement's _compose method
        result_py = self_py._compose(replacements_py, initial_poly_py)

        # Convert result back to FPolyElement
        return self.ring.from_python(result_py)

    def canonical_unit(self):
        domain = self.ring.domain
        return domain.canonical_unit(self.LC)

    def _div(self, fv):
        # Bootstrap remaining
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
            return qv[0], r
        else:
            return qv, r

    def _prem(self, g, x):
        # Bootstrap remaining
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
        # Bootstrap remaining
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
        # Bootstrap remaining
        f = self
        return f.pdiv(g, x)[0]

    def _pexquo(self, g, x):
        # Bootstrap remaining
        f = self
        q, r = f.pdiv(g, x)

        if r.is_zero:
            return q
        else:
            raise ExactQuotientFailed(f, g)

    def _subresultants(self, g, x):
        # Bootstrap remaining
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
        """Substitute variables with polynomials from subs_dict, using PyPolyElement implementation."""
        # Get the Python ring for conversions
        python_ring = self._get_python_ring()

        # Convert self to PyPolyElement
        self_py = python_ring.from_flint(self)

        # Convert subs_dict values to PyPolyElement if they are PolyElement

        # Call PyPolyElement's _subs method
        result_py = self_py._subs(subs_dict)

        # Convert the result back to FPolyElement
        return self.ring.from_python(result_py)

    def _evaluate(self, eval_dict):
        """Evaluate polynomial at scalar values from eval_dict, using PyPolyElement implementation."""
        # Get the Python ring for conversions
        python_ring = self._get_python_ring()

        # Convert self to PyPolyElement
        self_py = python_ring.from_flint(self)

        # Call PyPolyElement's _evaluate method
        result = self_py._evaluate(eval_dict)

        return result

    def _symmetrize(self):
        # Bootstrap remaining
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

    def half_gcdex(self, other):
        """Compute the half extended GCD using PyPolyElement implementation."""
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        # Convert or validate the other polynomial
        if isinstance(other, FPolyElement):
            other_py = python_ring.from_flint(other)
        elif isinstance(other, PyPolyElement):
            if not _supports_flint(other.ring.domain, other.ring.order):
                raise ValueError(f"Unsupported domain or order for FLINT: {other.ring.domain}, {other.ring.order}")
            other_py = other
        else:
            raise TypeError(f"Expected PolyElement, got {type(other)}")

        # Call the ring's dmp_half_gcdex method
        s_py, h_py = python_ring.dmp_half_gcdex(self_py, other_py)

        # Convert results back to FPolyElement
        if not _supports_flint(s_py.ring.domain, s_py.ring.order):
            raise ValueError(
                f"Resulting polynomial domain or order not supported by FLINT: {s_py.ring.domain}, {s_py.ring.order}")
        if not _supports_flint(h_py.ring.domain, h_py.ring.order):
            raise ValueError(
                f"Resulting polynomial domain or order not supported by FLINT: {h_py.ring.domain}, {h_py.ring.order}")

        s_flint = self.ring.from_python(s_py)
        h_flint = self.ring.from_python(h_py)

        return s_flint, h_flint

    def gcdex(self, other):
        """Compute the extended GCD using PyPolyElement implementation."""
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        # Convert or validate the other polynomial
        if isinstance(other, FPolyElement):
            other_py = python_ring.from_flint(other)
        elif isinstance(other, PyPolyElement):
            if not _supports_flint(other.ring.domain, other.ring.order):
                raise ValueError(f"Unsupported domain or order for FLINT: {other.ring.domain}, {other.ring.order}")
            other_py = other
        else:
            raise TypeError(f"Expected PolyElement, got {type(other)}")

        # Call the ring's dmp_gcdex method
        s_py, t_py, h_py = python_ring.dmp_gcdex(self_py, other_py)

        # Convert results back to FPolyElement
        if not _supports_flint(s_py.ring.domain, s_py.ring.order):
            raise ValueError(
                f"Resulting polynomial domain or order not supported by FLINT: {s_py.ring.domain}, {s_py.ring.order}")
        if not _supports_flint(t_py.ring.domain, t_py.ring.order):
            raise ValueError(
                f"Resulting polynomial domain or order not supported by FLINT: {t_py.ring.domain}, {t_py.ring.order}")
        if not _supports_flint(h_py.ring.domain, h_py.ring.order):
            raise ValueError(
                f"Resulting polynomial domain or order not supported by FLINT: {h_py.ring.domain}, {h_py.ring.order}")

        s_flint = self.ring.from_python(s_py)
        t_flint = self.ring.from_python(t_py)
        h_flint = self.ring.from_python(h_py)

        return s_flint, t_flint, h_flint

    def resultant(self, other):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        if isinstance(other, FPolyElement):
            other_py = python_ring.from_flint(other)
        elif isinstance(other, PyPolyElement):
            assert _supports_flint(other.ring.domain, other.ring.order)
            other_py = other

        result = self_py.resultant(other_py)

        return result

    def discriminant(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.discriminant()

        if isinstance(result, PyPolyElement):
            if _supports_flint(result.ring.domain, result.ring.order):
                return self.ring.from_python(result)
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result.ring.domain, result.ring.order)
        elif isinstance(result, FPolyElement):
            return result

        return result

    def decompose(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.decompose()
        result_flint = []

        for poly in result:
            if isinstance(poly, FPolyElement):
                result_flint.append(poly)
            elif isinstance(poly, PyPolyElement):
                assert _supports_flint(poly.ring.domain, poly.ring.order)
                poly_flint = self.ring.from_python(poly)
                result_flint.append(poly_flint)
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly)}")

        return result_flint

    def shift(self, a):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.shift(a)

        if isinstance(result_py, FPolyElement):
            return result_py
        elif isinstance(result_py, PyPolyElement):
            assert _supports_flint(result_py.ring.domain, result_py.ring.order)
            result_flint = self.ring.from_python(result_py)
            return result_flint
        else:
            raise TypeError(f"Expected PolyElement, got {type(result_py)}")

    def shift_list(self, a):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.shift_list(a)

        if isinstance(result_py, FPolyElement):
            return result_py
        elif isinstance(result_py, PyPolyElement):
            assert _supports_flint(result_py.ring.domain, result_py.ring.order)
            result_flint = self.ring.from_python(result_py)
            return result_flint
        else:
            raise TypeError(f"Expected PolyElement, got {type(result_py)}")

    def sturm(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.sturm()
        result_flint = []

        for poly in result:
            if isinstance(poly, PyPolyElement):
                if _supports_flint(poly.ring.domain, poly.ring.order):
                    result_flint.append(self.ring.from_python(poly))
            else:
                result_flint.append(poly)

        return result_flint

    def gff_list(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.gff_list()
        result_flint = []

        for poly in result:
            if isinstance(poly, PyPolyElement):
                if _supports_flint(poly.ring.domain, poly.ring.order):
                    result_flint.append(self.ring.from_python(poly))
            else:
                result_flint.append(poly)

        return result_flint

    def norm(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result = self_py.norm()

        if isinstance(result, PyPolyElement):
            if _supports_flint(result.ring.domain, result.ring.order):
                result_flint = self.ring.from_python(result)
                return result_flint
        else:
            return result

    def sqf_part(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        result_py = self_py.sqf_part()

        if isinstance(result_py, PyPolyElement):
            if _supports_flint(result_py.ring.domain, result_py.ring.order):
                result_flint = self.ring.from_python(result_py)
                return result_flint
            else:
                raise NotImplementedError("Unsupported domain or order by flint got %s, %s" % result_py.ring.domain, result_py.ring.order)
        elif isinstance(result_py, FPolyElement):
            return result_py
        else:
            return result_py

    def sqf_list(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        cont, pypoly_list = self_py.sqf_list()
        fpoly_list = []

        for poly_tuple in pypoly_list:
            if isinstance(poly_tuple[0], PyPolyElement):
                if _supports_flint(poly_tuple[0].ring.domain, poly_tuple[0].ring.order):
                    fpoly_tuple = self.ring.from_python(poly_tuple[0]), poly_tuple[1]
                    fpoly_list.append(fpoly_tuple)
            elif isinstance(poly_tuple[1], FPolyElement):
                fpoly_list.append(poly_tuple)
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly_tuple[0])}")

        return cont, fpoly_list

    def factor_list(self):
        python_ring = self._get_python_ring()
        self_py = python_ring.from_flint(self)

        cont, pypoly_list = self_py.factor_list()
        fpoly_list = []

        for poly_tuple in pypoly_list:
            if isinstance(poly_tuple[0], PyPolyElement):
                if _supports_flint(poly_tuple[0].ring.domain, poly_tuple[0].ring.order):
                    fpoly_tuple = self.ring.from_python(poly_tuple[0]), poly_tuple[1]
                    fpoly_list.append(fpoly_tuple)
            elif isinstance(poly_tuple[1], FPolyElement):
                fpoly_list.append(poly_tuple)
            else:
                raise TypeError(f"Expected PolyElement, got {type(poly_tuple[0])}")

        return cont, fpoly_list
