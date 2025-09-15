"""Sparse polynomial rings."""

from __future__ import annotations

from typing import (
    Generic,
    overload,
    Callable,
    Iterable,
    Iterator,
    TYPE_CHECKING,
    Mapping,
    cast,
    Sequence,
)

from operator import add, mul, lt, le, gt, ge
from functools import reduce
from types import GeneratorType

from sympy.external.gmpy import MPQ
from sympy.core.cache import cacheit
from sympy.core.expr import Expr
from sympy.core.intfunc import igcd
from sympy.core.symbol import Symbol, symbols as _symbols
from sympy.core.sympify import CantSympify, sympify
from sympy.ntheory.multinomial import multinomial_coefficients
from sympy.polys.compatibility import IPolys
from sympy.polys.constructor import construct_domain
from sympy.polys.densebasic import dup, dmp, dmp_to_dict, dup_from_dict, dmp_from_dict
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domain import Domain, Er, Es, Et
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


if TYPE_CHECKING:
    from typing import TypeIs
    from sympy.polys.fields import FracField
    from types import NotImplementedType

    _str = str


Mon = tuple[int, ...]


ninf = float("-inf")


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


def _parse_symbols(
    symbols: str | Expr | Sequence[str] | Sequence[Expr],
) -> Sequence[Expr]:
    if isinstance(symbols, str):
        return _symbols(symbols, seq=True) if symbols else ()
    elif isinstance(symbols, Expr):
        return (symbols,)
    elif is_sequence(symbols):
        if all(isinstance(s, str) for s in symbols):
            return _symbols(symbols)  # type: ignore
        elif all(isinstance(s, Expr) for s in symbols):
            return cast("Sequence[Expr]", symbols)

    raise GeneratorsError(
        "expected a string, Symbol or expression or a non-empty sequence of strings, Symbols or expressions"
    )


class PolyRing(DefaultPrinting, IPolys[Er], Generic[Er]):
    """Multivariate distributed polynomial ring."""

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
    monomial_mulpow: Callable[[Mon, Mon, int], Mon]
    monomial_ldiv: Callable[[Mon, Mon], Mon]
    monomial_div: Callable[[Mon, Mon], Mon]
    monomial_lcm: Callable[[Mon, Mon], Mon]
    monomial_gcd: Callable[[Mon, Mon], Mon]
    leading_expv: Callable[[PolyElement[Er]], Mon]
    zero_monom: Mon

    @overload
    def __new__(
        cls,
        symbols: str | Expr | Sequence[str] | Sequence[Expr],
        domain: Domain[Er],
        order: str | MonomialOrder | None = lex,
    ) -> PolyRing[Er]: ...

    @overload
    def __new__(
        cls,
        symbols: str | Expr | Sequence[str] | Sequence[Expr],
        domain: PolyRing[Es],
        order: str | MonomialOrder | None = lex,
    ) -> PolyRing[PolyElement[Es]]: ...

    def __new__(  # type: ignore
        cls,
        symbols: str | Expr | Sequence[str] | Sequence[Expr],
        domain: Domain[Er] | PolyRing[Es],
        order: str | MonomialOrder | None = lex,
    ) -> PolyRing[Er] | PolyRing[PolyElement[Es]]:
        # Create a new ring instance.
        symbols = tuple(_parse_symbols(symbols))
        ngens = len(symbols)
        domain = DomainOpt.preprocess(domain)
        morder = OrderOpt.preprocess(order)

        # Validate that symbols do not overlap with domain symbols
        if isinstance(domain, CompositeDomain) and set(symbols) & set(domain.symbols):
            raise GeneratorsError(
                "polynomial ring and it's ground domain share generators"
            )

        # Create and initialize instance
        obj = object.__new__(cls)
        obj._hash_tuple = (cls.__name__, symbols, ngens, domain, order)
        obj._hash = hash(obj._hash_tuple)
        obj.symbols = symbols
        obj.ngens = ngens
        obj.domain = domain
        obj.order = morder

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

    def _init_monomial_operations(self) -> None:
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

    def _init_leading_expv_function(self, order) -> None:
        # Initialize the leading exponent vector function.
        if order is lex:
            self.leading_expv = max
        else:
            self.leading_expv = lambda f: max(f, key=order)

    def _add_generator_attributes(self) -> None:
        """Add generator attributes for Symbol names."""
        for symbol, generator in zip(self.symbols, self.gens):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if not hasattr(self, name):
                    setattr(self, name, generator)

    # Pickle support
    def __getnewargs__(self) -> tuple[tuple[Expr, ...], Domain[Er], MonomialOrder]:
        return self.symbols, self.domain, self.order

    # Hash and equality
    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, other):
        return isinstance(other, PolyRing) and self._ring_equality(other)

    def __ne__(self, other):
        return not self == other

    @overload
    def __getitem__(self, key: int) -> PolyRing[Er]: ...
    @overload
    def __getitem__(self, key: slice) -> PolyRing[Er] | Domain[Er]: ...

    def __getitem__(self, key: slice | int) -> PolyRing[Er] | Domain[Er]:
        # Get a subring with subset of symbols.
        symbols = self.symbols[key]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

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
    @overload
    def clone(
        self,
        symbols: Expr | list[Expr] | tuple[Expr, ...] | None = None,
        domain: None = None,
        order: None = None,
    ) -> PolyRing[Er]: ...

    @overload
    def clone(
        self,
        symbols: Expr | list[Expr] | tuple[Expr, ...] | None = None,
        *,
        domain: Domain[Es],
        order: None = None,
    ) -> PolyRing[Es]: ...

    @overload
    def clone(
        self,
        symbols: Expr | list[Expr] | tuple[Expr, ...] | None = None,
        *,
        domain: PolyRing[Es],
        order: None = None,
    ) -> PolyRing[PolyElement[Es]]: ...

    def clone(
        self,
        symbols: Expr | list[Expr] | tuple[Expr, ...] | None = None,
        domain: PolyRing[Es] | Domain[Es] | None = None,
        order: str | MonomialOrder | None = None,
    ) -> PolyRing[Er] | PolyRing[Es] | PolyRing[PolyElement[Es]]:
        """Create a clone with modified parameters."""
        # Convert list to tuple for hashability
        if symbols is not None and isinstance(symbols, list):
            symbols = tuple(symbols)
        return self._clone(symbols, domain, order)

    @cacheit
    def _clone(
        self,
        symbols: Expr | tuple[Expr, ...] | None,
        domain: PolyRing[Es] | Domain[Et] | None,
        order: str | MonomialOrder | None,
    ) -> PolyRing[Er] | PolyRing[Et] | PolyRing[PolyElement[Es]]:
        return PolyRing(
            symbols or self.symbols, domain or self.domain, order or self.order
        )

    def compose(self, other: PolyRing[Er]) -> PolyRing[Er]:
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

    def to_ground(self: PolyRing[PolyElement[Es]]) -> PolyRing[Es]:
        """Convert to ground domain."""
        if isinstance(self.domain, CompositeDomain) or hasattr(self.domain, "domain"):
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

    def term_new(self, monom: Mon, coeff: int | Er) -> PolyElement[Er]:
        """Create a polynomial with a single term."""
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    # Polynomial creation from various formats
    def from_dict(
        self,
        element: Mapping[Mon, int | Er | Expr] | PolyElement[Er],
        orig_domain: Domain[Er] | None = None,
    ) -> PolyElement[Er]:
        """Create polynomial from dictionary of monomials to coefficients."""
        if not isinstance(element, dict):
            raise TypeError(
                "Input must be a dictionary mapping monomials to coefficients"
            )
        return self._from_dict_ground(element, orig_domain)

    def from_terms(
        self, element: Iterable[tuple[Mon, Er]], orig_domain: Domain[Er] | None = None
    ) -> PolyElement[Er]:
        """Create polynomial from sequence of (monomial, coefficient) pairs."""
        return self.from_dict(dict(element), orig_domain)

    def from_list(self, element: dmp[Er]) -> PolyElement[Er]:
        """Create polynomial from list(dense) representation."""
        poly_dict = dmp_to_dict(element, self.ngens - 1, self.domain)
        return self.from_dict(poly_dict)

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

    def index(self, gen: PolyElement[Er] | int | str | None) -> int:
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

    def _gen_index(self, gen: int | str) -> int:
        # Get generator index from int or string.
        if isinstance(gen, int):
            if 0 <= gen < self.ngens:
                return gen
            else:
                raise ValueError(f"invalid generator index: {gen}")
        else:
            try:
                return self.symbols.index(gen)
            except ValueError:
                raise ValueError(f"invalid generator: {gen}")

    def add_gens(self, symbols: Iterable[Symbol]) -> PolyRing[Er]:
        """Add new generators to the ring."""
        syms = set(self.symbols).union(set(symbols))
        return self.clone(symbols=list(syms))

    def drop(self, *gens: PolyElement[Er] | int | str) -> PolyRing[Er] | Domain[Er]:
        """Remove specified generators from the ring."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def drop_to_ground(
        self, *gens: PolyElement[Er] | int | str | None
    ) -> PolyRing[PolyElement[Er]] | PolyRing[Er]:
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

    def symmetric_poly(self, n: int) -> PolyElement[Er]:
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

    # Internal helper methods
    def _gens(self) -> tuple[PolyElement[Er], ...]:
        # Generate the polynomial generators.
        one = self.domain.one
        generators = []

        for i in range(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            generators.append(poly)

        return tuple(generators)

    def _ring_equality(self, other: PolyRing[Er]) -> bool:
        # Check equality of two polynomial rings.
        return (self.symbols, self.domain, self.ngens, self.order) == (
            other.symbols,
            other.domain,
            other.ngens,
            other.order,
        )

    def _from_dict_ground(
        self, element: Mapping[Mon, int | Er | Expr], orig_domain=None
    ) -> PolyElement[Er]:
        # Create polynomial from dictionary with ground domain conversion.
        poly = self.zero
        domain_new = self.domain_new

        for monom, coeff in element.items():
            if coeff:  # Skip zero coefficients
                coeff = domain_new(coeff, orig_domain)
                poly[monom] = coeff

        return poly


class PolyElement(
    DomainElement, DefaultPrinting, CantSympify, dict[tuple[int, ...], Er], Generic[Er]
):
    """Element of multivariate distributed polynomial ring."""

    def __init__(
        self, ring: PolyRing[Er], init: dict[Mon, Er] | Iterable[tuple[Mon, Er]]
    ):
        super().__init__(init)
        self.ring = ring
        # This check would be too slow to run every time:
        # self._check()

    def __getnewargs__(self):
        return (self.ring, list(self.iterterms()))

    _hash = None

    def __hash__(self) -> int:  # type: ignore
        # XXX: This computes a hash of a dictionary, but currently we don't
        # protect dictionary from being changed so any use site modifications
        # will make hashing go wrong. Use this feature with caution until we
        # figure out how to make a safe API without compromising speed of this
        # low-level class.
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.ring, frozenset(self.items())))
        return _hash

    def __ne__(self, other) -> bool:
        return not self == other

    def __pos__(self) -> PolyElement[Er]:
        return self

    def __lt__(self, other) -> bool:
        return self._cmp(other, lt)

    def __le__(self, other) -> bool:
        return self._cmp(other, le)

    def __gt__(self, other) -> bool:
        return self._cmp(other, gt)

    def __ge__(self, other) -> bool:
        return self._cmp(other, ge)

    def as_expr(self, *symbols: Expr) -> Expr:
        if not symbols:
            symbols = self.ring.symbols
        elif len(symbols) != self.ring.ngens:
            raise ValueError(
                "Wrong number of symbols, expected %s got %s"
                % (self.ring.ngens, len(symbols))
            )
        return expr_from_dict(self.as_expr_dict(), *symbols)

    @overload
    def __add__(self, other: PolyElement[Er] | Er | int, /) -> PolyElement[Er]: ...

    @overload
    def __add__(
        self, other: PolyElement[PolyElement[Er]], /
    ) -> PolyElement[PolyElement[Er]]: ...

    def __add__(
        self, other: PolyElement[Er] | Er | int | PolyElement[PolyElement[Er]], /
    ) -> PolyElement[Er] | PolyElement[PolyElement[Er]]:
        """Add two polynomials.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.rings import ring

        >>> _, x, y = ring('x, y', ZZ)
        >>> (x + y)**2 + (x - y)**2
        2*x**2 + 2*y**2

        """
        if self.ring.is_element(other):
            return self._add(other)

        if isinstance(other, PolyElement):
            domain = other.ring.domain
            if isinstance(domain, PolynomialRing) and domain.ring.is_element(self):
                return cast("PolyElement[PolyElement[Er]]", other)._add_ground(self)

        res = self._try_add_ground(other)
        if res is not NotImplemented:
            return res

        if isinstance(other, PolyElement):
            return other._try_add_ground(self)

        return NotImplemented

    def __radd__(self, other: Er | int) -> PolyElement[Er]:
        return self._try_add_ground(other)

    def _try_add_ground(self, other: object) -> PolyElement[Er] | NotImplementedType:
        ring = self.ring

        domain = ring.domain
        if domain.of_type(other):
            return self._add_ground(other)

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._add_ground(cp2)

    @overload
    def __sub__(self, other: PolyElement[Er] | Er | int, /) -> PolyElement[Er]: ...

    @overload
    def __sub__(
        self, other: PolyElement[PolyElement[Er]], /
    ) -> PolyElement[PolyElement[Er]]: ...

    def __sub__(
        self, other: PolyElement[Er] | Er | int | PolyElement[PolyElement[Er]], /
    ) -> PolyElement[Er] | PolyElement[PolyElement[Er]]:
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
        if self.ring.is_element(other):
            return self._sub(other)

        if isinstance(other, PolyElement):
            domain = other.ring.domain
            if isinstance(domain, PolynomialRing) and domain.ring.is_element(self):
                return cast("PolyElement[PolyElement[Er]]", other)._sub_ground(self)

        res = self._try_sub_ground(other)
        if res is not NotImplemented:
            return res

        if isinstance(other, PolyElement):
            return other._try_rsub_ground(self)

        return NotImplemented

    def __rsub__(self, other: Er | int) -> PolyElement[Er]:
        return self._try_rsub_ground(other)

    def _try_sub_ground(self, other: object) -> PolyElement[Er] | NotImplementedType:
        ring = self.ring
        domain = ring.domain

        if domain.of_type(other):
            return self._sub_ground(other)

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._sub_ground(cp2)

    def _try_rsub_ground(self, other: object) -> PolyElement[Er] | NotImplementedType:
        ring = self.ring
        domain = ring.domain

        if domain.of_type(other):
            return self._rsub_ground(other)

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._rsub_ground(cp2)

    @overload
    def __mul__(self, other: PolyElement[Er] | Er | int, /) -> PolyElement[Er]: ...

    @overload
    def __mul__(
        self, other: PolyElement[PolyElement[Er]], /
    ) -> PolyElement[PolyElement[Er]]: ...

    def __mul__(
        self, other: PolyElement[Er] | Er | int | PolyElement[PolyElement[Er]], /
    ) -> PolyElement[Er] | PolyElement[PolyElement[Er]]:
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
        if not self or not other:
            return self.ring.zero

        if self.ring.is_element(other):
            return self._mul(other)

        if isinstance(other, PolyElement):
            domain = other.ring.domain
            if isinstance(domain, PolynomialRing) and domain.ring.is_element(self):
                return cast("PolyElement[PolyElement[Er]]", other).mul_ground(self)

        res = self._try_mul_ground(other)
        if res is not NotImplemented:
            return res

        if isinstance(other, PolyElement):
            return other._try_mul_ground(self)

        return NotImplemented

    def __rmul__(self, other: Er | int) -> PolyElement[Er]:
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
        return self._try_mul_ground(other)

    def _try_mul_ground(self, other: object) -> PolyElement[Er] | NotImplementedType:
        ring = self.ring

        domain = ring.domain
        if domain.of_type(other):
            return self.mul_ground(other)

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

    def __divmod__(
        self, other: PolyElement[Er] | Er | int
    ) -> tuple[PolyElement[Er], PolyElement[Er]]:
        ring = self.ring
        if not other:
            raise ZeroDivisionError("polynomial division")
        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._divmod(other)
            elif (
                isinstance(ring.domain, PolynomialRing)
                and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rdivmod__(self)  # type: ignore
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._divmod_ground(cp2)

    def __rdivmod__(self, other: Er | int) -> tuple[PolyElement[Er], PolyElement[Er]]:
        ring = self.ring
        try:
            other_poly = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other_poly._divmod(self)

    def __mod__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring
        if not other:
            raise ZeroDivisionError("polynomial division")
        if isinstance(other, PolyElement):
            if other.ring == ring:
                return self._mod(other)
            elif (
                isinstance(ring.domain, PolynomialRing)
                and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rmod__(self)  # type: ignore
            else:
                return NotImplemented
        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._mod_ground(cp2)

    def __rmod__(self, other: Er | int) -> PolyElement[Er]:
        ring = self.ring
        try:
            other_poly = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other_poly._mod(self)

    def __floordiv__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(other):
            return self._floordiv(other)
        elif isinstance(other, PolyElement):
            if (
                isinstance(ring.domain, PolynomialRing)
                and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rtruediv__(self)  # type: ignore
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._floordiv_ground(other)

    def __rfloordiv__(self, other: Er | int) -> PolyElement[Er]:
        ring = self.ring
        try:
            other_poly = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other_poly._floordiv(self)

    def __truediv__(self, other: PolyElement[Er] | Er | int) -> PolyElement[Er]:
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif ring.is_element(other):
            return self._truediv(cast("Er", other))
        elif isinstance(other, PolyElement):
            if (
                isinstance(ring.domain, PolynomialRing)
                and ring.domain.ring == other.ring
            ):
                pass
            elif (
                isinstance(other.ring.domain, PolynomialRing)
                and other.ring.domain.ring == ring
            ):
                return other.__rtruediv__(self)  # type: ignore
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self._floordiv_ground(other)

    def __rtruediv__(self, other: Er | int) -> PolyElement[Er]:
        ring = self.ring
        try:
            other_poly = ring.ground_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other_poly._truediv(self)

    @property
    def is_generator(self) -> bool:
        return self in self.ring._gens_set

    @property
    def is_monomial(self) -> bool:
        return not self or (len(self) == 1 and self.LC == 1)

    @property
    def is_term(self) -> bool:
        return len(self) <= 1

    @property
    def is_negative(self) -> bool:
        return self.ring.domain.is_negative(self.LC)

    @property
    def is_positive(self) -> bool:
        return self.ring.domain.is_positive(self.LC)

    @property
    def is_nonnegative(self) -> bool:
        return self.ring.domain.is_nonnegative(self.LC)

    @property
    def is_nonpositive(self) -> bool:
        return self.ring.domain.is_nonpositive(self.LC)

    @property
    def is_monic(self) -> bool:
        return self.ring.domain.is_one(self.LC)

    @property
    def is_primitive(self) -> bool:
        return self.ring.domain.is_one(self.content())

    @property
    def is_linear(self) -> bool:
        return all(sum(monom) <= 1 for monom in self.itermonoms())

    @property
    def is_quadratic(self) -> bool:
        return all(sum(monom) <= 2 for monom in self.itermonoms())

    def _check(self) -> None:
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
            return cast("PolyElement[Es]", self)
        return self._change_ring(new_ring)

    def strip_zero(self) -> None:
        """Eliminate monomials with zero coefficient."""
        for monom, coeff in self.listterms():
            if not coeff:
                del self[monom]

    def almosteq(
        self, other: PolyElement[Er] | Er | int, tolerance: float | None = None
    ) -> bool:
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

    def sort_key(self) -> tuple[int, list[tuple[Mon, Er]]]:
        """Return a key for sorting polynomials."""
        return len(self), self.terms()

    def _drop(
        self, gen: PolyElement[Er] | int | str
    ) -> tuple[int, PolyRing[Er] | Domain[Er]]:
        ring = self.ring
        i = ring.index(gen)

        if ring.ngens == 1:
            return i, ring.domain
        else:
            new_ring = ring.drop(gen)
            return i, new_ring

    def _drop_multi(self, i: int) -> PolyElement[Er]:
        assert self.ring.ngens > 1
        return cast("PolyElement[Er]", self.drop(i))

    def drop(self, gen: PolyElement[Er] | int | str) -> PolyElement[Er] | Er:
        i, ring = self._drop(gen)

        if self.ring.ngens == 1:
            if self.is_ground:
                return self.coeff(1)
            else:
                raise ValueError(f"Cannot drop {gen}")
        else:
            if not isinstance(ring, PolyRing):
                raise TypeError("Ring after drop must be a PolyRing")

            poly = ring.zero

            for k, v in self.iterterms():
                if k[i] == 0:
                    K = list(k)
                    del K[i]
                    poly[tuple(K)] = v
                else:
                    raise ValueError(f"Cannot drop {gen}")

            return poly

    def drop_to_ground(
        self, gen: PolyElement[Er] | int | str | None
    ) -> PolyElement[PolyElement[Er]]:
        ring = self.ring
        if ring.ngens == 1:
            raise ValueError("Cannot drop only generator to ground")

        i = ring.index(gen)
        new_syms = list(ring.symbols)
        ground_sym = new_syms[i]
        del new_syms[i]

        new_ground = PolyRing([ground_sym], ring.domain, ring.order)
        new_ring = PolyRing(new_syms, new_ground, ring.order)

        poly = new_ring.zero
        gen = new_ground.gens[0]

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

    def degree(self, x: PolyElement[Er] | int | str | None = None) -> float:
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

    def degrees(self) -> tuple[float, ...]:
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not self:
            return (ninf,) * self.ring.ngens
        else:
            return self._degrees()

    def tail_degree(self, x: PolyElement[Er] | int | str | None = None) -> float:
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

    def tail_degrees(self) -> tuple[float, ...]:
        """
        A tuple containing tail degrees in all variables.

        Note that the degree of 0 is negative infinity (``float('-inf')``)

        """
        if not self:
            return (ninf,) * self.ring.ngens
        else:
            return tuple(map(min, list(zip(*self.itermonoms()))))

    def monic(self) -> PolyElement[Er]:
        """Divides all coefficients by the leading coefficient."""
        if not self:
            return self
        else:
            return self.quo_ground(self.LC)

    @overload
    def div(self, fv: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er]]: ...

    @overload
    def div(
        self, fv: Iterable[PolyElement[Er]]
    ) -> tuple[list[PolyElement[Er]], PolyElement[Er]]: ...

    def div(
        self, fv: PolyElement[Er] | Iterable[PolyElement[Er]]
    ) -> (
        tuple[PolyElement[Er], PolyElement[Er]]
        | tuple[list[PolyElement[Er]], PolyElement[Er]]
    ):
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

        if isinstance(fv, PolyElement):
            if fv.ring != ring:
                raise ValueError("self and f must have the same ring")
            if not fv:
                raise ZeroDivisionError("polynomial division")
            if not self:
                return (ring.zero, ring.zero)
            return self._div(fv)
        else:
            fv_list = list(fv)
            if not all(f.ring == ring for f in fv_list):
                raise ValueError("self and f must have the same ring")
            if not all(fv_list):
                raise ZeroDivisionError("polynomial division")
            if not self:
                return ([], ring.zero)
            return self._div_list(fv_list)

    def quo_ground(self, x: Er) -> PolyElement[Er]:
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError("polynomial division")
        if not self or x == domain.one:
            return self
        return self._quo_ground(x)

    def extract_ground(
        self, g: PolyElement[Er]
    ) -> tuple[Er, PolyElement[Er], PolyElement[Er]]:
        f = self
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def quo_term(self, term: tuple[Mon, Er]) -> PolyElement[Er]:
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

    def deflate(
        self, *G: PolyElement[Er]
    ) -> tuple[tuple[int, ...], list[PolyElement[Er]]]:
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

        J2 = tuple(J)

        if all(b == 1 for b in J2):
            return J2, polys

        return J2, self._deflate(J2, polys)

    def canonical_unit(self):
        domain = self.ring.domain
        return domain.canonical_unit(self.LC)

    def diff(self, x: int | str | PolyElement[Er]) -> PolyElement[Er]:
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

    def trunc_ground(self, p: Er) -> PolyElement[Er]:
        # XXX: This is not valid for all domains (e.g. GF(p))
        if self.ring.domain.is_ZZ:
            terms = []

            for monom, coeff in self.iterterms():
                coeff = coeff % p  # type: ignore

                if coeff > p // 2:  # type: ignore
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [(monom, coeff % p) for monom, coeff in self.iterterms()]  # type: ignore

        poly = self.new(terms)
        poly.strip_zero()
        return poly

    rem_ground = trunc_ground

    def lcm(self, g: PolyElement[Er]) -> PolyElement[Er]:
        f = self
        domain = f.ring.domain

        if not domain.is_Field:
            fc, f = f.primitive()
            gc, g = g.primitive()
            c = domain.lcm(fc, gc)
            h = (f * g).quo(f.gcd(g))
            return h.mul_ground(c)
        else:
            h = (f * g).quo(f.gcd(g))
            return h.monic()

    def coeff_wrt(self, x: int | str | PolyElement[Er], deg: int) -> PolyElement[Er]:
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
        monoms_list = [m[:i] + (0,) + m[i + 1 :] for m in monoms]
        return p.ring.from_dict(dict(zip(monoms_list, coeffs)))

    def compose(self, x, a=None):
        ring = self.ring
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

        return self._compose(replacements, initial_poly=poly)

    def __call__(self, *values):
        if 0 < len(values) <= self.ring.ngens:
            return self.evaluate(list(zip(self.ring.gens, values)))
        else:
            raise ValueError(
                "expected at least 1 and at most %s values, got %s"
                % (self.ring.ngens, len(values))
            )

    @overload
    def evaluate(
        self, values: list[tuple[PolyElement[Er], Er | int]]
    ) -> PolyElement[Er] | Er: ...

    @overload
    def evaluate(
        self, x: PolyElement[Er] | int | str, a: Er | int
    ) -> PolyElement[Er] | Er: ...

    def evaluate(self, *args, **kwargs) -> PolyElement[Er] | Er:
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
            return temp_result.set_ring(new_ring)  # type: ignore

    @overload
    def subs(self, values: list[tuple[Expr, Er | int]]) -> PolyElement[Er]: ...

    @overload
    def subs(self, x: PolyElement[Er] | int | str, a: Er | int) -> PolyElement[Er]: ...

    def subs(self, *args, **kwargs) -> PolyElement[Er]:
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

    def prem(
        self, g: PolyElement[Er], x: PolyElement[Er] | int | None = None
    ) -> PolyElement[Er]:
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
        x_index = self.ring.index(x)
        return self._prem(g, x_index)

    def pdiv(
        self, g: PolyElement[Er], x: PolyElement[Er] | int | None = None
    ) -> tuple[PolyElement[Er], PolyElement[Er]]:
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
        x_index = self.ring.index(x)
        return self._pdiv(g, x_index)

    def pquo(self, g: PolyElement[Er], x: PolyElement[Er] | int | None = None):
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
        x_index = self.ring.index(x)
        return self._pquo(g, x_index)

    def pexquo(self, g: PolyElement[Er], x: PolyElement[Er] | int | None = None):
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
        x_index = self.ring.index(x)
        return self._pexquo(g, x_index)

    def subresultants(self, g: PolyElement[Er], x: PolyElement[Er] | int | None = None):
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
        x_index = self.ring.index(x)
        return self._subresultants(g, x_index)

    def symmetrize(
        self,
    ) -> tuple[
        PolyElement[Er], PolyElement[Er], list[tuple[PolyElement[Er], PolyElement[Er]]]
    ]:
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

    def __eq__(self, other: object) -> bool:
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

    def __neg__(self) -> PolyElement[Er]:
        # Return (-1) * self in case of python-flint
        return self.new([(monom, -coeff) for monom, coeff in self.iterterms()])

    def _add(self, p2: PolyElement[Er]) -> PolyElement[Er]:
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

    def _add_ground(self, cp2: Er) -> PolyElement[Er]:
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

    def _sub(self, p2: PolyElement[Er]) -> PolyElement[Er]:
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

    def _sub_ground(self, cp2: Er) -> PolyElement[Er]:
        p = self.copy()
        if not cp2:
            return p
        ring = self.ring
        zm = ring.zero_monom
        v = self.get(zm, ring.domain.zero) - cp2
        if v:
            p[zm] = v
        else:
            del p[zm]
        return p

    def _rsub_ground(self, cp2: Er) -> PolyElement[Er]:
        return self.__neg__()._add_ground(cp2)

    def _mul(self, other: PolyElement[Er]) -> PolyElement[Er]:
        ring = self.ring
        p = ring.zero
        for exp1, v1 in self.iterterms():
            for exp2, v2 in other.iterterms():
                exp = ring.monomial_mul(exp1, exp2)
                v = v1 * v2
                p[exp] = p.get(exp, ring.domain.zero) + v
        p.strip_zero()
        return p

    def mul_ground(self, x: Er) -> PolyElement[Er]:
        if not x:
            return self.ring.zero

        terms = [(monom, coeff * x) for monom, coeff in self.iterterms()]
        return self.new(terms)

    def _pow_int(self, n: int) -> PolyElement[Er]:
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

    def _pow_generic(self, n: int) -> PolyElement[Er]:
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

    def _pow_multinomial(self, n: int) -> PolyElement[Er]:
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

    def _square(self) -> PolyElement[Er]:
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

    def _divmod(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er]]:
        return self.div(other)

    def _divmod_ground(self, x: Er) -> tuple[PolyElement[Er], PolyElement[Er]]:
        return self.quo_ground(x), self.rem_ground(x)

    def _floordiv(self, p2: PolyElement[Er]) -> PolyElement[Er]:
        return self.quo(p2)

    def _floordiv_ground(self, p2: Er) -> PolyElement[Er]:
        return self.quo_ground(p2)

    @overload
    def _truediv(self, p2: PolyElement[Er]) -> PolyElement[Er]: ...

    @overload
    def _truediv(self, p2: list[PolyElement[Er]]) -> list[PolyElement[Er]]: ...

    @overload
    def _truediv(self, p2: Er) -> PolyElement[Er]: ...

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

    @overload
    def rem(self, G: PolyElement[Er]) -> PolyElement[Er]: ...

    @overload
    def rem(self, G: list[PolyElement[Er]]) -> list[PolyElement[Er]]: ...

    def rem(
        self, G: PolyElement[Er] | list[PolyElement[Er]]
    ) -> PolyElement[Er] | list[PolyElement[Er]]:
        f = self
        if isinstance(G, PolyElement):
            return f._rem(G)
        else:
            if not all(G):
                raise ZeroDivisionError("polynomial division")
            return f._rem_list(G)

    @overload
    def quo(self, G: PolyElement[Er]) -> PolyElement[Er]: ...

    @overload
    def quo(self, G: list[PolyElement[Er]]) -> list[PolyElement[Er]]: ...

    def quo(self, G):
        return self.div(G)[0]

    @overload
    def exquo(self, G: list[PolyElement[Er]]) -> list[PolyElement[Er]]: ...

    @overload
    def exquo(self, G: PolyElement[Er]) -> PolyElement[Er]: ...

    def exquo(
        self, G: PolyElement[Er] | list[PolyElement[Er]]
    ) -> PolyElement[Er] | list[PolyElement[Er]]:
        q, r = self.div(G)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(self, G)

    def _iadd_monom(self, mc: tuple[Mon, Er]) -> PolyElement[Er]:
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

    def _iadd_poly_monom(
        self, p2: PolyElement[Er], mc: tuple[Mon, Er]
    ) -> PolyElement[Er]:
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

    def imul_num(self, c: Er | int) -> PolyElement[Er]:
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
        if self in self.ring._gens_set:
            return self * c
        if not c:
            self.clear()
            return self
        for exp in self:
            self[exp] *= c
        return self

    def _rem(self, g: PolyElement[Er]) -> PolyElement[Er]:
        """Remainder when dividing by a single polynomial g"""
        return self._rem_list([g])

    def _rem_list(self, G: list[PolyElement[Er]]) -> PolyElement[Er]:
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

    def _mod(self, other: PolyElement[Er]) -> PolyElement[Er]:
        return self.rem(other)

    def _mod_ground(self, x: Er) -> PolyElement[Er]:
        return self.rem_ground(x)

    @property
    def is_ground(self) -> bool:
        # Return self.flint_poly.is_constant() in case of python-flint
        return not self or (len(self) == 1 and self.ring.zero_monom in self)

    @property
    def is_zero(self) -> bool:
        # Return self.flint_poly.is_zero() in case of python-flint
        return not self

    @property
    def is_one(self) -> bool:
        # Return self.flint_poly.is_one() in case of python-flint
        return self == self.ring.one

    @property
    def is_squarefree(self) -> bool:
        if not self.ring.ngens:
            return True
        return self.ring.dmp_sqf_p(self)

    @property
    def is_irreducible(self) -> bool:
        if not self.ring.ngens:
            return True
        return self.ring.dmp_irreducible_p(self)

    @property
    def is_cyclotomic(self) -> bool:
        if self.ring.is_univariate:
            return self.ring.dup_cyclotomic_p(self)
        else:
            raise MultivariatePolynomialError("cyclotomic polynomial")

    @property
    def LC(self) -> Er:
        # Just use leafing_coefficient() in case of python-flint
        return self._get_coeff(self.leading_expv())

    @property
    def LM(self) -> Mon:
        # Use monomial(0) in case of python-flint
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom
        else:
            return expv

    @property
    def LT(self) -> tuple[Mon, Er]:
        # Use monomial(0) and leafing_coefficient() in case of python-flint
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
        # Use fmpz_mpoly.compose() or fmpz_mpoly.compose() in case of python-flint
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

    def _cmp(
        self,
        other: PolyElement[Er],
        op: Callable[
            [tuple[int, list[tuple[Mon, Er]]], tuple[int, list[tuple[Mon, Er]]]], bool
        ],
    ) -> bool:
        # We can override this for python-flint version
        # to use the native lt, le, gt, ge methods
        if self.ring.is_element(other):
            return op(self.sort_key(), other.sort_key())
        else:
            return NotImplemented

    def to_dense(self) -> dmp[Er]:
        return dmp_from_dict(self, self.ring.ngens - 1, self.ring.domain)

    def to_dup(self) -> dup[Er]:
        assert self.ring.ngens == 1
        return dup_from_dict(self, self.ring.domain)  # type: ignore

    def to_dict(self) -> dict[Mon, Er]:
        # Return a self.flint_poly.to_dict() in case of python-flint
        return dict(self)

    def str(self, printer, precedence, exp_pattern, mul_symbol) -> str:
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

    def _degree(self, i: int) -> int:
        return max(monom[i] for monom in self.itermonoms())

    def _degree_int(self, x: int | PolyElement[Er] | None = None) -> int:
        i = self.ring.index(x)

        if not self:
            return -1
        elif i < 0:
            return 0
        else:
            return self._degree(i)

    def _degrees(self) -> tuple[int, ...]:
        return tuple(map(max, list(zip(*self.itermonoms()))))

    def leading_expv(self) -> Mon | None:
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
        try:
            return self._leading_expv()
        except KeyError:
            return None

    def _leading_expv(self) -> Mon:
        if not self:
            raise KeyError
        else:
            return self.ring.leading_expv(self)

    def _get_coeff(self, expv) -> Er:
        return self.get(expv, self.ring.domain.zero)

    def const(self) -> Er:
        # Use
        """Returns the constant coefficient."""
        return self._get_coeff(self.ring.zero_monom)

    def coeff(self, element: PolyElement[Er] | int) -> Er:
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
            terms = list(cast("PolyElement[Er]", element).iterterms())
            if len(terms) == 1:
                monom, coeff = terms[0]
                if coeff == self.ring.domain.one:
                    return self._get_coeff(monom)

        raise ValueError("expected a monomial, got %s" % element)

    def leading_monom(self) -> PolyElement[Er]:
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

    def leading_term(self) -> PolyElement[Er]:
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

    def coeffs(self, order: _str | None = None) -> list[Er]:
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

    def monoms(self, order: _str | None = None) -> list[Mon]:
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

    def terms(self, order: _str | None = None) -> list[tuple[Mon, Er]]:
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

    def _sorted(
        self, seq: list[tuple[Mon, Er]], order: _str | None
    ) -> list[tuple[Mon, Er]]:
        if order is None:
            ordering = self.ring.order
        else:
            ordering = OrderOpt.preprocess(order)

        if ordering is lex:
            return sorted(seq, key=lambda monom: monom[0], reverse=True)
        else:
            return sorted(seq, key=lambda monom: ordering(monom[0]), reverse=True)

    def itercoeffs(self):
        """Iterator over coefficients of a polynomial."""
        return iter(self.values())

    def itermonoms(self):
        """Iterator over monomials of a polynomial."""
        return iter(self.keys())

    def iterterms(self) -> Iterator[tuple[Mon, Er]]:
        """Iterator over terms of a polynomial."""
        return iter(self.items())

    def listcoeffs(self) -> list[Er]:
        """Unordered list of polynomial coefficients."""
        return list(self.values())

    def listmonoms(self) -> list[Mon]:
        """Unordered list of polynomial monomials."""
        return list(self.keys())

    def listterms(self) -> list[tuple[Mon, Er]]:
        """Unordered list of polynomial terms."""
        return list(self.items())

    def content(self) -> Er:
        """Returns GCD of polynomial's coefficients."""
        # In the flint version, we will have to override
        # this to use the native content() method for ZZ
        # and use the pure python technique for other domains
        domain = self.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in self.itercoeffs():
            cont = gcd(cont, coeff)

        return cont

    def primitive(self) -> tuple[Er, PolyElement[Er]]:
        """Returns content and a primitive polynomial."""
        cont = self.content()
        if cont == self.ring.domain.zero:
            return (cont, self)
        return cont, self.quo_ground(cont)

    def mul_monom(self, monom: Mon) -> PolyElement[Er]:
        monomial_mul = self.ring.monomial_mul
        terms = [
            (monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in self.items()
        ]
        return self.new(terms)

    def mul_term(self, term: tuple[Mon, Er]) -> PolyElement[Er]:
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

    def _quo_ground(self, x: Er) -> PolyElement[Er]:
        domain = self.ring.domain
        if domain.is_Field:
            quo = domain.quo
            terms = [(monom, quo(coeff, x)) for monom, coeff in self.iterterms()]
        else:
            # XXX: This is not valid for all domains (e.g. GF(p))
            terms = [
                (monom, coeff // x)  # type: ignore
                for monom, coeff in self.iterterms()
                if not (coeff % x)  # type: ignore
            ]
        return self.new(terms)

    def _quo_term(self, term: tuple[Mon, Er]) -> PolyElement[Er]:
        term_div = self._term_div()
        terms = [term_div(t, term) for t in self.iterterms()]
        return self.new([t for t in terms if t is not None])

    def _deflate(
        self, J: tuple[int, ...], polys: list[PolyElement[Er]]
    ) -> list[PolyElement[Er]]:
        ring = self.ring
        H = []
        for p in polys:
            h = ring.zero
            for I, coeff in p.iterterms():
                N = [i // j for i, j in zip(I, J)]
                h[tuple(N)] = coeff
            H.append(h)
        return H

    def inflate(self, J: Sequence[int]) -> PolyElement[Er]:
        poly = self.ring.zero

        for I, coeff in self.iterterms():
            N = [i * j for i, j in zip(I, J)]
            poly[tuple(N)] = coeff

        return poly

    def gcd(self, other: PolyElement[Er]) -> PolyElement[Er]:
        return self.cofactors(other)[0]

    def _diff(self, i: int) -> PolyElement[Er]:
        # Use the native derivative() method in case of python-flint
        ring = self.ring
        m = ring.monomial_basis(i)
        g = ring.zero
        for expv, coeff in self.iterterms():
            if expv[i]:
                e = ring.monomial_ldiv(expv, m)
                g[e] = ring.domain_new(coeff * expv[i])
        return g

    def cofactors(
        self: PolyElement[Er], other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
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

    def _gcd_zero(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        one, zero = self.ring.one, self.ring.zero
        if other.is_nonnegative:
            return other, zero, one
        else:
            return -other, zero, -one

    def _gcd_monom(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
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

    def _gcd(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        ring = self.ring

        if ring.domain.is_QQ:
            return self._gcd_QQ(other)
        elif ring.domain.is_ZZ:
            return self._gcd_ZZ(other)
        else:  # TODO: don't use dense representation (port PRS algorithms)
            return ring.dmp_inner_gcd(self, other)

    def _gcd_ZZ(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        return heugcd(self, other)

    def _gcd_QQ(
        self, g: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
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

    def cancel(self, g: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er]]:
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

    # XXX: implement the same algorith for div from CLO
    # for python-flint
    def _div(self, fv: PolyElement[Er]) -> tuple[PolyElement[Er], PolyElement[Er]]:
        [q], r = self._div_list([fv])
        return q, r

    def _div_list(
        self, fv: list[PolyElement[Er]]
    ) -> tuple[list[PolyElement[Er]], PolyElement[Er]]:
        ring = self.ring
        s = len(fv)
        qv = [ring.zero for i in range(s)]
        p = self.copy()
        r = ring.zero
        term_div = self._term_div()
        expvs = [fx._leading_expv() for fx in fv]

        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p._leading_expv()
                term = term_div((expv, p[expv]), (expvs[i], fv[i][expvs[i]]))
                if term is not None:
                    expv1, c = term
                    qv[i] = qv[i]._iadd_monom((expv1, c))
                    p = p._iadd_poly_monom(fv[i], (expv1, -c))
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                expv = p._leading_expv()
                r = r._iadd_monom((expv, p[expv]))
                del p[expv]

        if expv == ring.zero_monom:
            r += p

        return qv, r

    # The following _p* and _subresultants methods can just be converted to pure python
    # methods in case of python-flint since their speeds don't exactly matter wrt the
    # flint version.
    def _prem(self, g: PolyElement[Er], x: int) -> PolyElement[Er]:
        f = self

        df = f._degree_int(x)
        dg = g._degree_int(x)

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

            dr = r._degree_int(x)

            if dr < dg:
                break

        c = lc_g**N

        return r * c

    def _pdiv(
        self, g: PolyElement[Er], x: int
    ) -> tuple[PolyElement[Er], PolyElement[Er]]:
        f = self

        df = f._degree_int(x)
        dg = g._degree_int(x)

        if dg < 0:
            raise ZeroDivisionError("polynomial division")

        q, r, dr = self.ring(x), f, df

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

            dr = r._degree_int(x)

            if dr < dg:
                break

        c = lc_g**N

        q = q * c
        r = r * c

        return q, r

    def _pquo(self, g: PolyElement[Er], x: int) -> PolyElement[Er]:
        f = self
        return f.pdiv(g, x)[0]

    def _pexquo(self, g: PolyElement[Er], x: int):
        f = self
        q, r = f.pdiv(g, x)

        if r.is_zero:
            return q
        else:
            raise ExactQuotientFailed(f, g)

    def _subresultants(self, g: PolyElement[Er], x: int):
        f = self

        n = f._degree_int(x)
        m = g._degree_int(x)

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

    def _subs(self, subs_dict: Mapping[int, Er]) -> PolyElement[Er]:
        ring = self.ring
        result_poly = ring.zero

        for monom, coeff in self.iterterms():
            new_coeff = coeff
            new_monom_list = list(monom)

            for i, val in subs_dict.items():
                exp = monom[i]
                if exp > 0:
                    new_coeff *= val**exp
                new_monom_list[i] = 0

            if new_coeff:
                new_monom = tuple(new_monom_list)
                if new_monom in result_poly:
                    result_poly[new_monom] += new_coeff
                    if not result_poly[new_monom]:
                        del result_poly[new_monom]
                else:
                    result_poly[new_monom] = new_coeff

        return result_poly

    def _evaluate(self, eval_dict: Mapping[int, Er]) -> Er:
        result = self.ring.domain.zero

        for monom, coeff in self.iterterms():
            monom_value = self.ring.domain.one
            for i, exp in enumerate(monom):
                if exp > 0:
                    monom_value *= eval_dict[i] ** exp

            result += coeff * monom_value

        return result

    def _symmetrize(
        self,
    ) -> tuple[
        PolyElement[Er], PolyElement[Er], list[tuple[PolyElement[Er], PolyElement[Er]]]
    ]:
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
                monom, coeff = cast("Mon", _monom), cast("Er", _coeff)
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

    def half_gcdex(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er]]:
        return self.ring.dmp_half_gcdex(self, other)

    def gcdex(
        self, other: PolyElement[Er]
    ) -> tuple[PolyElement[Er], PolyElement[Er], PolyElement[Er]]:
        return self.ring.dmp_gcdex(self, other)

    def resultant(self, other: PolyElement[Er]) -> PolyElement[Er] | Er:
        return self.ring.dmp_resultant(self, other)

    def discriminant(self) -> PolyElement[Er] | Er:
        return self.ring.dmp_discriminant(self)

    def decompose(self) -> list[PolyElement[Er]]:
        if self.ring.is_univariate:
            return self.ring.dup_decompose(self)
        else:
            raise MultivariatePolynomialError("polynomial decomposition")

    def shift(self, a: Er) -> PolyElement[Er]:
        if self.ring.is_univariate:
            return self.ring.dup_shift(self, a)
        else:
            raise MultivariatePolynomialError("shift: use shift_list instead")

    def shift_list(self, a: list[Er]) -> PolyElement[Er]:
        return self.ring.dmp_shift(self, a)

    def sturm(self) -> list[PolyElement[Er]]:
        if self.ring.is_univariate:
            return self.ring.dup_sturm(self)
        else:
            raise MultivariatePolynomialError("sturm sequence")

    def gff_list(self) -> list[tuple[PolyElement[Er], int]]:
        return self.ring.dmp_gff_list(self)

    def norm(self) -> PolyElement[MPQ]:
        # XXX: Only defined for AlgebraicField
        return self.ring.dmp_norm(self)

    def sqf_norm(self) -> tuple[list[int], PolyElement[Er], PolyElement[MPQ]]:
        # XXX: Only defined for AlgebraicField
        return self.ring.dmp_sqf_norm(self)

    def sqf_part(self) -> PolyElement[Er]:
        return self.ring.dmp_sqf_part(self)

    def sqf_list(
        self, all: bool = False
    ) -> tuple[Er, list[tuple[PolyElement[Er], int]]]:
        return self.ring.dmp_sqf_list(self, all=all)

    def factor_list(self) -> tuple[Er, list[tuple[PolyElement[Er], int]]]:
        return self.ring.dmp_factor_list(self)
