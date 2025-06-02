from sympy.core.expr import Expr
from sympy.polys.domains import Domain
from sympy.polys.polyoptions import (Domain as DomainOpt, Order as OrderOpt)
from sympy.polys.rings import PolyRing, PolyElement
from sympy.polys.orderings import lex, MonomialOrder
from sympy.polys.polyerrors import GeneratorsError

from typing import Any, Union, List, Tuple, Optional

def _parse_symbols(symbols):
    """Parse symbols from various input formats."""
    from sympy.core.symbol import symbols as _symbols
    from sympy.utilities.iterables import is_sequence

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


class PowerSeriesPolyRing:
    """A class for representing univariate power series ring."""

    symbols: Tuple[Expr, ...]
    domain: Domain
    ring: PolyRing
    gens: Tuple['PowerSeriesElement', ...]
    order: Union[MonomialOrder, str]
    prec: int

    def __init__(self, symbols: Union[str, List[Expr], Tuple[Expr, ...]],
                domain: Union[Domain, str],
                prec: int = 6,
                order: Union[MonomialOrder, str] = lex) -> None:
        symbols = _parse_symbols(symbols)
        domain = DomainOpt.preprocess(domain)
        processed_order = OrderOpt.preprocess(order)

        if prec < 1 or not isinstance(prec, int):
            raise ValueError("Precision must be a positive integer.")

        if domain.is_Composite and set(symbols) & set(domain.symbols):
            raise GeneratorsError("Power Series ring and it's ground domain share generators")
        elif len(symbols) != 1:
            raise GeneratorsError("Only univatiate power series rings are supported.")

        ring = PolyRing(symbols, domain)

        self.ring = ring
        self.domain = ring.domain
        self.prec = prec
        self.symbols = ring.symbols
        self.order = processed_order
        self.gens = tuple([self.from_poly(g) for g in ring.gens])

    @property
    def zero(self) -> 'PowerSeriesElement':
        """Return the zero element of the ring."""
        return self.from_poly(self.ring.zero)

    @property
    def one(self) -> 'PowerSeriesElement':
        """Return the one element of the ring."""
        return self.from_poly(self.ring.one)

    def __repr__(self) -> str:
        """Representation string for the power series polynomial ring."""
        return f"Power Series Ring in {', '.join(map(str, self.symbols))} over {self.domain} with {self.order} order"

    def __str__(self) -> str:
        """String representation of the power series polynomial ring."""
        return f"Power Series Ring in {', '.join(map(str, self.symbols))} over {self.domain} with {self.order} order"

    def __eq__(self, other: Any) -> bool:
        """Check if two power series polynomial rings are equal."""
        return (isinstance(other, PowerSeriesPolyRing) and
                self.domain == other.domain and
                self.symbols == other.symbols and
                self.prec == other.prec and
                self.order == other.order)

    def __hash__(self) -> int:
        """Hash function for the power series polynomial ring."""
        return hash((self.domain, self.symbols, self.prec, self.order))

    def from_poly(self, poly: PolyElement) -> 'PowerSeriesElement':
        """Convert a polynomial to a power series element."""
        if poly.ring.symbols != self.symbols:
            raise GeneratorsError("Polynomial symbols do not match the power series ring symbols.")
        return PowerSeriesElement(self, poly)


class PowerSeriesElement:
    """A class for representing univariate power series. """

    ring: PowerSeriesPolyRing
    poly: PolyElement
    domain: Domain
    symbols: Optional[Tuple[Expr, ...]]
    prec: int

    def __new__(cls, ring: PowerSeriesPolyRing, poly: PolyElement) -> 'PowerSeriesElement':
        obj = object.__new__(cls)
        obj.ring = ring
        obj.prec = ring.prec
        obj.poly = poly
        obj.domain = ring.domain
        obj.symbols = ring.symbols

        return obj

    def __repr__(self) -> str:
        """Representation string for the power series."""
        def format_term(coeff: Any, power: int, symbol: Expr) -> str:
            if power == 0:
                return str(abs(coeff))
            elif power == 1:
                return str(symbol) if abs(coeff) == 1 else f"{abs(coeff)}*{symbol}"
            else:
                return f"{symbol}**{power}" if abs(coeff) == 1 else f"{abs(coeff)}*{symbol}**{power}"

        symbol: Expr = self.ring.symbols[0]

        signs: List[str] = []
        expressions: List[str] = []

        for exp, coeff in self.terms():
            e: int = exp[0]
            if coeff != 0 and e < self.prec:
                if self.ring.domain.is_negative(coeff):
                    signs.append('-')
                else:
                    signs.append('+')
                expressions.append(format_term(coeff, e, symbol))

        # Building series representation
        if not expressions:
            series: str = f"O({symbol}**{self.prec})"
        else:
            if signs[0] == "+":
                result_parts = [expressions[0]]
            else:
                result_parts = ["-" + expressions[0]]

            for i in range(1, len(expressions)):
                result_parts.append(f"{signs[i]} {expressions[i]}")

            series = " ".join(result_parts)

        return f"{series} + O({symbol}**{self.prec})"

    def __eq__(self, other: Any) -> bool:
        """Equality check for two power series."""
        return (isinstance(other, PowerSeriesElement) and
                self.ring == other.ring and
                self.poly == other.poly)

    def terms(self) -> List[Tuple[Tuple[int, ...], Any]]:
        """Return an iterator over the terms of the power series."""
        return self.poly.terms()
