from sympy.core.expr import Expr
from sympy.polys.domains import Domain
from sympy.polys.domains.ring import Ring
from sympy.polys.fps_ring import PowerSeriesElement, PowerSeriesPolyRing
from sympy.polys.orderings import MonomialOrder, lex
from sympy.utilities import public

from typing import Any, Union, List, Tuple, Optional, Type


@public
class PowerSeriesRing:
    """A class for representing univariate power series ring."""

    is_PowerSeriesRing = is_PowerSeries = True
    has_assoc_Ring  = True
    has_assoc_Field = False
    DEFAULT_PRECISION = 6

    domain: Domain
    symbols: Union[str, Tuple[Expr, ...]]
    ring: PowerSeriesPolyRing
    gens: Tuple[PowerSeriesElement, ...]
    prec: int
    order: Optional[Union[MonomialOrder, str]]
    dtype: Type[PowerSeriesElement]

    def __init__(self, domain: Union[Domain, str],
                 symbols: Union[str, List[Expr], Tuple[Expr, ...]],
                 prec: Optional[int] = None,
                 order: Optional[Union[MonomialOrder, str]] = None) -> None:

        if prec is None:
            prec = self.DEFAULT_PRECISION
        if order is None:
            from sympy.polys.orderings import lex
            order = lex
        ring = PowerSeriesPolyRing(symbols, domain, prec, order)

        self.ring = ring
        self.dtype = PowerSeriesElement
        self.domain = ring.domain
        self.symbols = ring.symbols
        self.gens = ring.gens
        self.prec = ring.prec
        self.order = order

    @property
    def zero(self) -> PowerSeriesElement:
        return self.ring.zero

    @property
    def one(self) -> PowerSeriesElement:
        return self.ring.one

    def __repr__(self) -> str:
        return f'{self.domain}[[{", ".join(str(s) for s in self.symbols)}], {self.prec}]'

    def __str__(self) -> str:
        return f'{self.domain}[[{", ".join(str(s) for s in self.symbols)}], {self.prec}]'

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.ring, self.domain, self.symbols, self.prec, self.order))

    def __eq__(self, other: Any) -> bool:
        """Returns `True` if two domain are equivalent. """
        return (isinstance(other, PowerSeriesRing) and
                self.domain == other.domain and
                self.symbols == other.symbols and
                self.prec == other.prec and
                self.order == other.order)
