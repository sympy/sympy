from sympy.core.expr import Expr
from sympy.polys.domains import Domain
from sympy.polys.domains.ring import Ring
from sympy.polys.fps_ring import PowerSeriesElement, PowerSeriesPolyRing
from sympy.polys.orderings import LexOrder
from sympy.utilities import public


@public
class PowerSeriesRing:
    """A class for representing univariate power series ring."""

    is_PowerSeriesRing: bool = True
    is_PowerSeries: bool = True
    has_assoc_Ring: bool = True
    has_assoc_Field: bool = False
    DEFAULT_PRECISION: int = 6

    domain: Domain
    symbols: str | tuple[Expr, ...]
    ring: PowerSeriesPolyRing
    gens: tuple[PowerSeriesElement, ...]
    prec: int
    order: LexOrder | str | None
    dtype: type[PowerSeriesElement]

    def __init__(
        self,
        domain: Domain | str,
        symbols: str | list[Expr] | tuple[Expr, ...],
        prec: int | None = None,
        order: LexOrder | str | None = None
        ) -> None:

        if prec is None:
            prec = self.DEFAULT_PRECISION

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

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PowerSeriesRing)
            and self.ring == other.ring
        )
