from __future__ import annotations

from typing import cast, Literal, TYPE_CHECKING

from sympy.core import Atom, Basic

if TYPE_CHECKING:
    from sympy.matrices.dense import Matrix


class CartanType_generator():
    """
    Constructor for actually creating things
    """
    def __call__(self, *args: str | list[str | int]):
        c = args[0]
        if isinstance(c, list):
            letter, n = c[0], int(c[1])
        elif isinstance(c, str):
            letter, n = c[0], int(c[1:])
        else:
            raise TypeError("Argument must be a string (e.g. 'A3') or a list (e.g. ['A', 3])")

        if n < 0:
            raise ValueError("Lie algebra rank cannot be negative")
        if letter == "A":
            from . import type_a
            return type_a.TypeA(n)
        if letter == "B":
            from . import type_b
            return type_b.TypeB(n)

        if letter == "C":
            from . import type_c
            return type_c.TypeC(n)

        if letter == "D":
            from . import type_d
            return type_d.TypeD(n)

        if letter == "E":
            if n >= 6 and n <= 8:
                from . import type_e
                return type_e.TypeE(cast(Literal[6, 7, 8], n))

        if letter == "F":
            if n == 4:
                from . import type_f
                return type_f.TypeF(4)

        if letter == "G":
            if n == 2:
                from . import type_g
                return type_g.TypeG(2)

CartanType = CartanType_generator()


class Standard_Cartan(Atom):
    """
    Concrete base class for Cartan types such as A4, etc
    """
    n: int
    _series: str

    def __new__(cls, series: str, n: int):
        obj = Basic.__new__(cls)
        obj.n = n
        obj._series = series
        return obj

    def rank(self) -> int:
        """
        Returns the rank of the Lie algebra
        """
        return self.n

    @property
    def series(self) -> str:
        """
        Returns the type of the Lie algebra
        """
        return self._series

    # for simple_root type_e is list[int | float], the rest are list[int]
    def simple_root(self, i: int) -> list:
        raise NotImplementedError(
            "simple_root is implemented by subclasses of Standard_Cartan"
        )

    def positive_roots(self) -> dict[int, list[int]]:
        raise NotImplementedError(
            "positive_roots is implemented by subclasses of Standard_Cartan"
        )

    def cartan_matrix(self) -> Matrix:
        raise NotImplementedError(
            "cartan_matrix is implemented by subclasses of Standard_Cartan"
        )

    def dynkin_diagram(self) -> str:
        raise NotImplementedError(
            "dynkin_diagram is implemented by subclasses of Standard_Cartan"
        )
