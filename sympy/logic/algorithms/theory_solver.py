from __future__ import annotations

from typing import Any, Literal, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.assumptions.cnf import EncodedCNF

ConflictClause = list[int]
Model = Any


class TheorySolver(Protocol):
    """
    This class implements parameter X of a DPLL(X) engine.
    Each Theory  can be instantiated with a specialized solver
    called Solvert_T for a given Theory T.
    """
    @classmethod
    def from_encoded_cnf(
        cls, encoded_cnf: EncodedCNF
    ) -> tuple[TheorySolver, list[ConflictClause]]:
        raise NotImplementedError

    def assert_lit(
        self, literal: int
    ) -> tuple[Literal[False], ConflictClause] | None:
        raise NotImplementedError

    def check(self) -> tuple[Literal[True], Model] | tuple[Literal[False], ConflictClause]:
        raise NotImplementedError

    def reset(self) -> None:
        """
        Reset the (pre) state M to initial state M'.

        This operation is not needed to be comprehensive, it is enough
        for the state M' to be logically equivalent to the initial state.
        Or in other words, it can be more "simplified" version of the initial state.
        """

        raise NotImplementedError
