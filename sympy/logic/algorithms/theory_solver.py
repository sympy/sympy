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
        raise NotImplementedError
