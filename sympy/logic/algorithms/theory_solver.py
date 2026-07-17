from __future__ import annotations

from typing import Any, Literal, Optional, Protocol, Union, TYPE_CHECKING, runtime_checkable

if TYPE_CHECKING:
    from sympy.assumptions.cnf import EncodedCNF

ConflictClause = list[int]
Model = Any


@runtime_checkable
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
        ...

    def assert_lit(
        self, literal: int
    ) -> Optional[tuple[Literal[False], ConflictClause]]:
        ...

    def check(self) -> Union[tuple[Literal[True], Model], tuple[Literal[False], ConflictClause]]:
        ...
