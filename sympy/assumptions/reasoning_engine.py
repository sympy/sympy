from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.logic.algorithms.dpll2 import SATSolver, IpasirStatus

if TYPE_CHECKING:
    from sympy.assumptions.assume import AppliedPredicate
    from sympy.assumptions.cnf import CNF, EncodedCNF


class ReasoningEngine:
    def __init__(self, factbase: EncodedCNF) -> None:
        if {0} in factbase.data:
            raise ValueError("Inconsistent assumptions")

        self._factbase = factbase
        self._solver: SATSolver | None = None # TODO: Initialize sat solver here instead of in `create_query`.

    def create_query(self, prop: CNF, _prop: CNF) -> int:
        if self._solver is not None:
            raise NotImplementedError("Multiple queries is not implemented yet.")
        guarded = self._factbase.copy()
        selector = guarded.new_auxiliary_variable()
        _encode_with_selector(prop, guarded, selector)
        _encode_with_selector(_prop, guarded, -selector)
        self._encoding = guarded.encoding
        self._solver = SATSolver(guarded.data, guarded.variables,
                                 set(), guarded.symbols)

        if self._solver.propagate() == IpasirStatus.UNSATISFIABLE:
            raise ValueError("Inconsistent assumptions")

        return selector

    def lookup(self, pred: AppliedPredicate) -> bool | None:
        assert self._solver is not None
        lit = self._encoding.get(pred)
        if lit is None:
            return None

        return self.fixed(lit)

    def fixed(self, lit: int) -> bool | None:
        assert self._solver is not None
        return {1: True, -1: False, 0: None}[self._solver.fixed(lit)]

    def ask_query(self, selector: int) -> bool | None:
        assert self._solver is not None
        # TODO: Run additional checks to see which combination of the
        # assumptions, global_assumptions, and relevant_facts are inconsistent.
        if self._solver.solve() == IpasirStatus.UNSATISFIABLE:
            raise ValueError("Inconsistent assumptions")

        # The polarity of the selector literal corresponds to whether
        # prop or _prop is true in a given model: positive implies prop
        # is true while negative implies _prop.
        #
        # Thus, if the model found by the solver sets the selector to true,
        # then prop is true in that model (and vice versa). After finding the
        # initial model, the solver looks for a model with the opposing polarity.
        # Once it does that, it knows the satisfiability of both prop and _prop.

        selector_value: int = self._solver.val(selector)
        self._solver.assume(-selector_value)
        other_value_is_satisfiable = self._solver.solve() == IpasirStatus.SATISFIABLE

        if other_value_is_satisfiable:
            return None
        else:
            return selector_value > 0


def _encode_with_selector(prop: CNF, encoded: EncodedCNF,
                          activation_literal: int) -> None:
    """
    Add each clause in prop to encoded with an additional literal: the
    negated `activation_literal`.

    Each resulting clause encodes `activation_literal => c`, where `c`
    is an original clause. When `activation_literal` is True, the original
    clauses must hold; when False, the added literal satisfies every clause,
    deactivating them.
    """
    encoded.data.extend(
        (encoded.encode(clause) - {0}) | {-activation_literal}
        for clause in prop.clauses
    )
