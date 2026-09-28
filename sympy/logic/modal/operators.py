"""
Layer 3 (Operators): Modal and higher-order logical operators.
Subclasses SymPy objects.
"""
from __future__ import annotations

from typing import Any
from sympy.logic.boolalg import Boolean, BooleanFunction
from sympy.logic.modal.types import PredicateVariable

class ModalOperator(BooleanFunction):
    """Base class for modal operators (Box, Diamond)."""

    def __new__(cls, *args: Any, **kwargs: Any) -> 'ModalOperator':
        return super().__new__(cls, *args, **kwargs)

    @classmethod
    def eval(cls, arg: Any, modality: str = "") -> None:
        pass


class Box(ModalOperator):
    """
    Necessity operator (□).
    """
    def __new__(cls, arg: Any, modality: str = "") -> 'Box':
        obj = super().__new__(cls, arg)
        obj._modality = modality
        return obj

    @property
    def modality(self) -> str:
        return self._modality

    def is_well_typed(self) -> bool:
        # Simplistic well-typedness check: if it's a Boolean it's well-typed
        return isinstance(self.args[0], Boolean)


class Diamond(ModalOperator):
    """
    Possibility operator (◇). Equivalent to ~Box(~P).
    """
    def __new__(cls, arg: Any, modality: str = "") -> 'Diamond':
        obj = super().__new__(cls, arg)
        obj._modality = modality
        return obj

    @property
    def modality(self) -> str:
        return self._modality

    def is_well_typed(self) -> bool:
        return isinstance(self.args[0], Boolean)


class ProvabilityBox(Box):
    """Necessity in Provability logic (GL)."""
    def __new__(cls, arg: Any) -> 'ProvabilityBox':
        return super().__new__(cls, arg, modality="provability")


class AlethicBox(Box):
    """Necessity in Alethic logic (S5/T)."""
    def __new__(cls, arg: Any) -> 'AlethicBox':
        return super().__new__(cls, arg, modality="alethic")


class EpistemicBox(Box):
    """Necessity in Epistemic logic (K45)."""
    def __new__(cls, arg: Any) -> 'EpistemicBox':
        return super().__new__(cls, arg, modality="epistemic")


class DeonticBox(Box):
    """Necessity in Deontic logic (D)."""
    def __new__(cls, arg: Any) -> 'DeonticBox':
        return super().__new__(cls, arg, modality="deontic")


class TemporalBox(Box):
    """Necessity in Temporal logic (S4)."""
    def __new__(cls, arg: Any) -> 'TemporalBox':
        return super().__new__(cls, arg, modality="temporal")


class ForAllPredicates(BooleanFunction):
    """
    Second-order universal quantification over a PredicateVariable.
    """
    def __new__(cls, variable: PredicateVariable, formula: Boolean) -> 'ForAllPredicates':
        if not isinstance(variable, PredicateVariable):
            raise TypeError("variable must be a PredicateVariable")
        if not isinstance(formula, Boolean):
            raise TypeError("formula must be a Boolean")
        return super().__new__(cls, variable, formula)

    @property
    def variable(self) -> PredicateVariable:
        return self.args[0]

    @property
    def formula(self) -> Boolean:
        return self.args[1]

    def is_well_typed(self) -> bool:
        # A rough check: we assume formula is well-typed if it parses into a Boolean.
        return isinstance(self.variable, PredicateVariable) and isinstance(self.formula, Boolean)


class AgentBox(Box):
    """
    Necessity (Knowledge) operator for a specific agent in Multi-Agent Epistemic Logic (e.g., K_A).
    """
    def __new__(cls, agent: Any, arg: Any) -> 'AgentBox':
        from sympy.core.sympify import sympify
        agent = sympify(agent)
        # Must include all parameters in args for SymPy tree reconstruction (subs, etc)
        obj = BooleanFunction.__new__(cls, agent, arg)
        obj._modality = f"epistemic_{agent}"
        return obj

    @property
    def agent(self) -> Any:
        # SymPy string sympification creates Symbols, we return the string representation
        # for backwards compatibility with the manual property access
        return str(self.args[0])

    # We override the formula property since Box expects it to be self.args[0]
    @property
    def formula(self) -> Any:
        return self.args[1]

    def is_well_typed(self) -> bool:
        return isinstance(self.args[1], Boolean)


class CommonKnowledge(Box):
    """
    Common Knowledge operator in Multi-Agent Epistemic Logic (C_G).
    Represents that a group of agents know P, know they know P, etc.
    """
    def __new__(cls, group: Any, arg: Any) -> 'CommonKnowledge':
        from sympy.core.sympify import sympify
        group = sympify(group)
        obj = BooleanFunction.__new__(cls, group, arg)
        obj._modality = f"common_knowledge_{group}"
        return obj

    @property
    def group(self) -> Any:
        return str(self.args[0])

    @property
    def formula(self) -> Any:
        return self.args[1]

    def is_well_typed(self) -> bool:
        return isinstance(self.args[1], Boolean)


class Next(ModalOperator):
    """
    'Next' operator (X) for Linear Temporal Logic.
    """
    def __new__(cls, arg: Any) -> 'Next':
        obj = super().__new__(cls, arg)
        obj._modality = "temporal_next"
        return obj

    @property
    def modality(self) -> str:
        return self._modality

    def is_well_typed(self) -> bool:
        return isinstance(self.args[0], Boolean)


class Until(BooleanFunction):
    """
    'Until' operator (U) for Linear Temporal Logic.
    A U B means A holds until B holds (and B must eventually hold).
    """
    def __new__(cls, arg1: Any, arg2: Any) -> 'Until':
        return super().__new__(cls, arg1, arg2)

    @property
    def left(self) -> Boolean:
        return self.args[0]

    @property
    def right(self) -> Boolean:
        return self.args[1]

    def is_well_typed(self) -> bool:
        return isinstance(self.args[0], Boolean) and isinstance(self.args[1], Boolean)


class ExistsPredicates(BooleanFunction):
    """
    Second-order existential quantification over a PredicateVariable.
    """
    def __new__(cls, variable: PredicateVariable, formula: Boolean) -> 'ExistsPredicates':
        if not isinstance(variable, PredicateVariable):
            raise TypeError("variable must be a PredicateVariable")
        if not isinstance(formula, Boolean):
            raise TypeError("formula must be a Boolean")
        return super().__new__(cls, variable, formula)

    @property
    def variable(self) -> PredicateVariable:
        return self.args[0]

    @property
    def formula(self) -> Boolean:
        return self.args[1]

    def is_well_typed(self) -> bool:
        return isinstance(self.variable, PredicateVariable) and isinstance(self.formula, Boolean)
