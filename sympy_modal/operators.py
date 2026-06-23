"""
Layer 3 (Operators): Modal and higher-order logical operators.
Subclasses SymPy objects.
"""

from typing import Any
from sympy.logic.boolalg import Boolean, BooleanFunction
from sympy_modal.types import PredicateVariable, FunctionType

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
