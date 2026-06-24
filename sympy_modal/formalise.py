"""
Layer 5: Formalisation Interface.
Parses strings representing code, leverages named modal operators to explicitly identify
the modality, and falls back to inferring modality by analyzing axioms.
"""

from typing import Optional, Dict, Any, List, Union
from dataclasses import dataclass

from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Boolean

from sympy_modal.types import PredicateVariable, ModalPredicate, Universe, BoolType, FunctionType
from sympy_modal.operators import (
    Box, Diamond, ProvabilityBox, AlethicBox, EpistemicBox,
    DeonticBox, TemporalBox, ForAllPredicates, ExistsPredicates
)
from sympy_modal.frames import KripkeFrame, Axiom
from sympy_modal.errors import FormalisationError, AmbiguousModalityError

@dataclass
class ModalSignature:
    operators: List[str]
    frame: KripkeFrame

@dataclass
class ScopeResolution:
    readings: List[Boolean]
    is_ambiguous: bool

@dataclass
class QuantifierOrder:
    is_first_order: bool
    is_second_order: bool

class FormalisationInterface:
    def __init__(self) -> None:
        # A dictionary mapping named operator classes to their modalities and default frames
        self.operator_map = {
            ProvabilityBox: ("provability", KripkeFrame.GL()),
            AlethicBox: ("alethic", KripkeFrame.S5()),
            EpistemicBox: ("epistemic", KripkeFrame.K45()),
            DeonticBox: ("deontic", KripkeFrame.D()),
            TemporalBox: ("temporal", KripkeFrame.S4()),
        }

    def _get_local_dict(self) -> Dict[str, Any]:
        """Returns the local dictionary for sympy's parse_expr."""
        d = {
            "Box": Box,
            "Diamond": Diamond,
            "ProvabilityBox": ProvabilityBox,
            "AlethicBox": AlethicBox,
            "EpistemicBox": EpistemicBox,
            "DeonticBox": DeonticBox,
            "TemporalBox": TemporalBox,
            "ForAllPredicates": ForAllPredicates,
            "ExistsPredicates": ExistsPredicates,
            "PredicateVariable": PredicateVariable,
            "FunctionType": FunctionType,
            "Universe": Universe,
            "BoolType": BoolType,
            "Symbol": Symbol
        }
        return d

    def parse_code(self, code_str: str) -> Boolean:
        """Parses a string into a SymPy modal formula."""
        from sympy.parsing.sympy_parser import parse_expr
        from sympy import Integer, Implies
        try:
            global_dict = {'Integer': Integer, 'Implies': Implies}
            expr = parse_expr(code_str, local_dict=self._get_local_dict(), global_dict=global_dict)
            return expr
        except Exception as e:
            raise FormalisationError(f"Failed to parse formula: {e}")

    def resolve_modality(self, code_str: str) -> ModalSignature:
        """
        Identifies which operators are present and which frame they require.
        First looks for explicit named operators. If none, tries to infer from axioms.
        """
        formula = self.parse_code(code_str)

        # 1. Look for explicit named operators
        found_operators = set()
        for op_class in self.operator_map.keys():
            if formula.has(op_class):
                found_operators.add(op_class)

        if len(found_operators) > 1:
            names = [op.__name__ for op in found_operators]
            raise AmbiguousModalityError(
                f"Multiple distinct named modal operators found: {names}", names
            )

        if len(found_operators) == 1:
            op_class = found_operators.pop()
            modality, frame = self.operator_map[op_class]
            return ModalSignature([modality], frame)

        # 2. Try to infer from axioms if only generic Box/Diamond are present
        if formula.has(Box) or formula.has(Diamond):
            # Check if formula contains recognizable axioms
            if KripkeFrame.GL()._is_axiom_instance(formula):
                 return ModalSignature(["provability"], KripkeFrame.GL())
            if KripkeFrame.S4()._is_axiom_instance(formula):
                 return ModalSignature(["temporal"], KripkeFrame.S4())
            if KripkeFrame.S5()._is_axiom_instance(formula):
                 return ModalSignature(["alethic"], KripkeFrame.S5())
            if KripkeFrame.D()._is_axiom_instance(formula):
                 return ModalSignature(["deontic"], KripkeFrame.D())

            # Default fallback for generic modal operators if no specific axioms match
            return ModalSignature(["minimal"], KripkeFrame.K())

        # No modal operators
        return ModalSignature(["none"], KripkeFrame.K())

    def resolve_quantifier_scope(self, code_str: str) -> ScopeResolution:
        """
        Distinguishes de re from de dicto readings.
        For SymPy code strings, the scope is explicit in the syntax. We return the single explicit reading.
        """
        formula = self.parse_code(code_str)

        # In a fully natural language interface, this would generate multiple readings.
        # Since we parse strict SymPy code, the scope is unambiguously defined by the AST.
        # So we just return the explicit reading.
        return ScopeResolution(readings=[formula], is_ambiguous=False)

    def resolve_order(self, code_str: str) -> QuantifierOrder:
        """
        Identifies whether quantification is first-order or second-order.
        """
        formula = self.parse_code(code_str)

        is_second_order = formula.has(ForAllPredicates) or formula.has(ExistsPredicates)

        # Check for first-order quantification
        # SymPy doesn't actually have first-order quantifiers natively in boolalg, they are often in sympy.tensor or custom logical modules.
        # But for modal logic, we'll assume any free variables could be implicitly universally quantified if we had a first-order logic.
        # For simplicity, if it's not second order, we'll mark it first-order iff it has basic unbound Symbols not used in second-order.
        is_first_order = False
        if not is_second_order:
            if any(isinstance(a, Symbol) for a in formula.free_symbols):
                is_first_order = True

        return QuantifierOrder(is_first_order=is_first_order, is_second_order=is_second_order)

    def infer_frame(self, modal_signature: ModalSignature) -> KripkeFrame:
        """
        Selects the most conservative frame consistent with the detected modality.
        """
        return modal_signature.frame

    def formalise(self, code_str: str, frame: Optional[KripkeFrame] = None) -> Union[Boolean, FormalisationError]:
        """
        Composes the above into a complete formalisation.
        """
        try:
            formula = self.parse_code(code_str)
            if frame is None:
                sig = self.resolve_modality(code_str)
                frame = self.infer_frame(sig)

            # Additional checks could be added here to ensure the formula is well-typed
            # with respect to the modal extension.
            return formula
        except Exception as e:
            if isinstance(e, FormalisationError):
                return e
            return FormalisationError(f"Formalisation failed: {e}")
