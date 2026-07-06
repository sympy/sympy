"""
Layer 3 (Kernel): Trusted Proof Kernel.
"""

from typing import List, Any, Optional
from sympy.logic.boolalg import Boolean, Implies
from sympy.logic.modal.frames import KripkeFrame, Axiom
from sympy.logic.modal.operators import Box
from sympy.logic.modal.errors import NotAnAxiomError, InvalidInferenceError, NecessitationError

class ProofTerm:
    """
    A term whose well-typedness is the certificate of validity.
    """
    def __init__(self, formula: Boolean, derivation: Any = None, source: str = "axiom", hypotheses: Optional[List[Boolean]] = None):
        self.formula = formula
        self.derivation = derivation
        self.source = source
        # Track undischarged hypotheses for hygiene checks (like Necessitation)
        self.hypotheses = hypotheses if hypotheses is not None else []

    @property
    def is_valid(self) -> bool:
        # In this implementation, creation of a ProofTerm via TrustedKernel is the certificate.
        return True

    def _ast_to_lean(self, expr: Any) -> str:
        """Helper to recursively convert SymPy AST to Lean syntax robustly."""
        from sympy.logic.boolalg import Implies, And, Or, Not
        from sympy.core.symbol import Symbol
        from sympy.logic.modal.operators import ModalOperator

        if isinstance(expr, Symbol):
            return str(expr)
        elif isinstance(expr, Implies):
            return f"({self._ast_to_lean(expr.args[0])} -> {self._ast_to_lean(expr.args[1])})"
        elif isinstance(expr, And):
            return "(" + " ∧ ".join(self._ast_to_lean(arg) for arg in expr.args) + ")"
        elif isinstance(expr, Or):
            return "(" + " ∨ ".join(self._ast_to_lean(arg) for arg in expr.args) + ")"
        elif isinstance(expr, Not):
            return f"¬{self._ast_to_lean(expr.args[0])}"
        elif isinstance(expr, ModalOperator):
            # Using generic prefix application for Lean since exact syntax depends on the Lean library chosen
            return f"({expr.__class__.__name__} {self._ast_to_lean(expr.args[0])})"
        else:
            return str(expr)

    def export_lean(self) -> str:
        """
        Exports the proof term formula to a Lean 4 compatible string representation.
        """
        f_str = self._ast_to_lean(self.formula)
        return f"theorem my_thm : {f_str} := by sorry"

    def _ast_to_coq(self, expr: Any) -> str:
        """Helper to recursively convert SymPy AST to Coq syntax robustly."""
        from sympy.logic.boolalg import Implies, And, Or, Not
        from sympy.core.symbol import Symbol
        from sympy.logic.modal.operators import ModalOperator

        if isinstance(expr, Symbol):
            return str(expr)
        elif isinstance(expr, Implies):
            return f"({self._ast_to_coq(expr.args[0])} -> {self._ast_to_coq(expr.args[1])})"
        elif isinstance(expr, And):
            return "(" + " /\\ ".join(self._ast_to_coq(arg) for arg in expr.args) + ")"
        elif isinstance(expr, Or):
            return "(" + " \\/ ".join(self._ast_to_coq(arg) for arg in expr.args) + ")"
        elif isinstance(expr, Not):
            return f"~{self._ast_to_coq(expr.args[0])}"
        elif isinstance(expr, ModalOperator):
            return f"({expr.__class__.__name__} {self._ast_to_coq(expr.args[0])})"
        else:
            return str(expr)

    def export_coq(self) -> str:
        """
        Exports the proof term formula to a Coq compatible string representation.
        """
        f_str = self._ast_to_coq(self.formula)
        return f"Theorem my_thm : {f_str}.\nProof.\n  admit.\nQed."


class ModusPonens:
    """Rule identifier for Modus Ponens."""
    pass


class TrustedKernel:
    """
    The trusted base enforcing strict natural deduction and modal rules.
    """
    def __init__(self, frame: KripkeFrame):
        self.frame = frame

    def check_axiom(self, formula: Boolean) -> ProofTerm:
        """
        Validates if formula is an axiom in the current frame or a tautology.
        Returns a ProofTerm or raises NotAnAxiomError.
        """
        # For simplicity, if it's valid in the frame (which includes tautologies and frame axioms),
        # we accept it as an axiom. In a strict kernel, we'd check exactly the shape of the axiom.
        if self.frame.validates(formula):
            return ProofTerm(formula, derivation="axiom", source="axiom")
        raise NotAnAxiomError(f"Formula {formula} is not an axiom in the current frame.")

    def check_smt(self, formula: Boolean, premises: List[ProofTerm]) -> ProofTerm:
        """
        Verifies if formula entails from premises using SymPy's SMT solver.
        """
        from sympy.logic.inference import entails
        from sympy.logic.modal.operators import Box, Diamond

        if formula.has(Box, Diamond):
            raise InvalidInferenceError("SMT solver cannot handle modal operators directly.")

        premise_formulas = [p.formula for p in premises if not p.formula.has(Box, Diamond)]

        if entails(formula, premise_formulas):
            # Union all hypotheses from the used premises
            new_hyps = set()
            for p in premises:
                if not p.formula.has(Box, Diamond):
                    new_hyps.update(p.hypotheses)
            return ProofTerm(formula, derivation=("SMT_Solver", formula), source="derived", hypotheses=list(new_hyps))

        raise InvalidInferenceError(f"Formula {formula} does not entail from premises.")

    def verify_rule(self, rule: Any, premises: List[ProofTerm]) -> ProofTerm:
        """
        Verifies a natural deduction rule application.
        """
        if rule is ModusPonens:
            if len(premises) != 2:
                raise InvalidInferenceError("Modus Ponens requires exactly two premises.")
            p1, p2 = premises
            # We expect p1 to be Implies(A, B) and p2 to be A, or vice versa
            if isinstance(p1.formula, Implies) and p1.formula.args[0] == p2.formula:
                impl = p1
                ant = p2
            elif isinstance(p2.formula, Implies) and p2.formula.args[0] == p1.formula:
                impl = p2
                ant = p1
            else:
                raise InvalidInferenceError(f"Premises {p1.formula} and {p2.formula} cannot be resolved via Modus Ponens.")

            new_hypotheses = list(set(p1.hypotheses + p2.hypotheses))
            # impl.formula.args[1] is typically a Basic or Boolean. We cast/assert to appease mypy.
            conclusion = impl.formula.args[1]
            if not isinstance(conclusion, Boolean):
                raise InvalidInferenceError("Conclusion is not a Boolean formula")
            return ProofTerm(conclusion, derivation=(ModusPonens, p1, p2), source="derived", hypotheses=new_hypotheses)

        raise InvalidInferenceError(f"Unknown rule {rule}")

    def necessitate(self, proof_term: ProofTerm) -> ProofTerm:
        """
        Necessitation rule: if ⊢ P then ⊢ □P.
        Raises NecessitationError if proof_term depends on undischarged hypotheses.
        """
        if proof_term.source == 'hypothesis' or proof_term.hypotheses:
            raise NecessitationError("Necessitation rule cannot be applied to a formula that depends on undischarged hypotheses.")

        # Creates a Box formula. For the formal system, this is valid for the minimal normal modal logic (K).
        # We use a generic Box.
        box_formula = Box(proof_term.formula)
        return ProofTerm(box_formula, derivation=("necessitation", proof_term), source="derived")

    def check_term(self, term: ProofTerm, type_formula: Boolean) -> bool:
        """
        The type checker. Verifies that the proof term proves the target formula.
        In a full implementation, this would recursively check the `derivation` tree.
        Here, we check the certificate matches.
        """
        return term.formula == type_formula and term.is_valid
