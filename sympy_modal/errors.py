"""
Exceptions for the second-order modal logic extension.
"""

from typing import Any, List

class SymPyModalError(Exception):
    """Base class for all sympy_modal exceptions."""
    pass


class FrameViolationError(SymPyModalError):
    """
    Raised when an inference is attempted that is not valid in the current Kripke frame.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)


class NecessitationError(SymPyModalError):
    """
    Raised when the necessitation rule is applied improperly, e.g., to a formula
    that depends on undischarged hypotheses.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidInferenceError(SymPyModalError):
    """
    Raised when an invalid logical inference is attempted.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)


class FormalisationError(SymPyModalError):
    """
    Raised when a natural language or code string cannot be formalised into a modal formula.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)


class NotAnAxiomError(SymPyModalError):
    """
    Raised when a formula is claimed to be an axiom but is not recognised as one
    in the current frame or logic.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)


class AmbiguousModalityError(FormalisationError):
    """
    Raised when the formalisation interface cannot uniquely determine the modality.
    """
    def __init__(self, message: str, candidates: List[str]) -> None:
        super().__init__(message)
        self.candidates = candidates


class ProofFailure:
    """
    A precise account of why proof search failed, returned by ProofContext.prove.
    Not an exception, but used to indicate a failure to find a proof.
    """
    def __init__(self, formula: Any, obstacle: str, missing_axioms: List[Any]) -> None:
        self.formula = formula
        self.obstacle = obstacle
        self.missing_axioms = missing_axioms

    def __repr__(self) -> str:
        return f"ProofFailure(formula={self.formula}, obstacle={self.obstacle!r}, missing_axioms={self.missing_axioms})"
