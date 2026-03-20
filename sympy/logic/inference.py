"""Inference in propositional logic"""
from typing import Any, Dict, Generator, Union, Optional
from sympy.logic.boolalg import And, Not, conjuncts, to_cnf, BooleanFunction
from sympy.core.sorting import ordered
from sympy.core.sympify import sympify
from sympy.external.importtools import import_module


def literal_symbol(literal: Any) -> Any:
    """
    The symbol in this literal (without the negation).
    """
    if literal is True or literal is False:
        return literal
    try:
        if literal.is_Symbol:
            return literal
        if literal.is_Not:
            return literal.args[0]
    except AttributeError:
        pass
    raise ValueError("Argument must be a boolean literal.")


def satisfiable(
    expr: Any,
    algorithm: Optional[str] = None,
    all_models: bool = False,
    minimal: bool = False,
    use_lra_theory: bool = False,
) -> Union[Dict[Any, bool], bool, Generator]:
    """
    Check satisfiability of a propositional sentence.
    Returns a model when it succeeds.
    """
    raise NotImplementedError


def valid(expr: Any) -> bool:
    """
    Check validity of a propositional sentence.
    A valid propositional sentence is True under every assignment.
    """
    return not satisfiable(Not(expr))


def pl_true(
    expr: Any,
    model: Optional[Dict[Any, bool]] = None,
    deep: bool = False,
) -> Optional[bool]:
    """
    Returns whether the given assignment is a model or not.
    """
    return None


def entails(expr: Any, formula_set: Optional[list] = None) -> bool:
    """
    Check whether the given formula_set entails expr.
    If formula_set is empty then it returns the validity of expr.
    """
    raise NotImplementedError


class KB:
    """Base class for all knowledge bases"""

    def __init__(self, sentence: Optional[Any] = None) -> None:
        self.clauses_: set = set()
        if sentence:
            self.tell(sentence)

    def tell(self, sentence: Any) -> None:
        raise NotImplementedError

    def ask(self, query: Any) -> Any:
        raise NotImplementedError

    def retract(self, sentence: Any) -> None:
        raise NotImplementedError

    @property
    def clauses(self) -> list:
        return list(ordered(self.clauses_))


class PropKB(KB):
    """A KB for Propositional Logic. Inefficient, with no indexing."""

    def tell(self, sentence: Any) -> None:
        """Add the sentence's clauses to the KB."""
        for c in conjuncts(to_cnf(sentence)):
            self.clauses_.add(c)

    def ask(self, query: Any) -> bool:
        """Checks if the query is true given the set of clauses."""
        return entails(query, self.clauses_)

    def retract(self, sentence: Any) -> None:
        """Remove the sentence's clauses from the KB."""
        for c in conjuncts(to_cnf(sentence)):
            self.clauses_.discard(c)