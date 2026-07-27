"""
Interface between a DPLL(X) engine and a theory solver.

Glossary
========
atom
    An indivisible proposition, EncodedCNF.encoding maps each atom to a
    natural number, so we usually work with integers

literal
    An atom together with a polarity, encoded as a nonzero integer whose
    magnitude is the atom and whose sign is the polarity.

clause
    A disjunction of literals, represented as a list of integers.

formula
    A conjunction of clauses (CNF), represented as a
    list of clauses. This is what EncodedCNF.data holds. Usually denoted as F.

assignment
    The set of literals the engine currently set to true, usually written
    M. We call M partial assignment if some of its atoms are yet to have a
    value.

satisfaction
    M satisfies a formula F when every clause of F is true in M.

T-consistency
    An assignment M is T-consistent when the conjunction of its literals is
    satisfiable in the theory T. Consistency and T-consistency are not the
    same thing!

model
    M is a model of F if M satisfies F.

explanation
    A clause stating why an assignment is T-inconsistent. It is obtained
    by negating the literals responsible for the conflict.

References
==========

.. [1] R. Nieuwenhuis, A. Oliveras, C. Tinelli, Solving SAT and SAT Modulo
       Theories: From an Abstract Davis-Putnam-Logemann-Loveland Procedure
       to DPLL(T), Journal of the ACM 53(6), 2006, pp. 937-977.

"""

from __future__ import annotations

from typing import Any, Literal, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.assumptions.cnf import EncodedCNF

Explanation = list[int]
Model = Any


class TheorySolver(Protocol):
    """
    This class implements parameter X of a DPLL(X) engine.
    Each Theory  can be instantiated with a specialized solver
    called Solver_T for a given Theory T.
    """
    @classmethod
    def from_encoded_cnf(
        cls, encoded_cnf: EncodedCNF
    ) -> tuple[TheorySolver, list[Explanation]]:
        """
        Create a Theory Solver instance from an encoded CNF.

        Parameters
        ==========

        encoded_cnf : EncodedCNF
            The formula to solve. Only the atoms belonging to the theory T
            are handled by the solver; the rest are ignored.

        Returns
        =======

        (solver, conflicts)

        conflicts : list of conflict clauses
            Conflicts that are known before anything is asserted, for
            example a pair of atoms that can never hold at the same time.
            They should be added to the formula by the caller.
        """
        raise NotImplementedError

    def assert_lit(
        self, literal: int
    ) -> tuple[Literal[False], Explanation] | None:
        """
        Notify the Theory Solver that the literal has been set true.
        This method should also have its own weaker and very cheap check().

        Parameters
        ==========
        literal: The literal that has been set to be true.

        Returns
        =======
        (False, Explanation): The literal is not satisfiable, and an explanation.

        None: No conflict was detected, or the literal does not belong to
            the theory T.

        """
        raise NotImplementedError

    def check(self) -> tuple[Literal[True], Model] | tuple[Literal[False], Explanation]:
        """
        Fully check whether if the (partial) assignment M is T-consistent.
        This method should be comprehensive, although in the future, for the optimization
        reasons we may need weaker checks.

        Returns
        =======

        (True, Model): The assignment M is T-consistent, and a model of it.

        (False, Explanation): The assignment M is T-inconsistent, and a
            conflict clause explaining which of the asserted literals
            cannot hold together.
        """
        raise NotImplementedError

    def reset(self) -> None:
        """
        Reset the (partial) state M to initial state M'.

        This operation is not needed to be comprehensive, it is enough
        for the state M' to be logically equivalent to the initial state.
        Or in other words, it can be more "simplified" version of the initial state.

        The most simplistic approach would be to rebuild the solver.
        """
dd
        raise NotImplementedError
