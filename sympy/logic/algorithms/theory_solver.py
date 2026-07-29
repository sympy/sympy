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

explanation (conflict clause)
    A clause stating why an assignment is T-inconsistent. It is obtained
    by negating the literals responsible for the conflict.

References
==========

.. [1] R. Nieuwenhuis, A. Oliveras, C. Tinelli, Solving SAT and SAT Modulo
       Theories: From an Abstract Davis-Putnam-Logemann-Loveland Procedure
       to DPLL(T), Journal of the ACM 53(6), 2006, pp. 937-977.

"""

from __future__ import annotations

from typing import Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.assumptions.cnf import EncodedCNF

Literal = int
Clause = list[Literal]
Assignment = list[Literal]
TheoryLemma = Clause
Explanation = Clause
NO_DECISION: Literal = 0


class TheorySolver(Protocol):
    """
    This class implements parameter X of a DPLL(X) engine.
    Each Theory  can be instantiated with a specialized solver
    called Solver_T for a given Theory T.
    """

    # Whether if the CDCL(X) loop is lazy or eager
    # Currently, dpll2.py is lazy i.e it calls TheorySolver
    # during/after the search. Eager is when all the Theory cruft
    # is done before any SAT search is done.
    is_lazy: bool = True

    @classmethod
    def from_encoded_cnf(
        cls, formula: EncodedCNF
    ) -> tuple[TheorySolver, list[TheoryLemma]]:
        """
        Create a Theory Solver instance from a formula. Preprocessing also
        evaluates the atoms that T decides on its own, e.g. Q.gt(2, 1) is
        True and Q.gt(1, 2) is False. Each such atom is returned as a
        lemma.

        Parameters
        ==========

        formula : EncodedCNF
            The formula to solve. Only the atoms belonging to the theory T
            are handled by the solver; the rest are ignored.

        Returns
        =======

        (solver, lemmas)

        lemmas : list of clauses/lemmas
            Clauses valid in T, so the caller (DPLL(X) in this case) must combine them with F
            before solving.
        """
        raise NotImplementedError

    def notify_assignment(self, lit: Literal, fixed: bool = False) -> None:
        """
        Notify the Theory Solver that the literal has been set true.

        Parameters
        ==========
        lit: The literal that has been set to be true.

        fixed: Whether the assignment is permanent, so that it is never
            undone by a later notify_backtrack().
        """
        raise NotImplementedError

    def notify_new_level(self) -> None:
        """
        Notify DPLL(X) to make a new level. (similar to push() )
        """
        pass

    def notify_backtrack(self, to: int) -> None:
        """
        Notify DPLL(X) to backtrack. (similar to pop(), used in parallel with
        notify_new_level() to do backjumping)
        """
        pass

    def check_model(self, model: Assignment) -> bool:
        """
        Fully check whether if the assignment M is T-consistent.
        This method should be comprehensive, although in the future, for the optimization
        reasons we may need weaker checks.

        Returns
        =======

        True: The assignment M is T-consistent.

        False: The assignment M is T-inconsistent. The conflict clause
            explaining which of the asserted literals cannot hold together
            is then handed over by add_clause().
        """
        raise NotImplementedError

    def decide(self) -> Literal:
        """
        Ask the Theory Solver whether it has its own decision to branch off.
        Otherwise, DPLL(X) will handle it.
        """
        return NO_DECISION

    def propagate(self) -> list[Literal]:
        """
        Ask the Theory Solver whether it has its own propagation to make in the current
        assignment.
        """
        return []

    def provide_reason(self, lit: Literal) -> Explanation:
        """
        Ask the THeory Solver why propagate() was done.
        """
        raise NotImplementedError

    def has_clause(self) -> bool:
        """
        Whether a clause is waiting to be handed over to the engine.
        """
        return False

    def add_clause(self) -> Explanation:
        """
        Hand over one waiting clause. An empty clause means unsatisfiable.
        """
        raise NotImplementedError
