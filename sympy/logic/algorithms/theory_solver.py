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

Atom = int
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

    # Whether T only decides complete assignments
    is_lazy: bool = False
    # Whether DPLL(X) may discard the reasons given by provide_reason()
    are_reasons_forgettable: bool = False
    # The atoms that belong to the theory T.
    observed: set[Atom]

    # This is the only engine-internal method, i.e this is only used in internals of Theory Solver engine. All other methods are used between DPLL(X) and Theory Solvers
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

    # ------ Methods for communicating between DPLL(X) and Theory Solvers
    def notify_assignment(self, lits: Assignment) -> None:
        """
        DPLL(X) notifies that the literals have been set true. Batched, so a
        conflict noticed here is handed (should handed) over by has_clause()/add_clause().
        """
        raise NotImplementedError

    def notify_new_decision_level(self) -> None:
        """
        A push/pop pair (notify/backtrack are pus/pop respectively). DPLL(X) notifies that it has opened a new decision
        level, or that it has already backtracked to `to`. Levels are unnamed,
        so the Theory Solver counts them itself, and `to` is an absolute level
        whose every assertion above must be undone.
        """
        raise NotImplementedError

    def notify_backtrack(self, to: int) -> None:
        raise NotImplementedError

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
        Ask the Theory Solver whether if has its own decision to branch off.
        Otherwise, DPLL(X) will handle it.
        """
        return NO_DECISION

    def propagate(self) -> list[Literal]:
        """
        The literals that T has derived to be true. DPLL(X) assigns them itself,
        so returning nothing is always allowed.
        """
        return []

    def provide_reason(self, lit: Literal) -> Explanation:
        """
        Why a literal handed over by propagate() had to be true, as a clause
        valid in T.
        """
        raise NotImplementedError

    def has_clause(self) -> tuple[bool, bool]:
        """
        Whether a clause is waiting to be handed over to the engine, and
        whether the engine may forget it once added.
        """
        raise NotImplementedError

    def add_clause(self) -> Explanation:
        """
        Hand over one waiting clause. An empty clause means unsatisfiable.
        """
        raise NotImplementedError
