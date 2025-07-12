"""
sympy.logic.algorithms.euf_theory
=================================

A congruence-closure based solver for Equality with Uninterpreted Functions (EUF),
using union-find for equivalence classes and leveraging SymPy's known fact system
for predicate propagation.

Classes
-------
- EUFUnhandledInput: Exception for unsupported or contradictory input.
- EUFVariable: Represents a variable in the EUF solver, with union-find structure.
- EUFSolver: Main solver class; supports fact assignment, equality, disequality, and queries.

Doctests
--------
>>> from sympy.assumptions.ask import Q
>>> from sympy import symbols
>>> from sympy.logic.algorithms.euf_theory import EUFSolver
>>> x, y = symbols('x y')
>>> solver = EUFSolver()
>>> solver.assign_fact(solver.get_var(x), Q.prime, True)
>>> solver.query(Q.prime, x)
True
>>> solver.assign_fact(solver.get_var(y), Q.composite, True)
>>> solver = EUFSolver()
>>> solver.assign_fact(solver.get_var(x), Q.even, True)
>>> solver.query(Q.integer, x)
True
>>> solver.assign_fact(solver.get_var(y), Q.odd, True)
>>> solver = EUFSolver()
>>> solver.assign_fact(solver.get_var(x), Q.composite, True)
>>> solver.query(Q.prime, x)
False
>>> solver.query(Q.integer, x)
True
>>> solver.query(Q.positive, x)
True
"""

from sympy.assumptions.ask import Q
from sympy.assumptions.assume import AppliedPredicate
from collections import defaultdict
from sympy.assumptions.ask_generated import get_known_facts_dict

class EUFUnhandledInput(Exception):
    """Exception raised for unsupported or contradictory input in EUF solver."""
    pass

class EUFVariable:
    """
    Represents a variable in the EUF solver, supporting union-find operations.

    Attributes
    ----------
    parent : EUFVariable
        Parent pointer for union-find.
    rank : int
        Rank for union-by-rank optimization.
    predicates : defaultdict
        Mapping from predicate to polarity (True/False).
    symbol : Symbol
        The SymPy symbol this variable represents.

    Methods
    -------
    add_predicate(pred, polarity)
        Assigns a predicate with the given polarity, raising on contradiction.
    """
    def __init__(self, symbol):
        self.parent = self
        self.rank = 0
        self.predicates = defaultdict(bool)  # {Predicate: polarity}
        self.symbol = symbol

    def add_predicate(self, pred, polarity):
        """
        Assign a predicate with the given polarity to this variable's class.

        Raises
        ------
        ValueError
            If a contradictory predicate is already present.
        """
        if self.predicates.get(pred, None) == (not polarity):
            raise ValueError(f"Contradiction: {pred} and ~{pred}")
        self.predicates[pred] = polarity

class EUFSolver:
    """
    EUF theory solver implementing congruence closure and fact propagation.

    Attributes
    ----------
    variables : dict
        Maps symbols to EUFVariable instances.
    disequalities : set
        Set of disequality constraints (pairs of symbols).
    known_facts : dict
        Mapping from predicates to (implied, rejected) sets, from SymPy's fact system.

    Methods
    -------
    get_var(symbol)
        Returns the EUFVariable for a given symbol, creating it if necessary.
    find(var)
        Finds the root of the equivalence class for a variable.
    union(var1, var2)
        Merges the equivalence classes of two variables, propagating predicates.
    assign_fact(var, predicate, polarity)
        Assigns a predicate with polarity to a variable, propagating implications.
    _apply_known_facts(predicate, polarity, root)
        Propagates all implied and rejected facts for a predicate.
    add_equality(a, b)
        Declares that symbols a and b are equal.
    add_disequality(a, b)
        Declares that symbols a and b are disequal.
    check_consistency()
        Returns True if the current state is consistent, False otherwise.
    query(predicate, symbol)
        Queries the value of a predicate for a given symbol.
    from_encoded_cnf(enc_cnf)
        Constructs a solver from an encoded CNF (SymPy internal).
    """
    def __init__(self):
        self.variables = {}
        self.disequalities = set()
        self.known_facts = get_known_facts_dict()

    def get_var(self, symbol):
        """
        Get or create the EUFVariable for a given symbol.

        Parameters
        ----------
        symbol : Symbol
            The SymPy symbol.

        Returns
        -------
        EUFVariable
        """
        return self.variables.setdefault(symbol, EUFVariable(symbol))

    def find(self, var):
        """
        Find the root representative of the equivalence class for a variable.

        Parameters
        ----------
        var : EUFVariable

        Returns
        -------
        EUFVariable
        """
        if var.parent != var:
            var.parent = self.find(var.parent)
        return var.parent

    def union(self, var1, var2):
        """
        Merge the equivalence classes of two variables, propagating all predicates.

        Parameters
        ----------
        var1, var2 : EUFVariable

        Raises
        ------
        ValueError
            If merging would create a contradiction.
        """
        root1 = self.find(var1)
        root2 = self.find(var2)
        if root1 == root2:
            return

        # Union by rank
        if root1.rank < root2.rank:
            root1, root2 = root2, root1

        # Merge predicates and check for conflicts
        for pred, pol in root2.predicates.items():
            if pred in root1.predicates:
                if root1.predicates[pred] != pol:
                    raise ValueError(f"Conflict in merged predicates: {pred}")
            else:
                root1.predicates[pred] = pol

        root2.parent = root1
        if root1.rank == root2.rank:
            root1.rank += 1

    def assign_fact(self, var, predicate, polarity):
        """
        Assign a predicate with polarity to a variable and propagate implications.

        Parameters
        ----------
        var : EUFVariable
        predicate : Predicate
        polarity : bool

        Raises
        ------
        EUFUnhandledInput
            If a contradiction is encountered.
        """
        root = self.find(var)
        try:
            root.add_predicate(predicate, polarity)
        except ValueError as e:
            raise EUFUnhandledInput(str(e))

        # Apply known fact implications and rejections
        self._apply_known_facts(predicate, polarity, root)

    def _apply_known_facts(self, predicate, polarity, root):
        """
        Apply all implications and rejections for a predicate from known facts.

        Parameters
        ----------
        predicate : Predicate
        polarity : bool
        root : EUFVariable
        """
        if predicate in self.known_facts:
            accepted_facts = self.known_facts[predicate][0]
            rejected_facts = self.known_facts[predicate][1]
            for fact in accepted_facts:
                root.add_predicate(fact, polarity)
            for fact in rejected_facts:
                root.add_predicate(fact, not polarity)

    def add_equality(self, a, b):
        """
        Assert that two symbols are equal (merge their classes).

        Parameters
        ----------
        a, b : Symbol
        """
        var_a = self.get_var(a)
        var_b = self.get_var(b)
        self.union(var_a, var_b)

    def add_disequality(self, a, b):
        """
        Assert that two symbols are disequal.

        Parameters
        ----------
        a, b : Symbol
        """
        self.disequalities.add((a, b))

    def check_consistency(self):
        """
        Check if the current state is consistent (no disequality conflicts).

        Returns
        -------
        bool
        """
        for a, b in self.disequalities:
            if self.find(self.get_var(a)) == self.find(self.get_var(b)):
                return False
        return True

    def query(self, predicate, symbol):
        """
        Query the value of a predicate for a symbol.

        Parameters
        ----------
        predicate : Predicate
        symbol : Symbol

        Returns
        -------
        bool or None
            True/False if known, None if undetermined.

        Examples
        --------
        >>> from sympy.assumptions.ask import Q
        >>> from sympy import symbols
        >>> from sympy.logic.algorithms.euf_theory import EUFSolver
        >>> solver = EUFSolver()
        >>> x, y = symbols('x y')
        >>> solver.assign_fact(solver.get_var(x), Q.even, True)
        >>> solver.query(Q.integer, x)
        True
        """
        root = self.find(self.get_var(symbol))
        if predicate in root.predicates:
            return root.predicates[predicate]
        return None

    @classmethod
    def from_encoded_cnf(cls, enc_cnf):
        """
        Construct a solver from an encoded CNF (SymPy internal).

        Parameters
        ----------
        enc_cnf : EncodedCNF

        Returns
        -------
        EUFSolver
        """
        solver = cls()
        for lit, enc in enc_cnf.encoding.items():
            if isinstance(lit, AppliedPredicate):
                pred = lit.function
                args = lit.arguments
                if pred in (Q.eq, Q.ne):
                    a, b = args
                    if pred == Q.eq:
                        solver.add_equality(a, b) if enc > 0 else solver.add_disequality(a, b)
                    else:
                        solver.add_disequality(a, b) if enc > 0 else solver.add_equality(a, b)
                elif len(args) == 1:
                    solver.assign_fact(solver.get_var(args[0]), pred, enc > 0)
        return solver
