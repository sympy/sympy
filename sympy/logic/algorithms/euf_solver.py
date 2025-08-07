"""
EUF Theory Solver for DPLL(T) integrating EUFCongruenceClosure.

This solver implements the theory solver interface expected by SymPy's SAT solver (DPLL(X))
and uses the congruence closure algorithm described in:

Nieuwenhuis & Oliveras, "Congruence Closure with Integer Offsets"
https://www.cs.upc.edu/~oliveras/dpllt.pdf

API
---
- from_encoded_cnf(encoded_cnf): classmethod that creates solver with CNF encoding
- assert_lit(literal): asserts a SAT literal into the theory solver
- check(): returns (True, proof) if consistent, (False, conflict_clause) if conflict
- explanation(literal): returns a list of literals explaining the conflict or assertion
- reset(): resets internal state for reuse

Integration
-----------
Supply this class instance as `euf_theory` in SymPy's SATSolver constructor for DPLL(T).

Example
-------
>>> from sympy.logic.algorithms.euf_theory import EUFSolver
>>> from sympy.assumptions.cnf import CNF, EncodedCNF
>>> from sympy import symbols, Eq, Function
>>> f = Function('f')
>>> a, b, x = symbols('a b x')
>>> phi = Eq(a, b) & Eq(f(a), x)
>>> cnf = CNF.from_prop(phi)
>>> enc = EncodedCNF()
>>> enc.from_cnf(cnf)
>>> euf_solver, immediate_conflicts = EUFSolver.from_encoded_cnf(enc)
>>> # Now supply euf_solver to SATSolver as euf_theory
"""

from sympy.core.relational import Eq, Ne
from sympy.core.symbol import Symbol
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
from collections import defaultdict

class EUFSolver:
    """
    Equality with Uninterpreted Functions (EUF) theory solver for DPLL(T).

    Wraps EUFCongruenceClosure to provide incremental equality reasoning.
    """
    def __init__(self, enc_to_eq, eq_to_enc):
        self.enc_to_eq = enc_to_eq          # int -> Eq/Ne object
        self.eq_to_enc = eq_to_enc          # Eq/Ne object -> int encoding
        self.cc = EUFCongruenceClosure([])  # Core congruence closure engine
        self.pos_eqs = set()                # set of tuples (a,b) for asserted equalities
        self.diseqs = set()                 # set of tuples (a,b) for asserted disequalities
        self.lit_stack = []                 # ordered stack of asserted literals
        self.is_sat = True
        self.conflict_clause = None
        self.explanation_map = defaultdict(set)  # literal -> set of literals explaining it

    @classmethod
    def from_encoded_cnf(cls, encoded_cnf, testing=False):
        enc_to_eq = {}
        eq_to_enc = {}
        immediate_conflicts = []
        for atom, enc in encoded_cnf.encoding.items():
            if isinstance(atom, (Eq, Ne)):
                # Canonical ordering for consistency
                try:
                    a, b = sorted((atom.lhs, atom.rhs), key=str)
                except Exception:
                    a, b = atom.lhs, atom.rhs
                norm_eq = Eq(a, b) if isinstance(atom, Eq) else Ne(a, b)
                enc_to_eq[enc] = norm_eq
                eq_to_enc[norm_eq] = enc
            elif atom is True:
                immediate_conflicts.append([enc])
            elif atom is False:
                immediate_conflicts.append([-enc])
        solver = cls(enc_to_eq, eq_to_enc)
        return solver, immediate_conflicts

    def reset(self):
        """Reset the solver state for reuse."""
        self.cc = EUFCongruenceClosure([])
        self.pos_eqs.clear()
        self.diseqs.clear()
        self.lit_stack.clear()
        self.is_sat = True
        self.conflict_clause = None
        self.explanation_map.clear()

    def assert_lit(self, enc_lit):
        """
        Assert literal (encoded int). Negative means negated literal (disequality).
        Returns None if no conflict, else (False, conflict_clause).
        """
        abs_lit = abs(enc_lit)
        is_neg = enc_lit < 0
        if abs_lit not in self.enc_to_eq:
            # Non-EUF atoms or bool constants ignored here
            return None
        atom = self.enc_to_eq[abs_lit]
        self.lit_stack.append(enc_lit)

        try:
            a, b = sorted((atom.lhs, atom.rhs), key=str)
        except Exception:
            a, b = atom.lhs, atom.rhs

        canonical_eq = Eq(a, b)
        canonical_ne = Ne(a, b)

        if not is_neg:
            # Positive equality
            if (a, b) in self.diseqs or (b, a) in self.diseqs:
                # Conflict: contradictory equality and disequality asserted
                self.is_sat = False
                self.conflict_clause = [-self.eq_to_enc[canonical_eq], -self.eq_to_enc[canonical_ne]]
                return False, self.conflict_clause
            self.cc.add_equality(a, b)
            self.pos_eqs.add((a, b))
            self.explanation_map[enc_lit] = {enc_lit}
        else:
            # Disequality
            self.diseqs.add((a, b))
            if self.cc.are_equal(a, b):
                # Conflict: equality already known
                self.is_sat = False
                self.conflict_clause = [-self.eq_to_enc[canonical_eq], -self.eq_to_enc[canonical_ne]]
                return False, self.conflict_clause
            self.explanation_map[enc_lit] = {enc_lit}
        return None

    def check(self):
        """
        Check consistency of asserted literals.
        Returns (True, proof) if satisfiable, or (False, conflict_clause) if conflict.
        """
        if self.is_sat:
            return True, list(self.lit_stack)
        else:
            return False, self.conflict_clause

    def explanation(self, enc_lit):
        """
        Retrieve explanation of why literal is implied or conflict.
        Returns list of encoded literals explaining the conflict.
        """
        if self.is_sat:
            return []
        else:
            # Minimal explanation: return stored conflict clause
            return self.conflict_clause

