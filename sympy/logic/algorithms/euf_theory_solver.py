"""
Theory solver for EUF, to be used as the T in a DPLL(T) loop.

The congruence closure engine in ``euf_theory`` only ever hears about
equalities that hold.  This module turns it into a theory solver: it picks the
equality atoms out of an ``EncodedCNF``, feeds the ones the SAT solver assigns
to true into the closure, and reports a conflict as soon as the closure makes
the two sides of an asserted disequality congruent.  The proof of that
congruence, translated back into literals, is the conflict clause.

The interface (``from_encoded_cnf``, ``assert_lit``, ``check`` and ``reset``)
is the one ``LRASolver`` offers, so ``dpll2`` drives both solvers identically.

References
==========

.. [1] Nieuwenhuis, R., Oliveras, A., Tinelli, C.:
       Solving SAT and SAT Modulo Theories
       https://www.cs.upc.edu/~oliveras/espai/papers/dpllt.pdf

Classes
-------
    EUFTheorySolver: Wraps EUFCongruenceClosure for the DPLL(T) loop.
"""
from __future__ import annotations

from sympy.assumptions.ask import Q
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.symbol import Dummy
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure

# The constant a predicate atom is equated to when the atom is asserted true.
TRUE = Dummy('true')


class EUFTheorySolver:
    """
    Decides conjunctions of equalities and disequalities between ground terms.

    A predicate atom such as ``Q.positive(x)`` is reified as the equality
    ``Q.positive(x) = TRUE``, so congruence reaches predicates as well: once
    ``a = b`` holds, ``Q.positive(a)`` and ``Q.positive(b)`` have to agree.
    Anything that is not an applied predicate is left alone, as EUF says
    nothing about it.
    """

    def __init__(self, atom_id_to_equality):
        """
        Use the "from_encoded_cnf" method to create a new EUFTheorySolver.
        """
        # atom id -> (lhs, rhs, whether the atom asserts lhs = rhs when true)
        self.atom_id_to_equality = atom_id_to_equality
        self.reset()

    @staticmethod
    def from_encoded_cnf(encoded_cnf):
        """
        Create an EUFTheorySolver from an EncodedCNF object and a list of
        conflict clauses for the atoms that are decided without the closure.

        Parameters
        ==========

        encoded_cnf : EncodedCNF

        Returns
        =======

        (euf, conflicts)

        euf : EUFTheorySolver

        conflicts : list
            Contains a one-literal conflict clause for each atom that is
            already true or false on its own, such as ``Q.eq(x, x)``.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy.assumptions.ask import Q
        >>> from sympy.assumptions.cnf import CNF, EncodedCNF
        >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
        >>> cnf = CNF.from_prop(Q.eq(x, y))
        >>> enc = EncodedCNF()
        >>> enc.from_cnf(cnf)
        >>> euf, conflicts = EUFTheorySolver.from_encoded_cnf(enc)
        >>> conflicts
        []
        >>> euf.atom_id_to_equality
        {1: (x, y, True)}
        >>> euf.assert_lit(1) is None
        True
        >>> euf.check()[0]
        True
        """
        atom_id_to_equality = {}
        conflicts = []

        for prop, atom_id in encoded_cnf.encoding.items():
            if prop == True:
                conflicts.append([atom_id])
                continue
            if prop == False:
                conflicts.append([-atom_id])
                continue
            if not isinstance(prop, AppliedPredicate):
                continue
            if prop.function not in (Q.eq, Q.ne):
                atom_id_to_equality[atom_id] = (prop, TRUE, True)
                continue

            is_equality = prop.function == Q.eq
            if prop.lhs == prop.rhs:
                # x = x holds and x != x does not, whatever else is asserted.
                conflicts.append([atom_id if is_equality else -atom_id])
                continue
            atom_id_to_equality[atom_id] = (prop.lhs, prop.rhs, is_equality)

        return EUFTheorySolver(atom_id_to_equality), conflicts

    def reset(self):
        """
        Reset the solver to the state it had before anything was asserted.
        """
        self.cc = EUFCongruenceClosure([])
        # the equations merged so far, so that a proof can name its literals
        self.eq_literal = {}
        self.disequalities = []
        self.result = None

    def assert_lit(self, literal):
        """
        Assert a literal and update the closure accordingly.

        Literals over atoms that are not equalities are ignored.

        Parameters
        ==========

        literal : int
            A mapping of ids to equalities can be found in
            ``self.atom_id_to_equality``.

        Returns
        =======

        None or (False, explanation)

        explanation : list of ints
            A conflict clause that "explains" why the literals asserted so
            far cannot all hold.
        """
        atom = self.atom_id_to_equality.get(abs(literal))
        if atom is None:
            return None

        lhs, rhs, is_equality = atom
        if (literal > 0) == is_equality:
            self.eq_literal[Q.eq(lhs, rhs)] = literal
            self.cc.merge(lhs, rhs)
        else:
            self.disequalities.append((literal, lhs, rhs))

        if self.result is None:
            self.result = self._find_conflict()
        return self.result

    def check(self):
        """
        Report whether the literals asserted so far are consistent, with a
        conflict clause that "explains" why if they are not.

        Returns
        =======

        (True, classes) or (False, explanation)

        classes : dict
            Maps every term the closure knows to its class representative.

        explanation : list of ints
        """
        if self.result is not None:
            return self.result
        return True, dict(self.cc.representative)

    def _find_conflict(self):
        """
        A disequality is violated the moment the closure makes its two sides
        congruent.  There are few disequalities per branch compared to the
        work the closure itself does, so they are simply scanned.
        """
        for literal, lhs, rhs in self.disequalities:
            if not self.cc.are_congruent(lhs, rhs):
                continue
            # The equalities that prove lhs = rhs, together with the
            # disequality itself, are exactly what cannot hold at once.
            proof = self.cc.explain(lhs, rhs)
            return False, [-literal] + [-self.eq_literal[eq] for eq in proof]
        return None
