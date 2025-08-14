from collections import defaultdict
from sympy import Eq, Unequality, Not
from sympy.logic import boolalg
from sympy.assumptions.ask import Q, AppliedPredicate
from sympy.core.relational import Relational
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure

def _order_key(expr):
    return (expr.__class__.__name__, hash(expr))

def _ordered_pair(a, b):
    return (a, b) if _order_key(a) <= _order_key(b) else (b, a)

def _canonical_lit(lit):
    if isinstance(lit, Not):
        eq = lit.args[0]
        if not isinstance(eq, Eq):
            raise TypeError("Unsupported negated literal")
        return eq.lhs, eq.rhs, False
    elif isinstance(lit, Unequality):
        return lit.lhs, lit.rhs, False
    elif isinstance(lit, Eq):
        return lit.lhs, lit.rhs, True
    elif isinstance(lit, AppliedPredicate):
        if lit.function == Q.eq:
            return lit.arguments[0], lit.arguments[1], True
        elif lit.function == Q.ne:
            return lit.arguments[0], lit.arguments[1], False
        else:
            raise TypeError(f"Unsupported AppliedPredicate {lit}")
    else:
        raise TypeError(f"Unsupported literal type {type(lit)} for EUF")

def _canon_eq(lhs, rhs):
    return Eq(lhs, rhs) if _order_key(lhs) <= _order_key(rhs) else Eq(rhs, lhs)

class EUFTheorySolver:
    def __init__(self):
        self.interpretation_stack = []      # [(canonical_eq, ts)]
        self.timestamp_counter = 0
        self.literal_timestamp = {}
        self.disequalities_set = defaultdict(set)
        self.disequality_causes = {}
        self.positive_literal_list = defaultdict(list)
        self.negative_literal_list = defaultdict(list)
        self.literal_terms = {}            # canonical_eq -> (lhs, rhs, is_pos)
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}
        self._enc_to_lit = {}               # encoding int -> sympy literal
        self._lit_to_enc = {}               # sympy literal -> encoding int

    @classmethod
    def from_encoded_cnf(cls, encoded_cnf, testing_mode=False):
        solver = cls()
        conflict_clauses = []
        encoded_items = (
            sorted(encoded_cnf.encoding.items(), key=lambda x: str(x[0]))
            if testing_mode else encoded_cnf.encoding.items()
        )
        literals = set()
        for prop, enc in encoded_items:
            solver._enc_to_lit[enc] = prop
            solver._lit_to_enc[prop] = enc
            if isinstance(prop, boolalg.BooleanTrue):
                conflict_clauses.append([enc]); continue
            if isinstance(prop, boolalg.BooleanFalse):
                conflict_clauses.append([-enc]); continue
            if not isinstance(prop, (Eq, Unequality, Not, Relational, AppliedPredicate)):
                raise TypeError(f"Unsupported predicate in EUF: {prop}")
            lhs, rhs, is_positive = _canonical_lit(prop)
            # trivial
            if lhs == rhs:
                conflict_clauses.append([enc] if is_positive else [-enc]); continue
            if lhs.is_Number and rhs.is_Number:
                conflict_clauses.append([-enc] if is_positive else [enc]); continue
            literals.add(prop)
        solver.Initialize(literals)
        return solver, conflict_clauses

    def Initialize(self, literal_set):
        for lit in literal_set:
            lhs, rhs, is_positive = _canonical_lit(lit)
            ceq = _canon_eq(lhs, rhs)
            self.literal_terms[ceq] = (lhs, rhs, is_positive)
            lhs_c = self.cc._flatten(lhs)
            rhs_c = self.cc._flatten(rhs)
            if is_positive:
                self.positive_literal_list[lhs_c].append(ceq)
                self.positive_literal_list[rhs_c].append(ceq)
            else:
                self.negative_literal_list[lhs_c].append(ceq)
                self.negative_literal_list[rhs_c].append(ceq)

    def assert_lit(self, enc_lit):
        sign = enc_lit > 0
        enc = abs(enc_lit)
        if enc not in self._enc_to_lit:
            return None  # literal not relevant to theory
        lit = self._enc_to_lit[enc]
        if not sign:
            lit = Not(lit)
        try:
            self.SetTrue(lit)
        except Exception:
            return False, self.Explanation(lit)
        return None

    def SetTrue(self, literal):
        self.timestamp_counter += 1
        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))
        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)
        consequences = set()
        if is_positive:
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)
            k1, k2 = _ordered_pair(rep_lhs, rep_rhs)
            self.cc.history[(k1, k2)] = (lhs_c, rhs_c, ceq)
            self.cc.add_equality(lhs, rhs)
            for (da, db), _ in self.disequality_causes.items():
                if self.cc._find(da) == self.cc._find(db):
                    raise ValueError("EUF conflict: equality contradicts disequality")
            consequences |= self._propagate_positive(rep_lhs)
            consequences |= self._propagate_positive(rep_rhs)
        else:
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)
            if rep_lhs == rep_rhs:
                raise ValueError("EUF conflict: disequality within same class")
            self.disequalities_set[rep_lhs].add(rep_rhs)
            self.disequalities_set[rep_rhs].add(rep_lhs)
            self.disequality_causes[(rep_lhs, rep_rhs)] = ceq
            self.disequality_causes[(rep_rhs, rep_lhs)] = ceq
            consequences |= self._propagate_negative(rep_lhs)
            consequences |= self._propagate_negative(rep_rhs)
        self.interpretation_stack.append((ceq, self.timestamp_counter))
        self.literal_timestamp[ceq] = self.timestamp_counter
        return consequences

    def check(self):
        for (a, b), lit in self.disequality_causes.items():
            if self.cc.are_equal(a, b):
                return False, {lit}
        return True, {}

    def Explanation(self, literal):
        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))
        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)
        explanation_set = set()
        if is_positive:
            self._explain_path(lhs_c, rhs_c, explanation_set,
                               self.literal_timestamp.get(ceq, 0))
        else:
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)
            cause_ceq = self.disequality_causes.get((rep_lhs, rep_rhs))
            if cause_ceq is not None:
                explanation_set.add(cause_ceq)
            self._explain_path(lhs_c, rep_lhs, explanation_set,
                               self.literal_timestamp.get(ceq, 0))
            self._explain_path(rhs_c, rep_rhs, explanation_set,
                               self.literal_timestamp.get(ceq, 0))
        return explanation_set

    def _propagate_positive(self, rep):
        res = set()
        for ceq in self.positive_literal_list[rep]:
            if self.IsTrue(ceq):
                res.add(ceq)
        for ceq in self.negative_literal_list[rep]:
            if self.IsTrue(ceq):
                res.add(ceq)
        return res

    def _propagate_negative(self, rep):
        res = set()
        for ceq in self.negative_literal_list[rep]:
            if self.IsTrue(ceq):
                res.add(ceq)
        return res

    def IsTrue(self, literal):
        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))
        if is_positive:
            return self.cc.are_equal(lhs, rhs)
        else:
            lhs_c = self.cc._flatten(lhs)
            rhs_c = self.cc._flatten(rhs)
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)
            if rep_lhs == rep_rhs:
                return False
            if ceq in self.disequality_causes.values():
                return True
            return None

    def _explain_path(self, a, b, out_set, max_ts):
        rep_a = self.cc._find(a)
        rep_b = self.cc._find(b)
        if rep_a == rep_b:
            return
        k1, k2 = _ordered_pair(rep_a, rep_b)
        if (k1, k2) not in self.cc.history:
            return
        _, _, cause_ceq = self.cc.history[k1, k2]
        if cause_ceq is not None and self.literal_timestamp.get(cause_ceq, 0) <= max_ts:
            out_set.add(cause_ceq)
            lhs_term, rhs_term, _ = self.cc.history[k1, k2]
            self._explain_path(lhs_term, rhs_term, out_set, max_ts)

    def Backtrack(self, n):
        for _ in range(n):
            if not self.interpretation_stack:
                break
            lit, _ = self.interpretation_stack.pop()
            self.literal_timestamp.pop(lit, None)
        remaining = [lit for lit, _ in self.interpretation_stack]
        # Rebuild
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}
        self.disequalities_set.clear()
        self.disequality_causes.clear()
        for ceq in remaining:
            lhs, rhs, is_pos = self.literal_terms[ceq]
            self.SetTrue(Eq(lhs, rhs) if is_pos else Not(Eq(lhs, rhs)))



