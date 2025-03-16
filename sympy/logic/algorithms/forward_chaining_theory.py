from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent,
    Exclusive, Or, to_nnf, BooleanFunction)
from sympy.assumptions.facts import get_number_facts, get_composite_predicates
from collections import defaultdict
from sympy.core.cache import cacheit
from sympy.assumptions import AppliedPredicate, Predicate
from types import MappingProxyType
from functools import reduce
from sympy.assumptions.cnf import Literal

from sympy.logic.algorithms.rules_engine import RulesEngine
from sympy.assumptions.facts2 import (id_rules_dict, rules_dict, id_direct_dict,
                                      id_to_pred, direct_dict, pred_to_id, pred_id_to_bitvec,
                                      pred_id_neg_to_direct_implicants_bitset,
                                      implication_counts_by_lit)


rules_engine = RulesEngine(id_rules_dict, direct_dict)

class FCSolver():
    """
    Theory solver for SymPy's unary facts
    """
    def __init__(self, pred_to_enc = None, testing_mode=False):
        self.engine = rules_engine
        self.engine.reset_state()
        self.enc_to_decomposed_pred = None
        self.expr_count = 0
        self.print_vars = False
        self.theory_prop_enabled = False

        # expr_id to set of unassigned (inactive) variables
        self.unassigned_variables = defaultdict(set)

        expr_to_id = {}
        if pred_to_enc is not None:
            self.enc_to_decomposed_pred = {}
            for appliedPred, enc in pred_to_enc.items():
                if not isinstance(appliedPred, AppliedPredicate) and not isinstance(appliedPred, Not):
                    continue
                expr, pred, neg = self.decompose_AppliedPredicate(appliedPred)
                if expr not in expr_to_id:
                    expr_to_id[expr] = self.expr_count
                    self.expr_count += 1
                assert neg is False
                pred_id = self.predicate_to_pred_id(pred)
                expr_id = expr_to_id[expr]
                self.enc_to_decomposed_pred[enc] = expr_id, pred_id, neg
                self.unassigned_variables[expr_id].add(enc)

        if self.print_vars:
            self.expr_id_to_expr = {value:key for key, value in expr_to_id.items()}

        # disable when finished
        if True:
            self.expr_id_to_expr = {value: key for key, value in expr_to_id.items()}

        self.pred_to_enc = pred_to_enc
        #self.enc_to_pred = {v: k for k, v in pred_to_enc.items()} if pred_to_enc else None

        # dictionary from expr_id to set of active literals associated with expr
        self.active_literals = {expr_id: set() for expr_id in range(self.expr_count)}

        self.testing_mode = testing_mode
        self.conflict = False

        # dictionary of rules with only one antecedent
        # all literals imply themselves
        all_lits = set.union(*rules_dict.values(), {key[0] for key in rules_dict})
        self.direct = id_direct_dict

    @staticmethod
    def from_encoded_cnf(encoded_cnf, testing_mode=False):
        if testing_mode:
            # sort to reduce nondeterminism
            pred_to_enc = dict(sorted(encoded_cnf.encoding.items(), key=lambda x: str(x)))
        else:
            pred_to_enc = encoded_cnf.encoding.copy()

        fc = FCSolver(pred_to_enc)
        state = {expr_id : 0 for expr_id in range(fc.expr_count)}
        # state is a dictionary mapping expr_id to a bitset representing all of the direct
        # implications of the active literals
        return FCSolver(pred_to_enc), state

    def decompose_AppliedPredicate(self, appliedPredicate):
        neg = False
        if isinstance(appliedPredicate, Not):
            neg = True
            appliedPredicate = appliedPredicate.args[0]

        assert isinstance(appliedPredicate, AppliedPredicate)
        expr, pred = appliedPredicate.arg, appliedPredicate.function
        return expr, pred, neg

    def decompose_literal(self, literal):
        assert not (isinstance(literal, AppliedPredicate) or isinstance(literal, Not))

        expr_id, pred_id, neg =  self.enc_to_decomposed_pred[abs(literal)]
        neg = literal < 0
        return expr_id, pred_id, neg

    def literal_to_expr_id(self, literal):
        expr_id, _, _ = self.enc_to_decomposed_pred[abs(literal)]
        return expr_id

    def literal_to_pred_id(self, literal):
        _, pred_id, _ = self.enc_to_decomposed_pred[abs(literal)]
        return pred_id


    def predicate_to_pred_id(self, pred):
        return pred_to_id[pred]



    def lit_in_theory(self, literal):
        return abs(literal) in self.enc_to_decomposed_pred

    def activate_literal(self, new_literal, expr_id, pred_id, neg, state):
        #assert new_literal not in self.active_literals[expr_id]
        if state[expr_id] is None:
            state[expr_id] = 0
            for active_literal in self.active_literals[expr_id]:
                state[expr_id] |= self.get_direct_implicants_bitset_from_literal(active_literal)
        state[expr_id] |= self.get_direct_implicants_bitset_from_literal(new_literal)
        self.active_literals[expr_id].add(new_literal)
        self.unassigned_variables[expr_id].discard(abs(new_literal))

    def deactivate_literal(self, literal, state):
        expr_id, _, _ = self.decompose_literal(literal)
        self.active_literals[expr_id].remove(literal)
        state[expr_id] = None
        self.unassigned_variables[expr_id].add(abs(literal))

    def to_bitvector(self, pred_id, neg):
        return pred_id_to_bitvec[pred_id + neg*len(id_to_pred)]

    def get_direct_implicants_bitset_from_literal(self, literal):
        _, pred_id, neg = self.decompose_literal(literal)
        return self.get_direct_implicants_bitset(pred_id, neg)

    def get_direct_implicants_bitset(self, pred_id, neg):
        return pred_id_neg_to_direct_implicants_bitset[(pred_id, neg)]


    def immediate_conflict(self, expr_id, pred_id, neg, state):
        negated_pred_bv = self.to_bitvector(pred_id, not neg)
        # check if negated pred is in state
        return (negated_pred_bv & state[expr_id]) > 0

    def identify_directly_conflicting_literal(self, new_expr_id, new_pred_id, new_neg):
        # Check each active lit to find the one that's causing the conflict:
        # - Get the set of direct implicants for each active lit.
        # - One of the sets will contain the negated new lit.
        # - The associated active lit is the one responsible for conflcit

        negated_pred_bv = self.to_bitvector(new_pred_id, not new_neg)

        for active_lit in self.active_literals[new_expr_id]:
            bitset = self.get_direct_implicants_bitset_from_literal(active_lit)
            if bitset & negated_pred_bv:
                return active_lit

        # there existing a direct conflict is a precondition to this function
        assert False

    def assert_lit(self, literal, state, source_facts = None, force_assertion=False):
        """
        Adds new fact to asserted dictionary.
        Returns (False, conflict_clause) if there is a 2 literal conflict.

        literal : int
            A literal from a SAT solver.

        state : dict
            A dictionary mapping each `expr_id` to a bitset that represents all currently true direct implications of the literals asserted so far.

        source_fact : set of ints

        force_assertion : bool
            If set to true, will add the literal to the asserted dictionary regardless of whether
            there was a conflict.

        1. Checks if literal is a unary number predicate
        2. Decomposes literal
        3. Check for conflict clause
        4. Add to data stucture keeping track of asserted literals
            - there should be external bitset

        Returns:
            (False, conflict_clause) or (True, None)
        """
        assert isinstance(literal, int)
        if not self.lit_in_theory(literal):
            return True, None

        assert self.conflict is not True

        expr_id, pred_id, neg = self.decompose_literal(literal)

        if self.immediate_conflict(expr_id, pred_id, neg, state):
            conflicting_literal = self.identify_directly_conflicting_literal(expr_id, pred_id, neg)
            if self.print_vars:
               print(f"found conflict caused by: {self._unecode_literals([conflicting_literal, literal])}")
               print()
            conflict_clause = [-conflicting_literal, -literal]
            if force_assertion:
                self.activate_literal(literal, expr_id, pred_id, neg, state)
                # it might be smart to avoid recalculating state here

            return False, conflict_clause

        self.activate_literal(literal, expr_id, pred_id, neg, state)
        return True, None



    def reset_state(self):
        self.engine.reset_state()
        self.active_literals = {expr_id: set() for expr_id in range(self.expr_count)}
        self.conflict = False

    def sanity_check(self, initial_literals=[]):
        res = (True, None)
        for lit in initial_literals:
            res = self.assert_lit(lit)
            if not res[0]:
                return res
        return res

    def theory_prop(self, new_lit):
        if not self.theory_prop_enabled:
            return []
        if not self.lit_in_theory(new_lit):
            return []

        #expr_id, pred_id, neg = self.decompose_literal(new_lit)
        expr_id = self.literal_to_expr_id(new_lit)

        # if self.print_vars:
        #     new_pred = ~id_to_pred[pred_id] if neg else id_to_pred[pred_id]

        implicants = []
        direct_implicants_set = self.get_direct_implicants_bitset_from_literal(new_lit)
        for unassigned_variable in self.unassigned_variables[expr_id]:
            if not self.lit_in_theory(unassigned_variable):
                continue

            uv_pred_id = self.literal_to_pred_id(unassigned_variable)
            #print(cur_pred)

            # Try positive polarity first, then negative if needed
            for polarity in (False, True):
                literal_bitvec = self.to_bitvector(uv_pred_id, polarity)
                if literal_bitvec & direct_implicants_set:
                    imp = -unassigned_variable if polarity else unassigned_variable
                    implicants.append(imp)
                    break


        return implicants





        # For all unassigned literals that have the same expr id:
        # Check if their negation is in the set of direct implicants of new lit/

    def check(self, initial_literals, pre_encoded=True):
        assert len(initial_literals) > 0

        fc_lit_to_sat_lit = {}

        def sat_lit_to_fc_lit(sat_lit):
            _, pred_id, neg = self.decompose_literal(sat_lit)
            return pred_id, neg

        expr_id_to_fc_literal = defaultdict(list)

        for sat_lit in initial_literals:
            if not self.lit_in_theory(sat_lit):
                continue
            fc_lit = sat_lit_to_fc_lit(sat_lit)
            fc_lit = {fc_lit : {sat_lit}}
            expr_id, _, _ = self.decompose_literal(sat_lit)
            expr_id_to_fc_literal[expr_id].append(fc_lit)

        for expr_id, fc_lits in expr_id_to_fc_literal.items():
            fc_lits = reduce(lambda a, b: {**a, **b}, fc_lits)
            if len(fc_lits) <= 2:
                continue
            res = self.engine.check(fc_lits)
            if res[0] is False:
                if self.print_vars:
                    print(f"Found conflict caused by: {self._unecode_literals(res[1])}")
                    print()
                return res

        return True, None


    def _unecode_literals(self, literals):
        ret = []
        for lit in literals:
            if not self.lit_in_theory(lit):
                continue
            expr_id, pred_id, neg = self.decompose_literal(lit)

            new_expr = self.expr_id_to_expr[expr_id]
            new_pred = ~id_to_pred[pred_id](new_expr) if neg else id_to_pred[pred_id](new_expr)

            ret.append(new_pred)
        return ret

    def get_heuristic_multipliers(self, variable):
        if not self.lit_in_theory(variable):
            return 1, 1
        expr_id, pred_id, _ = self.enc_to_decomposed_pred[variable]
        positive_weighting = implication_counts_by_lit[(pred_id, False)]#*len(self.unassigned_variables[expr_id])
        negative_weighting = implication_counts_by_lit[(pred_id, True)]#*len(self.unassigned_variables[expr_id])
        return positive_weighting, negative_weighting

    # def starts_negated(self, variable):
    #     pred_id = self.literal_to_pred_id(variable)
    #     return should_pred_start_negated[pred_id]




