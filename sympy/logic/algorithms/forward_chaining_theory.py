from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent,
    Exclusive, Or, to_nnf, BooleanFunction)
from sympy.assumptions.facts import get_number_facts, get_composite_predicates
from collections import defaultdict
from sympy.core.cache import cacheit
from sympy.assumptions import AppliedPredicate, Predicate
from types import MappingProxyType
from functools import reduce
from sympy.assumptions.cnf import Literal


from sympy.core.symbol import Symbol
from sympy.strategies.core import switch
from collections.abc import Iterable


class BitSet:
    def __init__(self, value=0):
        """Initialize the BitSet with an integer value (default: 0)."""
        self.value = value

    def add(self, x):
        """Add an element to the bit set."""
        self.value |= (1 << x)

    def remove(self, x):
        """Remove an element from the bit set."""
        self.value &= ~(1 << x)

    def contains(self, x):
        """Check if an element is in the bit set."""
        return (self.value & (1 << x)) != 0

    def toggle(self, x):
        """Toggle (flip) an element in the bit set."""
        self.value ^= (1 << x)

    def size(self):
        """Return the number of elements in the bit set."""
        return bin(self.value).count('1')

    def to_list(self):
        """Convert the bit set to a list of integers."""
        return [i for i in range(self.value.bit_length()) if self.contains(i)]

    def clear(self):
        """Clear the bit set (reset to empty)."""
        self.value = 0

    def union(self, other):
        """Return a new BitSet with the union of two sets."""
        return BitSet(self.value | other.value)

    def intersection(self, other):
        """Return a new BitSet with the intersection of two sets."""
        return BitSet(self.value & other.value)

    def difference(self, other):
        """Return a new BitSet with the difference (self - other)."""
        return BitSet(self.value & ~other.value)

    def __repr__(self):
        """String representation of the BitSet."""
        return f"BitSet({self.to_list()})"

def _AppliedPredicate_to_Predicate(pred):
    if isinstance(pred, AppliedPredicate):
        return pred.function
    if isinstance(pred, Predicate):
        return pred
    if isinstance(pred, BooleanFunction):
        return pred.func(*[_AppliedPredicate_to_Predicate(arg) for arg in pred.args])

    raise ValueError(f"Object {pred} is not of a recognized type")

def _is_literal(expr):
    return isinstance(expr, Predicate) or (isinstance(expr.args[0], Predicate) and isinstance(expr, Not)), f"{expr} is not a literal"

def _add_rule(rules_dict, antecedent, implicant, remove_var = True):



    if remove_var:
        antecedent = _AppliedPredicate_to_Predicate(antecedent)
        implicant = _AppliedPredicate_to_Predicate(implicant)
    if isinstance(antecedent, Predicate) and isinstance(implicant, Predicate):
        rules_dict[(antecedent,)].add(implicant)
        rules_dict[(to_nnf(~implicant),)].add(to_nnf(~antecedent))  # contrapositive
        if len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) > 0:
            assert len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) == 0
        return

    antecedent = to_nnf(antecedent, simplify=False)
    implicant = to_nnf(implicant, simplify=False)
    if isinstance(antecedent, Or) or isinstance(implicant, And):
        if isinstance(implicant, And):
            antecedent, implicant = to_nnf(~implicant), to_nnf(~antecedent)

        for disjunct in antecedent.args:
            _add_rule(rules_dict, disjunct, implicant, remove_var=False)

        if len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) > 0:
            assert len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) == 0
        return

    assert _is_literal(antecedent) or (all(_is_literal(arg) for arg in antecedent.args) and isinstance(antecedent, And)), antecedent

    if  isinstance(antecedent, And) or isinstance(implicant, Or):
        if isinstance(implicant, Or):
            antecedent, implicant = to_nnf(~implicant), to_nnf(~antecedent)

        assert isinstance(antecedent, And)
        rules_dict[antecedent.args].add(implicant)

        assert len(antecedent.args) >= 2
        for i in range(len(antecedent.args)):
            new = And(*antecedent.args[:i] + antecedent.args[i+1:])
            assert isinstance(And(~implicant, new), And)
            rules_dict[And(~implicant, new).args].add(~antecedent.args[i])
        if len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) > 0:
            assert len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) == 0

    else:
        assert not isinstance(to_nnf(antecedent), And) and not isinstance(to_nnf(antecedent), Or)
        assert not isinstance(to_nnf(implicant), And) and not isinstance(to_nnf(implicant), Or)
        rules_dict[(antecedent,)].add(implicant)
        rules_dict[(to_nnf(~implicant),)].add(~antecedent)
        if len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) > 0:
            assert len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) == 0

    if len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) > 0:
        assert len([key for key in rules_dict.keys() if not isinstance(key, Iterable)]) == 0

from sympy.assumptions import Q
# add additional rules:
additional_rules = And(
    Implies(Q.nonzero, ~Q.zero),
    Implies(Q.nonpositive, ~Q.positive),
    Implies(Q.nonnegative, ~Q.negative),
    Implies(~Q.real, ~Q.nonzero),
    Implies(~Q.real, ~Q.nonpositive),
    Implies(~Q.real, ~Q.nonnegative),
)
#@cacheit
def facts_to_dictionary(x = None):

    rules_dict = defaultdict(set)

    facts = And(get_number_facts(), additional_rules)

    for fact in list(facts.args):
        if isinstance(fact, Implies):
            _add_rule(rules_dict, fact.args[0], fact.args[1])
        elif isinstance(fact, Equivalent):
            _add_rule(rules_dict, fact.args[0], fact.args[1])
            _add_rule(rules_dict, fact.args[1], fact.args[0])
        elif isinstance(fact, Not) and isinstance(fact.args[0], And):
            fact = fact.args[0]
            _add_rule(rules_dict, fact.args[0], ~fact.args[1])
        else:
            raise ValueError(f"{fact} is not of a recognized form")

    composite_predicate = get_composite_predicates()
    for superset, subsets in composite_predicate.items():
        _add_rule(rules_dict, subsets, superset)
        _add_rule(rules_dict, superset, subsets)
        # for subset in subsets.args:
        #     _add_rule(rules_dict, subset, superset)




    return rules_dict


from collections import deque


def transitive_closure(graph):
    closure = {node: set(neighbors) for node, neighbors in graph.items()}

    for start_node in graph:
        visited = set()
        queue = deque([start_node])

        while queue:
            node = queue.popleft()
            if node in visited or len(node) > 1:
                continue
            visited.add(node)

            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    queue.append((neighbor,))
                    closure[start_node].add(neighbor)

    return {node: set(neighbors) for node, neighbors in closure.items()}



from sympy.assumptions.ask_generated import get_known_facts_dict

rules_dict = facts_to_dictionary()
rules_dict = transitive_closure(rules_dict)
rules_dict = dict(sorted(rules_dict.items(), key=lambda item: str(item)))

# dictionary of rules with only one antecedent
# all literals imply themselves
all_lits = set.union(*rules_dict.values(), {key[0] for key in rules_dict})
direct_dict = defaultdict(set, {key[0]: val | {key[0]} for key, val in rules_dict.items() if len(key) == 1})
for lit in all_lits:
    direct_dict[lit].add(lit)

all_preds = [lit for lit in all_lits if not isinstance(lit, Not)]
all_preds = sorted(all_preds, key=lambda pred: str(pred))

# if we give more commonly used preds lower numbers, this will help
pred_to_id = {all_preds[i] : i for i in range(len(all_preds))}
def pred_to_id_neg_tup(pred):
    if isinstance(pred, Not):
        return -pred_to_id[pred.args[0]]#, True
    else:
        return pred_to_id[pred]#, False

id_to_pred = [pred for pred, _ in sorted(pred_to_id.items(), key=lambda x: x[1])]
num_preds = len(all_preds)

id_direct_dict = {pred_to_id_neg_tup(pred) : {pred_to_id_neg_tup(im) for im in imps } for pred, imps in direct_dict.items()}
id_rules_dict = {}
for ante, imps in rules_dict.items():
    id_rules_dict[tuple(pred_to_id_neg_tup(pred) for pred in ante)] = set(pred_to_id_neg_tup(imp) for imp in imps)

#id_rules_dict = { tuple(pred_to_id_neg_tup(pred) for pred in ante) : {pred_to_id_neg_tup(im) for im in imps} for ante, imps in rules_dict.items()}


# def pred_lit_set_to_bit_set(pred_lit_set):
#     i = pred_lit_set_to_int(pred_lit_set)
#     return BitSet(i)



def pred_lit_set_to_int(pred_lit_set):
    binary = []
    for i in range(2*len(all_preds)):
        neg = i >= len(all_preds)
        lit = id_to_pred[i % len(all_preds)]
        if neg:
            lit = ~lit

        in_set = (lit in pred_lit_set)*1
        binary.append(str(in_set))

    return int("".join(binary), 2)

direct_dict_bitset = {}
for lit in all_lits:
    direct_imps_bitset = pred_lit_set_to_int(direct_dict[lit])
    neg = False
    if isinstance(lit, Not):
        pred_id = pred_to_id[lit.args[0]]
        lit = -pred_to_id[lit.args[0]]
        neg = True
    else:
        pred_id = pred_to_id[lit]

    direct_dict_bitset[(pred_id, neg)] = direct_imps_bitset


pred_id_to_bitvec = [pred_lit_set_to_int({pred}) for pred in id_to_pred]
pred_id_to_bitvec += [pred_lit_set_to_int({~pred}) for pred in id_to_pred]
pred_id_neg_to_direct_implicants_bitset = direct_dict_bitset

#pred_id_direct_dict = {pred_to_id[pred] : 2  for pred, implications in direct_dict.items()}

# I want a mapping of pred number to an int representing the negated preds implied by it
# I need to give negated preds a seperate number

# I can create an int representing a signle pred or negated pred
# then I can AND that with the int representing the current state.
# If non-zero, then the pred is in the set.

direct_dict = dict(direct_dict)
direct_dict = MappingProxyType(direct_dict)


from collections.abc import Iterable

# class RuleTree:
#     def __init__(self, rules):


class RulesEngine:
    def __init__(self, dictionary):
        """Initialize the rules engine with a nested dictionary structure."""
        self.rule_tree = {}
        self.knowledge_base = {}

        for key, value in dictionary.items():
            self.add_rule(key, value)

        #self.original_keys = list(self.rule_tree.keys())

        self.rules = self.rule_tree.copy()
        self.direct = direct_dict

    def reset_state(self):
        self.rules = self.rule_tree.copy()
        self.knowledge_base = {}

        #assert list(self.rule_tree.keys()) == self.original_keys


    def add_rule(self, conditions, consequence):
        """
        Add a rule to the rules engine.

        :param conditions: A set of conditions (e.g., {"a", "b", "c", "d"})
        :param consequence: The resulting fact (e.g., "e")
        """
        rule_tree = self.rule_tree
        self.print_steps = False

        for condition in sorted(conditions, key=lambda x: str(x)):  # Sort to maintain consistency
            rule_tree = rule_tree.setdefault(condition, {})
        rule_tree["__result__"] = (consequence, conditions)   # Store the consequence

    def add_fact(self, fact, source_facts):
        """Add a fact to the knowledge base."""
        print(f"adding fact: {self._to_pred(fact)}")
        self.knowledge_base[fact] = source_facts

    def add_facts(self, facts):
        """Add a fact to the knowledge base."""
        if self.print_steps:
            print(f"adding facts: {[self._to_pred(f)for f in facts]}")
        self.knowledge_base.update(facts)

    def _to_pred(self, literal):
        pred_id, neg = literal.arg, literal.is_Not
        pred = id_to_pred[pred_id]
        return ~pred if neg else pred

    def trigger(self, trig_fact):
        if trig_fact not in self.rules or trig_fact not in self.knowledge_base:
            return False

        implicants, antecedents = [], []
        pending_facts = set()
        next = self.rules.pop(trig_fact)
        for rule_root in next:
            if rule_root != "__result__":
                self.rules.update({rule_root: next[rule_root]})
                pending_facts.update(next[rule_root].keys())
            else:
                implicants, antecedents = next["__result__"]
                if self.print_steps:
                    pred_implicants = [self._to_pred(imp) for imp in implicants]
                    pred_antecedents = [self._to_pred(ant) for ant in antecedents]
                    sources = set(source for ant in antecedents for source in self.knowledge_base[ant])
                    for imp in pred_implicants:
                        print(f"Derived {imp} from {pred_antecedents} which was sourced from {sources}")


        source_facts = set.union(*[self.knowledge_base[ant] for ant in antecedents], set())

        res = (True, None)
        for fact in implicants:
            if -fact in self.knowledge_base:
                if self.print_steps:
                    explanation = self.knowledge_base[trig_fact] | self.knowledge_base[-fact]
                    print(f"Found conflict triggered by {self._to_pred(trig_fact)} from {self._to_pred(fact)} caused by {explanation}")
                    #print(f"Found conflict caused by {explanation}")
                explanation = self.knowledge_base[trig_fact] | self.knowledge_base[-fact]
                # res = self._assert_lit(expr, fact, source_facts)
                #assert len(explanation) >= 2
                res = (False, [-lit for lit in explanation])
                break
            pending_facts.add(fact)
            self.add_fact(fact, source_facts)

        return res, source_facts, pending_facts

    # def traverse_rule_tree(self, rule_tree, current_facts):
    #     """Recursively traverse the rule tree to check for matching conditions."""
    #     if "__result__" in rule_tree:  # Consequence found
    #         return rule_tree["__result__"]
    #     for key, sub_tree in rule_tree.items():
    #         if key in current_facts:  # Move deeper if condition exists
    #             result = self.traverse_rule_tree(sub_tree, current_facts)
    #             if result:
    #                 return result
    #     return None

    def infer_facts(self):
        """Infer new facts based on existing knowledge."""
        new_facts = set()


        # Traverse from each known fact
        for fact in self.knowledge_base:
            result = self.traverse_rule_tree(self.rules.get(fact, {}), self.knowledge_base)
            if result and result not in self.knowledge_base:
                new_facts.add(result)

        # Update knowledge base
        self.knowledge_base.update(new_facts)
        return new_facts

    def run_inference_until_stable(self):
        """Keep running inference until no new facts are inferred."""
        while True:
            new_facts = self.infer_facts()
            if not new_facts:
                break

    def __repr__(self):
        """Return a string representation of the current knowledge base."""
        return f"Knowledge Base: {self.knowledge_base}"

rules_engine = RulesEngine(id_rules_dict)

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

        if pred_to_enc is not None:
            self.enc_to_decomposed_pred = {}
            expr_to_id = {}
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

        if self.print_vars:
            new_pred = ~id_to_pred[pred_id] if neg else id_to_pred[pred_id]

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
            return Literal(pred_id, neg)

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
            res = self._check_assignment(fc_lits)
            if res[0] is False:
                return res

        return True, None


    def _check_assignment(self, fc_lit_to_source_lits):

        queue = fc_lit_to_source_lits
        self.engine.add_facts(fc_lit_to_source_lits)

        while queue:
            pending_facts = set()
            for antecedent in queue:
                res = self.engine.trigger(antecedent)
                if not res:
                    continue
                results, source_facts, new_pending_facts = res

                if results[0] is False:
                    if self.testing_mode:
                        results = False, [self.engine._to_pred(lit) for lit in results[1]]
                    return results

                pending_facts.update(new_pending_facts)

            queue = pending_facts

        return True, None

    def _unecode_literals(self, literals):
        ret = []
        for lit in literals:
            if not self.lit_in_theory(lit):
                continue
            expr_id, pred_id, neg = self.decompose_literal(lit)

            new_pred = ~id_to_pred[pred_id] if neg else id_to_pred[pred_id]
            ret.append(new_pred)
        return ret



