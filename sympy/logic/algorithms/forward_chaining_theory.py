from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent,
    Exclusive, Or, to_nnf, BooleanFunction)
from sympy.assumptions.facts import get_number_facts, get_composite_predicates
from collections import defaultdict
from sympy.core.cache import cacheit
from sympy.assumptions import AppliedPredicate, Predicate
from types import MappingProxyType

from sympy.core.symbol import Symbol
from sympy.strategies.core import switch
from collections.abc import Iterable

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

# make sure all literals imply themselves
# for lit in set.union(*rules_dict.values(), {key[0] for key in rules_dict}).:
#     rules_dict[key].add(key[0])

rules_dict = transitive_closure(rules_dict)
rules_dict = dict(sorted(rules_dict.items(), key=lambda item: str(item)))

# bad_keys = [key for key in rules_dict.keys()
#             if not(isinstance(key, Predicate) or isinstance(key[0], Predicate))]
#assert len(bad_keys) == 0, bad_keys

# print(blah)
#dic = get_known_facts_dict()
from collections.abc import Iterable


class RulesEngine:
    def __init__(self, dictionary):
        """Initialize the rules engine with a nested dictionary structure."""
        self.rule_tree = {}
        self.knowledge_base = {}

        for key, value in dictionary.items():
            self.add_rule(key, value)

        #self.original_keys = list(self.rule_tree.keys())

        self.rules = self.rule_tree.copy()

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

        for condition in sorted(conditions, key=lambda x: str(x)):  # Sort to maintain consistency
            rule_tree = rule_tree.setdefault(condition, {})
        rule_tree["__result__"] = (consequence, conditions)   # Store the consequence

    def add_fact(self, fact, source_facts):
        """Add a fact to the knowledge base."""
        self.knowledge_base[fact] = source_facts

    def add_facts(self, facts):
        """Add a fact to the knowledge base."""
        self.knowledge_base.update(facts)

    def trigger(self, fact):
        if fact not in self.rules or fact not in self.knowledge_base:
            return False

        if "__result__" in self.rules[fact]:
            implicants, antecedents = self.rules[fact]["__result__"]
            #implicants = [imp for imp in implicants if imp not in self.knowledge_base]
        else:
            implicants, antecedents = [], []


        pending_facts = set()
        next = self.rules.pop(fact)
        for rule_root in next:
            if rule_root != "__result__":
                self.rules.update({rule_root: next[rule_root]})
                pending_facts.update(next[rule_root].keys())

        return implicants, antecedents, pending_facts

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

rules_engine = RulesEngine(rules_dict)

class FCSolver():
    """
    Theory solver for SymPy's unary facts
    """
    def __init__(self, pred_to_enc = None, testing_mode=False):
        self.engine = rules_engine
        self.engine.reset_state()
        self.pred_to_enc = pred_to_enc
        self.enc_to_pred = {v: k for k, v in pred_to_enc.items()} if pred_to_enc else None
        self.asserted = defaultdict(dict)
        self.testing_mode = testing_mode
        self.conflict = False

        # dictionary of rules with only one antecedent
        # all literals imply themselves
        all_lits = set.union(*rules_dict.values(), {key[0] for key in rules_dict})
        self.direct = defaultdict(set, {key[0]: val | {key[0]} for key, val in rules_dict.items() if len(key) == 1})
        for lit in all_lits:
            self.direct[lit].add(lit)

    @staticmethod
    def from_encoded_cnf(encoded_cnf, testing_mode=False):
        if testing_mode:
            # sort to reduce nondeterminism
            pred_to_enc = dict(sorted(encoded_cnf.encoding.items(), key=lambda x: str(x)))
        else:
            pred_to_enc = encoded_cnf.encoding.copy()


        return FCSolver(pred_to_enc)

    def decomposte_fact(self, fact):
        neg = False
        if self.enc_to_pred is not None:
            if fact < 0:
                neg = True
            fact = self.enc_to_pred[abs(fact)]
        elif isinstance(fact, Not):
            neg = True
            fact = fact.args[0]

        assert isinstance(fact, AppliedPredicate)
        expr, fact = fact.arg, fact.function
        fact = ~fact if neg else fact
        return expr, fact

    def unassert_lit(self, fact):
        expr, fact = self.decomposte_fact(fact)
        if fact not in self.asserted[expr]:
            return False
            assert fact in self.asserted[expr]
        del self.asserted[expr][fact]
        self.conflict = None


    def assert_lit(self, new_fact, source_facts = None, force_assertion=False):

        #print(f"Asserting {new_fact} in assert_lit")

        assert self.conflict is not True
        if source_facts is None:
            # initial facts are their own source facts.
            source_facts = {new_fact}

        if isinstance(new_fact, int) and not isinstance(self.enc_to_pred[abs(new_fact)], AppliedPredicate):
            return True, None

        expr, new_fact = self.decomposte_fact(new_fact)

        return self._assert_lit(expr, new_fact, source_facts, external=True, force_assertion=force_assertion)

    def _assert_lit(self, expr, new_fact, source_facts, external=False, force_assertion=False):
        assert type(source_facts) == set

        # note: each fact implies itself and is included in its list of implicants
        for implicant in self.direct[new_fact]:
            if ~implicant in self.asserted[expr]:
                # we have found a contradiciton: some literal and its negation are true
                conflict_clause = [lit for lit in source_facts | self.asserted[expr][~implicant]]
                if isinstance(conflict_clause[0], int):
                    conflict_clause = [-lit for lit in conflict_clause]
                else:
                    conflict_clause = [~lit for lit in conflict_clause]
                assert len(conflict_clause) <= 4
                if external:
                    assert len(conflict_clause) == 2
                #print(len(conflict_clause))
                if self.testing_mode:
                    conflict_clause = sorted(conflict_clause, key=lambda x: str(x))
                self.conflict= True
                if force_assertion:
                    self.asserted[expr][new_fact] = source_facts
                return False, conflict_clause

        self.asserted[expr][new_fact] = source_facts
        #print(f"{new_fact} asserted in _assert_lit()")
        return True, None




    def reset_state(self):
        self.engine.reset_state()
        self.asserted = defaultdict(dict)
        self.conflict = False

    def sanity_check(self, initial_literals=[]):
        res = (True, None)
        for lit in initial_literals:
            res = self.assert_lit(lit)
            if not res[0]:
                return res
        return res

    def check(self, initial_literals = [], res = None):

        # res = self.sanity_check(initial_literals)
        if res is not None and res[0] is False:
            return res

        for lit in []:
            if isinstance(lit, int) and not isinstance(self.enc_to_pred[abs(lit)], AppliedPredicate):
                continue
            # initial facts are their own source facts.
            res = self.assert_lit(lit, {lit})
            if res[0] is False:
                assert len(res[1]) == 2
                return res

        res = None
        assert len(self.asserted) > 0
        for expr in self.asserted:
            self.engine.reset_state()
            res = self._check_expr(expr)
            if res[0] is False:
                return res
        return res

    def _check_expr(self, expr):
        assignment = self.asserted[expr].copy()
        return self._check_assignment(expr)


    def _check_assignment(self, expr):
        queue = self.asserted[expr].copy()

        self.engine.add_facts(queue)

        # for fact in queue:
        #     self.engine.add_fact(fact, self.asserted[expr][fact])

        while queue:
            pending_facts = set()
            for antecedent in queue:
                res = self.engine.trigger(antecedent)
                if not res:
                    continue
                new_facts, antecedents, new_pending_facts = res
                source_facts = set.union(*[self.engine.knowledge_base[ant] for ant in antecedents], set())

                for fact in new_facts:
                    if ~fact in self.engine.knowledge_base:
                        explanation = self.engine.knowledge_base[antecedent] | self.engine.knowledge_base[~fact]
                    #res = self._assert_lit(expr, fact, source_facts)
                        assert len(explanation) >= 2
                        return False, [-lit for lit in explanation]
                    pending_facts.add(fact)
                    self.engine.add_fact(fact, source_facts)

                pending_facts.update(new_pending_facts)

            queue = pending_facts

        return True, None


