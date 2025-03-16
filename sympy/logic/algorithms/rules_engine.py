from collections.abc import Iterable
from sympy.assumptions.facts2 import id_to_pred

class RulesEngine:
    def __init__(self, full_rules, direct_implications, testing_mode=False):
        """Initialize the rules engine with a nested dictionary structure."""
        self.rule_tree = {}
        self.knowledge_base = {}
        self.encoded_literal_to_pred = {}
        self.testing_mode = testing_mode

        for key, value in full_rules.items():
            self.add_rule(key, value)

        #self.original_keys = list(self.rule_tree.keys())
        self.debugging_print_enabled = False

        self.rules = self.rule_tree.copy()
        self.rules_list = full_rules
        self.direct = direct_implications

    def reset_state(self):
        self.rules = self.rule_tree.copy()
        self.knowledge_base = {}

        #assert list(self.rule_tree.keys()) == self.original_keys


    def add_rule(self, conditions, consequence):
        """
        Add a rule to the rules engine.

        :param conditions: A set of conditions (e.g., {"a", "b", "c", "d"})
        :param consequence: The resulting facts (e.g., {"e"})
        """
        rule_tree = self.rule_tree
        self.print_steps = False

        for condition in sorted(conditions, key=lambda x: str(x)):  # Sort to maintain consistency
            rule_tree = rule_tree.setdefault(condition, {})
        rule_tree["__result__"] = (consequence, conditions)   # Store the consequence

    def _add_facts(self, facts):
        """Internal method to add facts."""
        if self.debugging_print_enabled:
            print(f"adding facts: {[self._to_pred(f)for f in facts]}")
        self.knowledge_base.update(facts)

    def add_facts(self, facts):
        """Add a fact to the knowledge base."""
        if self.debugging_print_enabled:
            self.encoded_literal_to_pred = {list(value)[0]: key for key, value in facts.items()}
        self._add_facts(facts)


    def _to_pred(self, literal):
        pred_id, neg = literal
        pred = id_to_pred[pred_id]
        return ~pred if neg else pred

    def get_negated_fact(self, fact):
        return fact[0], not fact[1]


    def create_explanation_if_conflict(self, fact, antecdents):
        neg_fact = self.get_negated_fact(fact)
        if neg_fact in self.knowledge_base:
            explanation = set.union(*(self.knowledge_base[cond] for cond in antecdents))
            explanation |=  self.knowledge_base[neg_fact]
            if self.debugging_print_enabled:
                print(
                    f"Found conflict triggered by {self._to_pred(fact)} caused by {self._unecode_literals(explanation)}")

            return False, [-lit for lit in explanation]

        return None

    def trigger_rules(self):
        new_facts = []
        for conditions, consequences in self.rules_list.items():
            if any(condition not in self.knowledge_base for condition in conditions):
                continue

            _new_facts = []
            for consequence in consequences:
                if consequence in self.knowledge_base:
                    continue

                res = self.create_explanation_if_conflict(consequence, conditions)
                if res:
                    return res

                _new_facts.append(consequence)


            if _new_facts:
                source_facts = set.union(*(self.knowledge_base[cond] for cond in conditions))
                for fact in _new_facts:
                    self.knowledge_base[fact] = source_facts
            new_facts.extend(_new_facts)

        return True, new_facts

    def trigger(self, trig_fact):
        if trig_fact not in self.rules or trig_fact not in self.knowledge_base:
            return False

        implicants, antecedents = [], []
        pending_facts = set()
        next_rules = self.rules.pop(trig_fact)
        for rule_root in next_rules:
            if rule_root != "__result__":
                self.rules[rule_root] = next_rules[rule_root]

        # there might be new facts that could be triggered
        pending_facts.update(next_rules.keys())
        pending_facts.discard("__result__")


        if "__result__" in next_rules:
            implicants, antecedents = next_rules["__result__"]
            if self.debugging_print_enabled:
                pred_implicants = [self._to_pred(imp) for imp in implicants]
                pred_antecedents = [self._to_pred(ant) for ant in antecedents]
                sources = set(source for ant in antecedents for source in self.knowledge_base[ant])
                for imp in pred_implicants:
                    print(f"Derived {imp} from {pred_antecedents} which was sourced from {self._unecode_literals(sources)}")



        source_facts =  set.union(*(self.knowledge_base[ant] for ant in antecedents)) if len(antecedents) > 0 else set()

        res = (True, None)
        for fact in implicants:
            neg_fact = self.get_negated_fact(fact)
            if neg_fact in self.knowledge_base:
                explanation = source_facts | self.knowledge_base[neg_fact]
                if self.debugging_print_enabled:
                    print(f"Found conflict triggered by {self._to_pred(trig_fact)} from {self._to_pred(fact)} caused by {self._unecode_literals(explanation)}")

                return (False, [-lit for lit in explanation]), source_facts, pending_facts

        pending_facts.update(implicants)
        new_facts = {imp: source_facts for imp in implicants}
        self._add_facts(new_facts)

        return res, source_facts, pending_facts


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

    def _unecode_literals(self, literals):
        ret = []
        for lit in literals:
            pred_id, neg = self.encoded_literal_to_pred[lit]
            new_pred = ~id_to_pred[pred_id] if neg else id_to_pred[pred_id]
            ret.append(new_pred)
        return ret

    def __repr__(self):
        """Return a string representation of the current knowledge base."""
        return f"Knowledge Base: {self.knowledge_base}"

    def get_direct_implications(self, fc_lit):
        return self.direct[fc_lit]


    def trigger_direct_implications(self, fc_lit_to_source_lits):

        pending_facts = set()
        for fc_lit in fc_lit_to_source_lits:
            for implication in self.get_direct_implications(fc_lit):
                if implication in self.knowledge_base:
                    continue
                self.knowledge_base[implication] = fc_lit_to_source_lits[fc_lit]
                pending_facts.add(implication)

        return pending_facts


    def _check_by_rules(self, fc_lit_to_source_lits):
        # loops through rules and continously applies them until no new facts are learned
        self.reset_state()
        queue = self.trigger_direct_implications(fc_lit_to_source_lits)

        while queue:
            res = self.trigger_rules()
            if res[0] is False:
                return res

            queue = res[1]

        return True, None



    def check(self, fc_lit_to_source_lits):
        return self._check_by_rules(fc_lit_to_source_lits)


