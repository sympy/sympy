from collections.abc import Iterable
class RulesEngine:
    def __init__(self, full_rules, direct_implications):
        """Initialize the rules engine with a nested dictionary structure."""
        self.rule_tree = {}
        self.knowledge_base = {}

        for key, value in full_rules.items():
            self.add_rule(key, value)

        #self.original_keys = list(self.rule_tree.keys())

        self.rules = self.rule_tree.copy()
        self.direct = direct_implications

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
