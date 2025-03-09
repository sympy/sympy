from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent,
    Exclusive, Or, to_nnf, BooleanFunction)
from sympy.assumptions.facts import get_number_facts, get_composite_predicates
from collections import defaultdict
from sympy.core.cache import cacheit
from sympy.assumptions import AppliedPredicate, Predicate

from sympy.core.symbol import Symbol
from sympy.strategies.core import switch


class RulesEngine:
    def __init__(self, dictionary):
        """Initialize the rules engine with a nested dictionary structure."""
        self.rules = {}
        self.knowledge_base = {}

        for key, value in dictionary.items():
            self.add_rule(key, value)

    def add_rule(self, conditions, consequence):
        """
        Add a rule to the rules engine.

        :param conditions: A set of conditions (e.g., {"a", "b", "c", "d"})
        :param consequence: The resulting fact (e.g., "e")
        """
        rule_tree = self.rules
        for condition in sorted(conditions):  # Sort to maintain consistency
            rule_tree = rule_tree.setdefault(condition, {})
        rule_tree["__result__"] = consequence  # Store the consequence

    def add_fact(self, fact, source_facts):
        """Add a fact to the knowledge base."""
        self.knowledge_base[fact] = source_facts

    def traverse_rule_tree(self, rule_tree, current_facts):
        """Recursively traverse the rule tree to check for matching conditions."""
        if "__result__" in rule_tree:  # Consequence found
            return rule_tree["__result__"]
        for key, sub_tree in rule_tree.items():
            if key in current_facts:  # Move deeper if condition exists
                result = self.traverse_rule_tree(sub_tree, current_facts)
                if result:
                    return result
        return None

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
        rules_dict[antecedent].add(implicant)
        rules_dict[to_nnf(~implicant)].add(to_nnf(~antecedent))  # contrapositive
        return

    antecedent = to_nnf(antecedent, simplify=False)
    implicant = to_nnf(implicant, simplify=False)
    if isinstance(antecedent, Or) or isinstance(implicant, And):
        if isinstance(implicant, And):
            antecedent, implicant = to_nnf(~implicant), to_nnf(~antecedent)

        for disjunct in antecedent.args:
            _add_rule(rules_dict, disjunct, implicant, remove_var=False)
        return

    assert _is_literal(antecedent) or (all(_is_literal(arg) for arg in antecedent.args) and isinstance(antecedent, And)), antecedent

    if  isinstance(antecedent, And) or isinstance(implicant, Or):
        if isinstance(implicant, Or):
            antecedent, implicant = to_nnf(~implicant), to_nnf(~antecedent)

        rules_dict[antecedent].add(implicant)

        assert len(antecedent.args) >= 2
        for i in range(len(antecedent.args)):
            new = And(*antecedent.args[:i] + antecedent.args[i+1:])
            rules_dict[And(~implicant, new)].add(~antecedent.args[i])

    else:
        rules_dict[antecedent].add(implicant)
        rules_dict[to_nnf(~implicant)].add(~antecedent)

#@cacheit
def facts_to_dictionary(x = None):

    rules_dict = defaultdict(set)

    facts = get_number_facts()
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

from sympy.assumptions.ask_generated import get_known_facts_dict

rules_dict = facts_to_dictionary()

bad_keys = [key for key in rules_dict.keys()
            if not(isinstance(key, Predicate) or isinstance(key.args[0], Predicate))]
assert len(bad_keys) == 0, bad_keys

# print(blah)
#dic = get_known_facts_dict()


def check_consistency(initial_literals):
  """
  Parameters
  ==========

  initial_literals (list): A list of initial literals (facts) to seed the knowledge base.

  Returns
  =======

  tuple:
      - (bool): `True` if the knowledge base is consistent, `False` if a contradiction is found.
      - (set or None): If inconsistent, a set of conflicting facts; otherwise, `None`.
  """
  # A dictionary keeping track of all of the literals known to be implied by
  # the initial list of literals. Each known literal is mapped to its "source
  # literals", a subset of `initial_literals` that implies that literal.
  knowledge_base = {}
  rules_engine = RulesEngine(rules_dict)

  def add_new_fact(new_fact, source_facts):
    assert type(source_facts) == set
    if ~new_fact in rules_engine.knowledge_base:
      # we have found a contradiciton: some literal and its negation are true
      return False, source_facts | knowledge_base[~new_fact]
    rules_engine.add_fact(new_fact, source_facts)
    return True, None

  for lit in initial_literals:
    # initial facts are their own source facts.
    res = add_new_fact(lit, {lit})
    if res[0] is False:
      return res

  queue = initial_literals
  while queue:
    pending_facts = set()
    for antecedent in queue:
      print(f"Checking {antecedent}")
      if antecedent not in rules_dict:
        print(f"\t{antecedent} not in rules")
        continue
      for implicant in rules_dict[antecedent]:
        if implicant in knowledge_base:
          print(f"\t{antecedent} already known")
          continue
        print(f"\tDeriving {implicant} from {antecedent}")
        source_facts = knowledge_base[antecedent]
        new_fact = implicant
        res = add_new_fact(new_fact, source_facts)
        if res[0] is False:
          return res
        pending_facts.add(new_fact)

    queue = pending_facts

  return True, None

from sympy import Q
print(check_consistency([~Q.positive, Q.real, ~Q.negative, ~Q.zero]))
#print(check_consistency([Q.positive, ~Q.real]))

