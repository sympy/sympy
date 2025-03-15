from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent,
    Exclusive, Or, to_nnf, BooleanFunction)
from sympy.assumptions.facts import get_number_facts, get_composite_predicates
from collections import defaultdict
from sympy.core.cache import cacheit
from sympy.assumptions import AppliedPredicate, Predicate
from types import MappingProxyType
from collections.abc import Iterable
from functools import reduce
from sympy.assumptions.cnf import Literal

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



# class RuleTree:
#     def __init__(self, rules):
