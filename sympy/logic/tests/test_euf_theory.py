from __future__ import annotations
from sympy.logic.algorithms.euf_theory import (EUFApp, EUFCongruenceClosure,
    EUFEquation)
from sympy.testing.pytest import raises
from sympy.core.random import choice, randint, sample, shuffle
from collections import defaultdict
import random


def fn(name):
    """Return a builder for applications of the uninterpreted function ``name``."""
    return lambda *args: EUFApp(name, args)


f, g, h = fn('f'), fn('g'), fn('h')
add, mul = fn('+'), fn('*')
x, y, z, w, a, b, c, d = 'x y z w a b c d'.split()


def test_basic_and_chain_equality():
    cc = EUFCongruenceClosure([EUFEquation(x, y), EUFEquation(y, z)])
    assert cc.are_congruent(x, y)
    assert cc.are_congruent(y, z)
    assert cc.are_congruent(x, z)
    assert not cc.are_congruent(x, w)


def test_unary_function_congruence():
    cc = EUFCongruenceClosure([])
    for t in (f(a), f(b), x):
        cc._flatten(t)
    cc.merge(a, b)
    cc.merge(f(a), x)
    assert cc.are_congruent(f(b), x)   # f(a) = x, a=b -> f(b) = x


def test_binary_congruence_and_propagation():
    cc = EUFCongruenceClosure([])
    for t in (g(a, c), g(b, d), x):
        cc._flatten(t)
    cc.merge(a, b)
    cc.merge(c, d)
    cc.merge(g(a, c), x)
    assert cc.are_congruent(g(b, d), x)  # g(a,c) = x; a=b, c=d -> g(b,d)=x


def test_currying_shares_partial_applications():
    # f(a, b) is stored as EUFApp(EUFApp(f, (a,)), (b,)), so the two spellings
    # are the same term and the partial application is a single constant.
    cc = EUFCongruenceClosure([])
    assert cc._flatten(g(a, b)) == cc._flatten(EUFApp(g(a), (b,)))
    assert cc._flatten(g(a)) == cc._flatten(EUFApp('g', (a,)))


def test_equal_functions_are_congruent():
    # currying puts the head in an argument position, so heads can be merged
    cc = EUFCongruenceClosure([EUFEquation('f', 'g')])
    assert cc.are_congruent(f(a), g(a))
    assert not cc.are_congruent(f(a), g(b))
    cc.merge(a, b)
    assert cc.are_congruent(f(a), g(b))


def test_equal_functions_propagate_to_every_arity():
    cc = EUFCongruenceClosure([EUFEquation('f', 'g')])
    assert cc.are_congruent(f(a), g(a))
    assert cc.are_congruent(f(a, b), g(a, b))
    assert cc.are_congruent(f(f(a)), g(g(a)))
    assert not cc.are_congruent(f(a), g(b))


def test_equal_functions_are_not_extensional():
    # f(t) = g(t) for every t in play does NOT give f = g back
    cc = EUFCongruenceClosure([EUFEquation(f(a), g(a)), EUFEquation(f(b), g(b))])
    assert not cc.are_congruent('f', 'g')


def test_partial_application_used_as_a_function():
    # g(a) is an ordinary term, so it can be equated to a head and applied
    cc = EUFCongruenceClosure([EUFEquation(g(a), 'f')])
    assert cc.are_congruent(g(a, b), f(b))


def test_explain_walks_the_function_side():
    # the congruence edge here holds because the heads are equal, so the proof
    # has to recurse into the function half of the application, not the argument
    eqs = [EUFEquation('f', 'g'), EUFEquation(f(a), x), EUFEquation(g(a), y)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, x, y)
    assert EUFEquation('f', 'g') in expl


def test_backtrack_retracts_function_equality():
    cc = EUFCongruenceClosure([])
    cc.merge('f', 'g')
    assert cc.are_congruent(f(a), g(a))
    cc.backtrack(1)
    assert not cc.are_congruent(f(a), g(a))
    assert not cc.are_congruent(f(a, b), g(a, b))


def test_congruence_is_not_injectivity():
    cc = EUFCongruenceClosure([EUFEquation(h(x, y), h(y, x))])
    assert cc.are_congruent(h(x, y), h(y, x))
    assert not cc.are_congruent(x, y)


def test_permuted_arguments_no_commutativity():
    cc = EUFCongruenceClosure([
        EUFEquation(h(x, y), h(y, x)),     # h(x,y) = h(y,x)
        EUFEquation(x, y)                  # x = y
    ])
    # Even without commutativity, if x=y, h(x,y)=h(y,x) by congruence
    assert cc.are_congruent(h(x, y), h(y, x))


def test_add_equality_registers_and_merges():
    cc = EUFCongruenceClosure([])

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))

    cc.merge(a, b)
    assert cc._find_repr(cc._flatten(a)) == cc._find_repr(cc._flatten(b))
    assert cc._find_repr(fa) == cc._find_repr(fb)


def test_distinct_heads_flatten_unique():
    cc = EUFCongruenceClosure([])
    assert cc._flatten(f(x)) != cc._flatten(g(x))


def test_flatten_application_and_cache():
    cc = EUFCongruenceClosure([])
    testf = fn('testf')
    ax = cc._flatten(testf(x))
    bx = cc._flatten(testf(x))
    assert ax == bx


def test_flatten_shares_partial_applications():
    cc = EUFCongruenceClosure([])
    cc._flatten(g(a, b))
    cc._flatten(g(a, c))
    # both terms went through one constant for g(a), so only three exist:
    # g(a), g(a)(b) and g(a)(c)
    assert cc._flatten(g(a)) == cc._flatten(EUFApp('g', (a,)))
    assert len(cc._app_to_const) == 3


def test_flatten_nested_spellings_agree():
    cc = EUFCongruenceClosure([])
    nary = cc._flatten(h(a, b, c))
    assert nary == cc._flatten(EUFApp(EUFApp(EUFApp('h', (a,)), (b,)), (c,)))
    assert nary == cc._flatten(EUFApp(h(a, b), (c,)))


def test_flatten_head_is_a_term_too():
    cc = EUFCongruenceClosure([])
    fx = cc._flatten(f(x))
    assert cc._app_to_const[EUFApp(cc._flatten('f'), (cc._flatten(x),))] == fx


def test_flatten_atom_and_application_never_collide():
    cc = EUFCongruenceClosure([])
    assert cc._flatten('f') != cc._flatten(f(a))
    assert cc._flatten(f(a)) != cc._flatten(a)


def test_flatten_is_stable():
    cc = EUFCongruenceClosure([])
    t = 't'
    assert cc._find_repr(cc._flatten(t)) == cc._flatten(t)
    assert cc._flatten(f(t)) == cc._flatten(f(t))


def test_process_pending_chain_merges():
    cc = EUFCongruenceClosure([])
    f1 = fn('alsof')
    x1, y1, z1 = 'x1', 'y1', 'z1'
    fx, fy, fz = (cc._flatten(f1(t)) for t in (x1, y1, z1))
    for p, q in ((x1, y1), (y1, z1)):
        cc.pending.append((cc._flatten(p), cc._flatten(q), EUFEquation(p, q)))
    cc._process_pending_unions()
    assert cc._find_repr(cc._flatten(x1)) == cc._find_repr(cc._flatten(y1)) \
        == cc._find_repr(cc._flatten(z1))
    assert cc._find_repr(fx) == cc._find_repr(fy) == cc._find_repr(fz)


def test_flatten_atoms_are_opaque():
    cc = EUFCongruenceClosure([])

    # the same atom always maps to the same constant
    assert cc._flatten('a') == cc._flatten('a')
    # distinct atoms map to distinct constants, whatever they are
    assert cc._flatten('a') != cc._flatten('b')
    assert cc._flatten(1) != cc._flatten('a')
    assert cc._flatten(1) == cc._flatten(1)


def test_use_list_merging_under_union():
    cc = EUFCongruenceClosure([])
    a1, b1, c1 = 'a1', 'b1', 'c1'
    f1 = fn('f1')
    # Register applications before any merging
    apps = [cc._flatten(f1(t)) for t in (a1, b1, c1)]
    cc.merge(a1, b1)
    cc.merge(b1, c1)
    # Test: all applications f1(x) for all class members x are congruent
    app_reps = {cc._find_repr(app) for app in apps}
    assert len(app_reps) == 1


def test_complex_deep_chaining():
    # Deep nesting of f
    depth = 190
    cc = EUFCongruenceClosure([])
    term_a = a
    term_b = b
    for _ in range(2 * depth):
        term_a = f(term_a)
        term_b = f(term_b)
    # declare both towers up front, then merge
    cc._flatten(term_a)
    cc._flatten(term_b)
    term_a = a
    term_b = b
    for _ in range(depth):
        term_a = f(term_a)
        term_b = f(term_b)
    cc._flatten(x)
    cc.merge(term_a, x)
    cc.merge(a, b)

    # All nestings over a and b should be equal to each other and to x
    for _ in range(depth):
        assert cc.are_congruent(term_a, term_b)
        term_a = f(term_a)
        term_b = f(term_b)


def test_long_chain_variables():
    names = ['a%s' % i for i in range(20)]
    eqs = [EUFEquation(names[i], names[i+1]) for i in range(len(names)-1)]
    cc = EUFCongruenceClosure(eqs)
    for i in range(len(names)):
        for j in range(len(names)):
            assert cc.are_congruent(names[i], names[j])


def test_composed_functions():
    cc = EUFCongruenceClosure([])
    for t in (f(a), f(b), g(c), h(b), g(f(b))):
        cc._flatten(t)
    cc.merge(a, b)
    cc.merge(f(a), c)
    cc.merge(g(c), h(b))
    # f(a) = c and a=b => f(b) = c
    assert cc.are_congruent(f(a), f(b))
    # g(c) = h(b) and c = f(a) = f(b)
    assert cc.are_congruent(g(f(b)), h(b))


def test_example_1():
    eq1 = EUFEquation(f(a), g(b))                  # f(a) = g(b)
    eq2 = EUFEquation(g(c), h(f(c), g(a)))         # g(c) = h(f(c), g(a))
    eq3 = EUFEquation(b, c)                        # b = c
    eq4 = EUFEquation(f(c), g(a))                  # f(c) = g(a)
    eq5 = EUFEquation(h(d, d), g(b))               # h(d, d) = g(b)
    eq6 = EUFEquation(g(a), d)                     # g(a) = d

    cc = EUFCongruenceClosure([eq1, eq2, eq3, eq4, eq5, eq6])

    # Assertions checking congruence closure identifies equalities properly
    assert cc.are_congruent(b, c)                  # b = c
    assert cc.are_congruent(g(a), d)               # g(a) = d
    assert cc.are_congruent(g(b), g(c))            # g(b) = g(c)
    assert cc.are_congruent(f(a), g(c))            # f(a) = g(c)


def test_example_2():
    eqs = [
        EUFEquation(f(a), g(b)),                   # f(a) = g(b)
        EUFEquation(g(b), h(c)),                   # g(b) = h(c)
        EUFEquation(h(c), f(d)),                   # h(c) = f(d)
        EUFEquation(a, b),                         # a = b
        EUFEquation(b, c),                         # b = c
        EUFEquation(c, d)                          # c = d
    ]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_congruent(g(b), h(c))            # g(a) = h(c)
    assert cc.are_congruent(f(a), h(c))            # f(a) = h(c)
    assert cc.are_congruent(a, d)                  # a = d


def test_compound_expression_propagation():
    # x = y => x*w + z = y*w + z (the caller decided + and * are uninterpreted)
    cc = EUFCongruenceClosure([])
    for t in (add(mul(x, w), z), add(mul(y, w), z)):
        cc._flatten(t)
    cc.merge(x, y)
    assert cc.are_congruent(add(mul(x, w), z), add(mul(y, w), z))


def test_compound_double_layer():
    v = 'v'
    cc = EUFCongruenceClosure([])
    expr1 = add(mul(x, v), z)
    expr2 = add(mul(y, v), w)
    for t in (expr1, expr2):
        cc._flatten(t)
    cc.merge(x, y)
    cc.merge(z, w)
    assert cc.are_congruent(expr1, expr2)


def test_compound_in_function_application():
    # Congruence: x*w + z = y*w + z => f(x*w + z) = f(y*w + z)
    cc = EUFCongruenceClosure([])
    for t in (f(add(mul(x, w), z)), f(add(mul(y, w), z))):
        cc._flatten(t)
    cc.merge(x, y)
    assert cc.are_congruent(f(add(mul(x, w), z)), f(add(mul(y, w), z)))


def test_are_congruent_on_unseen_terms():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    assert not cc.are_congruent(h(z), h(w))
    assert cc.are_congruent(f(a), f(b))
    cc.merge(z, w)
    assert cc.are_congruent(h(z), h(w))


# ---------------------------------------------------------------------------
# explain() -- classical proof-forest explanations (RTA'05 / Inf.Comput.'07).
# Explanations must be a subset of the input equations that alone re-proves
# the queried equality.  They are NOT required to be minimal.
# ---------------------------------------------------------------------------

def _check_explanation(cc, inputs, lhs, rhs):
    """Assert explain(lhs, rhs) is a sound subset-of-inputs explanation."""
    expl = cc.explain(lhs, rhs)
    assert expl is not None
    assert expl <= set(inputs)
    assert EUFCongruenceClosure(list(expl)).are_congruent(lhs, rhs)
    return expl


def test_explain_congruence_edge():
    # x = f(a) = f(b) = y needs the congruence edge f(a)-f(b), which must
    # recurse into the argument proof a = b.
    eqs = [EUFEquation(a, b), EUFEquation(f(a), x), EUFEquation(f(b), y)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, x, y)
    assert expl == set(eqs)


def test_explain_nested_congruence():
    # Two levels of congruence: a = b -> f(a) = f(b) -> g(f(a)) = g(f(b)).
    eqs = [EUFEquation(a, b), EUFEquation(g(f(a)), x), EUFEquation(g(f(b)), y)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, x, y)
    assert EUFEquation(a, b) in expl


def test_explain_ignores_irrelevant_inputs():
    # The z = w component is disjoint and must never leak into explanations.
    eqs = [EUFEquation(a, b), EUFEquation(b, c), EUFEquation(z, w)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, a, c)
    assert expl == {EUFEquation(a, b), EUFEquation(b, c)}


def test_explain_incremental_merges():
    # explain() interleaved with merges: answers must track the growing state.
    cc = EUFCongruenceClosure([])
    cc.merge(a, b)
    assert cc.explain(a, c) is None
    cc.merge(b, c)
    expl = _check_explanation(cc, [EUFEquation(a, b), EUFEquation(b, c)], a, c)
    assert expl == {EUFEquation(a, b), EUFEquation(b, c)}
    # A later query must not be affected by the earlier explain() call
    # (the auxiliary union-find is per-call state).
    assert cc.explain(a, b) == {EUFEquation(a, b)}


def test_explain_may_be_redundant_but_sound():
    # Example 10 of Nieuwenhuis & Oliveras (RTA'05): the proof forest can
    # yield a redundant explanation; it must still be sound and within inputs.
    a1, b1, c1 = 'a1', 'b1', 'c1'
    eqs = [EUFEquation(a1, b1), EUFEquation(a1, c1),
           EUFEquation(f(a1), a), EUFEquation(f(b1), b), EUFEquation(f(c1), c)]
    cc = EUFCongruenceClosure(eqs)
    _check_explanation(cc, eqs, a, c)


def test_explain_reflexive_and_disconnected():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    assert cc.explain(a, a) == set()
    assert cc.explain(f(a), f(a)) == set()
    assert cc.explain(a, z) is None


def test_explain_shortcut_is_level_bounded():
    v = ['n%s' % i for i in range(11)]
    shortcut = EUFEquation(v[0], v[10])
    eqs = [EUFEquation(v[i], v[i + 1]) for i in range(10)] + [shortcut]
    cc = EUFCongruenceClosure(eqs)
    assert _check_explanation(cc, eqs, v[0], v[10]) == set(eqs) - {shortcut}


def test_explain_is_stable_across_calls():
    eqs = [EUFEquation(a, b), EUFEquation(f(a), x), EUFEquation(f(b), y),
           EUFEquation(z, w)]
    cc = EUFCongruenceClosure(eqs)
    first = _check_explanation(cc, eqs, x, y)
    assert cc.explain(x, y) == first
    assert _check_explanation(cc, eqs, y, x) is not None
    assert cc.explain(z, w) == {EUFEquation(z, w)}


def test_explain_after_incremental_merges():
    eqs = [EUFEquation(f(a), x), EUFEquation(f(b), y), EUFEquation(a, b)]
    cc = EUFCongruenceClosure([])
    for eq in eqs:
        cc.merge(*eq)
    assert _check_explanation(cc, eqs, x, y) == set(eqs)


def test_explain_soundness_stress():
    rng = random.Random(20250821)
    s = ['s%s' % i for i in range(12)]
    eqs = []
    for _ in range(18):
        i, j = rng.sample(range(12), 2)
        if rng.random() < 0.3:
            eqs.append(EUFEquation(f(s[i]), f(s[j])))
        else:
            eqs.append(EUFEquation(s[i], s[j]))
    cc = EUFCongruenceClosure(eqs)
    for i in range(12):
        for j in range(i + 1, 12):
            if cc.are_congruent(s[i], s[j]):
                _check_explanation(cc, eqs, s[i], s[j])
            else:
                assert cc.explain(s[i], s[j]) is None


def test_any_head_is_an_uninterpreted_function():
    positive, negative = fn('positive'), fn('negative')
    eqs = [EUFEquation(a, b), EUFEquation(positive(a), x)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_congruent(positive(a), positive(b))
    assert cc.are_congruent(positive(b), x)
    assert not cc.are_congruent(negative(a), x)
    _check_explanation(cc, eqs, positive(b), x)


def test_trivial_and_duplicate_equalities():
    eqs = [EUFEquation(a, a), EUFEquation(a, b), EUFEquation(a, b),
           EUFEquation(b, a)]
    cc = EUFCongruenceClosure(eqs)
    cc.merge(a, a)
    cc.merge(f(a), f(a))
    assert cc.are_congruent(a, b)
    assert cc.explain(a, a) == set()
    assert cc.explain(a, b) <= set(eqs)


def _snapshot(cc):
    def kept(mapping):
        return {k: (set(v) if isinstance(v, set) else list(v))
                for k, v in mapping.items() if v}
    forest = {frozenset((child, parent)): cc.pf_label[child]
              for child, parent in cc.pf_parent.items()}
    return (dict(cc.representative), kept(cc.classlist), dict(cc.lookup_table),
            kept(cc.use_list), forest)


def test_backtrack_retracts_congruence():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    cc.merge(b, c)
    cc.merge(c, d)
    assert cc.are_congruent(a, d)
    cc.backtrack(2)
    assert cc.are_congruent(a, b)
    assert not cc.are_congruent(b, c)
    assert not cc.are_congruent(a, d)
    cc.backtrack(1)
    assert not cc.are_congruent(a, b)


def test_backtrack_retracts_propagated_congruence():
    cc = EUFCongruenceClosure([])
    for t in (g(f(a)), g(f(b))):
        cc._flatten(t)
    cc.merge(a, b)
    assert cc.are_congruent(g(f(a)), g(f(b)))
    cc.backtrack(1)
    assert not cc.are_congruent(f(a), f(b))
    assert not cc.are_congruent(g(f(a)), g(f(b)))


def test_backtrack_keeps_constant_identity():
    cc = EUFCongruenceClosure([])
    consts = {t: cc._flatten(t) for t in (f(a), g(b), h(a, b), x)}
    cc.merge(a, b)
    cc.merge(f(a), x)
    cc.backtrack(2)
    assert {t: cc._flatten(t) for t in consts} == consts
    for const in consts.values():
        assert const in cc.representative


def test_backtrack_is_deterministic():
    cc = EUFCongruenceClosure([])
    for t in (g(f(a)), g(f(b)), f(c), f(d), x, y):
        cc._flatten(t)
    cc.merge(c, d)
    cc.merge(a, b)
    cc.merge(g(f(a)), x)
    cc.backtrack(2)
    once = _snapshot(cc)
    assert cc.are_congruent(c, d)
    cc.merge(a, b)
    cc.merge(g(f(b)), y)
    cc.merge(x, y)
    cc.backtrack(3)
    assert _snapshot(cc) == once


def test_backtrack_counts_every_asserted_equation():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    cc.merge(b, c)
    cc.merge(a, c)
    assert len(cc._asserted) == 3
    cc.backtrack(1)
    assert cc.are_congruent(b, c)
    cc.backtrack(1)
    assert cc.are_congruent(a, b)
    assert not cc.are_congruent(b, c)
    cc.backtrack(1)
    assert not cc.are_congruent(a, b)
    assert len(cc._asserted) == 0


def test_backtrack_argument_validation():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    cc.backtrack(0)
    assert cc.are_congruent(a, b)
    raises(ValueError, lambda: cc.backtrack(2))
    raises(ValueError, lambda: cc.backtrack(-1))
    cc.backtrack(1)
    raises(ValueError, lambda: cc.backtrack(1))


def test_backtrack_handles_applications_registered_afterwards():
    cc = EUFCongruenceClosure([])
    cc.merge(a, b)
    assert cc.are_congruent(f(a), f(b))
    cc.backtrack(1)
    assert not cc.are_congruent(a, b)
    assert not cc.are_congruent(f(a), f(b))
    cc.merge(a, b)
    assert cc.are_congruent(f(a), f(b))


def test_backtrack_handles_application_merged_on_registration():
    cc = EUFCongruenceClosure([])
    cc._flatten(f(a))
    cc.merge(a, b)
    cc.merge(f(b), c)
    cc.backtrack(1)
    assert cc.are_congruent(a, b)
    assert cc.are_congruent(f(a), f(b))
    assert not cc.are_congruent(f(b), c)


def test_backtrack_splits_collapsed_applications():
    cc = EUFCongruenceClosure([])
    cc.merge(a, b)
    assert cc.are_congruent(g(f(a)), g(f(b)))
    cc.backtrack(1)
    assert not cc.are_congruent(a, b)
    assert not cc.are_congruent(g(f(a)), g(f(b)))


def test_backtrack_then_remerge():
    cc = EUFCongruenceClosure([])
    for t in (f(a), f(b), f(c)):
        cc._flatten(t)
    for _ in range(3):
        cc.merge(a, b)
        cc.merge(b, c)
        assert cc.are_congruent(f(a), f(c))
        cc.backtrack(2)
        assert not cc.are_congruent(a, b)
        assert not cc.are_congruent(f(a), f(b))


def test_explain_after_backtrack_uses_only_live_equations():
    v = ['p%s' % i for i in range(8)]
    base = [EUFEquation(v[i], v[i + 1]) for i in range(4)]
    cc = EUFCongruenceClosure(base)
    cc.explain(v[0], v[4])
    for eq in (EUFEquation(v[4], v[5]), EUFEquation(v[0], v[5])):
        cc.merge(*eq)
    cc.explain(v[0], v[5])
    cc.backtrack(2)
    assert not cc.are_congruent(v[0], v[5])
    live = base + [EUFEquation(v[4], v[6])]
    cc.merge(v[4], v[6])
    _check_explanation(cc, live, v[0], v[6])


def test_backtrack_random_differential():
    rng = random.Random(31337)
    s = ['q%s' % i for i in range(9)]
    cc = EUFCongruenceClosure([])
    live = []
    for _ in range(70):
        if live and rng.random() < 0.4:
            k = rng.randint(1, len(live))
            cc.backtrack(k)
            del live[-k:]
        else:
            i, j = rng.sample(range(9), 2)
            lhs, rhs = (s[i], s[j]) if rng.random() < 0.7 else (f(s[i]), f(s[j]))
            if cc.are_congruent(lhs, rhs):
                continue
            cc.merge(lhs, rhs)
            live.append(EUFEquation(lhs, rhs))
        assert len(cc._asserted) == len(live)
        _check_invariants(cc)
        ref = EUFCongruenceClosure(live)
        for i in range(9):
            for j in range(9):
                for term in (lambda u: u, f):
                    assert (cc.are_congruent(term(s[i]), term(s[j]))
                            is ref.are_congruent(term(s[i]), term(s[j])))


# ---------------------------------------------------------------------------
# Textbook problems.  These are the standard congruence closure exercises from
# Bradley & Manna, "The Calculus of Computation" (section 9.2) and the QF_UF
# regressions that ship with SMT solvers; they are small but each one fails
# unless the congruence rule (not just transitivity) is propagated correctly.
# ---------------------------------------------------------------------------

def test_self_applied_function():
    # g(a, b) = a  |-  g(g(a, b), b) = a
    cc = EUFCongruenceClosure([EUFEquation(g(a, b), a)])
    assert cc.are_congruent(g(g(a, b), b), a)
    assert cc.are_congruent(g(g(g(a, b), b), b), a)


def test_iterated_function_cycles():
    # f^3(a) = a together with f^5(a) = a entail f(a) = a.
    powers = [a]
    for _ in range(5):
        powers.append(f(powers[-1]))
    eqs = [EUFEquation(powers[3], a), EUFEquation(powers[5], a)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_congruent(f(a), a)
    assert cc.are_congruent(powers[2], a)
    _check_explanation(cc, eqs, f(a), a)


def test_iterated_function_cycle_is_not_entailed():
    # f^4(a) = a on its own says nothing about f(a).
    powers = [a]
    for _ in range(4):
        powers.append(f(powers[-1]))
    cc = EUFCongruenceClosure([EUFEquation(powers[4], a)])
    assert not cc.are_congruent(f(a), a)
    assert cc.are_congruent(f(powers[4]), f(a))


def test_arity_is_part_of_the_signature():
    # f(a) and f(a, b) share a head but must never be congruent.  Currying
    # makes f(a) the partial application inside f(a, b), which is why the
    # caller has to give each head a single arity.
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    assert not cc.are_congruent(f(a), f(a, b))
    assert not cc.are_congruent(f(a, b), f(b))
    assert cc.are_congruent(f(a, b), f(b, a))


def test_distinct_heads_never_merge():
    cc = EUFCongruenceClosure([EUFEquation(a, b)])
    assert not cc.are_congruent(f(a), g(b))
    assert cc.explain(f(a), g(b)) is None


def test_diamond_explanation_stays_linear():
    # The "diamonds" family used to benchmark proof-producing congruence
    # closure: each diamond offers two routes of equal length, so a proof of
    # the endpoints has to pick one route per diamond rather than explore the
    # 2**n combinations.
    n = 8
    v = ['dm%s' % i for i in range(3*n + 1)]
    eqs = []
    for i in range(n):
        lower, upper = v[3*i], v[3*(i + 1)]
        eqs += [EUFEquation(lower, v[3*i + 1]), EUFEquation(v[3*i + 1], upper),
                EUFEquation(lower, v[3*i + 2]), EUFEquation(v[3*i + 2], upper)]
    cc = EUFCongruenceClosure(eqs)
    assert len(_check_explanation(cc, eqs, v[0], v[3*n])) == 2*n
    _check_invariants(cc)


# ---------------------------------------------------------------------------
# Differential testing against a naive reference, plus internal invariants.
# The reference shares no code with the engine: it just saturates the merge
# rule over every subterm until nothing changes, so agreeing with it is real
# evidence and not the engine confirming itself.  It stays first-order, with
# no currying at all, which is only sound because every head in the pools
# below is used at a single arity.
# ---------------------------------------------------------------------------

def _subterms(term, acc):
    acc.add(term)
    if isinstance(term, EUFApp):
        _subterms(term.func, acc)
        for arg in term.args:
            _subterms(arg, acc)
    return acc


def _reference_congruent(eqs, terms):
    """Return ``congruent(p, q)`` computed by brute-force fixpoint."""
    universe = set()
    for eq in eqs:
        _subterms(eq.lhs, universe)
        _subterms(eq.rhs, universe)
    for term in terms:
        _subterms(term, universe)
    parent = {t: t for t in universe}

    def find(t):
        while parent[t] is not t:
            parent[t] = parent[parent[t]]
            t = parent[t]
        return t

    def union(u, v):
        root_u, root_v = find(u), find(v)
        if root_u is root_v:
            return False
        parent[root_u] = root_v
        return True

    for eq in eqs:
        union(eq.lhs, eq.rhs)
    apps = [t for t in universe if isinstance(t, EUFApp)]
    changed = True
    while changed:
        changed = False
        for i, u in enumerate(apps):
            for v in apps[i + 1:]:
                if u.func != v.func or len(u.args) != len(v.args):
                    continue
                if all(find(p) is find(q) for p, q in zip(u.args, v.args)):
                    changed |= union(u, v)
    return lambda p, q: find(p) is find(q)


def _check_invariants(cc):
    """Assert that the engine's internal state is self-consistent."""
    for const, rep in cc.representative.items():
        assert cc.representative[rep] == rep
        assert const in cc.classlist[rep]
    assert set(cc.classlist) == set(cc.representative.values())
    assert sum(len(m) for m in cc.classlist.values()) == len(cc.representative)

    # every application is still reachable through the signature it has now
    for app, const in cc._app_to_const.items():
        func, (arg,) = app
        key = EUFApp(cc._find_repr(func), (cc._find_repr(arg),))
        assert key in cc.lookup_table
        assert cc._find_repr(cc.lookup_table[key].rhs) == cc._find_repr(const)

    # use_list is filed under live representatives only
    for rep, eqs in cc.use_list.items():
        assert not eqs or cc._find_repr(rep) == rep

    # each class is spanned by an acyclic proof tree of |class| - 1 edges
    edges = defaultdict(int)
    for child, parent in cc.pf_parent.items():
        assert cc._find_repr(child) == cc._find_repr(parent)
        edges[cc._find_repr(child)] += 1
    for rep, members in cc.classlist.items():
        assert edges[rep] == len(members) - 1
    seen = set()
    for node in cc.pf_parent:
        walked = set()
        cursor = node
        while cursor in cc.pf_parent and cursor not in seen:
            assert cursor not in walked
            walked.add(cursor)
            seen.add(cursor)
            cursor = cc.pf_parent[cursor]


def _random_equations(pool, count):
    return [EUFEquation(*sample(pool, 2)) for _ in range(count)]


def _partition(cc, pool):
    """Group the pool of terms by the class each of them lands in."""
    classes = defaultdict(set)
    for term in pool:
        classes[cc._find_repr(cc._flatten(term))].add(term)
    return {frozenset(members) for members in classes.values()}


def test_random_closure_matches_reference():
    s = ['r%s' % i for i in range(6)]
    pool = list(s) + [f(t) for t in s] + [g(t, u) for t, u in zip(s, s[1:])]
    pool += [f(f(s[0])), f(g(s[0], s[1])), h(s[0], s[1], s[2])]
    for _ in range(15):
        eqs = _random_equations(pool, 8)
        cc = EUFCongruenceClosure(eqs)
        congruent = _reference_congruent(eqs, pool)
        for i, p in enumerate(pool):
            for q in pool[i + 1:]:
                assert cc.are_congruent(p, q) is congruent(p, q)
        _check_invariants(cc)


def test_random_explanations_are_sound():
    s = ['e%s' % i for i in range(6)]
    pool = list(s) + [f(t) for t in s] + [g(t, t) for t in s]
    for _ in range(10):
        eqs = _random_equations(pool, 8)
        cc = EUFCongruenceClosure(eqs)
        for i, p in enumerate(pool):
            for q in pool[i + 1:]:
                if cc.are_congruent(p, q):
                    _check_explanation(cc, eqs, p, q)
                else:
                    assert cc.explain(p, q) is None
        _check_invariants(cc)


def test_explain_interleaved_with_merge_and_backtrack():
    # explain() is not read-only: it grows the c-graph with extra edges.  The
    # closure reported afterwards must still match a freshly built engine.
    s = ['i%s' % i for i in range(6)]
    pool = list(s) + [f(t) for t in s] + [g(t, t) for t in s]
    cc = EUFCongruenceClosure([])
    live = []
    for _ in range(20):
        action = choice(['merge', 'merge', 'explain', 'backtrack'])
        if action == 'backtrack' and live:
            k = randint(1, len(live))
            cc.backtrack(k)
            del live[-k:]
        elif action == 'explain':
            p, q = sample(pool, 2)
            expl = cc.explain(p, q)
            if expl is None:
                assert not cc.are_congruent(p, q)
            else:
                assert expl <= set(live)
                assert EUFCongruenceClosure(list(expl)).are_congruent(p, q)
        else:
            p, q = sample(pool, 2)
            cc.merge(p, q)
            live.append(EUFEquation(p, q))
        _check_invariants(cc)
        assert _partition(cc, pool) == _partition(EUFCongruenceClosure(live), pool)


def test_closure_is_order_independent():
    # The closure is a fixpoint, so the order the equations arrive in cannot
    # change which terms end up congruent.
    s = ['o%s' % i for i in range(6)]
    pool = list(s) + [f(t) for t in s] + [g(t, t) for t in s]
    eqs = _random_equations(pool, 9)
    reference = _partition(EUFCongruenceClosure(eqs), pool)
    for _ in range(4):
        permuted = list(eqs)
        shuffle(permuted)
        assert _partition(EUFCongruenceClosure(permuted), pool) == reference
