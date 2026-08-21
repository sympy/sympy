from __future__ import annotations
from sympy.core.symbol import symbols, Symbol, Dummy
from sympy.assumptions.ask import Q
from sympy.core.numbers import Integer
from sympy.core.function import Function, Lambda
from sympy.core.relational import Eq
from sympy.logic.algorithms.euf_theory import (EUFCongruenceClosure,
    EUFUnhandledInput)
from sympy.testing.pytest import raises
import random

f, g, h = symbols('f g h', cls=Function)
x, y, z, w, a, b, c, d = symbols('x y z w a b c d')


def test_basic_and_chain_equality():
    cc = EUFCongruenceClosure([Q.eq(x, y), Q.eq(y, z)])
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


def test_lambda_curry_and_equivalent_application():
    lam = Lambda((x, y), x + 2*y)
    cc = EUFCongruenceClosure([
        Q.eq(lam(x,y), lam(y,x)),
        Q.eq(x, y)
    ])
    assert cc.are_congruent(lam(x,y), lam(y,x))


def test_permuted_arguments_no_commutativity():
    lam_h = Lambda((x, y), h(x, y))

    cc = EUFCongruenceClosure([
        Q.eq(lam_h(x, y), lam_h(y, x)),     # h(x,y) = h(y,x)
        Q.eq(x, y)                          # x = y
    ])
    # Even without commutativity, if x=y, h(x,y)=h(y,x) by congruence
    assert cc.are_congruent(lam_h(x, y), lam_h(y, x))  # h(x,y) = h(y,x)


def test_nested_lambdas_chain():
    lam = Lambda((x, y), x + y).curry()
    cc = EUFCongruenceClosure([
        Q.eq(lam(x)(y), lam(y)(x)),
        Q.eq(x, y)
    ])
    assert cc.are_congruent(lam(x)(y), lam(y)(x))


def test_add_equality_registers_and_merges():
    cc = EUFCongruenceClosure([])

    a, b = Symbol('a'), Symbol('b')
    f = Function('f')

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))

    cc.merge(a, b)
    assert cc._find_repr(a) == cc._find_repr(b)
    assert cc._find_repr(fa) == cc._find_repr(fb)


def test_mixed_lambdas_flatten_unique():
    cc = EUFCongruenceClosure([])
    lam1 = Lambda(x, x + 1).curry()
    lam2 = Lambda(x, x + 2).curry()
    d1 = cc._flatten(lam1(x))
    d2 = cc._flatten(lam2(x))
    assert d1 != d2


def test_flatten_application_and_cache():
    cc = EUFCongruenceClosure([])
    testf = Function('testf')
    ax = cc._flatten(testf(x))
    bx = cc._flatten(testf(x))
    assert ax == bx


def test_find_of_transformed_symbols():
    cc = EUFCongruenceClosure([])
    t = Symbol('t')
    assert cc._find_repr(cc._flatten(t)) == t
    t2 = Dummy('t2')
    assert cc._find_repr(cc._flatten(t2)) == t2


def test_process_pending_chain_merges():
    cc = EUFCongruenceClosure([])
    f1 = Function('alsof')
    x1, y1, z1 = symbols('x1 y1 z1')
    cc._register(x1)
    cc._register(y1)
    cc._register(z1)
    fx = cc._flatten(f1(x1))
    fy = cc._flatten(f1(y1))
    fz = cc._flatten(f1(z1))
    cc.use_list[x1].append((f1, (x1,), fx))
    cc.use_list[y1].append((f1, (y1,), fy))
    cc.use_list[z1].append((f1, (z1,), fz))
    cc.lookup_table[(f1, (x1,))] = (f1, (x1,), fx)
    cc.lookup_table[(f1, (y1,))] = (f1, (y1,), fy)
    cc.lookup_table[(f1, (z1,))] = (f1, (z1,), fz)
    cc.pending.append((x1, y1))
    cc.pending.append((y1, z1))
    cc.pending.append((fx, fy))
    cc.pending.append((fy, fz))
    cc._process_pending_unions()
    assert cc._find_repr(x1) == cc._find_repr(y1) == cc._find_repr(z1)
    assert cc._find_repr(fx) == cc._find_repr(fy) == cc._find_repr(fz)


def test_flatten_lambda_consistency_and_cache():
    cc = EUFCongruenceClosure([])
    t = Symbol('t')
    f1 = Function('ccf')
    lam = Lambda(t, f1(t))
    flat1 = cc._flatten(lam)
    flat2 = cc._flatten(lam)
    assert flat1 == flat2
    flat_app1 = cc._flatten(lam(t))
    flat_app2 = cc._flatten(lam(t))
    assert flat_app1 == flat_app2


def test_use_list_merging_under_union():
    cc = EUFCongruenceClosure([])
    a1, b1, c1 = Symbol('a1'), Symbol('b1'), Symbol('c1')
    f1 = Function('f1')
    # Register variables and applications before any merging
    apps = [cc._flatten(f1(t)) for t in (a1, b1, c1)]
    # Union a1 and b1
    cc._union(a1, b1)
    # Now add equality b1 = c1, so all three are merged
    cc._union(b1, c1)
    cc._process_pending_unions()
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
    vars = symbols('a0:20')
    eqs = [Q.eq(vars[i], vars[i+1]) for i in range(len(vars)-1)]
    cc = EUFCongruenceClosure(eqs)
    for i in range(len(vars)):
        for j in range(len(vars)):
            assert cc.are_congruent(vars[i], vars[j])


def test_composed_functions():
    f, g, h = symbols('f g h', cls=Function)
    a, b, c = symbols('a b c')
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
    lam_f = Lambda(symbols('x'), f('x'))
    lam_g = Lambda(symbols('x'), g('x'))
    lam_h = Lambda(x, Lambda(y, h(x, y)))
    eq1 = Q.eq(lam_f(a), lam_g(b))                   # f(a) = g(b)
    eq2 = Q.eq(lam_g(c), lam_h(lam_f(c))(lam_g(a)))  # g(c) = h(f(c), g(a))
    eq3 = Q.eq(b, c)                                 # b = c
    eq4 = Q.eq(lam_f(c), lam_g(a))                   # f(c) = g(a)
    eq5 = Q.eq(lam_h(d)(d), lam_g(b))                # h(d, d) = g(b)
    eq6 = Q.eq(lam_g(a), d)                          # g(a) = d

    cc = EUFCongruenceClosure([eq1, eq2, eq3, eq4, eq5, eq6])

    # Assertions checking congruence closure identifies equalities properly
    assert cc.are_congruent(b, c)                      # b = c
    assert cc.are_congruent(lam_g(a), d)               # g(a) = d
    assert cc.are_congruent(lam_g(b), lam_g(c))        # g(b) = g(c)
    assert cc.are_congruent(lam_f(a), lam_g(c))        # f(a) = g(c)


def test_example_2():
    lam_f = Lambda(x, f(x))
    lam_g = Lambda(x, g(x))
    lam_h = Lambda(x, h(x))
    eqs = [
        Q.eq(lam_f(a), lam_g(b)),                    # f(a) = g(b)
        Q.eq(lam_g(b), lam_h(c)),                    # g(b) = h(c)
        Q.eq(lam_h(c), lam_f(d)),                    # h(c) = f(d)
        Q.eq(a, b),                                  # a = b
        Q.eq(b, c),                                  # b = c
        Q.eq(c, d)                                   # c = d
    ]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_congruent(lam_g(b), lam_h(c))        # g(a) = h(c)
    assert cc.are_congruent(lam_f(a), lam_h(c))        # f(a) = h(c)
    assert cc.are_congruent(a, d)                      # a = d


def test_flatten_simple_atoms_and_numbers():
    cc = EUFCongruenceClosure([])

    a = Symbol('a')
    d = Dummy('d')
    n1 = Integer(1)
    n2 = Integer(2)

    # Symbols flatten to themselves
    assert cc._flatten(a) == a

    # Dummy flatten to themselves
    assert cc._flatten(d) == d

    # Same number always maps to same Dummy
    flat_n1_1 = cc._flatten(n1)
    flat_n1_2 = cc._flatten(Integer(1))
    flat_n2 = cc._flatten(n2)
    assert flat_n1_1 == flat_n1_2
    assert flat_n1_1 != flat_n2

    # Different atoms (non-symbol) map to distinct Dummies
    from sympy.logic.boolalg import BooleanTrue
    btrue = BooleanTrue()
    flat_btrue1 = cc._flatten(btrue)
    flat_btrue2 = cc._flatten(btrue)
    assert flat_btrue1 == flat_btrue2
    assert flat_btrue1 != d  # Different from other dummies


def test_compound_expression_propagation():
    x, y, z, w = symbols('x y z w')
    # x = y => x*w + z = y*w + z (Mul/Add treated as uninterpreted)
    cc = EUFCongruenceClosure([])
    for t in (x*w + z, y*w + z):
        cc._flatten(t)
    cc.merge(x, y)
    assert cc.are_congruent(x*w + z, y*w + z)


def test_compound_double_layer():
    x, y, z, w, v = symbols('x y z w v')
    cc = EUFCongruenceClosure([])
    expr1 = x*v + z
    expr2 = y*v + w
    for t in (expr1, expr2):
        cc._flatten(t)
    cc.merge(x, y)
    cc.merge(z, w)
    assert cc.are_congruent(expr1, expr2)


def test_mixed_equality_disequality():
    x, y, z = symbols('x y z')
    cc = EUFCongruenceClosure([])
    # x = y, so x*z + y = y*z + y (uninterpreted congruence)
    for t in (x*z + y, y*z + y):
        cc._flatten(t)
    cc.merge(x, y)
    assert cc.are_congruent(x*z + y, y*z + y)


def test_compound_in_function_application():
    from sympy import Function
    x, y, z, w = symbols('x y z w')
    f = Function('f')
    cc = EUFCongruenceClosure([])
    # Congruence: x*w + z = y*w + z => f(x*w + z) = f(y*w + z)
    for t in (f(x*w + z), f(y*w + z)):
        cc._flatten(t)
    cc.merge(x, y)
    assert cc.are_congruent(f(x*w + z), f(y*w + z))


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
    eqs = [Q.eq(a, b), Q.eq(f(a), x), Q.eq(f(b), y)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, x, y)
    assert expl == set(eqs)


def test_explain_nested_congruence():
    # Two levels of congruence: a = b -> f(a) = f(b) -> g(f(a)) = g(f(b)).
    eqs = [Q.eq(a, b), Q.eq(g(f(a)), x), Q.eq(g(f(b)), y)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, x, y)
    assert Q.eq(a, b) in expl


def test_explain_ignores_irrelevant_inputs():
    # The z = w component is disjoint and must never leak into explanations.
    eqs = [Q.eq(a, b), Q.eq(b, c), Q.eq(z, w)]
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, a, c)
    assert expl == {Q.eq(a, b), Q.eq(b, c)}


def test_explain_incremental_merges():
    # explain() interleaved with merges: answers must track the growing state.
    cc = EUFCongruenceClosure([])
    cc.merge(a, b)
    assert cc.explain(a, c) is None
    cc.merge(b, c)
    expl = _check_explanation(cc, [Q.eq(a, b), Q.eq(b, c)], a, c)
    assert expl == {Q.eq(a, b), Q.eq(b, c)}
    # A later query must not be affected by the earlier explain() call
    # (the auxiliary union-find is per-call state).
    assert cc.explain(a, b) == {Q.eq(a, b)}


def test_explain_may_be_redundant_but_sound():
    # Example 10 of Nieuwenhuis & Oliveras (RTA'05): the proof forest can
    # yield a redundant explanation; it must still be sound and within inputs.
    a1, b1, c1 = symbols('a1 b1 c1')
    eqs = [Q.eq(a1, b1), Q.eq(a1, c1),
           Q.eq(f(a1), a), Q.eq(f(b1), b), Q.eq(f(c1), c)]
    cc = EUFCongruenceClosure(eqs)
    _check_explanation(cc, eqs, a, c)


# ---------------------------------------------------------------------------
# Input validation and queries on unseen terms.
# ---------------------------------------------------------------------------

def test_rejects_non_equality_input():
    raises(EUFUnhandledInput, lambda: EUFCongruenceClosure([Q.positive(x)]))
    raises(EUFUnhandledInput, lambda: EUFCongruenceClosure([Q.ne(x, y)]))
    raises(EUFUnhandledInput, lambda: EUFCongruenceClosure([Eq(x, y)]))
    # A bad equation anywhere in the list is rejected, not silently dropped.
    raises(EUFUnhandledInput,
           lambda: EUFCongruenceClosure([Q.eq(x, y), Q.positive(z)]))


def test_are_congruent_on_unseen_terms():
    # Querying terms the engine has never flattened must answer, not raise.
    cc = EUFCongruenceClosure([Q.eq(a, b)])
    assert not cc.are_congruent(h(z), h(w))
    assert cc.are_congruent(f(a), f(b))     # congruence found on first query
    cc.merge(z, w)
    assert cc.are_congruent(h(z), h(w))     # closure catches up after merge


def test_get_canonical_form():
    cc = EUFCongruenceClosure([])
    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))
    # A constant that replaced no application is its own canonical form.
    assert cc._get_canonical_form(a) == a
    assert cc._get_canonical_form(fa) != cc._get_canonical_form(fb)
    cc.merge(a, b)
    # Arguments are rewritten to representatives, so both apps now agree.
    assert cc._get_canonical_form(fa) == cc._get_canonical_form(fb)


# ---------------------------------------------------------------------------
# explain() -- greedy (c-graph) path.
# ---------------------------------------------------------------------------

def test_explain_reflexive_and_disconnected():
    cc = EUFCongruenceClosure([Q.eq(a, b)])
    assert cc.explain(a, a) == set()
    assert cc.explain(f(a), f(a)) == set()
    assert cc.explain(a, z) is None          # z is unknown, hence not equal


def test_explain_prefers_shortcut_over_chain():
    # e0 = e10 is redundant when it arrives, so the proof forest drops it
    # while the c-graph keeps it.  Explaining e0 = e11 must take that
    # shortcut instead of walking the ten chain edges.
    v = symbols('e0:12')
    shortcut = Q.eq(v[0], v[10])
    eqs = ([Q.eq(v[i], v[i + 1]) for i in range(10)]
           + [shortcut, Q.eq(v[10], v[11])])
    cc = EUFCongruenceClosure(eqs)
    expl = _check_explanation(cc, eqs, v[0], v[11])
    assert expl == {shortcut, Q.eq(v[10], v[11])}


def test_explain_shortcut_is_level_bounded():
    # The same shortcut may not be used for e0 = e10 itself: at the moment
    # those two became equal the edge did not exist yet, so the greedy search
    # must ignore it and still return a sound explanation.
    v = symbols('n0:11')
    shortcut = Q.eq(v[0], v[10])
    eqs = [Q.eq(v[i], v[i + 1]) for i in range(10)] + [shortcut]
    cc = EUFCongruenceClosure(eqs)
    assert _check_explanation(cc, eqs, v[0], v[10]) == set(eqs) - {shortcut}


def test_explain_zero_fuel_falls_back_to_classical():
    # With no greedy fuel every congruence edge is expanded classically; the
    # answer must stay sound.
    eqs = [Q.eq(a, b), Q.eq(g(f(a)), x), Q.eq(g(f(b)), y)]
    cc = EUFCongruenceClosure(eqs)
    cc.greedy_fuel = 0
    assert _check_explanation(cc, eqs, x, y) == set(eqs)


def test_explain_is_stable_across_calls():
    # explain() mutates c-graph state (extra edges, seen pairs); repeating a
    # query must not change the answer or corrupt later ones.
    eqs = [Q.eq(a, b), Q.eq(f(a), x), Q.eq(f(b), y), Q.eq(z, w)]
    cc = EUFCongruenceClosure(eqs)
    first = _check_explanation(cc, eqs, x, y)
    assert cc.explain(x, y) == first
    assert _check_explanation(cc, eqs, y, x) is not None
    assert cc.explain(z, w) == {Q.eq(z, w)}


def test_extra_edges_stay_within_budget():
    v = symbols('m0:10')
    eqs = ([Q.eq(v[i], v[i + 1]) for i in range(9)]
           + [Q.eq(f(v[i]), g(v[i])) for i in range(10)])
    cc = EUFCongruenceClosure(eqs)
    for i in range(10):
        cc.explain(f(v[0]), f(v[i]))
    assert cc._n_edges_extra <= 2 * cc._n_edges_during_union


def test_explain_after_incremental_merges():
    # Labels produced by merge() must be usable as explanations too.
    eqs = [Q.eq(f(a), x), Q.eq(f(b), y), Q.eq(a, b)]
    cc = EUFCongruenceClosure([])
    for eq in eqs:
        cc.merge(eq.lhs, eq.rhs)
    assert _check_explanation(cc, eqs, x, y) == set(eqs)


def test_explain_soundness_stress():
    # Random equality graph: every derivable pair must get an explanation
    # that is a subset of the inputs and re-proves the equality on its own.
    rng = random.Random(20250821)
    s = symbols('s0:12')
    eqs = []
    for _ in range(18):
        i, j = rng.sample(range(12), 2)
        if rng.random() < 0.3:
            eqs.append(Q.eq(f(s[i]), f(s[j])))
        else:
            eqs.append(Q.eq(s[i], s[j]))
    cc = EUFCongruenceClosure(eqs)
    for i in range(12):
        for j in range(i + 1, 12):
            if cc.are_congruent(s[i], s[j]):
                _check_explanation(cc, eqs, s[i], s[j])
            else:
                assert cc.explain(s[i], s[j]) is None


def test_predicate_terms_are_uninterpreted_functions():
    # An AppliedPredicate appearing as a *term* is flattened like any other
    # application, so congruence applies to it.
    eqs = [Q.eq(a, b), Q.eq(Q.positive(a), x)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_congruent(Q.positive(a), Q.positive(b))
    assert cc.are_congruent(Q.positive(b), x)
    assert not cc.are_congruent(Q.negative(a), x)
    _check_explanation(cc, eqs, Q.positive(b), x)


def test_trivial_and_duplicate_equalities():
    # Self-equalities and repeats must be absorbed without corrupting state.
    eqs = [Q.eq(a, a), Q.eq(a, b), Q.eq(a, b), Q.eq(b, a)]
    cc = EUFCongruenceClosure(eqs)
    cc.merge(a, a)
    cc.merge(f(a), f(a))
    assert cc.are_congruent(a, b)
    assert cc.explain(a, a) == set()
    assert cc.explain(a, b) <= set(eqs)


def test_greedy_explanation_never_larger_than_classical():
    # explain() is not required to be minimal, but the greedy search must
    # never do worse than the plain proof-forest walk it replaces ([2] S.3).
    rng = random.Random(26)
    s = symbols('t0:12')
    eqs = []
    for _ in range(18):
        i, j = rng.sample(range(12), 2)
        if rng.random() < 0.3:
            eqs.append(Q.eq(f(s[i]), f(s[j])))
        else:
            eqs.append(Q.eq(s[i], s[j]))
    greedy = EUFCongruenceClosure(eqs)
    # A second engine, so the greedy c-graph state cannot skew the baseline.
    classical = EUFCongruenceClosure(eqs)
    shortened = 0
    for i in range(12):
        for j in range(i + 1, 12):
            if not greedy.are_congruent(s[i], s[j]):
                continue
            baseline = set()
            classical._explain_classical(classical._flatten(s[i]),
                                         classical._flatten(s[j]), baseline)
            assert EUFCongruenceClosure(list(baseline)).are_congruent(s[i], s[j])
            expl = _check_explanation(greedy, eqs, s[i], s[j])
            assert len(expl) <= len(baseline)
            shortened += len(expl) < len(baseline)
    # The fixture must actually exercise the shortening, not just tie.
    assert shortened > 0
