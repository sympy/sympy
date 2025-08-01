from sympy import Symbol, Dummy, Eq, Lambda, Integer
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
from sympy.abc import x, y, z, w, a, b, c, d
from sympy.core.function import Function
from sympy import symbols

f, g, h = symbols('f g h', cls=Function)

def test_basic_equality():
    cc = EUFCongruenceClosure([Eq(x, y), Eq(y, z)])
    assert cc.are_equal(x, y)
    assert cc.are_equal(y, z)
    assert cc.are_equal(x, z)
    assert not cc.are_equal(x, w)


def test_curried_lambda_congruence():
    lam_curried = Lambda((x, y), x + y).curry()
    assert lam_curried == Lambda(x, Lambda(y, x + y))
    cc = EUFCongruenceClosure([
        Eq(lam_curried(x)(y), lam_curried(y)(x)),
        Eq(x, y)
    ])
    # Since x=y, lam_curried(x)(y) == lam_curried(y)(x)
    a = lam_curried(x)(y)
    b = lam_curried(y)(x)
    assert cc.are_equal(a, b)


def test_nested_lambdas():
    lam = Lambda((x, y), x + y).curry()
    inner = lam(x)
    outer = inner(y)
    cc = EUFCongruenceClosure([
        Eq(lam(x)(y), lam(y)(x)),
        Eq(x, y),
        Eq(lam(x), lam(y))
    ])
    assert cc.are_equal(outer, lam(y)(x))
    assert cc.are_equal(lam(x), lam(y))


def test_inequality_of_different_functions():
    lam1 = Lambda(x, x + 1).curry()
    lam2 = Lambda(x, x + 2).curry()
    cc = EUFCongruenceClosure([
        Eq(lam1(x), Symbol('a')),
        Eq(lam2(x), Symbol('b'))
    ])
    assert not cc.are_equal(lam1(x), lam2(x))
    assert not cc.are_equal(Symbol('a'), Symbol('b'))


def test_identity_lambda_equality():
    lam = Lambda(x, x).curry()
    cc = EUFCongruenceClosure([Eq(lam, lam)])
    assert cc.are_equal(lam, lam)


def test_constants_and_dummies():
    a = Symbol('a')
    d = Dummy('d')
    cc = EUFCongruenceClosure([Eq(a, d)])
    assert cc.are_equal(a, d)


def test_trivial_symbol_equality():
    cc = EUFCongruenceClosure([Eq(x, y)])
    assert cc.are_equal(x, y)
    cc = EUFCongruenceClosure([Eq(x, x)])
    assert cc.are_equal(x, x)
    assert not cc.are_equal(x, y)


def test_unconnected_terms():
    cc = EUFCongruenceClosure([Eq(x, y)])
    assert not cc.are_equal(x, z)
    assert not cc.are_equal(y, z)


def test_example_1():
    lam_f = Lambda(symbols('x'), f('x'))
    lam_g = Lambda(symbols('x'), g('x'))
    lam_h = Lambda(x, Lambda(y, h(x, y)))
    eq1 = Eq(lam_f(a), lam_g(b))                   # f(a) = g(b)
    eq2 = Eq(lam_g(c), lam_h(lam_f(c))(lam_g(a)))  # g(c) = h(f(c), g(a))
    eq3 = Eq(b, c)                                 # b = c
    eq4 = Eq(lam_f(c), lam_g(a))                   # f(c) = g(a)
    eq5 = Eq(lam_h(d)(d), lam_g(b))                # h(d, d) = g(b)
    eq6 = Eq(lam_g(a), d)                          # g(a) = d

    cc = EUFCongruenceClosure([eq1, eq2, eq3, eq4, eq5, eq6])

    # Assertions checking congruence closure identifies equalities properly
    assert cc.are_equal(b, c)                      # b = c
    assert cc.are_equal(lam_g(a), d)               # g(a) = d
    assert cc.are_equal(lam_g(b), lam_g(c))        # g(b) = g(c)
    assert cc.are_equal(lam_f(a), lam_g(c))        # f(a) = g(c)


def test_example_2():
    lam_f = Lambda(x, f(x))
    lam_g = Lambda(x, g(x))
    lam_h = Lambda(x, h(x))
    eqs = [
        Eq(lam_f(a), lam_g(b)),                    # f(a) = g(b)
        Eq(lam_g(b), lam_h(c)),                    # g(b) = h(c)
        Eq(lam_h(c), lam_f(d)),                    # h(c) = f(d)
        Eq(a, b),                                  # a = b
        Eq(b, c),                                  # b = c
        Eq(c, d)                                   # c = d
    ]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(lam_g(b), lam_h(c))        # g(a) = h(c)
    assert cc.are_equal(lam_f(a), lam_h(c))        # f(a) = h(c)
    assert cc.are_equal(a, d)                      # a = d

def test_nested_applications_chaining():
    lam_f = Lambda(x, f(x))
    lam_g = Lambda(x, g(x))
    # Equate inner and outer applications via nested equalities
    eqs = [
        Eq(lam_f(lam_g(x)), lam_g(lam_f(x))),      # f(g(x)) = g(f(x))
        Eq(lam_g(y), lam_f(z)),                    # g(y) = f(z)
        Eq(x, y),                                  # x = y
        Eq(y, z)                                   # y = z
    ]
    cc = EUFCongruenceClosure(eqs)
    # All nested applications and all input vars should collapse together
    assert cc.are_equal(lam_g(y), lam_f(x))                  # g(y) = f(x)
    assert cc.are_equal(lam_f(x), lam_f(z))                  # f(x) = f(z)



def test_permuted_arguments_no_commutativity():
    lam_h = Lambda((x, y), h(x, y))

    cc = EUFCongruenceClosure([
        Eq(lam_h(x, y), lam_h(y, x)),     # h(x,y) = h(y,x)
        Eq(x, y)                          # x = y
    ])
    # Even without commutativity, if x=y, h(x,y)=h(y,x) by congruence
    assert cc.are_equal(lam_h(x, y), lam_h(y, x))  # h(x,y) = h(y,x)


def test_register_dummy_or_symbol_and_find_basic():
    cc = EUFCongruenceClosure([])

    a = Symbol('a')
    d = Dummy('d')

    # Register symbol and dummy
    cc._register_dummy_or_symbol(a)
    cc._register_dummy_or_symbol(d)

    # Find representatives immediately after registration
    assert cc._find(a) == a
    assert cc._find(d) == d

    # Re-registering does not cause error and maintains state
    cc._register_dummy_or_symbol(a)
    cc._register_dummy_or_symbol(d)


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


def test_flatten_lambdas_and_application_dummies():
    cc = EUFCongruenceClosure([])

    f = Function('f')
    x, y = Symbol('x'), Symbol('y')

    # Lambda with multiple variables gets curried and flattened
    lam = Lambda((x, y), x + y)
    lam_curried = lam.curry()
    flat_lam = cc._flatten(lam)
    flat_curried = cc._flatten(lam_curried)
    assert flat_lam == flat_curried  # Same flattening

    # Different lambdas with different bodies are different
    lam2 = Lambda(x, x*y)
    flat_lam2 = cc._flatten(lam2)
    assert flat_lam != flat_lam2


def test_find_with_unregistered_element_creates_singleton():
    cc = EUFCongruenceClosure([])

    s = Symbol('s')
    # s not registered yet
    rep1 = cc._find(s)
    rep2 = cc._find(s)
    assert rep1 == rep2 == s
    # Also registered now
    assert s in cc.rep
    assert s in cc.cls


def test_union_merges_classes_and_propagates_uses():
    cc = EUFCongruenceClosure([])

    a, b = Symbol('a'), Symbol('b')
    f = Function('f')

    cc._register_dummy_or_symbol(a)
    cc._register_dummy_or_symbol(b)

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))

    # Initially separate classes
    assert cc._find(a) != cc._find(b)
    assert cc._find(fa) != cc._find(fb)

    # Add UseLists manually to simulate function usage
    cc.use[a].append((f, (a,), fa))
    cc.use[b].append((f, (b,), fb))
    cc.lookup[(f, (a,))] = fa
    cc.lookup[(f, (b,))] = fb

    # Union a and b
    cc._union(a, b)

    # Then representatives should merge
    assert cc._find(a) == cc._find(b)


def test_process_pending_merges_until_fixpoint():
    cc = EUFCongruenceClosure([])

    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    f = Function('f')

    cc._register_dummy_or_symbol(x)
    cc._register_dummy_or_symbol(y)
    cc._register_dummy_or_symbol(z)

    fx = cc._flatten(f(x))
    fy = cc._flatten(f(y))
    fz = cc._flatten(f(z))

    # Initially disjoint
    assert cc._find(x) != cc._find(y) != cc._find(z)
    assert cc._find(fx) != cc._find(fy) != cc._find(fz)

    cc.use[x].append((f, (x,), fx))
    cc.use[y].append((f, (y,), fy))
    cc.use[z].append((f, (z,), fz))

    cc.lookup[(f, (x,))] = fx
    cc.lookup[(f, (y,))] = fy
    cc.lookup[(f, (z,))] = fz

    # Add equalities to pending, simulating chained merges
    cc.pending.append((x, y))
    cc.pending.append((y, z))
    cc.pending.append((fx, fy))
    cc.pending.append((fy, fz))

    cc._process_pending()

    # All merged now
    assert cc._find(x) == cc._find(y) == cc._find(z)
    assert cc._find(fx) == cc._find(fy) == cc._find(fz)


def test_add_equality_registers_and_merges():
    cc = EUFCongruenceClosure([])

    a, b = Symbol('a'), Symbol('b')
    f = Function('f')

    cc.add_equality(a, b)
    assert cc._find(a) == cc._find(b)

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))

    # Before adding equality of f(a), f(b) classes differ
    assert cc._find(fa) != cc._find(fb)


def test_flatten_cache_consistency():
    cc = EUFCongruenceClosure([])

    f = Function('f')
    x = Symbol('x')

    lam = Lambda(x, f(x))

    flat1 = cc._flatten(lam)
    flat2 = cc._flatten(lam)
    assert flat1 == flat2  # Same Flatten repeated

    flat_app1 = cc._flatten(lam(x))
    flat_app2 = cc._flatten(lam(x))
    assert flat_app1 == flat_app2  # Cached application flatten results


def test_union_merge_use_lists_correctness():
    cc = EUFCongruenceClosure([])

    a, b, c = Symbol('a'), Symbol('b'), Symbol('c')
    f = Function('f')

    cc._register_dummy_or_symbol(a)
    cc._register_dummy_or_symbol(b)
    cc._register_dummy_or_symbol(c)

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))
    fc = cc._flatten(f(c))

    # Setup UseList for a and b with different function applications
    cc.use[a].append((f, (a,), fa))
    cc.use[b].append((f, (b,), fb))
    cc.use[c].append((f, (c,), fc))

    cc.lookup[(f, (a,))] = fa
    cc.lookup[(f, (b,))] = fb
    cc.lookup[(f, (c,))] = fc

    # Union a, b; Use lists merge correctly, propagation triggers equality of f(a), f(b)
    cc._union(a, b)

    # After union, use list for the merged rep includes uses from both
    rep = cc._find(a)
    assert len(cc.use[rep]) >= 2