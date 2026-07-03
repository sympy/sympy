from __future__ import annotations
from sympy.core.symbol import symbols, Symbol, Dummy
from sympy.assumptions.ask import Q
from sympy.core.numbers import Integer
from sympy.core.function import Function, Lambda
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure

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
    cc.lookup_table[(f1, (x1,))] = fx
    cc.lookup_table[(f1, (y1,))] = fy
    cc.lookup_table[(f1, (z1,))] = fz
    cc.pending_unions.append((x1, y1))
    cc.pending_unions.append((y1, z1))
    cc.pending_unions.append((fx, fy))
    cc.pending_unions.append((fy, fz))
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
