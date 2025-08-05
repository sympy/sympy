from sympy.core.symbol import symbols, Symbol, Dummy
from sympy.core.relational import Eq
from sympy.core.numbers import Integer
from sympy.core.function import Function, Lambda
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
import random

f, g, h = symbols('f g h', cls=Function)
x, y, z, w, a, b, c, d = symbols('x y z w a b c d')


def test_basic_and_chain_equality():
    cc = EUFCongruenceClosure([Eq(x, y), Eq(y, z)])
    assert cc.are_equal(x, y)
    assert cc.are_equal(y, z)
    assert cc.are_equal(x, z)
    assert not cc.are_equal(x, w)


def test_unary_function_congruence():
    cc = EUFCongruenceClosure([Eq(a, b), Eq(f(a), x)])
    assert cc.are_equal(f(b), x)   # f(a) = x, a=b -> f(b) = x


def test_binary_congruence_and_propagation():
    cc = EUFCongruenceClosure([
        Eq(a, b),
        Eq(c, d),
        Eq(g(a, c), x)
    ])
    assert cc.are_equal(g(b, d), x)  # g(a,c) = x; a=b, c=d -> g(b,d)=x


def test_lambda_curry_and_equivalent_application():
    lam = Lambda((x, y), x + 2*y)
    lam_curried = lam.curry()
    cc = EUFCongruenceClosure([
        Eq(lam_curried(x)(y), lam_curried(y)(x)),
        Eq(x, y)
    ])
    assert cc.are_equal(lam_curried(x)(y), lam_curried(y)(x))


def test_permuted_arguments_no_commutativity():
    lam_h = Lambda((x, y), h(x, y))

    cc = EUFCongruenceClosure([
        Eq(lam_h(x, y), lam_h(y, x)),     # h(x,y) = h(y,x)
        Eq(x, y)                          # x = y
    ])
    # Even without commutativity, if x=y, h(x,y)=h(y,x) by congruence
    assert cc.are_equal(lam_h(x, y), lam_h(y, x))  # h(x,y) = h(y,x)


def test_nested_lambdas_chain():
    lam = Lambda((x, y), x + y).curry()
    cc = EUFCongruenceClosure([
        Eq(lam(x)(y), lam(y)(x)),
        Eq(x, y)
    ])
    assert cc.are_equal(lam(x)(y), lam(y)(x))


def test_add_equality_registers_and_merges():
    cc = EUFCongruenceClosure([])

    a, b = Symbol('a'), Symbol('b')
    f = Function('f')

    cc.add_equality(a, b)
    assert cc._find(a) == cc._find(b)

    fa = cc._flatten(f(a))
    fb = cc._flatten(f(b))

    assert cc._find(fa) == cc._find(fb)


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


def test_find_with_automatic_registration():
    cc = EUFCongruenceClosure([])
    t = Symbol('t')
    assert cc._find(t) == t
    t2 = Dummy('t2')
    assert cc._find(t2) == t2


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
    assert cc._find(x1) == cc._find(y1) == cc._find(z1)
    assert cc._find(fx) == cc._find(fy) == cc._find(fz)


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
    # Register variables and applications
    cc._register(a1)
    cc._register(b1)
    cc._register(c1)
    # Flatten applications
    fa = cc._flatten(f1(a1))
    fb = cc._flatten(f1(b1))
    fc = cc._flatten(f1(c1))
    # Union a1 and b1
    cc._union(a1, b1)
    # Now add equality b1 = c1, so all three are merged
    cc._union(b1, c1)
    rep = cc._find(a1)
    class_members = cc.classlist[rep]
    # Test: all applications f1(x) for all class members x are congruent
    app_reps = set(cc._find(cc._flatten(f1(x))) for x in class_members)
    assert len(app_reps) == 1


def test_complex_deep_chaining():
    # Very deep nesting of f
    depth = 190
    term_a = a
    term_b = b
    for _ in range(depth):
        term_a = f(term_a)
        term_b = f(term_b)
    eqs = [Eq(term_a, x), Eq(a, b)]
    cc = EUFCongruenceClosure(eqs)


    # All nestings over a and b should be equal to each other and to x
    for _ in range(depth):
        assert cc.are_equal(term_a, term_b)
        term_a = f(term_a)
        term_b = f(term_b)


def random_expr(vars, funcs, depth=3):
    """
    Recursively builds a random EUF expression (nested function or variable/integer).
    """
    if depth == 0 or (depth < 3 and random.random() < 0.5):
        # base case: variable or integer
        if random.random() < 0.7:
            return random.choice(vars)
        else:
            return Integer(random.randint(0, 10))
    else:
        func = random.choice(funcs)
        # unary for simplicity, but you can generalize to n-ary
        return func(random_expr(vars, funcs, depth-1))

def test_random_euf_equations(num_eqs=5, max_depth=4):
    vars = symbols('x y z')
    funcs = [Function(f'f{i}') for i in range(3)]
    eqs = []
    for _ in range(num_eqs):
        e1 = random_expr(vars, funcs, max_depth)
        e2 = random_expr(vars, funcs, max_depth)
        eqs.append(Eq(e1, e2))
    return eqs


def test_long_chain_variables():
    vars = symbols('a0:20')
    eqs = [Eq(vars[i], vars[i+1]) for i in range(len(vars)-1)]
    cc = EUFCongruenceClosure(eqs)
    for i in range(len(vars)):
        for j in range(len(vars)):
            assert cc.are_equal(vars[i], vars[j])


def test_composed_functions():
    f, g, h = symbols('f g h', cls=Function)
    a, b, c = symbols('a b c')
    eqs = [Eq(a, b), Eq(f(a), c), Eq(g(c), h(b))]
    cc = EUFCongruenceClosure(eqs)
    # f(a) = c and a=b => f(b) = c
    assert cc.are_equal(f(a), f(b))
    # g(c) = h(b) and c = f(a) = f(b)
    assert cc.are_equal(g(f(b)), h(b))

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
