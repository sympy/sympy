# from sympy import symbols, Eq
# from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure, EUFFunction


# def test_euf_basic_chain_equality():
#     """
#     Test simple equality chain propagation: x = y, y = z => x = z.
#     """

#     x, y, z = symbols("x y z")
#     eqs = [Eq(x, y), Eq(y, z)]
#     cc = EUFCongruenceClosure(eqs)
#     assert cc.are_equal(x, z)
#     assert cc.are_equal(y, z)
#     assert cc.are_equal(x, y)


# def test_euf_function_congruence_equal_args():
#     """
#     Test congruence: if x = y => f(x) = f(y);
#     and verify propagation of equality through function applications.
#     """

#     x, y = symbols("x y")
#     a, b = symbols("a b")
#     f = EUFFunction('f', 1)  # unary uninterpreted function
#     fx, fy = f(x), f(y)
#     eqs = [Eq(x, y), Eq(fx, a), Eq(fy, b)]
#     cc = EUFCongruenceClosure(eqs)
#     assert cc.are_equal(fx, fy)
#     assert cc.are_equal(a, b)


# def test_euf_nested_function_application():
#     """
#     Test nested function term propagation:
#     If f(x) = f(y), then g(f(x)) = g(f(y)) by congruence.
#     """

#     x, y = symbols("x y")
#     f = EUFFunction('f', 1)
#     g = EUFFunction('g', 1)
#     fx, fy = f(x), f(y)
#     gfx, gfy = g(fx), g(fy)
#     eqs = [Eq(x, y), Eq(gfx, gfy)]
#     cc = EUFCongruenceClosure(eqs)
#     assert cc.are_equal(fx, fy)
#     assert cc.are_equal(gfx, gfy)


# def test_euf_different_functions_do_not_propagate():
#     """
#     If f(x) != g(x), even when x is equal, f and g must not cause spurious merges.
#     """

#     x = symbols("x")
#     f = EUFFunction('f', 1)
#     g = EUFFunction('g', 1)
#     fx = f(x)
#     gx = g(x)
#     a, b = symbols("a b")
#     eqs = [Eq(fx, a), Eq(gx, b)]
#     cc = EUFCongruenceClosure(eqs)
#     assert not cc.are_equal(fx, gx)
#     assert not cc.are_equal(a, b)


# def test_euf_nested_equalities():
#     """
#     This test demonstrates congruence propagation through nested applications.

#     Given:
#         x = y
#         f(x) = a
#         f(y) = b
#         g(f(x)) = g(f(y))

#     The congruence closure should determine:
#         f(x) == f(y)        => True  (from x = y)
#         a == b              => True  (f(x) = a, f(y) = b, f(x) == f(y))
#         g(f(x)) == g(f(y))  => True  (by congruence)
#         x == z              => False (z is unconnected)
#     """

#     x, y, z = symbols("x y z")
#     a, b = symbols("a b")
#     f = EUFFunction('f', 1)
#     g = EUFFunction('g', 1)

#     fx, fy = f(x), f(y)
#     gfx, gfy = g(fx), g(fy)

#     eqs = [
#         Eq(x, y),
#         Eq(fx, a),
#         Eq(fy, b),
#         Eq(gfx, gfy),
#     ]

#     cc = EUFCongruenceClosure(eqs)

#     assert cc.are_equal(fx, fy)      # f(x) = f(y) from x = y
#     assert cc.are_equal(a, b)        # f(x) = a, f(y) = b => a = b
#     assert cc.are_equal(gfx, gfy)    # congruence: g(f(x)) = g(f(y))
#     assert not cc.are_equal(x, z)    # unrelated term z => different class



from sympy import Symbol, Dummy, Eq, Lambda
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
from sympy.abc import x, y, z, w

def test_basic_equality():
    cc = EUFCongruenceClosure([Eq(x, y), Eq(y, z)])
    assert cc.are_equal(x, y)
    assert cc.are_equal(y, z)
    assert cc.are_equal(x, z)
    assert not cc.are_equal(x, w)


def test_curried_lambda_congruence():
    lam = Lambda(x, y + x)
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
