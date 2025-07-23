from sympy import symbols, Eq
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure, EUFFunction


def test_euf_basic_chain_equality():
    """
    Test simple equality chain propagation: x = y, y = z ⇒ x = z.
    """

    x, y, z = symbols("x y z")
    eqs = [Eq(x, y), Eq(y, z)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(x, z)
    assert cc.are_equal(y, z)
    assert cc.are_equal(x, y)


def test_euf_function_congruence_equal_args():
    """
    Test congruence: if x = y ⇒ f(x) = f(y);
    and verify propagation of equality through function applications.
    """

    x, y = symbols("x y")
    a, b = symbols("a b")
    f = EUFFunction('f', 1)  # unary uninterpreted function
    fx, fy = f(x), f(y)
    eqs = [Eq(x, y), Eq(fx, a), Eq(fy, b)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(fx, fy)
    assert cc.are_equal(a, b)


def test_euf_nested_function_application():
    """
    Test nested function term propagation:
    If f(x) = f(y), then g(f(x)) = g(f(y)) by congruence.
    """

    x, y = symbols("x y")
    f = EUFFunction('f', 1)
    g = EUFFunction('g', 1)
    fx, fy = f(x), f(y)
    gfx, gfy = g(fx), g(fy)
    eqs = [Eq(x, y), Eq(gfx, gfy)]
    cc = EUFCongruenceClosure(eqs)
    assert cc.are_equal(fx, fy)
    assert cc.are_equal(gfx, gfy)


def test_euf_different_functions_do_not_propagate():
    """
    If f(x) ≠ g(x), even when x is equal, f and g must not cause spurious merges.
    """

    x = symbols("x")
    f = EUFFunction('f', 1)
    g = EUFFunction('g', 1)
    fx = f(x)
    gx = g(x)
    a, b = symbols("a b")
    eqs = [Eq(fx, a), Eq(gx, b)]
    cc = EUFCongruenceClosure(eqs)
    assert not cc.are_equal(fx, gx)
    assert not cc.are_equal(a, b)


def test_euf_nested_equalities():
    """
    This test demonstrates congruence propagation through nested applications.

    Given:
        x = y
        f(x) = a
        f(y) = b
        g(f(x)) = g(f(y))

    The congruence closure should determine:
        f(x) == f(y)        => True  (from x = y)
        a == b              => True  (f(x) = a, f(y) = b, f(x) == f(y))
        g(f(x)) == g(f(y))  => True  (by congruence)
        x == z              => False (z is unconnected)
    """

    x, y, z = symbols("x y z")
    a, b = symbols("a b")
    f = EUFFunction('f', 1)
    g = EUFFunction('g', 1)

    fx, fy = f(x), f(y)
    gfx, gfy = g(fx), g(fy)

    eqs = [
        Eq(x, y),
        Eq(fx, a),
        Eq(fy, b),
        Eq(gfx, gfy),
    ]

    cc = EUFCongruenceClosure(eqs)

    assert cc.are_equal(fx, fy)      # f(x) = f(y) from x = y
    assert cc.are_equal(a, b)        # f(x) = a, f(y) = b ⇒ a = b
    assert cc.are_equal(gfx, gfy)    # congruence: g(f(x)) = g(f(y))
    assert not cc.are_equal(x, z)    # unrelated term z ⇒ different class
