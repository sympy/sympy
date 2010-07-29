"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, Matrix, S
from sympy.integrals.prde import (prde_normal_denom, prde_special_denom,
    prde_linear_constraints, constant_system, prde_spde, prde_no_cancel_b_large,
    prde_no_cancel_b_small, limited_integrate_reduce, limited_integrate,
    is_log_deriv_k_t_radical, parametric_log_deriv_heu)

from sympy.abc import x, t, n

def test_prde_normal_denom():
    D = [Poly(1, x), Poly(1 + t**2, t)]
    fa = Poly(1, t)
    fd = Poly(x, t)
    G = [(Poly(t, t), Poly(1 + t**2, t)), (Poly(1, t), Poly(x + x*t**2, t))]
    assert prde_normal_denom(fa, fd, G, D, [x, t]) == \
        (Poly(x, t), (Poly(1, t), Poly(1, t)), [(Poly(x*t, t),
         Poly(t**2 + 1, t)), (Poly(1, t), Poly(t**2 + 1, t))], Poly(1, t))
    G = [(Poly(t, t), Poly(t**2 + 2*t + 1, t)), (Poly(x*t, t),
        Poly(t**2 + 2*t + 1, t)), (Poly(x*t**2, t), Poly(t**2 + 2*t + 1, t))]
    D = [Poly(1, x), Poly(t, t)]
    assert prde_normal_denom(Poly(x, t), Poly(1, t), G, D, [x, t]) == \
        (Poly(t + 1, t), (Poly((-1 + x)*t + x, t), Poly(1, t)), [(Poly(t, t),
        Poly(1, t)), (Poly(x*t, t), Poly(1, t)), (Poly(x*t**2, t),
        Poly(1, t))], Poly(t + 1, t))

def test_prde_special_denom():
    a = Poly(t + 1, t)
    ba = Poly(t**2, t)
    bd = Poly(1, t)
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    D = [Poly(1, x), Poly(t, t)]
    assert prde_special_denom(a, ba, bd, G, D, [x, t]) == \
        (Poly(t + 1, t), Poly(t**2, t), [(Poly(t, t), Poly(1, t)),
        (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))], Poly(1, t))
    G = [(Poly(t, t), Poly(1, t)), (Poly(1, t), Poly(t, t))]
    assert prde_special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), G, D, [x, t]) == \
        (Poly(1, t), Poly(t**2 - 1, t), [(Poly(t**2, t), Poly(1, t)),
        (Poly(1, t), Poly(1, t))], Poly(t, t))

def test_prde_linear_constraints():
    D = [Poly(1, x)]
    G = [(Poly(2*x**3 + 3*x + 1, x), Poly(x**2 - 1, x)), (Poly(1, x), Poly(x - 1, x)),
        (Poly(1, x), Poly(x + 1, x))]
    assert prde_linear_constraints(Poly(1, x), Poly(0, x), G, D, [x]) == \
        ((Poly(2*x, x), Poly(0, x), Poly(0, x)), Matrix([[1, 1, -1], [5, 1,  1]]))
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    D = [Poly(1, x), Poly(t, t)]
    assert prde_linear_constraints(Poly(t + 1, t), Poly(t**2, t), G, D, [x, t]) == \
        ((Poly(t, t), Poly(t**2, t), Poly(t**3, t)), Matrix())
    G = [(Poly(2*x, t), Poly(t, t)), (Poly(-x, t), Poly(t, t))]
    D = [Poly(1, x), Poly(1/x, t)]
    prde_linear_constraints(Poly(1, t), Poly(0, t), G, D, [x, t]) == \
        ((Poly(0, t), Poly(0, t)), Matrix([[2*x, -x]]))

def test_constant_system():
    A = Matrix([[-(x + 3)/(x - 1), (x + 1)/(x - 1), 1],
                [-x - 3, x + 1, x - 1],
                [2*(x + 3)/(x - 1), 0, 0]])
    u = Matrix([(x + 1)/(x - 1), x + 1, 0])
    assert constant_system(A, u, [Poly(1, x)], [x]) == \
        (Matrix([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 0],
                 [0, 0, 1]]), Matrix([0, 1, 0, 0]))

def test_prde_spde():
    D = [Poly(x, t), Poly(-x*t, t)]
    # TODO: when bound_degree() can handle this, test degree bound from that too
    assert prde_spde(Poly(t, t), Poly(-1/x, t), D, n, [Poly(1, x), Poly(1/x, t)], [x, t]) == \
        (Poly(t, t), Poly(0, t), [Poly(2*x, t), Poly(-x, t)],
        [Poly(-x**2, t), Poly(0, t)], n - 1)

def test_prde_no_cancel():
    # b large
    D = [Poly(1, x)]
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**2, x), Poly(1, x)], 2, D, [x]) == \
        ([Poly(x**2 - 2*x + 2, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                        [0, 1, 0, -1]]))
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**3, x), Poly(1, x)], 3, D, [x]) == \
        ([Poly(x**3 - 3*x**2 + 6*x - 6, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                                 [0, 1, 0, -1]]))
    # b small
    D = [Poly(1, x), Poly(t**3 + 1, t)]
    # XXX: Is there a better example of a monomial with D.degree() > 2?

    # My original q was t**4 + t + 1, but this solution implies q == t**4
    # (c1 = 4), with some of the ci for the original q equal to 0.
    G = [Poly(t**6, t), Poly(x*t**5, t), Poly(t**3, t), Poly(x*t**2, t), Poly(1 + x, t)]
    assert prde_no_cancel_b_small(Poly(x*t, t), G, 4, D, [x, t]) == \
        ([Poly(t**4/4 - x/12*t**3 + x**2/24*t**2 + (-S(11)/12 - x**3/24)*t + x/24, t),
        Poly(x/3*t**3 - x**2/6*t**2 + (-S(1)/3 + x**3/6)*t - x/6, t), Poly(t, t),
        Poly(0, t), Poly(0, t)], Matrix([[1, 0,      -1, 0, 0,  0,  0,  0,  0,  0],
                                         [0, 1, -S(1)/4, 0, 0,  0,  0,  0,  0,  0],
                                         [0, 0,       0, 0, 0,  0,  0,  0,  0,  0],
                                         [0, 0,       0, 1, 0,  0,  0,  0,  0,  0],
                                         [0, 0,       0, 0, 1,  0,  0,  0,  0,  0],
                                         [1, 0,       0, 0, 0, -1,  0,  0,  0,  0],
                                         [0, 1,       0, 0, 0,  0, -1,  0,  0,  0],
                                         [0, 0,       1, 0, 0,  0,  0, -1,  0,  0],
                                         [0, 0,       0, 1, 0,  0,  0,  0, -1,  0],
                                         [0, 0,       0, 0, 1,  0,  0,  0,  0, -1]]))

    # TODO: Add test for deg(b) <= 0 with b small

def test_limited_integrate_reduce():
    D = [Poly(1, x), Poly(1/x, t)]
    assert limited_integrate_reduce(Poly(x, t), Poly(t**2, t), [(Poly(x, t),
    Poly(t, t))], D, [x, t]) == \
        (Poly(t, t), Poly(-1/x, t), Poly(t, t), 1, (Poly(x, t), Poly(1, t)),
        [(Poly(-x*t, t), Poly(1, t))])

def test_limited_integrate():
    D = [Poly(1, x)]
    G = [(Poly(x, x), Poly(x + 1, x))]
    assert limited_integrate(Poly(-(1 + x + 5*x**2 - 3*x**3), x),
    Poly(1 - x - x**2 + x**3, x), G, D, [x]) == \
        ((Poly(x**2 - x + 2, x), Poly(x - 1, x)), [2])

def test_is_log_deriv_k_t_radical():
    pass
