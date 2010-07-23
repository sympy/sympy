"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly, Matrix
from sympy.integrals.prde import (prde_normal_denom, prde_special_denom,
    prde_linear_constraints, constant_system, par_spde,
    is_log_deriv_k_t_radical, parametric_log_deriv_heu)

from sympy.abc import x, t

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

def test_par_spde():
    D = [Poly(x, t), Poly(-x*t, t)]
    assert par_spde(Poly(t, t), Poly(-1/x, t), D, 2, [Poly(1, x), Poly(1/x, t)], [x, t]) == \
        (Poly(t, t), Poly(0, t), [Poly(2*x, t), Poly(-x, t)],
        [Poly(-x**2, t), Poly(0, t)], 1)

def test_is_log_deriv_k_t_radical():
    pass
