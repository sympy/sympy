"""Most of these tests come from the examples in Bronstein's book."""

from __future__ import division, print_function

from sympy import Poly, Matrix, S, symbols
from sympy.integrals.risch import DifferentialExtension
from sympy.integrals.prde import (prde_normal_denom, prde_special_denom,
    prde_linear_constraints, constant_system, prde_spde, prde_no_cancel_b_large,
    prde_no_cancel_b_small, limited_integrate_reduce, limited_integrate,
    is_deriv_k, is_log_deriv_k_t_radical, parametric_log_deriv_heu,
    is_log_deriv_k_t_radical_in_field)

from sympy.abc import x, t, n

t0, t1, t2, t3, k = symbols('t:4 k')


def test_prde_normal_denom():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    fa = Poly(1, t)
    fd = Poly(x, t)
    G = [(Poly(t, t), Poly(1 + t**2, t)), (Poly(1, t), Poly(x + x*t**2, t))]
    assert prde_normal_denom(fa, fd, G, DE) == \
        (Poly(x, t), (Poly(1, t), Poly(1, t)), [(Poly(x*t, t),
         Poly(t**2 + 1, t)), (Poly(1, t), Poly(t**2 + 1, t))], Poly(1, t))
    G = [(Poly(t, t), Poly(t**2 + 2*t + 1, t)), (Poly(x*t, t),
        Poly(t**2 + 2*t + 1, t)), (Poly(x*t**2, t), Poly(t**2 + 2*t + 1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_normal_denom(Poly(x, t), Poly(1, t), G, DE) == \
        (Poly(t + 1, t), (Poly((-1 + x)*t + x, t), Poly(1, t)), [(Poly(t, t),
        Poly(1, t)), (Poly(x*t, t), Poly(1, t)), (Poly(x*t**2, t),
        Poly(1, t))], Poly(t + 1, t))


def test_prde_special_denom():
    a = Poly(t + 1, t)
    ba = Poly(t**2, t)
    bd = Poly(1, t)
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_special_denom(a, ba, bd, G, DE) == \
        (Poly(t + 1, t), Poly(t**2, t), [(Poly(t, t), Poly(1, t)),
        (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))], Poly(1, t))
    G = [(Poly(t, t), Poly(1, t)), (Poly(1, t), Poly(t, t))]
    assert prde_special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), G, DE) == \
        (Poly(1, t), Poly(t**2 - 1, t), [(Poly(t**2, t), Poly(1, t)),
        (Poly(1, t), Poly(1, t))], Poly(t, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-2*x*t0, t0)]})
    DE.decrement_level()
    G = [(Poly(t, t), Poly(t**2, t)), (Poly(2*t, t), Poly(t, t))]
    assert prde_special_denom(Poly(5*x*t + 1, t), Poly(t**2 + 2*x**3*t, t), Poly(t**3 + 2, t), G, DE) == \
        (Poly(5*x*t + 1, t), Poly(0, t), [(Poly(t, t), Poly(t**2, t)),
        (Poly(2*t, t), Poly(t, t))], Poly(1, x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly((t**2 + 1)*2*x, t)]})
    G = [(Poly(t + x, t), Poly(t*x, t)), (Poly(2*t, t), Poly(x**2, x))]
    assert prde_special_denom(Poly(5*x*t + 1, t), Poly(t**2 + 2*x**3*t, t), Poly(t**3, t), G, DE) == \
        (Poly(5*x*t + 1, t), Poly(0, t), [(Poly(t + x, t), Poly(x*t, t)),
        (Poly(2*t, t, x), Poly(x**2, t, x))], Poly(1, t))
    assert prde_special_denom(Poly(t + 1, t), Poly(t**2, t), Poly(t**3, t), G, DE) == \
        (Poly(t + 1, t), Poly(0, t), [(Poly(t + x, t), Poly(x*t, t)), (Poly(2*t, t, x),
        Poly(x**2, t, x))], Poly(1, t))


def test_prde_linear_constraints():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    G = [(Poly(2*x**3 + 3*x + 1, x), Poly(x**2 - 1, x)), (Poly(1, x), Poly(x - 1, x)),
        (Poly(1, x), Poly(x + 1, x))]
    assert prde_linear_constraints(Poly(1, x), Poly(0, x), G, DE) == \
        ((Poly(2*x, x), Poly(0, x), Poly(0, x)), Matrix([[1, 1, -1], [5, 1, 1]]))
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_linear_constraints(Poly(t + 1, t), Poly(t**2, t), G, DE) == \
        ((Poly(t, t), Poly(t**2, t), Poly(t**3, t)), Matrix())
    G = [(Poly(2*x, t), Poly(t, t)), (Poly(-x, t), Poly(t, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    prde_linear_constraints(Poly(1, t), Poly(0, t), G, DE) == \
        ((Poly(0, t), Poly(0, t)), Matrix([[2*x, -x]]))


def test_constant_system():
    A = Matrix([[-(x + 3)/(x - 1), (x + 1)/(x - 1), 1],
                [-x - 3, x + 1, x - 1],
                [2*(x + 3)/(x - 1), 0, 0]])
    u = Matrix([(x + 1)/(x - 1), x + 1, 0])
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert constant_system(A, u, DE) == \
        (Matrix([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 0],
                 [0, 0, 1]]), Matrix([0, 1, 0, 0]))


def test_prde_spde():
    D = [Poly(x, t), Poly(-x*t, t)]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    # TODO: when bound_degree() can handle this, test degree bound from that too
    assert prde_spde(Poly(t, t), Poly(-1/x, t), D, n, DE) == \
        (Poly(t, t), Poly(0, t), [Poly(2*x, t), Poly(-x, t)],
        [Poly(-x**2, t), Poly(0, t)], n - 1)


def test_prde_no_cancel():
    # b large
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**2, x), Poly(1, x)], 2, DE) == \
        ([Poly(x**2 - 2*x + 2, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                        [0, 1, 0, -1]]))
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**3, x), Poly(1, x)], 3, DE) == \
        ([Poly(x**3 - 3*x**2 + 6*x - 6, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                                 [0, 1, 0, -1]]))
    # b small
    # XXX: Is there a better example of a monomial with D.degree() > 2?
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**3 + 1, t)]})

    # My original q was t**4 + t + 1, but this solution implies q == t**4
    # (c1 = 4), with some of the ci for the original q equal to 0.
    G = [Poly(t**6, t), Poly(x*t**5, t), Poly(t**3, t), Poly(x*t**2, t), Poly(1 + x, t)]
    assert prde_no_cancel_b_small(Poly(x*t, t), G, 4, DE) == \
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
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert limited_integrate_reduce(Poly(x, t), Poly(t**2, t), [(Poly(x, t),
    Poly(t, t))], DE) == \
        (Poly(t, t), Poly(-1/x, t), Poly(t, t), 1, (Poly(x, t), Poly(1, t)),
        [(Poly(-x*t, t), Poly(1, t))])


def test_limited_integrate():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    G = [(Poly(x, x), Poly(x + 1, x))]
    assert limited_integrate(Poly(-(1 + x + 5*x**2 - 3*x**3), x),
    Poly(1 - x - x**2 + x**3, x), G, DE) == \
        ((Poly(x**2 - x + 2, x), Poly(x - 1, x)), [2])
    G = [(Poly(1, x), Poly(x, x))]
    assert limited_integrate(Poly(5*x**2, x), Poly(3, x), G, DE) == \
        ((Poly(5*x**3/9, x), Poly(1, x)), [0])


def test_is_log_deriv_k_t_radical():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'E_K': [], 'L_K': [],
        'E_args': [], 'L_args': []})
    assert is_log_deriv_k_t_radical(Poly(2*x, x), Poly(1, x), DE) is None

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2*t1, t1), Poly(1/x, t2)],
        'L_K': [2], 'E_K': [1], 'L_args': [x], 'E_args': [2*x]})
    assert is_log_deriv_k_t_radical(Poly(x + t2/2, t2), Poly(1, t2), DE) == \
        ([(t1, 1), (x, 1)], t1*x, 2, 0)
    # TODO: Add more tests

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0), Poly(1/x, t)],
        'L_K': [2], 'E_K': [1], 'L_args': [x], 'E_args': [x]})
    assert is_log_deriv_k_t_radical(Poly(x + t/2 + 3, t), Poly(1, t), DE) == \
        ([(t0, 2), (x, 1)], x*t0**2, 2, 3)


def test_is_deriv_k():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)],
        'L_K': [1, 2], 'E_K': [], 'L_args': [x, x + 1], 'E_args': []})
    assert is_deriv_k(Poly(2*x**2 + 2*x, t2), Poly(1, t2), DE) == \
        ([(t1, 1), (t2, 1)], t1 + t2, 2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(t2, t2)],
        'L_K': [1], 'E_K': [2], 'L_args': [x], 'E_args': [x]})
    assert is_deriv_k(Poly(x**2*t2**3, t2), Poly(1, t2), DE) == \
        ([(x, 3), (t1, 2)], 2*t1 + 3*x, 1)
    # TODO: Add more tests, including ones with exponentials

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/x, t1)],
        'L_K': [1], 'E_K': [], 'L_args': [x**2], 'E_args': []})
    assert is_deriv_k(Poly(x, t1), Poly(1, t1), DE) == \
        ([(t1, S(1)/2)], t1/2, 1)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/(1 + x), t0)],
        'L_K': [1], 'E_K': [], 'L_args': [x**2 + 2*x + 1], 'E_args': []})
    assert is_deriv_k(Poly(1 + x, t0), Poly(1, t0), DE) == \
        ([(t0, S(1)/2)], t0/2, 1)


def test_is_log_deriv_k_t_radical_in_field():
    # NOTE: any potential constant factor in the second element of the result
    # doesn't matter, because it cancels in Da/a.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(5*t + 1, t), Poly(2*t*x, t), DE) == \
        (2, t*x**5)
    assert is_log_deriv_k_t_radical_in_field(Poly(2 + 3*t, t), Poly(5*x*t, t), DE) == \
        (5, x**3*t**2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-t/x**2, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(-(1 + 2*t), t),
    Poly(2*x**2 + 2*x**2*t, t), DE) == \
        (2, t + t**2)
    assert is_log_deriv_k_t_radical_in_field(Poly(-1, t), Poly(x**2, t), DE) == \
        (1, t)
    assert is_log_deriv_k_t_radical_in_field(Poly(1, t), Poly(2*x**2, t), DE) == \
        (2, 1/t)


def test_parametric_log_deriv():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert parametric_log_deriv_heu(Poly(5*t**2 + t - 6, t), Poly(2*x*t**2, t),
    Poly(-1, t), Poly(x*t**2, t), DE) == \
        (2, 6, t*x**5)
