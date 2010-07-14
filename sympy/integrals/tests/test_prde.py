"""Most of these tests come from the examples in Bronstein's book."""
from sympy import Poly
from sympy.integrals.prde import (prde_normal_denom, prde_special_denom,
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
        (Poly(t + 1, t), Poly(t**2, t), [Poly(t, t), Poly(t**2, t), Poly(t**3, t)], Poly(1, t))
    G = [(Poly(t, t), Poly(1, t)), (Poly(1, t), Poly(t, t))]
    assert prde_special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), G, D, [x, t])
        (Poly(1, t), Poly(t**2 - 1, t), [Poly(t**2, t), Poly(1, t)], Poly(t, t))

def test_is_log_deriv_k_t_radical():
    pass
