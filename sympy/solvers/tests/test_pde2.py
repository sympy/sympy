from sympy.core.function import Function
from sympy.core.relational import Eq
from sympy.core.symbol import symbols
from sympy.solvers.pde import classify_pde
from sympy.physics.units.quantities import Quantity

a, b, c, x, y = symbols('a b c x y')

def test_pde_classify_2nd_linear_constant_coeff_homogeneous():
    # Wave equation test
    x, t = symbols("x t")
    c = Quantity('c')
    u = Function('u')(x, t)
    wave_eq = Eq(u.diff(t, t) - c**2*u.diff(x, x), 0)
    # Verify classification matches the 2nd-order hint
    assert classify_pde(wave_eq, u) == ('2nd_linear_constant_coeff_homogeneous',)
    # Verify invalid order doesn't get classified
    non_linear_eq = Eq(u.diff(t)**2 - c**2*u.diff(x, x), 0)
    assert classify_pde(non_linear_eq, u) != ('2nd_linear_constant_coeff_homogeneous',)
