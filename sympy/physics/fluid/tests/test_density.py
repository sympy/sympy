from sympy.core import symbols
from sympy.physics.fluid.fluid_prop import density

def test_density():
    """
    performs test to check if density function
    is able to find correct density
    """
    m, v = symbols('m, v')
    assert density(m, v) == m/v
