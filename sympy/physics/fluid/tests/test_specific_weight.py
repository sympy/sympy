from sympy.core import symbols
from sympy.physics.fluid.fluid_prop import specific_weight

def test_specific_weight():
    """
    performs test to check if specific_weight function
    is able to find correct specific weight
    """
    m, v, g = symbols('m, v, g')
    assert specific_weight(m, v, g) == g*m/v
