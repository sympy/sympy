from sympy.physics.mechanics import Body
from sympy.physics.mechanics.joints import Joint
from sympy.testing.pytest import raises

def test_abstract_Joint():
    P = Body('P')
    C = Body('B')
    raises(TypeError, lambda: Joint('J', P, C))
