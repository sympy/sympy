from sympy.physics.continuum_mechanics.column import Column
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions import SingularityFunction

x = Symbol('x')
x_1 = Symbol('x_1')

def test_column():
    E = Symbol('E')
    E_1 = Symbol('E_1')
    A = Symbol('A')
    A_1 = Symbol('A_1')

    # Test the variables
    c = Column(10, E, A)
    assert c.length == 10
    assert c.elastic_modulus == E
    assert c.area == A
    assert c.variable == x

    # Test the str message
    assert str(c) == 'Column(10, E, A)'

    # Test the length setter
    c.length = 20
    assert c.length == 20

    # Test the elastic modulus setter
    c.elastic_modulus = E_1
    assert c.elastic_modulus == E_1

    # Test the area setter
    c.area = A_1
    assert c.area == A_1

    # Test the variable setter
    c.variable = x_1
    assert c.variable == x_1

    # Test applying supports
    c2 = Column(10, E, A)
    c2.apply_support(0)
    assert c2._bc_deflection == [(0, 0)]
    c2.apply_support(10)
    assert c2._bc_deflection == [(0, 0), (10, 0)]
    c2.apply_support(5)
    assert c2._bc_deflection == [(0, 0), (10, 0), (5, 0)]

    # Test applying a load
    c3 = Column(10, E, A)
    c3.apply_load(-10, 10, -1)
    assert c3.applied_loads == [(-10, 10, -1)]
    c3.apply_load(10, 5, -1)
    assert c3.applied_loads == [(-10, 10, -1), (10, 5, -1)]

    # Test the load equation
    p = c3.load
    q = -10 * SingularityFunction(x, 10, -1) + 10 * SingularityFunction(x, 5, -1)
    assert p == q

    # Test applying a symbolic load
    c4 = Column(10, E, A)
    F, G = symbols('F G')
    c4.apply_load(F, 0, -1)
    c4.apply_load(-G, 10, -1)
    assert c4.applied_loads == [(F, 0, -1), (-G, 10, -1)]

    # Test the load equation
    p = c4.load
    q = F * SingularityFunction(x, 0, -1) - G * SingularityFunction(x, 10, -1)
    assert p == q

    # Test the load equation for column with a support
    c5 = Column(10, E, A)
    R_0 = c5.apply_support(0)
    c5.apply_load(10, 0, -1)
    c5.apply_load(10, 5, -1)
    
    p = c5.load
    q = R_0 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 5, -1)
    assert p == q

test_column()