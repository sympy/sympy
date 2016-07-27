from sympy import Symbol
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.functions import SingularityFunction
from sympy.utilities.pytest import XFAIL, raises

x = Symbol('x')
y = Symbol('y')


def test_Beam():
    E = Symbol('E')
    E_1 = Symbol('E_1')
    I = Symbol('I')
    I_1 = Symbol('I_1')
    b = Beam(1, E, I)
    assert b.length == 1
    assert b.elastic_modulus == E
    assert b.second_moment == I
    assert b.variable == x

    # Test the length setter
    b.length = 4
    assert b.length == 4

    # Test the E setter
    b.elastic_modulus = E_1
    assert b.elastic_modulus == E_1

    # Test the I setter
    b.second_moment = I_1
    assert b.second_moment is I_1

    # Test the variable setter
    b.variable = y
    assert b.variable is y

    # Test for all boundary conditions.
    b.bc_moment = [(0, 4), (4, 0)]
    b.boundary_conditions['deflection'] = [(0, 2)]
    b.boundary_conditions['slope'] = [(0, 1)]
    assert b.boundary_conditions == {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}

    # Test for moment boundary condition method
    b.bc_moment.extend([(4, 3), (5, 0)])
    m_bcs = b.bc_moment
    assert m_bcs == [(0, 4), (4, 0), (4, 3), (5, 0)]

    # Test for slope boundary condition method
    b.boundary_conditions['slope'].extend([(4, 3), (5, 0)])
    s_bcs = b.boundary_conditions['slope']
    assert s_bcs == [(0, 1), (4, 3), (5, 0)]

    # Test for deflection boundary condition method
    b.boundary_conditions['deflection'].extend([(4, 3), (5, 0)])
    d_bcs = b.boundary_conditions['deflection']
    assert d_bcs == [(0, 2), (4, 3), (5, 0)]

    # Test for updated boundary conditions
    bcs_new = b.boundary_conditions
    assert bcs_new == {
        'deflection': [(0, 2), (4, 3), (5, 0)],
        'moment': [(0, 4), (4, 0), (4, 3), (5, 0)],
        'slope': [(0, 1), (4, 3), (5, 0)]}

    b1 = Beam(2, E, I)
    b1.apply_load(order=-2, value=-3, start=0)
    b1.apply_load(order=-1, value=4, start=2)
    b1.apply_load(order=2, value=-2, start=3)

    b1.bc_moment = [(0, 4), (4, 0)]
    b1.boundary_conditions['deflection'] = [(0, 2)]
    b1.boundary_conditions['slope'] = [(0, 1)]

    # Test for load distribution function.
    p = b1.load
    q = -3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1) - 2*SingularityFunction(x, 3, 2)
    assert p == q

    # Test for shear force distribution function
    p = b1.shear_force()
    q = -3*SingularityFunction(x, 0, -1) + 4*SingularityFunction(x, 2, 0) - 2*SingularityFunction(x, 3, 3)/3 - 71/24
    assert p == q

    # Test for bending moment distribution function
    p = b1.bending_moment()
    q = -71*x/24 - 3*SingularityFunction(x, 0, 0) + 4*SingularityFunction(x, 2, 1) - SingularityFunction(x, 3, 4)/6 + 7
    assert p == q

    # Test for slope distribution function
    p = b1.slope()
    q = -71*x**2/48 + 7*x - 3*SingularityFunction(x, 0, 1) + 2*SingularityFunction(x, 2, 2) - SingularityFunction(x, 3, 5)/30 + 1
    assert p == q/(E*I)

    # Test for deflection distribution function
    p = b1.deflection()
    q = -71*x**3/144 + 7*x**2/2 + x - 3*SingularityFunction(x, 0, 2)/2 + 2*SingularityFunction(x, 2, 3)/3 - SingularityFunction(x, 3, 6)/180 + 2
    assert p == q/(E*I)
