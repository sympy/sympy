from sympy import Symbol
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.functions import SingularityFunction

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
    b.bc_deflection = [(0, 2)]
    b.bc_slope = [(0, 1)]
    assert b.boundary_conditions == {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}

    # Test for moment boundary condition method
    b.bc_moment.extend([(4, 3), (5, 0)])
    m_bcs = b.bc_moment
    assert m_bcs == [(0, 4), (4, 0), (4, 3), (5, 0)]

    # Test for slope boundary condition method
    b.bc_slope.extend([(4, 3), (5, 0)])
    s_bcs = b.bc_slope
    assert s_bcs == [(0, 1), (4, 3), (5, 0)]

    # Test for deflection boundary condition method
    b.bc_deflection.extend([(4, 3), (5, 0)])
    d_bcs = b.bc_deflection
    assert d_bcs == [(0, 2), (4, 3), (5, 0)]

    # Test for updated boundary conditions
    bcs_new = b.boundary_conditions
    assert bcs_new == {
        'deflection': [(0, 2), (4, 3), (5, 0)],
        'moment': [(0, 4), (4, 0), (4, 3), (5, 0)],
        'slope': [(0, 1), (4, 3), (5, 0)]}

    b1 = Beam(2, E, I)
    b1.apply_load(-3, 0, -2)
    b1.apply_load(4, 2, -1)
    b1.apply_load(-2, 3, 2)

    b1.bc_moment = [(0, 4), (4, 0)]
    b1.bc_deflection = [(0, 2)]
    b1.bc_slope = [(0, 1)]

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

    # Test using symbols
    l = Symbol('l')
    w0 = Symbol('w0')
    w2 = Symbol('w2')
    a1 = Symbol('a1')
    c = Symbol('c')
    c1 = Symbol('c1')
    d = Symbol('d')
    e = Symbol('e')
    f = Symbol('f')

    b2 = Beam(l, E, I)

    b2.apply_load(w0, a1, 1)
    b2.apply_load(w2, c1, -1)

    b2.bc_deflection = [(c, d)]
    b2.bc_slope = [(e, f)]

    # Test for load distribution function.
    p = b2.load
    q = w0*SingularityFunction(x, a1, 1) + w2*SingularityFunction(x, c1, -1)
    assert p == q

    # Test for shear force distribution function
    p = b2.shear_force()
    q = w0*SingularityFunction(x, a1, 2)/2 + w2*SingularityFunction(x, c1, 0)
    assert p == q

    # Test for bending moment distribution function
    p = b2.bending_moment()
    q = -w0*SingularityFunction(x, a1, 3)/6 - w2*SingularityFunction(x, c1, 1)
    assert p == q

    # Test for slope distribution function
    p = b2.slope()
    q = (f + w0*SingularityFunction(e, a1, 4)/24 - w0*SingularityFunction(x, a1, 4)/24 + w2*SingularityFunction(e, c1, 2)/2 - w2*SingularityFunction(x, c1, 2)/2)
    assert p == q/(E*I)

    # Test for deflection distribution function
    p = b2.deflection()
    q = (-c*f - c*w0*SingularityFunction(e, a1, 4)/24 - c*w2*SingularityFunction(e, c1, 2)/2 + d + f*x + w0*x*SingularityFunction(e, a1, 4)/24 + w0*SingularityFunction(c, a1, 5)/120 - w0*SingularityFunction(x, a1, 5)/120 + w2*x*SingularityFunction(e, c1, 2)/2 + w2*SingularityFunction(c, c1, 3)/6 - w2*SingularityFunction(x, c1, 3)/6)
    assert p == q/(E*I)

    b3 = Beam(9, E, I)
    b3.apply_load(value=-2, start=2, order=2, end=3)
    p = b3.load
    q = - 2*SingularityFunction(x, 2, 2) + 2*SingularityFunction(x, 3, 0) + 2*SingularityFunction(x, 3, 2)
    assert p == q
