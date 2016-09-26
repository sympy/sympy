from sympy import Symbol, symbols
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.functions import SingularityFunction
from sympy.utilities.pytest import raises

x = Symbol('x')
y = Symbol('y')
R1, R2 = symbols('R1, R2')


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
    b.bc_deflection = [(0, 2)]
    b.bc_slope = [(0, 1)]
    assert b.boundary_conditions == {'deflection': [(0, 2)], 'slope': [(0, 1)]}

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
        'slope': [(0, 1), (4, 3), (5, 0)]}

    b1 = Beam(30, E, I)
    b1.apply_load(-8, 0, -1)
    b1.apply_load(R1, 10, -1)
    b1.apply_load(R2, 30, -1)
    b1.apply_load(120, 30, -2)
    b1.bc_deflection = [(10, 0), (30, 0)]
    b1.solve_for_reaction_loads(R1, R2)

    # Test for finding reaction forces
    p = b1.reaction_loads
    q = {R1: 6, R2: 2}
    assert p == q

    # Test for load distribution function.
    p = b1.load
    q = -8*SingularityFunction(x, 0, -1) + 6*SingularityFunction(x, 10, -1) + 120*SingularityFunction(x, 30, -2) + 2*SingularityFunction(x, 30, -1)
    assert p == q

    # Test for shear force distribution function
    p = b1.shear_force()
    q = -8*SingularityFunction(x, 0, 0) + 6*SingularityFunction(x, 10, 0) + 120*SingularityFunction(x, 30, -1) + 2*SingularityFunction(x, 30, 0)
    assert p == q

    # Test for bending moment distribution function
    p = b1.bending_moment()
    q = -8*SingularityFunction(x, 0, 1) + 6*SingularityFunction(x, 10, 1) + 120*SingularityFunction(x, 30, 0) + 2*SingularityFunction(x, 30, 1)
    assert p == q

    # Test for slope distribution function
    p = b1.slope()
    q = -4*SingularityFunction(x, 0, 2) + 3*SingularityFunction(x, 10, 2) + 120*SingularityFunction(x, 30, 1) + SingularityFunction(x, 30, 2) + 4000/3
    assert p == q/(E*I)

    # Test for deflection distribution function
    p = b1.deflection()
    q = 4000*x/3 - 4*SingularityFunction(x, 0, 3)/3 + SingularityFunction(x, 10, 3) + 60*SingularityFunction(x, 30, 2) + SingularityFunction(x, 30, 3)/3 - 12000
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
    q = w0*SingularityFunction(x, a1, 3)/6 + w2*SingularityFunction(x, c1, 1)
    assert p == q

    # Test for slope distribution function
    p = b2.slope()
    q = (f - w0*SingularityFunction(e, a1, 4)/24 + w0*SingularityFunction(x, a1, 4)/24 - w2*SingularityFunction(e, c1, 2)/2 + w2*SingularityFunction(x, c1, 2)/2)
    assert p == q/(E*I)

    # Test for deflection distribution function
    p = b2.deflection()
    q = (-c*f + c*w0*SingularityFunction(e, a1, 4)/24 + c*w2*SingularityFunction(e, c1, 2)/2 + d + f*x - w0*x*SingularityFunction(e, a1, 4)/24 - w0*SingularityFunction(c, a1, 5)/120 + w0*SingularityFunction(x, a1, 5)/120 - w2*x*SingularityFunction(e, c1, 2)/2 - w2*SingularityFunction(c, c1, 3)/6 + w2*SingularityFunction(x, c1, 3)/6)
    assert p == q/(E*I)

    b3 = Beam(9, E, I)
    b3.apply_load(value=-2, start=2, order=2, end=3)
    b3.bc_slope.append((0, 2))
    p = b3.load
    q = - 2*SingularityFunction(x, 2, 2) + 2*SingularityFunction(x, 3, 0) + 2*SingularityFunction(x, 3, 2)
    assert p == q

    p = b3.slope()
    q = -SingularityFunction(x, 2, 5)/30 + SingularityFunction(x, 3, 3)/3 + SingularityFunction(x, 3, 5)/30 + 2
    assert p == q/(E*I)

    p = b3.deflection()
    q = 2*x - SingularityFunction(x, 2, 6)/180 + SingularityFunction(x, 3, 4)/12 + SingularityFunction(x, 3, 6)/180
    assert p == q/(E*I)

    b4 = Beam(4, E, I)
    b4.apply_load(-3, 0, 0, end=3)

    p = b4.load
    q = -3*SingularityFunction(x, 0, 0) + 3*SingularityFunction(x, 3, 0)
    assert p == q

    p = b4.slope()
    q = -3*SingularityFunction(x, 0, 3)/6 + 3*SingularityFunction(x, 3, 3)/6
    assert p == q/(E*I)

    p = b4.deflection()
    q = -3*SingularityFunction(x, 0, 4)/24 + 3*SingularityFunction(x, 3, 4)/24
    assert p == q/(E*I)

    raises(ValueError, lambda: b4.apply_load(-3, 0, -1, end=3))
    with raises(TypeError):
        b4.variable = 1
