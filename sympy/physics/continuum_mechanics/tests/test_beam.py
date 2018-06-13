from sympy import Symbol, symbols, S
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.functions import SingularityFunction, Piecewise
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


def test_composite_beam():
    E = Symbol('E')
    I = Symbol('I')
    b1 = Beam(2, E, 1.5*I)
    b2 = Beam(2, E, I)
    b = b1.join(b2, "fixed")
    b.apply_load(-20, 0, -1)
    b.apply_load(80, 0, -2)
    b.apply_load(20, 4, -1)
    b.bc_slope = [(0, 0)]
    b.bc_deflection = [(0, 0)]

    assert b.length == 4
    assert b.second_moment == Piecewise((1.5*I, x <= 2), (I, x <= 4))
    assert b.slope().subs(x, 4) == 120.0/(E*I)
    assert b.slope().subs(x, 2) == 80.0/(E*I)
    assert int(b.deflection().subs(x, 4).args[0]) == 302  # Coefficient of 1/(E*I)

    l=symbols('l', positive=True)
    M1, A1, M2, A2, P = symbols('M1 A1 M2 A2 P')
    b1=Beam(l ,E,I)
    b2=Beam(2*l ,E,I)
    b=b1.join(b2,"hinge")
    b.apply_load(A1,0,-1)
    b.apply_load(M1,0,-2)
    b.apply_load(P,2*l,-1)
    b.apply_load(A2,3*l,-1)
    b.apply_load(M2,3*l,-2)
    b.bc_slope=[(0,0), (3*l, 0)]
    b.bc_deflection=[(0,0), (3*l, 0)]
    b._solve_hinge_beams(M1, A1, M2, A2)

    assert b.reaction_loads == {A1: -5*P/18, M1: 5*P*l/18, A2: -13*P/18, M2: -4*P*l/9}
    assert b.slope().subs(x, l) == 5*P*l**2/(36*E*I)
    assert b.slope().subs(x, 2*l) == -P*l**2/(12*E*I)
    assert b.deflection().subs(x, l) == 5*P*l**3/(54*E*I)
    assert b.deflection().subs(x, 2*l) == 11*P*l**3/(108*E*I)


def test_point_cflexure():
    E = Symbol('E')
    I = Symbol('I')
    b = Beam(10, E, I)
    b.apply_load(-4, 0, -1)
    b.apply_load(-46, 6, -1)
    b.apply_load(10, 2, -1)
    b.apply_load(20, 4, -1)
    b.apply_load(3, 6, 0)

    assert b.point_cflexure() == [S(10)/3]


def test_remove_load():
    E = Symbol('E')
    I = Symbol('I')
    b = Beam(4, E, I)

    try:
        b.remove_load(2, 1, -1)
    # As no load is applied on beam, ValueError should be returned.
    except ValueError:
        assert True
    else:
        assert False

    b.apply_load(-3, 0, -2)
    b.apply_load(4, 2, -1)
    b.apply_load(-2, 2, 2, end = 3)
    b.remove_load(-2, 2, 2, end = 3)
    assert b.load == -3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1)
    assert b.applied_loads == [(-3, 0, -2, None), (4, 2, -1, None)]

    try:
        b.remove_load(1, 2, -1)
    # As load of this magnitude was never applied at
    # this position, method should return a ValueError.
    except ValueError:
        assert True
    else:
        assert False

    b.remove_load(-3, 0, -2)
    b.remove_load(4, 2, -1)
    assert b.load == 0
    assert b.applied_loads == []
