from sympy.physics.continuum_mechanics.column import Column
from sympy.functions import SingularityFunction
from sympy.core.numbers import Rational
from sympy.core.symbol import (Symbol, symbols)


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
    assert c2._bc_extension == [0]
    c2.apply_support(10)
    assert c2._bc_extension == [0, 10]
    c2.apply_support(5)
    assert c2._bc_extension == [0, 10, 5]

    # Test applying a load
    c3 = Column(10, E, A)
    c3.apply_load(-10, 10, -1)
    assert c3.applied_loads == [(-10, 10, -1, None)]
    c3.apply_load(10, 5, -1)
    assert c3.applied_loads == [(-10, 10, -1, None), (10, 5, -1, None)]

    # Test the load equation
    p = c3.load
    q = -10 * SingularityFunction(x, 10, -1) + 10 * SingularityFunction(x, 5, -1)
    assert p == q

    # Test applying a symbolic load
    c4 = Column(10, E, A)
    F, G = symbols('F G')
    c4.apply_load(F, 0, -1)
    c4.apply_load(-G, 10, -1)
    assert c4.applied_loads == [(F, 0, -1, None), (-G, 10, -1, None)]

    # Test the load equation
    p = c4.load
    q = F * SingularityFunction(x, 0, -1) - G * SingularityFunction(x, 10, -1)
    assert p == q

    # Test the load equation for column with a support
    c5 = Column(10, E, A)
    c5.apply_support(0)
    c5.apply_load(10, 0, -1)
    c5.apply_load(10, 5, -1)

    p = c5.load
    R_0 = Symbol('R_0')
    q = R_0 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 5, -1)
    assert p == q

test_column()


def test_supports_using_manual_method():
    # Test the load equation for column with a support
    E = Symbol('E')
    A = Symbol('A')
    R_0 = Symbol('R_0')
    c5 = Column(10, E, A)
    c5.apply_load(R_0,0,-1)
    c5._bc_extension.append(0)
    c5.apply_load(10, 0, -1)
    c5.apply_load(10, 5, -1)

    p = c5.load
    q = R_0 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 0, -1) + 10 * SingularityFunction(x, 5, -1)
    assert p == q

    # Test for one fixed support at x=0
    E, A = symbols('E A')
    R_0 = Symbol('R_0')

    c = Column(5, E, A)
    c.apply_load(R_0,0,-1)
    c._bc_extension.append(0)
    c.apply_load(-1, 2.5, -1)
    c.apply_load(2, 5, -1)

    c.solve_for_reaction_loads(R_0)
    p = c.reaction_loads
    q = {R_0: -1}
    assert p == q

    # Test for one fixed support at x=L
    c1 = Column(5, E, A)
    R_5 = Symbol('R_5')
    c1.apply_load(R_5,5,-1)
    c1._bc_extension.append(5)
    c1.apply_load(-4, 2.5, -1)
    c1.apply_load(2, 5, -1)

    c1.solve_for_reaction_loads(R_5)
    p = c1.reaction_loads
    q = {R_5: 2}
    assert p == q

    # Test for two supports at ends
    c2 = Column(10, E, A)
    R_0 = Symbol('R_0')
    R_10 = Symbol('R_10')
    c2.apply_load(R_10,10,-1)
    c2.apply_load(R_0,0,-1)
    c2._bc_extension.append(0)
    c2._bc_extension.append(10)
    c2.apply_load(-1, 5, -1)

    c2.solve_for_reaction_loads(R_0,R_10)
    p = c2.reaction_loads
    q = {R_0: Rational(1,2), R_10: Rational(1,2)}
    assert p == q


def test_distributed_loads():
    E, A = symbols('E A')
    c = Column(10, E, A)

    # Test different load orders
    c.apply_load(5, 0, -1)
    c.apply_load(5, 2, 0, end=4)
    c.apply_load(5, 4, 1, end=6)

    p = c.applied_loads
    q = [(5, 0, -1, None), (5, 2, 0, 4), (5, 4, 1, 6)]
    assert p == q

    p = c.load
    q = (
        5 * SingularityFunction(x, 0, -1) +
        5 * SingularityFunction(x, 2, 0) -
        5 * SingularityFunction(x, 4, 0) +
        5 * SingularityFunction(x, 4, 1) -
        10 * SingularityFunction(x, 6, 0) -
        5 * SingularityFunction(x, 6, 1)
    )
    assert p == q

    # Test distributed loads symbolically
    E, A, F, L = symbols('E A F L')
    c1 = Column(L, E, A)

    c1.apply_load(F, 0, -1)
    c1.apply_load(F, L/4, 0, end=L/2)
    c1.apply_load(F, L/2, 1, end=L)

    # Test applied loads
    p = c1.applied_loads
    q = [(F, 0, -1, None), (F, L/4, 0, L/2), (F, L/2, 1, L)]
    assert p == q

    # Test symbolic load expression
    p = c1.load
    q = (
        F * SingularityFunction(x, 0, -1) +
        F * SingularityFunction(x, L/4, 0) -
        F * SingularityFunction(x, L/2, 0) +
        F * SingularityFunction(x, L/2, 1) -
        (L/2) * F * SingularityFunction(x, L, 0) -
        F * SingularityFunction(x, L, 1)
    )
    assert p == q

test_distributed_loads()


def test_remove_load():
    E, A = symbols('E A')
    c = Column(10, E, A)
    c.apply_load(10, 0, -1)
    c.apply_load(-10, 10, -1)

    # Test applied loads
    c.remove_load(10, 0, -1, None)
    p = c.applied_loads
    q = [(-10, 10, -1, None)]
    assert p == q

    # Test load equation
    c.apply_load(5, 0, 0, end=10)
    c.remove_load(5, 0, 0, end=10)
    p = c.load
    q = -10 * SingularityFunction(x, 10, -1)
    assert p == q

    # Test load equation for higher orders
    c.apply_load(5, 0, 1, end=5)
    c.remove_load(5, 0, 1, 5)
    p = c.load
    q = -10 * SingularityFunction(x, 10, -1)
    assert p == q

    # Test symbolically
    E, A, F, L = symbols('E A F L')
    c = Column(L, E, A)
    c.apply_load(F, 0, -1)
    c.apply_load(-F, L, -1)

    c.remove_load(F, 0, -1, None)

    p = c.applied_loads
    q = [(-F, L, -1, None)]
    assert p == q

    c.apply_load(F, 0, 0, end=L)
    c.remove_load(F, 0, 0, end=L)

    p = c.load
    q = -F * SingularityFunction(x, L, -1)
    assert p == q

    c.apply_load(F, 0, 1, end=L/2)
    c.remove_load(F, 0, 1, L/2)

    p = c.load
    q = -F * SingularityFunction(x, L, -1)
    assert p == q

    # Ramp load in opposite direction
    c1 = Column(10, E, A)
    c1.apply_load(1, 5, -1)
    c1.apply_load(10, 10, 1, end=0)
    c1.remove_load(10, 10, 1, end=0)

    p = c1.applied_loads
    q = [(1, 5, -1, None)]
    assert p == q

    p = c1.load
    q = SingularityFunction(x, 5, -1)
    assert p == q

test_remove_load()


def test_reactions_point_loads():
    E, A = symbols('E A')

    # Test for one fixed support at x=0
    c = Column(5, E, A)
    c.apply_support(0)
    c.apply_load(-1, 2.5, -1)
    c.apply_load(2, 5, -1)

    c.solve_for_reaction_loads()
    p = c.reaction_loads
    R_0 = Symbol('R_0')
    q = {R_0: -1}
    assert p == q

    # Test for one fixed support at x=L
    c1 = Column(5, E, A)
    c1.apply_support(5)
    c1.apply_load(-4, 2.5, -1)
    c1.apply_load(2, 5, -1)

    c1.solve_for_reaction_loads()
    p = c1.reaction_loads
    R_5 = Symbol('R_5')
    q = {R_5: 2}
    assert p == q

    # Test for two supports at ends
    c2 = Column(10, E, A)
    c2.apply_support(10) # Check if order of applying supports matters
    c2.apply_support(0)
    c2.apply_load(-1, 5, -1)

    c2.solve_for_reaction_loads()
    p = c2.reaction_loads
    R_10 = Symbol('R_10')
    q = {R_0: Rational(1,2), R_10: Rational(1,2)}
    assert p == q

    # Test for two supports, not at ends
    c3 = Column(10, E, A)
    c3.apply_support(2)
    c3.apply_support(8)
    c3.apply_load(-1, 5, -1)
    c3.apply_load(-1, 0, -1)
    c3.apply_load(-1, 10, -1)

    c3.solve_for_reaction_loads()
    p = c3.reaction_loads
    R_2, R_8 = symbols('R_2 R_8')
    q = {R_2: Rational(3,2), R_8: Rational(3,2)}
    assert p == q

    # Test for two supports at ends, unsymmetrical load
    c4 = Column(10, E, A)
    c4.apply_support(0)
    c4.apply_support(10)
    c4.apply_load(-1, 6, -1)

    c4.solve_for_reaction_loads()
    p = c4.reaction_loads
    q = {R_0: Rational(2,5), R_10: Rational(3,5)}
    assert p == q

    # Test for two supports, one at end and one not
    c5 = Column(10, E, A)
    c5.apply_support(10)
    c5.apply_support(2)
    c5.apply_load(-1, 0, -1)
    c5.apply_load(2, 5, -1)
    c5.apply_load(-3, 8, -1)

    c5.solve_for_reaction_loads()
    p = c5.reaction_loads
    q = {R_2: Rational(1,2), R_10: Rational(3,2)}
    assert p == q

    # Test for three supports
    c6 = Column(10, E, A)
    c6.apply_support(0)
    c6.apply_support(5)
    c6.apply_support(10)
    c6.apply_load(-1, 2, -1)
    c6.apply_load(-1, 8, -1)

    c6.solve_for_reaction_loads()
    p = c6.reaction_loads
    q = {R_0: Rational(3,5), R_5: Rational(4,5), R_10: Rational(3,5)}
    assert p == q

    # Test for three supports, force on one member
    c6 = Column(10, E, A)
    c6.apply_support(0)
    c6.apply_support(5)
    c6.apply_support(10)
    c6.apply_load(-1, 2, -1)

    c6.solve_for_reaction_loads()
    p = c6.reaction_loads
    q = {R_0: Rational(3,5), R_5: Rational(2,5), R_10: 0}
    assert p == q

    # Test using symbols
    L, F = symbols('L F', positive=True)

    c7 = Column(L, E, A)
    c7.apply_support(0)
    c7.apply_support(L)
    c7.apply_load(-F, L/2, -1)

    c7.solve_for_reaction_loads()
    p = c7.reaction_loads
    R_L = Symbol('R_L')
    q = {R_0: F/2, R_L: F/2}
    assert p == q

    # Test removing a load
    c8 = Column(10, E, A)
    c8.apply_support(0)
    c8.apply_support(10)
    c8.apply_load(-1, 0, -1)
    c8.apply_load(-1, 5, -1)
    c8.apply_load(-1, 10, -1)

    c8.solve_for_reaction_loads()
    p = c8.reaction_loads
    q = {R_0: Rational(3,2), R_10: Rational(3,2)}
    assert p == q

    c8.remove_load(-1, 10, -1)

    c8.solve_for_reaction_loads()
    p = c8.reaction_loads
    q = {R_0: Rational(3,2), R_10: Rational(1,2)}
    assert p == q

test_reactions_point_loads()


def test_reactions_higher_orders():
    E, A = symbols('E A')

    # Test UDE, one support
    c = Column(10, E, A)
    c.apply_support(0)
    c.apply_load(-1, 0, 0, end=10)

    c.solve_for_reaction_loads()
    p = c.reaction_loads
    R_0 = Symbol('R_0')
    q = {R_0: 10}
    assert p == q

    # Test UDE, two supports
    c1 = Column(10, E, A)
    c1.apply_support(0)
    c1.apply_support(10)
    c1.apply_load(1, 0, 0, end=5)
    c1.apply_load(2, 5, 0, end=10)

    c1.solve_for_reaction_loads()
    p = c1.reaction_loads
    R_10 = Symbol('R_10')
    q = {R_0: -Rational(25,4), R_10: -Rational(35,4)}
    assert p == q

    # Test ramp load, one support
    c2 = Column(10, E, A)
    c2.apply_support(10)
    c2.apply_load(1, 0, 1, end=5)

    c2.solve_for_reaction_loads()
    p = c2.reaction_loads
    q = {R_10: -Rational(25,2)}
    assert p == q

    # Test ramp loads, two supports
    c3 = Column(10, E, A)
    c3.apply_support(0)
    c3.apply_support(10)
    c3.apply_load(1, 0, 1, end=10)

    c3.solve_for_reaction_loads()
    p = c3.reaction_loads
    q = {R_0: -Rational(50,3), R_10: -Rational(100,3)}
    assert p == q

    # Test parabolic load
    c4 = Column(10, E, A)
    c4.apply_support(10)
    c4.apply_load(1, 0, 2, end=10)

    c4.solve_for_reaction_loads()
    p = c4.reaction_loads
    q = {R_10: -Rational(1000,3)}
    assert p == q

    # Test combination of loads
    c5 = Column(10, E, A)
    c5.apply_support(0)
    c5.apply_support(10)

    c5.apply_load(100, 0, -1)
    c5.apply_load(20, 0, 0, end=4)
    c5.apply_load(20, 4, 0, end=8)
    c5.apply_load(5, 4, 1, end=8)
    c5.apply_load(40, 8, 0, end=10)
    c5.apply_load(2, 8, 2, end=10)

    c5.solve_for_reaction_loads()
    p = c5.reaction_loads
    q = {R_0: -Rational(1088,5), R_10: -Rational(2516,15)}
    assert p == q

    # Test ramp load in opposite direction
    c6 = Column(8, 20000, 0.75)
    c6.apply_support(0)
    c6.apply_support(8)
    c6.apply_load(-100, 8, -1)
    c6.apply_load(-20, 0, 0, end = 8)
    c6.apply_load(-10, 4, 1, end = 0) # Ramp load, starts at x = 4
    c6.solve_for_reaction_loads()

    p = c6.load
    R_8 = Symbol('R_8')
    q = (R_0*SingularityFunction(x, 0, -1)
         + R_8*SingularityFunction(x, 8, -1)
         - 60*SingularityFunction(x, 0, 0)
         + 10*SingularityFunction(x, 0, 1)
         - 10*SingularityFunction(x, 4, 1)
         - 100*SingularityFunction(x, 8, -1)
         + 20*SingularityFunction(x, 8, 0))
    assert p == q

    p = c6.reaction_loads
    q = {R_0: Rational(440,3), R_10: Rational(580,3)}

test_reactions_higher_orders()


def test_telescope_hinge():
    E, A = symbols('E A')
    c = Column(10, E, A)
    c.apply_support(0)
    c.apply_support(10)
    c.apply_load(10, 5, -1)
    c.apply_telescope_hinge(7.5)

    # Test boundary conditions
    p = c._bc_hinge
    q = [7.5]
    assert p == q

    p = c._applied_hinges
    q = [Symbol('u_7.5')]
    assert p == q

    # Test load equations telescope hinge
    p = c.load
    R_0, R_10 = symbols('R_0, R_10')
    q = (
        R_0 * SingularityFunction(x, 0, -1) +
        10 * SingularityFunction(x, 5, -1) +
        E*A*Symbol('u_7.5') * SingularityFunction(x, 7.5, -2) +
        R_10 * SingularityFunction(x, 10, -1)
    )
    assert p == q

    # Test solution single telescope hinge
    c.solve_for_reaction_loads()
    p = c.reaction_loads
    q = {R_0: -10, R_10: 0}
    assert p == q

    p = c.hinge_extensions
    q = {Symbol('u_7.5'): 50/(E*A)}
    assert p == q

    # Test numeric solution, multiple forces
    c2 = Column(10, 20000, 0.5)
    c2.apply_support(0)
    c2.apply_support(10)
    c2.apply_telescope_hinge(4)
    c2.apply_load(-5, 3, -1)
    c2.apply_load(-10, 8, -1)

    c2.solve_for_reaction_loads()

    p = c2.reaction_loads
    q = {R_0: 5, R_10: 10}
    assert p == q

    p = c2.hinge_extensions
    u_4 = Symbol('u_4')
    q = {u_4: Rational(1, 2000)}
    assert p == q

test_telescope_hinge()


def test_equations():
    c = Column(10, 210000, 1)
    c.apply_support(0)
    c.apply_support(10)
    c.apply_load(5, 8, -1)
    R_0, R_10 = symbols("R_0 R_10")
    C_N, C_u = symbols("C_N C_u")

    # Test before solving the unkowns
    p = c.axial_force()
    q = (
        C_N
        - R_0*SingularityFunction(x, 0, 0)
        - R_10*SingularityFunction(x, 10, 0)
        - 5*SingularityFunction(x, 8, 0)
    )
    assert p == q

    p = c.extension()
    q = (
        C_N*x + C_u
        - R_0*SingularityFunction(x, 0, 1)/210000
        - R_10*SingularityFunction(x, 10, 1)/210000
        - SingularityFunction(x, 8, 1)/42000
    )
    assert p == q

    # Test after solving the unknowns
    c.solve_for_reaction_loads()

    p = c.axial_force()
    q = (
        SingularityFunction(x, 0, 0)
        - 5*SingularityFunction(x, 8, 0)
        + 4*SingularityFunction(x, 10, 0)
    )
    assert p == q

    p = c.extension()
    q = (
        SingularityFunction(x, 0, 1)/210000
        - SingularityFunction(x, 8, 1)/42000
        + SingularityFunction(x, 10, 1)/52500
    )
    assert p == q

test_equations()
