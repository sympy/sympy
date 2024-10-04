from sympy.physics.continuum_mechanics.structure2d import Structure2d
from sympy.core.function import expand
from sympy.core.numbers import (Rational, pi)
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.sets.sets import Interval
from sympy.simplify.simplify import simplify
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.functions import SingularityFunction, Piecewise, meijerg, Abs, log
from sympy.testing.pytest import raises
from sympy.physics.units import meter, newton, kilo, giga, milli
from sympy.physics.continuum_mechanics.beam import Beam3D
from sympy.geometry import Circle, Polygon, Point2D, Triangle
from sympy.core.sympify import sympify
from sympy.functions.elementary.trigonometric import sin, cos, atan2

x = Symbol('x')
y = Symbol('y')
R1, R2 = symbols('R1, R2')

def test_Structure2d_symbolic_1():
    ## symbolic test
    E = Symbol('E')
    I = Symbol('I')
    A = Symbol('A')

    s = Structure2d()

    s.add_member(0,0,2,0,E, I, A)
    s.add_member(2,0,4,0,E, I, A)

    assert s.members[0].x1 == 0
    assert s.members[0].y1 == 0
    assert s.members[0].x2 == 2
    assert s.members[0].y2 == 0
    assert s.members[0].E == E
    assert s.members[0].I == I
    assert s.members[0].A == A
    assert s.members[0].member_id == 0
    assert s.members[0].length == 2
    assert s.members[0].angle == 0
    assert s.members[0].angle_deg == 0
    assert s.members[0].member_loads == []
    assert s.members[0].member_eq == 0

    assert s.members[1].x1 == 2
    assert s.members[1].y1 == 0
    assert s.members[1].x2 == 4
    assert s.members[1].y2 == 0
    assert s.members[1].E == E
    assert s.members[1].I == I
    assert s.members[1].A == A
    assert s.members[1].member_id == 1
    assert s.members[1].length == 2
    assert s.members[1].angle == 0
    assert s.members[1].angle_deg == 0
    assert s.members[1].member_loads == []
    assert s.members[1].member_eq == 0

    Rv1,Rh1 = symbols('Rv1 Rh1')
    Rv1,Rh1 = s.apply_support(x=0,y=0,type='pin')

    Rv2 = symbols('Rv2')
    Rv2 = s.apply_support(x=4,y=0,type='roller')

    #check pin support
    assert s.supports[0].x == 0
    assert s.supports[0].y == 0
    assert s.supports[0].node_type == 'pin'

    order_explain = "order should be -1 for point loads -2 for moments and 0 for distributed loads"

    debug_str = f"expected 0 got {s.nodes[0].node_loads[1].global_angle}"

    assert s.nodes[0].node_loads[0].start_x == 0
    assert s.nodes[0].node_loads[0].start_y == 0
    assert s.nodes[0].node_loads[0].value == Symbol('R_v-0')
    assert s.nodes[0].node_loads[0].global_angle == pi*3/2
    assert s.nodes[0].node_loads[0].order == -1, order_explain
    assert s.nodes[0].node_loads[0].end_x == None
    assert s.nodes[0].node_loads[0].end_y == None
    assert s.nodes[0].node_loads[0].load_id == 0
    assert s.nodes[0].node_loads[0].x_component == 0
    assert s.nodes[0].node_loads[0].y_component == -1*Symbol('R_v-0')

    assert s.nodes[0].node_loads[1].start_x == 0
    assert s.nodes[0].node_loads[1].start_y == 0
    assert s.nodes[0].node_loads[1].value == Symbol('R_h-0')
    assert s.nodes[0].node_loads[1].global_angle == 0
    assert s.nodes[0].node_loads[1].order == -1, order_explain
    assert s.nodes[0].node_loads[1].end_x == None
    assert s.nodes[0].node_loads[1].end_y == None
    assert s.nodes[0].node_loads[1].load_id == 1
    assert s.nodes[0].node_loads[1].x_component == -1*Symbol('R_h-0')
    assert s.nodes[0].node_loads[1].y_component == 0

    # check that there are no other loads on the node
    assert len(s.nodes[0].node_loads) == 2

    #middle node should have no loads
    assert len(s.nodes[1].node_loads) == 0


    assert s.supports[1].x == 4
    assert s.supports[1].y == 0
    assert s.supports[1].node_type == 'roller'

    assert s.nodes[2].node_loads[0].start_x == 4
    assert s.nodes[2].node_loads[0].start_y == 0
    assert s.nodes[2].node_loads[0].value == Symbol('R_v-4')
    assert s.nodes[2].node_loads[0].global_angle == pi*3/2
    assert s.nodes[2].node_loads[0].order == -1, order_explain
    assert s.nodes[2].node_loads[0].end_x == None
    assert s.nodes[2].node_loads[0].end_y == None
    assert s.nodes[2].node_loads[0].load_id == 2
    assert s.nodes[2].node_loads[0].x_component == 0
    assert s.nodes[2].node_loads[0].y_component == -1*Symbol('R_v-4')

    # check that there are no other loads on the node
    assert len(s.nodes[2].node_loads) == 1

    #add loads
    s.apply_load(2,0,symbols('F'),global_angle=90,order=-1)

    assert s.nodes[1].node_loads[0].start_x == 2
    assert s.nodes[1].node_loads[0].start_y == 0
    assert s.nodes[1].node_loads[0].value == Symbol('F')
    assert s.nodes[1].node_loads[0].global_angle == pi/2
    assert s.nodes[1].node_loads[0].order == -1, order_explain
    assert s.nodes[1].node_loads[0].end_x == None
    assert s.nodes[1].node_loads[0].end_y == None
    assert s.nodes[1].node_loads[0].load_id == 3
    assert s.nodes[1].node_loads[0].x_component == 0
    assert s.nodes[1].node_loads[0].y_component == Symbol('F')

    s.solve_for_reaction_loads(Rv1, Rh1, Rv2)
    assert s.reaction_loads[Rv1] == Symbol('F') /2
    assert s.reaction_loads[Rh1] == 0
    assert s.reaction_loads[Rv2] == Symbol('F') / 2

    F = Symbol('F')
    Rv0 = Symbol('R_v-0')
    Rv4 = Symbol('R_v-4')

    assert s.load_qz == F*SingularityFunction(x, 2, -1) - Rv0*SingularityFunction(x, 0, -1) - Rv4*SingularityFunction(x, 4, -1)
    assert s.shear_force() == F*SingularityFunction(x, 0, 0)/2 - F*SingularityFunction(x, 2, 0) + F*SingularityFunction(x, 4, 0)/2
    assert s.bending_moment() == F*SingularityFunction(x, 0, 1)/2 - F*SingularityFunction(x, 2, 1) + F*SingularityFunction(x, 4, 1)/2

    print("Structure2d symbolic test 1: passed.")

def test_Structure2d_symbolic_2():
    ## symbolic test
    E = Symbol('E')
    I = Symbol('I')
    A = Symbol('A')

    s = Structure2d()

    s.add_member(0,0,2,0,E, I, A)
    s.add_member(2,0,4,0,E, I, A)

    Rv1,Rh1 = symbols('Rv1 Rh1')
    Rv1,Rh1 = s.apply_support(x=0,y=0,type='pin')

    Rv2 = symbols('Rv2')
    Rv2 = s.apply_support(x=4,y=0,type='roller')

    s.apply_load(2,0,symbols('F'),global_angle=45,order=-1)

    assert s.nodes[1].node_loads[0].start_x == 2
    assert s.nodes[1].node_loads[0].start_y == 0
    assert s.nodes[1].node_loads[0].value == Symbol('F')
    assert s.nodes[1].node_loads[0].global_angle == pi/4
    assert s.nodes[1].node_loads[0].order == -1
    assert s.nodes[1].node_loads[0].end_x == None
    assert s.nodes[1].node_loads[0].end_y == None
    assert s.nodes[1].node_loads[0].load_id == 3
    assert s.nodes[1].node_loads[0].x_component == Symbol('F')*cos(pi/4) * -1 #to the right is positive
    assert s.nodes[1].node_loads[0].y_component == Symbol('F')*sin(pi/4)

    s.solve_for_reaction_loads(Rv1, Rh1, Rv2)
    # assert s.reaction_loads[Rv1] == Symbol('F') /2
    # assert s.reaction_loads[Rh1] == 0
    # assert s.reaction_loads[Rv2] == Symbol('F') /2

    # F = Symbol('F')
    # Rv0 = Symbol('R_v-0')
    # Rv4 = Symbol('R_v-4')


    print("Structure2d symbolic test 2: passed.")




test_Structure2d_symbolic_2()

