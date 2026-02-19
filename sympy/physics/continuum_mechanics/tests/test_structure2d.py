from sympy.physics.continuum_mechanics.structure2d import Structure2d
from sympy.core.symbol import Symbol, symbols
from sympy.core.numbers import Rational
from sympy.functions import SingularityFunction
from sympy.simplify.simplify import simplify
# from sympy import pi, cos, sin, rad, sympify
import pytest
import math

x = Symbol('x')
F = Symbol('F')
# Rv_0, Rh_0, Rv_4 = symbols('R_v__0,__0 R_h__0,__0 R_v__4,__0')
Rv_0 = Symbol('R_v (x=0,y=0)')
Rh_0 = Symbol('R_h (x=0,y=0)')
Rv_4 = Symbol('R_v (x=4,y=0)')



# Test Data Definitions
symbolic_test_data = [
    {
        'id': 'symbolic_test_1',
        'angle': 270,
        'value': F,
        'reaction_loads': {Rv_0: -F / 2, Rv_4: -F / 2, Rh_0: 0},
        'bending_moment': {'symbolic': F * SingularityFunction(x, 0, 1) / 2 - F * SingularityFunction(x, 2, 1) + F * SingularityFunction(x, 4, 1) / 2}
    },
]


@pytest.fixture
def test_beam_fixture():
    E, I, A = symbols('E I A')
    s = Structure2d()
    s.add_member(0, 0, 4, 0, E, I, A)
    return s, E, I, A

@pytest.mark.parametrize(
    'test_data',
    symbolic_test_data,
    ids=[
        data['id'] if isinstance(data['id'], (str, float, int, type(None))) else str(data['id'])
        for data in symbolic_test_data
    ]
)

def test_structure2d_symbolic(test_beam_fixture, test_data):
    """Test symbolic calculations for a horizontal beam with a point load at the center."""

    # Extract test data
    s, E, I, A = test_beam_fixture
    angle = test_data['angle']
    value = test_data['value']
    reaction_loads = test_data['reaction_loads']
    bending_moment = test_data['bending_moment']

    # Assert member properties
    member = s.members[0]
    assert (member.x1, member.y1, member.x2, member.y2) == (0, 0, 4, 0)
    assert member.E == E
    assert member.I == I
    assert member.A == A
    assert member.length == 4
    assert member.angle == 0

    # Apply supports
    s.apply_support(x=0, y=0, type='pin')
    s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2, 0, value, global_angle=angle, order=-1)
    s.solve_for_reaction_loads()

    # Check reaction loads
    for load_symbol, expected_value in reaction_loads.items():
        assert simplify(s.reaction_loads[load_symbol] - expected_value) == 0

    # Check bending moment
    assert s.bending_moment() == bending_moment['symbolic']
#################################################################################################

numerical_test_data = [
    {
        'id': 'basic_numerical_test_1',
        'angle': 0,
        'value': 10,
        'reaction_loads': {Rh_0: -10, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_2',
        'angle': 30,
        'value': 10,
        'reaction_loads': {Rh_0: -8.66, Rv_0: 2.50, Rv_4: 2.50},
        'bending_moment': {0: 0, 2: -5, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_3',
        'angle': 60,
        'value': 10,
        'reaction_loads': {Rh_0: -5, Rv_0: 4.33, Rv_4: 4.33},
        'bending_moment': {0: 0, 2: -8.66, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_4',
        'angle': 90,
        'value': 10,
        'reaction_loads': {Rh_0: 0, Rv_0: 5, Rv_4: 5},
        'bending_moment': {0: 0, 2: -10, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_5',
        'angle': 120,
        'value': 10,
        'reaction_loads': {Rh_0: 5, Rv_0: 4.33, Rv_4: 4.33},
        'bending_moment': {0: 0, 2: -8.66, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_6',
        'angle': 150,
        'value': 10,
        'reaction_loads': {Rh_0: 8.66, Rv_0: 2.5, Rv_4: 2.5},
        'bending_moment': {0: 0, 2: -5, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_7',
        'angle': 180,
        'value': 10,
        'reaction_loads': {Rh_0: 10, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_8',
        'angle': 210,
        'value': 10,
        'reaction_loads': {Rh_0: 8.66, Rv_0: -2.5, Rv_4: -2.5},
        'bending_moment': {0: 0, 2: 5, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_9',
        'angle': 240,
        'value': 10,
        'reaction_loads': {Rh_0: 5, Rv_0: -4.33, Rv_4: -4.33},
        'bending_moment': {0: 0, 2: 8.66, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_10',
        'angle': 270,
        'value': 10,
        'reaction_loads': {Rh_0: 0, Rv_0: -5, Rv_4: -5},
        'bending_moment': {0: 0, 2: 10, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_11',
        'angle': 300,
        'value': 10,
        'reaction_loads': {Rh_0: -5, Rv_0: -4.33, Rv_4: -4.33},
        'bending_moment': {0: 0, 2: 8.66, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_12',
        'angle': 330,
        'value': 10,
        'reaction_loads': {Rh_0: -8.66, Rv_0: -2.5, Rv_4: -2.5},
        'bending_moment': {0: 0, 2: 5, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_13',
        'angle': 360,
        'value': 10,
        'reaction_loads': {Rh_0: -10, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
]

@pytest.mark.parametrize(
    'test_data',
    numerical_test_data,
    ids=[
        data['id'] if isinstance(data['id'], (str, float, int, type(None))) else str(data['id'])
        for data in numerical_test_data
    ]
)

def test_numerical_pointload(test_beam_fixture, test_data):
    """Test numerical calculations for a horizontal beam with a point load at the center. angles are 0-360 in 30 degree increments"""
    # Extract test data
    s, E, I, A = test_beam_fixture
    angle = test_data['angle']
    value = test_data['value']
    reaction_loads = test_data['reaction_loads']
    bending_moment = test_data['bending_moment']

    # Assert member properties
    member = s.members[0]
    assert (member.x1, member.y1, member.x2, member.y2) == (0, 0, 4, 0)
    assert member.E == E
    assert member.I == I
    assert member.A == A
    assert member.length == 4
    assert member.angle == 0

    # Apply supports
    s.apply_support(x=0, y=0, type='pin')
    s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2, 0, value, global_angle=angle, order=-1)
    s.solve_for_reaction_loads()

    # Check reaction loads SSSADASDAS
    for load_symbol, expected_value in reaction_loads.items():
        float_num = float(s.reaction_loads[load_symbol]) - float(expected_value)
        # assert math.isclose(float_num,0) == True
        assert math.isclose(float_num,0,rel_tol=1e-2,abs_tol=1e-2) == True, f"Expected: {expected_value}, Actual: {s.reaction_loads[load_symbol]}"

    # Check bending moment
    for loc, expected_bm in bending_moment.items():
        float_bm = float(s.bending_moment(loc)) - float(expected_bm)
        assert math.isclose(float_bm,0,rel_tol=1e-2,abs_tol=1e-2) == True
#####################################################################################################################

numerical_test_data_2 = [
    {
        'id': 'basic_numerical_test_1',
        'angle': 0,
        'value': 15,
        'reaction_loads': {Rh_0: -15, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_2',
        'angle': 30,
        'value': 15,
        'reaction_loads': {Rh_0: -12.99, Rv_0: 2.81, Rv_4: 4.69},
        'bending_moment': {0: 0, 2: -5.63, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_3',
        'angle': 60,
        'value': 15,
        'reaction_loads': {Rh_0: -7.5, Rv_0: 4.87, Rv_4: 8.12},
        'bending_moment': {0: 0, 2: -9.74, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_4',
        'angle': 90,
        'value': 15,
        'reaction_loads': {Rh_0: 0, Rv_0: 5.63, Rv_4: 9.38},
        'bending_moment': {0: 0, 2: -11.25, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_5',
        'angle': 120,
        'value': 15,
        'reaction_loads': {Rh_0: 7.5, Rv_0: 4.87, Rv_4: 8.12},
        'bending_moment': {0: 0, 2: -9.74, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_6',
        'angle': 150,
        'value': 15,
        'reaction_loads': {Rh_0: 12.99, Rv_0: 2.81, Rv_4: 4.69},
        'bending_moment': {0: 0, 2: -5.63, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_7',
        'angle': 180,
        'value': 15,
        'reaction_loads': {Rh_0: 15, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_8',
        'angle': 210,
        'value': 15,
        'reaction_loads': {Rh_0: 12.99, Rv_0: -2.81, Rv_4: -4.69},
        'bending_moment': {0: 0, 2: 5.63, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_9',
        'angle': 240,
        'value': 15,
        'reaction_loads': {Rh_0: 7.5, Rv_0: -4.87, Rv_4: -8.12},
        'bending_moment': {0: 0, 2: 9.74, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_10',
        'angle': 270,
        'value': 15,
        'reaction_loads': {Rh_0: 0, Rv_0: -5.63, Rv_4: -9.38},
        'bending_moment': {0: 0, 2: 11.25, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_11',
        'angle': 300,
        'value': 15,
        'reaction_loads': {Rh_0: -7.5, Rv_0: -4.87, Rv_4: -8.12},
        'bending_moment': {0: 0, 2: 9.74, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_12',
        'angle': 330,
        'value': 15,
        'reaction_loads': {Rh_0: -12.99, Rv_0: -2.81, Rv_4: -4.69},
        'bending_moment': {0: 0, 2: 5.63, 4: 0} # location: value
    },
    {
        'id': 'basic_numerical_test_13',
        'angle': 360,
        'value': 15,
        'reaction_loads': {Rh_0: -15, Rv_0: 0, Rv_4: 0},
        'bending_moment': {0: 0, 2: 0, 4: 0} # location: value
    },
]

@pytest.mark.parametrize(
    'test_data',
    numerical_test_data_2,
    ids=[
        data['id'] if isinstance(data['id'], (str, float, int, type(None))) else str(data['id'])
        for data in numerical_test_data_2
    ]
)

def test_numerical_distload(test_beam_fixture, test_data):
    """Test numerical calculations for a horizontal beam with a distributed load at the center 1m distance from both ends. angles are 0-360 in 30 degree increments"""
    # Extract test data
    s, E, I, A = test_beam_fixture
    angle = test_data['angle']
    value = test_data['value']
    reaction_loads = test_data['reaction_loads']
    bending_moment = test_data['bending_moment']

    # Assert member properties
    member = s.members[0]
    assert (member.x1, member.y1, member.x2, member.y2) == (0, 0, 4, 0)
    assert member.E == E
    assert member.I == I
    assert member.A == A
    assert member.length == 4
    assert member.angle == 0

    # Apply supports
    s.apply_support(x=0, y=0, type='pin')
    s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2,0,value,global_angle=angle,order=0,end_x=3,end_y=0)
    s.solve_for_reaction_loads()

    # Check reaction loads SSSADASDAS
    for load_symbol, expected_value in reaction_loads.items():
        float_num = float(s.reaction_loads[load_symbol]) - float(expected_value)
        # assert math.isclose(float_num,0) == True
        # print(reaction_loads)
        assert math.isclose(float_num,0,rel_tol=1e-2,abs_tol=1e-2) == True, f"Expected: {expected_value}, Actual: {s.reaction_loads[load_symbol]}"

    # Check bending moment
    for loc, expected_bm in bending_moment.items():
        float_bm = float(s.bending_moment(loc)) - float(expected_bm)
        assert math.isclose(float_bm,0,rel_tol=1e-2,abs_tol=1e-2) == True
#####################################################################################################################


# Test Data Definitions
# Rv_7, Rh_7, T_7 = symbols('R_v__7,__-1 R_h__7,__-1 T__7,__-1')
Rv_7 = Symbol('R_v (x=7,y=-1)')
Rh_7 = Symbol('R_h (x=7,y=-1)')
T_7 = Symbol('T (x=7,y=-1)')

numerical_test_data = [
    {
        'id': 'symbolic_test_1',
        'angle': 0,
        'value': 15,
        'reaction_loads': {Rv_7: -75, Rh_7: 3.75, T_7: -323.44},
        'bending_moment': {0: 0, 5: -120, 7: -171.54},

    },
]


@pytest.fixture
def test_beam_fixture_one_bend():
    E = 3e4
    I = 1
    A = 1e4
    s = Structure2d()
    s.add_member(0,0,3,4,E, I, A)
    s.add_member(3, 4, 7, -1, E, I, A)
    return s

@pytest.mark.parametrize(
    'test_data',
    numerical_test_data,
    ids=[
        data['id'] if isinstance(data['id'], (str, float, int, type(None))) else str(data['id'])
        for data in numerical_test_data
    ]
)
def test_structure2d_symbolic_onebend(test_beam_fixture_one_bend, test_data):
    """This one has a bend in the middle"""

    # Extract test data
    s = test_beam_fixture_one_bend
    # angle = test_data['angle']
    value = test_data['value']

    reaction_loads = test_data['reaction_loads']
    bending_moment = test_data['bending_moment']


    s.apply_load(start_x=1.5,start_y=2,value=value,global_angle=0,order=-1)
    s.apply_support(x=7, y=-1, type='fixed')

    s.apply_load(5,1.5,value/2,global_angle=s.members[1].angle_deg + 270,order=0,end_x=7,end_y=-1)
    s.apply_load(0,0,value*0.8,global_angle=270,order=0,end_x=3,end_y=4)
    s.solve_for_reaction_loads()
    # s.draw()
    # Check reaction loads
    for load_symbol, expected_value in reaction_loads.items():
        print(float(s.reaction_loads[load_symbol]), float(expected_value))
        # print(math.isclose(float(s.reaction_loads[load_symbol]),float(expected_value),rel_tol=1e-2,abs_tol=1e-2))
        assert math.isclose(float(s.reaction_loads[load_symbol]),float(expected_value),rel_tol=1e-2,abs_tol=1e-2) == True , f"Expected: {float(expected_value)}, Actual: {float(s.reaction_loads[load_symbol])}"

    # Check bending moment
    for loc, expected_bm in bending_moment.items():
        # float_bm = float(s.bending_moment(loc)) - float(expected_bm)
        print(float(s.bending_moment(loc)), float(expected_bm))
        print(math.isclose(float(s.bending_moment(loc)),float(expected_bm),rel_tol=1e-2,abs_tol=1e-2))

        assert math.isclose(float(s.bending_moment(loc)),float(expected_bm),rel_tol=1e-2,abs_tol=1e-2) == True


    """
    All test cases written in test_frame_problems(),
        test_distributed_loads(), test_hinge_problems(),
        test_joint_loads(), test_moment_loads()
        are tested and verified using Alex algorithm
        from the thesis https://oit.tudelft.nl/Macaulays-method/theses/alex.html
        and all validated testcases can be found here for reference https://github.com/saiudayagiri/GSoC/blob/main/Test%20cases%20using%20alex%20algorithm.ipynb
     """
def test_frame_problems():
    # Alex Example 1(https://oit.tudelft.nl/Macaulays-method/theses/Macaulay-2D/BEP%20v1.html)
    # 3 member structure with mixed loads
    # one fixed support
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 15, 0, -1)
    s.apply_load(2, 0, 16, 270, -1)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 6, 1.5)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -55

    assert result[Symbol('R_h (x=11,y=-1)')] == -15

    assert result[Symbol('T (x=11,y=-1)')] == -435

    expected_axial_force = -15*SingularityFunction(x, 0, 0) + 27*SingularityFunction(x, 4, 0) + 18*SingularityFunction(x, 4, 1)/5 - 18*SingularityFunction(x, Rational(13,2), 1)/5 - 74*SingularityFunction(x, 9, 0) + 53*SingularityFunction(x, 14, 0)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) - 16*SingularityFunction(x, 2, 0) - SingularityFunction(x, 4, 0) + 6*SingularityFunction(x, 4, 1)/5 + 24*SingularityFunction(x, Rational(13,2), 1)/5 + 32*SingularityFunction(x, 9, 0) + 435*SingularityFunction(x, 14, -1) + 21*SingularityFunction(x, 14, 0)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) - 16*SingularityFunction(x, 2, 1) - SingularityFunction(x, 4, 1) + 3*SingularityFunction(x, 4, 2)/5 + 12*SingularityFunction(x, Rational(13,2), 2)/5 + 32*SingularityFunction(x, 9, 1) + 435*SingularityFunction(x, 14, 0) + 21*SingularityFunction(x, 14, 1)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 15, 270, -1)
    s.apply_load(2, 0, 16, 0, -1)
    s.apply_load(0, 0, 6, 0, 0, 4, 0)
    s.apply_load(4, 0, 6, 0, 0, 6, 1.5)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -15

    assert result[Symbol('R_h (x=11,y=-1)')] == -55

    assert result[Symbol('T (x=11,y=-1)')] == Rational(-9875,100)

    expected_axial_force = -6*SingularityFunction(x, 0, 1) - 16*SingularityFunction(x, 2, 0) + 17*SingularityFunction(x, 4, 0) + 6*SingularityFunction(x, 4, 1)/5 + 24*SingularityFunction(x, Rational(13,2), 1)/5 - 10*SingularityFunction(x, 9, 0) + 45*SingularityFunction(x, 14, 0)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -15*SingularityFunction(x, 0, 0) - 21*SingularityFunction(x, 4, 0) - 18*SingularityFunction(x, 4, 1)/5 + 18*SingularityFunction(x, Rational(13,2), 1)/5 + 80*SingularityFunction(x, 9, 0) + 395*SingularityFunction(x, 14, -1)/4 - 35*SingularityFunction(x, 14, 0)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -15*SingularityFunction(x, 0, 1) - 21*SingularityFunction(x, 4, 1) - 9*SingularityFunction(x, 4, 2)/5 + 9*SingularityFunction(x, Rational(13,2), 2)/5 + 80*SingularityFunction(x, 9, 1) + 395*SingularityFunction(x, 14, 0)/4 - 35*SingularityFunction(x, 14, 1)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(4, 0, 15, 0, -1)
    s.apply_load(0, 0, 16, 270, -1)
    s.apply_load(4, 0, 6, 270, 0, 6, 1.5)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -31

    assert result[Symbol('R_h (x=11,y=-1)')] == -15

    assert result[Symbol('T (x=11,y=-1)')] == -251

    expected_axial_force = -12*SingularityFunction(x, 4, 0)/5 + 18*SingularityFunction(x, 4, 1)/5 - 18*SingularityFunction(x, Rational(13,2), 1)/5 - 202*SingularityFunction(x, 9, 0)/5 + 169*SingularityFunction(x, 14, 0)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -16*SingularityFunction(x, 0, 0) - 29*SingularityFunction(x, 4, 0)/5 - 24*SingularityFunction(x, 4, 1)/5 + 24*SingularityFunction(x, Rational(13,2), 1)/5 + 136*SingularityFunction(x, 9, 0)/5 + 251*SingularityFunction(x, 14, -1) + 33*SingularityFunction(x, 14, 0)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -16*SingularityFunction(x, 0, 1) - 29*SingularityFunction(x, 4, 1)/5 - 12*SingularityFunction(x, 4, 2)/5 + 12*SingularityFunction(x, Rational(13,2), 2)/5 + 136*SingularityFunction(x, 9, 1)/5 + 251*SingularityFunction(x, 14, 0) + 33*SingularityFunction(x, 14, 1)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(4, 0, 15, 0, -1)
    s.apply_load(0, 0, 16, 270, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -16

    assert result[Symbol('R_h (x=11,y=-1)')] == -15

    assert result[Symbol('T (x=11,y=-1)')] == -161

    expected_axial_force = -12*SingularityFunction(x, 4, 0)/5 - 97*SingularityFunction(x, 9, 0)/5 + 109*SingularityFunction(x, 14, 0)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -16*SingularityFunction(x, 0, 0) - 29*SingularityFunction(x, 4, 0)/5 + 121*SingularityFunction(x, 9, 0)/5 + 161*SingularityFunction(x, 14, -1) - 12*SingularityFunction(x, 14, 0)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -16*SingularityFunction(x, 0, 1) - 29*SingularityFunction(x, 4, 1)/5 + 121*SingularityFunction(x, 9, 1)/5 + 161*SingularityFunction(x, 14, 0) - 12*SingularityFunction(x, 14, 1)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # Alex Example 2(https://oit.tudelft.nl/Macaulays-method/theses/Macaulay-2D/BEP%20v2.html)
    # without the hinge
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 0, 0, 3, 4)
    s.apply_load(3, 4, 60, 0, 0, 6, 0)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 200

    assert result[Symbol('R_h (x=0,y=0)')] == -300

    assert result[Symbol('R_v (x=6,y=0)')] == -200

    assert result[Symbol('R_h (x=6,y=0)')] == -300

    expected_axial_force = 340*SingularityFunction(x, 0, 0) - 36*SingularityFunction(x, 0, 1) - 320*SingularityFunction(x, 5, 0) + 340*SingularityFunction(x, 10, 0) + 36*SingularityFunction(x, 10, 1)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 120*SingularityFunction(x, 0, 0) - 48*SingularityFunction(x, 0, 1) + 96*SingularityFunction(x, 5, 1) - 120*SingularityFunction(x, 10, 0) - 48*SingularityFunction(x, 10, 1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 120*SingularityFunction(x, 0, 1) - 24*SingularityFunction(x, 0, 2) + 48*SingularityFunction(x, 5, 2) - 120*SingularityFunction(x, 10, 1) - 24*SingularityFunction(x, 10, 2)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(3, 4, 60, 0, -1)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 40

    assert result[Symbol('R_h (x=0,y=0)')] == -30

    assert result[Symbol('R_v (x=6,y=0)')] == -40

    assert result[Symbol('R_h (x=6,y=0)')] == -30

    expected_axial_force = 50*SingularityFunction(x, 0, 0) - 100*SingularityFunction(x, 5, 0) + 50*SingularityFunction(x, 10, 0)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 0
    assert simplify(s.V.expand()) == simplify(expected_shear_force)

    expected_moment = 0
    assert simplify(s.M.expand()) == simplify(expected_moment)

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, 10**4, 10**4, 10**4)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(3, 4, 60, 270, -1)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == -30

    assert result[Symbol('R_h (x=0,y=0)')] == Rational(7920,427)

    assert result[Symbol('R_v (x=6,y=0)')] == -30

    assert result[Symbol('R_h (x=6,y=0)')] == Rational(-7920,427)

    expected_axial_force = -15000*SingularityFunction(x, 0, 0)/427 + 15000*SingularityFunction(x, 10, 0)/427
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 1350*SingularityFunction(x, 0, 0)/427 - 2700*SingularityFunction(x, 5, 0)/427 + 1350*SingularityFunction(x, 10, 0)/427
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 1350*SingularityFunction(x, 0, 1)/427 - 2700*SingularityFunction(x, 5, 1)/427 + 1350*SingularityFunction(x, 10, 1)/427
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 270, -1)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == -60

    assert result[Symbol('R_h (x=0,y=0)')] == 0

    assert result[Symbol('R_v (x=6,y=0)')] == 0

    assert result[Symbol('R_h (x=6,y=0)')] == 0

    expected_axial_force = 0
    assert simplify(s.N.expand()) == simplify(expected_axial_force)

    expected_shear_force = 0
    assert simplify(s.V.expand()) == simplify(expected_shear_force)

    expected_moment = 0
    assert simplify(s.M.expand()) == simplify(expected_moment)

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 0, -1)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 0

    assert result[Symbol('R_h (x=0,y=0)')] == -60

    assert result[Symbol('R_v (x=6,y=0)')] == 0

    assert result[Symbol('R_h (x=6,y=0)')] == 0

    expected_axial_force = 0
    assert simplify(s.N.expand()) == simplify(expected_axial_force)

    expected_shear_force = 0
    assert simplify(s.V.expand()) == simplify(expected_shear_force)

    expected_moment = 0
    assert simplify(s.M.expand()) == simplify(expected_moment)

    # Alex Example 3(https://oit.tudelft.nl/Macaulays-method/theses/Macaulay-2D/BEP%20v3.html)
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
    s.add_member(3, 4, 7, 1, E, I, A)
    s.apply_load(3, 4, 60, 270, -1)
    s.apply_load(3, 4, 18, 0, 0, 7, 1)
    s.apply_support(0, 0,"fixed")
    s.apply_support(7, 1,"roller")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == Rational(-44421,1120)

    assert result[Symbol('R_h (x=0,y=0)')] == -90

    assert result[Symbol('R_v (x=7,y=1)')] == Rational(-22779,1120)

    assert result[Symbol('T (x=0,y=0)')] == Rational(42021,160)

    expected_axial_force = 31179*SingularityFunction(x, 0, 0)/1400 + 30021*SingularityFunction(x, 5, 0)/800 - 72*SingularityFunction(x, 5, 1)/5 + 68337*SingularityFunction(x, 10, 0)/5600 + 72*SingularityFunction(x, 10, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -42021*SingularityFunction(x, 0, -1)/160 + 536463*SingularityFunction(x, 0, 0)/5600 - 929979*SingularityFunction(x, 5, 0)/5600 + 54*SingularityFunction(x, 5, 1)/5 + 22779*SingularityFunction(x, 10, 0)/1400 - 54*SingularityFunction(x, 10, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -42021*SingularityFunction(x, 0, 0)/160 + 536463*SingularityFunction(x, 0, 1)/5600 - 929979*SingularityFunction(x, 5, 1)/5600 + 27*SingularityFunction(x, 5, 2)/5 + 22779*SingularityFunction(x, 10, 1)/1400 - 27*SingularityFunction(x, 10, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, 10**2, 10**2)
    s.add_member(3, 4, 7, 1, E, I, A)
    s.apply_load(3, 4, 18, 0, 0, 7, 1)
    s.apply_support(0, 0,"fixed")
    s.apply_support(7, 1,"roller")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == Rational(5427,1120)

    assert result[Symbol('R_h (x=0,y=0)')] == -90

    assert result[Symbol('R_v (x=7,y=1)')] == Rational(-5427,1120)

    assert result[Symbol('T (x=0,y=0)')] == Rational(30573,160)

    expected_axial_force = 81027*SingularityFunction(x, 0, 0)/1400 + 8973*SingularityFunction(x, 5, 0)/800 - 72*SingularityFunction(x, 5, 1)/5 + 16281*SingularityFunction(x, 10, 0)/5600 + 72*SingularityFunction(x, 10, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -30573*SingularityFunction(x, 0, -1)/160 + 386919*SingularityFunction(x, 0, 0)/5600 - 711027*SingularityFunction(x, 5, 0)/5600 + 54*SingularityFunction(x, 5, 1)/5 + 5427*SingularityFunction(x, 10, 0)/1400 - 54*SingularityFunction(x, 10, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -30573*SingularityFunction(x, 0, 0)/160 + 386919*SingularityFunction(x, 0, 1)/5600 - 711027*SingularityFunction(x, 5, 1)/5600 + 27*SingularityFunction(x, 5, 2)/5 + 5427*SingularityFunction(x, 10, 1)/1400 - 27*SingularityFunction(x, 10, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 15, 0, -1)
    s.apply_load(8, 3, 16, 270, -1)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 8, 3)
    s.apply_load(8, 3, 6, 270, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -100

    assert result[Symbol('R_h (x=11,y=-1)')] == -15

    assert result[Symbol('T (x=11,y=-1)')] == -444

    expected_axial_force = -15*SingularityFunction(x, 0, 0) + 87*SingularityFunction(x, 4, 0)/5 + 18*SingularityFunction(x, 4, 1)/5 - 427*SingularityFunction(x, 9, 0)/5 - 42*SingularityFunction(x, 9, 1)/5 + 89*SingularityFunction(x, 14, 0) + 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) - 21*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 111*SingularityFunction(x, 9, 0)/5 + 6*SingularityFunction(x, 9, 1)/5 + 444*SingularityFunction(x, 14, -1) + 48*SingularityFunction(x, 14, 0) + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) - 21*SingularityFunction(x, 4, 1)/5 + 3*SingularityFunction(x, 4, 2)/5 + 111*SingularityFunction(x, 9, 1)/5 + 3*SingularityFunction(x, 9, 2)/5 + 444*SingularityFunction(x, 14, 0) + 48*SingularityFunction(x, 14, 1) + 9*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(8, 3, 15, 0, -1)
    s.apply_load(2, 0, 16, 270, -1)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 8, 3)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -70

    assert result[Symbol('R_h (x=11,y=-1)')] == -15

    assert result[Symbol('T (x=11,y=-1)')] == -450

    expected_axial_force = 24*SingularityFunction(x, 4, 0) + 18*SingularityFunction(x, 4, 1)/5 - 107*SingularityFunction(x, 9, 0) - 18*SingularityFunction(x, 9, 1)/5 + 65*SingularityFunction(x, 14, 0)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) - 16*SingularityFunction(x, 2, 0) + 8*SingularityFunction(x, 4, 0) + 6*SingularityFunction(x, 4, 1)/5 + 26*SingularityFunction(x, 9, 0) + 24*SingularityFunction(x, 9, 1)/5 + 450*SingularityFunction(x, 14, -1) + 30*SingularityFunction(x, 14, 0)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) - 16*SingularityFunction(x, 2, 1) + 8*SingularityFunction(x, 4, 1) + 3*SingularityFunction(x, 4, 2)/5 + 26*SingularityFunction(x, 9, 1) + 12*SingularityFunction(x, 9, 2)/5 + 450*SingularityFunction(x, 14, 0) + 30*SingularityFunction(x, 14, 1)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, 100, 100, 100)
    s.add_member(3, 4, 7, 1, 100, 100, 100)
    s.apply_load(3, 4, 18, 0, 0, 7, 1)
    s.apply_support(0, 0,"fixed")
    s.apply_support(7, 1,"roller")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == Rational(5427,1120)

    assert result[Symbol('R_h (x=0,y=0)')] == -90

    assert result[Symbol('R_v (x=7,y=1)')] == Rational(-5427,1120)

    assert result[Symbol('T (x=0,y=0)')] == Rational(30573,160)

    expected_axial_force = 81027*SingularityFunction(x, 0, 0)/1400 + 8973*SingularityFunction(x, 5, 0)/800 - 72*SingularityFunction(x, 5, 1)/5 + 16281*SingularityFunction(x, 10, 0)/5600 + 72*SingularityFunction(x, 10, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -30573*SingularityFunction(x, 0, -1)/160 + 386919*SingularityFunction(x, 0, 0)/5600 - 711027*SingularityFunction(x, 5, 0)/5600 + 54*SingularityFunction(x, 5, 1)/5 + 5427*SingularityFunction(x, 10, 0)/1400 - 54*SingularityFunction(x, 10, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -30573*SingularityFunction(x, 0, 0)/160 + 386919*SingularityFunction(x, 0, 1)/5600 - 711027*SingularityFunction(x, 5, 1)/5600 + 27*SingularityFunction(x, 5, 2)/5 + 5427*SingularityFunction(x, 10, 1)/1400 - 27*SingularityFunction(x, 10, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())


def test_distributed_loads():
    # pure horizontal UDL load through out the structure
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 0, 0, 3, 4)
    s.apply_load(3, 4, 60, 0, 0, 6, 0)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 200

    assert result[Symbol('R_h (x=0,y=0)')] == -300

    assert result[Symbol('R_v (x=6,y=0)')] == -200

    assert result[Symbol('R_h (x=6,y=0)')] == -300

    expected_axial_force = 340*SingularityFunction(x, 0, 0) - 36*SingularityFunction(x, 0, 1) - 320*SingularityFunction(x, 5, 0) + 340*SingularityFunction(x, 10, 0) + 36*SingularityFunction(x, 10, 1)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 120*SingularityFunction(x, 0, 0) - 48*SingularityFunction(x, 0, 1) + 96*SingularityFunction(x, 5, 1) - 120*SingularityFunction(x, 10, 0) - 48*SingularityFunction(x, 10, 1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 120*SingularityFunction(x, 0, 1) - 24*SingularityFunction(x, 0, 2) + 48*SingularityFunction(x, 5, 2) - 120*SingularityFunction(x, 10, 1) - 24*SingularityFunction(x, 10, 2)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    #pure vertical UDL load through out the structure
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 8, 3)
    s.apply_load(8, 3, 6, 270, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -84

    assert result[Symbol('R_h (x=11,y=-1)')] == 0

    assert result[Symbol('T (x=11,y=-1)')] == -411

    expected_axial_force = 72*SingularityFunction(x, 4, 0)/5 + 18*SingularityFunction(x, 4, 1)/5 - 378*SingularityFunction(x, 9, 0)/5 - 42*SingularityFunction(x, 9, 1)/5 + 336*SingularityFunction(x, 14, 0)/5 + 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) + 24*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 54*SingularityFunction(x, 9, 0)/5 + 6*SingularityFunction(x, 9, 1)/5 + 411*SingularityFunction(x, 14, -1) + 252*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) + 24*SingularityFunction(x, 4, 1)/5 + 3*SingularityFunction(x, 4, 2)/5 + 54*SingularityFunction(x, 9, 1)/5 + 3*SingularityFunction(x, 9, 2)/5 + 411*SingularityFunction(x, 14, 0) + 252*SingularityFunction(x, 14, 1)/5 + 9*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 8, 3)
    s.apply_load(8, 3, 6, 270, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -84

    assert result[Symbol('R_h (x=11,y=-1)')] == 0

    assert result[Symbol('T (x=11,y=-1)')] == -411

    expected_axial_force = 72*SingularityFunction(x, 4, 0)/5 + 18*SingularityFunction(x, 4, 1)/5 - 378*SingularityFunction(x, 9, 0)/5 - 42*SingularityFunction(x, 9, 1)/5 + 336*SingularityFunction(x, 14, 0)/5 + 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) + 24*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 54*SingularityFunction(x, 9, 0)/5 + 6*SingularityFunction(x, 9, 1)/5 + 411*SingularityFunction(x, 14, -1) + 252*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) + 24*SingularityFunction(x, 4, 1)/5 + 3*SingularityFunction(x, 4, 2)/5 + 54*SingularityFunction(x, 9, 1)/5 + 3*SingularityFunction(x, 9, 2)/5 + 411*SingularityFunction(x, 14, 0) + 252*SingularityFunction(x, 14, 1)/5 + 9*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 6, 0, 0, 4, 0)
    s.apply_load(4, 0, 6, 0, 0, 8, 3)
    s.apply_load(8, 3, 6, 0, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == 0

    assert result[Symbol('R_h (x=11,y=-1)')] == -84

    assert result[Symbol('T (x=11,y=-1)')] == 159

    expected_axial_force = -6*SingularityFunction(x, 0, 1) + 24*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 54*SingularityFunction(x, 9, 0)/5 + 6*SingularityFunction(x, 9, 1)/5 + 252*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -72*SingularityFunction(x, 4, 0)/5 - 18*SingularityFunction(x, 4, 1)/5 + 378*SingularityFunction(x, 9, 0)/5 + 42*SingularityFunction(x, 9, 1)/5 - 159*SingularityFunction(x, 14, -1) - 336*SingularityFunction(x, 14, 0)/5 - 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -72*SingularityFunction(x, 4, 1)/5 - 9*SingularityFunction(x, 4, 2)/5 + 378*SingularityFunction(x, 9, 1)/5 + 21*SingularityFunction(x, 9, 2)/5 - 159*SingularityFunction(x, 14, 0) - 336*SingularityFunction(x, 14, 1)/5 - 12*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # UDL horizontal loads on end members
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 6, 0, 0, 4, 0)
    s.apply_load(8, 3, 6, 0, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == 0

    assert result[Symbol('R_h (x=11,y=-1)')] == -54

    assert result[Symbol('T (x=11,y=-1)')] == 84

    expected_axial_force = -6*SingularityFunction(x, 0, 1) + 24*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1) + 24*SingularityFunction(x, 9, 0)/5 - 18*SingularityFunction(x, 9, 1)/5 + 162*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -72*SingularityFunction(x, 4, 0)/5 + 168*SingularityFunction(x, 9, 0)/5 + 24*SingularityFunction(x, 9, 1)/5 - 84*SingularityFunction(x, 14, -1) - 216*SingularityFunction(x, 14, 0)/5 - 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -72*SingularityFunction(x, 4, 1)/5 + 168*SingularityFunction(x, 9, 1)/5 + 12*SingularityFunction(x, 9, 2)/5 - 84*SingularityFunction(x, 14, 0) - 216*SingularityFunction(x, 14, 1)/5 - 12*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # UDL vertical loads on end members
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 6, 270, 0, 4, 0)
    s.apply_load(8, 3, 6, 270, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -54

    assert result[Symbol('R_h (x=11,y=-1)')] == 0

    assert result[Symbol('T (x=11,y=-1)')] == -261

    expected_axial_force = 72*SingularityFunction(x, 4, 0)/5 - 168*SingularityFunction(x, 9, 0)/5 - 24*SingularityFunction(x, 9, 1)/5 + 216*SingularityFunction(x, 14, 0)/5 + 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 0, 1) + 24*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1) + 24*SingularityFunction(x, 9, 0)/5 - 18*SingularityFunction(x, 9, 1)/5 + 261*SingularityFunction(x, 14, -1) + 162*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 0, 2) + 24*SingularityFunction(x, 4, 1)/5 + 3*SingularityFunction(x, 4, 2) + 24*SingularityFunction(x, 9, 1)/5 - 9*SingularityFunction(x, 9, 2)/5 + 261*SingularityFunction(x, 14, 0) + 162*SingularityFunction(x, 14, 1)/5 + 9*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    #UDL loads from middle of the members
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(2, 0, 6, 270, 0, 4, 0)
    s.apply_load(4, 0, 6, 270, 0, 6, 1.5)
    s.apply_load(9.5, 1, 6, 270, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == -42

    assert result[Symbol('R_h (x=11,y=-1)')] == 0

    assert result[Symbol('T (x=11,y=-1)')] == Rational(-19725,100)

    expected_axial_force = 36*SingularityFunction(x, 4, 0)/5 + 18*SingularityFunction(x, 4, 1)/5 - 18*SingularityFunction(x, Rational(13,2), 1)/5 - 189*SingularityFunction(x, 9, 0)/5 - 24*SingularityFunction(x, Rational(23,2), 1)/5 + 168*SingularityFunction(x, 14, 0)/5 + 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -6*SingularityFunction(x, 2, 1) + 12*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 24*SingularityFunction(x, Rational(13,2), 1)/5 + 27*SingularityFunction(x, 9, 0)/5 - 18*SingularityFunction(x, Rational(23,2), 1)/5 + 789*SingularityFunction(x, 14, -1)/4 + 126*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -3*SingularityFunction(x, 2, 2) + 12*SingularityFunction(x, 4, 1)/5 + 3*SingularityFunction(x, 4, 2)/5 + 12*SingularityFunction(x, Rational(13,2), 2)/5 + 27*SingularityFunction(x, 9, 1)/5 - 9*SingularityFunction(x, Rational(23,2), 2)/5 + 789*SingularityFunction(x, 14, 0)/4 + 126*SingularityFunction(x, 14, 1)/5 + 9*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    #UDL loads from middle of the members
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(2, 0, 6, 0, 0, 4, 0)
    s.apply_load(4, 0, 6, 0, 0, 6, 1.5)
    s.apply_load(9.5, 1, 6, 0, 0, 11, -1)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == 0

    assert result[Symbol('R_h (x=11,y=-1)')] == -42

    assert result[Symbol('T (x=11,y=-1)')] == Rational(5325,100)

    expected_axial_force = -6*SingularityFunction(x, 2, 1) + 12*SingularityFunction(x, 4, 0)/5 + 6*SingularityFunction(x, 4, 1)/5 + 24*SingularityFunction(x, Rational(13,2), 1)/5 + 27*SingularityFunction(x, 9, 0)/5 - 18*SingularityFunction(x, Rational(23,2), 1)/5 + 126*SingularityFunction(x, 14, 0)/5 + 18*SingularityFunction(x, 14, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -36*SingularityFunction(x, 4, 0)/5 - 18*SingularityFunction(x, 4, 1)/5 + 18*SingularityFunction(x, Rational(13,2), 1)/5 + 189*SingularityFunction(x, 9, 0)/5 + 24*SingularityFunction(x, Rational(23,2), 1)/5 - 213*SingularityFunction(x, 14, -1)/4 - 168*SingularityFunction(x, 14, 0)/5 - 24*SingularityFunction(x, 14, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -36*SingularityFunction(x, 4, 1)/5 - 9*SingularityFunction(x, 4, 2)/5 + 9*SingularityFunction(x, Rational(13,2), 2)/5 + 189*SingularityFunction(x, 9, 1)/5 + 12*SingularityFunction(x, Rational(23,2), 2)/5 - 213*SingularityFunction(x, 14, 0)/4 - 168*SingularityFunction(x, 14, 1)/5 - 12*SingularityFunction(x, 14, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())


def test_hinge_problems():
    # Alex Example 2(https://oit.tudelft.nl/Macaulays-method/theses/Macaulay-2D/BEP%20v2.html)
    #horizontal UDL with hinge
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 0, 0, 3, 4)
    s.apply_load(3, 4, 60, 0, 0, 6, 0)
    s.apply_rotation_hinge(3, 4)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 200

    assert result[Symbol('R_h (x=0,y=0)')] == -300

    assert result[Symbol('R_v (x=6,y=0)')] == -200

    assert result[Symbol('R_h (x=6,y=0)')] == -300

    expected_axial_force = 340*SingularityFunction(x, 0, 0) - 36*SingularityFunction(x, 0, 1) - 320*SingularityFunction(x, 5, 0) + 340*SingularityFunction(x, 10, 0) + 36*SingularityFunction(x, 10, 1)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 120*SingularityFunction(x, 0, 0) - 48*SingularityFunction(x, 0, 1) + 96*SingularityFunction(x, 5, 1) - 120*SingularityFunction(x, 10, 0) - 48*SingularityFunction(x, 10, 1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 120*SingularityFunction(x, 0, 1) - 24*SingularityFunction(x, 0, 2) + 48*SingularityFunction(x, 5, 2) - 120*SingularityFunction(x, 10, 1) - 24*SingularityFunction(x, 10, 2)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # vertical UDL with hinge
    E, I, A = symbols('E I A')
    s = Structure2d()
    s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 270, 0, 3, 4)
    s.apply_load(3, 4, 60, 270, 0, 6, 0)
    s.apply_rotation_hinge(3, 4)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == -300

    assert result[Symbol('R_h (x=0,y=0)')] == Rational(225,2)

    assert result[Symbol('R_v (x=6,y=0)')] == -300

    assert result[Symbol('R_h (x=6,y=0)')] == Rational(-225,2)

    expected_axial_force = -615*SingularityFunction(x, 0, 0)/2 + 48*SingularityFunction(x, 0, 1) - 96*SingularityFunction(x, 5, 1) + 615*SingularityFunction(x, 10, 0)/2 + 48*SingularityFunction(x, 10, 1)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 90*SingularityFunction(x, 0, 0) - 36*SingularityFunction(x, 0, 1) - 375*SingularityFunction(x, 5, -2)/4 + 180*SingularityFunction(x, 5, 0) + 90*SingularityFunction(x, 10, 0) + 36*SingularityFunction(x, 10, 1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 90*SingularityFunction(x, 0, 1) - 18*SingularityFunction(x, 0, 2) - 375*SingularityFunction(x, 5, -1)/4 + 180*SingularityFunction(x, 5, 1) + 90*SingularityFunction(x, 10, 1) + 18*SingularityFunction(x, 10, 2)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # point load on hinge
    s = Structure2d()
    s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(3, 4, 60, 270, -1)
    s.apply_rotation_hinge(3, 4)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == -30

    assert result[Symbol('R_h (x=0,y=0)')] == Rational(45,2)

    assert result[Symbol('R_v (x=6,y=0)')] == -30

    assert result[Symbol('R_h (x=6,y=0)')] == Rational(-45,2)

    expected_axial_force = -75*SingularityFunction(x, 0, 0)/2 + 75*SingularityFunction(x, 10, 0)/2
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 225*SingularityFunction(x, 5, -2)/4
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 225*SingularityFunction(x, 5, -1)/4
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    # vertical point load on hinge
    E, I, A = symbols('E I A')
    s = Structure2d()
    E = 10**4
    I = 10**4
    A=10**4
    s.add_member(0,0,4,3,E,I,A)
    s.add_member(4,3,7,-1,E,I,A)
    s.apply_support(0,0,"pin")
    s.apply_support(7,-1,"pin")
    s.apply_rotation_hinge(4,3)
    s.apply_load(4,3,375,270,-1)
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == -135

    assert result[Symbol('R_h (x=0,y=0)')] == 180

    assert result[Symbol('R_v (x=7,y=-1)')] == -240

    assert result[Symbol('R_h (x=7,y=-1)')] == -180

    expected_axial_force = -225*SingularityFunction(x, 0, 0) - 75*SingularityFunction(x, 5, 0) + 300*SingularityFunction(x, 10, 0)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 525*SingularityFunction(x, 5, -2)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 525*SingularityFunction(x, 5, -1)
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())


def test_joint_loads():
    # load on a joint
    E = 10**4
    I = 10**4
    A = 10**4
    s1 = Structure2d()
    s1.add_member(0, 0, 3, 4, E, I, A)
    s1.add_member(3, 4, 8, 16, E, I, A)
    s1.add_member(8, 16, 15, 40, E, I, A)
    s1.apply_support(0, 0, "pin")
    s1.apply_support(15, 40, "pin")
    s1.apply_load(3, 4, 150, 270, -1)
    s1.apply_rotation_hinge(3, 4)
    result = s1.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 120

    assert result[Symbol('R_h (x=0,y=0)')] == -90

    assert result[Symbol('R_v (x=15,y=40)')] == -270

    assert result[Symbol('R_h (x=15,y=40)')] == 90

    expected_axial_force = 150*SingularityFunction(x, 0, 0) + 1740*SingularityFunction(x, 5, 0)/13 + 36*SingularityFunction(x, 18, 0)/65 - 1422*SingularityFunction(x, 43, 0)/5
    assert simplify(s1.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 1733053*SingularityFunction(x, 5, -2)/130 - 270*SingularityFunction(x, 5, 0)/13 + 2052*SingularityFunction(x, 18, 0)/65 - 54*SingularityFunction(x, 43, 0)/5
    assert simplify(s1.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 1733053*SingularityFunction(x, 5, -1)/130 - 270*SingularityFunction(x, 5, 1)/13 + 2052*SingularityFunction(x, 18, 1)/65 - 54*SingularityFunction(x, 43, 1)/5
    assert simplify(s1.M.expand()) == simplify(expected_moment.expand())

    # horizontal load on a hinge
    E = 10**4
    I = 10**4
    A = 10**4
    s1 = Structure2d()
    s1.add_member(0, 0, 3, 4, E, I, A)
    s1.add_member(3, 4, 8, 16, E, I, A)
    s1.add_member(8, 16, 15, 40, E, I, A)
    s1.apply_support(0, 0, "pin")
    s1.apply_support(15, 40, "pin")
    s1.apply_load(3, 4, 150, 0, -1)
    s1.apply_rotation_hinge(3, 4)
    result = s1.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 360

    assert result[Symbol('R_h (x=0,y=0)')] == -270

    assert result[Symbol('R_v (x=15,y=40)')] == -360

    assert result[Symbol('R_h (x=15,y=40)')] == 120

    expected_axial_force = 450*SingularityFunction(x, 0, 0) - 930*SingularityFunction(x, 5, 0)/13 + 48*SingularityFunction(x, 18, 0)/65 - 1896*SingularityFunction(x, 43, 0)/5
    assert simplify(s1.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 2421779*SingularityFunction(x, 5, -2)/130 - 360*SingularityFunction(x, 5, 0)/13 + 2736*SingularityFunction(x, 18, 0)/65 - 72*SingularityFunction(x, 43, 0)/5
    assert simplify(s1.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 2421779*SingularityFunction(x, 5, -1)/130 - 360*SingularityFunction(x, 5, 1)/13 + 2736*SingularityFunction(x, 18, 1)/65 - 72*SingularityFunction(x, 43, 1)/5
    assert simplify(s1.M.expand()) == simplify(expected_moment.expand())

    # mix of point loads in different directions
    s2 = Structure2d()
    E = 10**4
    I = 10**4
    A = 10**4
    s2.add_member(0, 0, 5, 12, E, I, A)
    s2.add_member(5, 12, 14, 52, E, I, A)
    s2.apply_support(0, 0, "pin")
    s2.apply_support(14, 52, "pin")
    s2.apply_load(0, 0, 150, 0, -1)
    s2.apply_load(5, 12, 150, 270, -1)
    s2.apply_load(14, 52, 150, 180, -1)
    s2.apply_rotation_hinge(5, 12)
    result = s2.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == Rational(4050,23)

    assert result[Symbol('R_h (x=0,y=0)')] == Rational(-10275,46)

    assert result[Symbol('R_v (x=14,y=52)')] == Rational(-7500,23)

    assert result[Symbol('R_h (x=14,y=52)')] == Rational(10275,46)

    expected_axial_force = 8775*SingularityFunction(x, 0, 0)/46 + 3300*SingularityFunction(x, 13, 0)/23 - 15375*SingularityFunction(x, 54, 0)/46
    assert simplify(s2.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = 10001775*SingularityFunction(x, 13, -2)/1058
    assert simplify(s2.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = 10001775*SingularityFunction(x, 13, -1)/1058
    assert simplify(s2.M.expand()) == simplify(expected_moment.expand())


def test_moment_loads():
    # 3 member structure with mixed loads
    # one fixed support
    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 4, 0, E, I, A)
    s.add_member(4, 0, 8, 3, E, I, A)
    s.add_member(8, 3, 11, -1, E, I, A)
    s.apply_load(0, 0, 15, 0, -2)
    s.apply_load(2, 0, 16, 270, -2)
    s.apply_load(6, 1.5, -6, 270, -2)
    s.apply_load(0, 0, 6, 270, -2)
    s.apply_support(11,-1,"fixed")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=11,y=-1)')] == 0

    assert result[Symbol('R_h (x=11,y=-1)')] == 0

    assert result[Symbol('T (x=11,y=-1)')] == -10

    expected_axial_force = 0
    assert simplify(s.N.expand()) == simplify(expected_axial_force)

    expected_shear_force = -21*SingularityFunction(x, 0, -1) - 16*SingularityFunction(x, 2, -1) + 6*SingularityFunction(x, Rational(13,2), -1) + 10*SingularityFunction(x, 14, -1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -21*SingularityFunction(x, 0, 0) - 16*SingularityFunction(x, 2, 0) + 6*SingularityFunction(x, Rational(13,2), 0) + 10*SingularityFunction(x, 14, 0) + 21
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, I, A)
    s.add_member(3, 4, 6, 0, E, I, A)
    s.apply_load(0, 0, 60, 0, -2)
    s.apply_load(3, 4, 60, 0, 0, 6, 0)
    s.apply_rotation_hinge(3, 4)
    s.apply_support(0, 0,"pin")
    s.apply_support(6, 0,"pin")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == 100

    assert result[Symbol('R_h (x=0,y=0)')] == -75

    assert result[Symbol('R_v (x=6,y=0)')] == -100

    assert result[Symbol('R_h (x=6,y=0)')] == -225

    expected_axial_force = 125*SingularityFunction(x, 0, 0) - 160*SingularityFunction(x, 5, 0) - 36*SingularityFunction(x, 5, 1) + 215*SingularityFunction(x, 10, 0) + 36*SingularityFunction(x, 10, 1)
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -60*SingularityFunction(x, 0, -1) + 250*SingularityFunction(x, 5, -2) - 120*SingularityFunction(x, 5, 0) + 48*SingularityFunction(x, 5, 1) - 120*SingularityFunction(x, 10, 0) - 48*SingularityFunction(x, 10, 1)
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -60*SingularityFunction(x, 0, 0) + 250*SingularityFunction(x, 5, -1) - 120*SingularityFunction(x, 5, 1) + 24*SingularityFunction(x, 5, 2) - 120*SingularityFunction(x, 10, 1) - 24*SingularityFunction(x, 10, 2) + 60
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())

    s = Structure2d()
    E, I, A = symbols('E I A')
    s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
    s.add_member(3, 4, 7, 1, E, I, A)
    s.apply_load(3, 4, 60, 270, -2)
    s.apply_load(3, 4, 18, 0, 0, 7, 1)
    s.apply_support(0, 0,"fixed")
    s.apply_support(7, 1,"roller")
    result = s.solve_for_reaction_loads()

    assert result[Symbol('R_v (x=0,y=0)')] == Rational(-639,160)

    assert result[Symbol('R_h (x=0,y=0)')] == -90

    assert result[Symbol('R_v (x=7,y=1)')] == Rational(639,160)

    assert result[Symbol('T (x=0,y=0)')] == Rational(30873,160)

    expected_axial_force = 10161*SingularityFunction(x, 0, 0)/200 + 18873*SingularityFunction(x, 5, 0)/800 - 72*SingularityFunction(x, 5, 1)/5 - 1917*SingularityFunction(x, 10, 0)/800 + 72*SingularityFunction(x, 10, 1)/5
    assert simplify(s.N.expand()) == simplify(expected_axial_force.expand())

    expected_shear_force = -30873*SingularityFunction(x, 0, -1)/160 + 59517*SingularityFunction(x, 0, 0)/800 - 60*SingularityFunction(x, 5, -1) - 100161*SingularityFunction(x, 5, 0)/800 + 54*SingularityFunction(x, 5, 1)/5 - 639*SingularityFunction(x, 10, 0)/200 - 54*SingularityFunction(x, 10, 1)/5
    assert simplify(s.V.expand()) == simplify(expected_shear_force.expand())

    expected_moment = -30873*SingularityFunction(x, 0, 0)/160 + 59517*SingularityFunction(x, 0, 1)/800 - 60*SingularityFunction(x, 5, 0) - 100161*SingularityFunction(x, 5, 1)/800 + 27*SingularityFunction(x, 5, 2)/5 - 639*SingularityFunction(x, 10, 1)/200 - 27*SingularityFunction(x, 10, 2)/5
    assert simplify(s.M.expand()) == simplify(expected_moment.expand())
