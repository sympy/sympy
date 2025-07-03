from sympy.physics.continuum_mechanics.structure2d import Structure2d
from sympy.core.symbol import Symbol, symbols
from sympy.functions import SingularityFunction
from sympy.simplify.simplify import simplify
# from sympy import pi, cos, sin, rad, sympify
import pytest
import math

x = Symbol('x')
F = Symbol('F')
# Rv_0, Rh_0, Rv_4 = symbols('R_v__0,__0 R_h__0,__0 R_v__4,__0')
Rv_0 = Symbol('R_v__0,__0')
Rh_0 = Symbol('R_h__0,__0')
Rv_4 = Symbol('R_v__4,__0')



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
    Rv_1, Rh_1 = s.apply_support(x=0, y=0, type='pin')
    Rv_2 = s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2, 0, value, global_angle=angle, order=-1)
    s.solve_for_reaction_loads(Rv_1, Rh_1, Rv_2)

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
    Rv_1, Rh_1 = s.apply_support(x=0, y=0, type='pin')
    Rv_2 = s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2, 0, value, global_angle=angle, order=-1)
    s.solve_for_reaction_loads(Rv_1, Rh_1, Rv_2)

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
    Rv_1, Rh_1 = s.apply_support(x=0, y=0, type='pin')
    Rv_2 = s.apply_support(x=4, y=0, type='roller')

    # Check support properties
    supports = [(0, 0, 'pin'), (4, 0, 'roller')]
    for i, (x, y, node_type) in enumerate(supports):
        support = s.supports[i]
        assert support.x == x
        assert support.y == y
        assert support.node_type == node_type

    # Apply load and solve for reactions
    s.apply_load(2,0,value,global_angle=angle,order=0,end_x=3,end_y=0)
    s.solve_for_reaction_loads(Rv_1, Rh_1, Rv_2)

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
Rv_7 = Symbol('R_v__7,__-1')
Rh_7 = Symbol('R_h__7,__-1')
T_7 = Symbol('T__7,__-1')

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
    Rh1, Rv1, T1 = s.apply_support(x=7, y=-1, type='fixed')

    s.apply_load(5,1.5,value/2,global_angle=s.members[1].angle_deg + 270,order=0,end_x=7,end_y=-1)
    s.apply_load(0,0,value*0.8,global_angle=270,order=0,end_x=3,end_y=4)
    s.solve_for_reaction_loads(Rh1,Rv1,T1)
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
def test_structure2d_multiple_horizontal_reactions():
    """Test that the structure can handle multiple horizontal reactions."""
    E, I, A, Q = symbols('E I A Q')
    s = Structure2d()
    s.add_member(0, 0, 6, 0, E, I, A)

    # Apply supports
    T_1,Rv_1, Rh_1 = s.apply_support(x=0, y=0, type='fixed')
    Rv_2,Rh_2 = s.apply_support(x=6, y=0, type='pin')

    # Apply horizontal loads
    s.apply_load(0, 0, -Q, global_angle=0, order=0, end_x=6, end_y=0)
    s.apply_load(0, 0, -Q, global_angle=90, order=0, end_x=6, end_y=0)

    # Solve for reaction loads
    s.solve_for_reaction_loads(T_1,Rv_1, Rh_1, Rv_2, Rh_2)
    p= s.reaction_loads
    q = {
        Symbol('R_v__0,__0'): -15*Q/4,
        Symbol('R_v__6,__0'):  -9*Q/4,
        Symbol('T__0,__0'):    9*Q/2,
        Symbol('R_h__0,__0'):  3*Q,
        Symbol('R_h__6,__0'):  3*Q,
    }
    assert p==q
