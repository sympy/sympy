from sympy.physics.continuum_mechanics.structure2d import Structure2d
from sympy.core.symbol import Symbol, symbols
from sympy.functions import SingularityFunction
from sympy.simplify.simplify import simplify
from sympy import pi, cos, sin, rad
import pytest

x = Symbol('x')
F = Symbol('F')
Rv_0, Rh_0, Rv_4 = symbols('Rv_0 Rh_0 Rv_4')

# Test Data Definitions
symbolic_test_data = [
    {
        'id': 'symbolic_test_1',
        'angle': 270,
        'value': F,
        'reaction_loads': {Rv_0: -F / 2, Rv_4: -F / 2, Rh_0: 0},
        'shear_force': F * SingularityFunction(x, 0, 0) / 2 - F * SingularityFunction(x, 2, 0) + F * SingularityFunction(x, 4, 0) / 2,
        'bending_moment': {'symbolic': F * SingularityFunction(x, 0, 1) / 2 - F * SingularityFunction(x, 2, 1) + F * SingularityFunction(x, 4, 1) / 2}
    }
]

numerical_test_data = [
    {
        'id': 'numerical_test_1',
        'angle': 270,
        'value': 10,
        'reaction_loads': {Rh_0: 0.0, Rv_0: -5, Rv_4: -5},
        'shear_force': 5 * SingularityFunction(x, 0, 0) - 10 * SingularityFunction(x, 2, 0) + 5 * SingularityFunction(x, 4, 0),
        'bending_moment': {0: 0, 2: 10, 4: 0}
    }
]

@pytest.mark.parametrize('test_data', symbolic_test_data + numerical_test_data, ids=[data['id'] for data in (symbolic_test_data + numerical_test_data)])
def test_structure2d_basic(test_data):
    # Extract test data
    angle = test_data['angle']
    value = test_data['value']
    reaction_loads = test_data['reaction_loads']
    shear_force = test_data['shear_force']
    bending_moment = test_data['bending_moment']

    # Create structure
    E, I, A = symbols('E I A')
    s = Structure2d()
    s.add_member(0, 0, 2, 0, E, I, A)
    s.add_member(2, 0, 4, 0, E, I, A)

    # Assert member properties
    for i, (x1, y1, x2, y2) in enumerate([(0, 0, 2, 0), (2, 0, 4, 0)]):
        member = s.members[i]
        assert (member.x1, member.y1, member.x2, member.y2) == (x1, y1, x2, y2)
        assert member.E == E
        assert member.I == I
        assert member.A == A
        assert member.length == 2
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
        assert s.reaction_loads[load_symbol] - expected_value == 0

    # Check bending moment
    for loc, expected_bm in bending_moment.items():
        if loc == 'symbolic':
            assert s.bending_moment() == expected_bm
        else:
            computed_bm = s.bending_moment(loc)
            assert simplify(computed_bm) == simplify(expected_bm), f"Bending moment mismatch at x={loc}"
