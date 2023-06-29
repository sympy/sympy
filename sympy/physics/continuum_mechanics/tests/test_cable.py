from sympy.physics.continuum_mechanics.cable import Cable
from sympy.core.symbol import Symbol

def test_cable():
    c = Cable(('A', 0, 10), ('B', 10, 10))
    assert c.supports == {'A': [0, 10], 'B': [10, 10]}
    assert c.left_support == [0, 10]
    assert c.right_support == [10, 10]
    assert c.loads == {'distributed': {}, 'point_load': {}}
    assert c.loads_position == {}
    assert c.length == 0
    assert c.reaction_loads == {Symbol("R_A_x"): 0, Symbol("R_A_y"): 0, Symbol("R_B_x"): 0, Symbol("R_B_y"): 0}

    # tests for change_support method
    c.change_support('A', ('C', 12, 3))
    assert c.supports == {'B': [10, 10], 'C': [12, 3]}
    assert c.left_support == [10, 10]
    assert c.right_support == [12, 3]
    assert c.reaction_loads == {Symbol("R_B_x"): 0, Symbol("R_B_y"): 0, Symbol("R_C_x"): 0, Symbol("R_C_y"): 0}

    c.change_support('C', ('A', 0, 10))

    # tests for apply_load method for point loads
    c.apply_load(-1, ('X', 2, 5, 3, 30))
    c.apply_load(-1, ('Y', 5, 8, 5, 60))
    assert c.loads == {'distributed': {}, 'point_load': {'X': [3, 30], 'Y': [5, 60]}}
    assert c.loads_position == {'X': [2, 5], 'Y': [5, 8]}
    assert c.length == 0
    assert c.reaction_loads == {Symbol("R_A_x"): 0, Symbol("R_A_y"): 0, Symbol("R_B_x"): 0, Symbol("R_B_y"): 0}

    # tests for remove_loads method
    c.remove_loads('X')
    assert c.loads == {'distributed': {}, 'point_load': {'Y': [5, 60]}}
    assert c.loads_position == {'Y': [5, 8]}
    assert c.length == 0
    assert c.reaction_loads == {Symbol("R_A_x"): 0, Symbol("R_A_y"): 0, Symbol("R_B_x"): 0, Symbol("R_B_y"): 0}

    c.remove_loads('Y')

    #tests for apply_load method for distributed load
    c.apply_load(0, ('Z', 9))
    assert c.loads == {'distributed': {'Z': 9}, 'point_load': {}}
    assert c.loads_position == {}
    assert c.length == 0
    assert c.reaction_loads == {Symbol("R_A_x"): 0, Symbol("R_A_y"): 0, Symbol("R_B_x"): 0, Symbol("R_B_y"): 0}

    # tests for apply_length method
    c.apply_length(20)
    assert c.length == 20
