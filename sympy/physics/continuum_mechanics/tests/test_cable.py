from sympy.physics.continuum_mechanics.cable import Cable

def test_cable():
    c = Cable()
    assert c.supports == {}
    assert c.left_support == []
    assert c.right_support == []
    assert c.loads == {}
    assert c.loads_position == {}
    assert c.length == 0
    assert c.reaction_loads == {}

    # tests for add_supports method
    c.add_supports(('A', 0, 5), ('B', 10, 6))
    assert c.supports == {'A': [0, 5], 'B': [10, 6]}
    assert c.left_support == [0, 5]
    assert c.right_support == [10, 6]
    assert c.length == 0
    assert c.reaction_loads == {}

    # tests for remove_supports method
    c.remove_supports()
    assert c.supports == {}
    assert c.left_support == []
    assert c.right_support == []

    c.add_supports(('A', 0, 5), ('B', 10, 6))

    # tests for change_support method
    c.change_support('B', ('C', 10, 3))
    assert c.supports == {'A': [0, 5], 'C': [10, 3]}
    assert c.left_support == [0, 5]
    assert c.right_support == [10, 3]

    # tests for add_loads method
    c.add_loads(('X', 2, 5, 19, 94), ('Y', 5, 7, 23, 30), ('Z', 8, 3, 40, 45))
    assert c.loads == {'X': [19, 94], 'Y': [23, 30], 'Z': [40, 45]}
    assert c.loads_position == {'X': [2, 5], 'Y': [5, 7], 'Z': [8, 3]}
    assert c.length == 0
    assert c.reaction_loads == {}

    # tests for remove_loads method
    c.remove_loads('Y', 'Z')
    assert c.loads == {'X': [19, 94]}
    assert c.loads_position == {'X': [2, 5]}
    assert c.length == 0
    assert c.reaction_loads == {}
