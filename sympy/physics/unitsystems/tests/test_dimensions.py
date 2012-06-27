from sympy import sympify
from sympy.physics.unitsystems.dimensions import Dimension
from sympy.physics.unitsystems.mks import length, mass, time, velocity, energy
from sympy.utilities.pytest import raises

def test_dimension_definition():
    raises(TypeError, lambda: Dimension(["length", 1, 2]))

    assert dict(length) == {sympify('length'): 1}
    assert length.get('length') == 1
    assert length.get('time') is None
    assert length.get('time', 'def') == 'def'
    #assert energy.as_dict == dict(zip(energy.names, energy.powers))

    assert length.is_dimensionless is False

def test_dimension_symbols():
    assert length.symbol == 'L'
    #assert str(length) == 'L'

def test_dimension_operations():
    assert length == length + length
    assert length / length == 1
    assert -length == length
    assert length**2 == Dimension(length=2)
    assert velocity == length / time
    assert velocity == length * (1 / time)
    assert energy == mass * length**2 / time**2
    assert 1 / length == Dimension(length=-1)

    raises(TypeError, lambda: length + 1)
    raises(TypeError, lambda: length + time)
    raises(TypeError, lambda: length - 1)
    raises(TypeError, lambda: length - time)
