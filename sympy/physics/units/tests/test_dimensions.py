from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.mks import length, mass, time, velocity, energy

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
