from sympy import *
from sympy.physics.units import *

def test_units():
    assert (5*m/s * day) / km == 432
    assert foot / meter == Rational('0.3048')

    # Light from the sun needs about 8.3 minutes to reach earth
    t = (1*au / speed_of_light).evalf() / minute
    assert abs(t - 8.31) < 0.1
