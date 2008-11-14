from sympy import Rational
from sympy.physics.units import m, s, day, km, foot, meter, au, \
        speed_of_light, minute

def test_units():
    assert (5*m/s * day) / km == 432
    assert foot / meter == Rational('0.3048')

    # Light from the sun needs about 8.3 minutes to reach earth
    t = (1*au / speed_of_light).evalf() / minute
    assert abs(t - 8.31) < 0.1

    assert (m**2)**Rational(1,2) == m
    assert (m**Rational(1,2))**2 == m
