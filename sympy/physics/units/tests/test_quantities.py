from sympy.physics.units import meter

def test_quantity_subs_behavior():
    assert meter.subs(meter, meter) == meter
    assert meter.subs(meter, 5) == 5