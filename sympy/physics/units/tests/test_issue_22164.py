from sympy.physics.units.systems.si import SI
from sympy.physics.units.definitions.dimension_definitions import length
from sympy.physics.units.quantities import Quantity
from warnings import simplefilter

def test_deprecated_warning():
    simplefilter("error")
    dm = Quantity("dm")
    SI.set_quantity_dimension(dm, length)
    SI.set_quantity_scale_factor(dm, 1)

    bad_exp = Quantity("bad_exp")
    SI.set_quantity_dimension(bad_exp, length)
    SI.set_quantity_scale_factor(bad_exp, 1)

    expr = dm ** bad_exp

    quantity = SI._collect_factor_and_dimension(expr)
