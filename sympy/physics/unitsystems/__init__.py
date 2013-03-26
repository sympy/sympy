# List of ideas: https://github.com/sympy/sympy/wiki/Unit-systems

from sympy.physics.unitsystems.units import (Quantity as Q, unit_simplify,
                                             set_system, get_system)
# FIX: gather all systems in a submodule
from sympy.physics.unitsystems.mks import mks
from sympy.physics.unitsystems.mksa import mksa
from sympy.physics.unitsystems.astro import astro
