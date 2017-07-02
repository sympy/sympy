"""
Deprecated: use ``sympy.physics.units``
"""

# DEPRECATED: use `units`
from sympy.utilities.exceptions import SymPyDeprecationWarning

exec("from sympy.physics.units import *")

SymPyDeprecationWarning(
    feature ="sympy.physics.unitsystems",
    useinstead ="sympy.physics.units",
    deprecated_since_version ="1.1",
    issue=12856,
).warn()
