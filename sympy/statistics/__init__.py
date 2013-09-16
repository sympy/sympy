"""
SymPy statistics module
Deprecated
See sympy.stats
"""
from sympy.utilities.exceptions import SymPyDeprecationWarning
SymPyDeprecationWarning(
    feature="sympy.statistics",
    useinstead="sympy.stats",
    issue=3386,
    deprecated_since_version="0.7.2",
).warn()
from .distributions import Normal, Uniform
del SymPyDeprecationWarning
