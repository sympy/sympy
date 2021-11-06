from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.randtest",
    useinstead="Import from sympy.core.random",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.core.random import *  # noqa:F401
