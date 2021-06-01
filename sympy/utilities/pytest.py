from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.pytest",
    useinstead="Import from sympy.testing.pytest",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.pytest import *  # noqa:F401
