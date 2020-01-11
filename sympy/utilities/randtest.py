from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.randtest",
    useinstead="Import from sympy.testing.randtest",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.randtest import *  # noqa:F401
