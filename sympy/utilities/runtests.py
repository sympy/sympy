from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.runtests",
    useinstead="Import from sympy.testing.runtests",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.runtests import *  # noqa:F401
