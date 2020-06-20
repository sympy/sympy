from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.benchmarking",
    useinstead="Import from sympy.testing.benchmarking",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.benchmarking import *  # noqa:F401
