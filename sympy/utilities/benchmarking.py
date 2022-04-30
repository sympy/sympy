"""
.. deprecated:: 1.6

   sympy.utilities.benchmarking has been renamed to sympy.testing.benchmarking.
"""
from sympy.utilities.exceptions import sympy_deprecation_warning

sympy_deprecation_warning("The sympy.utilities.benchmarking submodule is deprecated. Use sympy.testing.benchmarking instead.",
    deprecated_since_version="1.6",
    active_deprecations_target="deprecated-sympy-utilities-submodules")

from sympy.testing.benchmarking import *  # noqa:F401
