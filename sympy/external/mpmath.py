#
# This module is to centralize sympy's direct access of mpmath functions
# so that there are not many places in the codebase that reach into what are
# arguably mpmath's internal interfaces. The functions referenced here are all
# quite small and so could just be copied from mpmath into sympy if it is a
# problem to depend on being able to import them from mpmath.
#
# More functions could be added here if needed but these are the ones that were
# widely used throughout the sympy codebase beforew this module was added.
#
from mpmath.libmp.libmpf import (
    prec_to_dps,
    dps_to_prec,
    to_str,
    repr_dps,
    from_float,
    to_rational,
)

__all__ = [
    "prec_to_dps",
    "dps_to_prec",
    "to_str",
    "repr_dps",
    "from_float",
    "to_rational",
]
