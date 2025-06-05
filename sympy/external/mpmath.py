#
# This module is to centralize sympy's direct access of mpmath functions
# so that there are not many places in the codebase that reach into what are
# arguably mpmath's internal interfaces. The functions referenced here are all
# quite small and so could just be copied from mpmath into sympy if it is a
# problem to depend on being able to import them from mpmath.
#
# More functions could be added here if needed but these are the ones that were
# widely used throughout the sympy codebase before this module was added.
#
from mpmath.libmp import (
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
    "_giant_steps_mpmath",
]


def _giant_steps_mpmath(start, target, n=2):
    """
    Return a list of integers ~=

    [start, n*start, ..., target/n^2, target/n, target]

    but conservatively rounded so that the quotient between two
    successive elements is actually slightly less than n.

    With n = 2, this describes suitable precision steps for a
    quadratically convergent algorithm such as Newton's method;
    with n = 3 steps for cubic convergence (Halley's method), etc.

        >>> from sympy.polys.ring_series import _giant_steps_mpmath
        >>> _giant_steps_mpmath(50,1000)
        [66, 128, 253, 502, 1000]
        >>> _giant_steps_mpmath(50,1000,4)
        [65, 252, 1000]

    """
    # This function is copied from mpmath to avoid depending on mpmath
    # internals. This implementation could possibly be improved for the
    # particular purpose here though.
    L = [target]
    while L[-1] > start*n:
        L = L + [L[-1]//n + 2]
    return L[::-1]
