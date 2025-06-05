#
# This module is to centralize sympy's direct access of mpmath functions
# so that there are not many places in the codebase that reach into what are
# arguably mpmath's internal interfaces.
#
# More functions could be added here if needed but these are the ones that were
# used in the sympy codebase when this module was added.
#
from mpmath.libmp import (
    MPZ,
    MPZ_ONE,
    ComplexResult,
    catalan_fixed,
    dps_to_prec,
    euler_fixed,
    fhalf,
    finf,
    fnan,
    fninf,
    fnone,
    fone,
    from_float,
    from_int,
    from_man_exp,
    from_rational,
    from_str,
    fzero,
    ifac,
    ifib,
    int_types,
    isqrt,
    mpc_abs,
    mpc_exp,
    mpc_pow,
    mpc_pow_int,
    mpc_pow_mpf,
    mpc_sqrt,
    mpf_abs,
    mpf_add,
    mpf_atan,
    mpf_atan2,
    mpf_ceil,
    mpf_cmp,
    mpf_cos,
    mpf_cosh_sinh,
    mpf_div,
    mpf_e,
    mpf_exp,
    mpf_floor,
    mpf_ge,
    mpf_gt,
    mpf_le,
    mpf_log,
    mpf_lt,
    mpf_mod,
    mpf_mul,
    mpf_neg,
    mpf_pi,
    mpf_pow,
    mpf_pow_int,
    mpf_shift,
    mpf_sin,
    mpf_sqrt,
    mpf_sub,
    mpf_tan,
    normalize,
    phi_fixed,
    prec_to_dps,
    repr_dps,
    round_nearest,
    sqrtrem,
    to_float,
    to_int,
    to_rational,
    to_str,
)
from mpmath.libmp.libintmath import giant_steps

__all__ = [
    "MPZ",
    "MPZ_ONE",
    "ComplexResult",
    "catalan_fixed",
    "dps_to_prec",
    "euler_fixed",
    "fhalf",
    "finf",
    "fnan",
    "fninf",
    "fnone",
    "fone",
    "from_float",
    "from_int",
    "from_man_exp",
    "from_rational",
    "from_str",
    "fzero",
    "giant_steps",
    "ifac",
    "ifib",
    "int_types",
    "isqrt",
    "mpc_abs",
    "mpc_exp",
    "mpc_pow",
    "mpc_pow_int",
    "mpc_pow_mpf",
    "mpc_sqrt",
    "mpf_abs",
    "mpf_add",
    "mpf_atan",
    "mpf_atan2",
    "mpf_ceil",
    "mpf_cmp",
    "mpf_cos",
    "mpf_cosh_sinh",
    "mpf_div",
    "mpf_e",
    "mpf_exp",
    "mpf_floor",
    "mpf_ge",
    "mpf_gt",
    "mpf_le",
    "mpf_log",
    "mpf_lt",
    "mpf_mod",
    "mpf_mul",
    "mpf_neg",
    "mpf_pi",
    "mpf_pow",
    "mpf_pow_int",
    "mpf_shift",
    "mpf_sin",
    "mpf_sqrt",
    "mpf_sub",
    "mpf_tan",
    "normalize",
    "phi_fixed",
    "prec_to_dps",
    "repr_dps",
    "round_nearest",
    "sqrtrem",
    "to_float",
    "to_int",
    "to_rational",
    "to_str",
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
