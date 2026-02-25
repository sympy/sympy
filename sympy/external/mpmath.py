#
# This module is to centralize sympy's direct access of mpmath functions
# so that there are not many places in the codebase that reach into what are
# arguably mpmath's internal interfaces.
#
# More functions could be added here if needed but these are the ones that were
# used in the sympy codebase when this module was added.
#
from __future__ import annotations
from functools import update_wrapper as _update_wrapper


from mpmath import (
    MPContext,
    MPIntervalContext,
    bernfrac,
    diff,
    eulernum,
    fac,
    findroot,
    inf,
    make_mpc,
    make_mpf,
    mp,
    mpc,
    mpf,
    mpi,
    ninf,
    sqrt,
    workprec,
)

NoConvergence = mp.NoConvergence

from mpmath.ctx_mp_python import mpnumeric
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
from mpmath.matrices.matrices import _matrix

__all__ = [
    "MPContext",
    "MPIntervalContext",
    "MPZ",
    "MPZ_ONE",
    "ComplexResult",
    "NoConvergence",
    "bernfrac",
    "catalan_fixed",
    "diff",
    "dps_to_prec",
    "euler_fixed",
    "eulernum",
    "fac",
    "findroot",
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
    "inf",
    "int_types",
    "isqrt",
    "make_mpf",
    "make_mpc",
    "mp",
    "mpc",
    "mpf",
    "mpi",
    "mpnumeric",
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
    "ninf",
    "normalize",
    "phi_fixed",
    "prec_to_dps",
    "repr_dps",
    "round_nearest",
    "sqrtrem",
    "sqrt",
    "to_float",
    "to_int",
    "to_rational",
    "to_str",
    "workprec",
    "_matrix",
]


def conserve_mpmath_dps(func):
    """Decorator to rest mpmath's global precision after a function call.

    It is not recommended to use this in new code which should instead use
    the :class:`local_workdps` or :class:`local_workprec` context managers.
    """
    import mpmath

    def func_wrapper(*args, **kwargs):
        dps = mpmath.mp.dps
        try:
            return func(*args, **kwargs)
        finally:
            mpmath.mp.dps = dps

    func_wrapper = _update_wrapper(func_wrapper, func)
    return func_wrapper


def _new_mpcontext() -> MPContext:
    # Note sure how much is really needed here...
    ctx = MPContext()
    iv = MPIntervalContext()
    ctx._iv = iv
    ctx.mpi = iv._mpi
    return ctx


class local_workprec:
    """Context manager to borrow an mpmath MPContext with given precision.

    >>> from sympy.external.mpmath import local_workprec
    >>> with local_workprec(10) as ctx:
    ...     print(ctx.sin(1))
    0.84

    The context must not be used outside the ``with`` block, as it will be
    returned to the context pool for reuse. Even just converting an mpmath
    object to a string or sympifying it will use its context so these
    operations should also be done within the ``with`` block.

    Unlike mpmath's ``workdps``, this context manager does not change the
    global precision, it only provides a local context with the requested
    precision for the duration of the ``with`` block. Context methods like
    ``ctx.sin`` must be used instead of global functions like ``mpmath.sin``.
    Likewise ``ctx.fadd(a, b)`` should be used rather than ``a + b``.

    See Also
    ========

    local_workdps
    """

    _contexts: list[MPContext] = []
    _ctx: MPContext

    def __init__(self, prec):
        self.prec = prec
        try:
            self._ctx = self._contexts.pop()
        except IndexError:
            self._ctx = _new_mpcontext()
        self._ctx.prec = prec

    def __enter__(self) -> MPContext:
        """Return the MPContext with the requested precision."""
        return self._ctx

    def __exit__(self, exc_type, exc_value, traceback):
        """Return the MPContext to the pool."""
        self._contexts.append(self._ctx)


class local_workdps:
    """Context manager to borrow an mpmath MPContext with given dps.

    See :class:`local_workprec` for usage.

    See Also
    ========

    local_workprec
    """

    _contexts: list[MPContext] = []
    _ctx: MPContext

    def __init__(self, dps):
        self.dps = dps
        try:
            self._ctx = self._contexts.pop()
        except IndexError:
            self._ctx = _new_mpcontext()
        self._ctx.dps = dps

    def __enter__(self) -> MPContext:
        """Return the MPContext with the requested dps."""
        return self._ctx

    def __exit__(self, exc_type, exc_value, traceback):
        """Return the MPContext to the pool."""
        self._contexts.append(self._ctx)
