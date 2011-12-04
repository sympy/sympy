"""Decorator to conserve mpmath precision."""

import functools
from sympy import mpmath

def CONSERVE_MPMATH_DPS(func):
    """After the function finishes, resets the value of mpmath.mp.dps to the value it had before the function was run."""
    def func_wrapper():
        dps = mpmath.mp.dps
        try:
            func()
        finally:
            mpmath.mp.dps = dps

    func_wrapper = functools.update_wrapper(func_wrapper, func)
    return func_wrapper
