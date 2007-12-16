"""Utilities for testing"""
from sympy import Basic

def REPR0(func):
    """decorator: run func under Basic.set_repr_level(0)"""
    def func_wrapper():
        lev = Basic.set_repr_level(0)
        try:
            func()
        finally:
            Basic.set_repr_level(lev)

    return func_wrapper

