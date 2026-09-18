from __future__ import annotations
from sympy.external.mpmath import (PythonMPContext, _constant, _mpc, _mpf,
    mpmath, mpnumeric, repr_dps)


def test_repr_dps_is_stable():
    assert repr_dps(13) == 6
    assert repr_dps(33) == 12
    assert repr_dps(53) == 17
    assert repr_dps(66) == 22


def test_mpmath_module_reexported():
    assert mpmath is __import__('mpmath')


def test_ctx_mp_python_names_reexported():
    assert PythonMPContext is __import__('mpmath').ctx_mp_python.PythonMPContext
    assert _mpf is __import__('mpmath').ctx_mp_python._mpf
    assert _mpc is __import__('mpmath').ctx_mp_python._mpc
    assert _constant is __import__('mpmath').ctx_mp_python._constant
    assert mpnumeric is __import__('mpmath').ctx_mp_python.mpnumeric
