from sympy import Symbol
from sympy.codegen.ffunctions import isign, dsign, cmplx, kind


def test_isign():
    x = Symbol('x', integer=True)
    assert isign(1, x) == isign(1, x)


def test_dsign():
    x = Symbol('x')
    assert dsign(1, x) == dsign(1, x)


def test_cmplx():
    x = Symbol('x')
    assert cmplx(1, x) == cmplx(1, x)


def test_kind():
    x = Symbol('x')
    assert kind(x) == kind(x)
