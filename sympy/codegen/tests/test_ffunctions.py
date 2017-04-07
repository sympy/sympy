from sympy import Symbol
from sympy.codegen.ffunctions import isign, dsign, cmplx, kind, literal_dp
from sympy.printing.fcode import fcode


def test_isign():
    x = Symbol('x', integer=True)
    assert isign(1, x) == isign(1, x)
    assert fcode(isign(1, x), standard=95, source_format='free') == 'isign(1, x)'


def test_dsign():
    x = Symbol('x')
    assert dsign(1, x) == dsign(1, x)
    assert fcode(dsign(literal_dp(1), x), standard=95, source_format='free') == 'dsign(1d0, x)'


def test_cmplx():
    x = Symbol('x')
    assert cmplx(1, x) == cmplx(1, x)


def test_kind():
    x = Symbol('x')
    assert kind(x) == kind(x)


def test_literal_dp():
    assert fcode(literal_dp(0), source_format='free') == '0d0'
