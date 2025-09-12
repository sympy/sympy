from hypothesis import given
from hypothesis import strategies as st
from hypothesis.strategies import composite

from sympy import ring, QQ
from sympy.abc import x
from sympy.external.gmpy import GROUND_TYPES
from sympy.functions.elementary.trigonometric import atan, asin, sin, cos, tan
from sympy.functions.elementary.hyperbolic import atanh, asinh, sinh, cosh, tanh
from sympy.functions.elementary.exponential import exp, log
from sympy.polys.series import power_series_ring
from sympy.polys.polyutils import expr_from_dict
from sympy.polys.densebasic import dup_to_dict
from sympy.polys.ring_series import (
    rs_log,
    rs_exp,
    rs_atan,
    rs_atanh,
    rs_asin,
    rs_asinh,
    rs_tan,
    rs_tanh,
    rs_sin,
    rs_sinh,
    rs_cos,
    rs_cosh,
)

from sympy.testing.pytest import slow

flint = False
if GROUND_TYPES == "flint":
    flint = True


def _dup_QQ():
    """This is a strategy of creating random dup_QQ elements."""
    elems = st.tuples(
        st.integers(min_value=-100, max_value=100),
        st.integers(min_value=1, max_value=100),
    ).map(lambda t: QQ(t[0], t[1]))

    return elems


@composite
def dup_zero_const(draw, min_size=3, max_size=25):
    """DUP zero constant strategy."""
    element_strategy = _dup_QQ()

    lst = draw(st.lists(element_strategy, min_size=min_size, max_size=max_size))
    lst[-1] = QQ.zero
    return lst


@composite
def dup_one_const(draw, min_size=3, max_size=25):
    """DUP one constant strategy."""
    element_strategy = _dup_QQ()

    lst = draw(st.lists(element_strategy, min_size=min_size, max_size=max_size))
    lst[-1] = QQ.one
    return lst


@given(f=dup_zero_const())
def test_rs_series_zero(f):
    Rs, _ = power_series_ring("x", QQ, 25)
    Rp, x = ring("x", QQ)

    s = Rs.from_list(f[::-1])
    p = Rp.from_list(f)

    assert Rs.to_dense((Rs.exp(s))) == (rs_exp(p, x, 25)).to_dense()

    assert Rs.to_dense((Rs.atan(s))) == (rs_atan(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.atanh(s))) == (rs_atanh(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.asin(s))) == (rs_asin(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.asinh(s))) == (rs_asinh(p, x, 25)).to_dense()

    assert Rs.to_dense((Rs.tan(s))) == (rs_tan(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.tanh(s))) == (rs_tanh(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.sin(s))) == (rs_sin(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.sinh(s))) == (rs_sinh(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.cos(s))) == (rs_cos(p, x, 25)).to_dense()
    assert Rs.to_dense((Rs.cosh(s))) == (rs_cosh(p, x, 25)).to_dense()


@given(f=dup_one_const())
def test_rs_series_one(f):
    Rs, _ = power_series_ring("x", QQ, 25)
    Rp, x = ring("x", QQ)

    s = Rs.from_list(f[::-1])
    p = Rp.from_list(f)

    assert Rs.to_dense((Rs.log(s))) == (rs_log(p, x, 25)).to_dense()


@slow
@given(f=dup_zero_const(max_size=5))
def test_global_series_zero(f):
    Rs, _ = power_series_ring("x", QQ, 5)

    s = Rs.from_list(f[::-1])
    e = expr_from_dict(dup_to_dict(f, QQ), x)
    assert Rs.to_expr((Rs.exp(s))) == exp(e).series(x, 0, 5)

    assert Rs.to_expr((Rs.atan(s))) == atan(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.atanh(s))) == atanh(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.asin(s))) == asin(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.asinh(s))) == asinh(e).series(x, 0, 5)

    assert Rs.to_expr((Rs.tan(s))) == tan(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.tanh(s))) == tanh(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.sin(s))) == sin(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.sinh(s))) == sinh(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.cos(s))) == cos(e).series(x, 0, 5)
    assert Rs.to_expr((Rs.cosh(s))) == cosh(e).series(x, 0, 5)


@given(f=dup_one_const(max_size=5))
def test_global_series_one(f):
    Rs, _ = power_series_ring("x", QQ, 5)

    s = Rs.from_list(f[::-1])
    e = expr_from_dict(dup_to_dict(f, QQ), x)

    assert Rs.to_expr((Rs.log(s))) == log(e).series(x, 0, 5)
