from hypothesis import given
from hypothesis import strategies as st
from hypothesis.strategies import composite

from sympy import ring, QQ
from sympy.polys.series import power_series_ring
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
)


@composite
def dup_zero_const(draw, min_size=3, max_size=25):
    """DUP zero constant strategy."""
    element_strategy = st.tuples(st.integers(), st.integers(min_value=1)).map(
        lambda t: QQ(t[0], t[1])
    )

    lst = draw(st.lists(element_strategy, min_size=min_size, max_size=max_size))
    lst[-1] = QQ.zero
    return lst


@composite
def dup_one_const(draw, min_size=3, max_size=25):
    """DUP zero constant strategy."""
    element_strategy = st.tuples(st.integers(), st.integers(min_value=1)).map(
        lambda t: QQ(t[0], t[1])
    )

    lst = draw(st.lists(element_strategy, min_size=min_size, max_size=max_size))
    lst[-1] = QQ.one
    return lst


@given(f=dup_zero_const())
def test_trancendental_function_zero(f):
    Rs = power_series_ring(QQ, 25)
    Rp, x = ring("x", QQ)

    s = Rs(f[::-1])
    p = Rp.from_list(f)

    assert Rs.to_dense(Rs.exp(s)) == (rs_exp(p, x, 25)).to_dense()

    assert Rs.to_dense(Rs.atan(s)) == (rs_atan(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.atanh(s)) == (rs_atanh(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.asin(s)) == (rs_asin(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.asinh(s)) == (rs_asinh(p, x, 25)).to_dense()

    assert Rs.to_dense(Rs.tan(s)) == (rs_tan(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.tanh(s)) == (rs_tanh(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.sin(s)) == (rs_sin(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.sinh(s)) == (rs_sinh(p, x, 25)).to_dense()
    assert Rs.to_dense(Rs.cos(s)) == (rs_cos(p, x, 25)).to_dense()
    # assert Rs.to_dense(Rs.cosh(s)) == (rs_cosh(p, x, 25)).to_dense()


@given(f=dup_one_const())
def test_trancendental_function_one(f):
    Rs = power_series_ring(QQ, 25)
    Rp, x = ring("x", QQ)

    s = Rs(f[::-1])
    p = Rp.from_list(f)

    assert Rs.to_dense(Rs.log(s)) == (rs_log(p, x, 25)).to_dense()
