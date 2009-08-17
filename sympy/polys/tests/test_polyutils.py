"""Tests for useful utilities for higher level polynomial classes. """

from sympy import S, I, sqrt, symbols

from sympy.polys.polyutils import (
    _analyze_power,
    _dict_from_basic_if_gens,
    _dict_from_basic_no_gens,
    dict_from_basic,
    basic_from_dict,
)

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
)

def test__analyze_power():
    x, y = symbols('x,y')

    assert _analyze_power(x, S(1)) == (x, S(1))
    assert _analyze_power(x, S(2)) == (x, S(2))
    assert _analyze_power(x, -S(1)) == (x**(-1), S(1))
    assert _analyze_power(x, -S(2)) == (x**(-1), S(2))

    assert _analyze_power(x, S(1)/3) == (x**(S(1)/3), S(1))
    assert _analyze_power(x, S(2)/3) == (x**(S(1)/3), S(2))
    assert _analyze_power(x, -S(1)/3) == (x**(-S(1)/3), S(1))
    assert _analyze_power(x, -S(2)/3) == (x**(-S(1)/3), S(2))

    assert _analyze_power(x, y) == (x**y, S(1))
    assert _analyze_power(x, -y) == (x**(-y), S(1))
    assert _analyze_power(x, 2*y) == (x**y, S(2))
    assert _analyze_power(x, -2*y) == (x**(-y), S(2))

    assert _analyze_power(x, y/3) == (x**(y/3), S(1))
    assert _analyze_power(x, -y/3) == (x**(-y/3), S(1))
    assert _analyze_power(x, 2*y/3) == (x**(y/3), S(2))
    assert _analyze_power(x, -2*y/3) == (x**(-y/3), S(2))

    assert _analyze_power(x, S(1.0)) == (x**S(1.0), S(1))
    assert _analyze_power(x, S(2.0)) == (x**S(2.0), S(1))
    assert _analyze_power(x, -S(1.0)) == (x**(-S(1.0)), S(1))
    assert _analyze_power(x, -S(2.0)) == (x**(-S(2.0)), S(1))

    assert _analyze_power(x, S(1.0)*y) == (x**(S(1.0)*y), S(1))
    assert _analyze_power(x, S(2.0)*y) == (x**(S(2.0)*y), S(1))
    assert _analyze_power(x, -S(1.0)*y) == (x**(-S(1.0)*y), S(1))
    assert _analyze_power(x, -S(2.0)*y) == (x**(-S(2.0)*y), S(1))

