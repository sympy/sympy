from sympy.concrete.guess import (
            find_simple_recurrence,
            rationalize
        )
from sympy import Function, var, sympify


def test_find_simple_recurrence():
    from mpmath import fib
    a = Function('a')
    n = var('n')
    assert find_simple_recurrence( [ fib(k) for k in range(12) ] ) == (
        -a(n) - a(n + 1) + a(n + 2) )

    f = Function('a')
    i = var('n')
    a = [1, 1, 1]
    for k in range(15): a.append(5*a[-1]-3*a[-2]+8*a[-3])
    assert find_simple_recurrence(a, A=f, N=i) == (
        -8*f(i) + 3*f(i + 1) - 5*f(i + 2) + f(i + 3) )


def test_rationalize():
    from mpmath import cos, pi, mpf
    assert rationalize( cos(pi/3) ) == sympify("1/2")
    assert rationalize( mpf("0.333333333333333") ) == sympify("1/3")
    assert rationalize( mpf("-0.333333333333333") ) == sympify("-1/3")
    assert rationalize( pi, maxcoeff = 250 ) == sympify("355/113")
