from sympy.sandbox.core import *

x = Symbol('x')
y = Symbol('y')

def test_split_add():
    assert (x).split('+') == [x]
    assert (3+x).split('+') == sorted([3, x])
    assert (3+x+y).split('+') == sorted([3, x, y])
    assert (1+4*x+y).split('+') == sorted([1, 4*x, y])
    assert (3*x).split('+') == [3*x]
    assert (3*x*y).split('+') == [3*x*y]

def test_split_mul():
    assert (x).split('*') == [x]
    assert (3+x).split('*') == [3+x]
    assert (3*x).split('*') == sorted([3, x])
    assert (3*x*y).split('*') == sorted([3, x, y])
    assert (3*x**2*y).split('*') == sorted([3, x**2, y])

def test_split_pow():
    assert (x).split('**') == [x, 1]
    assert (3+x).split('**') == [3+x, 1]
    assert (3*x).split('**') == [3*x, 1]
    assert (x**4).split('**') == [x, 4]
