# This testfile tests SymPy <-> NumPy compatibility

# Don't test any SymPy features here. Just pure interaction with NumPy.
# Always write regular SymPy tests for anything, that can be tested in pure
# Python (without numpy). Here we test everything, that a user may need when
# using SymPy with NumPy

try:
    from numpy import array, ndarray
except ImportError:
    #py.test will not execute any tests now
    disabled = True


import sys
sys.path.append("..")
from sympy import Rational, Symbol, list2numpy

def test_basics():
    one = Rational(1)
    zero = Rational(0)
    x = Symbol("x")
    assert array(1) == array(one)
    assert array([one]) == array([one])
    assert array([x]) == array([x])
    assert array(x) == array(Symbol("x"))
    assert array(one+x) == array(1+x)

    X = array([one, zero, zero])
    assert X == array([one, zero, zero])
    assert X == array([one, 0, 0])

def test_arrays():
    one = Rational(1)
    zero = Rational(0)
    X = array([one, zero, zero])
    #fails:
    #Y = one*X
    X = array([Symbol("a")+Rational(1,2)])
    Y = X+X
    #fails:
    #assert Y == numpy.array([1+2*sympy.Symbol("a")])
    #fails:
    #Y = Y + 1
    #fails:
    #Y = X-X

def test_conversion1():
    x = Symbol("x")
    a = list2numpy([x**2, x])
    #looks like an array? 
    assert isinstance(a, ndarray)
    assert a[0] == x**2
    assert a[1] == x
    assert len(a) == 2
    #yes, it's the array

def test_conversion2():
    x = Symbol("x")
    a = 2*list2numpy([x**2, x])
    b = list2numpy([2*x**2, 2*x])
    assert (a == b).all()

    one = Rational(1)
    zero = Rational(0)
    X = list2numpy([one, zero, zero])
    #Y = one*X
    X = list2numpy([Symbol("a")+Rational(1,2)])
    Y = X+X
    #fails:
    #assert Y == numpy.array([1+2*sympy.Symbol("a")])
    #fails:
    #Y = Y + 1
    #fails:
    #Y = X-X
