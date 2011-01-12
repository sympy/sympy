from sympy.physics.qbit import *
from sympy import symbols, Rational
from sympy.core.numbers import *
from sympy.functions.elementary import *
from sympy.physics.shor import *
import random
x, y = symbols('xy')

epsilon = .000001

def test_controlledMod():
    assert apply_gates(controlledMod(4, 2, 2)*Qbit(0,0,1,0,0,0,0,0)) ==\
    Qbit(0,0,1,0,0,0,0,0)
    assert apply_gates(controlledMod(5, 5, 7)*Qbit(0,0,1,0,0,0,0,0,0,0)) ==\
    Qbit(0,0,1,0,0,0,0,0,1,0)
    assert apply_gates(controlledMod(3, 2, 3)*Qbit(0,1,0,0,0,0)) ==\
    Qbit(0,1,0,0,0,1)

def test_continuedFrac():
    assert continuedFraction(3245, 10000) == [0,3,12,4,13]
    assert continuedFraction(1932, 2568) == [0, 1, 3, 26, 2]
    assert continuedFraction(6589, 2569) == [2, 1, 1, 3, 2, 1, 3, 1, 23]
    assert getr(513, 1024, 10) == 2
    assert getr(169, 1024, 11) == 6
    assert getr(314, 4096, 16) == 13
    assert arr(5,3) == [1,0,1]
    assert arr(4,3) == [1,0,0]
    assert arr(8,5) == [0,1,0,0,0]

