from sympy import Set, BinaryOperator
from sympy.algebras import (
    LeftQuasigroup, RightQuasigroup, Quasigroup
)
from sympy.testing.pytest import raises

A = Set('A')

class F_1(BinaryOperator):
    is_left_divisible=True
    domain = A**2
    codomain = A
f_1 = F_1()

class F_2(BinaryOperator):
    is_right_divisible=True
    domain = A**2
    codomain = A
f_2 = F_2()

def test_LeftQuasigroup():
    # LeftQuasigroup's operator must have left division
    raises(TypeError, lambda: LeftQuasigroup('Q', (A,), (f_2,)))

def test_RightQuasigroup():
    # RightQuasigroup's operator must have left division
    raises(TypeError, lambda: RightQuasigroup('Q', (A,), (f_1,)))

def test_Quasigroup():
    # Quasigroup's operator must have left and right division
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f_1,)))
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f_2,)))
