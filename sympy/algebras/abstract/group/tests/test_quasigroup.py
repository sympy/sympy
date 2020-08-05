from sympy import Set, BinaryOperator
from sympy.algebras import (
    LeftQuasigroup, RightQuasigroup, Quasigroup
)
from sympy.testing.pytest import raises

A = Set('A')
a, b, c = [A.element(i) for i in 'abc']

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

class G(F_1, F_2):
    pass
g = G()

def test_LeftQuasigroup():
    # LeftQuasigroup's operator must have left division
    raises(TypeError, lambda: LeftQuasigroup('Q', (A,), (f_2,)))

    Q = LeftQuasigroup('Q', (A,), (f_1,))
    Q_op = Q.operator
    Q_ld = Q.left_division
    assert Q_op(a, Q_ld(a, b), evaluate=True) == b
    assert Q_ld(a, Q_op(a, b), evaluate=True) == b

def test_RightQuasigroup():
    # RightQuasigroup's operator must have left division
    raises(TypeError, lambda: RightQuasigroup('Q', (A,), (f_1,)))

    Q = RightQuasigroup('Q', (A,), (f_2,))
    Q_op = Q.operator
    Q_rd = Q.right_division
    assert Q_op(Q_rd(a, b), b, evaluate=True) == a
    assert Q_rd(Q_op(a, b), b, evaluate=True) == a

def test_Quasigroup():
    # Quasigroup's operator must have left and right division
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f_1,)))
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f_2,)))

    Q = Quasigroup('Q', (A,), (g,))
    Q_op = Q.operator
    Q_ld, Q_rd = Q.left_division, Q.right_division
    assert Q_op(a, Q_ld(a, b), evaluate=True) == b
    assert Q_ld(a, Q_op(a, b), evaluate=True) == b
    assert Q_op(Q_rd(a, b), b, evaluate=True) == a
    assert Q_rd(Q_op(a, b), b, evaluate=True) == a
