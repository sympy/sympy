from sympy import Symbol, Integer
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.hilbert import HilbertSpace

x = Symbol('x')
y = Symbol('y')

def test_qexpr_new():
    q = QExpr(0)
    assert q.label == (0,)
    assert q.hilbert_space == HilbertSpace()
    assert q.is_commutative == False

    q = QExpr(0,1)
    assert q.label == (Integer(0),Integer(1))

    q = QExpr._new_rawargs(HilbertSpace(), Integer(0), Integer(1))
    assert q.label == (Integer(0),Integer(1))
    assert q.hilbert_space == HilbertSpace()


def test_qexpr_commutative():
    q1 = QExpr(x)
    q2 = QExpr(y)
    assert q1.is_commutative == False
    assert q2.is_commutative == False
    assert q1*q2 != q2*q1

    q = QExpr._new_rawargs(0,1,HilbertSpace())
    assert q.is_commutative == False


def test_qexpr_subs():
    q1 = QExpr(x,y)
    assert q1.subs(x, y) == QExpr(y,y)
    assert q1.subs({x:1,y:2}) == QExpr(1,2)


    