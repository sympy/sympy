from sympy import Symbol, Tuple
from sympy.physics.qexpr import QExpr
from sympy.physics.hilbert import HilbertSpace

x = Symbol('x')
y = Symbol('y')

def test_qexpr_new():
    q = QExpr(0)
    assert q.label == Tuple(0)
    assert q.hilbert_space == HilbertSpace()
    assert q.is_commutative == False

    q = QExpr((0,1))
    assert q.label == Tuple(0,1)

    q = QExpr(Tuple(0,1))
    assert q.label == Tuple(0,1)

    q = QExpr._new_rawargs(HilbertSpace(), Tuple(0,1))
    assert q.label == Tuple(0,1)
    assert q.hilbert_space == HilbertSpace()

def test_qexpr_commutative():
    q1 = QExpr(x)
    q2 = QExpr(y)
    assert q1*q2 != q2*q1

def test_qexpr_subs():
    q1 = QExpr((x,y))
    assert q1.subs(x, y) == QExpr((y,y))
    assert q1.subs({x:1,y:2}) == QExpr((1,2))

    