from sympy import Symbol, Function, sin, exp, cosh, Matrix, Number
import pickle
from sympy.utilities.pytest import XFAIL

def test_symbolpickling():
    x = Symbol("x")
    x_after = pickle.loads(pickle.dumps(x))
    assert x==x_after
    assert dir(x)==dir(x_after)

def test_numberpickling():
    for n in (Number(0), Number(1), Number(3), Number(100)/4, Number("1.1")):
        n_after = pickle.loads(pickle.dumps(n))
        assert n==n_after
        assert dir(n)==dir(n_after)

def test_functionpickling():
    for f in (sin, exp, cosh):
        f_after = pickle.loads(pickle.dumps(f))
        assert f==f_after
        assert dir(f)==dir(f_after)

@XFAIL
def test_dynamicfunctionpickling():
    f = Function("f")
    f_after = pickle.loads(pickle.dumps(f))
    assert f==f_after
    assert dir(f)==dir(f_after)

def test_expressionpickling():
    x = Symbol("x")
    y = Symbol("y")
    for expr in (x**2, x+y, x*y, x/y, x*sin(y)):
        expr_after = pickle.loads(pickle.dumps(expr))
        assert expr==expr_after
        assert dir(expr)==dir(expr_after)

def test_matrixpickling():
    m = Matrix([1,2,3])
    m_after = pickle.loads(pickle.dumps(m))
    assert m==m_after
    assert dir(m)==dir(m_after)
