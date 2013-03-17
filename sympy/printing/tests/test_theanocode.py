from sympy.printing.theanocode import theano_code
import sympy
import theano
ts = theano.scalar
tt = theano.tensor
sy = sympy

from sympy.abc import x, y, z, a, b, c

xt, yt, zt = map(ts.Scalar('floatX'), 'xyz')

def theq(a, b):
    """ theano equality """
    astr = theano.printing.debugprint(a, file='str')
    bstr = theano.printing.debugprint(b, file='str')

    if not astr == bstr:
        print
        print astr
        print bstr

    return astr == bstr

def test_symbol():
    xt = theano_code(x)
    assert isinstance(xt, ts.ScalarVariable)
    assert xt.name == x.name

def test_add():
    expr = x + y
    comp = theano_code(expr)
    assert comp.owner.op == theano.scalar.add

def test_trig():
    assert theq(theano_code(sympy.sin(x)), ts.sin(xt))
    assert theq(theano_code(sympy.tan(x)), ts.tan(xt))

def test_many():
    expr = sy.exp(x**2 + sy.cos(y)) * sy.log(2*z)
    comp = theano_code(expr)
    expected = ts.exp(xt**2 + ts.cos(yt)) * ts.log(2*zt)
    # assert theq(comp, expected)

def test_dtype():
    assert theano_code(x, {x: 'float32'}).type.dtype == 'float32'
    assert theano_code(x, {x: 'float64'}).type.dtype == 'float64'
    assert theano_code(x+1, {x: 'float32'}).type.dtype == 'float32'
    assert theano_code(x+y, {x: 'float64', y: 'float32'}).type.dtype == 'float64'

def test_MatrixSymbol():
    X = sympy.MatrixSymbol('X', 4, 5)
    Xt = theano_code(X)
    assert isinstance(Xt, tt.TensorVariable)
    assert Xt.broadcastable == (False, False)

def test_MatMul():
    X = sympy.MatrixSymbol('X', 4, 4)
    Y = sympy.MatrixSymbol('X', 4, 4)
    Z = sympy.MatrixSymbol('X', 4, 4)
    expr = X*Y*Z
    assert isinstance(theano_code(expr).owner.op, tt.Dot)
