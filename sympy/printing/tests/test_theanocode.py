from sympy.external import import_module

theano = import_module('theano')
if theano:
    ts = theano.scalar
    tt = theano.tensor
else:
    #bin/test will not execute any tests now
    disabled = True

import sympy
from sympy import S
sy = sympy
from sympy.abc import x, y, z, a, b, c
from sympy.printing.theanocode import theano_code, dim_handling, tensor_wrap

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

def test_MatAdd():
    X = sympy.MatrixSymbol('X', 4, 4)
    Y = sympy.MatrixSymbol('X', 4, 4)
    Z = sympy.MatrixSymbol('X', 4, 4)
    expr = X+Y+Z
    assert isinstance(theano_code(expr).owner.op, tt.Elemwise)

def test_symbols_are_created_once():
    expr = x**x
    comp = theano_code(expr)

    assert theq(comp, xt**xt)

def test_dim_handling():
    assert dim_handling([x], dim=2) == {x: (False, False)}
    assert dim_handling([x, y], dims={x: 1, y: 2}) == {x: (False, True),
                                                       y: (False, False)}
    assert dim_handling([x], broadcastable={x: (False,)}) == {x: (False,)}
    assert dim_handling([x], dims={'x': 1}, keys=['x']) == {x: (False,)}

def test_tensor_wrap():
    [Xt], Xtp1 = tensor_wrap([xt], [xt+1], dim=2)
    assert isinstance(Xt,   tt.TensorVariable)
    assert isinstance(Xtp1, tt.TensorVariable)
    assert Xtp1.type.broadcastable == (False, False)

def test_Rationals():
    assert theq(theano_code(sympy.Integer(2) / 3), ts.true_div(2, 3))
    assert theq(theano_code(S.Half), ts.true_div(1, 2))
