from sympy.external import import_module

theano = import_module('theano')
if theano:
    ts = theano.scalar
    tt = theano.tensor
    xt, yt, zt = [tt.scalar(name, 'floatX') for name in 'xyz']
else:
    #bin/test will not execute any tests now
    disabled = True

import sympy
from sympy import S
sy = sympy
from sympy.abc import x, y, z, a, b, c
from sympy.printing.theanocode import (theano_code, dim_handling,
        theano_function)

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
    assert isinstance(xt, (tt.TensorVariable, ts.ScalarVariable))
    assert xt.name == x.name

    assert theano_code(x, broadcastables={x: (False,)}).broadcastable == (False,)
    assert theano_code(x, broadcastables={x: (False,)}).name == x.name

def test_add():
    expr = x + y
    comp = theano_code(expr)
    assert comp.owner.op == theano.tensor.add

    comp = theano_code(expr, broadcastables={x: (False,), y: (False,)})
    assert comp.broadcastable == (False,)

    comp = theano_code(expr, broadcastables={x: (False, True), y: (False, False)})
    assert comp.broadcastable == (False, False)


def test_trig():
    assert theq(theano_code(sympy.sin(x)), tt.sin(xt))
    assert theq(theano_code(sympy.tan(x)), tt.tan(xt))

def test_many():
    expr = sy.exp(x**2 + sy.cos(y)) * sy.log(2*z)
    comp = theano_code(expr)
    expected = tt.exp(xt**2 + tt.cos(yt)) * tt.log(2*zt)
    # assert theq(comp, expected)

def test_dtype():
    assert theano_code(x, dtypes={x: 'float32'}).type.dtype == 'float32'
    assert theano_code(x, dtypes={x: 'float64'}).type.dtype == 'float64'
    assert theano_code(x+1, dtypes={x: 'float32'}).type.dtype == 'float32'
    assert theano_code(x+y, dtypes={x: 'float64', y: 'float32'}).type.dtype == 'float64'

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

def test_Transpose():
    X = sympy.MatrixSymbol('X', 4, 4)
    assert isinstance(theano_code(X.T).owner.op, tt.DimShuffle)

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
    assert dim_handling([x], broadcastables={x: (False,)}) == {x: (False,)}

def test_Rationals():
    assert theq(theano_code(sympy.Integer(2) / 3), tt.true_div(2, 3))
    assert theq(theano_code(S.Half), tt.true_div(1, 2))

def test_Integers():
    assert theano_code(sympy.Integer(3)) == 3

def test_factorial():
    n = sympy.Symbol('n')
    assert theano_code(sympy.factorial(n))

def test_Derivative():
    assert theq(theano_code(sy.Derivative(sy.sin(x), x, evaluate=False)),
                theano.grad(tt.sin(xt), xt))

def test_theano_function_simple():
    f = theano_function([x, y], [x+y])
    assert f(2, 3) == 5


def test_theano_function_numpy():
    import numpy as np
    f = theano_function([x, y], [x+y], dim=1)
    assert np.linalg.norm(f([1, 2], [3, 4]) - np.asarray([4, 6])) < 1e-9

    f = theano_function([x, y], [x+y], dtypes={x: 'float64', y: 'float64'},
                                     dim=1)
    xx = np.arange(3).astype('float64')
    yy = 2*np.arange(3).astype('float64')
    assert np.linalg.norm(f(xx, yy) - 3*np.arange(3)) < 1e-9
