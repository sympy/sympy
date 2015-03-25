from sympy.external import import_module
from sympy.utilities.pytest import raises, SKIP
from sympy.core.compatibility import range

theano = import_module('theano')
if theano:
    import numpy as np
    ts = theano.scalar
    tt = theano.tensor
    xt, yt, zt = [tt.scalar(name, 'floatX') for name in 'xyz']
else:
    #bin/test will not execute any tests now
    disabled = True

import sympy
from sympy import S
sy = sympy
from sympy.abc import x, y, z
from sympy.printing.theanocode import (theano_code, dim_handling,
        theano_function)

def fgraph_of(*exprs):
    """ Transform SymPy expressions into Theano Computation """
    outs = list(map(theano_code, exprs))
    ins = theano.gof.graph.inputs(outs)
    ins, outs = theano.gof.graph.clone(ins, outs)
    return theano.gof.FunctionGraph(ins, outs)

def theano_simplify(fgraph):
    """ Simplify a Theano Computation """
    mode = theano.compile.get_default_mode().excluding("fusion")
    fgraph = fgraph.clone()
    mode.optimizer.optimize(fgraph)
    return fgraph


def theq(a, b):
    """ theano equality """
    astr = theano.printing.debugprint(a, file='str')
    bstr = theano.printing.debugprint(b, file='str')

    if not astr == bstr:
        print()
        print(astr)
        print(bstr)

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
    simp = lambda expr: theano_simplify(fgraph_of(expr))
    assert theq(simp(theano_code(sy.Derivative(sy.sin(x), x, evaluate=False))),
                simp(theano.grad(tt.sin(xt), xt)))

def test_theano_function_simple():
    f = theano_function([x, y], [x+y])
    assert f(2, 3) == 5

def test_theano_function_numpy():
    f = theano_function([x, y], [x+y], dim=1,
                        dtypes={x: 'float64', y: 'float64'})
    assert np.linalg.norm(f([1, 2], [3, 4]) - np.asarray([4, 6])) < 1e-9

    f = theano_function([x, y], [x+y], dtypes={x: 'float64', y: 'float64'},
                                     dim=1)
    xx = np.arange(3).astype('float64')
    yy = 2*np.arange(3).astype('float64')
    assert np.linalg.norm(f(xx, yy) - 3*np.arange(3)) < 1e-9

def test_theano_function_kwargs():
    import numpy as np
    f = theano_function([x, y, z], [x+y], dim=1, on_unused_input='ignore',
            dtypes={x: 'float64', y: 'float64', z: 'float64'})
    assert np.linalg.norm(f([1, 2], [3, 4], [0, 0]) - np.asarray([4, 6])) < 1e-9

    f = theano_function([x, y, z], [x+y],
                        dtypes={x: 'float64', y: 'float64', z: 'float64'},
                        dim=1, on_unused_input='ignore')
    xx = np.arange(3).astype('float64')
    yy = 2*np.arange(3).astype('float64')
    zz = 2*np.arange(3).astype('float64')
    assert np.linalg.norm(f(xx, yy, zz) - 3*np.arange(3)) < 1e-9

def test_slice():
    assert theano_code(slice(1, 2, 3)) == slice(1, 2, 3)
    assert str(theano_code(slice(1, x, 3), dtypes={x: 'int32'})) ==\
           str(slice(1, xt, 3))

def test_MatrixSlice():
    n = sympy.Symbol('n', integer=True)
    X = sympy.MatrixSymbol('X', n, n)

    Y = X[1:2:3, 4:5:6]
    Yt = theano_code(Y)
    assert tuple(Yt.owner.op.idx_list) == (slice(1,2,3), slice(4,5,6))
    assert Yt.owner.inputs[0] == theano_code(X)

    k = sympy.Symbol('k')
    kt = theano_code(k, dtypes={k: 'int32'})
    start, stop, step = 4, k, 2
    Y = X[start:stop:step]
    Yt = theano_code(Y, dtypes={n: 'int32', k: 'int32'})
    # assert Yt.owner.op.idx_list[0].stop == kt

def test_BlockMatrix():
    n = sympy.Symbol('n', integer=True)
    A = sympy.MatrixSymbol('A', n, n)
    B = sympy.MatrixSymbol('B', n, n)
    C = sympy.MatrixSymbol('C', n, n)
    D = sympy.MatrixSymbol('D', n, n)
    At, Bt, Ct, Dt = map(theano_code, (A, B, C, D))
    Block = sympy.BlockMatrix([[A, B], [C, D]])
    Blockt = theano_code(Block)
    solutions = [tt.join(0, tt.join(1, At, Bt), tt.join(1, Ct, Dt)),
                 tt.join(1, tt.join(0, At, Ct), tt.join(0, Bt, Dt))]
    assert any(theq(Blockt, solution) for solution in solutions)

@SKIP
def test_BlockMatrix_Inverse_execution():
    k, n = 2, 4
    dtype = 'float32'
    A = sympy.MatrixSymbol('A', n, k)
    B = sympy.MatrixSymbol('B', n, n)
    inputs = A, B
    output = B.I*A

    cutsizes = {A: [(n//2, n//2), (k//2, k//2)],
                B: [(n//2, n//2), (n//2, n//2)]}
    cutinputs = [sympy.blockcut(i, *cutsizes[i]) for i in inputs]
    cutoutput = output.subs(dict(zip(inputs, cutinputs)))

    dtypes = dict(zip(inputs, [dtype]*len(inputs)))
    f = theano_function(inputs, [output], dtypes=dtypes, cache={})
    fblocked = theano_function(inputs, [sympy.block_collapse(cutoutput)],
                               dtypes=dtypes, cache={})

    ninputs = [np.random.rand(*x.shape).astype(dtype) for x in inputs]
    ninputs = [np.arange(n*k).reshape(A.shape).astype(dtype),
               np.eye(n).astype(dtype)]
    ninputs[1] += np.ones(B.shape)*1e-5

    assert np.allclose(f(*ninputs), fblocked(*ninputs), rtol=1e-5)

def test_DenseMatrix():
    t = sy.Symbol('theta')
    for MatrixType in [sy.Matrix, sy.ImmutableMatrix]:
        X = MatrixType([[sy.cos(t), -sy.sin(t)], [sy.sin(t), sy.cos(t)]])
        tX = theano_code(X)
        assert isinstance(tX, tt.TensorVariable)
        assert tX.owner.op == tt.join

def test_AppliedUndef():
    t = sy.Symbol('t')
    f = sy.Function('f')
    ft = theano_code(f(t))
    assert isinstance(ft, tt.TensorVariable)
    assert ft.name == 'f_t'

def test_bad_keyword_args_raise_error():
    raises(Exception, lambda : theano_function([x], [x+1], foobar=3))

def test_cache():
    sx = sy.Symbol('x')
    cache = {}
    tx = theano_code(sx, cache=cache)
    assert theano_code(sx, cache=cache) is tx
    assert theano_code(sx, cache={}) is not tx

def test_Piecewise():
    # A piecewise linear
    xt, yt = theano_code(x), theano_code(y)
    expr = sy.Piecewise((0, x<0), (x, x<2), (1, True))  # ___/III
    result = theano_code(expr)
    assert result.owner.op == tt.switch

    expected = tt.switch(xt<0, 0, tt.switch(xt<2, xt, 1))
    assert theq(result, expected)

    expr = sy.Piecewise((x, x < 0))
    result = theano_code(expr)
    expected = tt.switch(xt < 0, xt, np.nan)
    assert theq(result, expected)

def test_Relationals():
    xt, yt = theano_code(x), theano_code(y)
    assert theq(theano_code(x > y), xt > yt)
    assert theq(theano_code(x < y), xt < yt)
    assert theq(theano_code(x >= y), xt >= yt)
    assert theq(theano_code(x <= y), xt <= yt)
    
def test_complexfunctions():
    xt, yt = theano_code(x, dtypes={x:'complex128'}), theano_code(y, dtypes={y: 'complex128'})
    from sympy import conjugate
    from theano.tensor import as_tensor_variable as atv
    from theano.tensor import complex as cplx
    assert theq(theano_code(y*conjugate(x)), yt*(xt.conj()))
    assert theq(theano_code((1+2j)*x), xt*(atv(1.0)+atv(2.0)*cplx(0,1)))

def test_constantfunctions():
    tf = theano_function([],[1+1j])
    assert(tf()==1+1j)
