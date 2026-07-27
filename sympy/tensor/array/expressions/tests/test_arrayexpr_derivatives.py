from __future__ import annotations
from sympy import (exp, log, Matrix, Array, derive_by_array, ImmutableMatrix,
    ImmutableDenseNDimArray, Transpose, Inverse, Determinant, MatPow, ZeroMatrix,
    Trace, atan2)
from sympy.core.symbol import symbols
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.matrices.expressions.diagonal import DiagMatrix, DiagonalMatrix
from sympy.matrices.expressions.matexpr import MatrixExpr, MatrixSymbol, MatrixElement
from sympy.matrices.expressions.special import Identity, MatrixUnit
from sympy.testing.pytest import raises
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.tensor.array.expressions.array_expressions import ArraySymbol, ArrayTensorProduct, \
    PermuteDims, ArrayDiagonal, ArrayElementwiseApplyFunc, ArrayContraction, _permute_dims, \
    Reshape, ZeroArray
from sympy.tensor.array.expressions.arrayexpr_derivatives import array_derive

k = symbols("k")

I = Identity(k)
X = MatrixSymbol("X", k, k)
x = MatrixSymbol("x", k, 1)

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)

A1 = ArraySymbol("A", (3, 2, k))


def test_arrayexpr_derivatives1():

    res = array_derive(X, X)
    assert res == PermuteDims(ArrayTensorProduct(I, I), [0, 2, 1, 3])

    cg = ArrayTensorProduct(A, X, B)
    res = array_derive(cg, X)
    assert res == _permute_dims(
        ArrayTensorProduct(I, A, I, B),
        [0, 4, 2, 3, 1, 5, 6, 7])

    cg = ArrayContraction(X, (0, 1))
    res = array_derive(cg, X)
    assert res == ArrayContraction(ArrayTensorProduct(I, I), (1, 3))

    cg = ArrayDiagonal(X, (0, 1))
    res = array_derive(cg, X)
    assert res == ArrayDiagonal(ArrayTensorProduct(I, I), (1, 3))

    cg = ElementwiseApplyFunction(sin, X)
    res = array_derive(cg, X)
    assert res.dummy_eq(ArrayDiagonal(
        ArrayTensorProduct(
            ElementwiseApplyFunction(cos, X),
            I,
            I
        ), (0, 3), (1, 5)))

    cg = ArrayElementwiseApplyFunc(sin, X)
    res = array_derive(cg, X)
    assert res.dummy_eq(ArrayDiagonal(
        ArrayTensorProduct(
            I,
            I,
            ArrayElementwiseApplyFunc(cos, X)
        ), (1, 4), (3, 5)))

    res = array_derive(A1, A1)
    assert res == PermuteDims(
        ArrayTensorProduct(Identity(3), Identity(2), Identity(k)),
        [0, 2, 4, 1, 3, 5]
    )

    cg = ArrayElementwiseApplyFunc(sin, A1)
    res = array_derive(cg, A1)
    assert res.dummy_eq(ArrayDiagonal(
        ArrayTensorProduct(
            Identity(3), Identity(2), Identity(k),
            ArrayElementwiseApplyFunc(cos, A1)
        ), (1, 6), (3, 7), (5, 8)
    ))

    cg = Reshape(A, (k**2,))
    res = array_derive(cg, A)
    assert res == Reshape(PermuteDims(ArrayTensorProduct(I, I), [0, 2, 1, 3]), (k, k, k**2))


def test_array_derive_with_common_matrices():
    expr = Matrix([[k, k**2], [sin(k), exp(k)]])
    assert array_derive(expr, Matrix([k])) == Array([[[[1, 2*k], [cos(k), exp(k)]]]])

    I2 = Identity(2)
    expr = MatrixSymbol("M", 2, 2)
    assert array_derive(expr, expr) == PermuteDims(ArrayTensorProduct(I2, I2), [0, 2, 1, 3])


def test_array_derive_permutedims_non_matrix():
    # The PermuteDims derivative rule must not assume a rank-2 operand: it
    # prepends rank(x) axes, so those must be fixed and the permutation shifted
    # by rank(x).  Previously this was hardcoded for a matrix (rank 2) and
    # crashed for any other rank.  See the ArrayContraction/ArrayDiagonal rules,
    # which already handle this correctly.
    A3 = ArraySymbol("A", (2, 2, 2))
    expr = PermuteDims(A3, [1, 2, 0])
    res = array_derive(expr, A3)
    # the derivative of a rank-3 array w.r.t. itself has rank 6
    assert res.as_explicit() == derive_by_array(expr.as_explicit(), A3.as_explicit())


def test_array_derive_matrix_element_of_composite_parent():
    # The MatrixElement rule used to return MatrixUnit(...) for *any*
    # MatrixElement, i.e. it silently pretended ``expr.parent`` was ``x``.
    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)

    # fast path: element of x itself
    assert array_derive(M[1, 0], M) == MatrixUnit(2, 2, 1, 0)

    # (M*N)[0, 0] would auto-expand for concrete shapes, so build the
    # unevaluated MatrixElement directly:
    res = array_derive(MatrixElement(M*N, 0, 0), M)
    assert res.shape == (2, 2)
    Mexp = ImmutableMatrix(2, 2, symbols("m:4"))
    Nexp = ImmutableMatrix(2, 2, symbols("n:4"))
    repl = dict(zip(list(M), list(Mexp)))
    repl.update(dict(zip(list(N), list(Nexp))))
    expected = derive_by_array((Mexp*Nexp)[0, 0], ImmutableDenseNDimArray(Mexp))
    assert ImmutableDenseNDimArray(res.as_explicit().subs(repl)) == expected


def test_array_derive_zero_expressions():
    M = MatrixSymbol("M", 2, 2)
    v = ArraySymbol("v", (2,))
    # ZeroArray used to fall into the scalar Expr rule (wrong shape),
    # ZeroMatrix used to raise NotImplementedError:
    assert array_derive(ZeroArray(2, 2), M) == ZeroArray(2, 2, 2, 2)
    assert array_derive(ZeroArray(3,), v) == ZeroArray(2, 3)
    assert array_derive(ZeroMatrix(2, 2), M) == ZeroArray(2, 2, 2, 2)
    assert array_derive(ZeroMatrix(2, 3), v) == ZeroArray(2, 2, 3)


def test_array_derive_matrix_rules_non_matrix_x():
    # Transpose/Inverse/Determinant/MatPow rules used to hardcode rank-2 x
    # and crashed for any other rank of the differentiation variable.
    M = MatrixSymbol("M", 2, 2)
    v = ArraySymbol("v", (2,))
    w = ArraySymbol("w", (2, 2, 2))
    for x in (v, w):
        for expr, shape in [(Transpose(M), (2, 2)), (Inverse(M), (2, 2)),
                            (Determinant(M), ()), (MatPow(M, 2), (2, 2))]:
            res = array_derive(expr, x)
            assert res.shape == x.shape + shape
            assert ImmutableDenseNDimArray(res.as_explicit()) == \
                ImmutableDenseNDimArray.zeros(*(x.shape + shape))


def test_array_derive_matrix_rules_rank2_x_numeric():
    # regression: rank-2 x must still match derive_by_array
    M = MatrixSymbol("M", 2, 2)
    Mexp = ImmutableMatrix(2, 2, symbols("m:4"))
    for expr, expl in [(Transpose(M), ImmutableDenseNDimArray(Mexp.T)),
                       (Inverse(M), ImmutableDenseNDimArray(Mexp.inv())),
                       (Determinant(M), Mexp.det())]:
        res = array_derive(expr, M)
        got = ImmutableDenseNDimArray(res.as_explicit().subs(M, Mexp).doit())
        expected = derive_by_array(expl, ImmutableDenseNDimArray(Mexp))
        diff = (got - expected).applyfunc(lambda t: t.simplify())
        assert diff == ImmutableDenseNDimArray.zeros(*expected.shape)


def test_array_derive_ndimarray_symbolic_x():
    # The NDimArray rule used to return a permanently unevaluated Derivative
    # when x was symbolic; it now converts x to explicit form like the
    # MatrixBase rule does.
    M = MatrixSymbol("M", 2, 2)
    arr = ImmutableDenseNDimArray([[M[0, 0], 0], [0, 0]])
    res = array_derive(arr, M)
    assert ImmutableDenseNDimArray(res) == \
        ImmutableDenseNDimArray(array_derive(arr.tomatrix(), M))
    assert ImmutableDenseNDimArray(res) == derive_by_array(arr, M.as_explicit())
    # constant arrays give a ZeroArray instead of a stuck Derivative:
    assert array_derive(ImmutableDenseNDimArray([[1, 2], [3, 4]]), M) == ZeroArray(2, 2, 2, 2)
    # ArraySymbol x:
    v = ArraySymbol("v", (2,))
    arr = ImmutableDenseNDimArray([v[0], v[1]**2])
    assert ImmutableDenseNDimArray(array_derive(arr, v)) == \
        derive_by_array(arr, v.as_explicit())


def _check_numeric(expr, wrt, replmap, numreplmap={}):
    # Compare the array derivative against entrywise differentiation of the
    # explicit form.  ``replmap`` maps matrix symbols to explicit matrices;
    # ``numreplmap`` optionally substitutes numbers at the end (to keep the
    # simplification of expressions with matrix inverses cheap).
    res = array_derive(expr, wrt)
    got = res.as_explicit().subs(replmap)
    got = ImmutableDenseNDimArray(got).applyfunc(lambda t: t.doit(deep=True))
    if isinstance(expr, MatrixExpr):
        ground = expr.as_explicit().subs(replmap)
    else:
        ground = expr.subs(replmap).doit(deep=True)
    expected = derive_by_array(ground, ImmutableDenseNDimArray(replmap[wrt]))
    got = got.subs(numreplmap)
    expected = expected.subs(numreplmap)
    diff = (got - expected).applyfunc(lambda t: t.doit(deep=True).simplify())
    assert diff == ImmutableDenseNDimArray.zeros(*expected.shape), expr


def test_array_derive_matpow_negative_exponent():
    # The MatPow rule used to build a summation with an empty (hence zero)
    # summation range for negative exponents, so e.g. ``MatPow(X, -2)``
    # silently differentiated to zero:
    M = MatrixSymbol("M", 2, 2)
    Mexp = ImmutableMatrix(2, 2, symbols("m:4"))
    numrepl = dict(zip(symbols("m:4"), [3, 1, 2, 5]))  # invertible matrix
    for e in (-1, -2, -3, 2, 3):
        _check_numeric(MatPow(M, e), M, {M: Mexp}, numrepl)
    # MatPow(M, -1) has to agree with the Inverse rule:
    d1 = array_derive(MatPow(M, -1), M).as_explicit().subs(M, Mexp).subs(numrepl)
    d2 = array_derive(Inverse(M), M).as_explicit().subs(M, Mexp).subs(numrepl)
    diff = (ImmutableDenseNDimArray(d1) - ImmutableDenseNDimArray(d2)
            ).applyfunc(lambda t: t.doit().simplify())
    assert diff == ImmutableDenseNDimArray.zeros(2, 2, 2, 2)
    # products of inverses auto-canonicalize to negative powers:
    _check_numeric((Inverse(M)*Inverse(M)).doit(), M, {M: Mexp}, numrepl)
    # a zeroth power is the identity matrix, a constant:
    assert array_derive(MatPow(M, 0), M) == ZeroArray(2, 2, 2, 2)


def test_array_derive_function_of_matrix_scalar():
    # Scalar functions wrapping matrix scalars (e.g. traces) used to raise
    # NotImplementedError: ``Trace`` is not a ``MatrixExpr``, so the function
    # rule recursed into the generic ``Expr`` fallback.
    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    Mexp = ImmutableMatrix(2, 2, symbols("m:4"))
    Nexp = ImmutableMatrix(2, 2, symbols("n:4"))
    repl = {M: Mexp, N: Nexp}
    numrepl = dict(zip(symbols("m:4"), [3, 1, 2, 5]))
    _check_numeric(Trace(M), M, repl)
    _check_numeric(Trace(Inverse(M)), M, repl, numrepl)
    _check_numeric(sin(Trace(M)), M, repl)
    _check_numeric(exp(Trace(M*N)), M, repl)
    _check_numeric(log(Trace(M.T*M)), M, repl)
    _check_numeric(sin(Trace(M))*cos(Trace(N*M)), M, repl)
    _check_numeric(exp(Trace(M))*Trace(M*N), M, repl)
    # multi-argument functions use the full chain rule:
    _check_numeric(atan2(Trace(M), Trace(M*M)), M, repl)
    # a function of a constant scalar has a zero derivative:
    assert array_derive(sin(Trace(N)), M) == ZeroArray(2, 2)


def test_array_derive_diag_matrix():
    # DiagMatrix used to raise NotImplementedError because it is left
    # unchanged by convert_matrix_to_array and no dedicated rule existed.
    A2 = MatrixSymbol("A", 2, 2)
    B2 = MatrixSymbol("B", 2, 2)
    v = MatrixSymbol("v", 2, 1)
    vexp = ImmutableMatrix(2, 1, symbols("v:2"))
    Aexp = ImmutableMatrix(2, 2, symbols("a:4"))
    Bexp = ImmutableMatrix(2, 2, symbols("b:4"))
    repl = {A2: Aexp, B2: Bexp, v: vexp}
    _check_numeric(DiagMatrix(v), v, repl)
    _check_numeric(DiagMatrix(v.T), v, repl)   # row vector
    _check_numeric(DiagMatrix(A2*v), v, repl)  # composite argument
    _check_numeric(A2*DiagMatrix(v)*B2, v, repl)
    # differentiation with respect to an unrelated matrix gives zero:
    assert array_derive(DiagMatrix(v), A2) == ZeroArray(2, 2, 2, 2)


def test_array_derive_diagonal_matrix():
    A2 = MatrixSymbol("A", 2, 2)
    M = MatrixSymbol("M", 2, 2)
    Mexp = ImmutableMatrix(2, 2, symbols("m:4"))
    Aexp = ImmutableMatrix(2, 2, symbols("a:4"))
    repl = {M: Mexp, A2: Aexp}
    _check_numeric(DiagonalMatrix(M), M, repl)
    _check_numeric(DiagonalMatrix(A2*M), M, repl)
    assert array_derive(DiagonalMatrix(M), A2) == ZeroArray(2, 2, 2, 2)
    # non-square arguments are not supported (elements beyond the diagonal
    # length are zeroed, which the diagonalization cannot express):
    N23 = MatrixSymbol("N", 2, 3)
    raises(NotImplementedError, lambda: array_derive(DiagonalMatrix(N23), N23))


def test_array_derive_ndimarray_returns_result():
    # The NDimArray rule used to be missing its ``return``, so array_derive
    # silently gave ``None`` for any concrete Array.
    x = symbols("x")
    arr = Array([x, x**2, x**3])
    assert array_derive(arr, x) == derive_by_array(arr, x) == Array([1, 2*x, 3*x**2])
