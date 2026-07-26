from __future__ import annotations
import operator
from functools import reduce, singledispatch

from sympy.core.singleton import S
from sympy import MatrixBase, derive_by_array, Integer, Determinant, Function, MatPow, Dummy, Pow, Mul
from sympy.tensor.array import NDimArray
from sympy.core.expr import Expr
from sympy.matrices.expressions.hadamard import HadamardProduct
from sympy.matrices.expressions.inverse import Inverse
from sympy.matrices.expressions.matexpr import (MatrixExpr, MatrixSymbol, MatrixElement)
from sympy.matrices.expressions.special import Identity, OneMatrix, MatrixUnit, ZeroMatrix
from sympy.matrices.expressions.transpose import Transpose
from sympy.combinatorics.permutations import _af_invert
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.tensor.array.expressions.array_expressions import (
    _ArrayExpr, ZeroArray, ArraySymbol, ArrayTensorProduct, ArrayAdd,
    PermuteDims, ArrayDiagonal, ArrayElementwiseApplyFunc, get_ndim,
    get_shape, ArrayContraction, _array_tensor_product, _array_contraction,
    _array_diagonal, _array_add, _permute_dims, Reshape, ArraySum)
from sympy.tensor.array.expressions.from_matrix_to_array import convert_matrix_to_array


@singledispatch
def array_derive(expr, x):
    """
    Derivatives (gradients) for array expressions.
    """
    raise NotImplementedError(f"not implemented for type {type(expr)}")


@array_derive.register(Expr)
def _(expr: Expr, x: _ArrayExpr):
    if expr.free_symbols & x.free_symbols:
        if isinstance(expr, MatrixElement):
            if expr.parent == x and isinstance(x, MatrixSymbol):
                return MatrixUnit(x.shape[0], x.shape[1], expr.i, expr.j)
            # Differentiate the parent matrix and extract the (i, j) element
            # by contracting the two matrix axes with a matrix unit:
            dparent = array_derive(expr.parent, x)
            rank_x = len(get_shape(x))
            m, n = expr.parent.shape
            tp = _array_tensor_product(dparent, MatrixUnit(m, n, expr.i, expr.j))
            return _array_contraction(tp, (rank_x, rank_x + 2), (rank_x + 1, rank_x + 3))
        raise NotImplementedError("algorithm not implemented for this case")
    return ZeroArray(*x.shape)


@array_derive.register(Mul)
def _(expr: Mul, x: _ArrayExpr):
    args = expr.args
    return ArrayAdd.fromiter([
        _array_tensor_product(Mul.fromiter(args[:i]), array_derive(arg, x), Mul.fromiter(args[(i+1):]))
        for i, arg in enumerate(args)
    ])


@array_derive.register(Pow)
def _(expr: Pow, x: _ArrayExpr):
    return Pow._eval_derivative(expr, x)


@array_derive.register(Function)
def _(expr: Function, x: _ArrayExpr):
    if len(expr.args) != 1:
        raise NotImplementedError("only 1-parameter functions are supported")
    dexpr = array_derive(expr.args[0], x)
    return _array_tensor_product(expr.fdiff(), dexpr)


@array_derive.register(ArrayTensorProduct)
def _(expr: ArrayTensorProduct, x: Expr):
    args = expr.args
    addend_list = []
    for i, arg in enumerate(expr.args):
        darg = array_derive(arg, x)
        if darg == 0:
            continue
        args_prev = args[:i]
        args_succ = args[i+1:]
        shape_prev = reduce(operator.add, map(get_shape, args_prev), ())
        shape_succ = reduce(operator.add, map(get_shape, args_succ), ())
        addend = _array_tensor_product(*args_prev, darg, *args_succ)
        tot1 = len(get_shape(x))
        tot2 = tot1 + len(shape_prev)
        tot3 = tot2 + len(get_shape(arg))
        tot4 = tot3 + len(shape_succ)
        perm = list(range(tot1, tot2)) + \
               list(range(tot1)) + list(range(tot2, tot3)) + \
               list(range(tot3, tot4))
        addend = _permute_dims(addend, _af_invert(perm))
        addend_list.append(addend)
    if len(addend_list) == 1:
        return addend_list[0]
    elif len(addend_list) == 0:
        return S.Zero
    else:
        return _array_add(*addend_list)


@array_derive.register(ArraySymbol)
def _(expr: ArraySymbol, x: _ArrayExpr):
    if expr == x:
        return _permute_dims(
            ArrayTensorProduct.fromiter(Identity(i) for i in expr.shape),
            [2*i for i in range(len(expr.shape))] + [2*i+1 for i in range(len(expr.shape))]
        )
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(MatrixSymbol)
def _(expr: MatrixSymbol, x: _ArrayExpr):
    m, n = expr.shape
    if expr == x:
        return _permute_dims(
            _array_tensor_product(Identity(m), Identity(n)),
            [0, 2, 1, 3]
        )
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(Determinant)
def _(expr: Determinant, x: Expr):
    arg = expr.arg
    arg_inverse = arg.inv()
    darg = array_derive(arg, x)
    rank_x = len(get_shape(x))
    tp = _array_tensor_product(expr, arg_inverse, darg)
    # axes: (0, 1) = arg_inverse, (2, ..., rank_x + 1) = x,
    # (rank_x + 2, rank_x + 3) = matrix axes of darg
    tc = _array_contraction(tp, (0, rank_x + 3), (1, rank_x + 2))
    return tc


@array_derive.register(Identity)
def _(expr: Identity, x: _ArrayExpr):
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(OneMatrix)
def _(expr: OneMatrix, x: _ArrayExpr):
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(ZeroMatrix)
def _(expr: ZeroMatrix, x: _ArrayExpr):
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(ZeroArray)
def _(expr: ZeroArray, x: _ArrayExpr):
    return ZeroArray(*(get_shape(x) + expr.shape))


@array_derive.register(Transpose)
def _(expr: Transpose, x: Expr):
    # D(A.T, A) ==> (m,n,i,j) ==> D(A_ji, A_mn) = d_mj d_ni
    # D(B.T, A) ==> (m,n,i,j) ==> D(B_ji, A_mn)
    fd = array_derive(expr.arg, x)
    # swap the two matrix axes that follow the ``rank_x`` axes of ``x``:
    rank_x = len(get_shape(x))
    return _permute_dims(fd, list(range(rank_x)) + [rank_x + 1, rank_x])


@array_derive.register(Inverse)
def _(expr: Inverse, x: Expr):
    mat = expr.I
    dexpr = array_derive(mat, x)
    rank_x = len(get_shape(x))
    tp = _array_tensor_product(-expr, dexpr, expr)
    # axes: (0, 1) = -expr, (2, ..., rank_x + 1) = x,
    # (rank_x + 2, rank_x + 3) = matrix axes of dexpr, last two = expr
    mp = _array_contraction(tp, (1, rank_x + 2), (rank_x + 3, rank_x + 4))
    pp = _permute_dims(mp, list(range(1, rank_x + 1)) + [0, rank_x + 1])
    return pp


@array_derive.register(ElementwiseApplyFunction)
def _(expr: ElementwiseApplyFunction, x: Expr):
    assert get_ndim(expr) == 2
    assert get_ndim(x) == 2
    fdiff = expr._get_function_fdiff()
    dexpr = array_derive(expr.expr, x)
    tp = _array_tensor_product(
        ElementwiseApplyFunction(fdiff, expr.expr),
        dexpr
    )
    td = _array_diagonal(
        tp, (0, 4), (1, 5)
    )
    return td


@array_derive.register(ArrayElementwiseApplyFunc)
def _(expr: ArrayElementwiseApplyFunc, x: Expr):
    fdiff = expr._get_function_fdiff()
    subexpr = expr.expr
    dsubexpr = array_derive(subexpr, x)
    tp = _array_tensor_product(
        dsubexpr,
        ArrayElementwiseApplyFunc(fdiff, subexpr)
    )
    b = get_ndim(x)
    c = get_ndim(expr)
    diag_indices = [(b + i, b + c + i) for i in range(c)]
    return _array_diagonal(tp, *diag_indices)


@array_derive.register(MatrixExpr)
def _(expr: MatrixExpr, x: Expr):
    cg = convert_matrix_to_array(expr)
    if cg == expr:
        # Avoid infinite looping:
        raise NotImplementedError()
    return array_derive(cg, x)


@array_derive.register(MatPow)
def _(expr: MatPow, x: Expr):
    base = expr.base
    exponent = expr.exp
    dbase = array_derive(base, x)
    dexponent = array_derive(exponent, x)
    if not isinstance(dexponent, ZeroArray) or (dexponent == 0) == True:
        raise NotImplementedError()
    d = Dummy("d")
    rank_x = len(get_shape(x))
    tp = _array_tensor_product(base**d, dbase, base**(exponent-d-1))
    tc = _array_contraction(tp, (1, rank_x + 2), (rank_x + 3, rank_x + 4))
    pd = _permute_dims(tc, list(range(1, rank_x + 1)) + [0, rank_x + 1])
    if isinstance(pd, ZeroArray):
        # a sum of zero arrays is zero:
        return pd
    return ArraySum(pd,(d, 0, exponent-1))


@array_derive.register(HadamardProduct)
def _(expr: HadamardProduct, x: Expr):
    raise NotImplementedError()


@array_derive.register(ArrayContraction)
def _(expr: ArrayContraction, x: Expr):
    fd = array_derive(expr.expr, x)
    rank_x = len(get_shape(x))
    contraction_indices = expr.contraction_indices
    new_contraction_indices = [tuple(j + rank_x for j in i) for i in contraction_indices]
    return _array_contraction(fd, *new_contraction_indices)


@array_derive.register(ArrayDiagonal)
def _(expr: ArrayDiagonal, x: Expr):
    dsubexpr = array_derive(expr.expr, x)
    rank_x = len(get_shape(x))
    diag_indices = [[j + rank_x for j in i] for i in expr.diagonal_indices]
    return _array_diagonal(dsubexpr, *diag_indices)


@array_derive.register(ArrayAdd)
def _(expr: ArrayAdd, x: Expr):
    return _array_add(*[array_derive(arg, x) for arg in expr.args])


@array_derive.register(PermuteDims)
def _(expr: PermuteDims, x: Expr):
    de = array_derive(expr.expr, x)
    # Differentiating prepends the ``rank_x`` axes of ``x`` to the result, so
    # those must be left in place and the original permutation shifted past
    # them.  (This used to assume ``x`` was always a matrix, i.e. ``rank_x == 2``.)
    rank_x = len(get_shape(x))
    perm = list(range(rank_x)) + [i + rank_x for i in expr.permutation.array_form]
    return _permute_dims(de, perm)


@array_derive.register(Reshape)
def _(expr: Reshape, x: Expr):
    de = array_derive(expr.expr, x)
    return Reshape(de, get_shape(x) + expr.shape)


@array_derive.register(MatrixBase)
def _(expr: MatrixBase, x):
    if not set.intersection(expr.free_symbols, x.free_symbols):
        return ZeroArray(*x.shape, *expr.shape)
    if isinstance(x, MatrixExpr) and all(isinstance(i, (int, Integer)) for i in x.shape):
        x = x.as_explicit()
    if isinstance(x, MatrixBase):
        return derive_by_array(expr, x)
    raise NotImplementedError("could not determine derivative")


@array_derive.register(NDimArray)
def _(expr: NDimArray, x):
    if not set.intersection(expr.free_symbols, x.free_symbols):
        return ZeroArray(*get_shape(x), *expr.shape)
    if isinstance(x, (MatrixExpr, ArraySymbol)) and all(isinstance(i, (int, Integer)) for i in get_shape(x)):
        x = x.as_explicit()
    return derive_by_array(expr, x)


def matrix_derive(expr, x):
    from sympy.tensor.array.expressions.from_array_to_matrix import convert_array_to_matrix
    ce = convert_matrix_to_array(expr)
    dce = array_derive(ce, x)
    return convert_array_to_matrix(dce).doit()
