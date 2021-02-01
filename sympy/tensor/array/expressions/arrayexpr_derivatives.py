import operator
from functools import reduce, singledispatch

from sympy import Expr, Transpose, Identity, MatrixSymbol, S
from sympy.codegen.array_utils import CodegenArrayElementwiseAdd, CodegenArrayPermuteDims, CodegenArrayContraction, \
    get_shape, CodegenArrayTensorProduct
from sympy.combinatorics.permutations import _af_invert
from sympy.tensor.array.expressions.array_expressions import ZeroArray


@singledispatch
def array_derive(expr, x):
    raise NotImplementedError(f"not implemented for type {type(expr)}")


@array_derive.register(Expr)
def _(expr: Expr, x: Expr):
    return ZeroArray(*x.shape)


@array_derive.register(CodegenArrayTensorProduct)
def _(expr: CodegenArrayTensorProduct, x: Expr):
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
        addend = CodegenArrayTensorProduct(*args_prev, darg, *args_succ)
        tot1 = len(get_shape(x))
        tot2 = tot1 + len(shape_prev)
        tot3 = tot2 + len(get_shape(arg))
        tot4 = tot3 + len(shape_succ)
        perm = [i for i in range(tot1, tot2)] + \
               [i for i in range(tot1)] + [i for i in range(tot2, tot3)] + \
               [i for i in range(tot3, tot4)]
        addend = CodegenArrayPermuteDims(addend, _af_invert(perm))
        addend_list.append(addend)
    if len(addend_list) == 1:
        return addend_list[0]
    elif len(addend_list) == 0:
        return S.Zero
    else:
        return CodegenArrayElementwiseAdd(*addend_list)


@array_derive.register(MatrixSymbol)
def _(expr: MatrixSymbol, x: Expr):
    m, n = expr.shape
    if expr == x:
        return CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(Identity(m), Identity(n)),
            [0, 2, 1, 3]
        )
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(Identity)
def _(expr: Identity, x: Expr):
    return ZeroArray(*(x.shape + expr.shape))


@array_derive.register(Transpose)
def _(expr: Transpose, x: Expr):
    # D(A.T, A) ==> (m,n,i,j) ==> D(A_ji, A_mn) = d_mj d_ni
    # D(B.T, A) ==> (m,n,i,j) ==> D(B_ji, A_mn)
    fd = array_derive(expr.arg, x)
    return CodegenArrayPermuteDims(fd, [0, 1, 3, 2])


@array_derive.register(CodegenArrayContraction)
def _(expr: CodegenArrayContraction, x: Expr):
    fd = array_derive(expr.expr, x)
    rank_x = len(get_shape(x))
    contraction_indices = expr.contraction_indices
    new_contraction_indices = [tuple(j + rank_x for j in i) for i in contraction_indices]
    return CodegenArrayContraction(fd, *new_contraction_indices)


@array_derive.register(CodegenArrayElementwiseAdd)
def _(expr: CodegenArrayElementwiseAdd, x: Expr):
    return CodegenArrayElementwiseAdd(*[array_derive(arg, x) for arg in expr.args])


@array_derive.register(CodegenArrayPermuteDims)
def _(expr: CodegenArrayPermuteDims, x: Expr):
    de = array_derive(expr.expr, x)
    perm = [0, 1] + [i + 2 for i in expr.permutation.array_form]
    return CodegenArrayPermuteDims(de, perm)


def matrix_derive(expr, x):
    from sympy.codegen.array_utils import parse_matrix_expression, recognize_matrix_expression
    ce = parse_matrix_expression(expr)
    dce = array_derive(ce, x)
    return recognize_matrix_expression(dce).doit()
