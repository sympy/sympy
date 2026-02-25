from __future__ import annotations
from sympy import KroneckerProduct, Add
from sympy.core.function import Lambda
from sympy.core.mul import Mul
from sympy.core.numbers import Integer
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import (Dummy, symbols)
from sympy.matrices.expressions.hadamard import (HadamardPower, HadamardProduct)
from sympy.matrices.expressions.matadd import MatAdd
from sympy.matrices.expressions.matmul import MatMul
from sympy.matrices.expressions.matpow import MatPow
from sympy.matrices.expressions.trace import Trace
from sympy.matrices.expressions.transpose import Transpose
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array.expressions.array_expressions import \
    ArrayElementwiseApplyFunc, _array_tensor_product, _array_contraction, \
    _array_diagonal, _array_add, _permute_dims, Reshape, get_shape, _ArrayExpr, _CodegenArrayAbstract
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.core.basic import Basic


def _array_elementwise_apply_func(function, element):
    if not isinstance(element, (MatrixExpr, _ArrayExpr, _CodegenArrayAbstract)):
        return function(element)
    return ArrayElementwiseApplyFunc(function, element)


def convert_matrix_to_array(expr: Basic) -> Basic:
    if isinstance(expr, MatMul):
        args_nonmat = []
        args = []
        for arg in expr.args:
            if isinstance(arg, MatrixExpr):
                args.append(arg)
            else:
                args_nonmat.append(convert_matrix_to_array(arg))
        contractions = [(2*i+1, 2*i+2) for i in range(len(args)-1)]
        scalar = _array_tensor_product(*args_nonmat) if args_nonmat else S.One
        if scalar == 1:
            tprod = _array_tensor_product(
                *[convert_matrix_to_array(arg) for arg in args])
        else:
            tprod = _array_tensor_product(
                scalar,
                *[convert_matrix_to_array(arg) for arg in args])
        return _array_contraction(
                tprod,
                *contractions
        )
    elif isinstance(expr, (Add, MatAdd)):
        return _array_add(
                *[convert_matrix_to_array(arg) for arg in expr.args]
        )
    elif isinstance(expr, Transpose):
        return _permute_dims(
                convert_matrix_to_array(expr.args[0]), [1, 0]
        )
    elif isinstance(expr, Trace):
        inner_expr: MatrixExpr = convert_matrix_to_array(expr.arg) # type: ignore
        return _array_contraction(inner_expr, (0, len(inner_expr.shape) - 1))
    elif isinstance(expr, Mul):
        return _array_tensor_product(*[convert_matrix_to_array(i) for i in expr.args])
    elif isinstance(expr, Pow):
        base = convert_matrix_to_array(expr.base)
        if (expr.exp > 0) == True:
            base_conv = convert_matrix_to_array(base)
            if get_shape(base_conv) == ():
                return _array_elementwise_apply_func(lambda x: x**expr.exp, base_conv)
            else:
                return _array_tensor_product(*[base_conv for i in range(expr.exp)])
        elif get_shape(base) == ():
            d = Dummy("d")
            return _array_elementwise_apply_func(Lambda(d, d**expr.exp), base)
        else:
            return expr
    elif isinstance(expr, MatPow):
        base = convert_matrix_to_array(expr.base)
        if expr.exp.is_Integer != True:
            if get_shape(expr) == (1, 1):
                b = symbols("b", cls=Dummy)
                return _array_elementwise_apply_func(Lambda(b, b ** expr.exp), convert_matrix_to_array(base))
            return expr
        elif (expr.exp > 0) == True:
            return convert_matrix_to_array(MatMul.fromiter(base for i in range(expr.exp)))
        else:
            return expr
    elif isinstance(expr, HadamardProduct):
        tp = _array_tensor_product(*[convert_matrix_to_array(arg) for arg in expr.args])
        diag = [[2*i for i in range(len(expr.args))], [2*i+1 for i in range(len(expr.args))]]
        return _array_diagonal(tp, *diag)
    elif isinstance(expr, HadamardPower):
        base, exp = expr.args
        if isinstance(exp, Integer) and exp > 0:
            return convert_matrix_to_array(HadamardProduct.fromiter(base for i in range(exp)))
        else:
            d = Dummy("d")
            return _array_elementwise_apply_func(Lambda(d, d**exp), base)
    elif isinstance(expr, KroneckerProduct):
        kp_args = [convert_matrix_to_array(arg) for arg in expr.args]
        permutation = [2*i for i in range(len(kp_args))] + [2*i + 1 for i in range(len(kp_args))]
        return Reshape(_permute_dims(_array_tensor_product(*kp_args), permutation), expr.shape)
    else:
        return expr
