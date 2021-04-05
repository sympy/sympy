from sympy import Mul, Basic, MatMul, MatAdd, Transpose, Trace, Pow, \
    MatPow, symbols, Dummy, Lambda, HadamardProduct, HadamardPower, S
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array.expressions.array_expressions import ArrayDiagonal, ArrayTensorProduct, \
    PermuteDims, ArrayAdd, ArrayContraction, ArrayElementwiseApplyFunc


def convert_matrix_to_array(expr: MatrixExpr) -> Basic:
    if isinstance(expr, MatMul):
        args_nonmat = []
        args = []
        for arg in expr.args:
            if isinstance(arg, MatrixExpr):
                args.append(arg)
            else:
                args_nonmat.append(convert_matrix_to_array(arg))
        contractions = [(2*i+1, 2*i+2) for i in range(len(args)-1)]
        scalar = ArrayTensorProduct.fromiter(args_nonmat) if args_nonmat else S.One
        if scalar == 1:
            tprod = ArrayTensorProduct(
                *[convert_matrix_to_array(arg) for arg in args])
        else:
            tprod = ArrayTensorProduct(
                scalar,
                *[convert_matrix_to_array(arg) for arg in args])
        return ArrayContraction(
                tprod,
                *contractions
        )
    elif isinstance(expr, MatAdd):
        return ArrayAdd(
                *[convert_matrix_to_array(arg) for arg in expr.args]
        )
    elif isinstance(expr, Transpose):
        return PermuteDims(
                convert_matrix_to_array(expr.args[0]), [1, 0]
        )
    elif isinstance(expr, Trace):
        inner_expr = convert_matrix_to_array(expr.arg)
        return ArrayContraction(inner_expr, (0, len(inner_expr.shape) - 1))
    elif isinstance(expr, Mul):
        return ArrayTensorProduct.fromiter(convert_matrix_to_array(i) for i in expr.args)
    elif isinstance(expr, Pow):
        base = convert_matrix_to_array(expr.base)
        if (expr.exp > 0) == True:
            return ArrayTensorProduct.fromiter(base for i in range(expr.exp))
        else:
            return expr
    elif isinstance(expr, MatPow):
        base = convert_matrix_to_array(expr.base)
        if expr.exp.is_Integer != True:
            b = symbols("b", cls=Dummy)
            return ArrayElementwiseApplyFunc(Lambda(b, b**expr.exp), convert_matrix_to_array(base))
        elif (expr.exp > 0) == True:
            return convert_matrix_to_array(MatMul.fromiter(base for i in range(expr.exp)))
        else:
            return expr
    elif isinstance(expr, HadamardProduct):
        tp = ArrayTensorProduct.fromiter(expr.args)
        diag = [[2*i for i in range(len(expr.args))], [2*i+1 for i in range(len(expr.args))]]
        return ArrayDiagonal(tp, *diag)
    elif isinstance(expr, HadamardPower):
        base, exp = expr.args
        return convert_matrix_to_array(HadamardProduct.fromiter(base for i in range(exp)))
    else:
        return expr
