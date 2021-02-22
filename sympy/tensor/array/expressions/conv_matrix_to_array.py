from sympy import Mul, Basic, MatMul, MatAdd, Transpose, Trace, Pow, \
    MatPow, symbols, Dummy, Lambda, HadamardProduct, HadamardPower, S
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array.expressions.array_expressions import CodegenArrayDiagonal, CodegenArrayTensorProduct, \
    CodegenArrayPermuteDims, CodegenArrayElementwiseAdd, CodegenArrayContraction, ArrayElementwiseApplyFunc


def parse_matrix_expression(expr: MatrixExpr) -> Basic:
    if isinstance(expr, MatMul):
        args_nonmat = []
        args = []
        for arg in expr.args:
            if isinstance(arg, MatrixExpr):
                args.append(arg)
            else:
                args_nonmat.append(parse_matrix_expression(arg))
        contractions = [(2*i+1, 2*i+2) for i in range(len(args)-1)]
        scalar = CodegenArrayTensorProduct.fromiter(args_nonmat) if args_nonmat else S.One
        if scalar == 1:
            tprod = CodegenArrayTensorProduct(
                *[parse_matrix_expression(arg) for arg in args])
        else:
            tprod = CodegenArrayTensorProduct(
                scalar,
                *[parse_matrix_expression(arg) for arg in args])
        return CodegenArrayContraction(
                tprod,
                *contractions
        )
    elif isinstance(expr, MatAdd):
        return CodegenArrayElementwiseAdd(
                *[parse_matrix_expression(arg) for arg in expr.args]
        )
    elif isinstance(expr, Transpose):
        return CodegenArrayPermuteDims(
                parse_matrix_expression(expr.args[0]), [1, 0]
        )
    elif isinstance(expr, Trace):
        inner_expr = parse_matrix_expression(expr.arg)
        return CodegenArrayContraction(inner_expr, (0, len(inner_expr.shape) - 1))
    elif isinstance(expr, Mul):
        return CodegenArrayTensorProduct.fromiter(parse_matrix_expression(i) for i in expr.args)
    elif isinstance(expr, Pow):
        base = parse_matrix_expression(expr.base)
        if (expr.exp > 0) == True:
            return CodegenArrayTensorProduct.fromiter(base for i in range(expr.exp))
        else:
            return expr
    elif isinstance(expr, MatPow):
        base = parse_matrix_expression(expr.base)
        if expr.exp.is_Integer != True:
            b = symbols("b", cls=Dummy)
            return ArrayElementwiseApplyFunc(Lambda(b, b**expr.exp), parse_matrix_expression(base))
        elif (expr.exp > 0) == True:
            return parse_matrix_expression(MatMul.fromiter(base for i in range(expr.exp)))
        else:
            return expr
    elif isinstance(expr, HadamardProduct):
        tp = CodegenArrayTensorProduct.fromiter(expr.args)
        diag = [[2*i for i in range(len(expr.args))], [2*i+1 for i in range(len(expr.args))]]
        return CodegenArrayDiagonal(tp, *diag)
    elif isinstance(expr, HadamardPower):
        base, exp = expr.args
        return parse_matrix_expression(HadamardProduct.fromiter(base for i in range(exp)))
    else:
        return expr
