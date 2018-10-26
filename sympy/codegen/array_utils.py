from collections import defaultdict
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.mul import Mul
from sympy.matrices.expressions import MatMul
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array import NDimArray


class CodegenArrayContraction(Basic):
    """
    This class is meant to represent contractions of arrays in a form easily
    processable by the code printers.
    """
    def __new__(cls, expr, *indices):
        indices = [_sympify(ind) for ind in indices]
        expr = _sympify(expr)
        return Basic.__new__(cls, expr, *indices)

    @property
    def expr(self):
        return self.args[0]

    @property
    def contraction_indices(self):
        return self.args[1:]

    @staticmethod
    def from_summation(summation):
        from sympy import Indexed
        from sympy.matrices.expressions.matexpr import MatrixElement
        function = summation.function
        indices = summation.variables
        if function.is_Mul:
            pass
            args = []
            ranks = []
            total_rank = 0
            axes_contraction = defaultdict(list)
            for arg in function.args:
                loc_indices = None
                if isinstance(arg, Indexed):
                    args.append(arg.base)
                    loc_indices = arg.indices
                    ranks.append(len(arg.indices))
                elif isinstance(arg, MatrixElement):
                    args.append(arg.parent)
                    loc_indices = arg.indices
                    ranks.append(2)
                elif isinstance(arg, Matrix):
                    args.append(arg)
                    ranks.append(2)
                else:
                    args.append(arg)
                    ranks.append(0)
                for i, ind in enumerate(loc_indices):
                    if ind in indices:
                        axes_contraction[ind].append(total_rank + i)
                total_rank += ranks[-1]
            return CodegenArrayContraction(
                    CodegenArrayTensorProduct(*args),
                    *[tuple(v) for v in axes_contraction.values()]
                )

    @staticmethod
    def from_MatMul(expr):
        args_nonmat = []
        args = []
        contractions = []
        for arg in expr.args:
            if isinstance(arg, MatrixExpr):
                args.append(arg)
            else:
                args_nonmat.append(arg)
        contractions = [(2*i+1, 2*i+2) for i in range(len(args)-1)]
        return Mul.fromiter(args_nonmat)*CodegenArrayContraction(
                CodegenArrayTensorProduct(*args),
                *contractions
            )


class CodegenArrayTensorProduct(Basic):
    def __new__(cls, *args):
        args = [_sympify(arg) for arg in args]
        ranks = [cls.get_rank(arg) for arg in args]
        if len(args) == 1:
            return args[0]
        obj = Basic.__new__(cls, *args)
        obj._ranks = ranks
        return obj

    @property
    def ranks(self):
        return self._ranks[:]

    @classmethod
    def get_rank(cls, expr):
        if isinstance(expr, MatrixExpr):
            return 2
        if isinstance(expr, NDimArray):
            return expr.rank()
        return 0
