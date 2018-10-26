from collections import defaultdict
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.mul import Mul
from sympy.core.compatibility import accumulate, default_sort_key
from sympy.matrices.expressions import MatMul
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array import NDimArray


class CodegenArrayContraction(Basic):
    """
    This class is meant to represent contractions of arrays in a form easily
    processable by the code printers.
    """
    def __new__(cls, expr, *contraction_indices):
        contraction_indices = [_sympify(ind) for ind in contraction_indices]
        expr = _sympify(expr)
        obj = Basic.__new__(cls, expr, *contraction_indices)
        obj._mapping = cls._get_mapping_from_contraction_indices(expr, *contraction_indices)
        return obj

    @staticmethod
    def _get_mapping_from_contraction_indices(expr, *contraction_indices):
        if not isinstance(expr, CodegenArrayTensorProduct):
            return {}
        args = expr.args
        ranks = expr.ranks
        mapping = {}
        counter = 0
        for i, rank in enumerate(ranks):
            for j in range(rank):
                mapping[counter] = (i, j)
                counter += 1
        return mapping

    def _contraction_indices_to_contraction_tuples(self):
        mapping = self._mapping
        if not mapping:
            raise NotImplementedError
        return [[mapping[j] for j in i] for i in self.contraction_indices]

    @staticmethod
    def _contraction_tuple_to_contraction_indices(expr, contraction_tuples):
        # TODO: check that `expr` has `.ranks`:
        ranks = expr.ranks
        cumulative_ranks = [0] + list(accumulate(ranks))
        return [tuple(cumulative_ranks[j]+k for j, k in i) for i in contraction_tuples]

    @property
    def expr(self):
        return self.args[0]

    @property
    def contraction_indices(self):
        return self.args[1:]

    def _contraction_indices_to_components(self):
        expr = self.expr
        if not isinstance(expr, CodegenArrayTensorProduct):
            raise NotImplementedError("only for contractions of tensor products")
        contraction_indices = self.contraction_indices
        args = expr.args
        ranks = expr.ranks
        mapping = {}
        counter = 0
        for i, rank in enumerate(ranks):
            for j in range(rank):
                mapping[counter] = (i, j)
                counter += 1
        return mapping

    def sort_args_by_name(self):
        expr = self.expr
        if not isinstance(expr, CodegenArrayTensorProduct):
            return self
        args = expr.args
        sorted_data = sorted(enumerate(args), key=lambda x: default_sort_key(x[1]))
        pos_sorted, args_sorted = zip(*sorted_data)
        reordering_map = {i: pos_sorted.index(i) for i, arg in enumerate(args)}
        contraction_tuples = self._contraction_indices_to_contraction_tuples()
        contraction_tuples = [[(reordering_map[j], k) for j, k in i] for i in contraction_tuples]
        c_tp = CodegenArrayTensorProduct(*args_sorted)
        new_contr_indices = self._contraction_tuple_to_contraction_indices(
                c_tp,
                contraction_tuples
        )
        return CodegenArrayContraction(c_tp, *new_contr_indices)

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
