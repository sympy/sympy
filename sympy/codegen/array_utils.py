from collections import defaultdict
from sympy import Indexed, IndexedBase, Tuple
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.mul import Mul
from sympy.core.compatibility import accumulate, default_sort_key
from sympy.matrices.expressions import MatMul, Trace, Transpose
from sympy.matrices.expressions.matexpr import MatrixExpr, MatrixElement
from sympy.tensor.array import NDimArray


class CodegenArrayContraction(Basic):
    """
    This class is meant to represent contractions of arrays in a form easily
    processable by the code printers.
    """
    def __new__(cls, expr, *contraction_indices, **kwargs):
        contraction_indices = cls._sort_contraction_indices(contraction_indices)
        expr = _sympify(expr)
        obj = Basic.__new__(cls, expr, *contraction_indices)
        obj._mapping, obj._ranks = cls._get_mapping_from_contraction_indices(expr, *contraction_indices)
        free_indices = kwargs.get("free_indices", None)
        if "free_indices" in kwargs:
            free_indices = kwargs["free_indices"]
        else:
            free_indices = {i: i for i in range(sum(obj._ranks)) if all([i not in cind for cind in contraction_indices])}
        obj._free_indices = free_indices
        return obj

    @staticmethod
    def _sort_contraction_indices(contraction_indices):
        contraction_indices = [Tuple(*sorted(i)) for i in contraction_indices]
        contraction_indices.sort(key=lambda x: min(x))
        return contraction_indices

    @staticmethod
    def _get_mapping_from_contraction_indices(expr, *contraction_indices):
        if isinstance(expr, CodegenArrayTensorProduct):
            args = expr.args
            ranks = expr.ranks
        else:
            args = [expr]
            ranks = [CodegenArrayTensorProduct.get_rank(expr)]
        mapping = {}
        counter = 0
        for i, rank in enumerate(ranks):
            for j in range(rank):
                mapping[counter] = (i, j)
                counter += 1
        return mapping, ranks

    def _get_contraction_tuples(self):
        r"""
        Return tuples containing the argument index and position within the
        argument of the index position.

        Examples
        ========

        >>> from sympy import MatrixSymbol, MatrixExpr, Sum, Symbol
        >>> from sympy.abc import i, j, k, l, N
        >>> from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)

        >>> cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (1, 2))
        >>> cg._get_contraction_tuples()
        [[(0, 1), (1, 0)]]

        Here the contraction pair `(1, 2)` meaning that the 2nd and 3rd indices
        of the tensor product `A\otimes B` are contracted, has been transformed
        into `(0, 1)` and `(1, 0)`, identifying the same indices in a different
        notation. `(0, 1)` is the second index (1) of the first argument (i.e.
                0 or `A`). `(1, 0)` is the first index (i.e. 0) of the second
        argument (i.e. 1 or `B`).
        """
        mapping = self._mapping
        #if mapping is None:
            #raise NotImplementedError
        return [[mapping[j] for j in i] for i in self.contraction_indices]

    @staticmethod
    def _contraction_tuples_to_contraction_indices(expr, contraction_tuples):
        # TODO: check that `expr` has `.ranks`:
        ranks = expr.ranks
        cumulative_ranks = [0] + list(accumulate(ranks))
        return [tuple(cumulative_ranks[j]+k for j, k in i) for i in contraction_tuples]

    @property
    def free_indices(self):
        return dict(self._free_indices)

    @property
    def expr(self):
        return self.args[0]

    @property
    def ranks(self):
        return self._ranks[:]

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
        """
        Sort arguments in the tensor product so that their order is lexicographical.

        Examples
        ========

        >>> from sympy import MatrixSymbol, MatrixExpr, Sum, Symbol
        >>> from sympy.abc import i, j, k, l, N
        >>> from sympy.codegen.array_utils import CodegenArrayContraction
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)
        >>> C = MatrixSymbol("C", N, N)
        >>> D = MatrixSymbol("D", N, N)

        >>> cg = CodegenArrayContraction.from_MatMul(C*D*A*B)
        >>> cg
        CodegenArrayContraction(CodegenArrayTensorProduct(C, D, A, B), (1, 2), (3, 4), (5, 6))
        >>> cg.sort_args_by_name()
        CodegenArrayContraction(CodegenArrayTensorProduct(A, B, C, D), (0, 7), (1, 2), (5, 6))
        """
        expr = self.expr
        if not isinstance(expr, CodegenArrayTensorProduct):
            return self
        args = expr.args
        sorted_data = sorted(enumerate(args), key=lambda x: default_sort_key(x[1]))
        pos_sorted, args_sorted = zip(*sorted_data)
        reordering_map = {i: pos_sorted.index(i) for i, arg in enumerate(args)}
        contraction_tuples = self._get_contraction_tuples()
        contraction_tuples = [[(reordering_map[j], k) for j, k in i] for i in contraction_tuples]
        c_tp = CodegenArrayTensorProduct(*args_sorted)
        new_contr_indices = self._contraction_tuples_to_contraction_indices(
                c_tp,
                contraction_tuples
        )
        return CodegenArrayContraction(c_tp, *new_contr_indices)

    def _get_contraction_links(self):
        r"""
        Returns a dictionary of links between arguments in the tensor product
        being contracted.

        See the example for an explanation of the values.

        Examples
        ========

        >>> from sympy import MatrixSymbol, MatrixExpr, Sum, Symbol
        >>> from sympy.abc import i, j, k, l, N
        >>> from sympy.codegen.array_utils import CodegenArrayContraction
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)
        >>> C = MatrixSymbol("C", N, N)
        >>> D = MatrixSymbol("D", N, N)

        Matrix multiplications are pairwise contractions between neighboring
        matrices:

        `A_{ij} B_{jk} C_{kl} D_{lm}`

        >>> cg = CodegenArrayContraction.from_MatMul(A*B*C*D)
        >>> cg
        CodegenArrayContraction(CodegenArrayTensorProduct(A, B, C, D), (1, 2), (3, 4), (5, 6))
        >>> cg._get_contraction_links()
        {0: [(1, 0)], 1: [(0, 1), (2, 0)], 2: [(1, 1), (3, 0)], 3: [(2, 1)]}

        This dictionary is interpreted as follows: argument in position 0 (i.e.
        matrix `A`) is contracted to `(1, 0)`, that is argument in position 1
        (matrix `B`) on the first index slot of `B`, this is the contraction
        provided by the index `j` from `A`.

        The argument in position 1 (that is, matrix `B`) has two contractions,
        the ones provided by the indices `j` and `k`.  The link `(0, 1)` and
        `(2, 0)` respectively. `(0, 1)` is the index slot 1 (the 2nd) of
        argument in position 0 (that is, `A_{\ldot j}`), and so on.
        """
        contraction_tuples = self._get_contraction_tuples()
        dlinks = defaultdict(list)
        for links in contraction_tuples:
            if len(links) > 2:
                raise NotImplementedError("three or more axes contracted at the same time")
            (arg1, pos1), (arg2, pos2) = links
            dlinks[arg1].append((arg2, pos2))
            dlinks[arg2].append((arg1, pos1))
        return dict(dlinks)

    def _recognize_matrix_mul_lines(self, first_indices=None):
        r"""
        Recognize lines of matrix multiplications in the contractions of tensor
        products of two-dimensional array.  If the ``CodegenArrayContraction``
        object was created from a summation of indexed expressions, it will
        remember the starting and ending free indices.

        This can help perform the transformation expressed in mathematical
        notation as:

        `\sum_{j=0}^{N-1} A_{i,j} B_{j,k} \Longrightarrow \mathbf{A}\cdot \mathbf{B}`

        Optional parameter ``first_indices``: specify a list of free indices to
        use as the indices of the starting mat-mul lines in the expression.

        Examples
        ========

        >>> from sympy import MatrixSymbol, MatrixExpr, Sum, Symbol
        >>> from sympy.abc import i, j, k, l, N
        >>> from sympy.codegen.array_utils import CodegenArrayContraction
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)

        >>> expr = Sum(A[i, j]*B[j, k], (j, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg.free_indices
        {i: 0, k: 3}
        >>> cg._recognize_matrix_mul_lines()
        [(i, k, [A, B])]
        >>> cg._recognize_matrix_mul_lines(first_indices=[k])
        [(k, i, [B.T, A.T])]

        Transposition is detected:

        >>> expr = Sum(A[j, i]*B[j, k], (j, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg.free_indices
        {i: 1, k: 3}
        >>> cg._recognize_matrix_mul_lines()
        [(i, k, [A.T, B])]
        >>> cg._recognize_matrix_mul_lines(first_indices=[k])
        [(k, i, [B.T, A])]

        Detect the trace:

        >>> expr = Sum(A[i, i], (i, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg._recognize_matrix_mul_lines()
        [(None, None, [Trace(A)])]

        More complicated expressions:

        >>> expr = Sum(A[i, j]*B[k, j]*A[l, k], (j, 0, N-1), (k, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg._recognize_matrix_mul_lines()
        [(i, l, [A, B.T, A.T])]

        Expressions constructed from matrix expressions do not contain literal indices,
        the positions of free indices are returned instead:

        >>> expr = A*B
        >>> cg = CodegenArrayContraction.from_MatMul(expr)
        >>> cg.free_indices
        {0: 0, 3: 3}
        >>> cg._recognize_matrix_mul_lines()
        [(0, 3, [A, B])]
        """
        expr = self.expr
        if not isinstance(expr, CodegenArrayTensorProduct):
            args = [expr]
        else:
            args = expr.args
        free_indices = self.free_indices
        dlinks = self._get_contraction_links()
        # TODO: check that all elements have rank-2
        if free_indices:
            if first_indices:
                first_index = first_indices[0]
                starting_argind = free_indices.pop(first_index)
            else:
                first_index, starting_argind = min(free_indices.items(), key=lambda x: x[1])
                free_indices.pop(first_index)
            current_argind = starting_argind
        else:
            first_index = None
            starting_argind = 0
            current_argind = 0
        current_argind, current_pos = self._mapping[current_argind]
        starting_argind, starting_pos = self._mapping[starting_argind]
        matmul_args = []
        prev_argind = None
        prev_pos = 0
        last_index = None
        while True:
            elem = args[current_argind]
            if current_pos != prev_pos:
                elem = Transpose(elem)
            matmul_args.append(elem)
            if current_argind not in dlinks:
                break
            link_list = [(i, p) for i, p in dlinks.pop(current_argind) if i != prev_argind]
            if len(link_list) == 0:
                if free_indices:
                    last_index = [i for i, j in free_indices.items() if self._mapping[j] == (current_argind, 1-current_pos)][0]
                else:
                    last_index = None
                break
            if len(link_list) != 1:
                if len(link_list) != 2 and link_list[0][0] != link_list[1][0]:
                    raise NotImplementedError("not a matrix multiplication line")
            prev_argind = current_argind
            current_argind, current_pos = link_list[0]
            if current_argind == starting_argind:
                # This is a trace:
                matmul_args = [Trace(*matmul_args)]
                break
        return [(first_index, last_index, matmul_args)]

    @staticmethod
    def from_summation(summation):
        from sympy import Indexed
        from sympy.matrices.expressions.matexpr import MatrixElement
        function = summation.function
        indices = summation.variables
        free_indices = {}
        if function.is_Mul:
            function_args = function.args
        elif function.is_Add:
            pass
        elif isinstance(function, MatrixElement):
            function_args = [function]
        elif isinstance(function, Indexed):
            function_args = [function]
        else:
            return summation
        args = []
        ranks = []
        total_rank = 0
        axes_contraction = defaultdict(list)
        for arg in function_args:
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
            for i, ind in enumerate(loc_indices):
                if ind not in indices:
                    free_indices[ind] = total_rank + i
            total_rank += ranks[-1]
        return CodegenArrayContraction(
                CodegenArrayTensorProduct(*args),
                *[tuple(v) for v in axes_contraction.values()],
                free_indices=free_indices
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
        if isinstance(expr, (MatrixExpr, MatrixElement)):
            return 2
        if isinstance(expr, NDimArray):
            return expr.rank()
        if isinstance(expr, Indexed):
            return expr.rank
        if isinstance(expr, IndexedBase):
            shape = expr.shape
            if shape is None:
                return -1
            else:
                return len(shape)
        return 0
