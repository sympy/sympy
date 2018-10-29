from collections import defaultdict
from sympy import Indexed, IndexedBase, Tuple, Sum, Add
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.mul import Mul
from sympy.core.compatibility import accumulate, default_sort_key
from sympy.combinatorics import Permutation
from sympy.matrices.expressions import MatMul, Trace, Transpose, MatrixSymbol
from sympy.matrices.expressions.matexpr import MatrixExpr, MatrixElement
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.tensor.array import NDimArray


class CodegenArrayContraction(Basic):
    """
    This class is meant to represent contractions of arrays in a form easily
    processable by the code printers.
    """
    def __new__(cls, expr, *contraction_indices, **kwargs):
        contraction_indices = _sort_contraction_indices(contraction_indices)
        expr = _sympify(expr)

        if isinstance(expr, CodegenArrayContraction):
            return cls._flatten(expr, *contraction_indices)

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
    def _get_mapping_from_contraction_indices(expr, *contraction_indices):
        if isinstance(expr, CodegenArrayTensorProduct):
            args = expr.args
            ranks = expr.ranks
        else:
            args = [expr]
            ranks = [get_rank(expr)]
        mapping = {}
        counter = 0
        for i, rank in enumerate(ranks):
            for j in range(rank):
                mapping[counter] = (i, j)
                counter += 1
        return mapping, ranks

    @staticmethod
    def _flatten(expr, *outer_contraction_indices):
        inner_contraction_indices = expr.contraction_indices
        all_inner = [j for i in inner_contraction_indices for j in i]
        all_inner.sort()
        # TODO: add API for total rank and cumulative rank:
        total_rank = get_rank(expr)
        inner_rank = len(all_inner)
        outer_rank = total_rank - inner_rank
        shifts = [0 for i in range(outer_rank)]
        counter = 0
        pointer = 0
        for i in range(outer_rank):
            while pointer < inner_rank and counter >= all_inner[pointer]:
                counter += 1
                pointer += 1
            shifts[i] += pointer
            counter += 1
        outer_contraction_indices = tuple(tuple(shifts[j] + j for j in i) for i in outer_contraction_indices)
        contraction_indices = inner_contraction_indices + outer_contraction_indices
        return CodegenArrayContraction(expr.expr, *contraction_indices)

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
        {0: {1: (1, 0)}, 1: {0: (0, 1), 1: (2, 0)}, 2: {0: (1, 1), 1: (3, 0)}, 3: {0: (2, 1)}}

        This dictionary is interpreted as follows: argument in position 0 (i.e.
        matrix `A`) has its second index (i.e. 1) contracted to `(1, 0)`, that
        is argument in position 1 (matrix `B`) on the first index slot of `B`,
        this is the contraction provided by the index `j` from `A`.

        The argument in position 1 (that is, matrix `B`) has two contractions,
        the ones provided by the indices `j` and `k`, respectively the first
        and second indices (0 and 1 in the sub-dict).  The link `(0, 1)` and
        `(2, 0)` respectively. `(0, 1)` is the index slot 1 (the 2nd) of
        argument in position 0 (that is, `A_{\ldot j}`), and so on.
        """
        contraction_tuples = self._get_contraction_tuples()
        dlinks = defaultdict(dict)
        for links in contraction_tuples:
            if len(links) > 2:
                raise NotImplementedError("three or more axes contracted at the same time")
            (arg1, pos1), (arg2, pos2) = links
            dlinks[arg1][pos1] = (arg2, pos2)
            dlinks[arg2][pos2] = (arg1, pos1)
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
        >>> from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)
        >>> C = MatrixSymbol("C", N, N)
        >>> D = MatrixSymbol("D", N, N)

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

        Recognize some more complex traces:
        >>> expr = Sum(A[i, j]*B[j, i], (i, 0, N-1), (j, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg._recognize_matrix_mul_lines()
        [(None, None, [Trace(A*B)])]

        More complicated expressions:

        >>> expr = Sum(A[i, j]*B[k, j]*A[l, k], (j, 0, N-1), (k, 0, N-1))
        >>> cg = CodegenArrayContraction.from_summation(expr)
        >>> cg._recognize_matrix_mul_lines()
        [(i, l, [A, B.T, A.T])]

        Expressions constructed from matrix expressions do not contain literal
        indices, the positions of free indices are returned instead:

        >>> expr = A*B
        >>> cg = CodegenArrayContraction.from_MatMul(expr)
        >>> cg.free_indices
        {0: 0, 3: 3}
        >>> cg._recognize_matrix_mul_lines()
        [(0, 3, [A, B])]

        If more than one line of matrix multiplications is detected, return
        separate matrix multiplication factors:

        >>> cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B, C, D), (1, 2), (5, 6))
        >>> cg._recognize_matrix_mul_lines()
        [(0, 3, [A, B]), (4, 7, [C, D])]

        The two lines have free indices at axes 0, 3 and 4, 7, respectively.
        """
        expr = self.expr
        if not isinstance(expr, CodegenArrayTensorProduct):
            args = [expr]
        else:
            args = expr.args
        free_indices = self.free_indices
        dlinks = self._get_contraction_links()
        return_list = []
        while dlinks:
            if free_indices:
                if first_indices:
                    first_index = first_indices[0]
                    starting_argind = free_indices.pop(first_index)
                else:
                    first_index, starting_argind = min(free_indices.items(), key=lambda x: x[1])
                    free_indices.pop(first_index)
                starting_argind, starting_pos = self._mapping[starting_argind]
            else:
                first_index = None
                starting_argind = min(dlinks)
                starting_pos = 0  #max(dlinks[starting_argind])
            current_argind, current_pos = starting_argind, starting_pos
            matmul_args = []
            prev_argind = None
            prev_pos = None
            last_index = None
            while True:
                elem = args[current_argind]
                if current_pos == 1:
                    elem = Transpose(elem)
                matmul_args.append(elem)
                if current_argind not in dlinks:
                    break
                other_pos = 1 - current_pos
                link_dict = dlinks.pop(current_argind)
                if other_pos not in link_dict:
                    if free_indices:
                        last_index = [i for i, j in free_indices.items() if self._mapping[j] == (current_argind, other_pos)][0]
                    else:
                        last_index = None
                    break
                if len(link_dict) > 2:
                    raise NotImplementedError("not a matrix multiplication line")
                prev_argind = current_argind
                prev_pos = current_pos
                # Get the last element of `link_list` as the next link. The last
                # element is the correct start for trace expressions:
                current_argind, current_pos = link_dict[other_pos]
                if current_argind == starting_argind:
                    # This is a trace:
                    if len(matmul_args) > 1:
                        matmul_args = [Trace(MatMul(*matmul_args, check=True))]
                    else:
                        matmul_args = [Trace(*matmul_args)]
                    break
            dlinks.pop(starting_argind, None)
            free_indices.pop(last_index, None)
            return_list.append((first_index, last_index, matmul_args))
        return return_list

    def _recognize_addition_of_mul_lines(self):
        res = [recurse_expr(i) for i in expr.args]
        d = collections.defaultdict(list)
        for res_addend in res:
            scalar = 1
            for elem, indices in res_addend:
                if indices is None:
                    scalar = elem
                    continue
                indices = tuple(sorted(indices, key=default_sort_key))
                d[indices].append(scalar*remove_matelement(elem, *indices))
                scalar = 1
        return [(MatrixElement(Add.fromiter(v), *k), k) for k, v in d.items()]

    @staticmethod
    def from_summation(summation):
        expr, indices = CodegenArrayContraction._parse_indexed(summation)
        return expr

    @staticmethod
    def _parse_indexed(summation):
        from sympy import Indexed
        from sympy.matrices.expressions.matexpr import MatrixElement
        function = summation.function
        indices = summation.variables
        free_indices = {}
        if function.is_Mul:
            function_args = function.args
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
        indices_ret = list(free_indices)
        indices_ret.sort(key=lambda x: free_indices[x])
        return CodegenArrayContraction(
                CodegenArrayTensorProduct(*args),
                *[tuple(v) for v in axes_contraction.values()],
                free_indices=free_indices
            ), tuple(indices_ret)

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
    r"""
    Class to represent the tensor product of array-like objects.
    """
    def __new__(cls, *args):
        args = [_sympify(arg) for arg in args]
        args = cls._flatten(args)
        ranks = [get_rank(arg) for arg in args]

        if len(args) == 1:
            return args[0]

        # If there are contraction objects inside, transform the whole
        # expression into `CodegenArrayContraction`:
        contractions = {i: arg for i, arg in enumerate(args) if isinstance(arg, CodegenArrayContraction)}
        if contractions:
            cumulative_ranks = list(accumulate([0] + ranks))[:-1]
            tp = cls(*[arg.expr if isinstance(arg, CodegenArrayContraction) else arg for arg in args])
            contraction_indices = [tuple(cumulative_ranks[i] + k for k in j) for i, arg in contractions.items() for j in arg.contraction_indices]
            return CodegenArrayContraction(tp, *contraction_indices)

        obj = Basic.__new__(cls, *args)
        obj._ranks = ranks
        return obj

    @property
    def ranks(self):
        return self._ranks[:]

    @classmethod
    def _flatten(cls, args):
        args = [i for arg in args for i in (arg.args if isinstance(arg, cls) else [arg])]
        return args


class CodegenArrayElementwiseAdd(Basic):
    r"""
    Class for elementwise array additions.
    """
    def __new__(cls, *args):
        args = [_sympify(arg) for arg in args]
        obj = Basic.__new__(cls, *args)
        return obj


class CodegenArrayPermuteDims(Basic):
    r"""
    Class to represent permutation of axes of arrays.
    """
    def __new__(cls, expr, permutation):
        from sympy.combinatorics import Permutation
        expr = _sympify(expr)
        permutation = Permutation(permutation)
        plist = permutation.args[0]
        if plist == sorted(plist):
            return expr
        obj = Basic.__new__(cls, expr, permutation)
        return obj

    @property
    def expr(self):
        return self.args[0]

    @property
    def permutation(self):
        return self.args[1]


class CodegenArrayDiagonal(Basic):
    r"""
    Class to represent the diagonal operator.

    In a 2-dimensional array it returns the diagonal, this looks like the
    operation:

    `i \rightarrow A_{ii}`

    The diagonal over axes 1 and 2 (the second and third) of the tensor product
    of two 2-dimensional arrays `A \otimes B` is

    `\Big[ A_{ab} B_{cd} \Big]_{abcd} \rightarrow \Big[ A_{ai} B_{id} \Big]_{adi}`

    In this last example the array expression has been reduced from
    4-dimensional to 3-dimensional. Notice that no contraction has occurred,
    rather there is a new index `i` for the diagonal, contraction would have
    reduced the array to 2 dimensions.

    Notice that the diagonalized out dimensions are added as new dimensions at
    the end of the indices.
    """
    def __new__(cls, expr, *diagonal_indices):
        expr = _sympify(expr)
        diagonal_indices = [Tuple(*sorted(i)) for i in diagonal_indices]
        if isinstance(expr, CodegenArrayDiagonal):
            return cls._flatten(expr, *diagonal_indices)
        obj = Basic.__new__(cls, expr, *diagonal_indices)
        obj._mapping, obj._ranks = CodegenArrayContraction._get_mapping_from_contraction_indices(expr, *diagonal_indices)
        return obj

    @property
    def expr(self):
        return self.args[0]

    @property
    def diagonal_indices(self):
        return self.args[1:]

    @property
    def ranks(self):
        return self._ranks

    @staticmethod
    def _flatten(expr, *outer_diagonal_indices):
        inner_diagonal_indices = expr.diagonal_indices
        all_inner = [j for i in inner_diagonal_indices for j in i]
        all_inner.sort()
        # TODO: add API for total rank and cumulative rank:
        total_rank = get_rank(expr)
        inner_rank = len(all_inner)
        outer_rank = total_rank - inner_rank
        shifts = [0 for i in range(outer_rank)]
        counter = 0
        pointer = 0
        for i in range(outer_rank):
            while pointer < inner_rank and counter >= all_inner[pointer]:
                counter += 1
                pointer += 1
            shifts[i] += pointer
            counter += 1
        outer_diagonal_indices = tuple(tuple(shifts[j] + j for j in i) for i in outer_diagonal_indices)
        diagonal_indices = inner_diagonal_indices + outer_diagonal_indices
        return CodegenArrayDiagonal(expr.expr, *diagonal_indices)


def get_rank(expr):
    if isinstance(expr, (MatrixExpr, MatrixElement)):
        return 2
    if isinstance(expr, CodegenArrayContraction):
        return sum(expr.ranks)
    if isinstance(expr, CodegenArrayDiagonal):
        return sum(expr.ranks)
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


def _sort_contraction_indices(pairing_indices):
    pairing_indices = [Tuple(*sorted(i)) for i in pairing_indices]
    pairing_indices.sort(key=lambda x: min(x))
    return pairing_indices


def _codegen_array_parse(expr):
    if isinstance(expr, Sum):
        return CodegenArrayContraction._parse_indexed(expr)
    if isinstance(expr, Mul):
        args, indices = zip(*[_codegen_array_parse(arg) for arg in expr.args])
        ranks = [get_rank(arg) for arg in args]
        total_rank = 0
        axes_contraction = defaultdict(list)
        flattened_indices = []
        held_back = []
        kronecker_delta_repl = {}
        for arg in args:
            if not isinstance(arg, KroneckerDelta):
                continue
            # Diagonalize two indices:
            i, j = arg.indices
            kindices = frozenset(arg.indices)
            kronecker_delta_repl[i] = kindices
            kronecker_delta_repl[j] = kindices
        args = [arg for arg in args if not isinstance(arg, KroneckerDelta)]
        for arg, loc_indices in zip(args, indices):
            for i, ind in enumerate(loc_indices):
                ind = kronecker_delta_repl.get(ind, ind)
                if ind in flattened_indices:
                    other_pos = flattened_indices.index(ind)
                    if other_pos not in axes_contraction[ind]:
                        axes_contraction[ind].append(other_pos)
                    axes_contraction[ind].append(total_rank + i)
                    flattened_indices[other_pos] = None
                    held_back.append(ind)
                else:
                    flattened_indices.append(ind)
            total_rank += ranks[-1]
        tp = CodegenArrayTensorProduct(*args)
        flattened_indices = tuple(i for i in flattened_indices if i is not None)
        flattened_indices += tuple(held_back)
        if axes_contraction:
            return (CodegenArrayDiagonal(tp,
                *[Tuple(*v) for v in axes_contraction.values()]), flattened_indices)
        else:
            return tp, flattened_indices
    if isinstance(expr, MatrixElement):
        return expr.args[0], expr.args[1:]
    if isinstance(expr, Indexed):
        return expr.base, expr.indices
    if isinstance(expr, IndexedBase):
        raise NotImplementedError
    if isinstance(expr, KroneckerDelta):
        return expr, expr.indices
        raise NotImplementedError
    if isinstance(expr, Add):
        args, indices = zip(*[_codegen_array_parse(arg) for arg in expr.args])
        args = list(args)
        # Check if all indices are compatible. Otherwise expand the dimensions:
        index0set = set(indices[0])
        index0 = indices[0]
        for i in range(1, len(args)):
            if set(indices[i]) != index0set:
                raise NotImplementedError("indices must be the same")
            permutation = Permutation([index0.index(j) for j in indices[i]])
            # Perform index permutations:
            args[i] = CodegenArrayPermuteDims(args[i], permutation)
        return CodegenArrayElementwiseAdd(*args), index0
    raise NotImplementedError("could not recognize expression %s" % expr)


class _RecognizeMatAdd(list):
    def __repr__(self):
        return "_RecognizeMatAdd(%s)" % super(_RecognizeMatAdd, self).__repr__()

class _RecognizeMatMul(list):
    pass

def _recognize_matrix_expression(expr):
    # TODO: expr has to be a CodegenArray... type
    if isinstance(expr, CodegenArrayContraction):
        return expr._recognize_addition_of_mul_lines()
    elif isinstance(expr, CodegenArrayElementwiseAdd):
        add_args = _RecognizeMatAdd()
        for arg in expr.args:
            add_args.append(_recognize_matrix_expression(arg))
            #if isinstance(arg, MatrixSymbol):
                #add_args.append(arg)
        return add_args
    elif isinstance(expr, (MatrixSymbol, IndexedBase)):
        return expr
    elif isinstance(expr, CodegenArrayPermuteDims):
        if expr.permutation.args[0] == [1, 0]:
            return Trace(expr.expr)
        else:
            raise NotImplementedError
    elif isinstance(expr, CodegenArrayTensorProduct):
        pass
    raise NotImplementedError
