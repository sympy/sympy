import bisect
import itertools
import operator
from functools import reduce, singledispatch
from itertools import accumulate
from collections import defaultdict

from sympy import Indexed, IndexedBase, Tuple, Sum, Add, S, Integer, diagonalize_vector, DiagMatrix, ZeroMatrix, Pow, \
    MatPow, HadamardProduct, HadamardPower, tensorcontraction, tensorproduct, permutedims, tensordiagonal, Lambda, \
    Dummy, symbols, Function
from sympy.combinatorics import Permutation
from sympy.combinatorics.permutations import _af_invert
from sympy.core.basic import Basic
from sympy.core.compatibility import default_sort_key
from sympy.core.mul import Mul
from sympy.core.sympify import _sympify
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.matrices.common import MatrixCommon
from sympy.matrices.expressions import (MatAdd, MatMul, Trace, Transpose)
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.matrices.expressions.matexpr import MatrixExpr, MatrixElement
from sympy.tensor.array import NDimArray
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray, _ArrayExpr


class _CodegenArrayAbstract(Basic):

    @property
    def subranks(self):
        """
        Returns the ranks of the objects in the uppermost tensor product inside
        the current object.  In case no tensor products are contained, return
        the atomic ranks.

        Examples
        ========

        >>> from sympy.codegen.array_utils import CodegenArrayTensorProduct, CodegenArrayContraction
        >>> from sympy import MatrixSymbol
        >>> M = MatrixSymbol("M", 3, 3)
        >>> N = MatrixSymbol("N", 3, 3)
        >>> P = MatrixSymbol("P", 3, 3)

        Important: do not confuse the rank of the matrix with the rank of an array.

        >>> tp = CodegenArrayTensorProduct(M, N, P)
        >>> tp.subranks
        [2, 2, 2]

        >>> co = CodegenArrayContraction(tp, (1, 2), (3, 4))
        >>> co.subranks
        [2, 2, 2]
        """
        return self._subranks[:]

    def subrank(self):
        """
        The sum of ``subranks``.
        """
        return sum(self.subranks)

    @property
    def shape(self):
        return self._shape


class CodegenArrayContraction(_CodegenArrayAbstract):
    r"""
    This class is meant to represent contractions of arrays in a form easily
    processable by the code printers.
    """

    def __new__(cls, expr, *contraction_indices, **kwargs):
        contraction_indices = _sort_contraction_indices(contraction_indices)
        expr = _sympify(expr)

        if len(contraction_indices) == 0:
            return expr

        if isinstance(expr, CodegenArrayContraction):
            return cls._flatten(expr, *contraction_indices)

        if isinstance(expr, (ZeroArray, ZeroMatrix)):
            contraction_indices_flat = [j for i in contraction_indices for j in i]
            shape = [e for i, e in enumerate(expr.shape) if i not in contraction_indices_flat]
            return ZeroArray(*shape)

        if isinstance(expr, CodegenArrayPermuteDims):
            return cls._handle_nested_permute_dims(expr, *contraction_indices)

        if isinstance(expr, CodegenArrayTensorProduct):
            expr, contraction_indices = cls._sort_fully_contracted_args(expr, contraction_indices)
            expr, contraction_indices = cls._lower_contraction_to_addends(expr, contraction_indices)
            if len(contraction_indices) == 0:
                return expr

        if isinstance(expr, CodegenArrayDiagonal):
            return cls._handle_nested_diagonal(expr, *contraction_indices)

        if isinstance(expr, CodegenArrayElementwiseAdd):
            return CodegenArrayElementwiseAdd(*[CodegenArrayContraction(i, *contraction_indices) for i in expr.args])

        obj = Basic.__new__(cls, expr, *contraction_indices)
        obj._subranks = _get_subranks(expr)
        obj._mapping = _get_mapping_from_subranks(obj._subranks)

        free_indices_to_position = {i: i for i in range(sum(obj._subranks)) if all([i not in cind for cind in contraction_indices])}
        obj._free_indices_to_position = free_indices_to_position

        shape = expr.shape
        cls._validate(expr, *contraction_indices)
        if shape:
            shape = tuple(shp for i, shp in enumerate(shape) if not any(i in j for j in contraction_indices))
        obj._shape = shape
        return obj

    def __mul__(self, other):
        if other == 1:
            return self
        else:
            raise NotImplementedError("Product of N-dim arrays is not uniquely defined. Use another method.")

    def __rmul__(self, other):
        if other == 1:
            return self
        else:
            raise NotImplementedError("Product of N-dim arrays is not uniquely defined. Use another method.")

    @staticmethod
    def _validate(expr, *contraction_indices):
        shape = expr.shape
        if shape is None:
            return

        # Check that no contraction happens when the shape is mismatched:
        for i in contraction_indices:
            if len({shape[j] for j in i if shape[j] != -1}) != 1:
                raise ValueError("contracting indices of different dimensions")

    @classmethod
    def _push_indices_down(cls, contraction_indices, indices):
        flattened_contraction_indices = [j for i in contraction_indices for j in i]
        flattened_contraction_indices.sort()
        transform = _build_push_indices_down_func_transformation(flattened_contraction_indices)
        return _apply_recursively_over_nested_lists(transform, indices)

    @classmethod
    def _push_indices_up(cls, contraction_indices, indices):
        flattened_contraction_indices = [j for i in contraction_indices for j in i]
        flattened_contraction_indices.sort()
        transform = _build_push_indices_up_func_transformation(flattened_contraction_indices)
        return _apply_recursively_over_nested_lists(transform, indices)

    @classmethod
    def _lower_contraction_to_addends(cls, expr, contraction_indices):
        if isinstance(expr, CodegenArrayElementwiseAdd):
            raise NotImplementedError()
        if not isinstance(expr, CodegenArrayTensorProduct):
            return expr, contraction_indices
        subranks = expr.subranks
        cumranks = list(accumulate([0] + subranks))
        contraction_indices_remaining = []
        contraction_indices_args = [[] for i in expr.args]
        backshift = set([])
        for i, contraction_group in enumerate(contraction_indices):
            for j in range(len(expr.args)):
                if not isinstance(expr.args[j], CodegenArrayElementwiseAdd):
                    continue
                if all(cumranks[j] <= k < cumranks[j+1] for k in contraction_group):
                    contraction_indices_args[j].append([k - cumranks[j] for k in contraction_group])
                    backshift.update(contraction_group)
                    break
            else:
                contraction_indices_remaining.append(contraction_group)
        if len(contraction_indices_remaining) == len(contraction_indices):
            return expr, contraction_indices
        total_rank = get_rank(expr)
        shifts = list(accumulate([1 if i in backshift else 0 for i in range(total_rank)]))
        contraction_indices_remaining = [Tuple.fromiter(j - shifts[j] for j in i) for i in contraction_indices_remaining]
        ret = CodegenArrayTensorProduct(*[
            CodegenArrayContraction(arg, *contr) for arg, contr in zip(expr.args, contraction_indices_args)
        ])
        return ret, contraction_indices_remaining

    def split_multiple_contractions(self):
        """
        Recognize multiple contractions and attempt at rewriting them as paired-contractions.
        """
        from sympy import ask, Q

        contraction_indices = self.contraction_indices
        if isinstance(self.expr, CodegenArrayTensorProduct):
            args = list(self.expr.args)
        else:
            args = [self.expr]
        # TODO: unify API, best location in CodegenArrayTensorProduct
        subranks = [get_rank(i) for i in args]
        # TODO: unify API
        mapping = _get_mapping_from_subranks(subranks)
        reverse_mapping = {v:k for k, v in mapping.items()}
        new_contraction_indices = []
        for indl, links in enumerate(contraction_indices):
            if len(links) <= 2:
                new_contraction_indices.append(links)
                continue

            # Check multiple contractions:
            #
            # Examples:
            #
            # * `A_ij b_j0 C_jk` ===> `A*DiagMatrix(b)*C`
            #
            # Care for:
            # - matrix being diagonalized (i.e. `A_ii`)
            # - vectors being diagonalized (i.e. `a_i0`)

            # Also consider the case of diagonal matrices being contracted:
            current_dimension = self.expr.shape[links[0]]

            tuple_links = [mapping[i] for i in links]
            arg_indices, arg_positions = zip(*tuple_links)
            args_updates = {}
            if len(arg_indices) != len(set(arg_indices)):
                # Maybe trace should be supported?
                raise NotImplementedError
            not_vectors = []
            vectors = []
            for arg_ind, arg_pos in tuple_links:
                mat = args[arg_ind]
                other_arg_pos = 1-arg_pos
                other_arg_abs = reverse_mapping[arg_ind, other_arg_pos]
                if (((1 not in mat.shape) and (not ask(Q.diagonal(mat)))) or
                    ((current_dimension == 1) is True and mat.shape != (1, 1)) or
                    any([other_arg_abs in l for li, l in enumerate(contraction_indices) if li != indl])
                ):
                    not_vectors.append((arg_ind, arg_pos))
                    continue
                args_updates[arg_ind] = diagonalize_vector(mat)
                vectors.append((arg_ind, arg_pos))
                vectors.append((arg_ind, 1-arg_pos))
            if len(not_vectors) > 2:
                new_contraction_indices.append(links)
                continue
            if len(not_vectors) == 0:
                new_sequence = vectors[:1] + vectors[2:]
            elif len(not_vectors) == 1:
                new_sequence = not_vectors[:1] + vectors[:-1]
            else:
                new_sequence = not_vectors[:1] + vectors + not_vectors[1:]
            for i in range(0, len(new_sequence) - 1, 2):
                arg1, pos1 = new_sequence[i]
                arg2, pos2 = new_sequence[i+1]
                if arg1 == arg2:
                    raise NotImplementedError
                    continue
                abspos1 = reverse_mapping[arg1, pos1]
                abspos2 = reverse_mapping[arg2, pos2]
                new_contraction_indices.append((abspos1, abspos2))
            for ind, newarg in args_updates.items():
                args[ind] = newarg
        return CodegenArrayContraction(
            CodegenArrayTensorProduct(*args),
            *new_contraction_indices
        )

    def flatten_contraction_of_diagonal(self):
        if not isinstance(self.expr, CodegenArrayDiagonal):
            return self
        contraction_down = self.expr._push_indices_down(self.expr.diagonal_indices, self.contraction_indices)
        new_contraction_indices = []
        diagonal_indices = self.expr.diagonal_indices[:]
        for i in contraction_down:
            contraction_group = list(i)
            for j in i:
                diagonal_with = [k for k in diagonal_indices if j in k]
                contraction_group.extend([l for k in diagonal_with for l in k])
                diagonal_indices = [k for k in diagonal_indices if k not in diagonal_with]
            new_contraction_indices.append(sorted(set(contraction_group)))

        new_contraction_indices = CodegenArrayDiagonal._push_indices_up(diagonal_indices, new_contraction_indices)
        return CodegenArrayContraction(
            CodegenArrayDiagonal(
                self.expr.expr,
                *diagonal_indices
            ),
            *new_contraction_indices
        )

    @staticmethod
    def _get_free_indices_to_position_map(free_indices, contraction_indices):
        free_indices_to_position = {}
        flattened_contraction_indices = [j for i in contraction_indices for j in i]
        counter = 0
        for ind in free_indices:
            while counter in flattened_contraction_indices:
                counter += 1
            free_indices_to_position[ind] = counter
            counter += 1
        return free_indices_to_position

    @staticmethod
    def _get_index_shifts(expr):
        """
        Get the mapping of indices at the positions before the contraction
        occurs.

        Examples
        ========

        >>> from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
        >>> from sympy import MatrixSymbol
        >>> M = MatrixSymbol("M", 3, 3)
        >>> N = MatrixSymbol("N", 3, 3)
        >>> cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), [1, 2])
        >>> cg._get_index_shifts(cg)
        [0, 2]

        Indeed, ``cg`` after the contraction has two dimensions, 0 and 1. They
        need to be shifted by 0 and 2 to get the corresponding positions before
        the contraction (that is, 0 and 3).
        """
        inner_contraction_indices = expr.contraction_indices
        all_inner = [j for i in inner_contraction_indices for j in i]
        all_inner.sort()
        # TODO: add API for total rank and cumulative rank:
        total_rank = _get_subrank(expr)
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
        return shifts

    @staticmethod
    def _convert_outer_indices_to_inner_indices(expr, *outer_contraction_indices):
        shifts = CodegenArrayContraction._get_index_shifts(expr)
        outer_contraction_indices = tuple(tuple(shifts[j] + j for j in i) for i in outer_contraction_indices)
        return outer_contraction_indices

    @staticmethod
    def _flatten(expr, *outer_contraction_indices):
        inner_contraction_indices = expr.contraction_indices
        outer_contraction_indices = CodegenArrayContraction._convert_outer_indices_to_inner_indices(expr, *outer_contraction_indices)
        contraction_indices = inner_contraction_indices + outer_contraction_indices
        return CodegenArrayContraction(expr.expr, *contraction_indices)

    @classmethod
    def _handle_nested_permute_dims(cls, expr, *contraction_indices):
        permutation = expr.permutation
        plist = permutation.array_form
        new_contraction_indices = [tuple(permutation(j) for j in i) for i in contraction_indices]
        new_plist = [i for i in plist if all(i not in j for j in new_contraction_indices)]
        new_plist = cls._push_indices_up(new_contraction_indices, new_plist)
        return CodegenArrayPermuteDims(
            CodegenArrayContraction(expr.expr, *new_contraction_indices),
            Permutation(new_plist)
        )

    @classmethod
    def _handle_nested_diagonal(cls, expr: 'CodegenArrayDiagonal', *contraction_indices):
        diagonal_indices = list(expr.diagonal_indices)
        down_contraction_indices = expr._push_indices_down(expr.diagonal_indices, contraction_indices, get_rank(expr.expr))
        # Flatten diagonally contracted indices:
        down_contraction_indices = [[k for j in i for k in (j if isinstance(j, (tuple, Tuple)) else [j])] for i in down_contraction_indices]
        new_contraction_indices = []
        for contr_indgrp in down_contraction_indices:
            ind = contr_indgrp[:]
            for j, diag_indgrp in enumerate(diagonal_indices):
                if diag_indgrp is None:
                    continue
                if any(i in diag_indgrp for i in contr_indgrp):
                    ind.extend(diag_indgrp)
                    diagonal_indices[j] = None
            new_contraction_indices.append(sorted(set(ind)))

        new_diagonal_indices_down = [i for i in diagonal_indices if i is not None]
        new_diagonal_indices = CodegenArrayContraction._push_indices_up(new_contraction_indices, new_diagonal_indices_down)
        return CodegenArrayDiagonal(
            CodegenArrayContraction(expr.expr, *new_contraction_indices),
            *new_diagonal_indices
        )

    @classmethod
    def _sort_fully_contracted_args(cls, expr, contraction_indices):
        if expr.shape is None:
            return expr, contraction_indices
        cumul = list(accumulate([0] + expr.subranks))
        index_blocks = [list(range(cumul[i], cumul[i+1])) for i in range(len(expr.args))]
        contraction_indices_flat = {j for i in contraction_indices for j in i}
        fully_contracted = [all(j in contraction_indices_flat for j in range(cumul[i], cumul[i+1])) for i, arg in enumerate(expr.args)]
        new_pos = sorted(range(len(expr.args)), key=lambda x: (0, default_sort_key(expr.args[x])) if fully_contracted[x] else (1,))
        new_args = [expr.args[i] for i in new_pos]
        new_index_blocks_flat = [j for i in new_pos for j in index_blocks[i]]
        index_permutation_array_form = _af_invert(new_index_blocks_flat)
        new_contraction_indices = [tuple(index_permutation_array_form[j] for j in i) for i in contraction_indices]
        new_contraction_indices = _sort_contraction_indices(new_contraction_indices)
        return CodegenArrayTensorProduct(*new_args), new_contraction_indices

    def _get_contraction_tuples(self):
        r"""
        Return tuples containing the argument index and position within the
        argument of the index position.

        Examples
        ========

        >>> from sympy import MatrixSymbol
        >>> from sympy.abc import N
        >>> from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)

        >>> cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (1, 2))
        >>> cg._get_contraction_tuples()
        [[(0, 1), (1, 0)]]

        Notes
        =====

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
        # TODO: check that `expr` has `.subranks`:
        ranks = expr.subranks
        cumulative_ranks = [0] + list(accumulate(ranks))
        return [tuple(cumulative_ranks[j]+k for j, k in i) for i in contraction_tuples]

    @property
    def free_indices(self):
        return self._free_indices[:]

    @property
    def free_indices_to_position(self):
        return dict(self._free_indices_to_position)

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
        ranks = expr.subranks
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

        >>> from sympy import MatrixSymbol
        >>> from sympy.abc import N
        >>> from sympy.codegen.array_utils import parse_matrix_expression
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)
        >>> C = MatrixSymbol("C", N, N)
        >>> D = MatrixSymbol("D", N, N)

        >>> cg = parse_matrix_expression(C*D*A*B)
        >>> cg
        CodegenArrayContraction(CodegenArrayTensorProduct(A, D, C, B), (0, 3), (1, 6), (2, 5))
        >>> cg.sort_args_by_name()
        CodegenArrayContraction(CodegenArrayTensorProduct(A, D, B, C), (0, 3), (1, 4), (2, 7))
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

        >>> from sympy import MatrixSymbol
        >>> from sympy.abc import N
        >>> from sympy.codegen.array_utils import parse_matrix_expression
        >>> A = MatrixSymbol("A", N, N)
        >>> B = MatrixSymbol("B", N, N)
        >>> C = MatrixSymbol("C", N, N)
        >>> D = MatrixSymbol("D", N, N)

        Matrix multiplications are pairwise contractions between neighboring
        matrices:

        `A_{ij} B_{jk} C_{kl} D_{lm}`

        >>> cg = parse_matrix_expression(A*B*C*D)
        >>> cg
        CodegenArrayContraction(CodegenArrayTensorProduct(B, C, A, D), (0, 5), (1, 2), (3, 6))
        >>> cg._get_contraction_links()
        {0: {0: (2, 1), 1: (1, 0)}, 1: {0: (0, 1), 1: (3, 0)}, 2: {1: (0, 0)}, 3: {0: (1, 1)}}

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
        args, dlinks = _get_contraction_links([self], self.subranks, *self.contraction_indices)
        return dlinks

    def as_explicit(self):
        return tensorcontraction(self.expr.as_explicit(), *self.contraction_indices)


def get_shape(expr):
    if hasattr(expr, "shape"):
        return expr.shape
    return ()


class CodegenArrayTensorProduct(_CodegenArrayAbstract):
    r"""
    Class to represent the tensor product of array-like objects.
    """

    def __new__(cls, *args):
        args = [_sympify(arg) for arg in args]
        args = cls._flatten(args)
        ranks = [_get_subrank(arg) for arg in args]

        # Check if there are nested permutation and lift them up:
        permutation_cycles = []
        for i, arg in enumerate(args):
            if not isinstance(arg, CodegenArrayPermuteDims):
                continue
            permutation_cycles.extend([[k + sum(ranks[:i]) for k in j] for j in arg.permutation.cyclic_form])
            args[i] = arg.expr
        if permutation_cycles:
            return CodegenArrayPermuteDims(CodegenArrayTensorProduct(*args), Permutation(sum(ranks)-1)*Permutation(permutation_cycles))

        if len(args) == 1:
            return args[0]

        # If any object is a ZeroArray, return a ZeroArray:
        if any(isinstance(arg, (ZeroArray, ZeroMatrix)) for arg in args):
            shapes = reduce(operator.add, [get_shape(i) for i in args], ())
            return ZeroArray(*shapes)

        # If there are contraction objects inside, transform the whole
        # expression into `CodegenArrayContraction`:
        contractions = {i: arg for i, arg in enumerate(args) if isinstance(arg, CodegenArrayContraction)}
        if contractions:
            cumulative_ranks = list(accumulate([0] + ranks))[:-1]
            tp = cls(*[arg.expr if isinstance(arg, CodegenArrayContraction) else arg for arg in args])
            contraction_indices = [tuple(cumulative_ranks[i] + k for k in j) for i, arg in contractions.items() for j in arg.contraction_indices]
            return CodegenArrayContraction(tp, *contraction_indices)

        #newargs = [i for i in args if hasattr(i, "shape")]
        #coeff = reduce(lambda x, y: x*y, [i for i in args if not hasattr(i, "shape")], S.One)
        #newargs[0] *= coeff

        obj = Basic.__new__(cls, *args)
        obj._subranks = ranks
        shapes = [get_shape(i) for i in args]

        if any(i is None for i in shapes):
            obj._shape = None
        else:
            obj._shape = tuple(j for i in shapes for j in i)
        return obj

    @classmethod
    def _flatten(cls, args):
        args = [i for arg in args for i in (arg.args if isinstance(arg, cls) else [arg])]
        return args

    def as_explicit(self):
        return tensorproduct(*[arg.as_explicit() if hasattr(arg, "as_explicit") else arg for arg in self.args])


class CodegenArrayElementwiseAdd(_CodegenArrayAbstract):
    r"""
    Class for elementwise array additions.
    """

    def __new__(cls, *args):
        args = [_sympify(arg) for arg in args]
        ranks = [get_rank(arg) for arg in args]
        ranks = list(set(ranks))
        if len(ranks) != 1:
            raise ValueError("summing arrays of different ranks")
        shapes = [arg.shape for arg in args]
        if len({i for i in shapes if i is not None}) > 1:
            raise ValueError("mismatching shapes in addition")

        # Flatten:
        args = cls._flatten_args(args)

        args = [arg for arg in args if not isinstance(arg, (ZeroArray, ZeroMatrix))]
        if len(args) == 0:
            if any(i for i in shapes if i is None):
                raise NotImplementedError("cannot handle addition of ZeroMatrix/ZeroArray and undefined shape object")
            return ZeroArray(*shapes[0])
        elif len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj._subranks = ranks
        if any(i is None for i in shapes):
            obj._shape = None
        else:
            obj._shape = shapes[0]
        return obj

    @classmethod
    def _flatten_args(cls, args):
        new_args = []
        for arg in args:
            if isinstance(arg, CodegenArrayElementwiseAdd):
                new_args.extend(arg.args)
            else:
                new_args.append(arg)
        return new_args

    def as_explicit(self):
        return Add.fromiter([arg.as_explicit() for arg in self.args])


class CodegenArrayPermuteDims(_CodegenArrayAbstract):
    r"""
    Class to represent permutation of axes of arrays.

    Examples
    ========

    >>> from sympy.codegen.array_utils import CodegenArrayPermuteDims
    >>> from sympy import MatrixSymbol
    >>> M = MatrixSymbol("M", 3, 3)
    >>> cg = CodegenArrayPermuteDims(M, [1, 0])

    The object ``cg`` represents the transposition of ``M``, as the permutation
    ``[1, 0]`` will act on its indices by switching them:

    `M_{ij} \Rightarrow M_{ji}`

    This is evident when transforming back to matrix form:

    >>> from sympy.codegen.array_utils import recognize_matrix_expression
    >>> recognize_matrix_expression(cg)
    M.T

    >>> N = MatrixSymbol("N", 3, 2)
    >>> cg = CodegenArrayPermuteDims(N, [1, 0])
    >>> cg.shape
    (2, 3)

    Permutations of tensor products are simplified in order to achieve a
    standard form:

    >>> from sympy.codegen.array_utils import CodegenArrayTensorProduct
    >>> M = MatrixSymbol("M", 4, 5)
    >>> tp = CodegenArrayTensorProduct(M, N)
    >>> tp.shape
    (4, 5, 3, 2)
    >>> perm1 = CodegenArrayPermuteDims(tp, [2, 3, 1, 0])

    The args ``(M, N)`` have been sorted and the permutation has been
    simplified, the expression is equivalent:

    >>> perm1.expr.args
    (N, M)
    >>> perm1.shape
    (3, 2, 5, 4)
    >>> perm1.permutation
    (2 3)

    The permutation in its array form has been simplified from
    ``[2, 3, 1, 0]`` to ``[0, 1, 3, 2]``, as the arguments of the tensor
    product `M` and `N` have been switched:

    >>> perm1.permutation.array_form
    [0, 1, 3, 2]

    We can nest a second permutation:

    >>> perm2 = CodegenArrayPermuteDims(perm1, [1, 0, 2, 3])
    >>> perm2.shape
    (2, 3, 5, 4)
    >>> perm2.permutation.array_form
    [1, 0, 3, 2]
    """

    def __new__(cls, expr, permutation, nest_permutation=True):
        from sympy.combinatorics import Permutation
        expr = _sympify(expr)
        permutation = Permutation(permutation)
        permutation_size = permutation.size
        expr_rank = get_rank(expr)
        if permutation_size != expr_rank:
            raise ValueError("Permutation size must be the length of the shape of expr")
        if isinstance(expr, CodegenArrayPermuteDims):
            subexpr = expr.expr
            subperm = expr.permutation
            permutation = permutation * subperm
            expr = subexpr
        if isinstance(expr, CodegenArrayContraction):
            expr, permutation = cls._handle_nested_contraction(expr, permutation)
        if isinstance(expr, CodegenArrayTensorProduct):
            expr, permutation = cls._sort_components(expr, permutation)
        if isinstance(expr, (ZeroArray, ZeroMatrix)):
            return ZeroArray(*[expr.shape[i] for i in permutation.array_form])
        plist = permutation.array_form
        if plist == sorted(plist):
            return expr
        obj = Basic.__new__(cls, expr, permutation)
        obj._subranks = [get_rank(expr)]
        shape = expr.shape
        if shape is None:
            obj._shape = None
        else:
            obj._shape = tuple(shape[permutation(i)] for i in range(len(shape)))
        return obj

    @property
    def expr(self):
        return self.args[0]

    @property
    def permutation(self):
        return self.args[1]

    @classmethod
    def _sort_components(cls, expr, permutation):
        # Get the permutation in its image-form:
        perm_image_form = _af_invert(permutation.array_form)
        args = list(expr.args)
        # Starting index global position for every arg:
        cumul = list(accumulate([0] + expr.subranks))
        # Split `perm_image_form` into a list of list corresponding to the indices
        # of every argument:
        perm_image_form_in_components = [perm_image_form[cumul[i]:cumul[i+1]] for i in range(len(args))]
        # Create an index, target-position-key array:
        ps = [(i, sorted(comp)) for i, comp in enumerate(perm_image_form_in_components)]
        # Sort the array according to the target-position-key:
        # In this way, we define a canonical way to sort the arguments according
        # to the permutation.
        ps.sort(key=lambda x: x[1])
        # Read the inverse-permutation (i.e. image-form) of the args:
        perm_args_image_form = [i[0] for i in ps]
        # Apply the args-permutation to the `args`:
        args_sorted = [args[i] for i in perm_args_image_form]
        # Apply the args-permutation to the array-form of the permutation of the axes (of `expr`):
        perm_image_form_sorted_args = [perm_image_form_in_components[i] for i in perm_args_image_form]
        new_permutation = Permutation(_af_invert([j for i in perm_image_form_sorted_args for j in i]))
        return CodegenArrayTensorProduct(*args_sorted), new_permutation

    @classmethod
    def _handle_nested_contraction(cls, expr, permutation):
        if not isinstance(expr, CodegenArrayContraction):
            return expr, permutation
        if not isinstance(expr.expr, CodegenArrayTensorProduct):
            return expr, permutation
        args = expr.expr.args
        subranks = [get_rank(arg) for arg in expr.expr.args]

        contraction_indices = expr.contraction_indices
        contraction_indices_flat = [j for i in contraction_indices for j in i]
        cumul = list(accumulate([0] + subranks))

        # Spread the permutation in its array form across the args in the corresponding
        # tensor-product arguments with free indices:
        permutation_array_blocks_up = []
        image_form = _af_invert(permutation.array_form)
        counter = 0
        for i, e in enumerate(subranks):
            current = []
            for j in range(cumul[i], cumul[i+1]):
                if j in contraction_indices_flat:
                    continue
                current.append(image_form[counter])
                counter += 1
            permutation_array_blocks_up.append(current)

        # Get the map of axis repositioning for every argument of tensor-product:
        index_blocks = [[j for j in range(cumul[i], cumul[i+1])] for i, e in enumerate(expr.subranks)]
        index_blocks_up = expr._push_indices_up(expr.contraction_indices, index_blocks)
        inverse_permutation = permutation**(-1)
        index_blocks_up_permuted = [[inverse_permutation(j) for j in i if j is not None] for i in index_blocks_up]

        # Sorting key is a list of tuple, first element is the index of `args`, second element of
        # the tuple is the sorting key to sort `args` of the tensor product:
        sorting_keys = list(enumerate(index_blocks_up_permuted))
        sorting_keys.sort(key=lambda x: x[1])

        # Now we can get the permutation acting on the args in its image-form:
        new_perm_image_form = [i[0] for i in sorting_keys]
        # Apply the args-level permutation to various elements:
        new_index_blocks = [index_blocks[i] for i in new_perm_image_form]
        new_index_perm_array_form = _af_invert([j for i in new_index_blocks for j in i])
        new_args = [args[i] for i in new_perm_image_form]
        new_contraction_indices = [tuple(new_index_perm_array_form[j] for j in i) for i in contraction_indices]
        new_expr = CodegenArrayContraction(CodegenArrayTensorProduct(*new_args), *new_contraction_indices)
        new_permutation = Permutation(_af_invert([j for i in [permutation_array_blocks_up[k] for k in new_perm_image_form] for j in i]))
        return new_expr, new_permutation

    @classmethod
    def _check_permutation_mapping(cls, expr, permutation):
        subranks = expr.subranks
        index2arg = [i for i, arg in enumerate(expr.args) for j in range(expr.subranks[i])]
        permuted_indices = [permutation(i) for i in range(expr.subrank())]
        new_args = list(expr.args)
        arg_candidate_index = index2arg[permuted_indices[0]]
        current_indices = []
        new_permutation = []
        inserted_arg_cand_indices = set([])
        for i, idx in enumerate(permuted_indices):
            if index2arg[idx] != arg_candidate_index:
                new_permutation.extend(current_indices)
                current_indices = []
                arg_candidate_index = index2arg[idx]
            current_indices.append(idx)
            arg_candidate_rank = subranks[arg_candidate_index]
            if len(current_indices) == arg_candidate_rank:
                new_permutation.extend(sorted(current_indices))
                local_current_indices = [j - min(current_indices) for j in current_indices]
                i1 = index2arg[i]
                new_args[i1] = CodegenArrayPermuteDims(new_args[i1], Permutation(local_current_indices))
                inserted_arg_cand_indices.add(arg_candidate_index)
                current_indices = []
        new_permutation.extend(current_indices)

        # TODO: swap args positions in order to simplify the expression:
        # TODO: this should be in a function
        args_positions = list(range(len(new_args)))
        # Get possible shifts:
        maps = {}
        cumulative_subranks = [0] + list(accumulate(subranks))
        for i in range(0, len(subranks)):
            s = set([index2arg[new_permutation[j]] for j in range(cumulative_subranks[i], cumulative_subranks[i+1])])
            if len(s) != 1:
                continue
            elem = next(iter(s))
            if i != elem:
                maps[i] = elem

        # Find cycles in the map:
        lines = []
        current_line = []
        while maps:
            if len(current_line) == 0:
                k, v = maps.popitem()
                current_line.append(k)
            else:
                k = current_line[-1]
                if k not in maps:
                    current_line = []
                    continue
                v = maps.pop(k)
            if v in current_line:
                lines.append(current_line)
                current_line = []
                continue
            current_line.append(v)
        for line in lines:
            for i, e in enumerate(line):
                args_positions[line[(i + 1) % len(line)]] = e

        # TODO: function in order to permute the args:
        permutation_blocks = [[new_permutation[cumulative_subranks[i] + j] for j in range(e)] for i, e in enumerate(subranks)]
        new_args = [new_args[i] for i in args_positions]
        new_permutation_blocks = [permutation_blocks[i] for i in args_positions]
        new_permutation2 = [j for i in new_permutation_blocks for j in i]
        return CodegenArrayTensorProduct(*new_args), Permutation(new_permutation2)  # **(-1)

    @classmethod
    def _check_if_there_are_closed_cycles(cls, expr, permutation):
        args = list(expr.args)
        subranks = expr.subranks
        cyclic_form = permutation.cyclic_form
        cumulative_subranks = [0] + list(accumulate(subranks))
        cyclic_min = [min(i) for i in cyclic_form]
        cyclic_max = [max(i) for i in cyclic_form]
        cyclic_keep = []
        for i, cycle in enumerate(cyclic_form):
            flag = True
            for j in range(0, len(cumulative_subranks) - 1):
                if cyclic_min[i] >= cumulative_subranks[j] and cyclic_max[i] < cumulative_subranks[j+1]:
                    # Found a sinkable cycle.
                    args[j] = CodegenArrayPermuteDims(args[j], Permutation([[k - cumulative_subranks[j] for k in cyclic_form[i]]]))
                    flag = False
                    break
            if flag:
                cyclic_keep.append(cyclic_form[i])
        return CodegenArrayTensorProduct(*args), Permutation(cyclic_keep, size=permutation.size)

    def nest_permutation(self):
        r"""
        DEPRECATED.
        """
        ret = self._nest_permutation(self.expr, self.permutation)
        if ret is None:
            return self
        return ret

    @classmethod
    def _nest_permutation(cls, expr, permutation):
        if isinstance(expr, CodegenArrayTensorProduct):
            return CodegenArrayPermuteDims(*cls._check_if_there_are_closed_cycles(expr, permutation))
        elif isinstance(expr, CodegenArrayContraction):
            # Invert tree hierarchy: put the contraction above.
            cycles = permutation.cyclic_form
            newcycles = CodegenArrayContraction._convert_outer_indices_to_inner_indices(expr, *cycles)
            newpermutation = Permutation(newcycles)
            new_contr_indices = [tuple(newpermutation(j) for j in i) for i in expr.contraction_indices]
            return CodegenArrayContraction(CodegenArrayPermuteDims(expr.expr, newpermutation), *new_contr_indices)
        elif isinstance(expr, CodegenArrayElementwiseAdd):
            return CodegenArrayElementwiseAdd(*[CodegenArrayPermuteDims(arg, permutation) for arg in expr.args])
        return None

    def as_explicit(self):
        return permutedims(self.expr.as_explicit(), self.permutation)


def nest_permutation(expr):
    if isinstance(expr, CodegenArrayPermuteDims):
        return expr.nest_permutation()
    else:
        return expr


class CodegenArrayDiagonal(_CodegenArrayAbstract):
    r"""
    Class to represent the diagonal operator.

    Explanation
    ===========

    In a 2-dimensional array it returns the diagonal, this looks like the
    operation:

    `A_{ij} \rightarrow A_{ii}`

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
        if isinstance(expr, CodegenArrayElementwiseAdd):
            return CodegenArrayElementwiseAdd(*[CodegenArrayDiagonal(arg, *diagonal_indices) for arg in expr.args])
        if isinstance(expr, CodegenArrayDiagonal):
            return cls._flatten(expr, *diagonal_indices)
        if isinstance(expr, CodegenArrayPermuteDims):
            return cls._handle_nested_permutedims_in_diag(expr, *diagonal_indices)
        shape = expr.shape
        if shape is not None:
            cls._validate(expr, *diagonal_indices)
            # Get new shape:
            positions, shape = cls._get_positions_shape(shape, diagonal_indices)
        else:
            positions = None
        if len(diagonal_indices) == 0:
            return expr
        if isinstance(expr, (ZeroArray, ZeroMatrix)):
            return ZeroArray(*shape)
        obj = Basic.__new__(cls, expr, *diagonal_indices)
        obj._positions = positions
        obj._subranks = _get_subranks(expr)
        obj._shape = shape
        return obj

    @staticmethod
    def _validate(expr, *diagonal_indices):
        # Check that no diagonalization happens on indices with mismatched
        # dimensions:
        shape = expr.shape
        for i in diagonal_indices:
            if len({shape[j] for j in i}) != 1:
                raise ValueError("diagonalizing indices of different dimensions")
            if len(i) <= 1:
                raise ValueError("need at least two axes to diagonalize")

    @staticmethod
    def _remove_trivial_dimensions(shape, *diagonal_indices):
        return [tuple(j for j in i) for i in diagonal_indices if shape[i[0]] != 1]

    @property
    def expr(self):
        return self.args[0]

    @property
    def diagonal_indices(self):
        return self.args[1:]

    @staticmethod
    def _flatten(expr, *outer_diagonal_indices):
        inner_diagonal_indices = expr.diagonal_indices
        all_inner = [j for i in inner_diagonal_indices for j in i]
        all_inner.sort()
        # TODO: add API for total rank and cumulative rank:
        total_rank = _get_subrank(expr)
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

    @classmethod
    def _handle_nested_permutedims_in_diag(cls, expr: CodegenArrayPermuteDims, *diagonal_indices):
        back_diagonal_indices = [[expr.permutation(j) for j in i] for i in diagonal_indices]
        nondiag = [i for i in range(get_rank(expr)) if not any(i in j for j in diagonal_indices)]
        back_nondiag = [expr.permutation(i) for i in nondiag]
        remap = {e: i for i, e in enumerate(sorted(back_nondiag))}
        new_permutation1 = [remap[i] for i in back_nondiag]
        shift = len(new_permutation1)
        diag_block_perm = [i + shift for i in range(len(back_diagonal_indices))]
        new_permutation = new_permutation1 + diag_block_perm
        return CodegenArrayPermuteDims(
            CodegenArrayDiagonal(
                expr.expr,
                *back_diagonal_indices
            ),
            new_permutation
        )

    def _push_indices_down_nonstatic(self, indices):
        transform = lambda x: self._positions[x] if x < len(self._positions) else None
        return _apply_recursively_over_nested_lists(transform, indices)

    def _push_indices_up_nonstatic(self, indices):

        def transform(x):
            for i, e in enumerate(self._positions):
                if (isinstance(e, int) and x == e) or (isinstance(e, tuple) and x in e):
                    return i

        return _apply_recursively_over_nested_lists(transform, indices)

    @classmethod
    def _push_indices_down(cls, diagonal_indices, indices, rank):
        positions, shape = cls._get_positions_shape(range(rank), diagonal_indices)
        transform = lambda x: positions[x] if x < len(positions) else None
        return _apply_recursively_over_nested_lists(transform, indices)

    @classmethod
    def _push_indices_up(cls, diagonal_indices, indices, rank):
        positions, shape = cls._get_positions_shape(range(rank), diagonal_indices)

        def transform(x):
            for i, e in enumerate(positions):
                if (isinstance(e, int) and x == e) or (isinstance(e, tuple) and x in e):
                    return i

        return _apply_recursively_over_nested_lists(transform, indices)

    @classmethod
    def _get_positions_shape(cls, shape, diagonal_indices):
        data1 = tuple((i, shp) for i, shp in enumerate(shape) if not any(i in j for j in diagonal_indices))
        pos1, shp1 = zip(*data1) if data1 else ((), ())
        data2 = tuple((i, shape[i[0]]) for i in diagonal_indices)
        pos2, shp2 = zip(*data2) if data2 else ((), ())
        positions = pos1 + pos2
        shape = shp1 + shp2
        return positions, shape

    def as_explicit(self):
        return tensordiagonal(self.expr.as_explicit(), *self.diagonal_indices)


class ArrayElementwiseApplyFunc(_CodegenArrayAbstract):

    def __new__(cls, function, element):

        if not isinstance(function, Lambda):
            d = Dummy('d')
            function = Lambda(d, function(d))

        obj = _CodegenArrayAbstract.__new__(cls, function, element)
        obj._subranks = _get_subranks(element)
        return obj

    @property
    def function(self):
        return self.args[0]

    @property
    def expr(self):
        return self.args[1]

    @property
    def shape(self):
        return self.expr.shape

    def _get_function_fdiff(self):
        d = Dummy("d")
        function = self.function(d)
        fdiff = function.diff(d)
        if isinstance(fdiff, Function):
            fdiff = type(fdiff)
        else:
            fdiff = Lambda(d, fdiff)
        return fdiff


def get_rank(expr):
    if isinstance(expr, (MatrixExpr, MatrixElement)):
        return 2
    if isinstance(expr, _CodegenArrayAbstract):
        return len(expr.shape)
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
    if hasattr(expr, "shape"):
        return len(expr.shape)
    return 0


def _get_subrank(expr):
    if isinstance(expr, _CodegenArrayAbstract):
        return expr.subrank()
    return get_rank(expr)


def _get_subranks(expr):
    if isinstance(expr, _CodegenArrayAbstract):
        return expr.subranks
    else:
        return [get_rank(expr)]


def _get_mapping_from_subranks(subranks):
    mapping = {}
    counter = 0
    for i, rank in enumerate(subranks):
        for j in range(rank):
            mapping[counter] = (i, j)
            counter += 1
    return mapping


def _get_contraction_links(args, subranks, *contraction_indices):
    mapping = _get_mapping_from_subranks(subranks)
    contraction_tuples = [[mapping[j] for j in i] for i in contraction_indices]
    dlinks = defaultdict(dict)
    for links in contraction_tuples:
        if len(links) == 2:
            (arg1, pos1), (arg2, pos2) = links
            dlinks[arg1][pos1] = (arg2, pos2)
            dlinks[arg2][pos2] = (arg1, pos1)
            continue

    return args, dict(dlinks)


def _sort_contraction_indices(pairing_indices):
    pairing_indices = [Tuple(*sorted(i)) for i in pairing_indices]
    pairing_indices.sort(key=lambda x: min(x))
    return pairing_indices


def _get_diagonal_indices(flattened_indices):
    axes_contraction = defaultdict(list)
    for i, ind in enumerate(flattened_indices):
        if isinstance(ind, (int, Integer)):
            # If the indices is a number, there can be no diagonal operation:
            continue
        axes_contraction[ind].append(i)
    axes_contraction = {k: v for k, v in axes_contraction.items() if len(v) > 1}
    # Put the diagonalized indices at the end:
    ret_indices = [i for i in flattened_indices if i not in axes_contraction]
    diag_indices = list(axes_contraction)
    diag_indices.sort(key=lambda x: flattened_indices.index(x))
    diagonal_indices = [tuple(axes_contraction[i]) for i in diag_indices]
    ret_indices += diag_indices
    ret_indices = tuple(ret_indices)
    return diagonal_indices, ret_indices


def _get_argindex(subindices, ind):
    for i, sind in enumerate(subindices):
        if ind == sind:
            return i
        if isinstance(sind, (set, frozenset)) and ind in sind:
            return i
    raise IndexError("%s not found in %s" % (ind, subindices))


def _codegen_array_parse(expr):
    if isinstance(expr, Sum):
        function = expr.function
        summation_indices = expr.variables
        subexpr, subindices = _codegen_array_parse(function)
        # Check dimensional consistency:
        shape = subexpr.shape
        if shape:
            for ind, istart, iend in expr.limits:
                i = _get_argindex(subindices, ind)
                if istart != 0 or iend+1 != shape[i]:
                    raise ValueError("summation index and array dimension mismatch: %s" % ind)
        contraction_indices = []
        subindices = list(subindices)
        if isinstance(subexpr, CodegenArrayDiagonal):
            diagonal_indices = list(subexpr.diagonal_indices)
            dindices = subindices[-len(diagonal_indices):]
            subindices = subindices[:-len(diagonal_indices)]
            for index in summation_indices:
                if index in dindices:
                    position = dindices.index(index)
                    contraction_indices.append(diagonal_indices[position])
                    diagonal_indices[position] = None
            diagonal_indices = [i for i in diagonal_indices if i is not None]
            for i, ind in enumerate(subindices):
                if ind in summation_indices:
                    pass
            if diagonal_indices:
                subexpr = CodegenArrayDiagonal(subexpr.expr, *diagonal_indices)
            else:
                subexpr = subexpr.expr

        axes_contraction = defaultdict(list)
        for i, ind in enumerate(subindices):
            if ind in summation_indices:
                axes_contraction[ind].append(i)
                subindices[i] = None
        for k, v in axes_contraction.items():
            contraction_indices.append(tuple(v))
        free_indices = [i for i in subindices if i is not None]
        indices_ret = list(free_indices)
        indices_ret.sort(key=lambda x: free_indices.index(x))
        return CodegenArrayContraction(
                subexpr,
                *contraction_indices,
                free_indices=free_indices
            ), tuple(indices_ret)
    if isinstance(expr, Mul):
        args, indices = zip(*[_codegen_array_parse(arg) for arg in expr.args])
        # Check if there are KroneckerDelta objects:
        kronecker_delta_repl = {}
        for arg in args:
            if not isinstance(arg, KroneckerDelta):
                continue
            # Diagonalize two indices:
            i, j = arg.indices
            kindices = set(arg.indices)
            if i in kronecker_delta_repl:
                kindices.update(kronecker_delta_repl[i])
            if j in kronecker_delta_repl:
                kindices.update(kronecker_delta_repl[j])
            kindices = frozenset(kindices)
            for index in kindices:
                kronecker_delta_repl[index] = kindices
        # Remove KroneckerDelta objects, their relations should be handled by
        # CodegenArrayDiagonal:
        newargs = []
        newindices = []
        for arg, loc_indices in zip(args, indices):
            if isinstance(arg, KroneckerDelta):
                continue
            newargs.append(arg)
            newindices.append(loc_indices)
        flattened_indices = [kronecker_delta_repl.get(j, j) for i in newindices for j in i]
        diagonal_indices, ret_indices = _get_diagonal_indices(flattened_indices)
        tp = CodegenArrayTensorProduct(*newargs)
        if diagonal_indices:
            return (CodegenArrayDiagonal(tp, *diagonal_indices), ret_indices)
        else:
            return tp, ret_indices
    if isinstance(expr, MatrixElement):
        indices = expr.args[1:]
        diagonal_indices, ret_indices = _get_diagonal_indices(indices)
        if diagonal_indices:
            return (CodegenArrayDiagonal(expr.args[0], *diagonal_indices), ret_indices)
        else:
            return expr.args[0], ret_indices
    if isinstance(expr, Indexed):
        indices = expr.indices
        diagonal_indices, ret_indices = _get_diagonal_indices(indices)
        if diagonal_indices:
            return (CodegenArrayDiagonal(expr.base, *diagonal_indices), ret_indices)
        else:
            return expr.args[0], ret_indices
    if isinstance(expr, IndexedBase):
        raise NotImplementedError
    if isinstance(expr, KroneckerDelta):
        return expr, expr.indices
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
    return expr, ()


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


def parse_indexed_expression(expr, first_indices=None):
    r"""
    Parse indexed expression into a form useful for code generation.

    Examples
    ========

    >>> from sympy.codegen.array_utils import parse_indexed_expression
    >>> from sympy import MatrixSymbol, Sum, symbols

    >>> i, j, k, d = symbols("i j k d")
    >>> M = MatrixSymbol("M", d, d)
    >>> N = MatrixSymbol("N", d, d)

    Recognize the trace in summation form:

    >>> expr = Sum(M[i, i], (i, 0, d-1))
    >>> parse_indexed_expression(expr)
    CodegenArrayContraction(M, (0, 1))

    Recognize the extraction of the diagonal by using the same index `i` on
    both axes of the matrix:

    >>> expr = M[i, i]
    >>> parse_indexed_expression(expr)
    CodegenArrayDiagonal(M, (0, 1))

    This function can help perform the transformation expressed in two
    different mathematical notations as:

    `\sum_{j=0}^{N-1} A_{i,j} B_{j,k} \Longrightarrow \mathbf{A}\cdot \mathbf{B}`

    Recognize the matrix multiplication in summation form:

    >>> expr = Sum(M[i, j]*N[j, k], (j, 0, d-1))
    >>> parse_indexed_expression(expr)
    CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2))

    Specify that ``k`` has to be the starting index:

    >>> parse_indexed_expression(expr, first_indices=[k])
    CodegenArrayContraction(CodegenArrayTensorProduct(N, M), (0, 3))
    """

    result, indices = _codegen_array_parse(expr)
    if not first_indices:
        return result
    for i in first_indices:
        if i not in indices:
            first_indices.remove(i)
    first_indices.extend([i for i in indices if i not in first_indices])
    permutation = [first_indices.index(i) for i in indices]
    return CodegenArrayPermuteDims(result, permutation)


def _a2m_mul(*args):
    if all(not isinstance(i, _CodegenArrayAbstract) for i in args):
        return MatMul(*args).doit()
    else:
        return CodegenArrayContraction(
            CodegenArrayTensorProduct(*args),
            *[(2*i-1, 2*i) for i in range(1, len(args))]
        )


def _a2m_tensor_product(*args):
    scalars = []
    arrays = []
    for arg in args:
        if isinstance(arg, (MatrixExpr, _ArrayExpr, _CodegenArrayAbstract)):
            arrays.append(arg)
        else:
            scalars.append(arg)
    scalar = Mul.fromiter(scalars)
    if len(arrays) == 0:
        return scalar
    if scalar != 1:
        if isinstance(arrays[0], _CodegenArrayAbstract):
            arrays = [scalar] + arrays
        else:
            arrays[0] *= scalar
    return CodegenArrayTensorProduct(*arrays)


def _a2m_add(*args):
    if all(not isinstance(i, _CodegenArrayAbstract) for i in args):
        return MatAdd(*args).doit()
    else:
        return CodegenArrayElementwiseAdd(*args)


def _a2m_trace(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return CodegenArrayContraction(arg, (0, 1))
    else:
        return Trace(arg)


def _a2m_transpose(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return CodegenArrayPermuteDims(arg, [1, 0])
    else:
        return Transpose(arg).doit()


def _support_function_tp1_recognize(contraction_indices, args):
    subranks = [get_rank(i) for i in args]
    coeff = reduce(lambda x, y: x*y, [arg for arg, srank in zip(args, subranks) if srank == 0], S.One)
    mapping = _get_mapping_from_subranks(subranks)
    new_contraction_indices = list(contraction_indices)
    newargs = args[:]  # make a copy of the list
    removed = [None for i in newargs]
    cumul = list(accumulate([0] + [get_rank(arg) for arg in args]))
    new_perms = [list(range(cumul[i], cumul[i+1])) for i, arg in enumerate(args)]
    for pi, contraction_pair in enumerate(contraction_indices):
        if len(contraction_pair) != 2:
            continue
        i1, i2 = contraction_pair
        a1, e1 = mapping[i1]
        a2, e2 = mapping[i2]
        while removed[a1] is not None:
            a1, e1 = removed[a1]
        while removed[a2] is not None:
            a2, e2 = removed[a2]
        if a1 == a2:
            trace_arg = newargs[a1]
            newargs[a1] = Trace(trace_arg)._normalize()
            new_contraction_indices[pi] = None
            continue
        if not isinstance(newargs[a1], MatrixExpr) or not isinstance(newargs[a2], MatrixExpr):
            continue
        arg1 = newargs[a1]
        arg2 = newargs[a2]
        if (e1 == 1 and e2 == 1) or (e1 == 0 and e2 == 0):
            arg2 = Transpose(arg2)
        if e1 == 1:
            argnew = arg1*arg2
        else:
            argnew = arg2*arg1
        removed[a2] = a1, e1
        new_perms[a1][e1] = new_perms[a2][1 - e2]
        new_perms[a2] = None
        newargs[a1] = argnew
        newargs[a2] = None
        new_contraction_indices[pi] = None
    new_contraction_indices = [i for i in new_contraction_indices if i is not None]
    newargs2 = [arg for arg in newargs if arg is not None]
    if len(newargs2) == 0:
        return coeff
    tp = _a2m_tensor_product(*newargs2)
    tc = CodegenArrayContraction(tp, *new_contraction_indices)
    new_perms2 = CodegenArrayContraction._push_indices_up(contraction_indices, [i for i in new_perms if i is not None])
    permutation = _af_invert([j for i in new_perms2 for j in i if j is not None])
    if permutation == [1, 0] and len(newargs2) == 1:
        return Transpose(newargs2[0]).doit()
    tperm = CodegenArrayPermuteDims(tc, permutation)
    return tperm


def _array_diag2contr_diagmatrix(expr: CodegenArrayDiagonal):
    if isinstance(expr.expr, CodegenArrayTensorProduct):
        args = list(expr.expr.args)
        diag_indices = list(expr.diagonal_indices)
        mapping = _get_mapping_from_subranks([_get_subrank(arg) for arg in args])
        tuple_links = [[mapping[j] for j in i] for i in diag_indices]
        contr_indices = []
        total_rank = get_rank(expr)
        replaced = [False for arg in args]
        for i, (abs_pos, rel_pos) in enumerate(zip(diag_indices, tuple_links)):
            if len(abs_pos) != 2:
                continue
            (pos1_outer, pos1_inner), (pos2_outer, pos2_inner) = rel_pos
            arg1 = args[pos1_outer]
            arg2 = args[pos2_outer]
            if get_rank(arg1) != 2 or get_rank(arg2) != 2:
                if replaced[pos1_outer]:
                    diag_indices[i] = None
                if replaced[pos2_outer]:
                    diag_indices[i] = None
                continue
            pos1_in2 = 1 - pos1_inner
            pos2_in2 = 1 - pos2_inner
            if arg1.shape[pos1_in2] == 1:
                darg1 = DiagMatrix(arg1)
                args.append(darg1)
                contr_indices.append(((pos2_outer, pos2_inner), (len(args)-1, pos1_inner)))
                total_rank += 1
                diag_indices[i] = None
                args[pos1_outer] = OneArray(arg1.shape[pos1_in2])
                replaced[pos1_outer] = True
            elif arg2.shape[pos2_in2] == 1:
                darg2 = DiagMatrix(arg2)
                args.append(darg2)
                contr_indices.append(((pos1_outer, pos1_inner), (len(args)-1, pos2_inner)))
                total_rank += 1
                diag_indices[i] = None
                args[pos2_outer] = OneArray(arg2.shape[pos2_in2])
                replaced[pos2_outer] = True
        diag_indices_new = [i for i in diag_indices if i is not None]
        cumul = list(accumulate([0] + [get_rank(arg) for arg in args]))
        contr_indices2 = [tuple(cumul[a] + b for a, b in i) for i in contr_indices]
        tc = CodegenArrayContraction(
            CodegenArrayTensorProduct(*args), *contr_indices2
        )
        td = CodegenArrayDiagonal(tc, *diag_indices_new)
        return td
    return expr


@singledispatch
def array2matrix(expr):
    return expr


@array2matrix.register(ZeroArray)
def _(expr: ZeroArray):
    if get_rank(expr) == 2:
        return ZeroMatrix(*expr.shape)
    else:
        return expr


@array2matrix.register(CodegenArrayTensorProduct)
def _(expr: CodegenArrayTensorProduct):
    return _a2m_tensor_product(*[array2matrix(arg) for arg in expr.args])


@array2matrix.register(CodegenArrayContraction)
def _(expr: CodegenArrayContraction):
    expr = expr.flatten_contraction_of_diagonal()
    expr = expr.split_multiple_contractions()
    subexpr = expr.expr
    contraction_indices: Tuple[Tuple[int]] = expr.contraction_indices
    if isinstance(subexpr, CodegenArrayTensorProduct):
        newexpr = CodegenArrayContraction(array2matrix(subexpr), *contraction_indices)
        contraction_indices = newexpr.contraction_indices
        if any(i > 2 for i in newexpr.subranks):
            addends = CodegenArrayElementwiseAdd(*[_a2m_tensor_product(*j) for j in itertools.product(*[i.args if isinstance(i, CodegenArrayElementwiseAdd) else [i] for i in expr.expr.args])])
            newexpr = CodegenArrayContraction(addends, *contraction_indices)
        if isinstance(newexpr, CodegenArrayElementwiseAdd):
            ret = array2matrix(newexpr)
            return ret
        assert isinstance(newexpr, CodegenArrayContraction)
        ret = _support_function_tp1_recognize(contraction_indices, list(newexpr.expr.args))
        return ret
    elif not isinstance(subexpr, _CodegenArrayAbstract):
        ret = array2matrix(subexpr)
        if isinstance(ret, MatrixExpr):
            assert expr.contraction_indices == ((0, 1),)
            return _a2m_trace(ret)
        else:
            return CodegenArrayContraction(ret, *expr.contraction_indices)


@array2matrix.register(CodegenArrayDiagonal)
def _(expr: CodegenArrayDiagonal):
    expr2 = array2matrix(expr.expr)
    pexpr = _array_diag2contr_diagmatrix(CodegenArrayDiagonal(expr2, *expr.diagonal_indices))
    if expr == pexpr:
        return expr
    return array2matrix(pexpr)


@array2matrix.register(CodegenArrayPermuteDims)
def _(expr: CodegenArrayPermuteDims):
    if expr.permutation.array_form == [1, 0]:
        return _a2m_transpose(array2matrix(expr.expr))
    elif isinstance(expr.expr, CodegenArrayTensorProduct):
        ranks = expr.expr.subranks
        inv_permutation = expr.permutation**(-1)
        newrange = [inv_permutation(i) for i in range(sum(ranks))]
        newpos = []
        counter = 0
        for rank in ranks:
            newpos.append(newrange[counter:counter+rank])
            counter += rank
        newargs = []
        newperm = []
        scalars = []
        for pos, arg in zip(newpos, expr.expr.args):
            if len(pos) == 0:
                scalars.append(array2matrix(arg))
            elif pos == sorted(pos):
                newargs.append((array2matrix(arg), pos[0]))
                newperm.extend(pos)
            elif len(pos) == 2:
                newargs.append((_a2m_transpose(array2matrix(arg)), pos[0]))
                newperm.extend(reversed(pos))
            else:
                raise NotImplementedError()
        newargs = [i[0] for i in newargs]
        return CodegenArrayPermuteDims(_a2m_tensor_product(*scalars, *newargs), _af_invert(newperm))
    elif isinstance(expr.expr, CodegenArrayContraction):
        mat_mul_lines = array2matrix(expr.expr)
        if not isinstance(mat_mul_lines, CodegenArrayTensorProduct):
            flat_cyclic_form = [j for i in expr.permutation.cyclic_form for j in i]
            expr_shape = get_shape(expr)
            if all(expr_shape[i] == 1 for i in flat_cyclic_form):
                return mat_mul_lines
            return mat_mul_lines
        permutation = Permutation(2*len(mat_mul_lines.args)-1)*expr.permutation
        permuted = [permutation(i) for i in range(2*len(mat_mul_lines.args))]
        args_array = [None for i in mat_mul_lines.args]
        for i in range(len(mat_mul_lines.args)):
            p1 = permuted[2*i]
            p2 = permuted[2*i+1]
            if p1 // 2 != p2 // 2:
                return CodegenArrayPermuteDims(mat_mul_lines, permutation)
            pos = p1 // 2
            if p1 > p2:
                args_array[i] = _a2m_transpose(mat_mul_lines.args[pos])
            else:
                args_array[i] = mat_mul_lines.args[pos]
        return _a2m_tensor_product(*args_array)
    else:
        raise NotImplementedError()


@array2matrix.register(CodegenArrayElementwiseAdd)
def _(expr: CodegenArrayElementwiseAdd):
    addends = [array2matrix(arg) for arg in expr.args]
    return _a2m_add(*addends)


@array2matrix.register(ArrayElementwiseApplyFunc)
def _(expr: ArrayElementwiseApplyFunc):
    subexpr = array2matrix(expr.expr)
    if isinstance(subexpr, MatrixExpr):
        return ElementwiseApplyFunction(expr.function, subexpr)
    else:
        return ArrayElementwiseApplyFunc(expr.function, subexpr)


@singledispatch
def _remove_trivial_dims(expr):
    return expr, []


@_remove_trivial_dims.register(CodegenArrayTensorProduct)
def _(expr: CodegenArrayTensorProduct):
    # Recognize expressions like [x, y] with shape (k, 1, k, 1) as `x*y.T`.
    # The matrix expression has to be equivalent to the tensor product of the
    # matrices, with trivial dimensions (i.e. dim=1) dropped.
    # That is, add contractions over trivial dimensions:

    removed = []
    newargs = []
    cumul = list(accumulate([0] + [get_rank(arg) for arg in expr.args]))
    pending = None
    prev_i = None
    for i, arg in enumerate(expr.args):
        current_range = list(range(cumul[i], cumul[i+1]))
        if isinstance(arg, OneArray):
            removed.extend(current_range)
            continue
        if not isinstance(arg, (MatrixExpr, MatrixCommon)):
            rarg, rem = _remove_trivial_dims(arg)
            removed.extend(rem)
            newargs.append(rarg)
            continue
        elif getattr(arg, "is_Identity", False):
            if arg.shape == (1, 1):
                # Ignore identity matrices of shape (1, 1) - they are equivalent to scalar 1.
                removed.extend(current_range)
                continue
            k = arg.shape[0]
            if pending == k:
                # OK, there is already
                removed.extend(current_range)
                continue
            elif pending is None:
                newargs.append(arg)
                pending = k
                prev_i = i
            elif pending != k:
                pending = k
                prev_i = i
                newargs.append(arg)
        elif arg.shape == (1, 1):
            arg, _ = _remove_trivial_dims(arg)
            # Matrix is equivalent to scalar:
            if len(newargs) == 0:
                newargs.append(arg)
            elif 1 in get_shape(newargs[-1]):
                if newargs[-1].shape[1] == 1:
                    newargs[-1] = newargs[-1]*arg
                else:
                    newargs[-1] = arg*newargs[-1]
                removed.extend(current_range)
            else:
                newargs.append(arg)
        elif 1 in arg.shape:
            k = [i for i in arg.shape if i != 1][0]
            if pending is None:
                pending = k
                prev_i = i
                newargs.append(arg)
            elif pending == k:
                prev = newargs[-1]
                if prev.is_Identity:
                    removed.extend([cumul[prev_i], cumul[prev_i]+1])
                    newargs[-1] = arg
                    prev_i = i
                    continue
                if prev.shape[0] == 1:
                    d1 = cumul[prev_i]
                    prev = _a2m_transpose(prev)
                else:
                    d1 = cumul[prev_i] + 1
                if arg.shape[1] == 1:
                    d2 = cumul[i] + 1
                    arg = _a2m_transpose(arg)
                else:
                    d2 = cumul[i]
                newargs[-1] = prev*arg
                pending = None
                removed.extend([d1, d2])
            else:
                newargs.append(arg)
                pending = k
                prev_i = i
        else:
            newargs.append(arg)
            pending = None
    return _a2m_tensor_product(*newargs), sorted(removed)


@_remove_trivial_dims.register(CodegenArrayElementwiseAdd)
def _(expr: CodegenArrayElementwiseAdd):
    rec = [_remove_trivial_dims(arg) for arg in expr.args]
    newargs, removed = zip(*rec)
    if len(set(map(tuple, removed))) != 1:
        return expr, []
    return _a2m_add(*newargs), removed[0]


@_remove_trivial_dims.register(CodegenArrayPermuteDims)
def _(expr: CodegenArrayPermuteDims):
    subexpr, subremoved = _remove_trivial_dims(expr.expr)
    p = expr.permutation.array_form
    pinv = _af_invert(expr.permutation.array_form)
    shift = list(accumulate([1 if i in subremoved else 0 for i in range(len(p))]))
    premoved = [pinv[i] for i in subremoved]
    p2 = [e - shift[e] for i, e in enumerate(p) if e not in subremoved]
    # TODO: check if subremoved should be permuted as well...
    newexpr = CodegenArrayPermuteDims(subexpr, p2)
    if newexpr != expr:
        newexpr = array2matrix(newexpr)
    return newexpr, sorted(premoved)


@_remove_trivial_dims.register(CodegenArrayContraction)
def _(expr: CodegenArrayContraction):
    newexpr, removed = _remove_trivial_dims(expr.expr)
    new_contraction_indices = [tuple(j for j in i if j not in removed) for i in expr.contraction_indices]
    # Remove possible empty tuples "()":
    new_contraction_indices = [i for i in new_contraction_indices if i]
    return CodegenArrayContraction(newexpr, *new_contraction_indices), removed


@_remove_trivial_dims.register(CodegenArrayDiagonal)
def _(expr: CodegenArrayDiagonal):
    newexpr, removed = _remove_trivial_dims(expr.expr)
    new_diag_indices = [tuple(j for j in i if j not in removed) for i in expr.diagonal_indices]
    return CodegenArrayDiagonal(newexpr, *new_diag_indices), removed


@_remove_trivial_dims.register(ElementwiseApplyFunction)
def _(expr: ElementwiseApplyFunction):
    subexpr, removed = _remove_trivial_dims(expr.expr)
    if subexpr.shape == (1, 1):
        # TODO: move this to ElementwiseApplyFunction
        return expr.function(subexpr), removed + [0, 1]
    return ElementwiseApplyFunction(expr.function, subexpr)


@_remove_trivial_dims.register(ArrayElementwiseApplyFunc)
def _(expr: ArrayElementwiseApplyFunc):
    subexpr, removed = _remove_trivial_dims(expr.expr)
    return ArrayElementwiseApplyFunc(expr.function, subexpr), removed


def recognize_matrix_expression(expr):
    r"""
    Recognize matrix expressions in codegen objects.

    If more than one matrix multiplication line have been detected, return a
    list with the matrix expressions.

    Examples
    ========

    >>> from sympy import MatrixSymbol, Sum
    >>> from sympy.abc import i, j, k, l, N
    >>> from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
    >>> from sympy.codegen.array_utils import recognize_matrix_expression, parse_indexed_expression, parse_matrix_expression
    >>> A = MatrixSymbol("A", N, N)
    >>> B = MatrixSymbol("B", N, N)
    >>> C = MatrixSymbol("C", N, N)
    >>> D = MatrixSymbol("D", N, N)

    >>> expr = Sum(A[i, j]*B[j, k], (j, 0, N-1))
    >>> cg = parse_indexed_expression(expr)
    >>> recognize_matrix_expression(cg)
    A*B
    >>> cg = parse_indexed_expression(expr, first_indices=[k])
    >>> recognize_matrix_expression(cg)
    B.T*A.T

    Transposition is detected:

    >>> expr = Sum(A[j, i]*B[j, k], (j, 0, N-1))
    >>> cg = parse_indexed_expression(expr)
    >>> recognize_matrix_expression(cg)
    A.T*B
    >>> cg = parse_indexed_expression(expr, first_indices=[k])
    >>> recognize_matrix_expression(cg)
    B.T*A

    Detect the trace:

    >>> expr = Sum(A[i, i], (i, 0, N-1))
    >>> cg = parse_indexed_expression(expr)
    >>> recognize_matrix_expression(cg)
    Trace(A)

    Recognize some more complex traces:

    >>> expr = Sum(A[i, j]*B[j, i], (i, 0, N-1), (j, 0, N-1))
    >>> cg = parse_indexed_expression(expr)
    >>> recognize_matrix_expression(cg)
    Trace(A*B)

    More complicated expressions:

    >>> expr = Sum(A[i, j]*B[k, j]*A[l, k], (j, 0, N-1), (k, 0, N-1))
    >>> cg = parse_indexed_expression(expr)
    >>> recognize_matrix_expression(cg)
    A*B.T*A.T

    Expressions constructed from matrix expressions do not contain literal
    indices, the positions of free indices are returned instead:

    >>> expr = A*B
    >>> cg = parse_matrix_expression(expr)
    >>> recognize_matrix_expression(cg)
    A*B

    If more than one line of matrix multiplications is detected, return
    separate matrix multiplication factors embedded in a tensor product object:

    >>> cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B, C, D), (1, 2), (5, 6))
    >>> recognize_matrix_expression(cg)
    CodegenArrayTensorProduct(A*B, C*D)

    The two lines have free indices at axes 0, 3 and 4, 7, respectively.
    """
    rec = array2matrix(expr)
    rec, removed = _remove_trivial_dims(rec)
    return rec


def _apply_recursively_over_nested_lists(func, arr):
    if isinstance(arr, (tuple, list, Tuple)):
        return tuple(_apply_recursively_over_nested_lists(func, i) for i in arr)
    elif isinstance(arr, Tuple):
        return Tuple.fromiter(_apply_recursively_over_nested_lists(func, i) for i in arr)
    else:
        return func(arr)


def _build_push_indices_up_func_transformation(flattened_contraction_indices):
    shifts = {0: 0}
    i = 0
    cumulative = 0
    while i < len(flattened_contraction_indices):
        j = 1
        while i+j < len(flattened_contraction_indices):
            if flattened_contraction_indices[i] + j != flattened_contraction_indices[i+j]:
                break
            j += 1
        cumulative += j
        shifts[flattened_contraction_indices[i]] = cumulative
        i += j
    shift_keys = sorted(shifts.keys())

    def func(idx):
        return shifts[shift_keys[bisect.bisect_right(shift_keys, idx)-1]]

    def transform(j):
        if j in flattened_contraction_indices:
            return None
        else:
            return j - func(j)

    return transform


def _build_push_indices_down_func_transformation(flattened_contraction_indices):
    N = flattened_contraction_indices[-1]+2

    shifts = [i for i in range(N) if i not in flattened_contraction_indices]

    def transform(j):
        if j < len(shifts):
            return shifts[j]
        else:
            return j + shifts[-1] - len(shifts) + 1

    return transform
