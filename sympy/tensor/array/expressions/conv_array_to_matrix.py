import itertools
from collections import defaultdict, Counter
from typing import Tuple, Union, FrozenSet, Dict, List, Optional
from functools import singledispatch
from itertools import accumulate

from sympy import Trace, MatrixExpr, Transpose, DiagMatrix, Mul, ZeroMatrix, hadamard_product, S
from sympy.combinatorics.permutations import _af_invert, Permutation
from sympy.matrices.common import MatrixCommon
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.tensor.array.expressions.array_expressions import PermuteDims, ArrayDiagonal, \
    ArrayTensorProduct, OneArray, get_rank, _get_subrank, ZeroArray, ArrayContraction, \
    ArrayAdd, _CodegenArrayAbstract, get_shape, ArrayElementwiseApplyFunc, _ArrayExpr, _EditArrayContraction, _ArgE
from sympy.tensor.array.expressions.utils import _get_mapping_from_subranks


def _get_candidate_for_matmul_from_contraction(scan_indices: List[Optional[int]], remaining_args: List[_ArgE]) -> Tuple[Optional[_ArgE], bool, int]:

    scan_indices = [i for i in scan_indices if i is not None]
    if len(scan_indices) == 0:
        return None, False, -1

    transpose: bool = False
    candidate: Optional[_ArgE] = None
    candidate_index: int = -1
    for arg_with_ind2 in remaining_args:
        if not isinstance(arg_with_ind2.element, MatrixExpr):
            continue
        for index in scan_indices:
            if candidate_index != -1 and candidate_index != index:
                # A candidate index has already been selected, check
                # repetitions only for that index:
                continue
            if index in arg_with_ind2.indices:
                if set(arg_with_ind2.indices) == {index}:
                    # Index repeated twice in arg_with_ind2
                    candidate = None
                    break
                if candidate is None:
                    candidate = arg_with_ind2
                    candidate_index = index
                    transpose = (index == arg_with_ind2.indices[1])
                else:
                    # Index repeated more than twice, break
                    candidate = None
                    break
    return candidate, transpose, candidate_index


def _insert_candidate_into_editor(editor: _EditArrayContraction, arg_with_ind: _ArgE, candidate: _ArgE, transpose1: bool, transpose2: bool):
    other = candidate.element
    other_index: int
    if transpose2:
        other = Transpose(other)
        other_index = candidate.indices[0]
    else:
        other_index = candidate.indices[1]
    new_element = (Transpose(arg_with_ind.element) if transpose1 else arg_with_ind.element) * other
    editor.args_with_ind.remove(candidate)
    new_arge = _ArgE(new_element)
    return new_arge, other_index


def _support_function_tp1_recognize(contraction_indices, args):
    if len(contraction_indices) == 0:
        return _a2m_tensor_product(*args)

    ac = ArrayContraction(ArrayTensorProduct(*args), *contraction_indices)
    editor = _EditArrayContraction(ac)
    editor.track_permutation_start()

    while True:
        flag_stop: bool = True
        for i, arg_with_ind in enumerate(editor.args_with_ind):
            if not isinstance(arg_with_ind.element, MatrixExpr):
                continue

            first_index = arg_with_ind.indices[0]
            second_index = arg_with_ind.indices[1]

            first_frequency = editor.count_args_with_index(first_index)
            second_frequency = editor.count_args_with_index(second_index)

            if first_index is not None and first_frequency == 1 and first_index == second_index:
                flag_stop = False
                arg_with_ind.element = Trace(arg_with_ind.element)._normalize()
                arg_with_ind.indices = []
                break

            scan_indices = []
            if first_frequency == 2:
                scan_indices.append(first_index)
            if second_frequency == 2:
                scan_indices.append(second_index)

            candidate, transpose, found_index = _get_candidate_for_matmul_from_contraction(scan_indices, editor.args_with_ind[i+1:])
            if candidate is not None:
                flag_stop = False
                editor.track_permutation_merge(arg_with_ind, candidate)
                transpose1 = found_index == first_index
                new_arge, other_index = _insert_candidate_into_editor(editor, arg_with_ind, candidate, transpose1, transpose)
                if found_index == first_index:
                    new_arge.indices = [second_index, other_index]
                else:
                    new_arge.indices = [first_index, other_index]
                set_indices = set(new_arge.indices)
                if len(set_indices) == 1 and set_indices != {None}:
                    # This is a trace:
                    new_arge.element = Trace(new_arge.element)._normalize()
                    new_arge.indices = []
                editor.args_with_ind[i] = new_arge
                # TODO: is this break necessary?
                break

        if flag_stop:
            break

    editor.refresh_indices()
    return editor.to_array_contraction()


@singledispatch
def _array2matrix(expr):
    return expr


@_array2matrix.register(ZeroArray)
def _(expr: ZeroArray):
    if get_rank(expr) == 2:
        return ZeroMatrix(*expr.shape)
    else:
        return expr


@_array2matrix.register(ArrayTensorProduct)
def _(expr: ArrayTensorProduct):
    return _a2m_tensor_product(*[_array2matrix(arg) for arg in expr.args])


@_array2matrix.register(ArrayContraction)
def _(expr: ArrayContraction):
    expr = expr.flatten_contraction_of_diagonal()
    expr = expr.split_multiple_contractions()
    expr = identify_hadamard_products(expr)
    if not isinstance(expr, ArrayContraction):
        return _array2matrix(expr)
    subexpr = expr.expr
    contraction_indices: Tuple[Tuple[int]] = expr.contraction_indices
    if isinstance(subexpr, ArrayTensorProduct):
        newexpr = ArrayContraction(_array2matrix(subexpr), *contraction_indices)
        contraction_indices = newexpr.contraction_indices
        if any(i > 2 for i in newexpr.subranks):
            addends = ArrayAdd(*[_a2m_tensor_product(*j) for j in itertools.product(*[i.args if isinstance(i,
                                                                                                                             ArrayAdd) else [i] for i in expr.expr.args])])
            newexpr = ArrayContraction(addends, *contraction_indices)
        if isinstance(newexpr, ArrayAdd):
            ret = _array2matrix(newexpr)
            return ret
        assert isinstance(newexpr, ArrayContraction)
        ret = _support_function_tp1_recognize(contraction_indices, list(newexpr.expr.args))
        return ret
    elif not isinstance(subexpr, _CodegenArrayAbstract):
        ret = _array2matrix(subexpr)
        if isinstance(ret, MatrixExpr):
            assert expr.contraction_indices == ((0, 1),)
            return _a2m_trace(ret)
        else:
            return ArrayContraction(ret, *expr.contraction_indices)


@_array2matrix.register(ArrayDiagonal)
def _(expr: ArrayDiagonal):
    pexpr = ArrayDiagonal(_array2matrix(expr.expr), *expr.diagonal_indices)
    pexpr = identify_hadamard_products(pexpr)
    if isinstance(pexpr, ArrayDiagonal):
        pexpr = _array_diag2contr_diagmatrix(pexpr)
    if expr == pexpr:
        return expr
    return _array2matrix(pexpr)


@_array2matrix.register(PermuteDims)
def _(expr: PermuteDims):
    if expr.permutation.array_form == [1, 0]:
        return _a2m_transpose(_array2matrix(expr.expr))
    elif isinstance(expr.expr, ArrayTensorProduct):
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
                scalars.append(_array2matrix(arg))
            elif pos == sorted(pos):
                newargs.append((_array2matrix(arg), pos[0]))
                newperm.extend(pos)
            elif len(pos) == 2:
                newargs.append((_a2m_transpose(_array2matrix(arg)), pos[0]))
                newperm.extend(reversed(pos))
            else:
                raise NotImplementedError()
        newargs = [i[0] for i in newargs]
        return PermuteDims(_a2m_tensor_product(*scalars, *newargs), _af_invert(newperm))
    elif isinstance(expr.expr, ArrayContraction):
        mat_mul_lines = _array2matrix(expr.expr)
        if not isinstance(mat_mul_lines, ArrayTensorProduct):
            flat_cyclic_form = [j for i in expr.permutation.cyclic_form for j in i]
            expr_shape = get_shape(expr)
            if all(expr_shape[i] == 1 for i in flat_cyclic_form):
                return mat_mul_lines
            return mat_mul_lines
        # TODO: this assumes that all arguments are matrices, it may not be the case:
        permutation = Permutation(2*len(mat_mul_lines.args)-1)*expr.permutation
        permuted = [permutation(i) for i in range(2*len(mat_mul_lines.args))]
        args_array = [None for i in mat_mul_lines.args]
        for i in range(len(mat_mul_lines.args)):
            p1 = permuted[2*i]
            p2 = permuted[2*i+1]
            if p1 // 2 != p2 // 2:
                return PermuteDims(mat_mul_lines, permutation)
            pos = p1 // 2
            if p1 > p2:
                args_array[i] = _a2m_transpose(mat_mul_lines.args[pos])
            else:
                args_array[i] = mat_mul_lines.args[pos]
        return _a2m_tensor_product(*args_array)
    else:
        return expr


@_array2matrix.register(ArrayAdd)
def _(expr: ArrayAdd):
    addends = [_array2matrix(arg) for arg in expr.args]
    return _a2m_add(*addends)


@_array2matrix.register(ArrayElementwiseApplyFunc)
def _(expr: ArrayElementwiseApplyFunc):
    subexpr = _array2matrix(expr.expr)
    if isinstance(subexpr, MatrixExpr):
        return ElementwiseApplyFunction(expr.function, subexpr)
    else:
        return ArrayElementwiseApplyFunc(expr.function, subexpr)


@singledispatch
def _remove_trivial_dims(expr):
    return expr, []


@_remove_trivial_dims.register(ArrayTensorProduct)
def _(expr: ArrayTensorProduct):
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
            else:
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


@_remove_trivial_dims.register(ArrayAdd)
def _(expr: ArrayAdd):
    rec = [_remove_trivial_dims(arg) for arg in expr.args]
    newargs, removed = zip(*rec)
    if len(set(map(tuple, removed))) != 1:
        return expr, []
    return _a2m_add(*newargs), removed[0]


@_remove_trivial_dims.register(PermuteDims)
def _(expr: PermuteDims):
    subexpr, subremoved = _remove_trivial_dims(expr.expr)
    p = expr.permutation.array_form
    pinv = _af_invert(expr.permutation.array_form)
    shift = list(accumulate([1 if i in subremoved else 0 for i in range(len(p))]))
    premoved = [pinv[i] for i in subremoved]
    p2 = [e - shift[e] for i, e in enumerate(p) if e not in subremoved]
    # TODO: check if subremoved should be permuted as well...
    newexpr = PermuteDims(subexpr, p2)
    if newexpr != expr:
        newexpr = _array2matrix(newexpr)
    return newexpr, sorted(premoved)


@_remove_trivial_dims.register(ArrayContraction)
def _(expr: ArrayContraction):
    newexpr, removed = _remove_trivial_dims(expr.expr)
    shifts = list(accumulate([1 if i in removed else 0 for i in range(get_rank(expr.expr))]))
    new_contraction_indices = [tuple(j for j in i if j not in removed) for i in expr.contraction_indices]
    # Remove possible empty tuples "()":
    new_contraction_indices = [i for i in new_contraction_indices if len(i) > 0]
    contraction_indices_flat = [j for i in expr.contraction_indices for j in i]
    removed = [i for i in removed if i not in contraction_indices_flat]
    new_contraction_indices = [tuple(j - shifts[j] for j in i) for i in new_contraction_indices]
    # Shift removed:
    removed = ArrayContraction._push_indices_up(expr.contraction_indices, removed)
    return ArrayContraction(newexpr, *new_contraction_indices), list(removed)


@_remove_trivial_dims.register(ArrayDiagonal)
def _(expr: ArrayDiagonal):
    newexpr, removed = _remove_trivial_dims(expr.expr)
    shifts = list(accumulate([0] + [1 if i in removed else 0 for i in range(get_rank(expr.expr))]))
    new_diag_indices = [tuple(j for j in i if j not in removed) for i in expr.diagonal_indices]
    new_diag_indices = [tuple(j - shifts[j] for j in i) for i in new_diag_indices]
    rank = get_rank(expr.expr)
    removed = ArrayDiagonal._push_indices_up(expr.diagonal_indices, removed, rank)
    removed = sorted({i for i in removed})
    # If there are single axes to diagonalize remaining, it means that their
    # corresponding dimension has been removed, they no longer need diagonalization:
    new_diag_indices = [i for i in new_diag_indices if len(i) > 1]
    return ArrayDiagonal(newexpr, *new_diag_indices), removed


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


def convert_array_to_matrix(expr):
    r"""
    Recognize matrix expressions in codegen objects.

    If more than one matrix multiplication line have been detected, return a
    list with the matrix expressions.

    Examples
    ========

    >>> from sympy.tensor.array.expressions.conv_indexed_to_array import convert_indexed_to_array
    >>> from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct
    >>> from sympy import MatrixSymbol, Sum
    >>> from sympy.abc import i, j, k, l, N
    >>> from sympy.tensor.array.expressions.array_expressions import ArrayContraction
    >>> from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
    >>> from sympy.tensor.array.expressions.conv_array_to_matrix import convert_array_to_matrix
    >>> A = MatrixSymbol("A", N, N)
    >>> B = MatrixSymbol("B", N, N)
    >>> C = MatrixSymbol("C", N, N)
    >>> D = MatrixSymbol("D", N, N)

    >>> expr = Sum(A[i, j]*B[j, k], (j, 0, N-1))
    >>> cg = convert_indexed_to_array(expr)
    >>> convert_array_to_matrix(cg)
    A*B
    >>> cg = convert_indexed_to_array(expr, first_indices=[k])
    >>> convert_array_to_matrix(cg)
    B.T*A.T

    Transposition is detected:

    >>> expr = Sum(A[j, i]*B[j, k], (j, 0, N-1))
    >>> cg = convert_indexed_to_array(expr)
    >>> convert_array_to_matrix(cg)
    A.T*B
    >>> cg = convert_indexed_to_array(expr, first_indices=[k])
    >>> convert_array_to_matrix(cg)
    B.T*A

    Detect the trace:

    >>> expr = Sum(A[i, i], (i, 0, N-1))
    >>> cg = convert_indexed_to_array(expr)
    >>> convert_array_to_matrix(cg)
    Trace(A)

    Recognize some more complex traces:

    >>> expr = Sum(A[i, j]*B[j, i], (i, 0, N-1), (j, 0, N-1))
    >>> cg = convert_indexed_to_array(expr)
    >>> convert_array_to_matrix(cg)
    Trace(A*B)

    More complicated expressions:

    >>> expr = Sum(A[i, j]*B[k, j]*A[l, k], (j, 0, N-1), (k, 0, N-1))
    >>> cg = convert_indexed_to_array(expr)
    >>> convert_array_to_matrix(cg)
    A*B.T*A.T

    Expressions constructed from matrix expressions do not contain literal
    indices, the positions of free indices are returned instead:

    >>> expr = A*B
    >>> cg = convert_matrix_to_array(expr)
    >>> convert_array_to_matrix(cg)
    A*B

    If more than one line of matrix multiplications is detected, return
    separate matrix multiplication factors embedded in a tensor product object:

    >>> cg = ArrayContraction(ArrayTensorProduct(A, B, C, D), (1, 2), (5, 6))
    >>> convert_array_to_matrix(cg)
    ArrayTensorProduct(A*B, C*D)

    The two lines have free indices at axes 0, 3 and 4, 7, respectively.
    """
    rec = _array2matrix(expr)
    rec, removed = _remove_trivial_dims(rec)
    return rec


def _array_diag2contr_diagmatrix(expr: ArrayDiagonal):
    if isinstance(expr.expr, ArrayTensorProduct):
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
        tc = ArrayContraction(
            ArrayTensorProduct(*args), *contr_indices2
        )
        td = ArrayDiagonal(tc, *diag_indices_new)
        return td
    return expr


def _a2m_mul(*args):
    if all(not isinstance(i, _CodegenArrayAbstract) for i in args):
        from sympy import MatMul
        return MatMul(*args).doit()
    else:
        return ArrayContraction(
            ArrayTensorProduct(*args),
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
    return ArrayTensorProduct(*arrays)


def _a2m_add(*args):
    if all(not isinstance(i, _CodegenArrayAbstract) for i in args):
        from sympy import MatAdd
        return MatAdd(*args).doit()
    else:
        return ArrayAdd(*args)


def _a2m_trace(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return ArrayContraction(arg, (0, 1))
    else:
        from sympy import Trace
        return Trace(arg)


def _a2m_transpose(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return PermuteDims(arg, [1, 0])
    else:
        from sympy import Transpose
        return Transpose(arg).doit()


def identify_hadamard_products(expr: Union[ArrayContraction, ArrayDiagonal]):
    mapping = _get_mapping_from_subranks(expr.subranks)

    editor: _EditArrayContraction
    if isinstance(expr, ArrayContraction):
        editor = _EditArrayContraction(expr)
    elif isinstance(expr, ArrayDiagonal):
        if isinstance(expr.expr, ArrayContraction):
            editor = _EditArrayContraction(expr.expr)
            diagonalized = ArrayContraction._push_indices_down(expr.expr.contraction_indices, expr.diagonal_indices)
        elif isinstance(expr.expr, ArrayTensorProduct):
            editor = _EditArrayContraction(None)
            editor.args_with_ind = [_ArgE(arg) for i, arg in enumerate(expr.expr.args)]
            diagonalized = expr.diagonal_indices
        else:
            return expr

        # Trick: add diagonalized indices as negative indices into the editor object:
        for i, e in enumerate(diagonalized):
            for j in e:
                arg_pos, rel_pos = mapping[j]
                editor.args_with_ind[arg_pos].indices[rel_pos] = -1 - i

    map_contr_to_args: Dict[FrozenSet, List[_ArgE]] = defaultdict(list)
    map_ind_to_inds = defaultdict(int)
    for arg_with_ind in editor.args_with_ind:
        for ind in arg_with_ind.indices:
            map_ind_to_inds[ind] += 1
        if None in arg_with_ind.indices:
            continue
        map_contr_to_args[frozenset(arg_with_ind.indices)].append(arg_with_ind)

    k: FrozenSet[int]
    v: List[_ArgE]
    for k, v in map_contr_to_args.items():
        make_trace: bool = False
        if len(k) == 1 and next(iter(k)) >= 0 and sum([next(iter(k)) in i for i in map_contr_to_args]) == 1:
            # This is a trace: the arguments are fully contracted with only one
            # index, and the index isn't used anywhere else:
            make_trace = True
            first_element = S.One
        elif len(k) != 2:
            # Hadamard product only defined for matrices:
            continue
        if len(v) == 1:
            # Hadamard product with a single argument makes no sense:
            continue
        for ind in k:
            if map_ind_to_inds[ind] <= 2:
                # There is no other contraction, skip:
                continue

        def check_transpose(x):
            x = [i if i >= 0 else -1-i for i in x]
            return x == sorted(x)

        # Check if expression is a trace:
        if all([map_ind_to_inds[j] == len(v) and j >= 0 for j in k]) and all([j >= 0 for j in k]):
            # This is a trace
            make_trace = True
            first_element = v[0].element
            if not check_transpose(v[0].indices):
                first_element = first_element.T
            hadamard_factors = v[1:]
        else:
            hadamard_factors = v

        # This is a Hadamard product:

        hp = hadamard_product(*[i.element if check_transpose(i.indices) else Transpose(i.element) for i in hadamard_factors])
        hp_indices = v[0].indices
        if not check_transpose(hadamard_factors[0].indices):
            hp_indices = list(reversed(hp_indices))
        if make_trace:
            hp = Trace(first_element*hp.T)._normalize()
            hp_indices = []
        editor.insert_after(v[0], _ArgE(hp, hp_indices))
        for i in v:
            editor.args_with_ind.remove(i)

    # Count the ranks of the arguments:
    counter = 0
    # Create a collector for the new diagonal indices:
    diag_indices = defaultdict(list)

    count_index_freq = Counter()
    for arg_with_ind in editor.args_with_ind:
        count_index_freq.update(Counter(arg_with_ind.indices))

    free_index_count = count_index_freq[None]

    # Construct the inverse permutation:
    inv_perm1 = []
    inv_perm2 = []
    # Keep track of which diagonal indices have already been processed:
    done = set([])

    # Counter for the diagonal indices:
    counter4 = 0

    for arg_with_ind in editor.args_with_ind:
        # If some diagonalization axes have been removed, they should be
        # permuted in order to keep the permutation.
        # Add permutation here
        counter2 = 0  # counter for the indices
        for i in arg_with_ind.indices:
            if i is None:
                inv_perm1.append(counter4)
                counter2 += 1
                counter4 += 1
                continue
            if i >= 0:
                continue
            # Reconstruct the diagonal indices:
            diag_indices[-1 - i].append(counter + counter2)
            if count_index_freq[i] == 1 and i not in done:
                inv_perm1.append(free_index_count - 1 - i)
                done.add(i)
            elif i not in done:
                inv_perm2.append(free_index_count - 1 - i)
                done.add(i)
            counter2 += 1
        # Remove negative indices to restore a proper editor object:
        arg_with_ind.indices = [i if i is not None and i >= 0 else None for i in arg_with_ind.indices]
        counter += len([i for i in arg_with_ind.indices if i is None or i < 0])

    inverse_permutation = inv_perm1 + inv_perm2
    permutation = _af_invert(inverse_permutation)

    if isinstance(expr, ArrayContraction):
        return editor.to_array_contraction()
    else:
        # Get the diagonal indices after the detection of HadamardProduct in the expression:
        diag_indices_filtered = [tuple(v) for v in diag_indices.values() if len(v) > 1]

        expr1 = editor.to_array_contraction()
        expr2 = ArrayDiagonal(expr1, *diag_indices_filtered)
        expr3 = PermuteDims(expr2, permutation)
        return expr3
