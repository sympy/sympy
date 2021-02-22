import itertools
from typing import Tuple
from functools import reduce, singledispatch
from itertools import accumulate

from sympy import S, Trace, MatrixExpr, Transpose, DiagMatrix, Mul, ZeroMatrix
from sympy.combinatorics.permutations import _af_invert, Permutation
from sympy.matrices.common import MatrixCommon
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.tensor.array.expressions.array_expressions import CodegenArrayPermuteDims, CodegenArrayDiagonal, \
    CodegenArrayTensorProduct, OneArray, get_rank, _get_subrank, ZeroArray, CodegenArrayContraction, \
    CodegenArrayElementwiseAdd, _CodegenArrayAbstract, get_shape, ArrayElementwiseApplyFunc, _ArrayExpr
from sympy.tensor.array.expressions.utils import _get_mapping_from_subranks


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
            addends = CodegenArrayElementwiseAdd(*[_a2m_tensor_product(*j) for j in itertools.product(*[i.args if isinstance(i,
                                                                                                                             CodegenArrayElementwiseAdd) else [i] for i in expr.expr.args])])
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

    >>> from sympy.tensor.array.expressions.conv_indexed_to_array import parse_indexed_expression
    >>> from sympy.tensor.array.expressions.array_expressions import CodegenArrayTensorProduct
    >>> from sympy import MatrixSymbol, Sum
    >>> from sympy.abc import i, j, k, l, N
    >>> from sympy.tensor.array.expressions.array_expressions import CodegenArrayContraction
    >>> from sympy.tensor.array.expressions.conv_matrix_to_array import parse_matrix_expression
    >>> from sympy.tensor.array.expressions.conv_array_to_matrix import recognize_matrix_expression
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


def _a2m_mul(*args):
    if all(not isinstance(i, _CodegenArrayAbstract) for i in args):
        from sympy import MatMul
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
        from sympy import MatAdd
        return MatAdd(*args).doit()
    else:
        return CodegenArrayElementwiseAdd(*args)


def _a2m_trace(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return CodegenArrayContraction(arg, (0, 1))
    else:
        from sympy import Trace
        return Trace(arg)


def _a2m_transpose(arg):
    if isinstance(arg, _CodegenArrayAbstract):
        return CodegenArrayPermuteDims(arg, [1, 0])
    else:
        from sympy import Transpose
        return Transpose(arg).doit()
