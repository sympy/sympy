"""Module for differentiation using CSE."""

from sympy import cse, Matrix, SparseMatrix, Derivative, MatrixBase
from collections import Counter
import re
from sympy import nan, S


def _process_cse(repl, reduced):
    """
    The `process_cse` function is designed to postprocess the output of a common subexpression elimination (CSE)
    operation. Specifically, it removes any CSE replacement symbols from the arguments of `Derivative` terms in
    the expression. This is necessary to ensure that the forward Jacobian function correctly handles derivative terms.

    Parameters
    ==========

    repl : list of (Symbol, expression) pairs
         Replacement Symbols and relative Common Subexpressions that have been replaced during a CSE operation.

    reduced : list of SymPy expressions
         The reduced expressions with all the replacements from the repl list above.

    Returns
    =======

    p_repl : list of (Symbol, expression) pairs
        Processed replacement list, in the same format of the 'repl' input list.

    p_reduced : list of SymPy expressions
         Processed reduced list, in the same format of the 'reduced' input list.

    """

    def _traverse(node, repl_dict):
        if isinstance(node, Derivative):
            return _replace_all(node, repl_dict)
        if not node.args:
            return node
        new_args = [_traverse(arg, repl_dict) for arg in node.args]
        return node.func(*new_args)

    def _replace_all(node, repl_dict):
        result = node
        while True:
            fs = result.free_symbols
            sl_dict = {k: repl_dict[k] for k in fs if k in repl_dict}
            if not sl_dict:
                break
            result = result.xreplace(sl_dict)
        return result

    repl_dict = dict(repl)
    p_repl = [(rep_sym, _traverse(sub_exp, repl_dict)) for rep_sym, sub_exp in repl]
    p_reduced = [Matrix([_traverse(exp, repl_dict) for exp in red_exp]) for red_exp in reduced]

    return p_repl, p_reduced


def _check_nan(A, nan_idx, n):
    """
    Check if any element of Matrix A multiplied by zero doesn't return zero.
    For each element which doesn't return zero, update the nan_idx set
    to store either the row index (n=0) or column index (n=1).
    """

    zero = S.Zero

    if n == 0:
        for (i, k), value in A.items():
            if zero * value != zero:
                nan_idx.add(i)

    elif n == 1:
        for (k, j), value in A.items():
            if zero * value != zero:
                nan_idx.add(j)

    return nan_idx


def _dok_matmul_with_nan_handling(A, B, nan_idx_A, nan_idx_B, rows_A, cols_B):
    """
    Sparse matrix multiplication handling possible undefined results (e.g., 0 * infinity).
    It uses the indexes of elements in A and B that produce NaN values, stored in nan_idx_A and nan_idx_B.
    """

    C = Counter()
    C_int = Counter()

    # Use set operations to handle NaN propagation
    nan_rows = set(nan_idx_A)
    nan_cols = set(nan_idx_B)

    # Handle NaN propagation for A
    for i in nan_rows:
        for j in range(cols_B):
            C[(i, j)] = nan

    # Handle NaN propagation for B
    for j in nan_cols:
        for i in range(rows_A):
            C[(i, j)] = nan

    # Standard matrix multiplication to form C_int
    for (i, k), A_value in A.items():
        for (k2, j), B_value in B.items():
            if k == k2:
                C_int[(i, j)] += A_value * B_value

    # Combine C and C_int to form the final result C
    for key, value in C_int.items():
        C[key] = value

    return C


def _forward_jacobian_core(replacements, reduced_expr, wrt):
    """
    Core function for Jacobian matrix calculation through forward accumulation.
    Takes directly the output of a CSE operation, and an iterable of variables
    with respect to which to differentiate the reduced expression and returns the
    Jacobian matrix in DAG form.

    The function also returns a list of precomputed free symbols for each subexpression,
    which are useful in the substitution process.

    NOTE: The output Jacobian might not be always exactly in a DAG form, but close to it.

    Parameters
    ==========

    replacements : list
        A list of tuples containing the CSE replacement symbols and their corresponding
        expressions.

    reduced_expr : list
        The reduced expression after the CSE operation.

    wrt : Matrix
        The matrix of expressions with respect to which to differentiate the reduced
        expression.

    """

    if not isinstance(reduced_expr[0], MatrixBase):
        raise TypeError("``expr`` must be of matrix type")

    if not (reduced_expr[0].shape[0] == 1 or reduced_expr[0].shape[1] == 1):
        raise TypeError("``expr`` must be a row or a column matrix")

    if not isinstance(wrt, (MatrixBase, list, tuple)):
        raise TypeError("``wrt`` must be an iterable of variables")

    elif isinstance(wrt, (list, tuple)):
        wrt = Matrix(wrt)

    if not (wrt.shape[0] == 1 or wrt.shape[1] == 1):
        raise TypeError("``wrt`` must be a row or a column matrix")

    replacements, reduced_expr = _process_cse(replacements, reduced_expr)

    if replacements:
        rep_sym, sub_expr = map(Matrix, zip(*replacements))
    else:
        rep_sym, sub_expr = Matrix([]), Matrix([])

    l_sub, l_wrt, l_red = len(sub_expr), len(wrt), len(reduced_expr[0])

    f1 = {
        (i, j): diff_value
        for i, r in enumerate(reduced_expr[0])
        for j, w in enumerate(wrt)
        if (diff_value := r.diff(w)) != 0
    }

    if not replacements:
        return [], SparseMatrix(l_red, l_wrt, f1), []

    f2 = {
        (i, j): diff_value
        for i, (r, fs) in enumerate([(r, r.free_symbols) for r in reduced_expr[0]])
        for j, s in enumerate(rep_sym)
        if s in fs and (diff_value := r.diff(s)) != 0
    }

    precomputed_fs = [
        {symbol for symbol in s.free_symbols if re.compile(r'x\d+').fullmatch(symbol.name)}
        for s in sub_expr
    ]


    C = Counter({(0, j): diff_value for j, w in enumerate(wrt) if (diff_value := sub_expr[0].diff(w)) != 0})

    nan_idx_C = set()
    nan_idx_C = _check_nan(C, nan_idx_C, 1)

    for i in range(1, l_sub):
        Bi = {(i, j): diff_value for j in range(i + 1)
              if rep_sym[j] in precomputed_fs[i] and (diff_value := sub_expr[i].diff(rep_sym[j])) != 0}

        nan_idx_Bi = set()
        nan_idx_Bi = _check_nan(Bi, nan_idx_Bi, 0)

        Ai = Counter({(i, j): diff_value for j, w in enumerate(wrt)
                      if (diff_value := sub_expr[i].diff(w)) != 0})

        if Bi:
            Ci = _dok_matmul_with_nan_handling(Bi, C, nan_idx_Bi, nan_idx_C, 1, l_wrt)
            nan_idx_C = _check_nan(Ci, nan_idx_C, 1)

            Ci.update(Ai)
            C.update(Ci)
        else:
            C.update(Ai)

    nan_idx_f2 = set()
    nan_idx_f2 = _check_nan(f2, nan_idx_f2, 0)
    J = _dok_matmul_with_nan_handling(f2, C, nan_idx_f2, nan_idx_C, l_red, l_wrt)


    for (i, j), value in f1.items():
        J[(i, j)] += value

    J = SparseMatrix(l_red, l_wrt, J)
    J = reduced_expr[0].__class__(J)

    return replacements, J, precomputed_fs


def _forward_jacobian_norm_in_dag_out(expr, wrt):

    replacements, reduced_expr = cse(expr)
    replacements, J, _ = _forward_jacobian_core(replacements, reduced_expr, wrt)

    return replacements, J


def _forward_jacobian(expr, wrt):
    r"""
    Returns the Jacobian matrix produced using a forward accumulation
    algorithm.

    Explanation
    ===========

    Expressions often contain repeated subexpressions. Using a tree structure,
    these subexpressions are duplicated and differentiated multiple times,
    leading to inefficiency.

    Instead, if a data structure called a directed acyclic graph (DAG) is used
    then each of these repeated subexpressions will only exist a single time.
    This function uses a combination of representing the expression as a DAG and
    a forward accumulation algorithm (repeated application of the chain rule
    symbolically) to more efficiently calculate the Jacobian matrix of a target
    expression ``expr`` with respect to an expression or set of expressions
    ``wrt``.

    Note that this function is intended to improve performance when
    differentiating large expressions that contain many common subexpressions.
    For small and simple expressions it is likely less performant than using
    SymPy's standard differentiation functions and methods.

    NOTE: When Derivative terms are present in the expression, the CSE output
    is post-processed to remove any CSE replacement symbols from the arguments of those terms.
    Thus, in that case some derivatives might be repeated.

    Parameters
    ==========

    expr : Matrix
        The vector to be differentiated.

    wrt : Matrix, list, or tuple
        The vector with respect to which to do the differentiation. Can be a matrix or an iterable of variables.

    See Also
    ========

    Direct Acyclic Graph : https://en.wikipedia.org/wiki/Directed_acyclic_graph

    """

    replacements, reduced_expr = cse(expr)

    if replacements:
        rep_sym, _ = map(Matrix, zip(*replacements))
    else:
        rep_sym = Matrix([])

    l_wrt = len(wrt)
    l_red = len(reduced_expr[0])

    replacements, J, precomputed_fs = _forward_jacobian_core(replacements, reduced_expr, wrt)

    if not replacements: return J

    sub_rep = dict(replacements)
    for i, ik in enumerate(precomputed_fs):
        sub_dict = {j: sub_rep[j] for j in ik}
        sub_rep[rep_sym[i]] = sub_rep[rep_sym[i]].xreplace(sub_dict)

    J = {key: expr.xreplace(sub_rep) for key, expr in J.todok().items()}
    J = SparseMatrix(l_red, l_wrt, J)
    J = expr.__class__(J)

    return J
