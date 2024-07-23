"""Module for differentiation using CSE."""

from sympy import cse, Matrix, SparseMatrix, Derivative, MatrixBase
from collections import Counter
import re

def _postprocess(repl, reduced):
    """
    Postprocess the CSE output to remove any CSE replacement symbols from the arguments of Derivative terms.
    """

    repl_dict = dict(repl)

    p_repl = [(rep_sym, _traverse(sub_exp, repl_dict)) for rep_sym, sub_exp in repl]
    p_reduced = [Matrix([_traverse(exp, repl_dict) for exp in red_exp]) for red_exp in reduced]

    return p_repl, p_reduced


def _traverse(node, repl_dict):
    """
    Traverse the node in preorder fashion, and apply replace_all() if the node
    is the argument of a Derivative.
    """

    if isinstance(node, Derivative):
        return _replace_all(node, repl_dict)

    if not node.args:
        return node

    new_args = [_traverse(arg, repl_dict) for arg in node.args]
    return node.func(*new_args)


def _replace_all(node, repl_dict):
    """
    Bring the node to its form before the CSE operation, by iteratively substituting
    the CSE replacement symbols in the node.
    """

    result = node
    while True:
        fs = result.free_symbols
        sl_dict = {k: repl_dict[k] for k in fs if k in repl_dict}
        if not sl_dict:
            break
        result = result.xreplace(sl_dict)
    return result


def _dok_matrix_multiply(A, B):
    """
    Multiply two sparse matrices in dok format (i, k) and (k, j) to get a dictionary of keys (i, j).
    """
    result = Counter()
    for (i, k), A_value in A.items():
        for (k2, j), B_value in B.items():
            if k == k2:
                result[(i, j)] += A_value * B_value
    return result


def _forward_jacobian(expr, wrt):
    r"""
    Returns the Jacobian matrix produced using a forward accumulation
    algorithm.

    Explanation
    ===========

    Expressions often contain repeated subexpressions. If and expression is
    represented as a tree structure then multiple copies of these subexpressions
    will be present in the expanded form of the expression. During
    differentiation these repeated subexpressions will be repeatedly and
    differentiated multiple times, resulting in repeated and wasted work.

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

    # Check Inputs
    if not isinstance(wrt, (MatrixBase, list, tuple)):
        raise TypeError("``wrt`` must be an iterable of variables")

    elif isinstance(wrt, (list, tuple)):
        wrt = Matrix(wrt)

    if not isinstance(expr, MatrixBase):
        raise TypeError("``expr`` must be of matrix type")

    # Both ``wrt`` and ``expr`` can be a row or a column matrix, so we need to make
    # sure all valid combinations work, but everything else fails:
    if not (expr.shape[0] == 1 or expr.shape[1] == 1):
        raise TypeError("``expr`` must be a row or a column matrix")

    if not (wrt.shape[0] == 1 or wrt.shape[1] == 1):
        raise TypeError("``wrt`` must be a row or a column matrix")

    replacements, reduced_expr = cse(expr)
    replacements, reduced_expr = _postprocess(replacements, reduced_expr)

    if replacements:
        rep_sym, sub_expr = map(Matrix, zip(*replacements))
    else:
        rep_sym, sub_expr, reduced_expr = Matrix([]), Matrix([]), [expr]

    l_sub, l_wrt, l_red = len(sub_expr), len(wrt), len(reduced_expr[0])

    f1 = {
        (i, j): diff_value
        for i, r in enumerate(reduced_expr[0])
        for j, w in enumerate(wrt)
        if (diff_value := r.diff(w)) != 0
    }

    if not replacements:
        return SparseMatrix(l_red, l_wrt, f1)

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

    for i in range(1, l_sub):
        Bi = {(i, j): diff_value for j in range(i + 1)
              if rep_sym[j] in precomputed_fs[i] and (diff_value := sub_expr[i].diff(rep_sym[j])) != 0}

        Ai = Counter({(i, j): diff_value for j, w in enumerate(wrt)
                      if (diff_value := sub_expr[i].diff(w)) != 0})

        if Bi:
            Ci = _dok_matrix_multiply(Bi, C)
            Ci.update(Ai)  # Use Counter's update method to add Ai to Ci
            C.update(Ci)  # Update C with the result
        else:
            C.update(Ai)

    J = _dok_matrix_multiply(f2, C)

    for (i, j), value in f1.items():
        J[(i, j)] += value

    sub_rep = dict(replacements)
    for i, ik in enumerate(precomputed_fs):
        sub_dict = {j: sub_rep[j] for j in ik}
        sub_rep[rep_sym[i]] = sub_rep[rep_sym[i]].xreplace(sub_dict)

    J = {key: expr.xreplace(sub_rep) for key, expr in J.items()}
    J = SparseMatrix(l_red, l_wrt, J)
    J = expr.__class__(J)

    return J
