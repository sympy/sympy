"""Module for differentiation using CSE."""

from collections import Counter

from sympy.core.containers import Tuple
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.traversal import postorder_traversal
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.simplify.cse_main import CseExpr, cse
from sympy.utilities.iterables import numbered_symbols


def forward_jacobian(
    expr: ImmutableDenseMatrix,
    wrt: ImmutableDenseMatrix,
    as_cse_expr: bool = True,
):
    """Returns the Jacobian matrix produced using a forward accumulation
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

    Parameters
    ==========

    expr : ``ImmutableDenseMatrix``
        The vector to be differentiated.
    wrt : ``ImmutableDenseMatrix``
        The vector with respect to which to do the differentiation.
    as_cse_expr : ``bool``
        Influences the return type. If ``False``, then a matrix with fully-
        replaced SymPy expressions for entries will be returned. If ``True``,
        then the return type will be left as a ``CseExpr`` where the matrix's
        entries are reduced expressions containing replacements. The default is
        ``True``.

    See Also
    ========

    Direct Acyclic Graph : https://en.wikipedia.org/wiki/Directed_acyclic_graph

    """

    def add_to_cache(node):
        if node in expr_to_replacement_cache:
            replacement_symbol = expr_to_replacement_cache[node]
            return replacement_symbol, replacement_to_reduced_expr_cache[replacement_symbol]
        elif node in replacement_to_reduced_expr_cache:
            return node, replacement_to_reduced_expr_cache[node]
        elif isinstance(node, Tuple):
            return None, None
        elif not node.free_symbols:
            return node, node

        replacement_symbol = replacement_symbols.__next__()
        replaced_subexpr = node.xreplace(expr_to_replacement_cache)
        replacement_to_reduced_expr_cache[replacement_symbol] = replaced_subexpr
        expr_to_replacement_cache[node] = replacement_symbol
        return replacement_symbol, replaced_subexpr

    if not isinstance(expr, ImmutableDenseMatrix):
        msg = (
            'The forward Jacobian differentiation algorithm can only be used '
            'to differentiate a single matrix expression at a time.'
        )
        raise NotImplementedError(msg)
    elif expr.shape[1] != 1:
        msg = 'Can only compute the Jacobian for column matrices.'
        raise NotImplementedError(msg)
    elif not isinstance(wrt, ImmutableDenseMatrix) or wrt.shape[1] != 1:
        msg = (
            'The forward Jacobian differentiation algorithm can only compute '
            'Jacobians with respect to column matrices.'
        )
        raise NotImplementedError(msg)

    replacement_symbols = numbered_symbols(
        prefix='_z',
        cls=Symbol,
        exclude=expr.free_symbols,
    )

    expr_to_replacement_cache = {}
    replacement_to_reduced_expr_cache = {}

    replacements, reduced_exprs = cse(expr.args[2], replacement_symbols)
    for replacement_symbol, reduced_subexpr in replacements:
        replaced_subexpr = reduced_subexpr.xreplace(expr_to_replacement_cache)
        replacement_to_reduced_expr_cache[replacement_symbol] = replaced_subexpr
        expr_to_replacement_cache[reduced_subexpr] = replacement_symbol
        for node in postorder_traversal(reduced_subexpr):
            _ = add_to_cache(node)
    for reduced_expr in reduced_exprs:
        for node in reduced_expr:
            _ = add_to_cache(node)

    reduced_matrix = ImmutableDenseMatrix(reduced_exprs).xreplace(expr_to_replacement_cache)
    replacements = list(replacement_to_reduced_expr_cache.items())

    absolute_derivative_mapping = {}
    for i, wrt_symbol in enumerate(wrt.args[2]):
        absolute_derivative = [S.Zero] * len(wrt)
        absolute_derivative[i] = S.One
        absolute_derivative_mapping[wrt_symbol] = ImmutableDenseMatrix([absolute_derivative])

    zeros = ImmutableDenseMatrix.zeros(1, len(wrt))
    for symbol, subexpr in replacements:
        free_symbols = subexpr.free_symbols
        absolute_derivative = zeros
        for free_symbol in free_symbols:
            replacement_symbol, partial_derivative = add_to_cache(subexpr.diff(free_symbol))
            absolute_derivative += partial_derivative * absolute_derivative_mapping.get(free_symbol, zeros)
        absolute_derivative_mapping[symbol] = ImmutableDenseMatrix([[add_to_cache(a)[0] for a in absolute_derivative]])

    replaced_jacobian = ImmutableDenseMatrix.vstack(*[absolute_derivative_mapping.get(e, ImmutableDenseMatrix.zeros(*wrt.shape).T) for e in reduced_matrix])

    required_replacement_symbols = set()
    stack = [entry for entry in replaced_jacobian if entry.free_symbols]
    while stack:
        entry = stack.pop()
        if entry in required_replacement_symbols or entry in wrt:
            continue
        children = list(replacement_to_reduced_expr_cache.get(entry, entry).free_symbols)
        for child in children:
            if child not in required_replacement_symbols and child not in wrt:
                stack.append(child)
        required_replacement_symbols.add(entry)

    required_replacements_dense = {
        replacement_symbol: replaced_subexpr
        for replacement_symbol, replaced_subexpr in replacement_to_reduced_expr_cache.items()
        if replacement_symbol in required_replacement_symbols
    }

    counter = Counter(replaced_jacobian.free_symbols)
    for replaced_subexpr in required_replacements_dense.values():
        counter.update(replaced_subexpr.free_symbols)

    required_replacements = {}
    unrequired_replacements = {}
    for replacement_symbol, replaced_subexpr in required_replacements_dense.items():
        if isinstance(replaced_subexpr, Symbol) or counter[replacement_symbol] == 1:
            unrequired_replacements[replacement_symbol] = replaced_subexpr.xreplace(unrequired_replacements)
        else:
            required_replacements[replacement_symbol] = replaced_subexpr.xreplace(unrequired_replacements)

    cse_expr = CseExpr((list(required_replacements.items()), replaced_jacobian.xreplace(unrequired_replacements)))
    if as_cse_expr:
        return cse_expr
    return cse_expr.reduced_exprs.subs(reversed(cse_expr.replacements))
