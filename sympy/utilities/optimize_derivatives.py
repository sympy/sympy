"""
Utilities for joint CSE optimization of a scalar function together with
its gradient and Hessian.

This module provides a helper to compute a function, its first derivatives,
and second derivatives, and apply common subexpression elimination (CSE)
jointly across all of them.
"""

def generate_optimized_derivatives(expr, variables):
    """
    Compute and CSE-optimize a scalar expression together with its
    gradient and Hessian.

    Parameters
    ----------
    expr : Expr
        A scalar SymPy expression.
    variables : sequence of Symbol
        Variables with respect to which derivatives are taken.

    Returns
    -------
    replacements : list of (Symbol, Expr)
        Common subexpressions extracted by CSE.
    reduced : list of Expr
        Reduced expressions in the following layout:

            reduced[0]          -> f
            reduced[1 : n+1]    -> gradient
            reduced[n+1 :]      -> Hessian (row-major, flattened)

    Notes
    -----
    Applying CSE jointly across the function, gradient, and Hessian can
    reduce redundant computations compared to applying CSE separately.
    """
    from sympy.simplify.cse_main import cse  

    variables = list(variables)
    n = len(variables)

    # gradient
    grad = [expr.diff(v) for v in variables]

    # hessian
    hess_flat = [
        expr.diff(vi).diff(vj)
        for vi in variables
        for vj in variables
    ]

    # joint CSE
    all_exprs = [expr] + grad + hess_flat
    replacements, reduced = cse(all_exprs)

    return replacements, reduced


def count_operations(replacements, reduced):
    """
    Count total operations in CSE-optimized expressions.
    """
    from sympy import count_ops

    total = sum(count_ops(r) for _, r in replacements)
    total += sum(count_ops(r) for r in reduced)
    return total