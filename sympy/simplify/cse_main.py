""" Tools for doing common subexpression elimination.
"""

from sympy import Basic
from sympy.utilities.iterables import postorder_traversal, numbered_symbols

import cse_opts

# (preprocessor, postprocessor) pairs which are commonly useful. They should
# each take a sympy expression and return a possibly transformed expression.
# When used in the function `cse()`, the target expressions will be transformed
# by each of the preprocessor functions in order. After the common
# subexpressions are eliminated, each resulting expression will have the
# postprocessor functions transform them in *reverse* order in order to undo the
# transformation if necessary. This allows the algorithm to operate on
# a representation of the expressions that allows for more optimization
# opportunities.
# `None` can be used to specify no transformation for either the preprocessor or
# postprocessor.
cse_optimizations = list(cse_opts.default_optimizations)

def preprocess_for_cse(expr, optimizations):
    """ Preprocess an expression to optimize for common subexpression
    elimination.

    Parameters
    ----------
    expr : sympy expression
        The target expression to optimize.
    optimizations : list of (callable, callable) pairs
        The (preprocessor, postprocessor) pairs.

    Returns
    -------
    expr : sympy expression
        The transformed expression.
    """
    for pre, post in optimizations:
        if pre is not None:
            expr = pre(expr)
    return expr

def postprocess_for_cse(expr, optimizations):
    """ Postprocess an expression after common subexpression elimination to
    return the expression to canonical sympy form.

    Parameters
    ----------
    expr : sympy expression
        The target expression to transform.
    optimizations : list of (callable, callable) pairs, optional
        The (preprocessor, postprocessor) pairs.  The postprocessors will be
        applied in reversed order to undo the effects of the preprocessors
        correctly.

    Returns
    -------
    expr : sympy expression
        The transformed expression.
    """
    if optimizations is None:
        optimizations = cse_optimizations
    for pre, post in reversed(optimizations):
        if post is not None:
            expr = post(expr)
    return expr

def cse(exprs, symbols=None, optimizations=None):
    """ Perform common subexpression elimination on an expression.

    Parameters:

    exprs : list of sympy expressions, or a single sympy expression
        The expressions to reduce.
    symbols : infinite iterator yielding unique Symbols
        The symbols used to label the common subexpressions which are pulled
        out. The `numbered_symbols` generator is useful. The default is a stream
        of symbols of the form "x0", "x1", etc. This must be an infinite
        iterator.
    optimizations : list of (callable, callable) pairs, optional
        The (preprocessor, postprocessor) pairs. If not provided,
        `sympy.simplify.cse.cse_optimizations` is used.

    Returns:

    replacements : list of (Symbol, expression) pairs
        All of the common subexpressions that were replaced. Subexpressions
        earlier in this list might show up in subexpressions later in this list.
    reduced_exprs : list of sympy expressions
        The reduced expressions with all of the replacements above.
    """
    if symbols is None:
        symbols = numbered_symbols()
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)
    seen_subexp = set()
    to_eliminate = []

    if optimizations is None:
        # Pull out the default here just in case there are some weird
        # manipulations of the module-level list in some other thread.
        optimizations = list(cse_optimizations)

    # Handle the case if just one expression was passed.
    if isinstance(exprs, Basic):
        exprs = [exprs]
    # Preprocess the expressions to give us better optimization opportunities.
    exprs = [preprocess_for_cse(e, optimizations) for e in exprs]

    # Find all of the repeated subexpressions.
    for expr in exprs:
        for subtree in postorder_traversal(expr):
            if subtree.args == ():
                # Exclude atoms, since there is no point in renaming them.
                continue
            if (subtree.args != () and
                subtree in seen_subexp and
                subtree not in to_eliminate):
                to_eliminate.append(subtree)
            seen_subexp.add(subtree)

    # Substitute symbols for all of the repeated subexpressions.
    replacements = []
    reduced_exprs = list(exprs)
    for i, subtree in enumerate(to_eliminate):
        sym = symbols.next()
        replacements.append((sym, subtree))
        # Make the substitution in all of the target expressions.
        for j, expr in enumerate(reduced_exprs):
            reduced_exprs[j] = expr.subs(subtree, sym)
        # Make the substitution in all of the subsequent substitutions.
        # WARNING: modifying iterated list in-place! I think it's fine,
        # but there might be clearer alternatives.
        for j in range(i+1, len(to_eliminate)):
            to_eliminate[j] = to_eliminate[j].subs(subtree, sym)

    # Postprocess the expressions to return the expressions to canonical form.
    for i, (sym, subtree) in enumerate(replacements):
        subtree = postprocess_for_cse(subtree, optimizations)
        replacements[i] = (sym, subtree)
    reduced_exprs = [postprocess_for_cse(e, optimizations) for e in reduced_exprs]

    return replacements, reduced_exprs
