""" Tools for doing common subexpression elimination.
"""
import bisect
import difflib

from sympy import Basic, Mul, Add
from sympy.utilities.iterables import preorder_traversal, numbered_symbols

import cse_opts

# (preprocessor, postprocessor) pairs which are commonly useful. They should
# each take a sympy expression and return a possibly transformed expression.
# When used in the function ``cse()``, the target expressions will be transformed
# by each of the preprocessor functions in order. After the common
# subexpressions are eliminated, each resulting expression will have the
# postprocessor functions transform them in *reverse* order in order to undo the
# transformation if necessary. This allows the algorithm to operate on
# a representation of the expressions that allows for more optimization
# opportunities.
# ``None`` can be used to specify no transformation for either the preprocessor or
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
        out. The ``numbered_symbols`` generator is useful. The default is a stream
        of symbols of the form "x0", "x1", etc. This must be an infinite
        iterator.
    optimizations : list of (callable, callable) pairs, optional
        The (preprocessor, postprocessor) pairs. If not provided,
        ``sympy.simplify.cse.cse_optimizations`` is used.

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
    muls = set()
    adds = set()
    to_eliminate = []
    to_eliminate_ops_count = []

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
    def insert(subtree):
        '''This helper will insert the subtree into to_eliminate while
        maintaining the ordering by op count and will skip the insertion
        if subtree is already present.'''
        ops_count = subtree.count_ops()
        index_to_insert = bisect.bisect(to_eliminate_ops_count, ops_count)
        # all i up to this index have op count <= the current op count
        # so check that subtree is not yet present from this index down
        # (if necessary) to zero.
        for i in xrange(index_to_insert - 1, -1, -1):
            if to_eliminate_ops_count[i] == ops_count and \
               subtree == to_eliminate[i]:
                return # already have it
        to_eliminate_ops_count.insert(index_to_insert, ops_count)
        to_eliminate.insert(index_to_insert, subtree)

    for expr in exprs:
        pt = preorder_traversal(expr)
        for subtree in pt:
            if subtree.is_Atom:
                # Exclude atoms, since there is no point in renaming them.
                continue

            if subtree in seen_subexp:
                insert(subtree)
                pt.skip()
                continue

            if subtree.is_Mul:
                muls.add(subtree)
            elif subtree.is_Add:
                adds.add(subtree)

            seen_subexp.add(subtree)

    # process adds - any adds that weren't repeated might contain
    # subpatterns that are repeated, e.g. x+y+z and x+y have x+y in common
    adds = [set(a.args) for a in adds]
    for i in xrange(len(adds)):
        for j in xrange(i + 1, len(adds)):
            com = adds[i].intersection(adds[j])
            if len(com) > 1:
                insert(Add(*com))

                # remove this set of symbols so it doesn't appear again
                adds[i] = adds[i].difference(com)
                adds[j] = adds[j].difference(com)
                for k in xrange(j + 1, len(adds)):
                    if not com.difference(adds[k]):
                        adds[k] = adds[k].difference(com)

    # process muls - any muls that weren't repeated might contain
    # subpatterns that are repeated, e.g. x*y*z and x*y have x*y in common

    # use SequenceMatcher on the nc part to find the longest common expression
    # in common between the two nc parts
    sm = difflib.SequenceMatcher()

    muls = [a.args_cnc() for a in muls]
    for i in xrange(len(muls)):
        if muls[i][1]:
            sm.set_seq1(muls[i][1])
        for j in xrange(i + 1, len(muls)):
            # the commutative part in common
            ccom = muls[i][0].intersection(muls[j][0])

            # the non-commutative part in common
            if muls[i][1] and muls[j][1]:
                # see if there is any chance of an nc match
                ncom = set(muls[i][1]).intersection(set(muls[j][1]))
                if len(ccom) + len(ncom) < 2:
                    continue

                # now work harder to find the match
                sm.set_seq2(muls[j][1])
                i1, _, n = sm.find_longest_match(0, len(muls[i][1]),
                                                 0, len(muls[j][1]))
                ncom = muls[i][1][i1:i1 + n]
            else:
                ncom = []

            com = list(ccom) + ncom
            if len(com) < 2:
                continue

            insert(Mul(*com))

            # remove ccom from all if there was no ncom; to update the nc part
            # would require finding the subexpr and then replacing it with a
            # dummy to keep bounding nc symbols from being identified as a
            # subexpr, e.g. removing B*C from A*B*C*D might allow A*D to be
            # identified as a subexpr which would not be right.
            if not ncom:
                muls[i][0] = muls[i][0].difference(ccom)
                for k in xrange(j, len(muls)):
                    if not ccom.difference(muls[k][0]):
                        muls[k][0] = muls[k][0].difference(ccom)

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
        for j in range(i+1, len(to_eliminate)):
            to_eliminate[j] = to_eliminate[j].subs(subtree, sym)

    # Postprocess the expressions to return the expressions to canonical form.
    for i, (sym, subtree) in enumerate(replacements):
        subtree = postprocess_for_cse(subtree, optimizations)
        replacements[i] = (sym, subtree)
    reduced_exprs = [postprocess_for_cse(e, optimizations) for e in reduced_exprs]

    return replacements, reduced_exprs
