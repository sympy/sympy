""" Tools for doing common subexpression elimination.
"""
from __future__ import print_function, division

from sympy.core import Basic, Mul, Add, Pow, sympify, Symbol, Tuple
from sympy.core.singleton import S
from sympy.core.function import _coeff_isneg
from sympy.core.exprtools import factor_terms
from sympy.core.compatibility import iterable, range
from sympy.utilities.iterables import filter_symbols, \
    numbered_symbols, sift, topological_sort, ordered

from . import cse_opts

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


basic_optimizations = [(cse_opts.sub_pre, cse_opts.sub_post),
                       (factor_terms, None)]

# sometimes we want the output in a different format; non-trivial
# transformations can be put here for users
# ===============================================================


def reps_toposort(r):
    """Sort replacements `r` so (k1, v1) appears before (k2, v2)
    if k2 is in v1's free symbols. This orders items in the
    way that cse returns its results (hence, in order to use the
    replacements in a substitution option it would make sense
    to reverse the order).

    Examples
    ========

    >>> from sympy.simplify.cse_main import reps_toposort
    >>> from sympy.abc import x, y
    >>> from sympy import Eq
    >>> for l, r in reps_toposort([(x, y + 1), (y, 2)]):
    ...     print(Eq(l, r))
    ...
    Eq(y, 2)
    Eq(x, y + 1)

    """
    r = sympify(r)
    E = []
    for c1, (k1, v1) in enumerate(r):
        for c2, (k2, v2) in enumerate(r):
            if k1 in v2.free_symbols:
                E.append((c1, c2))
    return [r[i] for i in topological_sort((range(len(r)), E))]


def cse_separate(r, e):
    """Move expressions that are in the form (symbol, expr) out of the
    expressions and sort them into the replacements using the reps_toposort.

    Examples
    ========

    >>> from sympy.simplify.cse_main import cse_separate
    >>> from sympy.abc import x, y, z
    >>> from sympy import cos, exp, cse, Eq, symbols
    >>> x0, x1 = symbols('x:2')
    >>> eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    >>> cse([eq, Eq(x, z + 1), z - 2], postprocess=cse_separate) in [
    ... [[(x0, y + 1), (x, z + 1), (x1, x + 1)],
    ...  [x1 + exp(x1/x0) + cos(x0), z - 2]],
    ... [[(x1, y + 1), (x, z + 1), (x0, x + 1)],
    ...  [x0 + exp(x0/x1) + cos(x1), z - 2]]]
    ...
    True
    """
    d = sift(e, lambda w: w.is_Equality and w.lhs.is_Symbol)
    r = r + [w.args for w in d[True]]
    e = d[False]
    return [reps_toposort(r), e]

# ====end of cse postprocess idioms===========================


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
    for pre, post in reversed(optimizations):
        if post is not None:
            expr = post(expr)
    return expr


def opt_cse(exprs, order='canonical'):
    """Find optimization opportunities in Adds, Muls, Pows and negative
    coefficient Muls

    Parameters
    ----------
    exprs : list of sympy expressions
        The expressions to optimize.
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. For large
        expressions where speed is a concern, use the setting order='none'.

    Returns
    -------
    opt_subs : dictionary of expression substitutions
        The expression substitutions which can be useful to optimize CSE.

    Examples
    ========

    >>> from sympy.simplify.cse_main import opt_cse
    >>> from sympy.abc import x
    >>> opt_subs = opt_cse([x**-2])
    >>> print(opt_subs)
    {x**(-2): 1/(x**2)}
    """
    from sympy.matrices.expressions import MatAdd, MatMul, MatPow
    opt_subs = dict()

    adds = set()
    muls = set()

    seen_subexp = set()

    def _find_opts(expr):

        if not isinstance(expr, Basic):
            return

        if expr.is_Atom or expr.is_Order:
            return

        if iterable(expr):
            list(map(_find_opts, expr))
            return

        if expr in seen_subexp:
            return expr
        seen_subexp.add(expr)

        list(map(_find_opts, expr.args))

        if _coeff_isneg(expr):
            neg_expr = -expr
            if not neg_expr.is_Atom:
                opt_subs[expr] = Mul(S.NegativeOne, neg_expr, evaluate=False)
                seen_subexp.add(neg_expr)
                expr = neg_expr

        if isinstance(expr, (Mul, MatMul)):
            muls.add(expr)

        elif isinstance(expr, (Add, MatAdd)):
            adds.add(expr)

        elif isinstance(expr, (Pow, MatPow)):
            if _coeff_isneg(expr.exp):
                opt_subs[expr] = Pow(Pow(expr.base, -expr.exp), S.NegativeOne,
                                     evaluate=False)

    for e in exprs:
        if isinstance(e, Basic):
            _find_opts(e)

    ## Process Adds and commutative Muls

    def _match_common_args(Func, funcs):
        if order != 'none':
            funcs = list(ordered(funcs))
        else:
            funcs = sorted(funcs, key=lambda x: len(x.args))

        func_args = [set(e.args) for e in funcs]
        for i in range(len(func_args)):
            for j in range(i + 1, len(func_args)):
                com_args = func_args[i].intersection(func_args[j])
                if len(com_args) > 1:
                    com_func = Func(*com_args)

                    # for all sets, replace the common symbols by the function
                    # over them, to allow recursive matches

                    diff_i = func_args[i].difference(com_args)
                    func_args[i] = diff_i | set([com_func])
                    if diff_i:
                        opt_subs[funcs[i]] = Func(Func(*diff_i), com_func,
                                                  evaluate=False)

                    diff_j = func_args[j].difference(com_args)
                    func_args[j] = diff_j | set([com_func])
                    opt_subs[funcs[j]] = Func(Func(*diff_j), com_func,
                                              evaluate=False)

                    for k in range(j + 1, len(func_args)):
                        if not com_args.difference(func_args[k]):
                            diff_k = func_args[k].difference(com_args)
                            func_args[k] = diff_k | set([com_func])
                            opt_subs[funcs[k]] = Func(Func(*diff_k), com_func,
                                                      evaluate=False)

    # split muls into commutative
    comutative_muls = set()
    for m in muls:
        c, nc = m.args_cnc(cset=True)
        if c:
            c_mul = m.func(*c)
            if nc:
                if c_mul == 1:
                    new_obj = m.func(*nc)
                else:
                    new_obj = m.func(c_mul, m.func(*nc), evaluate=False)
                opt_subs[m] = new_obj
            if len(c) > 1:
                comutative_muls.add(c_mul)

    _match_common_args(Add, adds)
    _match_common_args(Mul, comutative_muls)

    return opt_subs


def tree_cse(exprs, symbols, opt_subs=None, order='canonical'):
    """Perform raw CSE on expression tree, taking opt_subs into account.

    Parameters
    ==========

    exprs : list of sympy expressions
        The expressions to reduce.
    symbols : infinite iterator yielding unique Symbols
        The symbols used to label the common subexpressions which are pulled
        out.
    opt_subs : dictionary of expression substitutions
        The expressions to be substituted before any CSE action is performed.
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. For large
        expressions where speed is a concern, use the setting order='none'.
    """
    from sympy.matrices.expressions import MatrixExpr, MatrixSymbol, MatMul, MatAdd

    if opt_subs is None:
        opt_subs = dict()

    ## Find repeated sub-expressions

    to_eliminate = set()

    seen_subexp = set()

    def _find_repeated(expr):
        if not isinstance(expr, Basic):
            return

        if expr.is_Atom or expr.is_Order:
            return

        if iterable(expr):
            args = expr

        else:
            if expr in seen_subexp:
                to_eliminate.add(expr)
                return

            seen_subexp.add(expr)

            if expr in opt_subs:
                expr = opt_subs[expr]

            args = expr.args

        list(map(_find_repeated, args))

    for e in exprs:
        if isinstance(e, Basic):
            _find_repeated(e)

    ## Rebuild tree

    replacements = []

    subs = dict()

    def _rebuild(expr):
        if not isinstance(expr, Basic):
            return expr

        if not expr.args:
            return expr

        if iterable(expr):
            new_args = [_rebuild(arg) for arg in expr]
            return expr.func(*new_args)

        if expr in subs:
            return subs[expr]

        orig_expr = expr
        if expr in opt_subs:
            expr = opt_subs[expr]

        # If enabled, parse Muls and Adds arguments by order to ensure
        # replacement order independent from hashes
        if order != 'none':
            if isinstance(expr, (Mul, MatMul)):
                c, nc = expr.args_cnc()
                if c == [1]:
                    args = nc
                else:
                    args = list(ordered(c)) + nc
            elif isinstance(expr, (Add, MatAdd)):
                args = list(ordered(expr.args))
            else:
                args = expr.args
        else:
            args = expr.args

        new_args = list(map(_rebuild, args))
        if new_args != args:
            new_expr = expr.func(*new_args)
        else:
            new_expr = expr

        if orig_expr in to_eliminate:
            try:
                sym = next(symbols)
            except StopIteration:
                raise ValueError("Symbols iterator ran out of symbols.")

            if isinstance(orig_expr, MatrixExpr):
                sym = MatrixSymbol(sym.name, orig_expr.rows,
                    orig_expr.cols)

            subs[orig_expr] = sym
            replacements.append((sym, new_expr))
            return sym

        else:
            return new_expr

    reduced_exprs = []
    for e in exprs:
        if isinstance(e, Basic):
            reduced_e = _rebuild(e)
        else:
            reduced_e = e
        reduced_exprs.append(reduced_e)

    return replacements, reduced_exprs


def cse(exprs, symbols=None, optimizations=None, postprocess=None,
        order='canonical'):
    """ Perform common subexpression elimination on an expression.

    Parameters
    ==========

    exprs : list of sympy expressions, or a single sympy expression
        The expressions to reduce.
    symbols : infinite iterator yielding unique Symbols
        The symbols used to label the common subexpressions which are pulled
        out. The ``numbered_symbols`` generator is useful. The default is a
        stream of symbols of the form "x0", "x1", etc. This must be an
        infinite iterator.
    optimizations : list of (callable, callable) pairs
        The (preprocessor, postprocessor) pairs of external optimization
        functions. Optionally 'basic' can be passed for a set of predefined
        basic optimizations. Such 'basic' optimizations were used by default
        in old implementation, however they can be really slow on larger
        expressions. Now, no pre or post optimizations are made by default.
    postprocess : a function which accepts the two return values of cse and
        returns the desired form of output from cse, e.g. if you want the
        replacements reversed the function might be the following lambda:
        lambda r, e: return reversed(r), e
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. If set to
        'canonical', arguments will be canonically ordered. If set to 'none',
        ordering will be faster but dependent on expressions hashes, thus
        machine dependent and variable. For large expressions where speed is a
        concern, use the setting order='none'.

    Returns
    =======

    replacements : list of (Symbol, expression) pairs
        All of the common subexpressions that were replaced. Subexpressions
        earlier in this list might show up in subexpressions later in this
        list.
    reduced_exprs : list of sympy expressions
        The reduced expressions with all of the replacements above.

    Examples
    ========

    >>> from sympy import cse, SparseMatrix
    >>> from sympy.abc import x, y, z, w
    >>> cse(((w + x + y + z)*(w + y + z))/(w + x)**3)
    ([(x0, y + z), (x1, w + x)], [(w + x0)*(x0 + x1)/x1**3])

    Note that currently, y + z will not get substituted if -y - z is used.

     >>> cse(((w + x + y + z)*(w - y - z))/(w + x)**3)
     ([(x0, w + x)], [(w - y - z)*(x0 + y + z)/x0**3])

    List of expressions with recursive substitutions:

    >>> m = SparseMatrix([x + y, x + y + z])
    >>> cse([(x+y)**2, x + y + z, y + z, x + z + y, m])
    ([(x0, x + y), (x1, x0 + z)], [x0**2, x1, y + z, x1, Matrix([
    [x0],
    [x1]])])

    Note: the type and mutability of input matrices is retained.

    >>> isinstance(_[1][-1], SparseMatrix)
    True
    """
    from sympy.matrices import (MatrixBase, Matrix, ImmutableMatrix,
                                SparseMatrix, ImmutableSparseMatrix)

    # Handle the case if just one expression was passed.
    if isinstance(exprs, (Basic, MatrixBase)):
        exprs = [exprs]

    copy = exprs
    temp = []
    for e in exprs:
        if isinstance(e, (Matrix, ImmutableMatrix)):
            temp.append(Tuple(*e._mat))
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            temp.append(Tuple(*e._smat.items()))
        else:
            temp.append(e)
    exprs = temp
    del temp

    if optimizations is None:
        optimizations = list()
    elif optimizations == 'basic':
        optimizations = basic_optimizations

    # Preprocess the expressions to give us better optimization opportunities.
    reduced_exprs = [preprocess_for_cse(e, optimizations) for e in exprs]

    excluded_symbols = set().union(*[expr.atoms(Symbol)
                                   for expr in reduced_exprs])

    if symbols is None:
        symbols = numbered_symbols()
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)

    symbols = filter_symbols(symbols, excluded_symbols)

    # Find other optimization opportunities.
    opt_subs = opt_cse(reduced_exprs, order)

    # Main CSE algorithm.
    replacements, reduced_exprs = tree_cse(reduced_exprs, symbols, opt_subs,
                                           order)

    # Postprocess the expressions to return the expressions to canonical form.
    exprs = copy
    for i, (sym, subtree) in enumerate(replacements):
        subtree = postprocess_for_cse(subtree, optimizations)
        replacements[i] = (sym, subtree)
    reduced_exprs = [postprocess_for_cse(e, optimizations)
                     for e in reduced_exprs]

    # Get the matrices back
    for i, e in enumerate(exprs):
        if isinstance(e, (Matrix, ImmutableMatrix)):
            reduced_exprs[i] = Matrix(e.rows, e.cols, reduced_exprs[i])
            if isinstance(e, ImmutableMatrix):
                reduced_exprs[i] = reduced_exprs[i].as_immutable()
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            m = SparseMatrix(e.rows, e.cols, {})
            for k, v in reduced_exprs[i]:
                m[k] = v
            if isinstance(e, ImmutableSparseMatrix):
                m = m.as_immutable()
            reduced_exprs[i] = m

    if postprocess is None:
        return replacements, reduced_exprs

    return postprocess(replacements, reduced_exprs)


class Marker:
    pass

# XXX: Should this go in sympy.utilities.iterables?
def shortest_repeated_subsequence(S, ignore=(Marker,)):
    """
    Given a sequence S, find the shortest repeated subsequence of length >= 2

    A subsequence is only considered a repeated subsequence if it can't be
    extended in every instance. For example, [a, b] is not a subsequence of
    [a, b, c, a, b, c] because every occurrence of [a, b] can be extended to
    [a, b, c].

    If no repeated subsequence is found, returns None.

    The input sequence should support iteration, indexing, and slicing.

    To use concurrently with multiple sequences, concatenate them separated by
    markers (custom markers can be pass in through the 'ignore' parameter').

    Examples
    ========

    In the examples here, we use strings. A more common example for SymPy
    would be the args of a noncommutative Mul.

    >>> from sympy.simplify.cse_main import shortest_repeated_subsequence
    >>> s = 'abcabc'
    >>> shortest_repeated_subsequence(s)
    'abc'
    >>> print(shortest_repeated_subsequence('abc'))
    None

    >>> from itertools import chain
    >>> from sympy.utilities.iterables import flatten
    >>> s1, s2, s3 = 'abc', 'babc', 'abca'
    >>> marker = '$'
    >>> S = ''.join(chain(*flatten(zip([s1, s2, s3], '$$$'))))[:-1]
    >>> S
    'abc$babc$abca'
    >>> shortest_repeated_subsequence(S, ignore=['$'])
    'abc'
    """
    from sympy import oo

    # This is based on the longest_common_substring algorithm, for instance,
    # at
    # https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring#Python3. We
    # build a matrix to keep track of common subsequence of S with itself. For instance, for
    # the sequence [a, b, c, a, b], we would have

    #        a  b  c  a  b
    #   [[0, 0, 0, 0, 0, 0],
    # a  [0, 0, 0, 0, 1, 0],
    # b  [0, 0, 0, 0, 0, 2],
    # c  [0, 0, 0, 0, 0, 0],
    # a  [0, 1, 0, 0, 0, 0],
    # b  [0, 0, 2, 0, 0, 0]]

    # (the leading row and column of 0s are used so we can always reference
    # the previous row and column). Each entry in the matrix represents the
    # longest subsequence ending in those two entries. We additionally keep
    # the diagonal zero, so that we don't consider the entire sequence to be a
    # common subsequence with itself.

    m = [[0] * (1 + len(S)) for i in range(1 + len(S))]

    shortest, x_shortest = oo, 0
    for x in range(1, 1 + len(S)):
        if S[x-1] in ignore:
            continue
        for y in range(1, 1 + len(S)):
            if x == y:
                continue
            if S[x - 1] == S[y - 1]:
                m[x][y] = (m[x-1][y-1] + 1)

    # Now we traverse the matrix from the bottom up, keeping track diagonal
    # "runs" (representing the same subsequence), and look for the shortest
    # one that is >= 2 in length. If there are ties, the earliest one will be
    # returned, but it would be easy to extend this to return all shortest
    # subsequences.

    seen = {}
    for row in reversed(m):
        seen = {i-1:seen[i]-1 for i in seen if seen[i]}
        for x in range(len(row)):
            if row[x] >= 2 and x not in seen:
                if row[x] < shortest:
                    shortest = row[x]
                    x_shortest = x
                seen[x] = row[x]

    if shortest == oo:
        return None
    return S[x_shortest - shortest: x_shortest]
