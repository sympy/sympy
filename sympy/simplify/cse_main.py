""" Tools for doing common subexpression elimination.
"""
from __future__ import print_function, division

import difflib

from sympy.core import Basic, Mul, Add, Pow, sympify, Tuple, Symbol
from sympy.core.singleton import S
from sympy.core.basic import preorder_traversal
from sympy.core.function import _coeff_isneg
from sympy.core.exprtools import factor_terms
from sympy.core.compatibility import iterable, xrange
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
    y == 2
    x == y + 1

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
    if optimizations is None:
        optimizations = cse_optimizations
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
    --------
    >>> from sympy.simplify.cse_main import opt_cse
    >>> from sympy.abc import x
    >>> opt_subs = opt_cse([x**-2])
    >>> print(opt_subs)
    {x**(-2): 1/(x**2)}
    """
    from sympy.matrices import Matrix

    opt_subs = dict()

    adds = set()
    muls = set()

    seen_subexp = set()

    def _find_opts(expr):

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

        if expr.is_Mul:
            muls.add(expr)

        elif expr.is_Add:
            adds.add(expr)

        elif expr.is_Pow:
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
        for i in xrange(len(func_args)):
            for j in xrange(i + 1, len(func_args)):
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

                    for k in xrange(j + 1, len(func_args)):
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
            c_mul = Mul(*c)
            if nc:
                opt_subs[m] = Mul(c_mul, Mul(*nc), evaluate=False)
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
    from sympy.matrices import Matrix

    if opt_subs is None:
        opt_subs = dict()

    ## Find repeated sub-expressions

    to_eliminate = set()

    seen_subexp = set()

    def _find_repeated(expr):
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

        if expr.is_Atom:
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
            if expr.is_Mul:
                c, nc = expr.args_cnc()
                args = list(ordered(c)) + nc
            elif expr.is_Add:
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
    """
    from sympy.matrices import Matrix

    # Handle the case if just one expression was passed.
    if isinstance(exprs, Basic):
        exprs = [exprs]

    if optimizations is None:
        optimizations = list()
    elif optimizations == 'basic':
        optimizations = basic_optimizations

    # Preprocess the expressions to give us better optimization opportunities.
    reduced_exprs = [preprocess_for_cse(e, optimizations) for e in exprs]

    excluded_symbols = set.union(*[expr.atoms(Symbol)
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
    for i, (sym, subtree) in enumerate(replacements):
        subtree = postprocess_for_cse(subtree, optimizations)
        replacements[i] = (sym, subtree)
    reduced_exprs = [postprocess_for_cse(e, optimizations)
                     for e in reduced_exprs]

    if isinstance(exprs, Matrix):
        reduced_exprs = [Matrix(exprs.rows, exprs.cols, reduced_exprs)]
    if postprocess is None:
        return replacements, reduced_exprs
    return postprocess(replacements, reduced_exprs)
