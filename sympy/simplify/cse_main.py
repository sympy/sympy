""" Tools for doing common subexpression elimination.
"""
from __future__ import print_function, division

from sympy.core import Basic, Mul, Add, Pow, sympify, Symbol, Tuple, igcd
from sympy.core.numbers import Integer
from sympy.core.singleton import S
from sympy.core.function import _coeff_isneg
from sympy.core.exprtools import factor_terms
from sympy.core.compatibility import iterable, range, as_int
from sympy.utilities.iterables import filter_symbols, \
    numbered_symbols, sift, topological_sort, ordered, subsets

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


def pairwise_most_common(sets):
    """Return a list of `(s, L)` tuples where `s` is the largest subset
    of elements that appear in pairs of sets given by `sets` and `L`
    is a list of tuples giving the indices of the pairs of sets in
    which those elements appeared. All `s` will be of the same length.

    Examples
    ========

    >>> from sympy.simplify.cse_main import pairwise_most_common
    >>> pairwise_most_common((
    ...     {1,2,3},
    ...     {1,3,5},
    ...     {1,2,3,4,5},
    ...     {1,2,3,6}))
    [({1, 3, 5}, [(1, 2)]), ({1, 2, 3}, [(0, 2), (0, 3), (2, 3)])]
    >>>
    """
    from sympy.utilities.iterables import subsets
    from collections import defaultdict
    most = -1
    for i, j in subsets(list(range(len(sets))), 2):
        com = sets[i] & sets[j]
        if com and len(com) > most:
            best = defaultdict(list)
            best_keys = []
            most = len(com)
        if len(com) == most:
            if com not in best_keys:
                best_keys.append(com)
            best[best_keys.index(com)].append((i,j))
    if most == -1:
        return []
    for k in range(len(best)):
        best_keys[k] = (best_keys[k], best[k])
    best_keys.sort(key=lambda x: len(x[1]))
    return best_keys


def opt_cse(exprs, order='canonical', verbose=False):
    """Find optimization opportunities in Adds, Muls, Pows and negative
    coefficient Muls

    Parameters
    ----------
    exprs : list of sympy expressions
        The expressions to optimize.
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. For large
        expressions where speed is a concern, use the setting order='none'.
    verbose : bool
        Print debug information (default=False)

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

        if Func is Mul:
            F = Pow
            meth = 'as_powers_dict'
            from sympy.core.add import _addsort as inplace_sorter
        elif Func is Add:
            F = Mul
            meth = 'as_coefficients_dict'
            from sympy.core.mul import _mulsort as inplace_sorter
        else:
            assert None  # expected Mul or Add

        # ----------------- helpers ---------------------------
        def ufunc(*args):
            # return a well formed unevaluated function from the args
            # SHARES Func, inplace_sorter
            args = list(args)
            inplace_sorter(args)
            return Func(*args, evaluate=False)

        def as_dict(e):
            # creates a dictionary of the expression using either
            # as_coefficients_dict or as_powers_dict, depending on Func
            # SHARES meth
            d = getattr(e, meth, lambda: {a: S.One for a in e.args})()
            for k in list(d.keys()):
                try:
                    as_int(d[k])
                except ValueError:
                    d[F(k, d.pop(k))] = S.One
            return d

        def from_dict(d):
            # build expression from dict from
            # as_coefficients_dict or as_powers_dict
            # SHARES F
            return ufunc(*[F(k, v) for k, v in d.items()])

        def update(k):
            # updates all of the info associated with k using
            # the com_dict: func_dicts, func_args, opt_subs
            # returns True if all values were updated, else None
            # SHARES com_dict, com_func, func_dicts, func_args,
            #        opt_subs, funcs, verbose
            for di in com_dict:
                # don't allow a sign to change
                if com_dict[di] > func_dicts[k][di]:
                    return
            # remove it
            if Func is Add:
                take = min(func_dicts[k][i] for i in com_dict)
                com_func_take = Mul(take, from_dict(com_dict), evaluate=False)
            else:
                take = igcd(*[func_dicts[k][i] for i in com_dict])
                com_func_take = Pow(from_dict(com_dict), take, evaluate=False)
            for di in com_dict:
                func_dicts[k][di] -= take*com_dict[di]
            # compute the remaining expression
            rem = from_dict(func_dicts[k])
            # reject hollow change, e.g extracting x + 1 from x + 3
            if Func is Add and rem and rem.is_Integer and 1 in com_dict:
                return
            if verbose:
                print('\nfunc %s (%s) \ncontains %s \nas %s \nleaving %s' %
                    (funcs[k], func_dicts[k], com_func, com_func_take, rem))
            # recompute the dict since some keys may now
            # have corresponding values of 0; one could
            # keep track of which ones went to zero but
            # this seems cleaner
            func_dicts[k] = as_dict(rem)
            # update associated info
            func_dicts[k][com_func] = take
            func_args[k] = set(func_dicts[k])
            # keep the constant separate from the remaining
            # part of the expression, e.g. 2*(a*b) rather than 2*a*b
            opt_subs[funcs[k]] = ufunc(rem, com_func_take)
            # everything was updated
            return True

        def get_copy(i):
            return [func_dicts[i].copy(), func_args[i].copy(), funcs[i], i]

        def restore(dafi):
            i = dafi.pop()
            func_dicts[i], func_args[i], funcs[i] = dafi

        # ----------------- end helpers -----------------------

        func_dicts = [as_dict(f) for f in funcs]
        func_args = [set(d) for d in func_dicts]
        while True:
            hit = pairwise_most_common(func_args)
            if not hit or len(hit[0][0]) <= 1:
                break
            changed = False
            for com_args, ij in hit:
                take = len(com_args)
                ALL = list(ordered(com_args))
                while take >= 2:
                    for com_args in subsets(ALL, take):
                        com_func = Func(*com_args)
                        com_dict = as_dict(com_func)
                        for i, j in ij:
                            dafi = None
                            if com_func != funcs[i]:
                                dafi = get_copy(i)
                                ch = update(i)
                                if not ch:
                                    restore(dafi)
                                    continue
                            if com_func != funcs[j]:
                                dafj = get_copy(j)
                                ch = update(j)
                                if not ch:
                                    if dafi is not None:
                                        restore(dafi)
                                    restore(dafj)
                                    continue
                            changed = True
                        if changed:
                            break
                    else:
                        take -= 1
                        continue
                    break
                else:
                    continue
                break
            if not changed:
                break

    # split muls into commutative
    commutative_muls = set()
    for m in muls:
        c, nc = m.args_cnc(cset=True)
        if c:
            c_mul = m.func(*c)
            if nc:
                opt_subs[m] = m.func(c_mul, m.func(*nc), evaluate=False)
            if len(c) > 1:
                commutative_muls.add(c_mul)

    _match_common_args(Add, adds)
    _match_common_args(Mul, commutative_muls)

    return opt_subs


def tree_cse(exprs, symbols, opt_subs=None, order='canonical', ignore=()):
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
    ignore : iterable of Symbols
        Substitutions containing any Symbol from ``ignore`` will be ignored.
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
                for ign in ignore:
                    if ign in expr.free_symbols:
                        break
                else:
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

    # don't allow hollow nesting
    # e.g if p = [b + 2*d + e + f, b + 2*d + f + g, a + c + d + f + g]
    # and R, C = cse(p) then
    #     R = [(x0, d + f), (x1, b + d)]
    #     C = [e + x0 + x1, g + x0 + x1, a + c + d + f + g]
    # but the args of C[-1] should not be `(a + c, d + f + g)`
    nested = [[i for i in f.args if isinstance(i, f.func)] for f in exprs]
    for i in range(len(exprs)):
        F = reduced_exprs[i].func
        if not (F is Mul or F is Add):
            continue
        nested = [a for a in exprs[i].args if isinstance(a, F)]
        args = []
        for a in reduced_exprs[i].args:
            if isinstance(a, F):
                for ai in a.args:
                    if isinstance(ai, F) and ai not in nested:
                        args.extend(ai.args)
                    else:
                        args.append(ai)
            else:
                args.append(a)
        reduced_exprs[i] = F(*args)

    return replacements, reduced_exprs


def cse(exprs, symbols=None, optimizations=None, postprocess=None,
        order='canonical', ignore=()):
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
    ignore : iterable of Symbols
        Substitutions containing any Symbol from ``ignore`` will be ignored.

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
    ([(x0, w + y + z)], [x0*(x + x0)/(w + x)**3])

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

    The user may disallow substitutions containing certain symbols:
    >>> cse([y**2*(x + 1), 3*y**2*(x + 1)], ignore=(y,))
    ([(x0, x + 1)], [x0*y**2, 3*x0*y**2])

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
                                           order, ignore)

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
