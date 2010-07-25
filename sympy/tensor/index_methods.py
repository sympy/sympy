"""Module with functions operating on Indexed, IndexedElement and Idx objects

    - Check shape conformance
    - Determine indices in resulting expression

    etc.
"""

from sympy.tensor.indexed import Idx, Indexed, IndexedElement


class IndexConformanceException(Exception):
    pass

def _remove_repeated(c_inds, nc_inds, return_dummies=False):
    """Removes summation indices from index sequences."""

    sum_index = {}
    for i in c_inds + nc_inds:
        if i in sum_index:
            sum_index[i] += 1
            assert sum_index[i] == 1, "Index %s repeated more than twice" % i
        else:
            sum_index[i] = 0

    c_inds = filter(lambda x: not sum_index[x], c_inds)
    nc_inds = filter(lambda x: not sum_index[x], nc_inds)

    if return_dummies:
        return c_inds, nc_inds, tuple([ i for i in sum_index if sum_index[i] ])
    else:
        return c_inds, nc_inds

def _get_indices_Mul(expr, return_dummies=False):
    """Determine the outer indices of a Mul object.

    returns all non-repeated indices.
    """

    junk, factors = expr.as_coeff_terms()
    inds = map(get_indices, factors)
    c_inds, nc_inds = zip(*inds)

    # sort commuting objects according to their indices
    c_inds = sorted(c_inds, key=hash)

    # flatten
    c_inds  = reduce(lambda x, y: x + y, c_inds)
    nc_inds = reduce(lambda x, y: x + y, nc_inds)

    return _remove_repeated(c_inds, nc_inds, return_dummies)


def _get_indices_Add(expr):
    """Determine outer indices of an Add object.

    In a sum, each term must have identical outer indices.  A valid expression
    could be (provided that x and y commutes):

        x(i)*y(j) - x(j)*y(i)

    But we do not allow expressions like:

        x(i)*y(j) - z(j)*z(j)

    FIXME: Add support for Numpy broadcasting
    """

    inds = map(get_indices, expr.args)
    if not all(map(lambda x: x == inds[0], inds[1:])):
        raise IndexConformanceException("Indices are not consistent: %s"%expr)
    return inds[0]

def get_indices(expr):
    """Determine the outer indices of expression `expr'

    By `outer' we mean indices that are not summation indices.  Returns two
    tuples with indices of commuting and non-commuting terms respectively.

    Examples
    ========

    >>> from sympy.tensor.index_methods import get_indices
    >>> from sympy import symbols
    >>> from sympy.tensor import Indexed, Idx
    >>> x, y, A = map(Indexed, ['x', 'y', 'A'])
    >>> i, j, a, z = symbols('i j a z', integer=True)


    The indices of the total expression is determined, Repeated indices imply a
    summation, for instance the trace of a matrix A:

    >>> get_indices(A(i, i))
    ((), ())

    In the case of many terms, the terms are required to have identical
    outer indices.  Else an IndexConformanceException is raised.

    >>> get_indices(x(i) + A(i, j)*y(j))
    ((i,), ())

    The concept of `outer' indices applies recursively, starting on the deepest
    level.  This implies that dummies inside parenthesis are assumed to be
    summed first, and there following expression is handled gracefully:

    >>> get_indices((x(i) + A(i, j)*y(j))*x(j))
    ((i, j), ())

    Note that if the indexed objects commute, the indices will be sorted so
    that the symbol of the stem should have no influence on the identified
    outer indices.

    >>> get_indices(x(i)*y(j)) == get_indices(x(j)*y(i))
    True

    But the order of indices on a particular Indexed should be intact:

    >>> get_indices(x(i, j))
    ((i, j), ())
    >>> get_indices(x(j, i))
    ((j, i), ())

    Exceptions
    ==========

    An IndexConformanceException means that the terms ar not compatible, e.g.

    >>> get_indices(x(i) + y(j))                #doctest: +SKIP
            (...)
    IndexConformanceException: Indices are not consistent: x(i) + y(j)


    """
    # We call ourself recursively to determine indices of sub expressions.

    # break recursion
    if isinstance(expr, IndexedElement):
        if expr.is_commutative:
            c = expr.indices
            nc = tuple()
        else:
            c = tuple()
            nc = expr.indices
        return _remove_repeated(c, nc)
    elif expr.is_Atom:
        return tuple(), tuple()

    # recurse via specialized functions
    else:
        if expr.is_Mul:
            return _get_indices_Mul(expr)
        elif expr.is_Add:
            return _get_indices_Add(expr)

        # this test is expensive, so it should be at the end
        elif not expr.has(IndexedElement):
            return tuple(), tuple()
        else:
            raise NotImplementedError(
                    "FIXME: No specialized handling of type %s"%type(expr))

def get_contraction_structure(expr):
    """Determine dummy indices of expression `expr' and describe structure of the expression

    By `dummy' we mean indices that are summation indices.

    The stucture of the expression is determined and returned as follows:

    1) The terms of a conforming summation are returned as a dict where the keys
    are summation indices and the values are terms for which the dummies are
    relevant.

    2) If there are nested Add objects, we recurse to determine summation
    indices for the the deeper terms. The resulting dict is returned as a value
    in the dictionary, and the Add expression is the corresponding key.

    Examples
    ========

    >>> from sympy.tensor.index_methods import get_contraction_structure
    >>> from sympy import symbols
    >>> from sympy.tensor import Indexed, Idx
    >>> x, y, A = map(Indexed, ['x', 'y', 'A'])
    >>> i, j, k, l = symbols('i j k l', integer=True)
    >>> get_contraction_structure(x(i)*y(i) + A(j, j))
    {(i,): set([x(i)*y(i)]), (j,): set([A(j, j)])}
    >>> get_contraction_structure(x(i)*y(j))
    {None: set([x(i)*y(j)])}

    A nested Add object is returned as a nested dictionary.  The term
    containing the parenthesis is used as the key, and it stores the dictionary
    resulting from a recursive call on the Add expression.

    >>> d = get_contraction_structure(x(i)*(y(i) + A(i, j)*x(j)))
    >>> sorted(d.keys())
    [(i,), (x(j)*A(i, j) + y(i))*x(i)]
    >>> d[(Idx(i),)]
    set([(x(j)*A(i, j) + y(i))*x(i)])
    >>> d[(x(j)*A(i, j) + y(i))*x(i)]
    [{None: set([y(i)]), (j,): set([x(j)*A(i, j)])}]

    Note that the presence of expressions among the dictinary keys indicates a
    factorization of the array contraction.  The summation in the deepest
    nested level must be calculated first so that the external contraction can access
    the resulting array with index j.

    """

    # We call ourself recursively to inspect sub expressions.

    if isinstance(expr, IndexedElement):
        c = expr.indices
        nc = tuple()
        junk, junk, key = _remove_repeated(c, nc, return_dummies=True)
        return {key or None: set([expr])}
    elif expr.is_Atom:
        return {None: set([expr])}
    elif expr.is_Mul:
        junk, junk, key = _get_indices_Mul(expr, return_dummies=True)
        result = {key or None: set([expr])}
        # recurse if we have any Add objects
        addfactors = filter(lambda x: x.is_Add, expr.args)
        if addfactors:
            result[expr] = []
            for factor in addfactors:
                d = get_contraction_structure(factor)
                result[expr].append(d)
        return result
    elif expr.is_Add:
        # Note: we just collect all terms with identical summation indices, We
        # do nothing to identify equivalent terms here, as this would require
        # substitutions or pattern matching in expressions of unknown
        # complexity.
        result = {}
        for term in expr.args:
            # recurse on every term
            d = get_contraction_structure(term)
            for key in d:
                if key in result:
                    result[key] |= d[key]
                else:
                    result[key] = d[key]
        return result

    # this test is expensive, so it should be at the end
    elif not expr.has(IndexedElement):
        return {None: set([expr])}
    else:
        raise NotImplementedError(
                "FIXME: No specialized handling of type %s"%type(expr))
