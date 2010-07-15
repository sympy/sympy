"""Module with functions operating on Indexed, IndexedElement and Idx objects

    - Check shape conformance
    - Determine indices in resulting expression

    etc.
"""

from sympy.tensor.indexed import Idx, Indexed, IndexedElement


class IndexConformanceException(Exception):
    pass

def _remove_repeated(c_inds, nc_inds):
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

    return c_inds, nc_inds

def _get_indices_Mul(expr):
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

    return _remove_repeated(c_inds, nc_inds)


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
        else:
            raise NotImplementedError(
                    "No specialized function for type %s"%type(expr))
