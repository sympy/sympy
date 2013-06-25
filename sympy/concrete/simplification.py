""" Change index / Reorder / Reverse order of limits of Sums and Products"""
from sympy.concrete import Product, Sum
from sympy import Subs
from sympy import Symbol
from sympy import Integer


class ChangeIndexError(NotImplementedError):
    """
    Exception raised when something other than Sum or Product
    is passed to the function change_index().
    """
    def __init__(self, expr, msg):
        super(ShiftError, self).__init__(
            "%s Shift could not be computed: %s." % (expr, msg))


def change_index(expr, new, var=None):
    """
    Change index of a Sum or Product.

    Changes of the index of form x --> ax + b are allowed. Here a is 1 or -1.
    New variable to be used after the change of index should also be specified.

    Usage
    =====

        change_index(expr, new, var=None) -> Here new is an expression of the form
        ax + b where a = 1 or -1. var is an optional argument. If var is given,
        ax + b is replaced by var.

    Examples
    ========

    >>> from sympy.concrete.simplification import change_index
    >>> from sympy import Sum
    >>> from sympy.abc import x, y, a, b
    >>> change_index(Sum(x, (x, a, b)), x + 1, y)
    Sum(y - 1, (y, a + 1, b + 1))
    >>> change_index(Sum(x, (x, a, b)), -x - 1)
    Sum(-x - 1, (x, -b - 1, -a - 1))
    >>> from sympy.abc import v
    >>> change_index(Sum(x, (x, a, b)), x + v)
    Sum(-v + x, (x, a + v, b + v))
    >>> change_index(Sum(x, (x, a, b)), -x - v)
    Sum(-v - x, (x, -b - v, -a - v))
    """
    limits = []

    for x in list(new.free_symbols):
        for i, limit in enumerate(expr.limits):
            if x == limit[0]:
                old = x

    invert = not (new + old).has(old)

    if var is None:
        var = old

    for limit in expr.limits:
        if limit[0] == old:
            if not invert:
                limits.append((var, limit[1]+ new - old, limit[2] + new - old))
            else:
                limits.append((var, (new + old) - limit[2], (new + old) - limit[1]))
        else:
            limits.append(limit)

    if not invert:
        function = Subs(expr.function, old, var - (new - old)).doit()
    else:
        function = Subs(expr.function, old, (new + old) - var).doit()

    if isinstance(expr, Sum):
        return Sum(function, *tuple(limits))
    elif isinstance(expr, Product):
        return Product(function, *tuple(limits))
    else:
        raise ChangeIndexError(expr, "change_index only implemented for Sum/Product.")


class ReorderError(NotImplementedError):
    """
    Exception raised when trying to reorder dependent limits.
    """
    def __init__(self, expr, msg):
        super(ReorderError, self).__init__(
            "%s could not be reordered: %s." % (expr, msg))


def reorder(expr, *arg):
    """
    Reorder limits in a expression like a Sum or a Product.

    Usage
    =====

        reorder(expr, *arg) -> limits in the expr is reordered according to the
        list of tuples given by arg. These tuples can be tuples of indices or tuples
        of index variables or can involve both.

    Examples
    ========

    >>> from sympy.concrete.simplification import reorder
    >>> from sympy import Sum
    >>> from sympy.abc import x, y, z, a, b, c, d, m, n
    >>> reorder(Sum(x*y, (x, a, b), (y, c, d)), (0, 1))
    Sum(x*y, (y, c, d), (x, a, b))
    >>> reorder(Sum(x*y + z, (x, a, b), (z, m, n), (y, c, d)), (2, 0), (0, 1))
    Sum(x*y + z, (z, m, n), (y, c, d), (x, a, b))
    >>> reorder(Sum(x*y, (x, a, b), (y, c, d)), (y, x))
    Sum(x*y, (y, c, d), (x, a, b))
    >>> reorder(Sum(x*y, (x, a, b), (y, c, d)), (y, 0))
    Sum(x*y, (y, c, d), (x, a, b))
    """
    temp = expr

    for r in arg:
        if len(r) != 2:
            raise ReorderError(r, "Invalid number of arguments.")

        index1 = r[0]
        index2 = r[1]

        if not isinstance(r[0], int):
            index1 = index(expr, r[0])
        if not isinstance(r[1], int):
            index2 = index(expr, r[1])
        if index1 == -1 or index2 == -1:
            raise ReorderError(r, "Instances of variable not equal to one.")

        temp = reorder_limit(temp, index1, index2)

    return temp


def index(expr, x):
    """
    Return the index of a limit variable.

    Usage
    =====

        index(expr, x) -> return the index of the limit variable x in the limits
        of expr.

    Examples
    ========

    >>> from sympy.concrete.simplification import index
    >>> from sympy.abc import x, y, a, b
    >>> from sympy import Sum
    >>> index(Sum(x*y, (x, 1, 5), (y, a, b)), x)
    0
    >>> index(Sum(x*y, (x, 1, 5), (y, a, b)), y)
    1
    """
    if isinstance(expr, Sum) or isinstance(expr, Product):
        variables = [limit[0] for limit in expr.limits]

        if variables.count(x) != 1:
            return -1
        else:
            return variables.index(x)


def reorder_limit(expr, x , y):
    """
    Reorder the two limits corresponds to indices x and y.

    Usage
    =====

        reorder_limit(expr, x, y) -> expr is either a Sum or a Product.
        x and y are integers corresponding to the indices of the limits
        which are to be reordered.

    Examples
    ========

    >>> from sympy.concrete.simplification import reorder_limit
    >>> from sympy.abc import x, y, a, b
    >>> from sympy import Sum
    >>> reorder_limit(Sum(x*y, (x, 1, 6), (y, 2, 8)), 0, 1)
    Sum(x*y, (y, 2, 8), (x, 1, 6))
    >>> reorder_limit(Sum(x*y, (x, a, b), (x, 2, 8)), 0, 1)
    Sum(x*y, (x, 2, 8), (x, a, b))
    """
    var = [limit[0] for limit in expr.limits]
    limits = []

    limit_x = expr.limits[x]
    limit_y = expr.limits[y]

    if (set(limit_x[1].free_symbols).intersection(set(var)) == set([])
        and set(limit_x[2].free_symbols).intersection(set(var)) == set([])
        and set(limit_y[1].free_symbols).intersection(set(var)) == set([])
        and set(limit_y[2].free_symbols).intersection(set(var)) == set([])):

            for i, limit in enumerate(expr.limits):
                if i == x:
                    limits.append(limit_y)
                elif i == y:
                    limits.append(limit_x)
                else:
                    limits.append(limit)

            if isinstance(expr, Sum):
                return Sum(expr.function, *limits)
            elif isinstance(expr, Product):
                return Product(expr.function, *limits)
            else:
                raise ReorderError(expr, "reorder only implemented for Sum/Product.")

    else:
        raise ReorderError(expr, "reorder only implemented for Sum/Product.")


class ReverseOrderError(NotImplementedError):
    """
    Exception raised when trying to reorder dependent limits.
    """
    def __init__(self, expr, msg):
        super(ReverseOrderError, self).__init__(
            "%s limits could not be reversed: %s." % (expr, msg))


def reverse_order(expr, *indices):
    """
    Reverse the order of a limit in a Sum.

    Usage
    =====

        reverse_order(expr, *indices) -> Reverse the limits in the expression
        specified by the indexes. Here indices can be indices of limits to
        be reversed or index variables corresponding to those limits.

    Examples
    ========

    >>> from sympy.concrete.simplification import reverse_order
    >>> from sympy import Sum
    >>> from sympy.abc import x, y, a, b
    >>> reverse_order(Sum(x, (x, 0, 3)), 0)
    Sum(-x, (x, 4, -1))
    >>> reverse_order(Sum(x*y, (x, 1, 5), (y, 0, 6)), 0, 1)
    Sum(x*y, (x, 6, 0), (y, 7, -1))
    >>> reverse_order(Sum(x, (x, a, b)), 0)
    Sum(-x, (x, b + 1, a - 1))
    >>> reverse_order(Sum(x, (x, a, b)), x)
    Sum(-x, (x, b + 1, a - 1))
    >>> reverse_order(Sum(x*y, (x, a, b), (y, 2, 5)), x, 1)
    Sum(x*y, (x, b + 1, a - 1), (y, 6, 1))
    >>> reverse_order(Sum(x*y, (x, a, b), (y, 2, 5)), y, x)
    Sum(x*y, (x, b + 1, a - 1), (y, 6, 1))

    References
    ==========

    .. [1] Michael Karr, "Summation in Finite Terms", Journal of the ACM,
        Volume 28 Issue 2, April 1981, Pages 305-350
        http://dl.acm.org/citation.cfm?doid=322248.322255
    """
    e = 1
    limits = []
    l_indices = list(indices)

    for i, indx in enumerate(l_indices):
        if not isinstance(indx, int):
            if index(expr, indx) == -1:
                raise ReverseOrderError(expr, "Instances of variable not equal to one.")
            l_indices[i] = index(expr, indx)

    if isinstance(expr, Sum):
        for i, limit in enumerate(expr.limits):
            l = limit
            if i in l_indices:
                if limit[1] != limit[2]:
                    e = -e
                    l = (limit[0], limit[2] + 1 , limit[1] - 1)

            limits.append(l)

        return Sum(e * expr.function, *limits)
    else:
        return expr
