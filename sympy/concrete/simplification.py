""" Change index / Reorder limits of Sums and Products"""
from sympy.concrete import Product, Sum
from sympy import Subs
from sympy import Symbol


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

    **Usage**

        change_index(expr, new, var) -> Here new is an expression of the form
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
    """
    limits = []
    old = list(new.free_symbols)[0]
    invert = not (new - old).is_number

    if var == None:
        var = old

    for limit in expr.limits:
        if limit[0] == old:
            if not invert:
                limits.append((var, limit[1]+ new - \
                    old, limit[2] + new - old))
            else:
                limits.append((var, (new + old) - limit[2], \
                    (new + old) - limit[1]))
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
        ChangeIndexError(expr, "Not Relevant / Implemented.")


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

    **Usage**
        reorder(expr, *arg) -> limits in the expr is reordered according to the
        list of tuples given by arg.

    Examples
    ========

    >>> from sympy.concrete.simplification import reorder
    >>> from sympy import Sum
    >>> from sympy.abc import x, y, z, a, b, c, d, m, n
    >>> reorder(Sum(x*y, (x, a, b), (y, c, d)), (0, 1))
    Sum(x*y, (y, c, d), (x, a, b))
    >>> reorder(Sum(x*y + z, (x, a, b), (z, m, n), (y, c, d)), (2, 0), (0, 1))
    Sum(x*y + z, (z, m, n), (y, c, d), (x, a, b))
    """
    temp = expr

    for r in arg:
        if len(r) != 2:
            ReorderError(r, "Invalid number of arguments.")

        temp = reorder_limit(temp, r[0], r[1])

    return temp


def reorder_limit(expr, x , y):

    var = [limit[0] for limit in expr.limits]
    limits = []

    limit_x = expr.limits[x]; limit_y = expr.limits[y]

    if set(limit_x[1].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_x[2].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_y[1].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_y[2].free_symbols).intersection(set(var)) == set([]):

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
                ReorderError(expr, "Not relevant / Not implemented.")

    else:
        ReorderError(expr, "Not relevant / Not implemented.")


def reverse_order(expr, *indexes):
    """
    Reverse the order of a limit in a Sum.

    Limits are simplified by the Karr's convention.

    **Usage**
        reverse_order(expr, *indexes) -> Reorder the limits in the expression
        specified by the indexes.

    Examples
    ========

    >>> from sympy.concrete.simplification import reverse_order
    >>> from sympy import Sum
    >>> from sympy.abc import x, y
    >>> reverse_order(Sum(x, (x, 0, 3)), 0)
    Sum(-x, (x, 1, 2))
    >>> reverse_order(Sum(x*y, (x, 1, 5), (y, 0, 6)), 0, 1)
    Sum(x*y, (x, 2, 4), (y, 1, 5))

    References
    ==========

    .. [1] http://dl.acm.org/citation.cfm?doid=322248.322255
    """
    e = 1
    limits = []

    if isinstance(expr, Sum):
        for i, limit in enumerate(expr.limits):
            l = limit
            if i in indexes:
                if isinstance(limit[1], Symbol) or isinstance(limit[2], Symbol):
                    e = e * -1
                    l = (limit[0], limit[1] + 1 , limit[2] - 1)
                else:
                    if limit[1] == limit[2] - 1:
                        e = 0
                        l = (limit[0], limit[2], limit[1])
                    elif limit[1] < limit[2] - 1:
                        e = e * -1
                        l = (limit[0], limit[1] + 1 , limit[2] - 1)
                    else:
                        l = (limit[0], limit[2], limit[1])

            limits.append(l)

        return Sum(expr.function * e, *limits)
