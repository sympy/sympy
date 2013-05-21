""" Change index / Reorder of Sums and Products"""
from sympy.concrete import Product, Sum
from sympy import Subs

class ShiftError(NotImplementedError):
    """
    Exception raised when something other than Sum or Product
    is passed to the function change_index().
    """
    def __init__(self, expr, msg):
        super(ShiftError, self).__init__(
            "%s Shift could not be computed: %s." % (expr, msg))


def change_index(expr, new, var=None):

    if isinstance(expr, Sum) or isinstance(expr, Product):
        return sum_prod_change_index(expr, new, var)
    else:
        ShiftError(expr, "Not Relevant / Implemented.")


def sum_prod_change_index(expr, new, var=None):

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


class ReorderError(NotImplementedError):
    """
    Exception raised when trying to reorder dependent limits.
    """
    def __init__(self, expr, msg):
        super(ReorderError, self).__init__(
            "%s could not be reordered: %s." % (expr, msg))


def reorder(expr, x , y):

    var = [limit[0] for limit in expr.limits]
    limits = []

    if var.count(x) != 1 or var.count(y) != 1:
        return expr

    index_x = var.index(x); index_y = var.index(y)
    limit_x = expr.limits[index_x]; limit_y = expr.limits[index_y]

    if set(limit_x[1].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_x[2].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_y[1].free_symbols).intersection(set(var)) == set([]) \
        and set(limit_y[2].free_symbols).intersection(set(var)) == set([]):

            for i, limit in enumerate(expr.limits):
                if i == index_x:
                    limits.append(limit_y)
                elif i == index_y:
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
