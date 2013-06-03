""" Change index / Reorder limits of Sums and Products"""
from sympy.concrete import Product, Sum
from sympy import Subs

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
    Change index of a Sum or Product. Changes of the form x --> ax + b
    are allowed. Here a is 1 or -1.
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
    Reorder limits in a expression lik a Sum or a Product.
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
