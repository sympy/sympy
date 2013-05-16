""" Change index of Sums and Products"""
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


def change_index(expr, old, new, var=None):

    if isinstance(expr, Sum) or isinstance(expr, Product):
        return sum_prod_change_index(expr, old, new, var)
    else:
        ShiftError(expr, "Not Relevant / Implemented.")


def sum_prod_change_index(expr, old, new, var=None):

    limits = []
    invert = not (new - old).is_number

    if var == None:
        var = old

    for i in range(len(expr.limits)):
        if expr.limits[i][0] == old:
            if not invert:
                limits.append((var, expr.limits[i][1]+ new - \
                    old, expr.limits[i][2] + new - old))
            else:
                limits.append((var, (new + old) - expr.limits[i][2], \
                    (new + old) - expr.limits[i][1]))
        else:
            limits.append(expr.limits[i])

    if not invert:
        function = Subs(expr.function, old, var - (new - old)).doit()
    else:
        function = Subs(expr.function, old, (new + old) - var).doit()

    if isinstance(expr, Sum):
        return Sum(function, *tuple(limits))
    elif isinstance(expr, Product):
        return Product(function, *tuple(limits))
