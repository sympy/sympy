from sympy import *


class Piecewise(Function):
    """
    > Container for a Piecewise Function
    -- Piecewise represents a piecewise function within the SymPy context.

    Example Usage:

    from sympy import Piecewise, oo, log, symbols
    x = symbols('x')
    p = Piecewise(x, [[(-oo, -1), 1], [(1, 2), x*x], [(3, oo), log(x)]])

    Explanation:
    - The first argument is the variable of the piecewise function.
    - The second argument is a Python list of sublists.
    - The first entry in the sublist is a tuple with the range of a function.
    - The second entry is the actual function over the specified range.
    """
    nargs=1

    def __new__(cls, *args, **options):
        if cls is Function:
            if len(args) == 1 and isinstance(args[0], str):
                return FunctionClass(Function, *args)
            else:
                print args
                print type(args[0])
                raise Exception("You need to specify exactly one string")
                args = map(Basic.sympify, args)
        return Basic.__new__(cls, *args, **options)

    def fdiff(self, arg=1):
        """ Returns the differentiated piecewise function """
        e = self.args[1]

        for i in range(0, len(self.args[1])):
            e[i][1] = diff(e[i][1], self.args[0])

    def canonize(self, arg, elems):
        arg = Basic.sympify(arg)
        for x in self.elems:
            cond = x[0]
            if isinstance(cond, (Basic.Number, Basic.Number)):
                if cond[0] <= arg and arg <= cond[1]:
                    return x[1]
            else:
                cond = Basic.sympify(cond)
                if isinstance(cond, Bool):
                    if cond:
                        return x[1]
