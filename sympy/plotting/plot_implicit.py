"""Implicit plotting module for SymPy

The module implements a data series called ImplicitSeries which is used by
``Plot`` class to plot implicit plots for different backends. The module
implements plotting using interval arithmetic.

See Also
========
sympy.plotting.plot

References
==========
- Jeffrey Allen Tupper. Reliable Two-Dimensional Graphing Methods for
Mathematical Formulae with Two Free Variables.

- Jeffrey Allen Tupper. Graphing Equations with Generalized Interval
Arithmetic. Master's thesis. University of Toronto, 1996

"""

from plot import BaseSeries, Plot
from experimental_lambdify import experimental_lambdify, vectorized_lambdify
from intervalmath import interval
from sympy.core.relational import (Equality, GreaterThan, LessThan,
                Relational, StrictLessThan, StrictGreaterThan)
from sympy import Eq, Tuple, sympify, Expr
from sympy.external import import_module
from sympy.core.compatibility import set_union
from sympy.logic.boolalg import BooleanFunction

np = import_module('numpy')


class ImplicitSeries(BaseSeries):
    """ Representation for Implicit plot """
    is_implicit = True

    def __init__(self, expr, var_start_end_x, var_start_end_y):
        super(ImplicitSeries, self).__init__()
        self.expr = sympify(expr)
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])
        self.get_points = self.get_meshes

    def __str__(self):
        return ('Implicit equation: %s for '
                '%s over %s and %s over %s') % (
                str(self.expr),
                str(self.var_x),
                str((self.start_x, self.end_x)),
                str(self.var_y),
                str((self.start_y, self.end_y)))

    def get_meshes(self):
        WIDTH = 512   #TODO: Add as an attribute which can be changed
        HEIGHT = 512  #TODO: Add as an attribute which can be changed
        is_equality = isinstance(self.expr, Equality)
        func = experimental_lambdify((self.var_x, self.var_y), self.expr, use_interval=True)
        xinterval = interval(self.start_x, self.end_x)
        yinterval = interval(self.start_y, self.end_y)
        #TODO: Have a fall back algorithm for functions that are not implemented
        try:
            temp = func(xinterval, yinterval)
        except AttributeError:
            raise NotImplementedError("Some functions in the expression %s are "
                                    "not supported. Kindly report this as a bug" % str(self.expr))
        #contour array, acts like a bitmap
        contour = np.zeros((WIDTH, HEIGHT))
        k = 5 #Best case by trial and error
        interval_list = []

        xsample = np.linspace(self.start_x, self.end_x, WIDTH / 2**k + 1)
        ysample = np.linspace(self.start_y, self.end_y, HEIGHT / 2**k + 1)

        xinter = [interval(x1, x2) for x1, x2 in zip(xsample[:-1], xsample[1:])]
        yinter = [interval(y1, y2) for y1, y2 in zip(ysample[:-1], ysample[1:])]
        interval_list = [[x, y] for x in xinter for y in yinter]
        plot_list = []

        #recursive call refinepixels which subdivides the intervals which are neither
        # True nor False according to the expression.
        def refine_pixels(interval_list):
            temp_interval_list = []
            plot_list = []
            for intervals in interval_list:

                #Convert the array indices to x and y values
                intervalx = intervals[0]
                intervaly = intervals[1]
                func_eval = func(intervalx, intervaly)
                #The expression is valid in the interval. Change the contour array
                #values to 1.
                if func_eval[1] is False or func_eval[0] is False:
                    pass
                elif func_eval == (True, True) and not is_equality:
                    plot_list.append([intervalx, intervaly])
                elif func_eval[1] is None or func_eval[0] is None \
                        or is_equality:
                    #Subdivide
                    avgx = intervalx.mid
                    avgy = intervaly.mid
                    a = interval(intervalx.start, avgx)
                    b = interval(avgx, intervalx.end)
                    c = interval(intervaly.start, avgy)
                    d = interval(avgy, intervaly.end)
                    temp_interval_list.append([a, c])
                    temp_interval_list.append([a, d])
                    temp_interval_list.append([b, c])
                    temp_interval_list.append([b, d])
            return temp_interval_list, plot_list

        while k >= 0 and len(interval_list):
            interval_list, plot_list_temp = refine_pixels(interval_list)
            plot_list.extend(plot_list_temp)
            k = k - 1
        #Check whether the expression represents an equality
        #If it represents an equality, then none of the intervals
        #would have satisfied the expression due to floating point
        #differences. Add all the undecided values to the plot.
        if isinstance(self.expr, (Equality, GreaterThan, LessThan)):
            for intervals in interval_list:
                intervalx = intervals[0]
                intervaly = intervals[1]
                func_eval = func(intervalx, intervaly)
                if func_eval[1] and func_eval[0] is not False:
                    plot_list.append([intervalx, intervaly])
        xvals, yvals = _matplotlib_list(plot_list)
        return xvals, yvals


def plot_implicit(expr, var_start_end_x, var_start_end_y, **kwargs):
    """A plot function to plot implicit equations / inequalities.

    The input arguments are:
    expr : The equation / inequality that is to be plotted.
    var_start_end_x: A tuple of length 3, with the first element representing
    the variable and the next two elements representing the range
    var_start_end_y: A tuple of length 3, with the first element representing
    the variable and the next two elements representing the range

    Examples:
    ---------

    Plot expressions:
    >>> from sympy import plot_implicit, cos, sin, symbols, Eq
    >>> x, y = symbols('x y')
    >>> p1 = plot_implicit(Eq(y, x ** 2), (x, -5, 5), (y, -5, 5)) #doctest: +SKIP
    >>> p2 = plot_implicit(Eq(x ** 2 + y ** 2, 3), (x, -3, 3), (y, -3, 3)) #doctest: +SKIP
    >>> p3 = plot_implicit(y ** 2 < x ** 3 - x, (x, -4, 4), (y, -4, 4)) #doctest: +SKIP
    >>> p4 = plot_implicit(y > sin(x), (x, -5, 5), (y, -2, 2)) #doctest: +SKIP
    """
    #TODO: Add a global variable show = False for test runner
    assert isinstance(expr, Expr)

    is_equal = False #Represents whether the expression contains an equality,
                     #GreaterThan or LessThan
    arg_list = []
    def arg_expand(bool_expr):
        """
        Recursively expands the arguments of an Boolean Function
        """
        for arg in bool_expr.args:
            if isinstance(arg, BooleanFunction):
                arg_expand(arg)
            elif isinstance(arg, Relational):
                arg_list.append(arg)

    if isinstance(expr, BooleanFunction):
        arg_expand(expr)

    #Check whether there is an equality in the expression provided.
        if any(isinstance(e, (Equality, GreaterThan, LessThan))
                            for e in arg_list):
            is_equal = True

    elif not isinstance(expr, Relational):
        expr = Eq(expr, 0)
    free_symbols = expr.free_symbols
    range_symbols = set([var_start_end_x[0], var_start_end_y[0]])
    symbols = set_union(free_symbols, range_symbols)
    if len(symbols) > 2:
        raise NotImplementedError("Implicit plotting is not implemented for "
                                    "more than 2 variables")
    series_argument = ImplicitSeries(expr, var_start_end_x, var_start_end_y)
    show = kwargs.pop('show', True)
    p = Plot(series_argument, **kwargs)
    if show:
        p.show()
    return p

def _matplotlib_list(interval_list):
    """
    Returns lists for matplotlib ``fill`` command from a list of bounding
    rectangular intervals
    """
    xlist = []
    ylist = []
    for intervals in interval_list:
        intervalx = intervals[0]
        intervaly = intervals[1]
        xlist.extend([intervalx.start, intervalx.start, intervalx.end, intervalx.end, None])
        ylist.extend([intervaly.start, intervaly.end, intervaly.end, intervaly.start, None])
    return xlist, ylist
