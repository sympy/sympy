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
from experimental_lambdify import experimental_lambdify
from intervalmath import interval
from sympy.core.relational import Equality, GreaterThan, LessThan
from sympy.external import import_module
from sympy import sympify, Expr
from sympy.core.compatibility import set_union

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
        #TODO: Have a fallback algorithm for functions that are not implemented
        try:
            temp = func(xinterval, yinterval)
        except AttributeError:
            raise NotImplementedError("Some functions in the expression %s are \
                                    not supported. Kindly report this as a bug" % str(self.expr))
        #contour array, acts like a bitmap
        contour = np.zeros((WIDTH, HEIGHT))
        k = 6
        interval_list = []

        xsample = [x for x in xrange(0, WIDTH + 1, 2 ** k)]
        ysample = [y for y in xrange(0, HEIGHT + 1, 2 ** k)]
        xinter = [interval(x1, x2) for x1, x2 in zip(xsample[:-1], xsample[1:])]
        yinter = [interval(y1, y2) for y1, y2 in zip(ysample[:-1], ysample[1:])]
        interval_list = [[x, y] for x in xinter for y in yinter]

        #recursive call refinepixels which subdivides the intervals which are neither
        # True nor False according to the expression.
        def refine_pixels(interval_list):
            temp_interval_list = []
            for intervals in interval_list:
                xindexa, xindexb = intervals[0].start, intervals[0].end
                yindexa, yindexb = intervals[1].start, intervals[1].end

                #Convert the array indices to x and y values
                xa = self.start_x + xindexa * (self.end_x - self.start_x) / WIDTH
                xb = self.start_x + xindexb * (self.end_x - self.start_x) / WIDTH
                ya = self.start_y + yindexa * (self.end_y - self.start_y) / HEIGHT
                yb = self.start_y + yindexb * (self.end_y - self.start_y) / HEIGHT
                intervalx = interval(xa, xb)
                intervaly = interval(ya, yb)
                func_eval = func(intervalx, intervaly)
                #The expression is valid in the interval. Change the contour array
                #values to 1.
                if func_eval[1] is False or func_eval[0] is False:
                    pass
                elif func_eval == (True, True) and not is_equality:
                    contour[yindexa:yindexb, xindexa:xindexb] = 1
                elif func_eval[1] is None or func_eval[0] is None \
                        or is_equality:
                    #Subdivide
                    avgx_index = (xindexa + xindexb) // 2
                    avgy_index = (yindexa + yindexb) // 2
                    a = interval(xindexa, avgx_index)
                    b = interval(avgx_index, xindexb)
                    c = interval(yindexa, avgy_index)
                    d = interval(avgy_index, yindexb)
                    temp_interval_list.append([a, c])
                    temp_interval_list.append([a, d])
                    temp_interval_list.append([b, c])
                    temp_interval_list.append([b, d])
            return temp_interval_list
        while k >= 0 and len(interval_list):
            interval_list = refine_pixels(interval_list)
            k = k - 1
        #Check whether the expression represents an equality
        #If it represents an equality, then none of the intervals
        #would have satisfied the expression due to floating point
        #differences. Add all the undecided values to the plot.
        if isinstance(self.expr, (Equality, GreaterThan, LessThan)):
            for intervals in interval_list:
                xindexa, xindexb = intervals[0].start, intervals[0].end
                yindexa, yindexb = intervals[1].start, intervals[1].end
                xa = self.start_x + xindexa * (self.end_x - self.start_x) / WIDTH
                xb = self.start_x + xindexb * (self.end_x - self.start_x) / WIDTH
                ya = self.start_y + yindexa * (self.end_y - self.start_y) / HEIGHT
                yb = self.start_y + yindexb * (self.end_y - self.start_y) / HEIGHT
                intervalx = interval(xa, xb)
                intervaly = interval(ya, yb)
                func_eval = func(intervalx, intervaly)
                if func_eval[1] and func_eval[0] is not False:
                    contour[yindexa:yindexb, xindexa:xindexb] = 1
        xvals = np.linspace(self.start_x, self.end_x, WIDTH)
        yvals = np.linspace(self.start_y, self.end_y, WIDTH)

        return xvals, yvals, contour


def plot_implicit(expr, var_start_end_x, var_start_end_y, **kwargs):
    """A plot function to plot implicit equations / inequations.

    The input arguments are:
    expr : The equation / inequation that is to be plotted.
    var_start_end_x: A tuple of length 3, with the first element representing
    the variable and the next two elements representing the range
    var_start_end_y: A tuple of length 3, with the first element representing
    the variable and the next two elements representing the range

    Examples:
    ---------

    Plot expressions:
    >>> from sympy import plot_implicit, cos, sin, symbols
    >>> x, y = symbols('x y u v')
    >>> p1 = plot_implicit(Eq(y, x ** 2), (x, -5, 5), (y, -5, 5), show=False)
    >>> p2 = plot_implicit(Eq(x ** 2 + y ** 2, 3), (x, -3, 3), (y, -3, 3), show=False)
    >>> p3 = plot_implicit(y ** 2 < x ** 3 - x, (x, -4, 4), (y, -4, 4), show=False)
    >>> p4 = plot_implicit(y > sin(x), (x, -5, 5), (y, -2, 2), show=False)
"""

    assert isinstance(expr, Expr)
    free_symbols = expr.free_symbols
    range_symbols = set([var_start_end_x[0], var_start_end_y[0]])
    symbols = set_union(free_symbols, range_symbols)
    if len(symbols) > 2:
        raise NotImplementedError("Implicit plotting is not implemented for \
                                    more than 2 variables")
    series_argument = ImplicitSeries(expr, var_start_end_x, var_start_end_y)
    show = kwargs.pop('show', True)
    p = Plot(series_argument, **kwargs)
    if show:
        p.show()
    return p
