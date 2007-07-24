from pyglet.gl import *
from plot_function import PlotFunction, PlotFunctionRegistry, get_vars, count_vars, vrange
from polar2d import PolarFunction2d
from polar3d import PolarFunction3d

class PolarFunction(PlotFunction):

    def __new__(cls, f, intervals, options):
        d = max( [count_vars(f), len(intervals)] )
        if d == 1:
            return object.__new__(PolarFunction2d, f, intervals, options)
        elif d == 2:
            return object.__new__(PolarFunction3d, f, intervals, options)
        else:
            raise ValueError("Cannot plot a polar function with %i variables." % d)

PlotFunctionRegistry.register('polar', PolarFunction)
