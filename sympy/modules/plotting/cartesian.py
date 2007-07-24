from pyglet.gl import *
from plot_function import PlotFunction, PlotFunctionRegistry, count_vars
from cartesian2d import CartesianFunction2d
from cartesian3d import CartesianFunction3d

class CartesianFunction(PlotFunction):

    def __new__(cls, f, intervals, options):
        d = max( [count_vars(f), len(intervals)] )
        if d == 1:
            return object.__new__(CartesianFunction2d, f, intervals, options)
        elif d == 2:
            return object.__new__(CartesianFunction3d, f, intervals, options)
        else:
            raise ValueError("Cannot plot a cartesian function with %i variables." % d)

PlotFunctionRegistry.register('cartesian', CartesianFunction)
