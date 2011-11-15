"""Plotting module for Sympy.

A plot is represented by the Plot class that contains a reference to the
backend and a list of the data series to be plotted. The data series are
instances of classes meant to simplify getting points and meshes from sympy
expressions. `plot_backends` is a dictionary with all the backends.

This module gives only the essential. For all the fancy stuff use directly
the backend. You can get the backend wrapper for every plot from the _backend
attribute. Moreover the data series classes have various useful methods like
get_points, get_segments, get_meshes, etc, that may be useful if you wish to
use another plotting library.

Especially if you need publication ready graphs and this module is not enough
for you - just get the _backend attribute and add whatever you want directly to
it. In the case of matplotlib (the common way to graph data in python) just
copy _backend.fig which is the figure and _backend.ax which is the axis and
work on them as you would on any other matplotlib objects.

Simplicity of code takes much greater importance than performance. Don't use it
if you care at all about performance. That is especially true about aesthetics.
The color aesthetic for example creates an array with colors for each line
segment even if the the color is constant. process_series is called stupidly
often. A new backend instance is initialized every time you call show() and the
old one is left to the garbage collector.
"""

import warnings
from inspect import getargspec
from sympy import sympify, Expr, Tuple
from sympy.external import import_module

from experimental_lambdify import experimental_lambdify

np = import_module('numpy')

# Backend specific imports - matplotlib
matplotlib = import_module('matplotlib',
    __import__kwargs={'fromlist':['pyplot', 'cm', 'collections']})
if matplotlib:
    plt = matplotlib.pyplot
    cm = matplotlib.cm
    LineCollection = matplotlib.collections.LineCollection
    mpl_toolkits = import_module('mpl_toolkits',
            __import__kwargs={'fromlist':['mplot3d']})
    Axes3D = mpl_toolkits.mplot3d.Axes3D
    art3d = mpl_toolkits.mplot3d.art3d

# Backend specific imports - textplot
from sympy.plotting.textplot import textplot


##############################################################################
# The public interface
##############################################################################

class Plot(object):
    """The central class of the plotting module.

     For interactive work the function plot is better suited.

     This class permits the plotting of sympy expressions using numerous
    backends (matplotlib, textplot, the old pyglet module for sympy, Google
    charts api, etc).

     The figure can contain an arbitrary number of plots of sympy expressions,
    lists of coordinates of points, etc. Plot has a private attribute _series that
    contains all data series to be plotted (expressions for lines or surfaces,
    lists of points, etc (all subclasses of BaseSeries)). Those data series are
    instances of classes private to the module. They are listed below.

     The customization of the figure is on two levels. Global options that
    concern the figure as a whole (eg title, xlabel, scale, etc) and
    per-data series options (eg name) and aesthetics (eg. color, point shape,
    line type, etc.).

     The difference between options and aesthetics is that an aesthetic can be
    a function of the coordinates (or parameters in a parametric plot). The
    supported values for an aesthetic are:
    - None (the backend uses default values)
    - a constant
    - a function of one variable (the first coordinate or parameter)
    - a function of two variables (the first and second coordinate or
    parameters)
    - a function of three variables (only in nonparametric 3D plots)
     If the plot is parametric and the arity of the aesthetic function permits
    it the aesthetic is calculated over parameters and not over coordinates.
    If the arity does not permit calculation over parameters the calculation is
    done over coordinates.

     Only cartesian coordinates are suported for the moment, but you can use
    the paremetric plots to plot in polar, spherical and cylindrical
    coordinates.

    The arguments for the constructor Plot can be any of the following:
    - [x coords], [y coords] - 2D list plot
        implemented in ListSeries
    - expr, (var, start, end) - 2D function plot
        implemented in LineOver1DRangeSeries
    - expr_x, expr_y, (var, start, end) - 2D parametric plot
        implemented in Parametric2DLineSeries
    - expr_x, expr_y, expr_z, (var, start, end) - 3D parametric line plot
        implemented in Parametric3DLineSeries
    - expr, (var_x, start, end), (var_y, start, end) - surface plot
        implemented in SurfaceOver2DRangeSeries
    - expr_x, expr_y, expr_z, (var_u, start, end), (var_v, start, end) - 3D
      parametric surface plot
        implemented in ParametricSurfaceSeries
    - subclasses of BaseSeries
    - tuples containing any of the aforementioned options

    Any global option can be specified as a keyword argument.

    The global options for a figure are:
    (None leaves the default from the backend)
    title : str
    xlabel : str
    ylabel : str
    legend : bool
    xscale : {'linear', 'log'}
    yscale : {'linear', 'log'}
    axis : bool
    axis_center : tuple of two floats or {'center'}
    xlim : tuple of two floats
    ylim : tuple of two floats
    aspect_ratio : tuple of two floats or {'auto'}
    autoscale : bool
    margin : float in [0, 1[

    The per data series options are:
    None in the base series. See below for options for subclasses.

    The aesthetics for the data series are:
    None in the base series. See below for options for subclasses.

    Some data series support additional aesthetics or options:
    is_line:
      ListSeries, LineOver1DRangeSeries,
      Parametric2DLineSeries, Parametric3DLineSeries
        aesthetics:
          line_color : float
        options:
          label : str
    is_3Dsurface:
      SurfaceOver2DRangeSeries, ParametricSurfaceSeries
        aesthetics:
          surface_color : float
        options:
          none
    is_contour:
      ContourSeries
        aesthetics:
          none
        options:
          none

    Examples
    --------
    >>> # Creating Plot
    >>> # Creating Plot of multiple data series
    >>> # Setting global option (title)
    >>> # Setting per-data series option (label)
    >>> # Creating legend (global option)
    >>> # Changing aesthetics (color set to a constant)
    >>> # Changing aesthetics (color proportional to length)
    >>> # Changing backend
    >>> # Appending and extending plots
    """

    def __init__(self, *args, **kwargs):
        super(Plot, self).__init__()

        #  Options for the graph as a whole. The __setattr__ method is overridden so
        # every key from this dictionary is exposed as an attribute.
        #  The possible values for each option are described in the docstring of
        # Plot. They are based purely on convention, no checking is done.
        #  The rationale for using a dictionary for the options and aesthetics
        # here and in the Series classes is that it is easier to list them
        # later or to loop over them (compared to class attributes).
        _options = {
            'title': None,
            'xlabel': None,
            'ylabel': None,
            'aspect_ratio': 'auto',
            'xlim': None,
            'ylim': None,
            'axis_center': (0,0),
            'axis': True,
            'xscale': 'linear',
            'yscale': 'linear',
            'legend': False,
            'autoscale': True,
            'margin': 0,
        }
        # The __setattr__ of Plot is overridden and it needs _options to be
        # already defined. So we use the parent class __setattr__.
        super(Plot, self).__setattr__('_options', _options)

        # Contains the data objects to be plotted. The backend should be smart
        # enough to iterate over this list and do the appropriate actions for
        # every subclass of BaseSeries. The iteration is done in process_series
        # in BaseBackend and the appropriate action is done in
        # _series_representation in the backend instance.
        self._series = []

        # The arguments are either:
        #  - arguments for the Series subclass or instances
        # of subclasses of BaseSeries
        #  - list of tuples containing the above
        args = sympify(args)
        if len(args) != 0:
            if not all([isinstance(a, Tuple) for a in args]):
                self._series.append(Series(*args))
            else:
                for a in args:
                    self._series.append(Series(*a))

        # The backend type. On every show() a new backend instance is created
        # in self._backend which is tightly coupled to the Plot instance
        # (thanks to the parent attribute of the backend).
        self.backend = DefaultBackend

        # The keyword arguments should only contain options for the plot.
        for key, val in kwargs.iteritems():
            if key in self._options:
                setattr(self, key, val)

    def __setattr__(self, name, val):
        if name in self._options:
            self._options[name] = val
        else:
            super(Plot, self).__setattr__(name, val)

    def __getattr__(self, name):
        if name in self._options:
            return self._options[name]
        else:
            raise AttributeError('The backend has no such attribute ' + name)

    def show(self):
        if hasattr(self, '_backend'):
            self._backend.close()
        self._backend = self.backend(self)
        self._backend.show()

    def __str__(self):
        series_strs = [('[%d]: ' % i) + str(s)
                          for i, s in enumerate(self._series)]
        return 'Plot object containing:\n' + '\n'.join(series_strs)

    def __getitem__(self, index):
        return self._series[index]

    def __setitem__(self, index, *args):
        self._series[index] = Series(*args)

    def __delitem__(self, index):
        del self._series[index]

    def append(self, *args):
        """Adds one more graph to the figure."""
        if len(args) == 1 and isinstance(args[0], BaseSeries):
            self._series.append(*args)
        else:
            self._series.append(Series(*args))

    def extend(self, plot):
        """Adds the series from another plot."""
        self._series.extend(plot._series)


def plot(*args, **kwargs):
    """Convenient interface for the *new* Plot class.

    This function has a more relaxed input requirements than the Plot class
    and shows the figure immediately. As such is better suited for interactive
    work.

    It returns a Plot object with a matplotlib backend.
    For the more advanced options, detailed documentation and other types of
    plots see the Plot class from sympy.plotting.plot.

    The new Plot class is not as of yet imported by 'from sympy import *'.
    For compability with older versions of sympy 'from sympy import *' still
    imports the old Plot class which now resides in sympy.plotting.pygletplot.

    All keyword arguments that Plot supports are also supported by plot (like
    title, etc).

    There is also the 'show' argument that defaults to True (immediately
    showing the plot).

    Examples:
    ---------

    Plot lists:
    >> listx = range(10)
    >> listy = [x**2 for x in listx]
    >> p0 = plot(listx, listy)

    Plot expressions:
    >> p1 = plot(x**2) # with default [-10,+10] range
    >> p2 = plot(x**2, (0, 5)) # it finds the free variable itself
    >> p3 = plot(x**2, (x, 0, 5)) # fully explicit

    Fully implicit examples (finding the free variable and using default
    range). For the explicit versions just add the tuples with ranges:
    >> p4 = plot(x**2) # cartesian line
    >> p5 = plot(cos(u), sin(u)) # parametric line
    >> p6 = plot(cos(u), sin(u), u) # parametric line in 3d
    >> p7 = plot(x**2 + y**2) # cartesian surface
    >> p8 = plot(u, v, u+v) # parametric surface

    Multiple plots per figure:
    >> p9 = plot((x**2), (cos(u), sin(u))) # cartesian and parametric lines

    Set title or other options:
    >> p10 = plot(x**2, title='second order polynomial')
    """
    # The idea of this function is to permint more ambiguous input than Plot
    # and Series. It uses default values to make the input understandable for
    # Plot. It also shows the plots immediately.

    # We check if the user wants to plot one graph or multiple graps and then
    # add the necessary ranges and free variables if they are not explicitly
    # stated. All tuples are converted to lists for ease of work.

    # The code assumes that the only use for tuples is to give the range of the
    # plot either as (x, -4, +2) or just (-4, +2).

    args = sympify(args)
    default_range = Tuple(-10, 10)

    if all([isinstance(a, Tuple) for a in args]):
        list_of_plots = map(list, args)
    else:
        list_of_plots = [list(args), ]

    for pl_index, pl in enumerate(list_of_plots):
        if any([isinstance(e, Expr) for e in pl]):
            free_vars = set.union(*[expr.free_symbols for expr in pl
                                                if isinstance(expr, Expr)])
            while len([t for t in pl if isinstance(t, Tuple)]) < len(free_vars):
                pl.append(default_range)
            for index, item in enumerate(pl):
                if isinstance(item, Expr):
                    pass
                elif len(item) == 3:
                    free_vars.discard(item[0])
                elif len(item) == 2:
                    pl[index] = Tuple(free_vars.pop(), item[0], item[1])
        list_of_plots[pl_index] = Tuple(*pl)

    p = Plot(*list_of_plots, **kwargs)
    if 'show' in kwargs and kwargs['show'] == False:
        pass
    else:
        p.show()
    return p


##############################################################################
# Data Series
##############################################################################
#TODO nb_of_points can (should) became an option
#TODO more general way to calculate aesthetics (see get_color_array)

### The base class for all series
class BaseSeries(object):
    """Base class for the data objects containing stuff to be plotted.

     The backend should check if it supports the data series that it's given.
    (eg TextBackend supports only LineOver1DRange).
    It's the backend responsibility to know how to use the class of
    data series that it's given.

     Some data series classes are grouped (using a class attribute like is_2Dline)
    according to the api they present (based only on convention). The backend is
    not obliged to use that api (eg. The LineOver1DRange belongs to the
    is_2Dline group and presents the get_points method, but the
    TextBackend does not use the get_points method).
    """

    # Some flags follow. The rationale for using flags instead of checking base
    # classes is that setting multiple flags is simpler than multiple
    # inheritance.

    is_2Dline = False
    # Some of the backends expect:
    #  - get_points returning 1D np.arrays list_x, list_y
    #  - get_segments returning np.array (done in Line2DBaseSeries)
    #  - get_color_array returning 1D np.array (done in Line2DBaseSeries)
    # with the colors calculated at the points from get_points

    is_3Dline = False
    # Some of the backends expect:
    #  - get_points returning 1D np.arrays list_x, list_y, list_y
    #  - get_segments returning np.array (done in Line2DBaseSeries)
    #  - get_color_array returning 1D np.array (done in Line2DBaseSeries)
    # with the colors calculated at the points from get_points

    is_3Dsurface = False
    # Some of the backends expect:
    #   - get_meshes returning mesh_x, mesh_y, mesh_z (2D np.arrays)
    #   - get_points an alias for get_meshes

    is_contour = False
    # Some of the backends expect:
    #   - get_meshes returning mesh_x, mesh_y, mesh_z (2D np.arrays)
    #   - get_points an alias for get_meshes

    is_parametric = False
    # The calculation of aesthetics expects:
    #   - get_parameter_points returning one or two np.arrays (1D or 2D)
    # used for calculation aesthetics

    def __init__(self):
        super(BaseSeries, self).__init__()

        _aesthetics = {
        }

        _options = {
        }

        super(BaseSeries, self).__setattr__('_aesthetics', _aesthetics)
        super(BaseSeries, self).__setattr__('_options', _options)

    def __setattr__(self, name, val):
        if name in self._aesthetics:
            self._aesthetics[name] = val
        elif name in self._options:
            self._options[name] = val
        else:
            object.__setattr__(self, name, val)

    def __getattr__(self, name):
        if name in self._options:
            return self._options[name]
        elif name in self._aesthetics:
            return self._aesthetics[name]
        else:
            raise AttributeError('The backend has no such attribute ' + name)

    def _add_aesthetics(self, additional_aes):
        self._aesthetics.update(additional_aes)

    def _add_options(self, additional_ops):
        self._options.update(additional_ops)

    @property
    def is_3D(self):
        flags3D = [
            self.is_3Dline,
            self.is_3Dsurface
        ]
        return any(flags3D)

    @property
    def is_line(self):
        flagslines = [
            self.is_2Dline,
            self.is_3Dline
        ]
        return any(flagslines)


### 2D lines
class Line2DBaseSeries(BaseSeries):
    """A base class for 2D lines.

    - adding the label option
    - making is_2Dline true
    - defining get_segments and get_color_array
    """

    is_2Dline = True

    _dim = 2

    def __init__(self):
        super(Line2DBaseSeries, self).__init__()
        self._add_options({'label' : None})
        self._add_aesthetics({'line_color' : None})

    def get_segments(self):
        points = self.get_points()
        points = np.array(points).T.reshape(-1, 1, self._dim)
        return np.concatenate([points[:-1], points[1:]], axis=1)

    def get_color_array(self):
        c = self._aesthetics['line_color']
        if hasattr(c, '__call__'):
            f = np.vectorize(c)
            arity = len(getargspec(c)[0])
            if arity == 1 and self.is_parametric:
                x = self.get_parameter_points()
                return f(x)
            else:
                variables = self.get_points()
                if arity == 1:
                    return f(variables[0])
                elif arity == 2:
                    return f(*variables[:2])
                else: # only if the line is 3D (otherwise raises an error)
                    return f(*variables)
        else:
            return c*np.ones(self.nb_of_points)


class List2DSeries(Line2DBaseSeries):
    """Representation for a line consisting of list of points."""

    def __init__(self, list_x, list_y):
        super(List2DSeries, self).__init__()
        self.list_x = np.array(list_x)
        self.list_y = np.array(list_y)
        self.label = 'list'

    def __str__(self):
        return 'list plot'

    def get_points(self):
        return (self.list_x, self.list_y)


class LineOver1DRangeSeries(Line2DBaseSeries):
    """Representation for a line consisting of a sympy expression over a range."""

    def __init__(self, expr, var_start_end):
        super(LineOver1DRangeSeries, self).__init__()
        self.nb_of_points = 100
        self.expr = sympify(expr)
        self.label = str(self.expr)
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])

    def __str__(self):
        return 'cartesian line: %s for %s over %s' % (
                str(self.expr),
                str(self.var),
                str((self.start, self.end)))

    def get_points(self):
        list_x = np.linspace(self.start, self.end, num=self.nb_of_points)
        f = vectorized_lambdify([self.var], self.expr)
        list_y = f(list_x)
        return (list_x, list_y)


class Parametric2DLineSeries(Line2DBaseSeries):
    """Representation for a line consisting of two parametric sympy expressions
    over a range."""

    is_parametric = True

    def __init__(self, expr_x, expr_y, var_start_end):
        super(Parametric2DLineSeries, self).__init__()
        self.nb_of_points = 100
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.label = "(%s, %s)" % (str(self.expr_x), str(self.expr_y))
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])

    def __str__(self):
        return 'parametric cartesian line: (%s, %s) for %s over %s' % (
                str(self.expr_x),
                str(self.expr_y),
                str(self.var),
                str((self.start, self.end)))

    def get_parameter_points(self):
        return np.linspace(self.start, self.end, num=self.nb_of_points)

    def get_points(self):
        param = self.get_parameter_points()
        fx = vectorized_lambdify([self.var], self.expr_x)
        fy = vectorized_lambdify([self.var], self.expr_y)
        list_x = fx(param)
        list_y = fy(param)
        return (list_x, list_y)


### 3D lines
class Line3DBaseSeries(Line2DBaseSeries):
    """A base class for 3D lines.

    Most of the stuff is derived from Line2DBaseSeries."""

    is_2Dline = False
    is_3Dline = True
    _dim = 3

    def __init__(self):
        super(Line3DBaseSeries, self).__init__()


class Parametric3DLineSeries(Line3DBaseSeries):
    """Representation for a 3D line consisting of two parametric sympy
    expressions and a range."""

    def __init__(self, expr_x, expr_y, expr_z, var_start_end):
        super(Parametric3DLineSeries, self).__init__()
        self.nb_of_points = 100
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.expr_z = sympify(expr_z)
        self.label = "(%s, %s)" % (str(self.expr_x), str(self.expr_y))
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])

    def __str__(self):
        return '3D parametric cartesian line: (%s, %s, %s) for %s over %s' % (
                str(self.expr_x),
                str(self.expr_y),
                str(self.expr_z),
                str(self.var),
                str((self.start, self.end)))

    def get_parameter_points(self):
        return np.linspace(self.start, self.end, num=self.nb_of_points)

    def get_points(self):
        param = self.get_parameter_points()
        fx = vectorized_lambdify([self.var], self.expr_x)
        fy = vectorized_lambdify([self.var], self.expr_y)
        fz = vectorized_lambdify([self.var], self.expr_z)
        list_x = fx(param)
        list_y = fy(param)
        list_z = fz(param)
        return (list_x, list_y, list_z)


### Surfaces
class SurfaceBaseSeries(BaseSeries):
    """A base class for 3D surfaces."""

    is_3Dsurface = True

    def __init__(self):
        super(SurfaceBaseSeries, self).__init__()
        self._add_aesthetics({'surface_color' : None})

        self.get_points = self.get_meshes

    def get_color_array(self):
        c = self._aesthetics['surface_color']
        if hasattr(c, '__call__'):
            f = np.vectorize(c)
            arity = len(getargspec(c)[0])
            if self.is_parametric:
                variables = self.get_parameter_points()
                if arity == 1:
                    return f(variables[0])
                else:
                    return f(*variables)
            else:
                variables = self.get_points()
                if arity == 1:
                    return f(variables[0])
                elif arity == 2:
                    return f(*variables[:2])
                else:
                    return f(*variables)
        else:
            return c*np.ones(self.nb_of_points)


class SurfaceOver2DRangeSeries(SurfaceBaseSeries):
    """Representation for a 3D surface consisting of a sympy expression and 2D
    range."""
    def __init__(self, expr, var_start_end_x, var_start_end_y):
        super(SurfaceOver2DRangeSeries, self).__init__()
        self.nb_of_points_x = 50
        self.nb_of_points_y = 50
        self.expr = sympify(expr)
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])

    def __str__(self):
        return ('cartesian surface: %s for'
                ' %s over %s and %s over %s') % (
                str(self.expr),
                str(self.var_x),
                str((self.start_x, self.end_x)),
                str(self.var_y),
                str((self.start_y, self.end_y)))

    def get_meshes(self):
        mesh_x, mesh_y = np.meshgrid(np.linspace(self.start_x, self.end_x,
                                                 num=self.nb_of_points_x),
                                     np.linspace(self.start_y, self.end_y,
                                                 num=self.nb_of_points_y))
        f = vectorized_lambdify((self.var_x, self.var_y), self.expr)
        return (mesh_x, mesh_y, f(mesh_x, mesh_y))


class ParametricSurfaceSeries(SurfaceBaseSeries):
    """Representation for a 3D surface consisting of three parametric sympy
    expressions and a range."""
    def __init__(self, expr_x, expr_y, expr_z, var_start_end_u, var_start_end_v):
        super(ParametricSurfaceSeries, self).__init__()
        self.nb_of_points_u = 50
        self.nb_of_points_v = 50
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.expr_z = sympify(expr_z)
        self.var_u = sympify(var_start_end_u[0])
        self.start_u = float(var_start_end_u[1])
        self.end_u = float(var_start_end_u[2])
        self.var_v = sympify(var_start_end_v[0])
        self.start_v = float(var_start_end_v[1])
        self.end_v = float(var_start_end_v[2])

    def __str__(self):
        return ('parametric cartesian surface: (%s, %s, %s) for'
                ' %s over %s and %s over %s') % (
                str(self.expr_x),
                str(self.expr_y),
                str(self.expr_z),
                str(self.var_u),
                str((self.start_u, self.end_u)),
                str(self.var_v),
                str((self.start_v, self.end_v)))

    def get_meshes(self):
        mesh_u, mesh_v = np.meshgrid(np.linspace(self.start_u, self.end_u,
                                                 num=self.nb_of_points_u),
                                     np.linspace(self.start_v, self.end_v,
                                                 num=self.nb_of_points_v))
        fx = vectorized_lambdify((self.var_u, self.var_v), self.expr_x)
        fy = vectorized_lambdify((self.var_u, self.var_v), self.expr_y)
        fz = vectorized_lambdify((self.var_u, self.var_v), self.expr_z)
        return (fx(mesh_u, mesh_v), fy(mesh_u, mesh_v), fz(mesh_u, mesh_v))


### Contours
class ContourSeries(BaseSeries):
    """Representation for a contour plot."""
    #The code is mostly repetition of SurfaceOver2DRange.

    is_contour = True

    def __init__(self, expr, var_start_end_x, var_start_end_y):
        super(ContourSeries, self).__init__()
        self.nb_of_points_x = 50
        self.nb_of_points_y = 50
        self.expr = sympify(expr)
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])

        self.get_points = self.get_meshes

    def __str__(self):
        return ('contour: %s for '
                '%s over %s and %s over %s') % (
                str(self.expr),
                str(self.var_x),
                str((self.start_x, self.end_x)),
                str(self.var_y),
                str((self.start_y, self.end_y)))

    def get_meshes(self):
        mesh_x, mesh_y = np.meshgrid(np.linspace(self.start_x, self.end_x,
                                                 num=self.nb_of_points_x),
                                     np.linspace(self.start_y, self.end_y,
                                                 num=self.nb_of_points_y))
        f = vectorized_lambdify((self.var_x, self.var_y), self.expr)
        return (mesh_x, mesh_y, f(mesh_x, mesh_y))


### Factory class
class Series(BaseSeries):
    """ Construct a data series representation that makes sense for the given
    arguments. It does not work for contour plot for the moment.
    """
    # overloading would be great here :(
    # or the new function signature stuff from PEP 362
    def __new__(cls, *args):
        if isinstance(args[0], BaseSeries):
            return args[0]
        if (len(args) == 2
            and isinstance(args[0], list)
            and isinstance(args[1], list)):
            return List2DSeries(*args)

        args = sympify(args)
        if (len(args) == 2
              and isinstance(args[0], Expr)
              and isinstance(args[1], Tuple)
              and len(args[1]) == 3):
            inst = LineOver1DRangeSeries(*args)
        elif (len(args) == 3
              and isinstance(args[0], Expr)
              and isinstance(args[1], Expr)
              and isinstance(args[2], Tuple)
              and len(args[2]) == 3):
            inst = Parametric2DLineSeries(*args)
        elif (len(args) == 3
              and isinstance(args[0], Expr)
              and isinstance(args[1], Tuple)
              and isinstance(args[2], Tuple)
              and len(args[1]) == 3
              and len(args[2]) == 3):
            inst = SurfaceOver2DRangeSeries(*args)
        elif (len(args) == 4
              and isinstance(args[0], Expr)
              and isinstance(args[1], Expr)
              and isinstance(args[2], Expr)
              and isinstance(args[3], Tuple)
              and len(args[3]) == 3):
            inst = Parametric3DLineSeries(*args)
        elif (len(args) == 5
              and isinstance(args[0], Expr)
              and isinstance(args[1], Expr)
              and isinstance(args[2], Expr)
              and isinstance(args[3], Tuple)
              and isinstance(args[4], Tuple)
              and len(args[3]) == 3
              and len(args[4]) == 3):
            inst = ParametricSurfaceSeries(*args)
        else:
            raise ValueError('The supplied argument do not correspond to a'
                             ' valid plotable object.')
        return inst


##############################################################################
# Backends
##############################################################################

class BaseBackend(object):
    def __init__(self, parent):
        super(BaseBackend, self).__init__()
        self.parent = parent
        self._dict_of_series = {} # Each backend may need to keep stuff about
        # each data series. This dict is meant for this. It's populated using
        # _series_representation during process_series. If a backend is not
        # supporting certain data series _series_representation is the best
        # place to raise an error. Another possible place is the show function.
        # The class may create helper classes if it needs complicated
        # representations or it can just return None if it does not need any.

    def process_options(self):
        for opt, val in self.parent._options.iteritems():
            if val != None:
                getattr(self, 'set_glo_'+opt)(val)

    def process_series(self):
        for s in self.parent._series:
            if not s in self._dict_of_series:
                self._dict_of_series[s] = self._series_representation(s)
            for opt, val in s._options.iteritems():
                if val != None:
                    getattr(self, 'set_ser_opt_'+opt)(s, val)
            for aes, val in s._aesthetics.iteritems():
                if val != None:
                    getattr(self, 'set_ser_aes_'+aes)(s, val)

    def __getattr__(self, name):
        if name.startswith('set_glo_'):
            warnings.warn('The global option setter ' + name +
                          ' is not implemented in the backend. ' +
                          'The options is not available.')
            return self._dummy
        elif name.startswith('set_ser_opt_'):
            warnings.warn('The series option setter ' + name +
                          ' is not implemented in the backend. ' +
                          'The options is not available.')
            return self._dummy
        elif name.startswith('set_ser_aes_'):
            warnings.warn('The series aesthetic setter ' + name +
                          ' is not implemented in the backend. ' +
                          'The aesthetic is not available.')
            return self._dummy
        else:
            raise AttributeError('The backend has no such attribute ' + name)

    def _dummy(self, *args):
        return None

    def show(self):
        pass

    def close(self):
        pass

    def _series_representation(self):
        return None


class MatplotlibBackend(BaseBackend):
    def __init__(self, parent):
        super(MatplotlibBackend, self).__init__(parent)
        are_3D = [s.is_3D for s in self.parent._series]
        if any(are_3D) and not all(are_3D):
            raise ValueError('The matplotlib backend can not mix 2D and 3D.')
        elif not any(are_3D):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.spines['left'].set_position('zero')
            self.ax.spines['right'].set_color('none')
            self.ax.spines['bottom'].set_position('zero')
            self.ax.spines['top'].set_color('none')
            self.ax.spines['left'].set_smart_bounds(True)
            self.ax.spines['bottom'].set_smart_bounds(True)
            self.ax.xaxis.set_ticks_position('bottom')
            self.ax.yaxis.set_ticks_position('left')
        elif all(are_3D):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            #blah blah

    def _series_representation(self, s):
        if s.is_2Dline:
            lc = LineCollection(s.get_segments())
            self.ax.add_collection(lc)
            return lc
        elif s.is_contour:
            return self.ax.contour(*s.get_meshes())
        elif s.is_3Dline:
            # TODO too complicated, I blame matplotlib
            lc = art3d.Line3DCollection(s.get_segments())
            self.ax.add_collection(lc)
            x, y, z = s.get_points()
            self.ax.set_xlim((min(x), max(x)))
            self.ax.set_ylim((min(y), max(y)))
            self.ax.set_zlim((min(z), max(z)))
            return lc
        elif s.is_3Dsurface:
            x, y, z = s.get_meshes()
            return self.ax.plot_surface(x, y, z, cmap=cm.jet,
                                        rstride=1, cstride=1,
                                        linewidth = 0.1)
        else:
            raise ValueError('The matplotlib backend supports only '
                             'is_2Dline, is_3Dline and is_contour objects.')

    def show(self):
        self.process_series()
        self.process_options()
        #TODO after fixing https://github.com/ipython/ipython/issues/1255
        # you can uncomment the next line and remove the pyplot.show() call
        #self.fig.show()
        plt.show()

    def close(self):
        plt.close(self.fig)

    def set_glo_title(self, val):
        self.ax.set_title(val)

    def set_glo_xlabel(self, val):
        self.ax.set_xlabel(val, position=(1,0))

    def set_glo_ylabel(self, val):
        self.ax.set_ylabel(val, position=(0,1))

    def set_glo_aspect_ratio(self, val):
        if isinstance(val, tuple):
            self.ax.set_aspect(float(val[1])/val[0])
        elif val == 'auto':
            self.ax.set_aspect('auto')
        else:
            raise ValueError('Bad value for the aspect arg. Use tuple or \'auto\'.')

    def set_glo_xlim(self, val):
        if val:
            self.ax.set_xlim(val)

    def set_glo_ylim(self, val):
        if val:
            self.ax.set_ylim(val)

    def set_glo_axis_center(self, val):
        if isinstance(self.ax, Axes3D):
            warnings.warn('axis_center is not supported in 3D matplotlib backend.')
        elif val == 'center':
            self.ax.spines['left'].set_position('center')
            self.ax.spines['bottom'].set_position('center')
        else:
            self.ax.spines['left'].set_position(('data', val[0]))
            self.ax.spines['bottom'].set_position(('data', val[1]))

    def set_glo_axis(self, val):
        if val:
            self.ax.set_axis_on()
        else:
            self.ax.set_axis_off()

    def set_glo_xscale(self, val):
        if isinstance(self.ax, Axes3D):
            warnings.warn('xscale is not supported in 3D matplotlib backend.')
        else:
            self.ax.set_xscale(val)
            #XXX In matplotlib xscale resets xlim, so we must set xlim again.
            self.set_glo_xlim(self.parent.xlim)

    def set_glo_yscale(self, val):
        if isinstance(self.ax, Axes3D):
            warnings.warn('yscale is not supported in 3D matplotlib backend.')
        else:
           self.ax.set_yscale(val)
           #XXX In matplotlib yscale resets ylim, so we must set ylim again.
           self.set_glo_ylim(self.parent.ylim)

    def set_glo_legend(self, val):
        if val:
            self.ax.legend()
            self.ax.legend_.set_visible(val)
        elif hasattr(self.ax, 'ledend_'):
            self.ax.legend_.set_visible(val)

    def set_glo_autoscale(self, val):
        self.ax.set_autoscale_on(val)

    def set_glo_margin(self, val):
        self.ax.set_xmargin(val)
        self.ax.set_ymargin(val)

    def set_ser_opt_label(self, series, val):
        self._dict_of_series[series].set_label(val)

    def set_ser_aes_line_color(self, series, val):
        if isinstance(val, (float,int)) or callable(val):
            color_array = series.get_color_array()
            self._dict_of_series[series].set_array(color_array[:-1])
        else:
            self._dict_of_series[series].set_color(val)

    def set_ser_aes_surface_color(self, series, val):
        if isinstance(val, (float,int)) or callable(val):
            color_array = np.array(series.get_color_array()[:-1,:-1])
            color_array.shape = color_array.shape[0] * color_array.shape[1]
            self._dict_of_series[series].set_array(color_array)
        else:
            self._dict_of_series[series].set_color(val)



class TextBackend(BaseBackend):
    def __init__(self, parent):
        super(TextBackend, self).__init__(parent)

    def show(self):
        if len(self.parent._series) != 1:
            raise ValueError('The TextBackend supports only one graph per Plot.')
        elif not isinstance(self.parent._series[0], LineOver1DRangeSeries):
            raise ValueError('The TextBackend supports only expressions over a 1D range')
        else:
            ser = self.parent._series[0]
            textplot(ser.expr, ser.start, ser.end)


class DefaultBackend(BaseBackend):
    def __new__(cls, parent):
        if matplotlib:
            return MatplotlibBackend(parent)
        else:
            return TextBackend(parent)


plot_backends = {
        'matplotlib' : MatplotlibBackend,
        'text' : TextBackend,
        }

##############################################################################
# A helper for making vectorized functions
##############################################################################
# There are so many levels of evalf, ufunc and vectorize
# here that this is most probably hundred times slower
# than The Right Solution TM.

def vectorized_lambdify(args, expr):
    return np.vectorize(experimental_lambdify(args, expr))
