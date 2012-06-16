"""Plotting module for Sympy.

A plot is represented by the ``Plot`` class that contains a reference to the
backend and a list of the data series to be plotted. The data series are
instances of classes meant to simplify getting points and meshes from sympy
expressions. ``plot_backends`` is a dictionary with all the backends.

This module gives only the essential. For all the fancy stuff use directly
the backend. You can get the backend wrapper for every plot from the
``_backend`` attribute. Moreover the data series classes have various useful
methods like ``get_points``, ``get_segments``, ``get_meshes``, etc, that may
be useful if you wish to use another plotting library.

Especially if you need publication ready graphs and this module is not enough
for you - just get the ``_backend`` attribute and add whatever you want
directly to it. In the case of matplotlib (the common way to graph data in
python) just copy ``_backend.fig`` which is the figure and ``_backend.ax``
which is the axis and work on them as you would on any other matplotlib object.

Simplicity of code takes much greater importance than performance. Don't use it
if you care at all about performance. A new backend instance is initialized
every time you call ``show()`` and the old one is left to the garbage collector.
"""

from inspect import getargspec
from itertools import repeat, izip
from sympy import sympify, Expr, Tuple, Dummy
from sympy.external import import_module
from sympy.core.compatibility import set_union

from experimental_lambdify import vectorized_lambdify

#TODO probably all of the imports after this line can be put inside function to
# speed up the `from sympy import *` command.
np = import_module('numpy')

# Backend specific imports - matplotlib
matplotlib = import_module('matplotlib',
    __import__kwargs={'fromlist':['pyplot', 'cm', 'collections']},
    min_module_version='1.0.0')
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

    For interactive work the function ``plot`` is better suited.

    This class permits the plotting of sympy expressions using numerous
    backends (matplotlib, textplot, the old pyglet module for sympy, Google
    charts api, etc).

    The figure can contain an arbitrary number of plots of sympy expressions,
    lists of coordinates of points, etc. Plot has a private attribute _series that
    contains all data series to be plotted (expressions for lines or surfaces,
    lists of points, etc (all subclasses of BaseSeries)). Those data series are
    instances of classes not imported by ``from sympy import *``.

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
    Their implementation depends on the backend so they may not work in some
    backends.

    If the plot is parametric and the arity of the aesthetic function permits
    it the aesthetic is calculated over parameters and not over coordinates.
    If the arity does not permit calculation over parameters the calculation is
    done over coordinates.

    Only cartesian coordinates are suported for the moment, but you can use
    the paremetric plots to plot in polar, spherical and cylindrical
    coordinates.

    The arguments for the constructor Plot must be subclasses of BaseSeries.

    Any global option can be specified as a keyword argument.

    The global options for a figure are:
    title : str
    xlabel : str
    ylabel : str
    legend : bool
    xscale : {'linear', 'log'}
    yscale : {'linear', 'log'}
    axis : bool
    axis_center : tuple of two floats or {'center', 'auto'}
    xlim : tuple of two floats
    ylim : tuple of two floats
    aspect_ratio : tuple of two floats or {'auto'}
    autoscale : bool
    margin : float in [0, 1[

    The per data series options and aesthetics are:
    There are none in the base series. See below for options for subclasses.

    Some data series support additional aesthetics or options:
    is_line:
      ListSeries, LineOver1DRangeSeries,
      Parametric2DLineSeries, Parametric3DLineSeries
        aesthetics:
          line_color : float
        options:
          label : str
          steps : bool
          integers_only : bool
    is_3Dsurface:
      SurfaceOver2DRangeSeries, ParametricSurfaceSeries
        aesthetics:
          surface_color : float
    """

    def __init__(self, *args, **kwargs):
        super(Plot, self).__init__()

        #  Options for the graph as a whole.
        #  The possible values for each option are described in the docstring of
        # Plot. They are based purely on convention, no checking is done.
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.aspect_ratio = 'auto'
        self.xlim = None
        self.ylim = None
        self.axis_center = 'auto'
        self.axis = True
        self.xscale = 'linear'
        self.yscale = 'linear'
        self.legend = False
        self.autoscale = True
        self.margin = 0

        # Contains the data objects to be plotted. The backend should be smart
        # enough to iterate over this list.
        self._series = []
        self._series.extend(args)

        # The backend type. On every show() a new backend instance is created
        # in self._backend which is tightly coupled to the Plot instance
        # (thanks to the parent attribute of the backend).
        self.backend = DefaultBackend

        # The keyword arguments should only contain options for the plot.
        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)

    def show(self):
        # TODO move this to the backend (also for save)
        if hasattr(self, '_backend'):
            self._backend.close()
        self._backend = self.backend(self)
        self._backend.show()

    def save(self, path):
        if hasattr(self, '_backend'):
            self._backend.close()
        self._backend = self.backend(self)
        self._backend.save(path)

    def __str__(self):
        series_strs = [('[%d]: ' % i) + str(s)
                          for i, s in enumerate(self._series)]
        return 'Plot object containing:\n' + '\n'.join(series_strs)

    def __getitem__(self, index):
        return self._series[index]

    def __setitem__(self, index, *args):
        # XXX for ease of use here we permit also arguments for plot()
        if len(args)==1 and isinstance(args[0], BaseSeries):
            self._series[index] = args
        else:
            p = plot(*args, **{'show':False})
            self.extend(p)

    def __delitem__(self, index):
        del self._series[index]

    def append(self, *args):
        """Adds one more graph to the figure."""
        if len(args) == 1 and isinstance(args[0], BaseSeries):
            self._series.append(*args)
        else:
            self._series.append(Series(*args))

    def extend(self, arg):
        """Adds the series from another plot or a list of series."""
        if isinstance(arg, Plot):
            self._series.extend(arg._series)
        else:
            self._series.extend(arg)


def plot(*args, **kwargs):
    """A plot function for interactive use.

    It implements many heuristics in order to guess what the user wants on
    incomplete input.

    There is also the 'show' argument that defaults to True (immediately
    showing the plot).

    The input arguments can be:
      - lists with coordinates for ploting a line in 2D or 3D
      - the expressions and variable lists with ranges in order to plot any of
        the following: 2d line, 2d parametric line, 3d parametric line,
        surface, parametric surface
         - if the variable lists do not provide ranges a default range is used
         - if the variables are not provided, the free variables are
           automatically supplied in the order they are sorted (e.g. x, y, z)
         - if neither variables nor ranges are provided, both are guessed
         - if multiple expressions are provided in a list all of them are
           plotted
      - an instance of BaseSeries() subclass
      - another Plot() instance
      - tuples containing any of the above mentioned options, for plotting them
        together

    Examples:
    ---------

    Plot expressions:
    >>> from sympy import plot, cos, sin, symbols
    >>> x,y,u,v = symbols('x y u v')
    >>> p1 = plot(x**2, show=False) # with default [-10,+10] range
    >>> p2 = plot(x**2, (0, 5), show=False) # it finds the free variable itself
    >>> p3 = plot(x**2, (x, 0, 5), show=False) # fully explicit

    Fully implicit examples (finding the free variable and using default
    range). For the explicit versions just add the tuples with ranges:
    >>> p4 = plot(x**2, show=False) # cartesian line
    >>> p5 = plot(cos(u), sin(u), show=False) # parametric line
    >>> p6 = plot(cos(u), sin(u), u, show=False) # parametric line in 3d
    >>> p7 = plot(x**2 + y**2, show=False) # cartesian surface
    >>> p8 = plot(u, v, u+v, show=False) # parametric surface

    Multiple plots per figure:
    >>> p9 = plot((x**2, ), (cos(u), sin(u)), show=False) # cartesian and parametric lines

    Set title or other options:
    >>> p10 = plot(x**2, title='second order polynomial', show=False)

    Plot a list of expressions:
    >>> p11 = plot([x, x**2, x**3], show=False)
    >>> p12 = plot([x, x**2, x**3], (0,2), show=False) # explicit range
    >>> p13 = plot([x*y, -x*y], show=False) # list of surfaces

    And you can even plot a Plot or a Series object:
    >>> a = plot(x, show=False)
    >>> p14 = plot(a, show=False) # plotting a plot object
    >>> p15 = plot(a[0], show=False) # plotting a series object

    """

    plot_arguments = [p for p in args if isinstance(p, Plot)]
    series_arguments = [p for p in args if isinstance(p, BaseSeries)]
    args = sympify([np for np in args if not isinstance(np, (Plot, BaseSeries))])


    # Are the arguments for only one plot or are they tuples with arguments
    # for many plots. If it is the latter make the Tuples into list so they are
    # mutable.
    # args = (x, (x, 10, 20)) vs args = ((x, (x, 10, 20)), (y, (y, 20, 30)))
    if all([isinstance(a, Tuple) for a in args]):
        list_of_plots = map(list, args)
    else:
        list_of_plots = [list(args), ]

    # Check for arguments containing lists of expressions for the same ranges.
    # args = (x**2, (x, 10, 20)) vs args = ([x**3, x**2], (x, 10, 20))
    list_arguments = [p for p in list_of_plots
                              if any(isinstance(a, list) for a in p)]
    list_of_plots = [p for p in list_of_plots
                             if not any(isinstance(a, list) for a in p)]

    def expand(plot):
        """
        >> expand(([1,2,3],(1,2),(4,5)))
        [[1, (1, 2), (4, 5)], [2, (1, 2), (4, 5)], [3, (1, 2), (4, 5)]]
        """
        lists = [i for i in plot if isinstance(i, list)]
        not_lists = [i for i in plot if not isinstance(i, list)]
        return [list(l) + not_lists for l in zip(*lists)]

    for a in list_arguments:
        list_of_plots.extend(expand(a))

    def add_variables_and_ranges(plot):
        """make sure all limits are in the form (symbol, a, b)"""

        # find where the limits begin and expressions end
        for i in range(len(plot)):
            if isinstance(plot[i], Tuple):
                break
        else:
            i = len(plot) + 1
        exprs = list(plot[:i])
        assert all(isinstance(e, Expr) for e in exprs)
        assert all(isinstance(t, Tuple) for t in plot[i:])

        ranges = set([i for i in plot[i:] if isinstance(i, Tuple) and len(i) > 1])
        range_variables = set([t[0] for t in ranges if len(t) == 3])
        expr_free = set_union(*[e.free_symbols for e in exprs if isinstance(e, Expr)])

        default_range = Tuple(-10, 10)

        # unambiguous cases for limits
        #   no ranges
        if not ranges:
            plot = exprs + [Tuple(e) + default_range for e in expr_free or [Dummy()]]

        #   all limits of length 3
        elif all(len(i) == 3 for i in ranges):
            pass

        #   all ranges the same
        elif len(ranges) == 1:
            range1 = ranges.pop()
            if len(range1) == 2:
                plot = exprs + [Tuple(x) + range1 for x in expr_free]

        #   ranges cover free variables of expression
        elif expr_free < range_variables:
            plot = exprs + [i if len(i) == 3 else Tuple(Dummy()) + i for i in ranges]

        #   ranges cover all but 1 free variable
        elif len(expr_free - range_variables) == 1:
            x = (expr_free - range_variables).pop()
            ranges = list(ranges)
            for i, ri in enumerate(ranges):
                if len(ri) == 2:
                    ranges[i] = Tuple(x) + ri
                    break
            else:
                ranges.append(Tuple(x) + default_range)
            for i, ri in enumerate(ranges):
                if len(ri) == 2:
                    ranges[i] = Tuple(Dummy()) + ri
            plot = exprs + ranges

        #   all implicit ranges
        elif all(len(i) == 2 for i in ranges):
            more = len(ranges) - len(expr_free)
            all_free = list(expr_free) + [Dummy() for i in range(more)]
            ranges = list(ranges)
            ranges.extend(-more*[default_range])
            plot = exprs + [Tuple(all_free[i]) + ri for i, ri in enumerate(ranges)]

        else:
            raise ValueError('erroneous or unanticipated range input')

        return plot

    list_of_plots = [Tuple(*add_variables_and_ranges(pl)) for pl in list_of_plots]

    series_arguments.extend([OverloadedSeriesFactory(*pl) for pl in list_of_plots])

    show = kwargs.pop('show', True)
    p = Plot(*series_arguments, **kwargs)
    for plot_argument in plot_arguments:
        p.extend(plot_argument)
    if show:
        p.show()
    return p


##############################################################################
# Data Series
##############################################################################
#TODO more general way to cacollectionulate aesthetics (see get_color_array)

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
    # with the colors cacollectionulated at the points from get_points

    is_3Dline = False
    # Some of the backends expect:
    #  - get_points returning 1D np.arrays list_x, list_y, list_y
    #  - get_segments returning np.array (done in Line2DBaseSeries)
    #  - get_color_array returning 1D np.array (done in Line2DBaseSeries)
    # with the colors cacollectionulated at the points from get_points

    is_3Dsurface = False
    # Some of the backends expect:
    #   - get_meshes returning mesh_x, mesh_y, mesh_z (2D np.arrays)
    #   - get_points an alias for get_meshes

    is_contour = False
    # Some of the backends expect:
    #   - get_meshes returning mesh_x, mesh_y, mesh_z (2D np.arrays)
    #   - get_points an alias for get_meshes

    is_parametric = False
    # The cacollectionulation of aesthetics expects:
    #   - get_parameter_points returning one or two np.arrays (1D or 2D)
    # used for cacollectionulation aesthetics

    def __init__(self):
        super(BaseSeries, self).__init__()

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

    - adding the label, steps and only_integers options
    - making is_2Dline true
    - defining get_segments and get_color_array
    """

    is_2Dline = True

    _dim = 2

    def __init__(self):
        super(Line2DBaseSeries, self).__init__()
        self.label = None
        self.steps = False
        self.only_integers = False
        self.line_color = None

    def get_segments(self):
        points = self.get_points()
        if self.steps == True:
            x = np.array((points[0],points[0])).T.flatten()[1:]
            y = np.array((points[1],points[1])).T.flatten()[:-1]
            points = (x, y)
        points = np.array(points).T.reshape(-1, 1, self._dim)
        return np.concatenate([points[:-1], points[1:]], axis=1)

    def get_color_array(self):
        c = self.line_color
        if hasattr(c, '__call__'):
            f = np.vectorize(c)
            arity = len(getargspec(c)[0])
            if arity == 1 and self.is_parametric:
                x = self.get_parameter_points()
                return f(centers_of_segments(x))
            else:
                variables = map(centers_of_segments, self.get_points())
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
        self.nb_of_points = 300
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
        if self.only_integers == True:
            list_x = np.linspace(int(self.start), int(self.end),
                    num=int(self.end)-int(self.start)+1)
        else:
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
        self.nb_of_points = 300
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
        self.nb_of_points = 300
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
        self.surface_color = None

    def get_color_array(self):
        c = self.surface_color
        if callable(c):
            f = np.vectorize(c)
            arity = len(getargspec(c)[0])
            if self.is_parametric:
                variables = map(centers_of_faces,self.get_parameter_meshes())
                if arity == 1:
                    return f(variables[0])
                elif arity == 2:
                    return f(*variables)
            variables = map(centers_of_faces, self.get_meshes())
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

    is_parametric = True

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

    def get_parameter_meshes(self):
        return np.meshgrid(np.linspace(self.start_u, self.end_u,
                                       num=self.nb_of_points_u),
                           np.linspace(self.start_v, self.end_v,
                                       num=self.nb_of_points_v))

    def get_meshes(self):
        mesh_u, mesh_v = self.get_parameter_meshes()
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
class OverloadedSeriesFactory(BaseSeries):
    """ Construct a data series representation that makes sense for the given
    arguments. It works for a small subset of all possible plots.
    """
    # overloading would be great here :(
    # or the new function signature stuff from PEP 362
    def __new__(cls, *args):
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

    def process_series(self):
        parent = self.parent

        for s in self.parent._series:
            # Create the collections
            if s.is_2Dline:
                collection = LineCollection(s.get_segments())
                self.ax.add_collection(collection)
            elif s.is_contour:
                self.ax.contour(*s.get_meshes())
            elif s.is_3Dline:
                # TODO too complicated, I blame matplotlib
                collection = art3d.Line3DCollection(s.get_segments())
                self.ax.add_collection(collection)
                x, y, z = s.get_points()
                self.ax.set_xlim((min(x), max(x)))
                self.ax.set_ylim((min(y), max(y)))
                self.ax.set_zlim((min(z), max(z)))
            elif s.is_3Dsurface:
                x, y, z = s.get_meshes()
                collection = self.ax.plot_surface(x, y, z, cmap=cm.jet,
                                                  rstride=1, cstride=1,
                                                  linewidth = 0.1)
            else:
                raise ValueError('The matplotlib backend supports only '
                                 'is_2Dline, is_3Dline, is_3Dsurface and '
                                 'is_contour objects.')

            # Customise the collections with the corresponding per-series
            # options.
            if hasattr(s, 'label'):
                collection.set_label(s.label)
            if s.is_line and s.line_color:
                if isinstance(s.line_color, (float,int)) or callable(s.line_color):
                    color_array = s.get_color_array()
                    collection.set_array(color_array)
                else:
                    collection.set_color(s.line_color)
            if s.is_3Dsurface and s.surface_color:
                if matplotlib.__version__ < "1.2.0": #TODO in the distant future remove this check
                    warn.waring('The version of matplotlib is too old to use surface coloring.')
                elif isinstance(s.surface_color, (float,int)) or callable(s.surface_color):
                    color_array = s.get_color_array()
                    color_array = color_array.reshape(color_array.size)
                    collection.set_array(color_array)
                else:
                    collection.set_color(s.surface_color)

        # Set global options.
        # TODO The 3D stuff
        # XXX The order of those is important.
        if parent.xscale and not isinstance(self.ax, Axes3D):
            self.ax.set_xscale(parent.xscale)
        if parent.yscale and  not isinstance(self.ax, Axes3D):
            self.ax.set_yscale(parent.yscale)
        if parent.xlim:
            self.ax.set_xlim(parent.xlim)
        if parent.ylim:
            self.ax.set_ylim(parent.ylim)
        self.ax.set_autoscale_on(parent.autoscale)
        if parent.axis_center:
            val = parent.axis_center
            if isinstance(self.ax, Axes3D):
                pass
            elif val == 'center':
                self.ax.spines['left'].set_position('center')
                self.ax.spines['bottom'].set_position('center')
            elif val == 'auto':
                xl, xh = self.ax.get_xlim()
                yl, yh = self.ax.get_ylim()
                pos_left = ('data', 0) if xl*xh <= 0 else 'center'
                pos_bottom = ('data', 0) if yl*yh <= 0 else 'center'
                self.ax.spines['left'].set_position(pos_left)
                self.ax.spines['bottom'].set_position(pos_bottom)
            else:
                self.ax.spines['left'].set_position(('data', val[0]))
                self.ax.spines['bottom'].set_position(('data', val[1]))
        if not parent.axis:
            self.ax.set_axis_off()
        if parent.legend:
            self.ax.legend()
            self.ax.legend_.set_visible(parent.legend)
        if parent.margin:
            self.ax.set_xmargin(parent.margin)
            self.ax.set_ymargin(parent.margin)
        if parent.title:
            self.ax.set_title(parent.title)
        if parent.xlabel:
            self.ax.set_xlabel(parent.xlabel, position=(1,0))
        if parent.ylabel:
            self.ax.set_ylabel(parent.ylabel, position=(0,1))

    def show(self):
        self.process_series()
        #TODO after fixing https://github.com/ipython/ipython/issues/1255
        # you can uncomment the next line and remove the pyplot.show() call
        #self.fig.show()
        plt.show()

    def save(self, path):
        self.process_series()
        self.fig.savefig(path)

    def close(self):
        plt.close(self.fig)


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

    def close(self):
        pass


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
# Finding the centers of line segments or mesh faces
##############################################################################

def centers_of_segments(array):
    return np.average(np.vstack((array[:-1], array[1:])), 0)

def centers_of_faces(array):
    return np.average(np.dstack((array[:-1, :-1],
                                 array[1: , :-1],
                                 array[:-1, 1: ],
                                 array[:-1, :-1],
                                 )), 2)
