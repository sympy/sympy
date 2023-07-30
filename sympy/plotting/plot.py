"""Plotting module for SymPy.

A plot is represented by the ``Plot`` class that contains a reference to the
backend and a list of the data series to be plotted. The data series are
instances of classes meant to simplify getting points and meshes from SymPy
expressions. ``plot_backends`` is a dictionary with all the backends.

This module gives only the essential. For all the fancy stuff use directly
the backend. You can get the backend wrapper for every plot from the
``_backend`` attribute. Moreover the data series classes have various useful
methods like ``get_points``, ``get_meshes``, etc, that may
be useful if you wish to use another plotting library.

Especially if you need publication ready graphs and this module is not enough
for you - just get the ``_backend`` attribute and add whatever you want
directly to it. In the case of matplotlib (the common way to graph data in
python) just copy ``_backend.fig`` which is the figure and ``_backend.ax``
which is the axis and work on them as you would on any other matplotlib object.

Simplicity of code takes much greater importance than performance. Do not use it
if you care at all about performance. A new backend instance is initialized
every time you call ``show()`` and the old one is left to the garbage collector.
"""


from collections.abc import Callable


from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.expr import Expr
from sympy.core.function import arity, Function
from sympy.core.symbol import (Dummy, Symbol)
from sympy.core.sympify import sympify
from sympy.external import import_module
from sympy.printing.latex import latex
from sympy.utilities.exceptions import sympy_deprecation_warning
from sympy.utilities.iterables import is_sequence
from .experimental_lambdify import (vectorized_lambdify, lambdify)

# N.B.
# When changing the minimum module version for matplotlib, please change
# the same in the `SymPyDocTestFinder`` in `sympy/testing/runtests.py`

# Backend specific imports - textplot
from sympy.plotting.textplot import textplot

# Global variable
# Set to False when running tests / doctests so that the plots don't show.
_show = True


def unset_show():
    """
    Disable show(). For use in the tests.
    """
    global _show
    _show = False

def _str_or_latex(label):
    if isinstance(label, Basic):
        return latex(label, mode='inline')
    return str(label)

def _create_generic_data_series(**kwargs):
    keywords = ["annotations", "markers", "fill", "rectangles"]
    series = []
    for kw in keywords:
        dictionaries = kwargs.pop(kw, [])
        if dictionaries is None:
            dictionaries = []
        if isinstance(dictionaries, dict):
            dictionaries = [dictionaries]
        for d in dictionaries:
            args = d.pop("args", [])
            series.append(GenericDataSeries(kw, *args, **d))
    return series

##############################################################################
# The public interface
##############################################################################

def _deprecation_msg_m_a_r_f(attr):
    sympy_deprecation_warning(
        f"The `{attr}` property is deprecated. The `{attr}` keyword "
        "argument should be passed to a plotting function, which generates "
        "the appropriate data series. If needed, index the plot object to "
        "retrieve a specific data series.",
        deprecated_since_version="1.13",
        active_deprecations_target="deprecated-markers-annotations-fill-rectangles",
        stacklevel=4)

class Plot:
    """Base class for all backends. A backend represents the plotting library,
    which implements the necessary functionalities in order to use SymPy
    plotting functions.

    For interactive work the function :func:`plot()` is better suited.

    This class permits the plotting of SymPy expressions using numerous
    backends (:external:mod:`matplotlib`, textplot, the old pyglet module for SymPy, Google
    charts api, etc).

    The figure can contain an arbitrary number of plots of SymPy expressions,
    lists of coordinates of points, etc. Plot has a private attribute _series that
    contains all data series to be plotted (expressions for lines or surfaces,
    lists of points, etc (all subclasses of BaseSeries)). Those data series are
    instances of classes not imported by ``from sympy import *``.

    The customization of the figure is on two levels. Global options that
    concern the figure as a whole (e.g. title, xlabel, scale, etc) and
    per-data series options (e.g. name) and aesthetics (e.g. color, point shape,
    line type, etc.).

    The difference between options and aesthetics is that an aesthetic can be
    a function of the coordinates (or parameters in a parametric plot). The
    supported values for an aesthetic are:

    - None (the backend uses default values)
    - a constant
    - a function of one variable (the first coordinate or parameter)
    - a function of two variables (the first and second coordinate or parameters)
    - a function of three variables (only in nonparametric 3D plots)

    Their implementation depends on the backend so they may not work in some
    backends.

    If the plot is parametric and the arity of the aesthetic function permits
    it the aesthetic is calculated over parameters and not over coordinates.
    If the arity does not permit calculation over parameters the calculation is
    done over coordinates.

    Only cartesian coordinates are supported for the moment, but you can use
    the parametric plots to plot in polar, spherical and cylindrical
    coordinates.

    The arguments for the constructor Plot must be subclasses of BaseSeries.

    Any global option can be specified as a keyword argument.

    The global options for a figure are:

    - title : str
    - xlabel : str or Symbol
    - ylabel : str or Symbol
    - zlabel : str or Symbol
    - legend : bool
    - xscale : {'linear', 'log'}
    - yscale : {'linear', 'log'}
    - axis : bool
    - axis_center : tuple of two floats or {'center', 'auto'}
    - xlim : tuple of two floats
    - ylim : tuple of two floats
    - aspect_ratio : tuple of two floats or {'auto'}
    - autoscale : bool
    - margin : float in [0, 1]
    - backend : {'default', 'matplotlib', 'text'} or a subclass of BaseBackend
    - size : optional tuple of two floats, (width, height); default: None

    The per data series options and aesthetics are:
    There are none in the base series. See below for options for subclasses.

    Some data series support additional aesthetics or options:

    :class:`~.LineOver1DRangeSeries`, :class:`~.Parametric2DLineSeries`, and
    :class:`~.Parametric3DLineSeries` support the following:

    Aesthetics:

    - line_color : string, or float, or function, optional
        Specifies the color for the plot, which depends on the backend being
        used.

        For example, if ``MatplotlibBackend`` is being used, then
        Matplotlib string colors are acceptable (``"red"``, ``"r"``,
        ``"cyan"``, ``"c"``, ...).
        Alternatively, we can use a float number, 0 < color < 1, wrapped in a
        string (for example, ``line_color="0.5"``) to specify grayscale colors.
        Alternatively, We can specify a function returning a single
        float value: this will be used to apply a color-loop (for example,
        ``line_color=lambda x: math.cos(x)``).

        Note that by setting line_color, it would be applied simultaneously
        to all the series.

    Options:

    - label : str
    - steps : bool
    - integers_only : bool

    :class:`~.SurfaceOver2DRangeSeries` and :class:`~.ParametricSurfaceSeries`
    support the following:

    Aesthetics:

    - surface_color : function which returns a float.

    Notes
    =====

    How the plotting module works:

    1. Whenever a plotting function is called, the provided expressions are
       processed and a list of instances of the :class:`BaseSeries` class is
       created, containing the necessary information to plot the expressions
       (e.g. the expression, ranges, series name, ...). Eventually, these
       objects will generate the numerical data to be plotted.
    2. A subclass of :class:`~.Plot` class is instantiaed (referred to as
       backend, from now on), which stores the list of series and the main
       attributes of the plot (e.g. axis labels, title, ...).
       The backend implements the logic to generate the actual figure with
       some plotting library.
    3. When the ``show`` command is executed, series are processed one by one
       to generate numerical data and add it to the figure. The backend is also
       going to set the axis labels, title, ..., according to the values stored
       in the Plot instance.

    The backend should check if it supports the data series that it is given
    (e.g. :class:`TextBackend` supports only :class:`LineOver1DRangeSeries`).

    It is the backend responsibility to know how to use the class of data series
    that it's given. Note that the current implementation of the ``*Series``
    classes is "matplotlib-centric": the numerical data returned by the
    ``get_points`` and ``get_meshes`` methods is meant to be used directly by
    Matplotlib. Therefore, the new backend will have to pre-process the
    numerical data to make it compatible with the chosen plotting library.
    Keep in mind that future SymPy versions may improve the ``*Series`` classes
    in order to return numerical data "non-matplotlib-centric", hence if you code
    a new backend you have the responsibility to check if its working on each
    SymPy release.

    Please explore the :class:`MatplotlibBackend` source code to understand
    how a backend should be coded.

    In order to be used by SymPy plotting functions, a backend must implement
    the following methods:

    * show(self): used to loop over the data series, generate the numerical
        data, plot it and set the axis labels, title, ...
    * save(self, path): used to save the current plot to the specified file
        path.
    * close(self): used to close the current plot backend (note: some plotting
        library does not support this functionality. In that case, just raise a
        warning).
    """

    def __init__(self, *args,
        title=None, xlabel=None, ylabel=None, zlabel=None, aspect_ratio='auto',
        xlim=None, ylim=None, axis_center='auto', axis=True,
        xscale='linear', yscale='linear', legend=False, autoscale=True,
        margin=0, annotations=None, markers=None, rectangles=None,
        fill=None, backend='default', size=None, **kwargs):

        # Options for the graph as a whole.
        # The possible values for each option are described in the docstring of
        # Plot. They are based purely on convention, no checking is done.
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.aspect_ratio = aspect_ratio
        self.axis_center = axis_center
        self.axis = axis
        self.xscale = xscale
        self.yscale = yscale
        self.legend = legend
        self.autoscale = autoscale
        self.margin = margin
        self._annotations = annotations
        self._markers = markers
        self._rectangles = rectangles
        self._fill = fill

        # Contains the data objects to be plotted. The backend should be smart
        # enough to iterate over this list.
        self._series = []
        self._series.extend(args)
        self._series.extend(_create_generic_data_series(
            annotations=annotations, markers=markers, rectangles=rectangles,
            fill=fill))

        is_real = \
            lambda lim: all(getattr(i, 'is_real', True) for i in lim)
        is_finite = \
            lambda lim: all(getattr(i, 'is_finite', True) for i in lim)

        # reduce code repetition
        def check_and_set(t_name, t):
            if t:
                if not is_real(t):
                    raise ValueError(
                    "All numbers from {}={} must be real".format(t_name, t))
                if not is_finite(t):
                    raise ValueError(
                    "All numbers from {}={} must be finite".format(t_name, t))
                setattr(self, t_name, (float(t[0]), float(t[1])))

        self.xlim = None
        check_and_set("xlim", xlim)
        self.ylim = None
        check_and_set("ylim", ylim)
        self.size = None
        check_and_set("size", size)

    @property
    def _backend(self):
        return self

    @property
    def backend(self):
        return type(self)

    def __str__(self):
        series_strs = [('[%d]: ' % i) + str(s)
                       for i, s in enumerate(self._series)]
        return 'Plot object containing:\n' + '\n'.join(series_strs)

    def __getitem__(self, index):
        return self._series[index]

    def __setitem__(self, index, *args):
        if len(args) == 1 and isinstance(args[0], BaseSeries):
            self._series[index] = args

    def __delitem__(self, index):
        del self._series[index]

    def append(self, arg):
        """Adds an element from a plot's series to an existing plot.

        Examples
        ========

        Consider two ``Plot`` objects, ``p1`` and ``p2``. To add the
        second plot's first series object to the first, use the
        ``append`` method, like so:

        .. plot::
           :format: doctest
           :include-source: True

           >>> from sympy import symbols
           >>> from sympy.plotting import plot
           >>> x = symbols('x')
           >>> p1 = plot(x*x, show=False)
           >>> p2 = plot(x, show=False)
           >>> p1.append(p2[0])
           >>> p1
           Plot object containing:
           [0]: cartesian line: x**2 for x over (-10.0, 10.0)
           [1]: cartesian line: x for x over (-10.0, 10.0)
           >>> p1.show()

        See Also
        ========

        extend

        """
        if isinstance(arg, BaseSeries):
            self._series.append(arg)
        else:
            raise TypeError('Must specify element of plot to append.')

    def extend(self, arg):
        """Adds all series from another plot.

        Examples
        ========

        Consider two ``Plot`` objects, ``p1`` and ``p2``. To add the
        second plot to the first, use the ``extend`` method, like so:

        .. plot::
           :format: doctest
           :include-source: True

           >>> from sympy import symbols
           >>> from sympy.plotting import plot
           >>> x = symbols('x')
           >>> p1 = plot(x**2, show=False)
           >>> p2 = plot(x, -x, show=False)
           >>> p1.extend(p2)
           >>> p1
           Plot object containing:
           [0]: cartesian line: x**2 for x over (-10.0, 10.0)
           [1]: cartesian line: x for x over (-10.0, 10.0)
           [2]: cartesian line: -x for x over (-10.0, 10.0)
           >>> p1.show()

        """
        if isinstance(arg, Plot):
            self._series.extend(arg._series)
        elif is_sequence(arg):
            self._series.extend(arg)
        else:
            raise TypeError('Expecting Plot or sequence of BaseSeries')

    def show(self):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    def close(self):
        raise NotImplementedError

    # deprecations

    @property
    def markers(self):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("markers")
        return self._markers

    @markers.setter
    def markers(self, v):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("markers")
        self._series.extend(_create_generic_data_series(markers=v))
        self._markers = v

    @property
    def annotations(self):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("annotations")
        return self._annotations

    @annotations.setter
    def annotations(self, v):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("annotations")
        self._series.extend(_create_generic_data_series(annotations=v))
        self._annotations = v

    @property
    def rectangles(self):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("rectangles")
        return self._rectangles

    @rectangles.setter
    def rectangles(self, v):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("rectangles")
        self._series.extend(_create_generic_data_series(rectangles=v))
        self._rectangles = v

    @property
    def fill(self):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("fill")
        return self._fill

    @fill.setter
    def fill(self, v):
        """.. deprecated:: 1.13"""
        _deprecation_msg_m_a_r_f("fill")
        self._series.extend(_create_generic_data_series(fill=v))
        self._fill = v


class PlotGrid:
    """This class helps to plot subplots from already created SymPy plots
    in a single figure.

    Examples
    ========

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> from sympy import symbols
        >>> from sympy.plotting import plot, plot3d, PlotGrid
        >>> x, y = symbols('x, y')
        >>> p1 = plot(x, x**2, x**3, (x, -5, 5))
        >>> p2 = plot((x**2, (x, -6, 6)), (x, (x, -5, 5)))
        >>> p3 = plot(x**3, (x, -5, 5))
        >>> p4 = plot3d(x*y, (x, -5, 5), (y, -5, 5))

    Plotting vertically in a single line:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> PlotGrid(2, 1, p1, p2)
        PlotGrid object containing:
        Plot[0]:Plot object containing:
        [0]: cartesian line: x for x over (-5.0, 5.0)
        [1]: cartesian line: x**2 for x over (-5.0, 5.0)
        [2]: cartesian line: x**3 for x over (-5.0, 5.0)
        Plot[1]:Plot object containing:
        [0]: cartesian line: x**2 for x over (-6.0, 6.0)
        [1]: cartesian line: x for x over (-5.0, 5.0)

    Plotting horizontally in a single line:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> PlotGrid(1, 3, p2, p3, p4)
        PlotGrid object containing:
        Plot[0]:Plot object containing:
        [0]: cartesian line: x**2 for x over (-6.0, 6.0)
        [1]: cartesian line: x for x over (-5.0, 5.0)
        Plot[1]:Plot object containing:
        [0]: cartesian line: x**3 for x over (-5.0, 5.0)
        Plot[2]:Plot object containing:
        [0]: cartesian surface: x*y for x over (-5.0, 5.0) and y over (-5.0, 5.0)

    Plotting in a grid form:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> PlotGrid(2, 2, p1, p2, p3, p4)
        PlotGrid object containing:
        Plot[0]:Plot object containing:
        [0]: cartesian line: x for x over (-5.0, 5.0)
        [1]: cartesian line: x**2 for x over (-5.0, 5.0)
        [2]: cartesian line: x**3 for x over (-5.0, 5.0)
        Plot[1]:Plot object containing:
        [0]: cartesian line: x**2 for x over (-6.0, 6.0)
        [1]: cartesian line: x for x over (-5.0, 5.0)
        Plot[2]:Plot object containing:
        [0]: cartesian line: x**3 for x over (-5.0, 5.0)
        Plot[3]:Plot object containing:
        [0]: cartesian surface: x*y for x over (-5.0, 5.0) and y over (-5.0, 5.0)

    """
    def __init__(self, nrows, ncolumns, *args, show=True, size=None, **kwargs):
        """
        Parameters
        ==========

        nrows :
            The number of rows that should be in the grid of the
            required subplot.
        ncolumns :
            The number of columns that should be in the grid
            of the required subplot.

        nrows and ncolumns together define the required grid.

        Arguments
        =========

        A list of predefined plot objects entered in a row-wise sequence
        i.e. plot objects which are to be in the top row of the required
        grid are written first, then the second row objects and so on

        Keyword arguments
        =================

        show : Boolean
            The default value is set to ``True``. Set show to ``False`` and
            the function will not display the subplot. The returned instance
            of the ``PlotGrid`` class can then be used to save or display the
            plot by calling the ``save()`` and ``show()`` methods
            respectively.
        size : (float, float), optional
            A tuple in the form (width, height) in inches to specify the size of
            the overall figure. The default value is set to ``None``, meaning
            the size will be set by the default backend.
        """
        self.matplotlib = import_module('matplotlib',
            import_kwargs={'fromlist': ['pyplot', 'cm', 'collections']},
            min_module_version='1.1.0', catch=(RuntimeError,))
        self.nrows = nrows
        self.ncolumns = ncolumns
        self._series = []
        self._fig = None
        self.args = args
        for arg in args:
            self._series.append(arg._series)
        self.size = size
        if show and self.matplotlib:
            self.show()

    def _create_figure(self):
        gs = self.matplotlib.gridspec.GridSpec(self.nrows, self.ncolumns)
        mapping = {}
        c = 0
        for i in range(self.nrows):
            for j in range(self.ncolumns):
                if c < len(self.args):
                    mapping[gs[i, j]] = self.args[c]
                c += 1

        kw = {} if not self.size else {"figsize": self.size}
        self._fig = self.matplotlib.pyplot.figure(**kw)
        for spec, p in mapping.items():
            kw = ({"projection": "3d"} if (len(p._series) > 0 and
                p._series[0].is_3D) else {})
            cur_ax = self._fig.add_subplot(spec, **kw)
            p._plotgrid_fig = self._fig
            p._plotgrid_ax = cur_ax
            p.process_series()

    @property
    def fig(self):
        if not self._fig:
            self._create_figure()
        return self._fig

    @property
    def _backend(self):
        return self

    def close(self):
        self.matplotlib.pyplot.close(self.fig)

    def show(self):
        if _show:
            self.fig.tight_layout()
            self.matplotlib.pyplot.show()
        else:
            self.close()

    def save(self, path):
        self.fig.savefig(path)

    def __str__(self):
        plot_strs = [('Plot[%d]:' % i) + str(plot)
                      for i, plot in enumerate(self.args)]

        return 'PlotGrid object containing:\n' + '\n'.join(plot_strs)


def plot_factory(*args, **kwargs):
    backend = kwargs.pop("backend", "default")
    if isinstance(backend, str):
        if backend == "default":
            matplotlib = import_module('matplotlib',
                min_module_version='1.1.0', catch=(RuntimeError,))
            if matplotlib:
                return MatplotlibBackend(*args, **kwargs)
            return TextBackend(*args, **kwargs)
        return plot_backends[backend](*args, **kwargs)
    elif (type(backend) == type) and issubclass(backend, Plot):
        return backend(*args, **kwargs)
    else:
        raise TypeError("backend must be either a string or a subclass of ``Plot``.")

##############################################################################
# Data Series
##############################################################################
#TODO more general way to calculate aesthetics (see get_color_array)

def _set_discretization_points(kwargs, pt):
    """Allow the use of the keyword arguments ``n, n1, n2`` to
    specify the number of discretization points in one and two
    directions, while keeping back-compatibility with older keyword arguments
    like, ``nb_of_points, nb_of_points_*, points``.

    Parameters
    ==========

    kwargs : dict
        Dictionary of keyword arguments passed into a plotting function.
    pt : type
        The type of the series, which indicates the kind of plot we are
        trying to create.
    """
    deprecated_keywords = {
        "nb_of_points": "n",
        "nb_of_points_x": "n1",
        "nb_of_points_y": "n2",
        "nb_of_points_u": "n1",
        "nb_of_points_v": "n2",
        "points": "n"
    }
    for k, v in deprecated_keywords.items():
        if k in kwargs.keys():
            kwargs[v] = kwargs.pop(k)

    if pt in [LineOver1DRangeSeries, Parametric2DLineSeries,
        Parametric3DLineSeries]:
        if "n1" in kwargs.keys():
            kwargs["n"] = kwargs["n1"]
    elif pt in [
        SurfaceOver2DRangeSeries, ContourSeries, ParametricSurfaceSeries]:
        if "n" in kwargs.keys():
            kwargs["n1"] = kwargs["n2"] = kwargs["n"]
    return kwargs

### The base class for all series
class BaseSeries:
    """Base class for the data objects containing stuff to be plotted.

    Explanation
    ===========

    The backend should check if it supports the data series that is given.
    (e.g. TextBackend supports only LineOver1DRangeSeries).
    It is the backend responsibility to know how to use the class of
    data series that is given.

    Some data series classes are grouped (using a class attribute like is_2Dline)
    according to the api they present (based only on convention). The backend is
    not obliged to use that api (e.g. LineOver1DRangeSeries belongs to the
    is_2Dline group and presents the get_points method, but the
    TextBackend does not use the get_points method).
    """

    # Some flags follow. The rationale for using flags instead of checking base
    # classes is that setting multiple flags is simpler than multiple
    # inheritance.

    is_2Dline = False
    # Some of the backends expect:
    #  - get_points returning 1D np.arrays list_x, list_y
    #  - get_color_array returning 1D np.array (done in Line2DBaseSeries)
    # with the colors calculated at the points from get_points

    is_3Dline = False
    # Some of the backends expect:
    #  - get_points returning 1D np.arrays list_x, list_y, list_y
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

    is_implicit = False
    # Some of the backends expect:
    #   - get_meshes returning mesh_x (1D array), mesh_y(1D array,
    #     mesh_z (2D np.arrays)
    #   - get_points an alias for get_meshes
    # Different from is_contour as the colormap in backend will be
    # different

    is_parametric = False
    # The calculation of aesthetics expects:
    #   - get_parameter_points returning one or two np.arrays (1D or 2D)
    # used for calculation aesthetics

    is_generic = False
    # Represent generic user-provided numerical data

    def __init__(self):
        super().__init__()

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
        super().__init__()
        self.label = None
        self.steps = False
        self.only_integers = False
        self.line_color = None

    def get_data(self):
        """ Return lists of coordinates for plotting the line.

        Returns
        =======
            x : list
                List of x-coordinates

            y : list
                List of y-coordinates

            z : list
                List of z-coordinates in case of Parametric3DLineSeries
        """
        np = import_module('numpy')
        points = self._get_data_helper()
        if self.steps is True:
            if len(points) == 2:
                x = np.array((points[0], points[0])).T.flatten()[1:]
                y = np.array((points[1], points[1])).T.flatten()[:-1]
                points = (x, y)
            else:
                x = np.repeat(points[0], 3)[2:]
                y = np.repeat(points[1], 3)[:-2]
                z = np.repeat(points[2], 3)[1:-1]
                points = (x, y, z)
        return points

    def get_segments(self):
        sympy_deprecation_warning(
            """
            The Line2DBaseSeries.get_segments() method is deprecated.

            Instead, use the MatplotlibBackend.get_segments() method, or use
            The get_points() or get_data() methods.
            """,
            deprecated_since_version="1.9",
            active_deprecations_target="deprecated-get-segments")

        np = import_module('numpy')
        points = type(self).get_data(self)
        points = np.ma.array(points).T.reshape(-1, 1, self._dim)
        return np.ma.concatenate([points[:-1], points[1:]], axis=1)

    def get_color_array(self):
        np = import_module('numpy')
        c = self.line_color
        if hasattr(c, '__call__'):
            f = np.vectorize(c)
            nargs = arity(c)
            if nargs == 1 and self.is_parametric:
                x = self.get_parameter_points()
                return f(centers_of_segments(x))
            else:
                variables = list(map(centers_of_segments, self.get_points()))
                if nargs == 1:
                    return f(variables[0])
                elif nargs == 2:
                    return f(*variables[:2])
                else:  # only if the line is 3D (otherwise raises an error)
                    return f(*variables)
        else:
            return c*np.ones(self.nb_of_points)


class List2DSeries(Line2DBaseSeries):
    """Representation for a line consisting of list of points."""

    def __init__(self, list_x, list_y):
        np = import_module('numpy')
        super().__init__()
        self.list_x = np.array(list_x)
        self.list_y = np.array(list_y)
        self.label = 'list'

    def __str__(self):
        return 'list plot'

    def get_points(self):
        return (self.list_x, self.list_y)


class LineOver1DRangeSeries(Line2DBaseSeries):
    """Representation for a line consisting of a SymPy expression over a range."""

    def __init__(self, expr, var_start_end, **kwargs):
        super().__init__()
        self.expr = sympify(expr)
        self.label = kwargs.get('label', None) or self.expr
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])
        self.nb_of_points = kwargs.get('n', 300)
        self.adaptive = kwargs.get('adaptive', True)
        self.depth = kwargs.get('depth', 12)
        self.line_color = kwargs.get('line_color', None)
        self.xscale = kwargs.get('xscale', 'linear')

    def __str__(self):
        return 'cartesian line: %s for %s over %s' % (
            str(self.expr), str(self.var), str((self.start, self.end)))

    def get_points(self):
        """Return lists of coordinates for plotting. Depending on the
        ``adaptive`` option, this function will either use an adaptive algorithm
        or it will uniformly sample the expression over the provided range.

        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.

        Returns
        =======
            x : list
                List of x-coordinates

            y : list
                List of y-coordinates
        """
        return self._get_data_helper()

    def _get_data_helper(self):
        if self.only_integers or not self.adaptive:
            return self._uniform_sampling()
        return self._adaptive_sampling()

    def _adaptive_sampling(self):
        """The adaptive sampling is done by recursively checking if three
        points are almost collinear. If they are not collinear, then more
        points are added between those points.

        References
        ==========

        .. [1] Adaptive polygonal approximation of parametric curves,
               Luiz Henrique de Figueiredo.
        """
        f = lambdify([self.var], self.expr)
        x_coords = []
        y_coords = []
        np = import_module('numpy')
        def sample(p, q, depth):
            """ Samples recursively if three points are almost collinear.
            For depth < 6, points are added irrespective of whether they
            satisfy the collinearity condition or not. The maximum depth
            allowed is 12.
            """
            # Randomly sample to avoid aliasing.
            random = 0.45 + np.random.rand() * 0.1
            if self.xscale == 'log':
                xnew = 10**(np.log10(p[0]) + random * (np.log10(q[0]) -
                                                        np.log10(p[0])))
            else:
                xnew = p[0] + random * (q[0] - p[0])
            ynew = f(xnew)
            new_point = np.array([xnew, ynew])

            # Maximum depth
            if depth > self.depth:
                x_coords.append(q[0])
                y_coords.append(q[1])

            # Sample to depth of 6 (whether the line is flat or not)
            # without using linspace (to avoid aliasing).
            elif depth < 6:
                sample(p, new_point, depth + 1)
                sample(new_point, q, depth + 1)

            # Sample ten points if complex values are encountered
            # at both ends. If there is a real value in between, then
            # sample those points further.
            elif p[1] is None and q[1] is None:
                if self.xscale == 'log':
                    xarray = np.logspace(p[0], q[0], 10)
                else:
                    xarray = np.linspace(p[0], q[0], 10)
                yarray = list(map(f, xarray))
                if not all(y is None for y in yarray):
                    for i in range(len(yarray) - 1):
                        if not (yarray[i] is None and yarray[i + 1] is None):
                            sample([xarray[i], yarray[i]],
                                [xarray[i + 1], yarray[i + 1]], depth + 1)

            # Sample further if one of the end points in None (i.e. a
            # complex value) or the three points are not almost collinear.
            elif (p[1] is None or q[1] is None or new_point[1] is None
                    or not flat(p, new_point, q)):
                sample(p, new_point, depth + 1)
                sample(new_point, q, depth + 1)
            else:
                x_coords.append(q[0])
                y_coords.append(q[1])

        f_start = f(self.start)
        f_end = f(self.end)
        x_coords.append(self.start)
        y_coords.append(f_start)
        sample(np.array([self.start, f_start]),
                np.array([self.end, f_end]), 0)

        return (x_coords, y_coords)

    def _uniform_sampling(self):
        np = import_module('numpy')
        if self.only_integers is True:
            if self.xscale == 'log':
                list_x = np.logspace(int(self.start), int(self.end),
                        num=int(self.end) - int(self.start) + 1)
            else:
                list_x = np.linspace(int(self.start), int(self.end),
                    num=int(self.end) - int(self.start) + 1)
        else:
            if self.xscale == 'log':
                list_x = np.logspace(self.start, self.end, num=self.nb_of_points)
            else:
                list_x = np.linspace(self.start, self.end, num=self.nb_of_points)
        f = vectorized_lambdify([self.var], self.expr)
        list_y = f(list_x)
        return (list_x, list_y)


class Parametric2DLineSeries(Line2DBaseSeries):
    """Representation for a line consisting of two parametric SymPy expressions
    over a range."""

    is_parametric = True

    def __init__(self, expr_x, expr_y, var_start_end, **kwargs):
        super().__init__()
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.label = kwargs.get('label', None) or \
                            Tuple(self.expr_x, self.expr_y)
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])
        self.nb_of_points = kwargs.get('n', 300)
        self.adaptive = kwargs.get('adaptive', True)
        self.depth = kwargs.get('depth', 12)
        self.line_color = kwargs.get('line_color', None)

    def __str__(self):
        return 'parametric cartesian line: (%s, %s) for %s over %s' % (
            str(self.expr_x), str(self.expr_y), str(self.var),
            str((self.start, self.end)))

    def get_parameter_points(self):
        np = import_module('numpy')
        return np.linspace(self.start, self.end, num=self.nb_of_points)

    def get_points(self):
        """ Return lists of coordinates for plotting. Depending on the
        ``adaptive`` option, this function will either use an adaptive algorithm
        or it will uniformly sample the expression over the provided range.

        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.

        Returns
        =======
            x : list
                List of x-coordinates

            y : list
                List of y-coordinates

        """
        return self._get_data_helper()[:-1]

    def _get_data_helper(self):
        if not self.adaptive:
            return self._uniform_sampling()
        return self._adaptive_sampling()

    def _uniform_sampling(self):
        param = self.get_parameter_points()
        fx = vectorized_lambdify([self.var], self.expr_x)
        fy = vectorized_lambdify([self.var], self.expr_y)
        list_x = fx(param)
        list_y = fy(param)
        return (list_x, list_y, param)

    def _adaptive_sampling(self):
        """The adaptive sampling is done by recursively checking if three
        points are almost collinear. If they are not collinear, then more
        points are added between those points.

        References
        ==========

        .. [1] Adaptive polygonal approximation of parametric curves,
            Luiz Henrique de Figueiredo.
        """
        f_x = lambdify([self.var], self.expr_x)
        f_y = lambdify([self.var], self.expr_y)
        x_coords = []
        y_coords = []
        param = []

        def sample(param_p, param_q, p, q, depth):
            """ Samples recursively if three points are almost collinear.
            For depth < 6, points are added irrespective of whether they
            satisfy the collinearity condition or not. The maximum depth
            allowed is 12.
            """
            # Randomly sample to avoid aliasing.
            np = import_module('numpy')
            random = 0.45 + np.random.rand() * 0.1
            param_new = param_p + random * (param_q - param_p)
            xnew = f_x(param_new)
            ynew = f_y(param_new)
            new_point = np.array([xnew, ynew])

            # Maximum depth
            if depth > self.depth:
                x_coords.append(q[0])
                y_coords.append(q[1])
                param.append(param_p)

            # Sample irrespective of whether the line is flat till the
            # depth of 6. We are not using linspace to avoid aliasing.
            elif depth < 6:
                sample(param_p, param_new, p, new_point, depth + 1)
                sample(param_new, param_q, new_point, q, depth + 1)

            # Sample ten points if complex values are encountered
            # at both ends. If there is a real value in between, then
            # sample those points further.
            elif ((p[0] is None and q[1] is None) or
                    (p[1] is None and q[1] is None)):
                param_array = np.linspace(param_p, param_q, 10)
                x_array = list(map(f_x, param_array))
                y_array = list(map(f_y, param_array))
                if not all(x is None and y is None
                           for x, y in zip(x_array, y_array)):
                    for i in range(len(y_array) - 1):
                        if ((x_array[i] is not None and y_array[i] is not None) or
                                (x_array[i + 1] is not None and y_array[i + 1] is not None)):
                            point_a = [x_array[i], y_array[i]]
                            point_b = [x_array[i + 1], y_array[i + 1]]
                            sample(param_array[i], param_array[i], point_a,
                                   point_b, depth + 1)

            # Sample further if one of the end points in None (i.e. a complex
            # value) or the three points are not almost collinear.
            elif (p[0] is None or p[1] is None
                    or q[1] is None or q[0] is None
                    or not flat(p, new_point, q)):
                sample(param_p, param_new, p, new_point, depth + 1)
                sample(param_new, param_q, new_point, q, depth + 1)
            else:
                x_coords.append(q[0])
                y_coords.append(q[1])
                param.append(param_p)

        f_start_x = f_x(self.start)
        f_start_y = f_y(self.start)
        start = [f_start_x, f_start_y]
        f_end_x = f_x(self.end)
        f_end_y = f_y(self.end)
        end = [f_end_x, f_end_y]
        x_coords.append(f_start_x)
        y_coords.append(f_start_y)
        param.append(start)
        sample(self.start, self.end, start, end, 0)

        return x_coords, y_coords, param


### 3D lines
class Line3DBaseSeries(Line2DBaseSeries):
    """A base class for 3D lines.

    Most of the stuff is derived from Line2DBaseSeries."""

    is_2Dline = False
    is_3Dline = True
    _dim = 3

    def __init__(self):
        super().__init__()


class Parametric3DLineSeries(Line3DBaseSeries):
    """Representation for a 3D line consisting of three parametric SymPy
    expressions and a range."""

    is_parametric = True

    def __init__(self, expr_x, expr_y, expr_z, var_start_end, **kwargs):
        super().__init__()
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.expr_z = sympify(expr_z)
        self.label = kwargs.get('label', None) or \
                        Tuple(self.expr_x, self.expr_y)
        self.var = sympify(var_start_end[0])
        self.start = float(var_start_end[1])
        self.end = float(var_start_end[2])
        self.nb_of_points = kwargs.get('n', 300)
        self.line_color = kwargs.get('line_color', None)
        self._xlim = None
        self._ylim = None
        self._zlim = None

    def __str__(self):
        return '3D parametric cartesian line: (%s, %s, %s) for %s over %s' % (
            str(self.expr_x), str(self.expr_y), str(self.expr_z),
            str(self.var), str((self.start, self.end)))

    def get_parameter_points(self):
        np = import_module('numpy')
        return np.linspace(self.start, self.end, num=self.nb_of_points)

    def get_points(self):
        """Return the x,y,z coordinates for plotting the line.
        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.
        """
        return self.get_data()[:-1]

    def get_data(self):
        """Return coordinates for plotting the line.

        Returns
        =======

        x: np.ndarray
        y: np.ndarray
        z: np.ndarray
        param : np.ndarray
        """
        np = import_module('numpy')
        param = self.get_parameter_points()
        fx = vectorized_lambdify([self.var], self.expr_x)
        fy = vectorized_lambdify([self.var], self.expr_y)
        fz = vectorized_lambdify([self.var], self.expr_z)

        list_x = fx(param)
        list_y = fy(param)
        list_z = fz(param)

        list_x = np.array(list_x, dtype=np.float64)
        list_y = np.array(list_y, dtype=np.float64)
        list_z = np.array(list_z, dtype=np.float64)

        list_x = np.ma.masked_invalid(list_x)
        list_y = np.ma.masked_invalid(list_y)
        list_z = np.ma.masked_invalid(list_z)

        self._xlim = (np.amin(list_x), np.amax(list_x))
        self._ylim = (np.amin(list_y), np.amax(list_y))
        self._zlim = (np.amin(list_z), np.amax(list_z))
        return list_x, list_y, list_z, param


### Surfaces
class SurfaceBaseSeries(BaseSeries):
    """A base class for 3D surfaces."""

    is_3Dsurface = True

    def __init__(self):
        super().__init__()
        self.surface_color = None

    def get_color_array(self):
        np = import_module('numpy')
        c = self.surface_color
        if isinstance(c, Callable):
            f = np.vectorize(c)
            nargs = arity(c)
            if self.is_parametric:
                variables = list(map(centers_of_faces, self.get_parameter_meshes()))
                if nargs == 1:
                    return f(variables[0])
                elif nargs == 2:
                    return f(*variables)
            variables = list(map(centers_of_faces, self.get_meshes()))
            if nargs == 1:
                return f(variables[0])
            elif nargs == 2:
                return f(*variables[:2])
            else:
                return f(*variables)
        else:
            if isinstance(self, SurfaceOver2DRangeSeries):
                return c*np.ones(min(self.nb_of_points_x, self.nb_of_points_y))
            else:
                return c*np.ones(min(self.nb_of_points_u, self.nb_of_points_v))


class SurfaceOver2DRangeSeries(SurfaceBaseSeries):
    """Representation for a 3D surface consisting of a SymPy expression and 2D
    range."""
    def __init__(self, expr, var_start_end_x, var_start_end_y, **kwargs):
        super().__init__()
        self.expr = sympify(expr)
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])
        self.nb_of_points_x = kwargs.get('n1', 50)
        self.nb_of_points_y = kwargs.get('n2', 50)
        self.surface_color = kwargs.get('surface_color', None)

        self._xlim = (self.start_x, self.end_x)
        self._ylim = (self.start_y, self.end_y)

    def __str__(self):
        return ('cartesian surface: %s for'
                ' %s over %s and %s over %s') % (
                    str(self.expr),
                    str(self.var_x),
                    str((self.start_x, self.end_x)),
                    str(self.var_y),
                    str((self.start_y, self.end_y)))

    def get_meshes(self):
        """Return the x,y,z coordinates for plotting the surface.
        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.
        """
        return self.get_data()

    def get_data(self):
        """Return arrays of coordinates for plotting.

        Returns
        =======
        mesh_x : np.ndarray
            Discretized x-domain.
        mesh_y : np.ndarray
            Discretized y-domain.
        mesh_z : np.ndarray
            Results of the evaluation.
        """
        np = import_module('numpy')
        mesh_x, mesh_y = np.meshgrid(np.linspace(self.start_x, self.end_x,
                                                 num=self.nb_of_points_x),
                                     np.linspace(self.start_y, self.end_y,
                                                 num=self.nb_of_points_y))
        f = vectorized_lambdify((self.var_x, self.var_y), self.expr)
        mesh_z = f(mesh_x, mesh_y)
        mesh_z = np.array(mesh_z, dtype=np.float64)
        mesh_z = np.ma.masked_invalid(mesh_z)
        self._zlim = (np.amin(mesh_z), np.amax(mesh_z))
        return mesh_x, mesh_y, mesh_z


class ParametricSurfaceSeries(SurfaceBaseSeries):
    """Representation for a 3D surface consisting of three parametric SymPy
    expressions and a range."""

    is_parametric = True

    def __init__(
        self, expr_x, expr_y, expr_z, var_start_end_u, var_start_end_v,
            **kwargs):
        super().__init__()
        self.expr_x = sympify(expr_x)
        self.expr_y = sympify(expr_y)
        self.expr_z = sympify(expr_z)
        self.var_u = sympify(var_start_end_u[0])
        self.start_u = float(var_start_end_u[1])
        self.end_u = float(var_start_end_u[2])
        self.var_v = sympify(var_start_end_v[0])
        self.start_v = float(var_start_end_v[1])
        self.end_v = float(var_start_end_v[2])
        self.nb_of_points_u = kwargs.get('n1', 50)
        self.nb_of_points_v = kwargs.get('n2', 50)
        self.surface_color = kwargs.get('surface_color', None)

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
        np = import_module('numpy')
        return np.meshgrid(np.linspace(self.start_u, self.end_u,
                                       num=self.nb_of_points_u),
                           np.linspace(self.start_v, self.end_v,
                                       num=self.nb_of_points_v))

    def get_meshes(self):
        """Return the x,y,z coordinates for plotting the surface.
        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.
        """
        return self.get_data()[:3]

    def get_data(self):
        """Return arrays of coordinates for plotting.

        Returns
        =======
        mesh_x : np.ndarray
            Discretized x-domain.
        mesh_y : np.ndarray
            Discretized y-domain.
        mesh_z : np.ndarray
            Results of the evaluation.
        mesh_u : np.ndarray
            Discretized u range.
        mesh_v : np.ndarray
            Discretized v range.
        """
        np = import_module('numpy')

        mesh_u, mesh_v = self.get_parameter_meshes()
        fx = vectorized_lambdify((self.var_u, self.var_v), self.expr_x)
        fy = vectorized_lambdify((self.var_u, self.var_v), self.expr_y)
        fz = vectorized_lambdify((self.var_u, self.var_v), self.expr_z)

        mesh_x = fx(mesh_u, mesh_v)
        mesh_y = fy(mesh_u, mesh_v)
        mesh_z = fz(mesh_u, mesh_v)

        mesh_x = np.array(mesh_x, dtype=np.float64)
        mesh_y = np.array(mesh_y, dtype=np.float64)
        mesh_z = np.array(mesh_z, dtype=np.float64)

        mesh_x = np.ma.masked_invalid(mesh_x)
        mesh_y = np.ma.masked_invalid(mesh_y)
        mesh_z = np.ma.masked_invalid(mesh_z)

        self._xlim = (np.amin(mesh_x), np.amax(mesh_x))
        self._ylim = (np.amin(mesh_y), np.amax(mesh_y))
        self._zlim = (np.amin(mesh_z), np.amax(mesh_z))

        return mesh_x, mesh_y, mesh_z, mesh_u, mesh_v


### Contours
class ContourSeries(BaseSeries):
    """Representation for a contour plot."""
    # The code is mostly repetition of SurfaceOver2DRange.
    # Presently used in contour_plot function

    is_contour = True

    def __init__(self, expr, var_start_end_x, var_start_end_y):
        super().__init__()
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

        self._xlim = (self.start_x, self.end_x)
        self._ylim = (self.start_y, self.end_y)

    def __str__(self):
        return ('contour: %s for '
                '%s over %s and %s over %s') % (
                    str(self.expr),
                    str(self.var_x),
                    str((self.start_x, self.end_x)),
                    str(self.var_y),
                    str((self.start_y, self.end_y)))

    def get_meshes(self):
        """Return the x,y,z coordinates for plotting the surface.
        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.
        """
        return self.get_data()[:3]

    def get_data(self):
        """Return arrays of coordinates for plotting.

        Returns
        =======
        mesh_x : np.ndarray
            Discretized x-domain.
        mesh_y : np.ndarray
            Discretized y-domain.
        mesh_z : np.ndarray
            Results of the evaluation.
        """
        np = import_module('numpy')
        mesh_x, mesh_y = np.meshgrid(np.linspace(self.start_x, self.end_x,
                                                 num=self.nb_of_points_x),
                                     np.linspace(self.start_y, self.end_y,
                                                 num=self.nb_of_points_y))
        f = vectorized_lambdify((self.var_x, self.var_y), self.expr)
        return (mesh_x, mesh_y, f(mesh_x, mesh_y))


class GenericDataSeries(BaseSeries):
    """Represents generic numerical data.

    Notes
    =====
    This class serves the purpose of back-compatibility with the "markers,
    annotations, fill, rectangles" keyword arguments that represent
    user-provided numerical data. In particular, it solves the problem of
    combining together two or more plot-objects with the ``extend`` or
    ``append`` methods: user-provided numerical data is also taken into
    consideration because it is stored in this series class.

    Also note that the current implementation is far from optimal, as each
    keyword argument is stored into an attribute in the ``Plot`` class, which
    requires a hard-coded if-statement in the ``MatplotlibBackend`` class.
    The implementation suggests that it is ok to add attributes and
    if-statements to provide more and more functionalities for user-provided
    numerical data (e.g. adding horizontal lines, or vertical lines, or bar
    plots, etc). However, in doing so one would reinvent the wheel: plotting
    libraries (like Matplotlib) already implements the necessary API.

    Instead of adding more keyword arguments and attributes, users interested
    in adding custom numerical data to a plot should retrieve the figure
    created by this plotting module. For example, this code:

    .. plot::
       :context: close-figs
       :include-source: True

       from sympy import Symbol, plot, cos
       x = Symbol("x")
       p = plot(cos(x), markers=[{"args": [[0, 1, 2], [0, 1, -1], "*"]}])

    Becomes:

    .. plot::
       :context: close-figs
       :include-source: True

       p = plot(cos(x), backend="matplotlib")
       fig, ax = p._backend.fig, p._backend.ax[0]
       ax.plot([0, 1, 2], [0, 1, -1], "*")
       fig

    Which is far better in terms of readibility. Also, it gives access to the
    full plotting library capabilities, without the need to reinvent the wheel.
    """
    is_generic = True

    def __init__(self, tp, *args, **kwargs):
        self.type = tp
        self.args = args
        self.rendering_kw = kwargs

    def get_data(self):
        return self.args

##############################################################################
# Backends
##############################################################################


# Don't have to check for the success of importing matplotlib in each case;
# we will only be using this backend if we can successfully import matploblib
class MatplotlibBackend(Plot):
    """ This class implements the functionalities to use Matplotlib with SymPy
    plotting functions.
    """

    def __init__(self, *series, **kwargs):
        super().__init__(*series, **kwargs)
        self.matplotlib = import_module('matplotlib',
            import_kwargs={'fromlist': ['pyplot', 'cm', 'collections']},
            min_module_version='1.1.0', catch=(RuntimeError,))
        self.plt = self.matplotlib.pyplot
        self.cm = self.matplotlib.cm
        self.LineCollection = self.matplotlib.collections.LineCollection
        self.aspect = kwargs.get('aspect_ratio', 'auto')
        if self.aspect != 'auto':
            self.aspect = float(self.aspect[1]) / self.aspect[0]
        # PlotGrid can provide its figure and axes to be populated with
        # the data from the series.
        self._plotgrid_fig = kwargs.pop("fig", None)
        self._plotgrid_ax = kwargs.pop("ax", None)

    def _create_figure(self):
        def set_spines(ax):
            ax.spines['left'].set_position('zero')
            ax.spines['right'].set_color('none')
            ax.spines['bottom'].set_position('zero')
            ax.spines['top'].set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

        if self._plotgrid_fig is not None:
            self.fig = self._plotgrid_fig
            self.ax = self._plotgrid_ax
            if not any(s.is_3D for s in self._series):
                set_spines(self.ax)
        else:
            self.fig = self.plt.figure(figsize=self.size)
            if any(s.is_3D for s in self._series):
                self.ax = self.fig.add_subplot(1, 1, 1, projection="3d")
            else:
                self.ax = self.fig.add_subplot(1, 1, 1)
                set_spines(self.ax)

    @staticmethod
    def get_segments(x, y, z=None):
        """ Convert two list of coordinates to a list of segments to be used
        with Matplotlib's :external:class:`~matplotlib.collections.LineCollection`.

        Parameters
        ==========
            x : list
                List of x-coordinates

            y : list
                List of y-coordinates

            z : list
                List of z-coordinates for a 3D line.
        """
        np = import_module('numpy')
        if z is not None:
            dim = 3
            points = (x, y, z)
        else:
            dim = 2
            points = (x, y)
        points = np.ma.array(points).T.reshape(-1, 1, dim)
        return np.ma.concatenate([points[:-1], points[1:]], axis=1)

    def _process_series(self, series, ax):
        np = import_module('numpy')
        mpl_toolkits = import_module(
            'mpl_toolkits', import_kwargs={'fromlist': ['mplot3d']})

        # XXX Workaround for matplotlib issue
        # https://github.com/matplotlib/matplotlib/issues/17130
        xlims, ylims, zlims = [], [], []

        for s in series:
            # Create the collections
            if s.is_2Dline:
                if s.is_parametric:
                    x, y, param = s.get_data()
                else:
                    x, y = s.get_data()
                if (isinstance(s.line_color, (int, float)) or
                        callable(s.line_color)):
                    segments = self.get_segments(x, y)
                    collection = self.LineCollection(segments)
                    collection.set_array(s.get_color_array())
                    ax.add_collection(collection)
                else:
                    lbl = _str_or_latex(s.label)
                    line, = ax.plot(x, y, label=lbl, color=s.line_color)
            elif s.is_contour:
                ax.contour(*s.get_data())
            elif s.is_3Dline:
                x, y, z, param = s.get_data()
                if (isinstance(s.line_color, (int, float)) or
                        callable(s.line_color)):
                    art3d = mpl_toolkits.mplot3d.art3d
                    segments = self.get_segments(x, y, z)
                    collection = art3d.Line3DCollection(segments)
                    collection.set_array(s.get_color_array())
                    ax.add_collection(collection)
                else:
                    lbl = _str_or_latex(s.label)
                    ax.plot(x, y, z, label=lbl, color=s.line_color)

                xlims.append(s._xlim)
                ylims.append(s._ylim)
                zlims.append(s._zlim)
            elif s.is_3Dsurface:
                if s.is_parametric:
                    x, y, z, u, v = s.get_data()
                else:
                    x, y, z = s.get_data()
                collection = ax.plot_surface(x, y, z,
                    cmap=getattr(self.cm, 'viridis', self.cm.jet),
                    rstride=1, cstride=1, linewidth=0.1)
                if isinstance(s.surface_color, (float, int, Callable)):
                    color_array = s.get_color_array()
                    color_array = color_array.reshape(color_array.size)
                    collection.set_array(color_array)
                else:
                    collection.set_color(s.surface_color)

                xlims.append(s._xlim)
                ylims.append(s._ylim)
                zlims.append(s._zlim)
            elif s.is_implicit:
                points = s.get_data()
                if len(points) == 2:
                    # interval math plotting
                    x, y = _matplotlib_list(points[0])
                    ax.fill(x, y, facecolor=s.line_color, edgecolor='None')
                else:
                    # use contourf or contour depending on whether it is
                    # an inequality or equality.
                    # XXX: ``contour`` plots multiple lines. Should be fixed.
                    ListedColormap = self.matplotlib.colors.ListedColormap
                    colormap = ListedColormap(["white", s.line_color])
                    xarray, yarray, zarray, plot_type = points
                    if plot_type == 'contour':
                        ax.contour(xarray, yarray, zarray, cmap=colormap)
                    else:
                        ax.contourf(xarray, yarray, zarray, cmap=colormap)
            elif s.is_generic:
                if s.type == "markers":
                    # s.rendering_kw["color"] = s.line_color
                    ax.plot(*s.args, **s.rendering_kw)
                elif s.type == "annotations":
                    ax.annotate(*s.args, **s.rendering_kw)
                elif s.type == "fill":
                    # s.rendering_kw["color"] = s.line_color
                    ax.fill_between(*s.args, **s.rendering_kw)
                elif s.type == "rectangles":
                    # s.rendering_kw["color"] = s.line_color
                    ax.add_patch(
                        self.matplotlib.patches.Rectangle(
                            *s.args, **s.rendering_kw))
            else:
                raise NotImplementedError(
                    '{} is not supported in the SymPy plotting module '
                    'with matplotlib backend. Please report this issue.'
                    .format(ax))

        Axes3D = mpl_toolkits.mplot3d.Axes3D
        if not isinstance(ax, Axes3D):
            ax.autoscale_view(
                scalex=ax.get_autoscalex_on(),
                scaley=ax.get_autoscaley_on())
        else:
            # XXX Workaround for matplotlib issue
            # https://github.com/matplotlib/matplotlib/issues/17130
            if xlims:
                xlims = np.array(xlims)
                xlim = (np.amin(xlims[:, 0]), np.amax(xlims[:, 1]))
                ax.set_xlim(xlim)
            else:
                ax.set_xlim([0, 1])

            if ylims:
                ylims = np.array(ylims)
                ylim = (np.amin(ylims[:, 0]), np.amax(ylims[:, 1]))
                ax.set_ylim(ylim)
            else:
                ax.set_ylim([0, 1])

            if zlims:
                zlims = np.array(zlims)
                zlim = (np.amin(zlims[:, 0]), np.amax(zlims[:, 1]))
                ax.set_zlim(zlim)
            else:
                ax.set_zlim([0, 1])

        # Set global options.
        # TODO The 3D stuff
        # XXX The order of those is important.
        if self.xscale and not isinstance(ax, Axes3D):
            ax.set_xscale(self.xscale)
        if self.yscale and not isinstance(ax, Axes3D):
            ax.set_yscale(self.yscale)
        if not isinstance(ax, Axes3D) or self.matplotlib.__version__ >= '1.2.0':  # XXX in the distant future remove this check
            ax.set_autoscale_on(self.autoscale)
        if self.axis_center:
            val = self.axis_center
            if isinstance(ax, Axes3D):
                pass
            elif val == 'center':
                ax.spines['left'].set_position('center')
                ax.spines['bottom'].set_position('center')
            elif val == 'auto':
                xl, xh = ax.get_xlim()
                yl, yh = ax.get_ylim()
                pos_left = ('data', 0) if xl*xh <= 0 else 'center'
                pos_bottom = ('data', 0) if yl*yh <= 0 else 'center'
                ax.spines['left'].set_position(pos_left)
                ax.spines['bottom'].set_position(pos_bottom)
            else:
                ax.spines['left'].set_position(('data', val[0]))
                ax.spines['bottom'].set_position(('data', val[1]))
        if not self.axis:
            ax.set_axis_off()
        if self.legend:
            if ax.legend():
                ax.legend_.set_visible(self.legend)
        if self.margin:
            ax.set_xmargin(self.margin)
            ax.set_ymargin(self.margin)
        if self.title:
            ax.set_title(self.title)
        if self.xlabel:
            xlbl = _str_or_latex(self.xlabel)
            ax.set_xlabel(xlbl, position=(1, 0))
        if self.ylabel:
            ylbl = _str_or_latex(self.ylabel)
            ax.set_ylabel(ylbl, position=(0, 1))
        if isinstance(ax, Axes3D) and self.zlabel:
            zlbl = _str_or_latex(self.zlabel)
            ax.set_zlabel(zlbl, position=(0, 1))

        # xlim and ylim should always be set at last so that plot limits
        # doesn't get altered during the process.
        if self.xlim:
            ax.set_xlim(self.xlim)
        if self.ylim:
            ax.set_ylim(self.ylim)
        self.ax.set_aspect(self.aspect)


    def process_series(self):
        """
        Iterates over every ``Plot`` object and further calls
        _process_series()
        """
        self._create_figure()
        self._process_series(self._series, self.ax)

    def show(self):
        self.process_series()
        #TODO after fixing https://github.com/ipython/ipython/issues/1255
        # you can uncomment the next line and remove the pyplot.show() call
        #self.fig.show()
        if _show:
            self.fig.tight_layout()
            self.plt.show()
        else:
            self.close()

    def save(self, path):
        self.process_series()
        self.fig.savefig(path)

    def close(self):
        self.plt.close(self.fig)


class TextBackend(Plot):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def show(self):
        if not _show:
            return
        if len(self._series) != 1:
            raise ValueError(
                'The TextBackend supports only one graph per Plot.')
        elif not isinstance(self._series[0], LineOver1DRangeSeries):
            raise ValueError(
                'The TextBackend supports only expressions over a 1D range')
        else:
            ser = self._series[0]
            textplot(ser.expr, ser.start, ser.end)

    def close(self):
        pass


plot_backends = {
    'matplotlib': MatplotlibBackend,
    'text': TextBackend,
}


##############################################################################
# Finding the centers of line segments or mesh faces
##############################################################################

def centers_of_segments(array):
    np = import_module('numpy')
    return np.mean(np.vstack((array[:-1], array[1:])), 0)


def centers_of_faces(array):
    np = import_module('numpy')
    return np.mean(np.dstack((array[:-1, :-1],
                             array[1:, :-1],
                             array[:-1, 1:],
                             array[:-1, :-1],
                             )), 2)


def flat(x, y, z, eps=1e-3):
    """Checks whether three points are almost collinear"""
    np = import_module('numpy')
    # Workaround plotting piecewise (#8577):
    #   workaround for `lambdify` in `.experimental_lambdify` fails
    #   to return numerical values in some cases. Lower-level fix
    #   in `lambdify` is possible.
    vector_a = (x - y).astype(np.float64)
    vector_b = (z - y).astype(np.float64)
    dot_product = np.dot(vector_a, vector_b)
    vector_a_norm = np.linalg.norm(vector_a)
    vector_b_norm = np.linalg.norm(vector_b)
    cos_theta = dot_product / (vector_a_norm * vector_b_norm)
    return abs(cos_theta + 1) < eps


def _matplotlib_list(interval_list):
    """
    Returns lists for matplotlib ``fill`` command from a list of bounding
    rectangular intervals
    """
    xlist = []
    ylist = []
    if len(interval_list):
        for intervals in interval_list:
            intervalx = intervals[0]
            intervaly = intervals[1]
            xlist.extend([intervalx.start, intervalx.start,
                          intervalx.end, intervalx.end, None])
            ylist.extend([intervaly.start, intervaly.end,
                          intervaly.end, intervaly.start, None])
    else:
        #XXX Ugly hack. Matplotlib does not accept empty lists for ``fill``
        xlist.extend((None, None, None, None))
        ylist.extend((None, None, None, None))
    return xlist, ylist


####New API for plotting module ####

# TODO: Add color arrays for plots.
# TODO: Add more plotting options for 3d plots.
# TODO: Adaptive sampling for 3D plots.

def plot(*args, show=True, **kwargs):
    """Plots a function of a single variable as a curve.

    Parameters
    ==========

    args :
        The first argument is the expression representing the function
        of single variable to be plotted.

        The last argument is a 3-tuple denoting the range of the free
        variable. e.g. ``(x, 0, 5)``

        Typical usage examples are in the following:

        - Plotting a single expression with a single range.
            ``plot(expr, range, **kwargs)``
        - Plotting a single expression with the default range (-10, 10).
            ``plot(expr, **kwargs)``
        - Plotting multiple expressions with a single range.
            ``plot(expr1, expr2, ..., range, **kwargs)``
        - Plotting multiple expressions with multiple ranges.
            ``plot((expr1, range1), (expr2, range2), ..., **kwargs)``

        It is best practice to specify range explicitly because default
        range may change in the future if a more advanced default range
        detection algorithm is implemented.

    show : bool, optional
        The default value is set to ``True``. Set show to ``False`` and
        the function will not display the plot. The returned instance of
        the ``Plot`` class can then be used to save or display the plot
        by calling the ``save()`` and ``show()`` methods respectively.

    line_color : string, or float, or function, optional
        Specifies the color for the plot.
        See ``Plot`` to see how to set color for the plots.
        Note that by setting ``line_color``, it would be applied simultaneously
        to all the series.

    title : str, optional
        Title of the plot. It is set to the latex representation of
        the expression, if the plot has only one expression.

    label : str, optional
        The label of the expression in the plot. It will be used when
        called with ``legend``. Default is the name of the expression.
        e.g. ``sin(x)``

    xlabel : str or expression, optional
        Label for the x-axis.

    ylabel : str or expression, optional
        Label for the y-axis.

    xscale : 'linear' or 'log', optional
        Sets the scaling of the x-axis.

    yscale : 'linear' or 'log', optional
        Sets the scaling of the y-axis.

    axis_center : (float, float), optional
        Tuple of two floats denoting the coordinates of the center or
        {'center', 'auto'}

    xlim : (float, float), optional
        Denotes the x-axis limits, ``(min, max)```.

    ylim : (float, float), optional
        Denotes the y-axis limits, ``(min, max)```.

    annotations : list, optional
        A list of dictionaries specifying the type of annotation
        required. The keys in the dictionary should be equivalent
        to the arguments of the :external:mod:`matplotlib`'s
        :external:meth:`~matplotlib.axes.Axes.annotate` method.

    markers : list, optional
        A list of dictionaries specifying the type the markers required.
        The keys in the dictionary should be equivalent to the arguments
        of the :external:mod:`matplotlib`'s :external:func:`~matplotlib.pyplot.plot()` function
        along with the marker related keyworded arguments.

    rectangles : list, optional
        A list of dictionaries specifying the dimensions of the
        rectangles to be plotted. The keys in the dictionary should be
        equivalent to the arguments of the :external:mod:`matplotlib`'s
        :external:class:`~matplotlib.patches.Rectangle` class.

    fill : dict, optional
        A dictionary specifying the type of color filling required in
        the plot. The keys in the dictionary should be equivalent to the
        arguments of the :external:mod:`matplotlib`'s
        :external:meth:`~matplotlib.axes.Axes.fill_between` method.

    adaptive : bool, optional
        The default value is set to ``True``. Set adaptive to ``False``
        and specify ``n`` if uniform sampling is required.

        The plotting uses an adaptive algorithm which samples
        recursively to accurately plot. The adaptive algorithm uses a
        random point near the midpoint of two points that has to be
        further sampled. Hence the same plots can appear slightly
        different.

    depth : int, optional
        Recursion depth of the adaptive algorithm. A depth of value
        `n` samples a maximum of `2^{n}` points.

        If the ``adaptive`` flag is set to ``False``, this will be
        ignored.

    n : int, optional
        Used when the ``adaptive`` is set to ``False``. The function
        is uniformly sampled at ``n`` number of points. If the ``adaptive``
        flag is set to ``True``, this will be ignored.
        This keyword argument replaces ``nb_of_points``, which should be
        considered deprecated.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of
        the overall figure. The default value is set to ``None``, meaning
        the size will be set by the default backend.

    Examples
    ========

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> from sympy import symbols
       >>> from sympy.plotting import plot
       >>> x = symbols('x')

    Single Plot

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot(x**2, (x, -5, 5))
       Plot object containing:
       [0]: cartesian line: x**2 for x over (-5.0, 5.0)

    Multiple plots with single range.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot(x, x**2, x**3, (x, -5, 5))
       Plot object containing:
       [0]: cartesian line: x for x over (-5.0, 5.0)
       [1]: cartesian line: x**2 for x over (-5.0, 5.0)
       [2]: cartesian line: x**3 for x over (-5.0, 5.0)

    Multiple plots with different ranges.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot((x**2, (x, -6, 6)), (x, (x, -5, 5)))
       Plot object containing:
       [0]: cartesian line: x**2 for x over (-6.0, 6.0)
       [1]: cartesian line: x for x over (-5.0, 5.0)

    No adaptive sampling.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot(x**2, adaptive=False, n=400)
       Plot object containing:
       [0]: cartesian line: x**2 for x over (-10.0, 10.0)

    See Also
    ========

    Plot, LineOver1DRangeSeries

    """
    kwargs = _set_discretization_points(kwargs, LineOver1DRangeSeries)
    args = list(map(sympify, args))
    free = set()
    for a in args:
        if isinstance(a, Expr):
            free |= a.free_symbols
            if len(free) > 1:
                raise ValueError(
                    'The same variable should be used in all '
                    'univariate expressions being plotted.')
    x = free.pop() if free else Symbol('x')
    kwargs.setdefault('xlabel', x)
    kwargs.setdefault('ylabel', Function('f')(x))
    series = []
    plot_expr = check_arguments(args, 1, 1)
    series = [LineOver1DRangeSeries(*arg, **kwargs) for arg in plot_expr]

    plots = plot_factory(*series, **kwargs)
    if show:
        plots.show()
    return plots


def plot_parametric(*args, show=True, **kwargs):
    """
    Plots a 2D parametric curve.

    Parameters
    ==========

    args
        Common specifications are:

        - Plotting a single parametric curve with a range
            ``plot_parametric((expr_x, expr_y), range)``
        - Plotting multiple parametric curves with the same range
            ``plot_parametric((expr_x, expr_y), ..., range)``
        - Plotting multiple parametric curves with different ranges
            ``plot_parametric((expr_x, expr_y, range), ...)``

        ``expr_x`` is the expression representing $x$ component of the
        parametric function.

        ``expr_y`` is the expression representing $y$ component of the
        parametric function.

        ``range`` is a 3-tuple denoting the parameter symbol, start and
        stop. For example, ``(u, 0, 5)``.

        If the range is not specified, then a default range of (-10, 10)
        is used.

        However, if the arguments are specified as
        ``(expr_x, expr_y, range), ...``, you must specify the ranges
        for each expressions manually.

        Default range may change in the future if a more advanced
        algorithm is implemented.

    adaptive : bool, optional
        Specifies whether to use the adaptive sampling or not.

        The default value is set to ``True``. Set adaptive to ``False``
        and specify ``n`` if uniform sampling is required.

    depth :  int, optional
        The recursion depth of the adaptive algorithm. A depth of
        value $n$ samples a maximum of $2^n$ points.

    n : int, optional
        Used when the ``adaptive`` flag is set to ``False``. Specifies the
        number of the points used for the uniform sampling.
        This keyword argument replaces ``nb_of_points``, which should be
        considered deprecated.

    line_color : string, or float, or function, optional
        Specifies the color for the plot.
        See ``Plot`` to see how to set color for the plots.
        Note that by setting ``line_color``, it would be applied simultaneously
        to all the series.

    label : str, optional
        The label of the expression in the plot. It will be used when
        called with ``legend``. Default is the name of the expression.
        e.g. ``sin(x)``

    xlabel : str, optional
        Label for the x-axis.

    ylabel : str, optional
        Label for the y-axis.

    xscale : 'linear' or 'log', optional
        Sets the scaling of the x-axis.

    yscale : 'linear' or 'log', optional
        Sets the scaling of the y-axis.

    axis_center : (float, float), optional
        Tuple of two floats denoting the coordinates of the center or
        {'center', 'auto'}

    xlim : (float, float), optional
        Denotes the x-axis limits, ``(min, max)```.

    ylim : (float, float), optional
        Denotes the y-axis limits, ``(min, max)```.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of
        the overall figure. The default value is set to ``None``, meaning
        the size will be set by the default backend.

    Examples
    ========

    .. plot::
       :context: reset
       :format: doctest
       :include-source: True

       >>> from sympy import plot_parametric, symbols, cos, sin
       >>> u = symbols('u')

    A parametric plot with a single expression:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot_parametric((cos(u), sin(u)), (u, -5, 5))
       Plot object containing:
       [0]: parametric cartesian line: (cos(u), sin(u)) for u over (-5.0, 5.0)

    A parametric plot with multiple expressions with the same range:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot_parametric((cos(u), sin(u)), (u, cos(u)), (u, -10, 10))
       Plot object containing:
       [0]: parametric cartesian line: (cos(u), sin(u)) for u over (-10.0, 10.0)
       [1]: parametric cartesian line: (u, cos(u)) for u over (-10.0, 10.0)

    A parametric plot with multiple expressions with different ranges
    for each curve:

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot_parametric((cos(u), sin(u), (u, -5, 5)),
       ...     (cos(u), u, (u, -5, 5)))
       Plot object containing:
       [0]: parametric cartesian line: (cos(u), sin(u)) for u over (-5.0, 5.0)
       [1]: parametric cartesian line: (cos(u), u) for u over (-5.0, 5.0)

    Notes
    =====

    The plotting uses an adaptive algorithm which samples recursively to
    accurately plot the curve. The adaptive algorithm uses a random point
    near the midpoint of two points that has to be further sampled.
    Hence, repeating the same plot command can give slightly different
    results because of the random sampling.

    If there are multiple plots, then the same optional arguments are
    applied to all the plots drawn in the same canvas. If you want to
    set these options separately, you can index the returned ``Plot``
    object and set it.

    For example, when you specify ``line_color`` once, it would be
    applied simultaneously to both series.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> from sympy import pi
        >>> expr1 = (u, cos(2*pi*u)/2 + 1/2)
        >>> expr2 = (u, sin(2*pi*u)/2 + 1/2)
        >>> p = plot_parametric(expr1, expr2, (u, 0, 1), line_color='blue')

    If you want to specify the line color for the specific series, you
    should index each item and apply the property manually.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

        >>> p[0].line_color = 'red'
        >>> p.show()

    See Also
    ========

    Plot, Parametric2DLineSeries
    """
    kwargs = _set_discretization_points(kwargs, Parametric2DLineSeries)
    args = list(map(sympify, args))
    series = []
    plot_expr = check_arguments(args, 2, 1)
    series = [Parametric2DLineSeries(*arg, **kwargs) for arg in plot_expr]
    plots = plot_factory(*series, **kwargs)
    if show:
        plots.show()
    return plots


def plot3d_parametric_line(*args, show=True, **kwargs):
    """
    Plots a 3D parametric line plot.

    Usage
    =====

    Single plot:

    ``plot3d_parametric_line(expr_x, expr_y, expr_z, range, **kwargs)``

    If the range is not specified, then a default range of (-10, 10) is used.

    Multiple plots.

    ``plot3d_parametric_line((expr_x, expr_y, expr_z, range), ..., **kwargs)``

    Ranges have to be specified for every expression.

    Default range may change in the future if a more advanced default range
    detection algorithm is implemented.

    Arguments
    =========

    expr_x : Expression representing the function along x.

    expr_y : Expression representing the function along y.

    expr_z : Expression representing the function along z.

    range : (:class:`~.Symbol`, float, float)
        A 3-tuple denoting the range of the parameter variable, e.g., (u, 0, 5).

    Keyword Arguments
    =================

    Arguments for ``Parametric3DLineSeries`` class.

    n : int
        The range is uniformly sampled at ``n`` number of points.
        This keyword argument replaces ``nb_of_points``, which should be
        considered deprecated.

    Aesthetics:

    line_color : string, or float, or function, optional
        Specifies the color for the plot.
        See ``Plot`` to see how to set color for the plots.
        Note that by setting ``line_color``, it would be applied simultaneously
        to all the series.

    label : str
        The label to the plot. It will be used when called with ``legend=True``
        to denote the function with the given label in the plot.

    If there are multiple plots, then the same series arguments are applied to
    all the plots. If you want to set these options separately, you can index
    the returned ``Plot`` object and set it.

    Arguments for ``Plot`` class.

    title : str
        Title of the plot.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of
        the overall figure. The default value is set to ``None``, meaning
        the size will be set by the default backend.

    Examples
    ========

    .. plot::
       :context: reset
       :format: doctest
       :include-source: True

       >>> from sympy import symbols, cos, sin
       >>> from sympy.plotting import plot3d_parametric_line
       >>> u = symbols('u')

    Single plot.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d_parametric_line(cos(u), sin(u), u, (u, -5, 5))
       Plot object containing:
       [0]: 3D parametric cartesian line: (cos(u), sin(u), u) for u over (-5.0, 5.0)


    Multiple plots.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d_parametric_line((cos(u), sin(u), u, (u, -5, 5)),
       ...     (sin(u), u**2, u, (u, -5, 5)))
       Plot object containing:
       [0]: 3D parametric cartesian line: (cos(u), sin(u), u) for u over (-5.0, 5.0)
       [1]: 3D parametric cartesian line: (sin(u), u**2, u) for u over (-5.0, 5.0)


    See Also
    ========

    Plot, Parametric3DLineSeries

    """
    kwargs = _set_discretization_points(kwargs, Parametric3DLineSeries)
    args = list(map(sympify, args))
    series = []
    plot_expr = check_arguments(args, 3, 1)
    series = [Parametric3DLineSeries(*arg, **kwargs) for arg in plot_expr]
    kwargs.setdefault("xlabel", "x")
    kwargs.setdefault("ylabel", "y")
    kwargs.setdefault("zlabel", "z")
    plots = plot_factory(*series, **kwargs)
    if show:
        plots.show()
    return plots


def plot3d(*args, show=True, **kwargs):
    """
    Plots a 3D surface plot.

    Usage
    =====

    Single plot

    ``plot3d(expr, range_x, range_y, **kwargs)``

    If the ranges are not specified, then a default range of (-10, 10) is used.

    Multiple plot with the same range.

    ``plot3d(expr1, expr2, range_x, range_y, **kwargs)``

    If the ranges are not specified, then a default range of (-10, 10) is used.

    Multiple plots with different ranges.

    ``plot3d((expr1, range_x, range_y), (expr2, range_x, range_y), ..., **kwargs)``

    Ranges have to be specified for every expression.

    Default range may change in the future if a more advanced default range
    detection algorithm is implemented.

    Arguments
    =========

    expr : Expression representing the function along x.

    range_x : (:class:`~.Symbol`, float, float)
        A 3-tuple denoting the range of the x variable, e.g. (x, 0, 5).

    range_y : (:class:`~.Symbol`, float, float)
        A 3-tuple denoting the range of the y variable, e.g. (y, 0, 5).

    Keyword Arguments
    =================

    Arguments for ``SurfaceOver2DRangeSeries`` class:

    n1 : int
        The x range is sampled uniformly at ``n1`` of points.
        This keyword argument replaces ``nb_of_points_x``, which should be
        considered deprecated.

    n2 : int
        The y range is sampled uniformly at ``n2`` of points.
        This keyword argument replaces ``nb_of_points_y``, which should be
        considered deprecated.

    Aesthetics:

    surface_color : Function which returns a float
        Specifies the color for the surface of the plot.
        See :class:`~.Plot` for more details.

    If there are multiple plots, then the same series arguments are applied to
    all the plots. If you want to set these options separately, you can index
    the returned ``Plot`` object and set it.

    Arguments for ``Plot`` class:

    title : str
        Title of the plot.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of the
        overall figure. The default value is set to ``None``, meaning the size will
        be set by the default backend.

    Examples
    ========

    .. plot::
       :context: reset
       :format: doctest
       :include-source: True

       >>> from sympy import symbols
       >>> from sympy.plotting import plot3d
       >>> x, y = symbols('x y')

    Single plot

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d(x*y, (x, -5, 5), (y, -5, 5))
       Plot object containing:
       [0]: cartesian surface: x*y for x over (-5.0, 5.0) and y over (-5.0, 5.0)


    Multiple plots with same range

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d(x*y, -x*y, (x, -5, 5), (y, -5, 5))
       Plot object containing:
       [0]: cartesian surface: x*y for x over (-5.0, 5.0) and y over (-5.0, 5.0)
       [1]: cartesian surface: -x*y for x over (-5.0, 5.0) and y over (-5.0, 5.0)


    Multiple plots with different ranges.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d((x**2 + y**2, (x, -5, 5), (y, -5, 5)),
       ...     (x*y, (x, -3, 3), (y, -3, 3)))
       Plot object containing:
       [0]: cartesian surface: x**2 + y**2 for x over (-5.0, 5.0) and y over (-5.0, 5.0)
       [1]: cartesian surface: x*y for x over (-3.0, 3.0) and y over (-3.0, 3.0)


    See Also
    ========

    Plot, SurfaceOver2DRangeSeries

    """

    kwargs = _set_discretization_points(kwargs, SurfaceOver2DRangeSeries)
    args = list(map(sympify, args))
    series = []
    plot_expr = check_arguments(args, 1, 2)
    series = [SurfaceOver2DRangeSeries(*arg, **kwargs) for arg in plot_expr]
    kwargs.setdefault("xlabel", series[0].var_x)
    kwargs.setdefault("ylabel", series[0].var_y)
    kwargs.setdefault("zlabel", Function('f')(series[0].var_x, series[0].var_y))
    plots = plot_factory(*series, **kwargs)
    if show:
        plots.show()
    return plots


def plot3d_parametric_surface(*args, show=True, **kwargs):
    """
    Plots a 3D parametric surface plot.

    Explanation
    ===========

    Single plot.

    ``plot3d_parametric_surface(expr_x, expr_y, expr_z, range_u, range_v, **kwargs)``

    If the ranges is not specified, then a default range of (-10, 10) is used.

    Multiple plots.

    ``plot3d_parametric_surface((expr_x, expr_y, expr_z, range_u, range_v), ..., **kwargs)``

    Ranges have to be specified for every expression.

    Default range may change in the future if a more advanced default range
    detection algorithm is implemented.

    Arguments
    =========

    expr_x : Expression representing the function along ``x``.

    expr_y : Expression representing the function along ``y``.

    expr_z : Expression representing the function along ``z``.

    range_u : (:class:`~.Symbol`, float, float)
        A 3-tuple denoting the range of the u variable, e.g. (u, 0, 5).

    range_v : (:class:`~.Symbol`, float, float)
        A 3-tuple denoting the range of the v variable, e.g. (v, 0, 5).

    Keyword Arguments
    =================

    Arguments for ``ParametricSurfaceSeries`` class:

    n1 : int
        The ``u`` range is sampled uniformly at ``n1`` of points.
        This keyword argument replaces ``nb_of_points_u``, which should be
        considered deprecated.

    n2 : int
        The ``v`` range is sampled uniformly at ``n2`` of points.
        This keyword argument replaces ``nb_of_points_v``, which should be
        considered deprecated.

    Aesthetics:

    surface_color : Function which returns a float
        Specifies the color for the surface of the plot. See
        :class:`~Plot` for more details.

    If there are multiple plots, then the same series arguments are applied for
    all the plots. If you want to set these options separately, you can index
    the returned ``Plot`` object and set it.


    Arguments for ``Plot`` class:

    title : str
        Title of the plot.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of the
        overall figure. The default value is set to ``None``, meaning the size will
        be set by the default backend.

    Examples
    ========

    .. plot::
       :context: reset
       :format: doctest
       :include-source: True

       >>> from sympy import symbols, cos, sin
       >>> from sympy.plotting import plot3d_parametric_surface
       >>> u, v = symbols('u v')

    Single plot.

    .. plot::
       :context: close-figs
       :format: doctest
       :include-source: True

       >>> plot3d_parametric_surface(cos(u + v), sin(u - v), u - v,
       ...     (u, -5, 5), (v, -5, 5))
       Plot object containing:
       [0]: parametric cartesian surface: (cos(u + v), sin(u - v), u - v) for u over (-5.0, 5.0) and v over (-5.0, 5.0)


    See Also
    ========

    Plot, ParametricSurfaceSeries

    """

    kwargs = _set_discretization_points(kwargs, ParametricSurfaceSeries)
    args = list(map(sympify, args))
    series = []
    plot_expr = check_arguments(args, 3, 2)
    series = [ParametricSurfaceSeries(*arg, **kwargs) for arg in plot_expr]
    kwargs.setdefault("xlabel", "x")
    kwargs.setdefault("ylabel", "y")
    kwargs.setdefault("zlabel", "z")
    plots = plot_factory(*series, **kwargs)
    if show:
        plots.show()
    return plots

def plot_contour(*args, show=True, **kwargs):
    """
    Draws contour plot of a function

    Usage
    =====

    Single plot

    ``plot_contour(expr, range_x, range_y, **kwargs)``

    If the ranges are not specified, then a default range of (-10, 10) is used.

    Multiple plot with the same range.

    ``plot_contour(expr1, expr2, range_x, range_y, **kwargs)``

    If the ranges are not specified, then a default range of (-10, 10) is used.

    Multiple plots with different ranges.

    ``plot_contour((expr1, range_x, range_y), (expr2, range_x, range_y), ..., **kwargs)``

    Ranges have to be specified for every expression.

    Default range may change in the future if a more advanced default range
    detection algorithm is implemented.

    Arguments
    =========

    expr : Expression representing the function along x.

    range_x : (:class:`Symbol`, float, float)
        A 3-tuple denoting the range of the x variable, e.g. (x, 0, 5).

    range_y : (:class:`Symbol`, float, float)
        A 3-tuple denoting the range of the y variable, e.g. (y, 0, 5).

    Keyword Arguments
    =================

    Arguments for ``ContourSeries`` class:

    n1 : int
        The x range is sampled uniformly at ``n1`` of points.
        This keyword argument replaces ``nb_of_points_x``, which should be
        considered deprecated.

    n2 : int
        The y range is sampled uniformly at ``n2`` of points.
        This keyword argument replaces ``nb_of_points_y``, which should be
        considered deprecated.

    Aesthetics:

    surface_color : Function which returns a float
        Specifies the color for the surface of the plot. See
        :class:`sympy.plotting.Plot` for more details.

    If there are multiple plots, then the same series arguments are applied to
    all the plots. If you want to set these options separately, you can index
    the returned ``Plot`` object and set it.

    Arguments for ``Plot`` class:

    title : str
        Title of the plot.

    size : (float, float), optional
        A tuple in the form (width, height) in inches to specify the size of
        the overall figure. The default value is set to ``None``, meaning
        the size will be set by the default backend.

    See Also
    ========

    Plot, ContourSeries

    """

    kwargs = _set_discretization_points(kwargs, ContourSeries)
    args = list(map(sympify, args))
    plot_expr = check_arguments(args, 1, 2)
    series = [ContourSeries(*arg) for arg in plot_expr]
    plot_contours = plot_factory(*series, **kwargs)
    if len(plot_expr[0].free_symbols) > 2:
        raise ValueError('Contour Plot cannot Plot for more than two variables.')
    if show:
        plot_contours.show()
    return plot_contours

def check_arguments(args, expr_len, nb_of_free_symbols):
    """
    Checks the arguments and converts into tuples of the
    form (exprs, ranges).

    Examples
    ========

    .. plot::
       :context: reset
       :format: doctest
       :include-source: True

       >>> from sympy import cos, sin, symbols
       >>> from sympy.plotting.plot import check_arguments
       >>> x = symbols('x')
       >>> check_arguments([cos(x), sin(x)], 2, 1)
           [(cos(x), sin(x), (x, -10, 10))]

       >>> check_arguments([x, x**2], 1, 1)
           [(x, (x, -10, 10)), (x**2, (x, -10, 10))]
    """
    if not args:
        return []
    if expr_len > 1 and isinstance(args[0], Expr):
        # Multiple expressions same range.
        # The arguments are tuples when the expression length is
        # greater than 1.
        if len(args) < expr_len:
            raise ValueError("len(args) should not be less than expr_len")
        for i in range(len(args)):
            if isinstance(args[i], Tuple):
                break
        else:
            i = len(args) + 1

        exprs = Tuple(*args[:i])
        free_symbols = list(set().union(*[e.free_symbols for e in exprs]))
        if len(args) == expr_len + nb_of_free_symbols:
            #Ranges given
            plots = [exprs + Tuple(*args[expr_len:])]
        else:
            default_range = Tuple(-10, 10)
            ranges = []
            for symbol in free_symbols:
                ranges.append(Tuple(symbol) + default_range)

            for i in range(len(free_symbols) - nb_of_free_symbols):
                ranges.append(Tuple(Dummy()) + default_range)
            plots = [exprs + Tuple(*ranges)]
        return plots

    if isinstance(args[0], Expr) or (isinstance(args[0], Tuple) and
                                     len(args[0]) == expr_len and
                                     expr_len != 3):
        # Cannot handle expressions with number of expression = 3. It is
        # not possible to differentiate between expressions and ranges.
        #Series of plots with same range
        for i in range(len(args)):
            if isinstance(args[i], Tuple) and len(args[i]) != expr_len:
                break
            if not isinstance(args[i], Tuple):
                args[i] = Tuple(args[i])
        else:
            i = len(args) + 1

        exprs = args[:i]
        assert all(isinstance(e, Expr) for expr in exprs for e in expr)
        free_symbols = list(set().union(*[e.free_symbols for expr in exprs
                                        for e in expr]))

        if len(free_symbols) > nb_of_free_symbols:
            raise ValueError("The number of free_symbols in the expression "
                             "is greater than %d" % nb_of_free_symbols)
        if len(args) == i + nb_of_free_symbols and isinstance(args[i], Tuple):
            ranges = Tuple(*list(args[
                           i:i + nb_of_free_symbols]))
            plots = [expr + ranges for expr in exprs]
            return plots
        else:
            # Use default ranges.
            default_range = Tuple(-10, 10)
            ranges = []
            for symbol in free_symbols:
                ranges.append(Tuple(symbol) + default_range)

            for i in range(nb_of_free_symbols - len(free_symbols)):
                ranges.append(Tuple(Dummy()) + default_range)
            ranges = Tuple(*ranges)
            plots = [expr + ranges for expr in exprs]
            return plots

    elif isinstance(args[0], Tuple) and len(args[0]) == expr_len + nb_of_free_symbols:
        # Multiple plots with different ranges.
        for arg in args:
            for i in range(expr_len):
                if not isinstance(arg[i], Expr):
                    raise ValueError("Expected an expression, given %s" %
                                     str(arg[i]))
            for i in range(nb_of_free_symbols):
                if not len(arg[i + expr_len]) == 3:
                    raise ValueError("The ranges should be a tuple of "
                                     "length 3, got %s" % str(arg[i + expr_len]))
        return args
