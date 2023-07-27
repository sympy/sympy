### The base class for all series
from collections.abc import Callable
from sympy.core.containers import Tuple
from sympy.core.function import arity
from sympy.core.relational import (Equality, GreaterThan, StrictGreaterThan,
    LessThan, StrictLessThan)
from sympy.core.sympify import sympify
from sympy.external import import_module
from sympy.utilities.exceptions import sympy_deprecation_warning
from .experimental_lambdify import (vectorized_lambdify, lambdify,
    experimental_lambdify)
from .intervalmath import interval
import warnings


__all__ = [
    "List2DSeries", "LineOver1DRangeSeries", "Parametric2DLineSeries",
    "Parametric3DLineSeries", "SurfaceOver2DRangeSeries", "ContourSeries",
    "ParametricSurfaceSeries", "ImplicitSeries", "GenericDataSeries",
    # for back-compatibility
    "BaseSeries", "Line2DBaseSeries", "Line3DBaseSeries", "SurfaceBaseSeries",
    "centers_of_segments", "centers_of_faces", "flat"
]


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


class ImplicitSeries(BaseSeries):
    """ Representation for Implicit plot """
    is_implicit = True

    def __init__(self, expr, var_start_end_x, var_start_end_y,
            has_equality, use_interval_math, depth, nb_of_points,
            line_color):
        super().__init__()
        self.expr = sympify(expr)
        self.label = self.expr
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])
        self.get_points = self.get_raster
        self.has_equality = has_equality  # If the expression has equality, i.e.
                                         #Eq, Greaterthan, LessThan.
        self.nb_of_points = nb_of_points
        self.use_interval_math = use_interval_math
        self.depth = 4 + depth
        self.line_color = line_color

    def __str__(self):
        return ('Implicit equation: %s for '
                '%s over %s and %s over %s') % (
                    str(self.expr),
                    str(self.var_x),
                    str((self.start_x, self.end_x)),
                    str(self.var_y),
                    str((self.start_y, self.end_y)))

    def get_data(self):
        """Returns numerical data.

        Returns
        =======

        If the series is evaluated with the `adaptive=True` it returns:

        interval_list : list
            List of bounding rectangular intervals to be postprocessed and
            eventually used with Matplotlib's ``fill`` command.
        dummy : str
            A string containing ``"fill"``.

        Otherwise, it returns 2D numpy arrays to be used with Matplotlib's
        ``contour`` or ``contourf`` commands:

        x_array : np.ndarray
        y_array : np.ndarray
        z_array : np.ndarray
        plot_type : str
            A string specifying which plot command to use, ``"contour"``
            or ``"contourf"``.
        """
        func = experimental_lambdify((self.var_x, self.var_y), self.expr,
                                    use_interval=True)
        xinterval = interval(self.start_x, self.end_x)
        yinterval = interval(self.start_y, self.end_y)
        try:
            func(xinterval, yinterval)
        except AttributeError:
            # XXX: AttributeError("'list' object has no attribute 'is_real'")
            # That needs fixing somehow - we shouldn't be catching
            # AttributeError here.
            if self.use_interval_math:
                warnings.warn("Adaptive meshing could not be applied to the"
                            " expression. Using uniform meshing.", stacklevel=7)
            self.use_interval_math = False

        if self.use_interval_math:
            return self._get_raster_interval(func)
        else:
            return self._get_meshes_grid()

    def get_raster(self):
        """Returns numerical data.

        This function is available for back-compatibility purposes. Consider
        using ``get_data()`` instead.
        """
        return self.get_data()

    def _get_raster_interval(self, func):
        """ Uses interval math to adaptively mesh and obtain the plot"""
        k = self.depth
        interval_list = []
        #Create initial 32 divisions
        np = import_module('numpy')
        xsample = np.linspace(self.start_x, self.end_x, 33)
        ysample = np.linspace(self.start_y, self.end_y, 33)

        #Add a small jitter so that there are no false positives for equality.
        # Ex: y==x becomes True for x interval(1, 2) and y interval(1, 2)
        #which will draw a rectangle.
        jitterx = (np.random.rand(
            len(xsample)) * 2 - 1) * (self.end_x - self.start_x) / 2**20
        jittery = (np.random.rand(
            len(ysample)) * 2 - 1) * (self.end_y - self.start_y) / 2**20
        xsample += jitterx
        ysample += jittery

        xinter = [interval(x1, x2) for x1, x2 in zip(xsample[:-1],
                           xsample[1:])]
        yinter = [interval(y1, y2) for y1, y2 in zip(ysample[:-1],
                           ysample[1:])]
        interval_list = [[x, y] for x in xinter for y in yinter]
        plot_list = []

        #recursive call refinepixels which subdivides the intervals which are
        #neither True nor False according to the expression.
        def refine_pixels(interval_list):
            """ Evaluates the intervals and subdivides the interval if the
            expression is partially satisfied."""
            temp_interval_list = []
            plot_list = []
            for intervals in interval_list:

                #Convert the array indices to x and y values
                intervalx = intervals[0]
                intervaly = intervals[1]
                func_eval = func(intervalx, intervaly)
                #The expression is valid in the interval. Change the contour
                #array values to 1.
                if func_eval[1] is False or func_eval[0] is False:
                    pass
                elif func_eval == (True, True):
                    plot_list.append([intervalx, intervaly])
                elif func_eval[1] is None or func_eval[0] is None:
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
        if self.has_equality:
            for intervals in interval_list:
                intervalx = intervals[0]
                intervaly = intervals[1]
                func_eval = func(intervalx, intervaly)
                if func_eval[1] and func_eval[0] is not False:
                    plot_list.append([intervalx, intervaly])
        return plot_list, 'fill'

    def _get_meshes_grid(self):
        """Generates the mesh for generating a contour.

        In the case of equality, ``contour`` function of matplotlib can
        be used. In other cases, matplotlib's ``contourf`` is used.
        """
        equal = False
        if isinstance(self.expr, Equality):
            expr = self.expr.lhs - self.expr.rhs
            equal = True

        elif isinstance(self.expr, (GreaterThan, StrictGreaterThan)):
            expr = self.expr.lhs - self.expr.rhs

        elif isinstance(self.expr, (LessThan, StrictLessThan)):
            expr = self.expr.rhs - self.expr.lhs
        else:
            raise NotImplementedError("The expression is not supported for "
                                    "plotting in uniform meshed plot.")
        np = import_module('numpy')
        xarray = np.linspace(self.start_x, self.end_x, self.nb_of_points)
        yarray = np.linspace(self.start_y, self.end_y, self.nb_of_points)
        x_grid, y_grid = np.meshgrid(xarray, yarray)

        func = vectorized_lambdify((self.var_x, self.var_y), expr)
        z_grid = func(x_grid, y_grid)
        z_grid[np.ma.where(z_grid < 0)] = -1
        z_grid[np.ma.where(z_grid > 0)] = 1
        if equal:
            return xarray, yarray, z_grid, 'contour'
        else:
            return xarray, yarray, z_grid, 'contourf'


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
