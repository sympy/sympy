from sympy import (sympify, Dummy, Matrix, Basic, Expr, solve, diff, factorial,
        Function)
from sympy.core.compatibility import reduce

# TODO order the imports and make them explicit
# TODO issue 2070: all the stuff about .args and rebuilding
# TODO maybe a common class Field makes sense

class Manifold(Basic):
    """Object representing a mathematical manifold.

    The only role that this object plays is to keep a list of all patches
    defined on the manifold. It does not provide any means to study the
    topological characteristics of the manifold that it represents.
    """
    def __init__(self, name, dim):
        super(Manifold, self).__init__()
        self.name = name
        # TODO What is the purpose of this name, besides printing?
        self.dim = dim
        self.patches = []
        # The patches list is necessary if a Patch instance needs to enumerate
        # other Patch instance on the same manifold.


class Patch(Basic):
    """Object representing a patch on a manifold.

    On a manifold one can have many patches that do not always include the
    whole manifold. On these patches coordinate charts can be defined that
    permit the parametrization of any point on the patch in terms of a tuple
    of real numbers (the coordinates).

    This object serves as a container/parent for all coordinate system charts
    that can be defined on the patch it represents.

    Examples:
    =========

    Define a Manifold and a Patch on that Manifold:
    >>> from sympy.differential_geometry import Manifold, Patch
    >>> m = Manifold('M', 3)
    >>> p = Patch('P', m)
    >>> p in m.patches
    True

    """
    # Contains a reference to the parent manifold in order to be able to access
    # other patches.
    def __init__(self, name, manifold):
        super(Patch, self).__init__()
        self.name = name
        # TODO What is the purpose of this name, besides printing?
        self.manifold = manifold
        self.manifold.patches.append(self)
        self.coord_systems = []
        # The list of coordinate systems is necessary for an instance of
        # CoordSystem to enumerate other coord systems on the patch.


class CoordSystem(Basic):
    """Contains all coordinate transformation logic.

    Examples:
    =========

    Define a Manifold and a Patch, and then define two coord systems on that
    patch:
    >>> from sympy import symbols, sin, cos, pi
    >>> from sympy.differential_geometry import Manifold, Patch, CoordSystem
    >>> x, y, r, theta = symbols('x, y, r, theta')
    >>> m = Manifold('M', 2)
    >>> p = Patch('P', m)
    >>> rect = CoordSystem('rect', p)
    >>> polar = CoordSystem('polar', p)
    >>> rect in p.coord_systems
    True

    Connect the coordinate systems. An inverse transformation is automatically
    found by `solve` when possible:
    >>> polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])
    >>> polar.coord_transform_to(rect, [0, 2])
    [0]
    [0]
    >>> polar.coord_transform_to(rect, [2, pi/2])
    [0]
    [2]
    >>> rect.coord_transform_to(polar, [1, 1])
    [sqrt(2)]
    [   pi/4]

    Define a point using coordinates in one of the coordinate systems:
    >>> p = polar.point([1, 3*pi/4])
    >>> rect.point_to_coords(p)
    [-sqrt(2)/2]
    [ sqrt(2)/2]

    Define a basis scalar field (i.e. a coordinate function), that takes a
    point and returns its coordinates. It is an instance of `ScalarField`.
    >>> rect.coord_function(0)(p)
    -sqrt(2)/2
    >>> rect.coord_function(1)(p)
    sqrt(2)/2

    Define a basis vector field (i.e. a unit vector field along the coordinate
    line). Vectors are also differential operators on scalar fields. It is an
    instance of `VectorField`.
    >>> v_x = rect.base_vector(0)
    >>> x = rect.coord_function(0)
    >>> v_x(x)(p)
    1
    >>> v_x(v_x(x))(p)
    0

    """
    #  Contains a reference to the parent patch in order to be able to access
    # other coordinate system charts.
    def __init__(self, name, patch):
        super(CoordSystem, self).__init__()
        self.name = name
        # TODO What is the purpose of this name, besides printing?
        self.patch = patch
        self.patch.coord_systems.append(self)
        self.transforms = {}
        # All the coordinate transformation logic is in this dictionary in the
        # form of:
        #  key = other coordinate system
        #  value = tuple of  # TODO make these Lambda instances
        #          - list of `Dummy` coordinates in this coordinate system
        #          - list of expressions as a function of the Dummies giving
        #          the coordinates in another coordinate system

    def connect_to(self, to_sys, from_coords, to_exprs, inverse=True, fill_in_gaps=True):
        """Register the transformation used to switch to another coordinate system.

        Arguments:
        ==========
        - to_sys - another instance of `CoordSystem`
        - from_coords - list of symbols in terms of which `to_exprs` is given
        - to_exprs - list of the expressions of the new coordinate tuple
        - inverse - try to deduce and register the inverse transformation
        - fill_in_gaps - try to deduce other transformation that are made
        possible by composing the present transformation with other already
        registered transformation
        """
        from_coords, to_exprs = dummyfy(from_coords, to_exprs)
        self.transforms[to_sys] = Matrix(from_coords), Matrix(to_exprs)

        if inverse:
            to_sys.transforms[self] = self._inv_transf(from_coords, to_exprs)

        if fill_in_gaps:
            self._fill_gaps_in_transformations()

    @staticmethod
    def _inv_transf(from_coords, to_exprs):
        # TODO, check for results, get solve to return results in definite
        # format instead of wondering dict/tuple/whatever.
        # As it is at the moment this is an ugly hack for changing the format
        inv_from = [i.as_dummy() for i in from_coords]
        inv_to = solve([t[0]-t[1] for t in zip(inv_from, to_exprs)], list(from_coords))
        if isinstance(inv_to, dict):
            inv_to = [inv_to[fc] for fc in from_coords]
        else:
            inv_to = inv_to[0]
        return Matrix(inv_from), Matrix(inv_to)

    @staticmethod
    def _fill_gaps_in_transformations():
        # TODO
        pass

    def coord_transform_to(self, to_sys, coords):
        """Transform `coords` to coord system `to_sys`.

        See the docstring of `CoordSystem` for examples."""
        coords = Matrix(coords)
        if self != to_sys:
            transf = self.transforms[to_sys]
            coords = transf[1].subs(zip(transf[0], coords))
        return coords

    def coord_function(self, coord_index):
        """Return a `ScalarField` that takes a point and returns one of the coords.

        Takes a point and returns its coordinate in this coordinate system.

        See the docstring of `CoordSystem` for examples."""
        args = [Dummy() for i in range(self.patch.manifold.dim)]
        result = args[coord_index]
        return ScalarField(self, args, result)

    def coord_functions(self):
        """Returns a list of all coordinate functions.

        For more details see the coord_function method of this class."""
        return [self.coord_function(i) for i in range(self.patch.manifold.dim)]

    def base_vector(self, coord_index):
        """Return a basis VectorField.

        The basis vector field for this coordinate system. It is also an
        operator on scalar fields.

        See the docstring of `CoordSystem` for examples."""
        args = [Dummy() for i in range(self.patch.manifold.dim)]
        result = [0,] * self.patch.manifold.dim
        result[coord_index] = 1
        return VectorField(self, args, result)

    def point(self, coords):
        """Create a `Point` with coordinates given in this coord system.

        See the docstring of `CoordSystem` for examples."""
        return Point(self, coords)

    def point_to_coords(self, point):
        """Calculate the coordinates of a point in this coord system.

        See the docstring of `CoordSystem` for examples."""
        return point.coords(self)


class Point(Basic):
    """Point in a Manifold object.

    To define a point you must supply coordinates and a coordinate system.

    The usage of this object after its definition is independent of the
    coordinate system that was used in order to define it, however due to
    limitations in the simplification routines you can arrive at complicated
    expressions if you use inappropriate coordinate systems.

    Examples:
    =========

    Define the boilerplate Manifold, Patch and coordinate systems:
    >>> from sympy import symbols, sin, cos, pi
    >>> from sympy.differential_geometry import (
    ...        Manifold, Patch, CoordSystem, Point)
    >>> x, y, r, theta = symbols('x, y, r, theta')
    >>> m = Manifold('M', 2)
    >>> p = Patch('P', m)
    >>> rect = CoordSystem('rect', p)
    >>> polar = CoordSystem('polar', p)
    >>> polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])

    Define a point using coordinates from one of the coordinate systems:
    >>> p = Point(polar, [r, 3*pi/4])
    >>> p.coords()
    [     r]
    [3*pi/4]
    >>> p.coords(rect)
    [-sqrt(2)*r/2]
    [ sqrt(2)*r/2]

    """
    def __init__(self, coord_sys, coords):
        super(Point, self).__init__()
        self._coord_sys = coord_sys
        self._coords = Matrix(coords)

    def coords(self, to_sys=None):
        """Coordinates of the point in a given coordinate system.

        If `to_sys` is None it returns the coordinates in the system in
        which the point was defined."""
        if to_sys:
            return self._coord_sys.coord_transform_to(to_sys, self._coords)
        else:
            return self._coords


class ScalarField(Expr):
    """Scalar Field over a Manifold.

    A scalar field takes a point as an argument and returns a scalar.

    To define a scalar field you need to choose a coordinate system and define
    the scalar field in terms of that coordinate system.

    The use of the scalar field after its definition is independent of the
    coordinate system in which it was defined, however due to limitations in
    the simplification routines you may arrive at more complicated
    expression if you use unappropriate coordinate systems.

    Examples:
    =========

    Define boilerplate Manifold, Patch and coordinate systems:
    >>> from sympy import symbols, sin, cos, pi, Function
    >>> from sympy.differential_geometry import (
    ...        Manifold, Patch, CoordSystem, Point, ScalarField)
    >>> x, y, r, theta = symbols('x, y, r, theta')
    >>> m = Manifold('M', 2)
    >>> p = Patch('P', m)
    >>> rect = CoordSystem('rect', p)
    >>> polar = CoordSystem('polar', p)
    >>> polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])

    Points to be used as arguments for the filed:
    >>> pointA = polar.point([r, 0])

    Examples of fields:
    >>> field1 = ScalarField(rect, [x, y], x**2+y**2)
    >>> field1(pointA)
    r**2
    >>> (1+5*field1)(pointA)
    5*r**2 + 1

    >>> g = Function('g')
    >>> field2 = ScalarField(polar, [r, theta], g(theta-pi))
    >>> field2(pointA)
    g(-pi)

    """
    def __init__(self, coord_sys, coords, expr):
        super(ScalarField, self).__init__()
        self._coord_sys = coord_sys
        coords, (expr, ) = dummyfy(coords, (expr, ))
        self._coords = coords
        self._expr = expr
        self._args = self._coord_sys, self._coords, self._expr

    def __call__(self, point):
        # XXX see the if statement
        #this is necessary in order to construct vector fields from stuff
        # like scalar field * vector field, however it need more justification
        # It is the method used in Functional Differential Geometry
        if not isinstance(point, Point):
            return self
        coords = point.coords(self._coord_sys)
        return self._expr.subs(zip(self._coords, coords))


class VectorField(Expr):
    """Vector Field over a Manifold.

    A vector field is an operator taking a scalar field and returning a
    directional derivative (which is also a scalar field).

    To define a vector field you need to choose a coordinate system and define
    the vector field in terms of that coordinate system.

    The use of the vector field after its definition is independent of the
    coordinate system in which it was defined, however due to limitations in
    the simplification routines you may arrive at more complicated
    expression if you use unappropriate coordinate systems.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy import symbols, sin, cos, pi, Function
    >>> from sympy.differential_geometry.rn import R2, R2_p, R2_r
    >>> from sympy.differential_geometry import ScalarField, VectorField
    >>> x, y, r, theta = symbols('x, y, r, theta')
    >>> x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0')

    Points to be used as arguments for the field:
    >>> point_p = R2_p.point([r0, theta0])
    >>> point_r = R2_r.point([x0, y0])

    Scalar field to operate on:
    >>> g = Function('g')
    >>> s_field = g(R2.x, R2.y)
    >>> s_field(point_r)
    g(x0, y0)
    >>> s_field(point_p)
    g(r0*cos(theta0), r0*sin(theta0))

    Vector field:
    >>> v = VectorField(R2_r, [x, y], [1, 1])
    >>> v(s_field)(point_r)
    Derivative(g(x0, y0), x0) + Derivative(g(x0, y0), y0)
    >>> v(s_field)(point_p) # TODO this is correct but unusable
    Subs(Derivative(g(_x, r0*sin(theta0)), _x), (_x,), (r0*cos(theta0),)) + Subs(Derivative(g(r0*cos(theta0), _y), _y), (_y,), (r0*sin(theta0),))

    """
    def __init__(self, coord_sys, coords, components):
        super(VectorField, self).__init__()
        self._coord_sys = coord_sys
        coords, components = dummyfy(coords, components)
        self._coords = coords
        self._components = components
        self._args = self._coord_sys, self._coords, self._components

    def __call__(self, scalar_field):
        scalar_expr = scalar_field(Point(self._coord_sys, self._coords))
        diff_scalar_expr = (diff(scalar_expr, c) for c in self._coords)
        projected = sum(e[0]*e[1] for e in zip(diff_scalar_expr, self._components))
        # TODO This is the simplest to write, however is it the smartest
        # routine for applying vector field on a scalar field?
        # TODO Document the nontrivial jump-through-hoops that is done wrt to
        # coordinate systems here
        return ScalarField(self._coord_sys, self._coords, projected)


###############################################################################
# Integral curves on vector fields
###############################################################################
def intcurve_series(vector_field, param, start_point, n=6, coord_sys=None):
    """Return the series expansion for an integral curve of the field.

    Integral curve is a function `gamma` taking a parameter in R to a point
    in the manifold. It verifies the equation:

    `vector_field(f)(gamma(param)) = diff(f(gamma(t)), t)`

    for any value `t` for the parameter and any scalar field `f`.

    This function returns a series expansion of `gamma(t)` in terms of the
    coordinate system `coord_sys`. The equations and expansions are necessarily
    done in coordinate-system-dependent way as there is no other way to
    represent movement between points on the manifold (i.e. there is no such
    thing as a difference of points for a general manifold).

    Arguments:
    ==========

    - vector_field - the vector field for which an integral curve will be given
    - param - the argument of the function `gamma` from R to the curve
    - start_point - the point which coresponds to `gamma(0)`
    - n - the order to which to expand
    - coord_sys - the coordinate system in which to expand

    See Also: intcurve_diffequ

    Examples:
    =========

    Use the predefined R2 manifold:
    >>> from sympy.abc import t, x, y
    >>> from sympy.differential_geometry.rn import R2, R2_p, R2_r
    >>> from sympy.differential_geometry import intcurve_series

    Specify a starting point and a vector field:
    >>> start_point = R2_r.point([x, y])
    >>> vector_field = R2_r.d_dx

    Calculate the series:
    >>> series = intcurve_series(vector_field, t, start_point, n=3)
    >>> series[0]
    [x]
    [y]
    >>> series[1]
    [t]
    [0]
    >>> series[2]
    [0]
    [0]

    The series in the polar coordinate system: #TODO check by hand for correctness
    >>> series = intcurve_series(vector_field, t, start_point, 3, R2_p)
    >>> series[0]
    [sqrt(x**2 + y**2)]
    [      atan2(y, x)]
    >>> series[1]
    [t*x/sqrt(x**2 + y**2)]
    [   -t*y/(x**2 + y**2)]
    >>> series[2]
    [t**2*(-x**2/(x**2 + y**2)**(3/2) + 1/sqrt(x**2 + y**2))/2]
    [                                t**2*x*y/(x**2 + y**2)**2]

    """
    def iter_vfield(scalar_field, i):
        """Return `vector_field` called `i` times on `scalar_field`."""
        return reduce(lambda s, v: v(s), [vector_field,]*i, scalar_field)
    def taylor_terms_per_coord(coord_function):
        """Return the series for one of the coordinates."""
        return [param**i*iter_vfield(coord_function, i)(start_point)/factorial(i)
                for i in range(n)]
    coord_sys = coord_sys if coord_sys else start_point._coord_sys
    coord_functions = coord_sys.coord_functions()
    taylor_terms = zip(*[taylor_terms_per_coord(f) for f in coord_functions])
    return [Matrix(t) for t in taylor_terms]


def intcurve_diffequ(vector_field, param, start_point, coord_sys=None):
    """Return the differential equation for an integral curve of the field.

    Integral curve is a function `gamma` taking a parameter in R to a point
    in the manifold. It verifies the equation:

    `vector_field(f)(gamma(param)) = diff(f(gamma(t)), t)`

    for any value `t` for the parameter and any scalar field `f`.

    This function returns the differential equation of `gamma(t)` in terms of the
    coordinate system `coord_sys`. The equations and expansions are necessarily
    done in coordinate-system-dependent way as there is no other way to
    represent movement between points on the manifold (i.e. there is no such
    thing as a difference of points for a general manifold).

    Arguments:
    ==========

    - vector_field - the vector field for which an integral curve will be given
    - param - the argument of the function `gamma` from R to the curve
    - start_point - the point which coresponds to `gamma(0)`
    - coord_sys - the coordinate system in which to give the equations

    Returns:
    ========
    a tuple of (equations, initial conditions)

    See Also: intcurve_series

    Examples:
    =========

    Use the predefined R2 manifold:
    >>> from sympy.abc import t
    >>> from sympy.differential_geometry.rn import R2, R2_p, R2_r
    >>> from sympy.differential_geometry import intcurve_diffequ

    Specify a starting point and a vector field:
    >>> start_point = R2_r.point([1, 0])
    >>> vector_field = -R2.y*R2.d_dx + R2.x*R2.d_dy

    Get the equation:
    >>> equations, init_cond = intcurve_diffequ(vector_field, t, start_point)
    >>> equations
    [f_1(t) + Derivative(f_0(t), t), -f_0(t) + Derivative(f_1(t), t)]
    >>> init_cond
    [f_0(0) - 1, f_1(0)]

    The series in the polar coordinate system: #TODO check by hand for correctness
    >>> equations, init_cond = intcurve_diffequ(vector_field, t, start_point, R2_p)
    >>> #TODO correct but too complicated >>> equations
    >>> init_cond
    [f_0(0) - 1, f_1(0)]

    """
    coord_sys = coord_sys if coord_sys else start_point._coord_sys
    # TODO on the next line the argument in range is hilariously long.
    gammas = [Function('f_%d'%i)(param) for i in range(start_point._coord_sys.patch.manifold.dim)]
    arbitrary_p = Point(coord_sys, gammas)
    coord_functions = coord_sys.coord_functions()
    equations = [diff(cf(arbitrary_p), param) - vector_field(cf)(arbitrary_p)
                 for cf in coord_functions]
    init_cond = [cf(arbitrary_p).subs(param,0) - cf(start_point)
                 for cf in coord_functions]
    return equations, init_cond


###############################################################################
# Helpers
###############################################################################
def dummyfy(args, exprs):
    # TODO Is this a good idea?
    d_args = Matrix([s.as_dummy() for s in args])
    d_exprs = Matrix([sympify(expr).subs(zip(args, d_args)) for expr in exprs])
    return d_args, d_exprs
