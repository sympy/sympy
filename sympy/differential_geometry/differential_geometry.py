from sympy import sympify, Dummy, Matrix, Basic, Expr, solve, diff

# TODO order the imports and make them explicit
# TODO issue 2070: all the stuff about .args and rebuilding
# TODO better dosctrings and more distributed doctests
# TODO maybe a common class Field makes sense

class Manifold(Basic):
    """Object representing a mathematical manifold.

    The only role that it plays is to keep a list of all patches defined on the
    manifold.
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

    It serves as a container/parent for all coordinate system charts that can
    be defined on this patch.

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
    found by `solve`:
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

    >>> rect.coord_function(0)(p)
    -sqrt(2)/2
    >>> rect.coord_function(1)(p)
    sqrt(2)/2

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
        """Transform `coords` to coord system `to_sys`."""
        coords = Matrix(coords)
        if self != to_sys:
            transf = self.transforms[to_sys]
            coords = transf[1].subs(zip(transf[0], coords))
        return coords

    def coord_function(self, coord_index):
        """Return a ScalarField that takes a point and returns one of the coords."""
        args = [Dummy() for i in range(self.patch.manifold.dim)]
        result = args[coord_index]
        return ScalarField(self, args, result)

    def base_vector(self, coord_index):
        """Return a basis VectorField."""
        args = [Dummy() for i in range(self.patch.manifold.dim)]
        result = [0,] * self.patch.manifold.dim
        result[coord_index] = 1
        return VectorField(self, args, result)

    def point(self, coords):
        """Create a `Point` with coordinates given in this coord system."""
        return Point(self, coords)

    def point_to_coords(self, point):
        """Calculate the coordinates of a point in this coord system."""
        return point.coords(self)


class Point(Basic):
    """Point in a Manifold object.

    To define a point you must supply coordinates and a coordinate system.

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
        coords = point.coords(self._coord_sys)
        return self._expr.subs(zip(self._coords, coords))


class VectorField(Expr):
    """Vector Field over a Manifold.

    A vector field is an operator taking a scalar field and returning a
    directional derivative (which is also a scalar field).

    To define a vector field you need to choose a coordinate system and define
    the vector field in terms of that coordinate system.

    Examples:
    =========

    Use the predefined R2 manifold
    >>> from sympy import symbols, sin, cos, pi, Function
    >>> from sympy.differential_geometry.Rn import R2, R2_p, R2_r
    >>> from sympy.differential_geometry import ScalarField, VectorField
    >>> x, y, r, theta = symbols('x, y, r, theta')
    >>> x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0')

    Points to be used as arguments for the filed:
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
        return ScalarField(self._coord_sys, self._coords, projected)


def dummyfy(args, exprs):
    # TODO Is this a good idea?
    d_args = Matrix([s.as_dummy() for s in args])
    d_exprs = Matrix([sympify(expr).subs(zip(args, d_args)) for expr in exprs])
    return d_args, d_exprs
