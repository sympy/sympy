from sympy.matrices import Matrix
from sympy.core import Basic, Expr, Dummy, Function, sympify, diff, Mul, Add
from sympy.core.numbers import Zero
from sympy.solvers import solve
from sympy.functions import factorial
from sympy.simplify import simplify
from sympy.core.compatibility import reduce, permutations
from sympy.combinatorics import Permutation

# TODO you are a bit excessive in the use of Dummies
# TODO dummy point, literal field

class Manifold(Basic):
    """Object representing a mathematical manifold.

    The only role that this object plays is to keep a list of all patches
    defined on the manifold. It does not provide any means to study the
    topological characteristics of the manifold that it represents.

    """
    def __init__(self, name, dim):
        super(Manifold, self).__init__()
        self.name = name
        self.dim = dim
        self.patches = []
        # The patches list is necessary if a Patch instance needs to enumerate
        # other Patch instance on the same manifold.

    def _latex(self, printer, *args):
        return r'\mathbb{%s}' % self.name


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
    >>> from sympy.diffgeom import Manifold, Patch
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
        self.manifold = manifold
        self.manifold.patches.append(self)
        self.coord_systems = []
        # The list of coordinate systems is necessary for an instance of
        # CoordSystem to enumerate other coord systems on the patch.

    @property
    def dim(self):
        return self.manifold.dim

    def _latex(self, printer, *args):
        return r'\mathrm{%s}_{%s}' % (self.name, self.manifold._latex(printer, *args))


class CoordSystem(Basic):
    """Contains all coordinate transformation logic.

    Examples:
    =========

    Define a Manifold and a Patch, and then define two coord systems on that
    patch:
    >>> from sympy import symbols, sin, cos, pi
    >>> from sympy.diffgeom import Manifold, Patch, CoordSystem
    >>> r, theta = symbols('r, theta')
    >>> m = Manifold('M', 2)
    >>> patch = Patch('P', m)
    >>> rect = CoordSystem('rect', patch)
    >>> polar = CoordSystem('polar', patch)
    >>> rect in patch.coord_systems
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

    Calculate the jacobian of the polar to cartesian transformation:
    >>> polar.jacobian(rect, [r, theta])
    [cos(theta), -r*sin(theta)]
    [sin(theta),  r*cos(theta)]

    Define a point using coordinates in one of the coordinate systems:
    >>> p = polar.point([1, 3*pi/4])
    >>> rect.point_to_coords(p)
    [-sqrt(2)/2]
    [ sqrt(2)/2]

    Define a basis scalar field (i.e. a coordinate function), that takes a
    point and returns its coordinates. It is an instance of `BaseScalarField`.
    >>> rect.coord_function(0)(p)
    -sqrt(2)/2
    >>> rect.coord_function(1)(p)
    sqrt(2)/2

    Define a basis vector field (i.e. a unit vector field along the coordinate
    line). Vectors are also differential operators on scalar fields. It is an
    instance of `BaseVectorField`.
    >>> v_x = rect.base_vector(0)
    >>> x = rect.coord_function(0)
    >>> v_x(x)(p)
    1
    >>> v_x(v_x(x))(p)
    0

    Define a basis oneform field:
    >>> dx = rect.base_oneform(0)
    >>> dx(v_x)(p)
    1

    If you provide a list of names the fields will print nicely:
    - without provided names:
    >>> x, v_x, dx
    (rect_0, e_rect_0, drect_0)

    - with provided names
    >>> rect = CoordSystem('rect', patch, ['x', 'y'])
    >>> rect.coord_function(0), rect.base_vector(0), rect.base_oneform(0)
    (x, e_x, dx)

    """
    #  Contains a reference to the parent patch in order to be able to access
    # other coordinate system charts.
    def __init__(self, name, patch, names=None):
        super(CoordSystem, self).__init__()
        self.name = name
        if not names:
            names = ['%s_%d'%(name, i) for i in range(patch.dim)]
        self._names = names
        self.patch = patch
        self._args = self.name, self.patch
        # names is not in args because it is related only to printing, not to
        # identifying the CoordSystem instance.
        self.patch.coord_systems.append(self)
        self.transforms = {}
        # All the coordinate transformation logic is in this dictionary in the
        # form of:
        #  key = other coordinate system
        #  value = tuple of  # TODO make these Lambda instances
        #          - list of `Dummy` coordinates in this coordinate system
        #          - list of expressions as a function of the Dummies giving
        #          the coordinates in another coordinate system
        self._dummies = [Dummy(str(n)) for n in (names if names else range(self.patch.dim))]
        self._dummy = Dummy()

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

    def jacobian(self, to_sys, coords):
        """Return the jacobian matrix of a transformation."""
        return self.coord_transform_to(to_sys, coords).jacobian(coords)

    def coord_function(self, coord_index):
        """Return a `ScalarField` that takes a point and returns one of the coords.

        Takes a point and returns its coordinate in this coordinate system.

        See the docstring of `CoordSystem` for examples."""
        return BaseScalarField(self, coord_index)

    def coord_functions(self):
        """Returns a list of all coordinate functions.

        For more details see the coord_function method of this class."""
        return [self.coord_function(i) for i in range(self.dim)]

    def base_vector(self, coord_index):
        """Return a basis VectorField.

        The basis vector field for this coordinate system. It is also an
        operator on scalar fields.

        See the docstring of `CoordSystem` for examples."""
        return BaseVectorField(self, coord_index)

    def base_vectors(self):
        """Returns a list of all base vectors.

        For more details see the base_vector method of this class."""
        return [self.base_vector(i) for i in range(self.dim)]

    def base_oneform(self, coord_index):
        """Return a basis OneFormField.

        The basis one-form field for this coordinate system. It is also an
        operator on vector fields.

        See the docstring of `CoordSystem` for examples."""
        return Differential(self.coord_function(coord_index))

    def base_oneforms(self):
        """Returns a list of all base oneforms.

        For more details see the base_oneform method of this class."""
        return [self.base_oneform(i) for i in range(self.dim)]

    def point(self, coords):
        """Create a `Point` with coordinates given in this coord system.

        See the docstring of `CoordSystem` for examples."""
        return Point(self, coords)

    def point_to_coords(self, point):
        """Calculate the coordinates of a point in this coord system.

        See the docstring of `CoordSystem` for examples."""
        return point.coords(self)

    @property
    def dim(self):
        return self.patch.dim

    def _latex(self, printer, *args):
        return r'\mathrm{%s}^{\mathrm{%s}}_{%s}' % (
                self.name, self.patch.name, self.patch.manifold._latex(printer, *args))

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
    >>> from sympy.diffgeom import (
    ...        Manifold, Patch, CoordSystem, Point)
    >>> r, theta = symbols('r, theta')
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
        self._args = self._coord_sys, self._coords

    def coords(self, to_sys=None):
        """Coordinates of the point in a given coordinate system.

        If `to_sys` is None it returns the coordinates in the system in
        which the point was defined."""
        if to_sys:
            return self._coord_sys.coord_transform_to(to_sys, self._coords)
        else:
            return self._coords

    @property
    def free_symbols(self):
        return self._coords.free_symbols


class BaseScalarField(Expr):
    """Base Scalar Field over a Manifold for a given Coordinate System.

    A scalar field takes a point as an argument and returns a scalar.

    A base scalar field of a coordinate system takes a point and returns one of
    the coordinates of that point in the coordinate system in question.

    To define a scalar field you need to choose the coordinate system and the
    index of the coordinate.

    The use of the scalar field after its definition is independent of the
    coordinate system in which it was defined, however due to limitations in
    the simplification routines you may arrive at more complicated
    expression if you use unappropriate coordinate systems.

    You can build complicated scalar fields by just building up SymPy
    expressions containing ``BaseScalarField`` instances.

    Examples
    ========

    Define boilerplate Manifold, Patch and coordinate systems:
    >>> from sympy import symbols, sin, cos, pi, Function
    >>> from sympy.diffgeom import (
    ...        Manifold, Patch, CoordSystem, Point, BaseScalarField)
    >>> r0, theta0 = symbols('r0, theta0')
    >>> m = Manifold('M', 2)
    >>> p = Patch('P', m)
    >>> rect = CoordSystem('rect', p)
    >>> polar = CoordSystem('polar', p)
    >>> polar.connect_to(rect, [r0, theta0], [r0*cos(theta0), r0*sin(theta0)])

    Point to be used as an argument for the filed:
    >>> point = polar.point([r0, 0])

    Examples of fields:
    >>> fx = BaseScalarField(rect, 0)
    >>> fy = BaseScalarField(rect, 1)
    >>> (fx**2+fy**2)(point)
    r0**2

    >>> g = Function('g')
    >>> ftheta = BaseScalarField(polar, 1)
    >>> fg = g(ftheta-pi)
    >>> fg(point)
    g(-pi)

    """
    def __init__(self, coord_sys, index):
        super(BaseScalarField, self).__init__()
        self._coord_sys = coord_sys
        self._index = index
        self._args = self._coord_sys, self._index

    def __call__(self, *args):
        point = args[0]
        if not isinstance(point, Point):
            # We want the recursive calling mechanics for all fields hence
            # we need to return ScalarField itself if the arg is not a Point.
            return self
        coords = point.coords(self._coord_sys)
        # XXX Calling doit  is necessary with all the Subs expressions
        # XXX Calling simplify is necessary with all the trig expressions
        return simplify(coords[self._index]).doit()

    # Workaround for limitations on the content of args
    free_symbols=set([])
    def doit(self):
        return self


class BaseVectorField(Expr):
    r"""Vector Field over a Manifold.

    A vector field is an operator taking a scalar field and returning a
    directional derivative (which is also a scalar field).

    A base vector field is the same type of operator, however the derivation is
    specifically done wrt a chosen coordinate.

    To define a base vector field you need to choose the coordinate system and
    the index of the coordinate.

    The use of the vector field after its definition is independent of the
    coordinate system in which it was defined, however due to limitations in
    the simplification routines you may arrive at more complicated
    expression if you use unappropriate coordinate systems.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy import symbols, pi, Function
    >>> from sympy.diffgeom.Rn import R2, R2_p, R2_r
    >>> from sympy.diffgeom import BaseVectorField
    >>> from sympy import pprint
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
    >>> v = BaseVectorField(R2_r, 1)
    >>> pprint(v(s_field))
    /  d              \|
    |-----(g(x, xi_2))||
    \dxi_2            /|xi_2=y
    >>> pprint(v(s_field)(point_r).doit())
     d
    ---(g(x0, y0))
    dy0
    >>> pprint(v(s_field)(point_p).doit())
    /  d                           \|
    |-----(g(r0*cos(theta0), xi_2))||
    \dxi_2                         /|xi_2=r0*sin(theta0)

    """
    def __init__(self, coord_sys, index):
        super(BaseVectorField, self).__init__()
        self._coord_sys = coord_sys
        self._index = index
        self._args = self._coord_sys, self._index

    def __call__(self, scalar_field):
        base_scalars = list(scalar_field.atoms(BaseScalarField))

        # First step: e_x(x+r**2) -> e_x(x) + 2*r*e_x(r)
        d_var = self._coord_sys._dummy
        # TODO: you need a real dummy function for the next line
        d_funcs = [Function('_#_%s' % i)(d_var) for i, b in enumerate(base_scalars)]
        d_result = scalar_field.subs(zip(base_scalars, d_funcs))
        d_result = d_result.diff(d_var)

        # Second step: e_x(x) -> 1 and e_x(r) -> cos(atan2(x, y))
        coords = self._coord_sys._dummies
        d_funcs_deriv = [f.diff(d_var) for f in d_funcs]
        d_funcs_deriv_sub = []
        for b in base_scalars:
            jac = self._coord_sys.jacobian(b._coord_sys, coords)
            d_funcs_deriv_sub.append(jac[b._index, self._index])
        d_result = d_result.subs(zip(d_funcs_deriv, d_funcs_deriv_sub))

        # Remove the dummies
        result = d_result.subs(zip(d_funcs, base_scalars))
        result = result.subs(zip(coords, self._coord_sys.coord_functions()))
        return result.doit() # XXX doit for the Subs instances


class Commutator(Expr):
    r"""Commutator of two vector fields.

    The commutator of two vector fields `v_1` and `v_2` is defined as the
    vector field `[v_1, v_2]` that evaluated on each scalar field `f` is equal
    to `v_1(v_2(f)) - v_2(v_1(f))`.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import Commutator
    >>> from sympy import pprint
    >>> from sympy.simplify import simplify

    Vector fields:
    >>> e_x, e_y, e_r = R2.e_x, R2.e_y, R2.e_r
    >>> c_xy = Commutator(e_x, e_y)
    >>> c_xr = Commutator(e_x, e_r)

    >>> c_xy(R2.x + R2.y**2)
    0

    """
    # TODO simplify fails with an error
    #>>> pprint(simplify(c_xr(R2.y**2).doit()))
    #              -1
    #     / 2    2\
    #-2*y*\x  + y /  *cos(theta)*y

    #"""
    def __init__(self, v1, v2):
        super(Commutator, self).__init__()
        self._args = (v1, v2)
        self._v1 = v1
        self._v2 = v2

    def __call__(self, scalar_field):
        return self._v1(self._v2(scalar_field)) - self._v2(self._v1(scalar_field))


class Differential(Expr):
    """Return the differential of a form field.

    The differential of a form (i.e. the exterior derivative) has a complicated
    definition in the general case.

    The differential `df` of the 0-form `f` is defined for any vector field `v`
    as `df(v) = v(f)`.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy import Function
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import Differential
    >>> from sympy import pprint

    Scalar field (0-forms):
    >>> g = Function('g')
    >>> s_field = g(R2.x, R2.y)

    Vector fields:
    >>> e_x, e_y, = R2.e_x, R2.e_y

    Differentials:
    >>> dg = Differential(s_field)
    >>> dg
    d(g(x, y))
    >>> pprint(dg(e_x))
    /  d              \|
    |-----(g(xi_1, y))||
    \dxi_1            /|xi_1=x
    >>> pprint(dg(e_y))
    /  d              \|
    |-----(g(x, xi_2))||
    \dxi_2            /|xi_2=y

    Applying the exterior derivative operator twice always results in:
    >>> Differential(dg)
    0

    """
    def __new__(cls, form_field):
        if isinstance(form_field, Differential):
            return sympify(0)
        else:
            return super(Differential, cls).__new__(cls, form_field)

    def __init__(self, form_field):
        super(Differential, self).__init__()
        self._form_field = form_field
        self._args = (self._form_field, )

    def __call__(self, *vector_fields):
        k = len(vector_fields)
        if k==1:
            return vector_fields[0](self._form_field)
        else:
            f = self._form_field
            v = vector_fields
            ret = 0
            for i in range(k):
                t = v[i](f(*v[:i]+v[i+1:]))
                ret += (-1)**i*t
                for j in range(i+1,k):
                    c = Commutator(v[i], v[j])
                    t = f(*(c,)+v[:i]+v[i+1:j]+v[j+1:])
                    ret += (-1)**(i+j)*t
            return ret


class TensorProduct(Expr):
    """Tensor product of forms.

    The tensor product permits the creation of multilinear functionals (i.e.
    higher order forms) out of lower order forms (e.g. 1-forms). However, the
    higher forms thus created lack the interesting features provided by the
    other type of product, the wedge product.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy import Function
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import TensorProduct
    >>> from sympy import pprint

    >>> TensorProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y)
    1
    >>> TensorProduct(R2.dx, R2.dy)(R2.e_y, R2.e_x)
    0
    >>> TensorProduct(R2.dx, R2.x*R2.dy)(R2.x*R2.e_x, R2.e_y)
    x**2

    You can nest tensor products.
    >>> tp1 = TensorProduct(R2.dx, R2.dy)
    >>> TensorProduct(tp1, R2.dx)(R2.e_x, R2.e_y, R2.e_x)
    1

    """
    def __init__(self, *args):
        super(TensorProduct, self).__init__()
        self._args = args

    def __call__(self, *v_fields):
        orders = [order_of_form(f) for f in self._args]
        indices = [sum(orders[:i+1]) for i in range(len(orders)-1)]
        v_fields = [v_fields[i:j] for i, j in zip([0]+indices, indices+[None])]
        return Mul(*[t[0](*t[1]) for t in zip(self._args, v_fields)])

    def _latex(self, printer, *args):
        elements = [printer._print(a) for a in self.args]
        return r'\otimes'.join(elements)


class WedgeProduct(TensorProduct):
    """Wedge product of forms.

    In the context of integration only completely antisymetric forms make
    sense. The wedge product permits the creation of such forms.

    Examples:
    =========

    Use the predefined R2 manifold, setup some boilerplate.
    >>> from sympy import Function
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import WedgeProduct
    >>> from sympy import pprint

    >>> WedgeProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y)
    1
    >>> WedgeProduct(R2.dx, R2.dy)(R2.e_y, R2.e_x)
    -1
    >>> WedgeProduct(R2.dx, R2.x*R2.dy)(R2.x*R2.e_x, R2.e_y)
    x**2

    You can nest wedge products.
    >>> wp1 = WedgeProduct(R2.dx, R2.dy)
    >>> WedgeProduct(wp1, R2.dx)(R2.e_x, R2.e_y, R2.e_x)
    0

    """
    # TODO the caclulation of signatures is slow
    # TODO you do not need all these permutations (neither the prefactor)
    def __call__(self, *vector_fields):
        orders = (order_of_form(e) for e in self.args)
        mul = 1/Mul(*(factorial(o) for o in orders))
        perms = permutations(vector_fields)
        perms_par = (Permutation(p).signature() for p in permutations(range(len(vector_fields))))
        tensor_prod = TensorProduct(*self.args)
        return mul*Add(*[tensor_prod(*p[0])*p[1] for p in zip(perms, perms_par)])


class LieDerivative(Expr):
    """Lie derivative operator wrt a vector field.

    The transport operator that defines the Lie derivative is the pushforward
    of the field to be derived along the integral curve of the field wrt which
    one derives.

    Examples:
    =========

    >>> #TODO
    """
    def __new__(cls, expr, v_field):
        expr_form_ord = order_of_form(expr)
        if expr_form_ord>0:
            return super(LieDerivative, cls).__new__(cls, v_field, expr)
        if arg.atoms(BaseVectorField):
            return Commutator(v_field, expr)
        else:
            return v_field(expr)

    def __init__(self, expr, v_field):
        super(LieDerivative, self).__init__()
        self._v_field = v_field
        self._expr = expr
        self._args = (self._expr, self._v_field)

    def __call__(self, *args):
        v = self._v_field
        expr = self._expr
        lead_term = v(expr(*args))
        rest = Add(*[Mul(*args[:i] + (Commutator[v, args[i]],) + args[i+1:]) for i in range(len(args))])
        return lead_term - rest


###############################################################################
# Integral curves on vector fields
###############################################################################
def intcurve_series(vector_field, param, start_point, n=6, coord_sys=None, coeffs=False):
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
    - coeffs (default False) - if True return a list of elements of the expansion

    See Also: intcurve_diffequ

    Examples:
    =========

    Use the predefined R2 manifold:
    >>> from sympy.abc import t, x, y
    >>> from sympy.diffgeom.Rn import R2, R2_p, R2_r
    >>> from sympy.diffgeom import intcurve_series

    Specify a starting point and a vector field:
    >>> start_point = R2_r.point([x, y])
    >>> vector_field = R2_r.e_x

    Calculate the series:
    >>> intcurve_series(vector_field, t, start_point, n=3)
    [t + x]
    [    y]

    Or get the elements of the expansion in a list:
    >>> series = intcurve_series(vector_field, t, start_point, n=3, coeffs=True)
    >>> series[0]
    [x]
    [y]
    >>> series[1]
    [t]
    [0]
    >>> series[2]
    [0]
    [0]

    The series in the polar coordinate system:
    >>> series = intcurve_series(vector_field, t, start_point, n=3, coord_sys=R2_p, coeffs=True)
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
    taylor_terms = [taylor_terms_per_coord(f) for f in coord_functions]
    if coeffs:
        return [Matrix(t) for t in zip(*taylor_terms)]
    else:
        return Matrix([sum(c) for c in taylor_terms])


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
    >>> from sympy.diffgeom.Rn import R2, R2_p, R2_r
    >>> from sympy.diffgeom import intcurve_diffequ

    Specify a starting point and a vector field:
    >>> start_point = R2_r.point([0, 1])
    >>> vector_field = -R2.y*R2.e_x + R2.x*R2.e_y

    Get the equation:
    >>> equations, init_cond = intcurve_diffequ(vector_field, t, start_point)
    >>> equations
    [f_1(t) + Derivative(f_0(t), t), -f_0(t) + Derivative(f_1(t), t)]
    >>> init_cond
    [f_0(0), f_1(0) - 1]

    The series in the polar coordinate system:
    >>> equations, init_cond = intcurve_diffequ(vector_field, t, start_point, R2_p)
    >>> equations
    [Derivative(f_0(t), t), Derivative(f_1(t), t) - 1]
    >>> init_cond
    [f_0(0) - 1, f_1(0) - pi/2]

    """
    coord_sys = coord_sys if coord_sys else start_point._coord_sys
    gammas = [Function('f_%d'%i)(param) for i in range(start_point._coord_sys.dim)]
    arbitrary_p = Point(coord_sys, gammas)
    coord_functions = coord_sys.coord_functions()
    equations = [simplify(diff(cf(arbitrary_p), param) - vector_field(cf)(arbitrary_p))
                 for cf in coord_functions]
    init_cond = [simplify(cf(arbitrary_p).subs(param,0) - cf(start_point))
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


def order_of_form(expr):
    # TODO move some of this to class methods.
    if isinstance(expr, Add):
        orders = [order_of_form(e) for e in expr.args]
        if len(set(orders)) != 1:
            raise ValueError('Misformed expression containing form fields of varying order.')
        return order[0]
    elif isinstance(expr, Mul):
        orders = [order_of_form(e) for e in expr.args]
        not_zero = [o for o in orders if o != 0]
        if len(not_zero) != 1:
            raise ValueError('Misformed expression containing multiplication between forms.')
        return 0 if not not_zero else not_zero[0]
    elif isinstance(expr, Differential):
        return order_of_form(*expr.args) + 1
    elif isinstance(expr, TensorProduct):
        return sum(order_of_form(a) for a in expr.args)
    else:
        return 0


###############################################################################
# Coordinate-dependent functions
###############################################################################
def twoform_to_matrix(expr):
    """Return the matrix representing the twoform.

    For the twoform `w` return the matrix `M` such that `M[i,j]=w(e_i, e_j)`,
    where `e_i` is the i-th base vector field for the coordinate system in
    which the expression of `w` is given.

    Examples:
    =========

    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import twoform_to_matrix, TensorProduct
    >>> TP = TensorProduct
    >>> twoform_to_matrix(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [1, 0]
    [0, 1]
    >>> twoform_to_matrix(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [x, 0]
    [0, 1]
    >>> twoform_to_matrix(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy) - TP(R2.dx, R2.dy)/2)
    [   1, 0]
    [-1/2, 1]

    """
    coord_sys = expr.atoms(CoordSystem)
    if len(coord_sys) != 1:
        raise ValueError('The input expression concerns more than one '
                         'coordinate systems, hence there is no unambiguous '
                         'way to choose a coordinate system for the matrix.')
    coord_sys = coord_sys.pop()
    vectors = coord_sys.base_vectors()
    expr = expr.expand()
    matrix_content = [[expr(v1, v2) for v1 in vectors]
                                    for v2 in vectors]
    return Matrix(matrix_content)


def metric_to_Christoffel_1st(expr):
    """Return the nested list of Christoffel symbols for the given metric.

    This returns the Christoffel symbol of first kind.

    Examples:
    =========

    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import metric_to_Christoffel_1st, TensorProduct
    >>> TP = TensorProduct
    >>> metric_to_Christoffel_1st(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    >>> metric_to_Christoffel_1st(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[[1/2, 0], [0, 0]], [[0, 0], [0, 0]]]

    """
    matrix = twoform_to_matrix(expr)
    coord_sys = expr.atoms(CoordSystem).pop()
    deriv_matrices = [matrix.applyfunc(lambda a: d(a)) for d in coord_sys.base_vectors()]
    indices = range(coord_sys.dim)
    christoffel = [[[(deriv_matrices[k][i,j] + deriv_matrices[j][i,k] - deriv_matrices[i][j,k])/2
                     for k in indices]
                     for j in indices]
                     for i in indices]
    return christoffel


def metric_to_Christoffel_2nd(expr):
    """Return the nested list of Christoffel symbols for the given metric.

    This returns the Christoffel symbol of second kind.

    Examples:
    =========

    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import metric_to_Christoffel_2nd, TensorProduct
    >>> TP = TensorProduct
    >>> metric_to_Christoffel_2nd(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    >>> metric_to_Christoffel_2nd(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[[x**(-1)/2, 0], [0, 0]], [[0, 0], [0, 0]]]

    """
    ch_1st = metric_to_Christoffel_1st(expr)
    coord_sys = expr.atoms(CoordSystem).pop()
    indices = range(coord_sys.dim)
    # XXX workaround, inverting a matrix does not work if it contains non
    # symbols
    #matrix = twoform_to_matrix(expr).inv()
    matrix = twoform_to_matrix(expr)
    s_fields = set()
    for e in matrix:
        s_fields.update(e.atoms(BaseScalarField))
    s_fields = list(s_fields)
    dums = coord_sys._dummies
    matrix = matrix.subs(zip(s_fields, dums)).inv().subs(zip(dums, s_fields))
    # XXX end of workaround
    christoffel = [[[Add(*[matrix[i,l]*ch_1st[l][j][k] for l in indices])
                     for k in indices]
                     for j in indices]
                     for i in indices]
    return christoffel


def metric_to_Riemann_components(expr):
    """Return the components of the Riemann tensor expressed in a given basis.

    Given a metric it calculates the components of the Riemann tensor in the
    canonical basis of the coordinate system in which the metric expression is
    given.

    Examples:
    =========

    >>> from sympy import pprint, exp
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import metric_to_Riemann_components, TensorProduct
    >>> TP = TensorProduct
    >>> metric_to_Riemann_components(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[[[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]]

    >>> non_trivial_metric = exp(2*R2.r)*TP(R2.dr, R2.dr) + R2.r**2*TP(R2.dtheta, R2.dtheta)
    >>> non_trivial_metric
    exp(2*r)*TensorProduct(dr, dr) + r**2*TensorProduct(dtheta, dtheta)
    >>> riemann = metric_to_Riemann_components(non_trivial_metric)
    >>> riemann[0]
    [[[0, 0], [0, 0]], [[0, -exp(-2*r)*r + 2*r*exp(-2*r)], [exp(-2*r)*r - 2*r*exp(-2*r), 0]]]
    >>> riemann[1]
    [[[0, -r**(-1)], [r**(-1), 0]], [[0, 0], [0, 0]]]

    """
    ch_2nd = metric_to_Christoffel_2nd(expr)
    coord_sys = expr.atoms(CoordSystem).pop()
    indices = range(coord_sys.dim)
    deriv_ch = [[[[d(ch_2nd[i][j][k])
                     for d in coord_sys.base_vectors()]
                     for k in indices]
                     for j in indices]
                     for i in indices]
    riemann_a = [[[[deriv_ch[rho][sig][nu][mu] - deriv_ch[rho][sig][mu][nu]
                     for nu in indices]
                     for mu in indices]
                     for sig in indices]
                     for rho in indices]
    riemann_b = [[[[Add(*[ch_2nd[rho][l][mu]*ch_2nd[l][sig][nu] - ch_2nd[rho][l][nu]*ch_2nd[l][sig][mu] for l in indices])
                     for nu in indices]
                     for mu in indices]
                     for sig in indices]
                     for rho in indices]
    riemann = [[[[riemann_a[rho][sig][mu][nu] + riemann_b[rho][sig][mu][nu]
                     for nu in indices]
                     for mu in indices]
                     for sig in indices]
                     for rho in indices]
    return riemann


def metric_to_Ricci_components(expr):
    """Return the components of the Ricci tensor expressed in a given basis.

    Given a metric it calculates the components of the Ricci tensor in the
    canonical basis of the coordinate system in which the metric expression is
    given.

    Examples:
    =========

    >>> from sympy import pprint, exp
    >>> from sympy.diffgeom.Rn import R2
    >>> from sympy.diffgeom import metric_to_Ricci_components, TensorProduct
    >>> TP = TensorProduct
    >>> metric_to_Ricci_components(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    [[0, 0], [0, 0]]

    >>> non_trivial_metric = exp(2*R2.r)*TP(R2.dr, R2.dr) + R2.r**2*TP(R2.dtheta, R2.dtheta)
    >>> non_trivial_metric
    exp(2*r)*TensorProduct(dr, dr) + r**2*TensorProduct(dtheta, dtheta)
    >>> metric_to_Ricci_components(non_trivial_metric) #TODO why is this not simpler
    [[r**(-1), 0], [0, -exp(-2*r)*r + 2*r*exp(-2*r)]]

    """
    riemann = metric_to_Riemann_components(expr)
    coord_sys = expr.atoms(CoordSystem).pop()
    indices = range(coord_sys.dim)
    ricci = [[Add(*[riemann[k][i][k][j] for k in indices])
                     for j in indices]
                     for i in indices]
    return ricci
