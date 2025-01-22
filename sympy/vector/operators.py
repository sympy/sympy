import collections
from sympy.core.expr import Expr
from sympy.core import sympify, S, preorder_traversal
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.vector import Vector, VectorMul, VectorAdd, Cross, Dot
from sympy.core.function import Derivative
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.functions.elementary.trigonometric import sin


def _get_coord_systems(expr):
    g = preorder_traversal(expr)
    ret = set()
    for i in g:
        if isinstance(i, CoordSys3D):
            ret.add(i)
            g.skip()
    return frozenset(ret)


def _split_mul_args_wrt_coordsys(expr):
    d = collections.defaultdict(lambda: S.One)
    for i in expr.args:
        d[_get_coord_systems(i)] *= i
    return list(d.values())


class Gradient(Expr):
    """
    Represents unevaluated Gradient.

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, Gradient
    >>> R = CoordSys3D('R')
    >>> s = R.x*R.y*R.z
    >>> Gradient(s)
    Gradient(R.x*R.y*R.z)

    """

    def __new__(cls, expr):
        expr = sympify(expr)
        obj = Expr.__new__(cls, expr)
        obj._expr = expr
        return obj

    def doit(self, **hints):
        return gradient(self._expr, doit=True)


class Divergence(Expr):
    """
    Represents unevaluated Divergence.

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, Divergence
    >>> R = CoordSys3D('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Divergence(v)
    Divergence(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k)

    """

    def __new__(cls, expr):
        expr = sympify(expr)
        obj = Expr.__new__(cls, expr)
        obj._expr = expr
        return obj

    def doit(self, **hints):
        return divergence(self._expr, doit=True)


class Curl(Expr):
    """
    Represents unevaluated Curl.

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, Curl
    >>> R = CoordSys3D('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Curl(v)
    Curl(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k)

    """

    def __new__(cls, expr):
        expr = sympify(expr)
        obj = Expr.__new__(cls, expr)
        obj._expr = expr
        return obj

    def doit(self, **hints):
        return curl(self._expr, doit=True)


def curl(vect, doit=True):
    """
    Returns the curl of a vector field computed wrt the base scalars
    of the given coordinate system.

    Parameters
    ==========

    vect : Vector
        The vector operand

    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, curl
    >>> R = CoordSys3D('R')
    >>> v1 = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> curl(v1)
    0
    >>> v2 = R.x*R.y*R.z*R.i
    >>> curl(v2)
    R.x*R.y*R.j + (-R.x*R.z)*R.k

    """

    coord_sys = _get_coord_systems(vect)

    if len(coord_sys) == 0:
        return Vector.zero
    elif len(coord_sys) == 1:
        coord_sys = next(iter(coord_sys))
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        h1, h2, h3 = coord_sys.lame_coefficients()
        vectx = vect.dot(i)
        vecty = vect.dot(j)
        vectz = vect.dot(k)
        outvec = Vector.zero
        outvec += (Derivative(vectz * h3, y) -
                   Derivative(vecty * h2, z)) * i / (h2 * h3)
        outvec += (Derivative(vectx * h1, z) -
                   Derivative(vectz * h3, x)) * j / (h1 * h3)
        outvec += (Derivative(vecty * h2, x) -
                   Derivative(vectx * h1, y)) * k / (h2 * h1)

        if doit:
            return outvec.doit()
        return outvec
    else:
        if isinstance(vect, (Add, VectorAdd)):
            from sympy.vector import express
            try:
                cs = next(iter(coord_sys))
                args = [express(i, cs, variables=True) for i in vect.args]
            except ValueError:
                args = vect.args
            return VectorAdd.fromiter(curl(i, doit=doit) for i in args)
        elif isinstance(vect, (Mul, VectorMul)):
            vector = [i for i in vect.args if isinstance(i, (Vector, Cross, Gradient))][0]
            scalar = Mul.fromiter(i for i in vect.args if not isinstance(i, (Vector, Cross, Gradient)))
            res = Cross(gradient(scalar), vector).doit() + scalar*curl(vector, doit=doit)
            if doit:
                return res.doit()
            return res
        elif isinstance(vect, (Cross, Curl, Gradient)):
            return Curl(vect)
        else:
            raise ValueError("Invalid argument for curl")


def divergence(vect, doit=True):
    """
    Returns the divergence of a vector field computed wrt the base
    scalars of the given coordinate system.

    Parameters
    ==========

    vector : Vector
        The vector operand

    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, divergence
    >>> R = CoordSys3D('R')
    >>> v1 = R.x*R.y*R.z * (R.i+R.j+R.k)

    >>> divergence(v1)
    R.x*R.y + R.x*R.z + R.y*R.z
    >>> v2 = 2*R.y*R.z*R.j
    >>> divergence(v2)
    2*R.z

    """
    coord_sys = _get_coord_systems(vect)
    if len(coord_sys) == 0:
        return S.Zero
    elif len(coord_sys) == 1:
        if isinstance(vect, (Cross, Curl, Gradient)):
            return Divergence(vect)
        # TODO: is case of many coord systems, this gets a random one:
        coord_sys = next(iter(coord_sys))
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        h1, h2, h3 = coord_sys.lame_coefficients()
        vx = _diff_conditional(vect.dot(i), x, h2, h3) \
             / (h1 * h2 * h3)
        vy = _diff_conditional(vect.dot(j), y, h3, h1) \
             / (h1 * h2 * h3)
        vz = _diff_conditional(vect.dot(k), z, h1, h2) \
             / (h1 * h2 * h3)
        res = vx + vy + vz
        if doit:
            return res.doit()
        return res
    else:
        if isinstance(vect, (Add, VectorAdd)):
            return Add.fromiter(divergence(i, doit=doit) for i in vect.args)
        elif isinstance(vect, (Mul, VectorMul)):
            vector = [i for i in vect.args if isinstance(i, (Vector, Cross, Gradient))][0]
            scalar = Mul.fromiter(i for i in vect.args if not isinstance(i, (Vector, Cross, Gradient)))
            res = Dot(vector, gradient(scalar)) + scalar*divergence(vector, doit=doit)
            if doit:
                return res.doit()
            return res
        elif isinstance(vect, (Cross, Curl, Gradient)):
            return Divergence(vect)
        else:
            raise ValueError("Invalid argument for divergence")


def gradient(scalar_field, doit=True):
    """
    Returns the vector gradient of a scalar field computed wrt the
    base scalars of the given coordinate system.
    
    Parameters
    ==========
    scalar_field : SymPy Expr
        The scalar field to compute the gradient of
    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances
        
    Returns
    =======
    Vector
        The gradient vector
    """
    from sympy import S
    
    # Handle zero and constant cases
    if scalar_field == 0 or scalar_field == S.Zero:
        return Vector.zero
    if scalar_field.is_constant():
        return Vector.zero
        
    coord_sys = _get_coord_systems(scalar_field)
    if len(coord_sys) == 0:
        return Vector.zero
    elif len(coord_sys) == 1:
        coord_sys = next(iter(coord_sys))
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        h1, h2, h3 = coord_sys.lame_coefficients()
        
        # Handle coordinate transformations
        if coord_sys.transformation == 'cylindrical':
            # r, θ, z components
            gx = Derivative(scalar_field, x) / h1
            gy = Derivative(scalar_field, y) / (x * h2)  # Note the 1/r term
            gz = Derivative(scalar_field, z) / h3
        elif coord_sys.transformation == 'spherical':
            # r, θ, φ components
            gx = Derivative(scalar_field, x) / h1
            gy = Derivative(scalar_field, y) / (x * h2)  # 1/r term
            gz = Derivative(scalar_field, z) / (x * sin(y) * h3)  # 1/(r*sin(θ)) term
        else:  # Cartesian coordinates
            gx = Derivative(scalar_field, x) / h1
            gy = Derivative(scalar_field, y) / h2
            gz = Derivative(scalar_field, z) / h3
            
        res = gx * i + gy * j + gz * k
        if doit:
            return res.doit()
        return res
    else:
        raise ValueError("Gradient can only handle one coordinate system.")

def directional_derivative(vect1, vect2, doit=True):
    """
    Returns the directional derivative of vect1 along vect2,
    accounting for curvilinear coordinate systems.
    
    Parameters
    ==========
    vect1 : Vector
        The vector to differentiate
    vect2 : Vector
        The direction vector
    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances
        
    Returns
    =======
    Vector
        The directional derivative
    """
    coord_sys = _get_coord_systems(vect1)
    if len(coord_sys) != 1:
        raise ValueError("Directional derivative is only supported for a single coordinate system.")
    coord_sys = next(iter(coord_sys))
    
    # Handle zero vector cases
    if vect1 == Vector.zero or vect2 == Vector.zero:
        return Vector.zero
        
    base_vectors = coord_sys.base_vectors()
    base_scalars = coord_sys.base_scalars()
    h = coord_sys.lame_coefficients()
    
    # Express vect2 in terms of the coordinate system
    vect2 = vect2.express(coord_sys)
    
    # Initialize result vector
    result = Vector.zero
    
    if coord_sys.transformation == 'cylindrical':
        r = base_scalars[0]  # radial coordinate
        # Add centripetal terms for cylindrical coordinates
        v_theta = vect1.dot(base_vectors[1])
        result += (-v_theta**2/r) * base_vectors[0]
    elif coord_sys.transformation == 'spherical':
        # Add appropriate geometric terms for spherical coordinates
        pass
        
    # Add standard derivative terms
    for i, base_vec in enumerate(base_vectors):
        for j, base_vec2 in enumerate(base_vectors):
            term = vect2.dot(base_vec2) * Derivative(vect1.dot(base_vec), base_scalars[j])
            term /= h[j] * h[i]
            result += term * base_vec
            
    if doit:
        return result.doit()
    return result


class Laplacian(Expr):
    """
    Represents unevaluated Laplacian.

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, Laplacian
    >>> R = CoordSys3D('R')
    >>> v = 3*R.x**3*R.y**2*R.z**3
    >>> Laplacian(v)
    Laplacian(3*R.x**3*R.y**2*R.z**3)

    """

    def __new__(cls, expr):
        expr = sympify(expr)
        obj = Expr.__new__(cls, expr)
        obj._expr = expr
        return obj

    def doit(self, **hints):
        from sympy.vector.functions import laplacian
        return laplacian(self._expr)


def _diff_conditional(expr, base_scalar, coeff_1, coeff_2):
    """
    First re-expresses expr in the system that base_scalar belongs to.
    If base_scalar appears in the re-expressed form, differentiates
    it wrt base_scalar.
    Else, returns 0
    """
    from sympy.vector.functions import express
    new_expr = express(expr, base_scalar.system, variables=True)
    arg = coeff_1 * coeff_2 * new_expr
    return Derivative(arg, base_scalar) if arg else S.Zero
