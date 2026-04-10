from __future__ import annotations
import collections
from sympy.core.expr import Expr
from sympy.core import sympify, S, preorder_traversal
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.dyadic import Dyadic, DyadicAdd
from sympy.vector.vector import Vector, VectorMul, VectorAdd, Cross, Dot
from sympy.core.function import Derivative
from sympy.core.add import Add
from sympy.core.mul import Mul


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


def _christoffel_symbol_2nd_kind(h, q, i, j, k):
    """
    Compute the Christoffel symbol of the second kind, Gamma^{k}_{ij},
    for an orthogonal system.
    """
    if (i == j) and (i == k):
        return Derivative(h[i], q[i]) / h[i]
    if (i == k) and (i != j):
        return Derivative(h[i], q[j]) / h[i]
    if (j == k) and (i != j):
        return Derivative(h[j], q[i]) / h[j]
    if (i == j) and (i != k):
        return -Derivative(h[i], q[k]) * h[i] / h[k]**2
    return S.Zero


def gradient(field, doit=True):
    """
    Returns the gradient of the field computed wrt the base scalars
    of the given coordinate system.

    Parameters
    ==========

    field : Expr, Vector
        The scalar-valued or vector-valued function representing the field.

    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances

    Examples
    ========

    >>> from sympy import Function
    >>> from sympy.vector import CoordSys3D, gradient

    The gradient of a scalar field is a vector field:

    >>> R = CoordSys3D('R')
    >>> s1 = R.x*R.y*R.z
    >>> gradient(s1)
    R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> s2 = 5*R.x**2*R.z
    >>> gradient(s2)
    10*R.x*R.z*R.i + 5*R.x**2*R.k

    The gradient of a vector field is a dyadic (second order tensor):

    >>> C = CoordSys3D("C", transformation="cylindrical")
    >>> e_r, e_theta, e_z = C.base_vectors()
    >>> r, theta, z = C.base_scalars()
    >>> v_r = Function('v_r')(r, theta, z)
    >>> v_theta = Function('v_theta')(r, theta, z)
    >>> v_z = Function('v_z')(r, theta, z)
    >>> v = v_r * e_r + v_theta * e_theta + v_z * e_z
    >>> gradient(v)     # doctest: +NORMALIZE_WHITESPACE
    (Derivative(v_r(C.r, C.theta, C.z), C.r))*(C.i|C.i)
    + (Derivative(v_theta(C.r, C.theta, C.z), C.r))*(C.i|C.j)
    + (Derivative(v_z(C.r, C.theta, C.z), C.r))*(C.i|C.k)
    + (-v_theta(C.r, C.theta, C.z)/C.r + Derivative(v_r(C.r, C.theta, C.z), C.theta)/C.r)*(C.j|C.i)
    + (v_r(C.r, C.theta, C.z)/C.r + Derivative(v_theta(C.r, C.theta, C.z), C.theta)/C.r)*(C.j|C.j)
    + (Derivative(v_z(C.r, C.theta, C.z), C.theta)/C.r)*(C.j|C.k)
    + (Derivative(v_r(C.r, C.theta, C.z), C.z))*(C.k|C.i)
    + (Derivative(v_theta(C.r, C.theta, C.z), C.z))*(C.k|C.j)
    + (Derivative(v_z(C.r, C.theta, C.z), C.z))*(C.k|C.k)

    """
    coord_sys = _get_coord_systems(field)
    if len(coord_sys) == 0:
        return Vector.zero

    if not isinstance(field, (Vector, Dyadic)):
        # scalar field
        if len(coord_sys) > 1:
            if isinstance(field, Add):
                return VectorAdd.fromiter(gradient(i) for i in field.args)
            if isinstance(field, Mul):
                s = _split_mul_args_wrt_coordsys(field)
                return VectorAdd.fromiter(field / i * gradient(i) for i in s)
            return Gradient(field)

        coord_sys = next(iter(coord_sys))
        h1, h2, h3 = coord_sys.lame_coefficients()
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        vx = Derivative(field, x) / h1
        vy = Derivative(field, y) / h2
        vz = Derivative(field, z) / h3

        if doit:
            return (vx * i + vy * j + vz * k).doit()
        return vx * i + vy * j + vz * k

    elif isinstance(field, Vector):
        # vector field
        if len(coord_sys) > 1:
            if isinstance(field, VectorAdd):
                return DyadicAdd.fromiter(gradient(i) for i in field.args)
            if isinstance(field, VectorMul):
                s = _split_mul_args_wrt_coordsys(field)
                return DyadicAdd.fromiter(field / i * gradient(i) for i in s)
            return Gradient(field)

        coord_sys = next(iter(coord_sys))
        e = coord_sys.base_vectors()
        q = coord_sys.base_scalars()
        h = coord_sys.lame_coefficients()
        A = [field & u for u in e]
        s = Dyadic.zero

        for i in range(3):
            for j in range(3):
                component = Derivative(A[j], q[i]) / h[i] - \
                    A[j] / (h[i] * h[j]) * Derivative(h[j], q[i])

                for k in range(3):
                    Gamma = _christoffel_symbol_2nd_kind(h, q, i, k, j)
                    component += h[j] / (h[i] * h[k]) * A[k] * Gamma

                s += component * (e[i] | e[j])
        if doit:
            return s.doit()
        return s

    return Gradient(field)


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
