import collections
from sympy.core.expr import Expr
from sympy.core import sympify, S, preorder_traversal
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.vector import Vector, VectorMul, VectorAdd, Dot
from sympy.vector.scalar import BaseScalar
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.core.function import Derivative
from sympy import Add, Mul


def _get_coord_sys_from_expr(expr, coord_sys=None):
    """
    expr : expression
        The coordinate system is extracted from this parameter.
    """
    if coord_sys is not None:
        SymPyDeprecationWarning(
            feature="coord_sys parameter",
            useinstead="do not use it",
            deprecated_since_version="1.1",
            issue=12884,
        ).warn()

    try:
        return expr.atoms(CoordSys3D)
    except:
        return set([])


def _get_coord_systems(expr):
    g = preorder_traversal(expr)
    ret = set([])
    for i in g:
        if isinstance(i, CoordSys3D):
            ret.add(i)
            g.skip()
    return ret


def _split_mul_args_wrt_coordsys(expr):
    d = collections.defaultdict(lambda: S.One)
    for i in expr.args:
        d[frozenset(_get_coord_systems(i))] *= i
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

    def doit(self, **kwargs):
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

    def doit(self, **kwargs):
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

    def doit(self, **kwargs):
        return curl(self._expr, doit=True)


def curl(vect, coord_sys=None, doit=True):
    """
    Returns the curl of a vector field computed wrt the base scalars
    of the given coordinate system.

    Parameters
    ==========

    vect : Vector
        The vector operand

    coord_sys : CoordSys3D
        The coordinate system to calculate the gradient in.
        Deprecated since version 1.1

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

    coord_sys = _get_coord_sys_from_expr(vect, coord_sys)

    if len(coord_sys) == 0:
        return Vector.zero
    elif len(coord_sys) == 1:
        coord_sys = coord_sys.pop()
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
        # TODO: use some of the vector calculus properties for this:
        coord_sys = coord_sys.pop()  # get one random coord_sys
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        h1, h2, h3 = coord_sys.lame_coefficients()

        from .functions import express
        vectx = express(vect.dot(i), coord_sys, variables=True)
        vecty = express(vect.dot(j), coord_sys, variables=True)
        vectz = express(vect.dot(k), coord_sys, variables=True)

        # This is a repetition of previous code, it will be removed as soon as
        # we have some better algorithm to deal with this case:
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


def divergence(vect, coord_sys=None, doit=True):
    """
    Returns the divergence of a vector field computed wrt the base
    scalars of the given coordinate system.

    Parameters
    ==========

    vector : Vector
        The vector operand

    coord_sys : CoordSys3D
        The coordinate system to calculate the gradient in
        Deprecated since version 1.1

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
    coord_sys = _get_coord_sys_from_expr(vect, coord_sys)

    if len(coord_sys) == 0:
        return S.Zero
    else:
        # TODO: is case of many coord systems, this gets a random one:
        coord_sys = coord_sys.pop()
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        h1, h2, h3 = coord_sys.lame_coefficients()
        vx = _diff_conditional(vect.dot(i), x, h2, h3) \
             / (h1 * h2 * h3)
        vy = _diff_conditional(vect.dot(j), y, h3, h1) \
             / (h1 * h2 * h3)
        vz = _diff_conditional(vect.dot(k), z, h1, h2) \
             / (h1 * h2 * h3)
        if doit:
            return (vx + vy + vz).doit()
        return vx + vy + vz


def gradient(scalar_field, coord_sys=None, doit=True):
    """
    Returns the vector gradient of a scalar field computed wrt the
    base scalars of the given coordinate system.

    Parameters
    ==========

    scalar_field : SymPy Expr
        The scalar field to compute the gradient of

    coord_sys : CoordSys3D
        The coordinate system to calculate the gradient in
        Deprecated since version 1.1

    doit : bool
        If True, the result is returned after calling .doit() on
        each component. Else, the returned expression contains
        Derivative instances

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, gradient
    >>> R = CoordSys3D('R')
    >>> s1 = R.x*R.y*R.z
    >>> gradient(s1)
    R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> s2 = 5*R.x**2*R.z
    >>> gradient(s2)
    10*R.x*R.z*R.i + 5*R.x**2*R.k

    """
    coord_sys = _get_coord_sys_from_expr(scalar_field, coord_sys)

    if len(coord_sys) == 0:
        return Vector.zero
    elif len(coord_sys) == 1:
        coord_sys = coord_sys.pop()
        h1, h2, h3 = coord_sys.lame_coefficients()
        i, j, k = coord_sys.base_vectors()
        x, y, z = coord_sys.base_scalars()
        vx = Derivative(scalar_field, x) / h1
        vy = Derivative(scalar_field, y) / h2
        vz = Derivative(scalar_field, z) / h3

        if doit:
            return (vx * i + vy * j + vz * k).doit()
        return vx * i + vy * j + vz * k
    else:
        if isinstance(scalar_field, (Add, VectorAdd)):
            return VectorAdd.fromiter(gradient(i) for i in scalar_field.args)
        if isinstance(scalar_field, (Mul, VectorMul)):
            s = _split_mul_args_wrt_coordsys(scalar_field)
            return VectorAdd.fromiter(scalar_field / i * gradient(i) for i in s)
        return Gradient(scalar_field)


def _diff_conditional(expr, base_scalar, coeff_1, coeff_2):
    """
    First re-expresses expr in the system that base_scalar belongs to.
    If base_scalar appears in the re-expressed form, differentiates
    it wrt base_scalar.
    Else, returns S(0)
    """
    from sympy.vector.functions import express
    new_expr = express(expr, base_scalar.system, variables=True)
    if base_scalar in new_expr.atoms(BaseScalar):
        return Derivative(coeff_1 * coeff_2 * new_expr, base_scalar)
    return S(0)
