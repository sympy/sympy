from sympy.core.expr import Expr
from sympy.core import  sympify, S
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.vector import Vector
from sympy.vector.scalar import BaseScalar
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.core.function import Derivative


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
        coord_sys = list(expr.atoms(CoordSys3D))
        if len(coord_sys) == 1:
            return coord_sys[0]
        else:
            return None
    except:
        return None


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

    if coord_sys is None:
        return Vector.zero
    else:
        h1, h2, h3 = coord_sys.lame_coefficients()
        from sympy.vector.functions import express
        vectx = express(vect.dot(coord_sys.i), coord_sys, variables=True)
        vecty = express(vect.dot(coord_sys.j), coord_sys, variables=True)
        vectz = express(vect.dot(coord_sys.k), coord_sys, variables=True)
        outvec = Vector.zero
        outvec += (Derivative(vectz * h3, coord_sys.y) -
                   Derivative(vecty * h2, coord_sys.z)) * coord_sys.i / (h2 * h3)
        outvec += (Derivative(vectx * h1, coord_sys.z) -
                   Derivative(vectz * h3, coord_sys.x)) * coord_sys.j / (h1 * h3)
        outvec += (Derivative(vecty * h2, coord_sys.x) -
                   Derivative(vectx * h1, coord_sys.y)) * coord_sys.k / (h2 * h1)

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

    if coord_sys is None:
        return S.Zero
    else:
        h1, h2, h3 = coord_sys.lame_coefficients()
        vx = _diff_conditional(vect.dot(coord_sys.i), coord_sys.x, h2, h3) \
             / (h1 * h2 * h3)
        vy = _diff_conditional(vect.dot(coord_sys.j), coord_sys.y, h3, h1) \
             / (h1 * h2 * h3)
        vz = _diff_conditional(vect.dot(coord_sys.k), coord_sys.z, h1, h2) \
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

    if coord_sys is None:
        return Vector.zero
    else:
        from sympy.vector.functions import express
        h1, h2, h3 = coord_sys.lame_coefficients()
        scalar_field = express(scalar_field, coord_sys,
                                   variables=True)
        vx = Derivative(scalar_field, coord_sys.x) / h1
        vy = Derivative(scalar_field, coord_sys.y) / h2
        vz = Derivative(scalar_field, coord_sys.z) / h3

        if doit:
            return (vx * coord_sys.i + vy * coord_sys.j + vz * coord_sys.k).doit()
        return vx * coord_sys.i + vy * coord_sys.j + vz * coord_sys.k


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
