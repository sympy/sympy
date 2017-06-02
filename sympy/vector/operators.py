from sympy.core import Expr
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.vector import Vector
from sympy.utilities.exceptions import SymPyDeprecationWarning


def _get_coord_sys_from_expr(expr, coord_sys=None):
    """

    expr : vector or scalar
        From this parameter is obtained information about
        coordinate system.
    coord_sys : CoordSysCartesian
        Deprecated since version 1.1. This parameter is unnecessary.

    """
    if coord_sys is not None:
        SymPyDeprecationWarning(
            feature="coord_sys parameter",
            useinstead="do not use it",
            deprecated_since_version="1.1"
        ).warn()
    return list(expr.atoms(CoordSysCartesian))[0]


class Gradient(Expr):
    """
    Represents unevaluated Gradient.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Gradient
    >>> R = CoordSysCartesian('R')
    >>> s = R.x*R.y*R.z
    >>> Gradient(s)
    Gradient(R.x*R.y*R.z)

    """

    def __new__(cls, scalar):
        obj = Expr.__new__(cls, scalar)
        obj._scalar = scalar
        return obj

    def doit(self):
        return gradient(self._scalar)


class Divergence(Expr):
    """
    Represents unevaluated Divergence.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Divergence
    >>> R = CoordSysCartesian('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Divergence(v)
    Divergence(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k)

    """

    def __new__(cls, vect):
        obj = Expr.__new__(cls, vect)
        obj._vect = vect
        return obj

    def doit(self):
        return divergence(self._vect)


class Curl(Expr):
    """
    Represents unevaluated Curl.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Curl
    >>> R = CoordSysCartesian('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Curl(v)
    Curl(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k)

    """

    def __new__(cls, vect):
        obj = Expr.__new__(cls, vect)
        obj._vect = vect
        return obj

    def doit(self):
        return curl(self._vect)


def gradient(scalar, coord_sys=None):
    """
    Returns the vector gradient of a scalar field computed wrt the
    base scalars of the given coordinate system.

    Parameters
    ==========

    scalar : SymPy Expr
        The scalar field to compute the gradient of

    coord_sys : CoordSysCartesian
        The coordinate system to calculate the gradient in
        Deprecated since version 1.1

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, gradient
    >>> R = CoordSysCartesian('R')
    >>> s1 = R.x*R.y*R.z
    >>> gradient(s1)
    R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> s2 = 5*R.x**2*R.z
    >>> gradient(s2)
    10*R.x*R.z*R.i + 5*R.x**2*R.k

    """

    return _get_coord_sys_from_expr(scalar, coord_sys).delop(scalar).doit()


def divergence(vect, coord_sys=None):
    """
    Returns the divergence of a vector field computed wrt the base
    scalars of the given coordinate system.

    Parameters
    ==========

    vect : Vector
        The vector operand

    coord_sys : CoordSysCartesian
        The coordinate system to calculate the gradient in
        Deprecated since version 1.1

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, divergence
    >>> R = CoordSysCartesian('R')
    >>> v1 = R.x*R.y*R.z * (R.i+R.j+R.k)
    >>> divergence(v1)
    R.x*R.y + R.x*R.z + R.y*R.z
    >>> v2 = 2*R.y*R.z*R.j
    >>> divergence(v2)
    2*R.z

    """

    return _get_coord_sys_from_expr(vect, coord_sys).delop.dot(vect).doit()


def curl(vect, coord_sys=None):
    """
    Returns the curl of a vector field computed wrt the base scalars
    of the given coordinate system.

    Parameters
    ==========

    vect : Vector
        The vector operand

    coord_sys : CoordSysCartesian
        The coordinate system to calculate the gradient in.
        Deprecated since version 1.1

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, curl
    >>> R = CoordSysCartesian('R')
    >>> v1 = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> curl(v1)
    0
    >>> v2 = R.x*R.y*R.z*R.i
    >>> curl(v2)
    R.x*R.y*R.j + (-R.x*R.z)*R.k

    """

    return _get_coord_sys_from_expr(vect, coord_sys).delop.cross(vect).doit()
