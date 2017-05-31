from sympy.core import Expr
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.vector import Vector
from sympy.vector.deloperator import Del


def _get_coord_sys_from_expr(expr):
    return list(expr.atoms(CoordSysCartesian))[0]


class Gradient(Expr):
    """
    Represents unevaluated Gradient.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Gradient
    >>> R = CoordSysCartesian('R')
    >>> s = R.x*R.y*R.z
    >>> Gradient(s, R)
    Gradient(R.x*R.y*R.z, R)

    """

    def __new__(cls, scalar):
        obj = Expr.__new__(cls, scalar)
        obj._scalar = scalar
        return obj

    def gradient(self):
        return Del(_get_coord_sys_from_expr(self._scalar)).gradient(self._scalar)

    def doit(self):
        return self.gradient().doit()


class Divergence(Expr):
    """
    Represents unevaluated Divergence.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Gradient
    >>> R = CoordSysCartesian('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Divergence(v, R)
    Divergence(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k, R)

    """

    def __new__(cls, vect):
        obj = Expr.__new__(cls, vect)
        obj._vect = vect
        return obj

    def divergence(self):
        return Del(_get_coord_sys_from_expr(self._vect)).dot(self._vect)

    def doit(self):
        return self.divergence().doit()


class Curl(Expr):
    """
    Represents unevaluated Curl.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Curl
    >>> R = CoordSysCartesian('R')
    >>> v = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> Curl(v, R)
    Curl(R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k, R)

    """

    def __new__(cls, vect):
        obj = Expr.__new__(cls, vect)
        obj._vect = vect
        return obj

    def curl(self):
        return Del(_get_coord_sys_from_expr(self._vect)).cross(self._vect)

    def doit(self):
        return self.curl().doit()


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

    return _get_coord_sys_from_expr(scalar).delop(scalar).doit()


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

    return _get_coord_sys_from_expr(vect).delop.dot(vect).doit()


def curl(vect, coord_sys=None):
    """
    Returns the curl of a vector field computed wrt the base scalars
    of the given coordinate system.

    Parameters
    ==========

    vect : Vector
        The vector operand

    coord_sys : CoordSysCartesian
        The coordinate system to calculate the gradient in

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

    return _get_coord_sys_from_expr(vect).delop.cross(vect).doit()
