from sympy.core import Expr
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.vector import Vector
from sympy.vector.deloperator import Del


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

    def __new__(cls, scalar, coord_sys):
        obj = Expr.__new__(cls, scalar, coord_sys)
        obj._scalar = scalar
        obj._coord_sys = coord_sys
        return obj

    def doit(self):
        return Del(self._coord_sys).gradient(self._scalar).doit()


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

    def __new__(cls, vect, coord_sys):
        obj = Expr.__new__(cls, vect, coord_sys)
        obj._vect = vect
        obj._coord_sys = coord_sys
        return obj

    def doit(self):
        return Del(self._coord_sys).dot(self._vect).doit()


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

    def __new__(cls, vect, coord_sys):
        obj = Expr.__new__(cls, vect, coord_sys)
        obj._vect = vect
        obj._coord_sys = coord_sys
        return obj

    def doit(self):
        return Del(self._coord_sys).cross(self._vect).doit()


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
    >>> gradient(s1, R)
    R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    >>> s2 = 5*R.x**2*R.z
    >>> gradient(s2, R)
    10*R.x*R.z*R.i + 5*R.x**2*R.k

    """

    if coord_sys is None:
        coord_sys = list(scalar.atoms(CoordSysCartesian))[0]
        return Gradient(scalar, coord_sys).doit()
    else:
        return Gradient(scalar, coord_sys).doit()


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
    >>> divergence(v2, R)
    2*R.z

    """

    if coord_sys is None:
        coord_sys = vect._sys
        return Curl(vect, coord_sys).doit()
    else:
        return Curl(vect, coord_sys).doit()


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
    >>> curl(v2, R)
    R.x*R.y*R.j + (-R.x*R.z)*R.k

    """

    if coord_sys is None:
        coord_sys = vect._sys
        return Curl(vect, coord_sys).doit()
    else:
        return Curl(vect, coord_sys).doit()
