from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.core import Basic
from sympy.vector.vector import Vector
from sympy.vector.deloperator import Del


class Gradient(Basic):
    """
    Represents unevaluated Gradient.
    
    Examples
    ========
    
    >>> from sympy.vector import CoordSysCartesian, Gradient
    >>> R = CoordSysCartesian('R')
    >>> s1 = R.x*R.y*R.z
    >>> Gradient(s1, R)
    (Derivative(R.x*R.y*R.z, R.x))*R.i git + (Derivative(R.x*R.y*R.z, R.y))*R.j + (Derivative(R.x*R.y*R.z, R.z))*R.k
    """

    def __new__(cls, scalar, coord_sys):
        return Del(coord_sys).gradient(scalar)


class Divergence(Basic):
    """
    Represents unevaluated Divergence.

    Examples
    ========
    
    >>> from sympy.vector import CoordSysCartesian, Gradient
    >>> R = CoordSysCartesian('R')
    >>> v = R.x*R.y*R.z*R.i + R.j + R.k
    >>> Divergence(v, R)
    Derivative(R.x*R.y*R.z, R.x)
    
    """

    def __new__(cls, vect, coord_sys):
        return Del(coord_sys).dot(vect)


class Curl(Basic):
    """
    Represents unevaluated Curl.

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian, Curl
    >>> R = CoordSysCartesian('R')
    >>> s1 = R.x*R.y*R.z*R.i + R.j
    Curl(s1, R)
    (Derivative(0, R.y) - Derivative(1, R.z))*R.i + (-Derivative(0, R.x) + 
    Derivative(R.x*R.y*R.z, R.z))*R.j + (Derivative(1, R.x) - Derivative(R.x*R.y*R.z, R.y))*R.k    
    
    """

    def __new__(cls, vect, coord_sys):
        return Del(coord_sys).cross(vect)


def gradient(scalar, coord_sys):
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

    return Gradient(scalar, coord_sys).doit()


def divergence(vect, coord_sys):
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

    return Divergence(vect, coord_sys).doit()


def curl(vect, coord_sys):
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

    return Curl(vect, coord_sys).doit()
