from sympy.vector.coordsysrect import CoordSys3D
from sympy.integrals import Integral

class ParametricIntegral(Integral):
    """
    Represents integral of a scalar or vector field
    over a Parametric Region

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, ParametricRegion
    >>> from sympy.abc import t

    >>> C = CoordSys3D('C)
    >>> curve = ParametricRegion(t, (t, t**2), {t: (-1, 1)})
    >>> ParametricIntegral(C.x, curve)
    0
    >>> semisphere = ParametricRegion((phi, theta), (2*sin(phi)*cos(theta), 2*sin(phi)*sin(theta), 2*cos(phi))
                            {theta: (0, 2*pi), phi: (0, pi/2)})
    >>> ParametricIntegral(C.z, semisphere)
    8*pi

    >>> ParametricIntegral(C.j - C.k, ParametricRegion((r, theta), (r*cos(theta), r*sin(theta))))
    ParametricIntegral(C.j - C.k, ParametricRegion((r, theta), (r*cos(theta), r*sin(theta))))

    """

    def __int__(self, field, parametricregion):
        pass
