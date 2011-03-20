from sympy import *
from volumetricEntity import Geometry3dEntity

class Point3d(Geometry3dEntity):
 """
    A point in Euclidean N-space defined by a sequence of values. Can be
    constructed from a sequence of points.

    Examples:
    ======
        >>> from sympy.sterometry import Point3d
        >>> Point(10,10,20)
        Point(10, 10, 20)
 """
 def __new__(cls, *args, **kwargs):
        if isinstance(args[0], (tuple, list, set)):
            coords = tuple([sympify(x) for x in args[0]])
        else:
            coords = tuple([sympify(x) for x in args])

        if len(coords) != 3:
            raise NotImplementedError("Only three dimensional points currently supported")

        return Geometry3dEntity.__new__(cls, *coords)
