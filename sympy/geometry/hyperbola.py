"""
Hyperbolic Geometrical Entity.
Contains
* Hyperbola
"""
from sympy.geometry.entity import GeometrySet
class Hyperbola(GeometrySet):
    """
    A Hyperbolic GeometryEntity:
    A Hyperbola is the locus of a point in a plane which moves in such a way that the ratio of it's distance from
    a fixed point in the same plane to its distance from a fixed line is always constant which is always greater
    than unity .

    Parameters
    ==========
        eccentricity: Eccentricity of Hyperbola is always one(>1).

    """

