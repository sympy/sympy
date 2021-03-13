"""
Hyperbolic Geometrical Entity.
Contains
* Hyperbola
"""
from sympy import Point, Line
from sympy.geometry.entity import GeometrySet,GeometryEntity
class Hyperbola(GeometrySet):
    """
    A Hyperbolic GeometryEntity:
    A Hyperbola is the locus of a point in a plane which moves in such a way that the ratio of it's distance from
    a fixed point in the same plane to its distance from a fixed line is always constant which is always greater
    than unity .

    Parameters
    ==========
        focus=Point
        directrix=Line
        eccentricity: Eccentricity of Hyperbola is always one(>1).

    Attributes
    ==========
    focus
    directrix
    focal length
    eccentricity
    axis of symmetry(Tranverse and conjugate axis)
    Centre
    focal chord
    """

    def __new__(cls, focus=None,directrix=None, **kwargs):
        if focus:
            focus=Point(focus,dim=2)
        directrix=Line(directrix)
        if directrix.contains(focus):
            raise ValueError("The focus must not be on the directrix")
        return GeometryEntity.__new__(cls,focus,directrix,**kwargs)





