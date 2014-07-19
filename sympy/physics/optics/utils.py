"""
**Contains**

* refraction_angle
"""

from __future__ import division

__all__ = ['refraction_angle', 'deviation']

from sympy import Symbol, sympify, sqrt, Matrix, acos
from sympy.geometry.line3d import Ray3D
from sympy.geometry.util import intersection
from sympy.geometry.plane import Plane
from .medium import Medium
from .waves import TWave


def refraction_angle(incident, medium1, medium2, normal=None, plane=None):
    """
    This function calculates transmitted vector after refraction.
    `medium1` and `medium2` can be `Medium` or any sympifiable object.

    If `incident` is an object of `Ray3D`, `normal` also has to be an instance
    of `Ray3D` in order to get the output as a `Ray3D`. Please note that if
    plane of separation is not provided and normal is an instance of `Ray3D`,
    normal will be assumed to be intersecting incident ray at the plane of
    separation. This will not be the case when `normal` is a `Matrix` or
    any other sequence.
    If `incident` is an instance of `Ray3D` and `plane` has not been provided
    and `normal` is not `Ray3D`, output will be a `Matrix`.

    Parameters
    ==========

    incident : Matrix, Ray3D, tuple or list
        Incident vector
    medium1 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 1 or its refractive index
    medium2 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 2 or its refractive index
    normal : Matrix, Ray3D, tuple or list
        Normal vector
    plane : Plane
        Plane of separation of the two media.

    Examples
    ========

    >>> from sympy.physics.optics import refraction_angle
    >>> from sympy.geometry import Point3D, Ray3D, Plane
    >>> from sympy.matrices import Matrix
    >>> from sympy import symbols
    >>> n = Matrix([0, 0, 1])
    >>> P = Plane(Point3D(0, 0, 0), normal_vector=[0, 0, 1])
    >>> r1 = Ray3D(Point3D(-1, -1, 1), Point3D(0, 0, 0))
    >>> refraction_angle(r1, 1, 1, n)
    Matrix([
    [ 1],
    [ 1],
    [-1]])
    >>> refraction_angle(r1, 1, 1, plane=P)
    Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))

    With different index of refraction of the two media

    >>> n1, n2 = symbols('n1, n2')
    >>> refraction_angle(r1, n1, n2, n)
    Matrix([
    [                                n1/n2],
    [                                n1/n2],
    [-sqrt(3)*sqrt(-2*n1**2/(3*n2**2) + 1)]])
    >>> refraction_angle(r1, n1, n2, plane=P)
    Ray3D(Point3D(0, 0, 0), Point3D(n1/n2, n1/n2, -sqrt(3)*sqrt(-2*n1**2/(3*n2**2) + 1)))

    """
    # A flag to check whether to return Ray3D or not
    return_ray = False

    if plane is not None and normal is not None:
        raise ValueError("Either plane or normal is acceptable.")

    if not isinstance(incident, Matrix):
        if type(incident) == type(()) or type(incident) == type([]):
            _incident = Matrix(incident)
        elif isinstance(incident, Ray3D):
            _incident = Matrix(incident.direction_ratio)
        else:
            raise TypeError("incident should be a Matrix, Ray3D, tuple or list")
    else:
        _incident = incident

    # If plane is provided, get direction ratios of the normal
    # to the plane from the plane else go with `normal` param.
    if plane is not None:
        if not isinstance(plane, Plane):
            raise TypeError("plane should be an instance of geometry.plane.Plane")
        # If we have the plane, we can get the intersection
        # point of incident ray and the plane and thus return
        # an instance of Ray3D.
        return_ray = True
        intersection_pt = plane.intersection(incident)
        _normal = plane.normal_vector
        if isinstance(_normal, list):
            _normal = Matrix(_normal)
    else:
        if not isinstance(normal, Matrix):
            if type(normal) == type(()) or type(normal) == type([]):
                _normal = Matrix(normal)
            elif isinstance(normal, Ray3D):
                _normal = Matrix(normal.direction_ratio)
                if isinstance(incident, Ray3D):
                    intersection_pt = intersection(incident, normal)
                    if len(intersection_pt) == 0:
                        raise ValueError("Normal isn't concurrent with the incident ray.")
                    else:
                        return_ray = True
                        intersection_pt = intersection_pt[0]
            else:
                raise TypeError("Normal should be a Matrix, Ray3D, tuple or list")
        else:
            _normal = normal

    n1, n2 = None, None

    if isinstance(medium1, Medium):
        n1 = medium1.refractive_index
    else:
        n1 = sympify(medium1)

    if isinstance(medium2, Medium):
        n2 = medium2.refractive_index
    else:
        n2 = sympify(medium2)

    eta = n1/n2  # Relative index of refraction
    # Calculating magnitude of the vectors
    mag_incident = sqrt(sum([i**2 for i in _incident]))
    mag_normal = sqrt(sum([i**2 for i in _normal]))
    # Converting vectors to unit vectors by dividing
    # them with their magnitudes
    _incident /= mag_incident
    _normal /= mag_normal
    c1 = -_incident.dot(_normal)  # cos(angle_of_incidence)
    cs2 = 1 - eta**2*(1 - c1**2)  # cos(angle_of_refraction)**2
    if cs2.is_negative:  # This is the case of total internal reflection(TIR).
        return 0
    drs = eta*_incident + (eta*c1 - sqrt(cs2))*_normal
    # Multiplying unit vector by its magnitude
    drs = drs*mag_incident
    if not return_ray:
        return drs
    else:
        return Ray3D(intersection_pt, direction_ratio=drs)


def deviation(incident, medium1, medium2, normal=None, plane=None):
    """
    This function calculates the angle of deviation of a ray
    due to refraction.

    Parameters
    ==========

    incident : Matrix, Ray3D, tuple or list
        Incident vector
    medium1 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 1 or its refractive index
    medium2 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 2 or its refractive index
    normal : Matrix, Ray3D, tuple or list
        Normal vector
    plane : Plane
        Plane of separation of the two media.

    Examples
    ========

    >>> from sympy.physics.optics import deviation
    >>> from sympy.geometry import Point3D, Ray3D, Plane
    >>> from sympy.matrices import Matrix
    >>> from sympy import symbols
    >>> n = Matrix([0, 0, 1])
    >>> P = Plane(Point3D(0, 0, 0), normal_vector=[0, 0, 1])
    >>> r1 = Ray3D(Point3D(-1, -1, 1), Point3D(0, 0, 0))
    >>> deviation(r1, 1, 1, n)
    0

    """
    refracted = refraction_angle(incident,
                                 medium1,
                                 medium2,
                                 normal=normal,
                                 plane=plane)
    if refracted != 0:
        if isinstance(refracted, Ray3D):
            refracted = Matrix(refracted.direction_ratio)

        if not isinstance(incident, Matrix):
            if type(incident) == type(()) or type(incident) == type([]):
                _incident = Matrix(incident)
            elif isinstance(incident, Ray3D):
                _incident = Matrix(incident.direction_ratio)
            else:
                raise TypeError("incident should be a Matrix, Ray3D, tuple or list")
        else:
            _incident = incident

        if not isinstance(normal, Matrix):
            if type(normal) == type(()) or type(normal) == type([]):
                _normal = Matrix(normal)
            elif isinstance(normal, Ray3D):
                _normal = Matrix(normal.direction_ratio)
            else:
                raise TypeError("normal should be a Matrix, Ray3D, tuple or list")
        else:
            _normal = normal

        mag_incident = sqrt(sum([i**2 for i in _incident]))
        mag_normal = sqrt(sum([i**2 for i in _normal]))
        mag_refracted = sqrt(sum([i**2 for i in refracted]))
        _incident /= mag_incident
        _normal /= mag_normal
        refracted /= mag_refracted
        i = acos(_incident.dot(_normal))
        r = acos(refracted.dot(_normal))
        return i - r
