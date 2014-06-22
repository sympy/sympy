"""
**Contains**

* snellslaw
"""

from __future__ import division

__all__ = ['snellslaw']

from sympy import Symbol, sympify, sqrt, Matrix
from sympy.physics.units import c, u0, e0
from sympy.physics.optics import Medium


def snellslaw(incident, normal, medium1, medium2):
    """
    This function calculates transmitted vector after refraction.
    `medium1` and `medium2` can be `Medium` or any sympifiable object.

    Parameters
    ==========

    incident : Matrix, tuple or list
        Incident vector
    normal : Matrix, tuple or list
        Normal vector
    medium1 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 1 or its refractive index
    medium2 : sympy.physics.optics.medium.Medium or sympifiable
        Medium 2 or its refractive index

    Examples
    ========


    """
    if not isinstance(incident, Matrix):
        if type(incident) == type(()) or type(incident) == type([]):
            incident = Matrix(incident)
        else:
            raise TypeError("Incident vector should be a matrix, tuple or list")
    else:
        incident = incident

    if not isinstance(normal, Matrix):
        if type(normal) == type(()) or type(normal) == type([]):
            normal = Matrix(normal)
        else:
            raise TypeError("Normal vector should be a matrix, tuple or list")
    else:
        normal = normal

    m1, m2 = None, None
    n1, n2 = None, None

    # If m1 and m2 are instances of Medium, assign them to m1 and m2
    # respectively otherwise m1 and m2 will be treated as refrative
    # indices of the mediums and will be assigned to n1 and n2 respectively.

    if isinstance(medium1, Medium):
        m1 = medium1
        n1 = medium1.refractive_index
    else:
        n1 = sympify(medium1)

    if isinstance(medium2, Medium):
        m2 = medium2
        n2 = medium2.refractive_index
    else:
        n2 = sympify(medium2)

    eta = n1/n2  # Relative index of refraction
    c1 = -incident.dot(normal)  # cos(angle_of_incidence)
    cs2 = 1 - eta**2*(1 - c1**2)  # cos(angle_of_refraction)**2
    if cs2 < 0:  # This is the case of total internal refraction(TIR).
        return 0
    return eta*incident + (eta*c1 - sqrt(cs2))*normal
