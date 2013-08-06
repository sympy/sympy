##############################################################################
from sympy import diff, integrate, S
from sympy.physics.mechanics import MovingRefFrame


def is_conservative(field):
    """
    Check if a field is conservative

    Paramaters
    ==========

    field : Vector
        The field to check for conservative property

    Examples
    ========

    >>> from sympy.physics.mechanics import MovingRefFrame, is_conservative
    >>> R = MovingRefFrame('R')
    >>> is_conservative(R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z)
    True
    >>> is_conservative(R[2] * R.y)
    False

    """

    #Field is conservative irrespective of frame
    #Take the first frame in the result of the
    #separate() method
    frame = field.separate().keys()[0]
    return field.curl(frame) == S(0)


def is_solenoidal(field):
    """
    Check if a field is solenoidal

    Paramaters
    ==========

    field : Vector
        The field to check for solenoidal property

    Examples
    ========

    >>> from sympy.physics.mechanics import MovingRefFrame, is_solenoidal
    >>> R = MovingRefFrame('R')
    >>> is_solenoidal(R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z)
    True
    >>> is_solenoidal(R[1] * R.y)
    False

    """

    #Field is solenoidal irrespective of frame
    #Take the first frame in the result of the
    #separate() method
    frame = field.separate().keys()[0]
    return field.divergence(frame) == S(0)


def scalar_potential(field, frame):
    """
    Returns the scalar potential function of a field in a given frame
    (without the added integration constant)

    Parameters
    ==========

    field : Vector
        The vector field whose scalar potential function is to be
        calculated

    frame : MovingRefFrame
        The frame to do the calculation in

    Examples
    ========

    """

    #Check whether field is conservative
    if not is_conservative(field):
        raise ValueError("Field is not conservative")
    #Express the field exntirely in frame
    #Susbitute coordinate variables also
    field = field.express(frame, variables=True)
    #Make a list of dimensions of the frame
    dims = frame.base_vectors
    #Calculate scalar potential function
    temp_function = integrate(field.dot(dimensions[0]), frame[0])
    for i, dim in enumerate(dimensions[1:]):
        partial_diff = diff(temp_function, frame[i+1])
        partial_diff = field.dot(dim) - partial_diff
        temp_function += integrate(partial_diff, frame[i+1])
    return temp_function


def scalar_potential_difference(field, frame, position1, position2):
    """
    Calculate the scalar potential difference between two points
    in a certain frame, wrt a given field.

    If a scalar field is provided, its values at the two points are
    considered. If a conservative vector field is provided, the values
    of its scalar potential function at the two points are used.

    Returns potential of position 2 - potential of position 1

    Parameters
    ==========

    field : Vector/scalar field
        The field to calculate wrt

    frame : MovingRefFrame
        The frame to do the calculations in

    position1 : Vector
        The position vector of the initial point in given frame

    position2 : Vector
        The position vector of the final point in given frame

    Examples
    ========
    
    """

    if isinstance(field, Vector):
        #Get the scalar potential function
        scalar_fn = scalar_potential(field, frame)
    else:
        #Field is a scalar
        scalar_fn = field
    #Express positions in required frame
    position1 = position1.express(frame, variables=True)
    position2 = position2.express(frame, variables=True)
    #Get the two positions as substitution dicts for coordinate variables
    subs_dict1 = {}
    subs_dict2 = {}
    for i, x in enumerate(frame):
        subs_dict1[frame[i]] = position1.dot(x)
        subs_dict2[frame[i]] = position2.dot(x)
    return scalar_fn.subs(subs_dict2) - scalar_fn.subs(subs_dict1)
