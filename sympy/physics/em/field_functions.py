from sympy import diff, integrate, S
from sympy.physics.mechanics import MovingRefFrame, Vector
from sympy.physics.mechanics.core import _check_vector, _check_frame


def curl(vect, frame):
    """ The curl of a vector field in given frame """

    _check_vector(vect)
    if vect == 0:
        return 0
    vectx = vect.dot(frame.x)
    vecty = vect.dot(frame.y)
    vectz = vect.dot(frame.z)
    outvec = 0
    outvec += (diff(vectz, frame[1]) - diff(vecty, frame[2])) * frame.x
    outvec += (diff(vectx, frame[2]) - diff(vectz, frame[0])) * frame.y
    outvec += (diff(vecty, frame[0]) - diff(vectx, frame[1])) * frame.z
    return outvec


def divergence(vect, frame):
    """ The divergence of a vector field in given frame """

    _check_vector(vect)
    if vect == 0:
        return 0
    vectx = vect.dot(frame.x)
    vecty = vect.dot(frame.y)
    vectz = vect.dot(frame.z)
    out = 0
    out += diff(vectx, frame[0])
    out += diff(vecty, frame[1])
    out += diff(vectz, frame[2])
    return out


def separate(vect):
    
    _check_vector(vect)
    if vect == 0:
        return {}
    components = {}
    for x in vect.args:
        components[x[1]] = 0
        components[x[1]] += x[0][0] * x[1].x
        components[x[1]] += x[0][1] * x[1].y
        components[x[1]] += x[0][2] * x[1].z
    return components


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
    if field == 0:
        return True
    frame = separate(field).keys()[0]
    return curl(field, frame) == 0


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
    if field == 0:
        return True
    frame = separate(field).keys()[0]
    return divergence(field, frame) == 0


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

    >>> from sympy.physics.mechanics import MovingRefFrame
    >>> from sympy.physics.em import scalar_potential, gradient
    >>> scalar_potential(R.z, R) == R[2]
    True
    >>> scalar_field = 2*R[0]**2*R[1]*R[2]
    >>> grad_field = gradient(scalar_field, R)
    >>> scalar_potential(grad_field, R)
    2*R[0]**2*R[1]*R[2]

    """

    #Check whether field is conservative
    if not is_conservative(field):
        raise ValueError("Field is not conservative")
    if field == 0:
        return 0
    #Express the field exntirely in frame
    #Subsitute coordinate variables also
    _check_frame(frame)
    field = frame.express(field)
    #Make a list of dimensions of the frame
    dimensions = [x for x in frame]
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

    _check_frame(frame)
    if isinstance(field, Vector):
        #Get the scalar potential function
        scalar_fn = scalar_potential(field, frame)
    else:
        #Field is a scalar
        scalar_fn = field
    #Express positions in required frame
    position1 = frame.express(position1)
    position2 = frame.express(position2)
    #Get the two positions as substitution dicts for coordinate variables
    subs_dict1 = {}
    subs_dict2 = {}
    for i, x in enumerate(frame):
        subs_dict1[frame[i]] = x.dot(position1)
        subs_dict2[frame[i]] = x.dot(position2)
    return scalar_fn.subs(subs_dict2) - scalar_fn.subs(subs_dict1)

def gradient(scalar, frame):
    """
    The vector gradient of a scalar field in the given frame

    Parameters
    ==========

    scalar : sympyfiable
        The scalar field to take the gradient of

    frame : MovingRefFrame
        The frame to calculate the gradient in

    Examples
    ========

    """

    _check_frame(frame)
    outvec = 0
    for i, x in enumerate(frame.base_vectors):
        outvec += diff(scalar, frame[i]) * x
    return outvec
