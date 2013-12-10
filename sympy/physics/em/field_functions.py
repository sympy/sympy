from sympy import diff, integrate, S
from sympy.physics.mechanics import ReferenceFrame, Vector, express
from sympy.physics.mechanics.essential import _check_frame, _check_vector


def curl(vect, frame):
    """
    The curl of a vector field in the given frame

    Parameters
    ==========

    vect : Vector
        The vector operand

    frame : ReferenceFrame
        The reference frame to calculate the curl in

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import curl
    >>> R = ReferenceFrame('R')
    >>> v1 = R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z
    >>> curl(v1, R)
    0
    >>> v2 = R[0]*R[1]*R[2]*R.x
    >>> curl(v2, R)
    R_x*R_y*R.y - R_x*R_z*R.z

    """

    _check_vector(vect)
    if vect == 0:
        return Vector(0)
    vectx = vect.dot(frame.x)
    vecty = vect.dot(frame.y)
    vectz = vect.dot(frame.z)
    outvec = Vector(0)
    outvec += (diff(vectz, frame[1]) - diff(vecty, frame[2])) * frame.x
    outvec += (diff(vectx, frame[2]) - diff(vectz, frame[0])) * frame.y
    outvec += (diff(vecty, frame[0]) - diff(vectx, frame[1])) * frame.z
    return outvec


def divergence(vect, frame):
    """
    The divergence of a vector field in the given frame

    Parameters
    ==========

    vect : Vector
        The vector operand

    frame : ReferenceFrame
        The reference frame to calculate the divergence in

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import divergence
    >>> R = ReferenceFrame('R')
    >>> v1 = R[0]*R[1]*R[2] * (R.x+R.y+R.z)
    >>> divergence(v1, R)
    R_x*R_y + R_x*R_z + R_y*R_z
    >>> v2 = 2*R[1]*R[2]*R.y
    >>> divergence(v2, R)
    2*R_z

    """

    _check_vector(vect)
    if vect == 0:
        return S(0)
    vectx = vect.dot(frame.x)
    vecty = vect.dot(frame.y)
    vectz = vect.dot(frame.z)
    out = S(0)
    out += diff(vectx, frame[0])
    out += diff(vecty, frame[1])
    out += diff(vectz, frame[2])
    return out


def gradient(scalar, frame):
    """
    The vector gradient of a scalar field in the given frame

    Parameters
    ==========

    scalar : sympyfiable
        The scalar field to take the gradient of

    frame : ReferenceFrame
        The frame to calculate the gradient in

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import gradient
    >>> R = ReferenceFrame('R')
    >>> s1 = R[0]*R[1]*R[2]
    >>> gradient(s1, R)
    R_y*R_z*R.x + R_x*R_z*R.y + R_x*R_y*R.z
    >>> s2 = 5*R[0]**2*R[2]
    >>> gradient(s2, R)
    10*R_x*R_z*R.x + 5*R_x**2*R.z

    """

    _check_frame(frame)
    outvec = Vector(0)
    for i, x in enumerate(frame):
        outvec += diff(scalar, frame[i]) * x
    return outvec


def is_conservative(field):
    """
    Checks if a field is conservative

    Paramaters
    ==========

    field : Vector
        The field to check for conservative property

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import is_conservative
    >>> R = ReferenceFrame('R')
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
    frame = _separate(field).keys()[0]
    return curl(field, frame) == 0


def is_solenoidal(field):
    """
    Checks if a field is solenoidal

    Paramaters
    ==========

    field : Vector
        The field to check for solenoidal property

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import is_solenoidal
    >>> R = ReferenceFrame('R')
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
    frame = _separate(field).keys()[0]
    return divergence(field, frame) == 0


def scalar_potential(field, frame):
    """
    The scalar potential function of a field in a given frame
    (without the added integration constant)

    Parameters
    ==========

    field : Vector
        The vector field whose scalar potential function is to be
        calculated

    frame : ReferenceFrame
        The frame to do the calculation in

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.em import scalar_potential, gradient
    >>> R = ReferenceFrame('R')
    >>> scalar_potential(R.z, R) == R[2]
    True
    >>> scalar_field = 2*R[0]**2*R[1]*R[2]
    >>> grad_field = gradient(scalar_field, R)
    >>> scalar_potential(grad_field, R)
    2*R_x**2*R_y*R_z

    """

    #Check whether field is conservative
    if not is_conservative(field):
        raise ValueError("Field is not conservative")
    if field == 0:
        return Vector(0)
    #Express the field exntirely in frame
    #Subsitute coordinate variables also
    _check_frame(frame)
    field = express(field, frame, variables=True)
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
    The scalar potential difference between two points in a certain
    frame, wrt a given field.

    If a scalar field is provided, its values at the two points are
    considered. If a conservative vector field is provided, the values
    of its scalar potential function at the two points are used.

    Returns (potential at position 2) - (potential at position 1)

    Parameters
    ==========

    field : Vector/sympyfiable
        The field to calculate wrt

    frame : ReferenceFrame
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
    position1 = express(position1, frame, variables=True)
    position2 = express(position2, frame, variables=True)
    #Get the two positions as substitution dicts for coordinate variables
    subs_dict1 = {}
    subs_dict2 = {}
    for i, x in enumerate(frame):
        subs_dict1[frame[i]] = x.dot(position1)
        subs_dict2[frame[i]] = x.dot(position2)
    return scalar_fn.subs(subs_dict2) - scalar_fn.subs(subs_dict1)


def _separate(vect):
    """
    The components of a vector in different reference frames,
    as per its definition

    Returns a dict mapping each frame to the corresponding
    component vector.
    """
    _check_vector(vect)
    if vect == 0:
        return {}
    components = {}
    for x in vect.args:
        components[x[1]] = Vector([x])
    return components
