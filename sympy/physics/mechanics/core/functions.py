from sympy import Symbol, diff, integrate, solve, sympify
from sympy.vector import BaseScalar, Vector, VectMul, VectAdd


def get_motion_pos(position=0, timevar, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the position vector as a function of time.

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========
    
    position : Vector/VectMul/VectAdd
        Position vector of an object/frame as a function of time

    timevar : Symbol
        The Symbol to be used as the time variable

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """
    
    if position != 0:
        if not position.is_vector:
            raise ValueError("position must be a vector")
    else:
        return [0, 0, 0]
    if not timevar.is_Symbol:
        raise ValueError("time variable must be a Symbol")
    position = frame.express(position)
    vel = diff(position, timevar)
    acc = diff(vel, timevar)
    return [acc, vel, position]


def get_motion_vel(velocity=0, position=0, timevar, timevalue=0, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the velocity and a boundary condition(position
    vector at time = timevalue).

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========

    velocity : Vector/VectMul/VectAdd
        Velocity vector of an object/frame as a function of time

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue

    timevar : Symbol
        The Symbol to be used as the time variable

    timevalue : Number
        The numeric value of time at which the given boundary condition
        has been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    if velocity != 0:
        if not velocity.is_vector:
            raise ValueError("velocity must be a vector")
    if position != 0:
        if not position.is_vector:
            raise ValueError("position must be a vector")
    if not timevar.is_Symbol:
        raise ValueError("time variable must be a Symbol")
    timevalue = sympify(timevalue)
    if timevar in timevalue.atoms():
        raise ValueError("timevalue must be independent of time")
    return _process_vector_differential(velocity, position, timevar, timevalue, frame)


def get_motion_acc(acceleration=0, velocity=0, position=0, timevar, timevalue1=0, timevalue2=0, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the acceleration and two boundary conditions
    (velocity vector at time = timevalue1 and position vector at time = timevalue2).

    Returns a list of [acceleration, velocity, position]

    Parameters
    ==========

    acceleration : Vector/VectMul/VectAdd
        Acceleration of the object/frame as a function of time

    velocity : Vector/VectMul/VectAdd
        Boundary condition of velocity at time = timevalue1

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue2

    timevar : Symbol
        The Symbol to be used as the time variable

    timevalue1, timevalue2 : Number
        The numeric values of times at which the given boundary conditions
        have been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    if acceleration != 0:
        if not velocity.is_vector:
            raise ValueError("acceleration must be a vector")
    if velocity != 0:
        if not position.is_vector:
            raise ValueError("velocity must be a vector")
    if position != 0:
        if not velocity.is_vector:
            raise ValueError("position must be a vector")
    if not timevar.is_Symbol:
        raise ValueError("time variable must be a Symbol")
    timevalue1 = sympify(timevalue1)
    timevalue2 = sympify(timevalue2)
    if timevar in timevalue1.atoms() or timevar in timevalue2.atoms():
        raise ValueError("time values must be independent of time")
    vel = _process_vector_differential(acceleration, velocity, timevar, timevalue1, frame)[2]
    pos = _process_vector_differential(vel, position, timevar, timevalue2, frame)[2]
    return [acceleration, vel, pos]


def _process_vector_differential(vectdiff, condition, variable, valueofvar, frame):
    """
    Helper function for get_motion methods

    Returns a list of - derivative, vectdiff and integral

    """

    #Make sure boundary condition is time independent and expressed entirely
    #in one frame
    if condition != 0:
        if type(condition) == Vector:
            if condition.coord_sys != frame:
                raise ValueError("Boundary condition must be defined in given frame")
        else:
            for x in condition.atoms():
                if type(x) == Vector or type(x) == BaseScalar:
                    if x.coord_sys != frame:
                        raise ValueError("Boundary condition must be defined in given frame")
                elif x == variable:
                    raise ValueError("Boundary condition must be time-independent")
    #Special case of vectdiff == 0
    if vectdiff == 0:
        return [0, 0, condition]
    #Express vectdiff completely in condition's frame
    vectdiff = frame.express(vectdiff)
    #Find derivative of vectdiff
    vectdiff2 = diff(vectdiff, timevar)
    #Integrate and use boundary condition
    vectdiff0 = 0
    for dim in frame.base_vectors:
        function1 = vectdiff.dot(dim)
        vectdiff0 += _integrate_boundary(function1, variable, valueofvar, dim.dot(condition))
    #Return list
    return [vectdiff2, vectdiff, vectdiff0]


def _integrate_boundary(expr, var, valueofvar, value):
    """
    Returns indefinite integratal of expr wrt var, using the boundary
    condition of expr's value being 'value' at var = valueofvar.
    """

    CoI = Symbol('CoI')
    expr = integrate(expr, var) + CoI
    return expr.subs({CoI : solve(expr.subs({var : valueofvar}) - value, CoI)[0]})
