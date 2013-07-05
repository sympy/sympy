from sympy import Symbol, diff, integrate, solve, sympify
from sympy.vector import BaseScalar, Vector, VectMul, VectAdd


def get_motion_pos(position=0, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the position vector as a function of time.

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========
    
    position : Vector/VectMul/VectAdd
        Position vector of an object/frame as a function of time

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """
    
    if position != 0:
        _check_vector(position)
    else:
        return [0, 0, 0]
    vel = frame.dt(position)
    acc = frame.dt(vel)
    return [acc, vel, position]


def get_motion_vel(velocity=0, position=0, timevalue=0, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the velocity and a boundary
    condition(position vector at time = timevalue).

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========

    velocity : Vector/VectMul/VectAdd
        Velocity vector of an object/frame as a function of time

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue

    timevalue : Number
        The numeric value of time at which the given boundary condition
        has been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    _check_vector(velocity)
    _check_vector(position)
    timevalue = sympify(timevalue)
    if frame.time in timevalue.atoms():
        raise ValueError("timevalue must be independent of time")
    return _process_vector_differential(velocity, position, frame.time,
                                        timevalue, frame)


def get_motion_acc(acceleration=0, velocity=0, position=0, timevalue1=0,
                   timevalue2=0, frame):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the acceleration and two boundary
    conditions-
    velocity vector at time = timevalue1 and position vector at time = timevalue2.

    Returns a list of [acceleration, velocity, position]

    Parameters
    ==========

    acceleration : Vector/VectMul/VectAdd
        Acceleration of the object/frame as a function of time

    velocity : Vector/VectMul/VectAdd
        Boundary condition of velocity at time = timevalue1

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue2

    timevalue1, timevalue2 : Number
        The numeric values of times at which the given boundary conditions
        have been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    _check_vector(acceleration)
    _check_vector(velocity)
    _check_vector(position)
    timevalue1 = sympify(timevalue1)
    timevalue2 = sympify(timevalue2)
    if frame.time in timevalue1.atoms() or frame.time in timevalue2.atoms():
        raise ValueError("time values must be independent of time")
    vel = _process_vector_differential(acceleration, velocity, frame.time,
                                       timevalue1, frame)[2]
    pos = _process_vector_differential(vel, position, frame.time,
                                       timevalue2, frame)[2]
    return [acceleration, vel, pos]


def _process_vector_differential(vectdiff, condition, variable, valueofvar, frame):
    """
    Helper function for get_motion methods

    Returns a list of - derivative, vectdiff and integral

    """

    #Make sure boundary condition is independent of 'variable'
    if condition != 0:
        condition = frame.express(condition)
        if variable in condition.atoms():
            raise ValueError("Boundary condition must be independent of " + \
                             str(variable))
    #Special case of vectdiff == 0
    if vectdiff == 0:
        return [0, 0, condition]
    #Express vectdiff completely in condition's frame to give vectdiff1
    vectdiff1 = frame.express(vectdiff)
    #Find derivative of vectdiff
    vectdiff2 = frame.dt(vectdiff)
    #Integrate and use boundary condition
    vectdiff0 = 0
    for dim in frame.base_vectors:
        function1 = vectdiff1.dot(dim)
        vectdiff0 += _integrate_boundary(function1, variable, valueofvar,
                                         dim.dot(condition))
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


def _check_vector(test_vect):
    """
    Helper to check whether an instance is a vector
    """

    test_vect = sympify(test_vect)
    if test_vect != 0:
        if not test_vect.is_Vector:
            raise TypeError(str(test_vect) + " should be a vector.")
