from __future__ import print_function, division

__all__ = ['cross',
           'dot',
           'express',
           'time_derivative',
           'outer',
           'inertia',
           'mechanics_printing',
           'mprint',
           'msprint',
           'mpprint',
           'mlatex',
           'kinematic_equations',
           'inertia_of_point_mass',
           'get_motion_params',
           'partial_velocity',
           'linear_momentum',
           'angular_momentum',
           'kinetic_energy',
           'potential_energy',
           'Lagrangian']

from sympy.physics.mechanics.essential import (Vector, Dyadic, ReferenceFrame,
                                               CoordinateSym,
                                               MechanicsStrPrinter,
                                               MechanicsPrettyPrinter,
                                               MechanicsLatexPrinter,
                                               dynamicsymbols,
                                               _check_frame, _check_vector)
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.rigidbody import RigidBody
from sympy.physics.mechanics.point import Point
from sympy import sympify, solve, diff, sin, cos, Matrix, Symbol, integrate, \
     trigsimp
from sympy.core.basic import S


def cross(vec1, vec2):
    """Cross product convenience wrapper for Vector.cross(): \n"""
    if not isinstance(vec1, (Vector, Dyadic)):
        raise TypeError('Cross product is between two vectors')
    return vec1 ^ vec2
cross.__doc__ += Vector.cross.__doc__


def dot(vec1, vec2):
    """Dot product convenience wrapper for Vector.dot(): \n"""
    if not isinstance(vec1, (Vector, Dyadic)):
        raise TypeError('Dot product is between two vectors')
    return vec1 & vec2
dot.__doc__ += Vector.dot.__doc__


def express(expr, frame, frame2=None, variables=False):
    """
    Global function for 'express' functionality.

    Re-expresses a Vector, scalar(sympyfiable) or Dyadic in given frame.

    Refer to the local methods of Vector and Dyadic for details.
    If 'variables' is True, then the coordinate variables (CoordinateSym
    instances) of other frames present in the vector/scalar field or
    dyadic expression are also substituted in terms of the base scalars of
    this frame.

    Parameters
    ==========

    expr : Vector/Dyadic/scalar(sympyfiable)
        The expression to re-express in ReferenceFrame 'frame'

    frame: ReferenceFrame
        The reference frame to express expr in

    frame2 : ReferenceFrame
        The other frame required for re-expression(only for Dyadic expr)

    variables : boolean
        Specifies whether to substitute the coordinate variables present
        in expr, in terms of those of frame

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, outer, dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q = dynamicsymbols('q')
    >>> B = N.orientnew('B', 'Axis', [q, N.z])
    >>> d = outer(N.x, N.x)
    >>> from sympy.physics.mechanics import express
    >>> express(d, B, N)
    cos(q)*(B.x|N.x) - sin(q)*(B.y|N.x)
    >>> express(B.x, N)
    cos(q)*N.x + sin(q)*N.y
    >>> express(N[0], B, variables=True)
    B_x*cos(q(t)) - B_y*sin(q(t))

    """

    _check_frame(frame)

    if expr == 0:
        return S(0)

    if isinstance(expr, Vector):
        #Given expr is a Vector
        if variables:
            #If variables attribute is True, substitute
            #the coordinate variables in the Vector
            frame_list = [x[-1] for x in expr.args]
            subs_dict = {}
            for f in frame_list:
                subs_dict.update(f.variable_map(frame))
            expr = expr.subs(subs_dict)
        #Re-express in this frame
        outvec = Vector([])
        for i, v in enumerate(expr.args):
            if v[1] != frame:
                temp = frame.dcm(v[1]) * v[0]
                if Vector.simp:
                    temp = temp.applyfunc(lambda x: \
                                          trigsimp(x, method='fu'))
                outvec += Vector([(temp, frame)])
            else:
                outvec += Vector([v])
        return outvec

    if isinstance(expr, Dyadic):
        if frame2 is None:
            frame2 = frame
        _check_frame(frame2)
        ol = Dyadic(0)
        for i, v in enumerate(expr.args):
            ol += express(v[0], frame, variables=variables) * \
                  (express(v[1], frame, variables=variables) | \
                   express(v[2], frame2, variables=variables))
        return ol

    else:
        if variables:
            #Given expr is a scalar field
            frame_set = set([])
            expr = sympify(expr)
            #Subsitute all the coordinate variables
            for x in expr.atoms():
                if isinstance(x, CoordinateSym)and x.frame != frame:
                    frame_set.add(x.frame)
            subs_dict = {}
            for f in frame_set:
                subs_dict.update(f.variable_map(frame))
            return expr.subs(subs_dict)
        return expr


def time_derivative(expr, frame, order=1):
    """
    Calculate the time derivative of a vector/scalar field function
    or dyadic expression in given frame.

    References
    ==========

    http://en.wikipedia.org/wiki/Rotating_reference_frame#Time_derivatives_in_the_two_frames

    Parameters
    ==========

    expr : Vector/Dyadic/sympifyable
        The expression whose time derivative is to be calculated

    frame : ReferenceFrame
        The reference frame to calculate the time derivative in

    order : integer
        The order of the derivative to be calculated

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, Vector, dynamicsymbols
    >>> from sympy import Symbol
    >>> q1 = Symbol('q1')
    >>> u1 = dynamicsymbols('u1')
    >>> N = ReferenceFrame('N')
    >>> A = N.orientnew('A', 'Axis', [q1, N.x])
    >>> v = u1 * N.x
    >>> A.set_ang_vel(N, 10*A.x)
    >>> from sympy.physics.mechanics import time_derivative
    >>> time_derivative(v, N)
    u1'*N.x
    >>> time_derivative(u1*A[0], N)
    N_x*Derivative(u1(t), t)
    >>> B = N.orientnew('B', 'Axis', [u1, N.z])
    >>> from sympy.physics.mechanics import outer
    >>> d = outer(N.x, N.x)
    >>> time_derivative(d, B)
    - u1'*(N.y|N.x) - u1'*(N.x|N.y)

    """

    t = dynamicsymbols._t
    _check_frame(frame)

    if order == 0:
        return expr
    if order%1 != 0 or order < 0:
        raise ValueError("Unsupported value of order entered")

    if isinstance(expr, Vector):
        outvec = Vector(0)
        for i, v in enumerate(expr.args):
            if v[1] == frame:
                outvec += Vector([(express(v[0], frame, \
                                           variables=True).diff(t), frame)])
            else:
                outvec += time_derivative(Vector([v]), v[1]) + \
                          (v[1].ang_vel_in(frame) ^ Vector([v]))
        return time_derivative(outvec, frame, order - 1)

    if isinstance(expr, Dyadic):
        ol = Dyadic(0)
        for i, v in enumerate(expr.args):
            ol += (v[0].diff(t) * (v[1] | v[2]))
            ol += (v[0] * (time_derivative(v[1], frame) | v[2]))
            ol += (v[0] * (v[1] | time_derivative(v[2], frame)))
        return time_derivative(ol, frame, order - 1)

    else:
        return diff(express(expr, frame, variables=True), t, order)


def outer(vec1, vec2):
    """Outer product convenience wrapper for Vector.outer():\n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Outer product is between two Vectors')
    return vec1 | vec2
outer.__doc__ += Vector.outer.__doc__


def inertia(frame, ixx, iyy, izz, ixy=0, iyz=0, izx=0):
    """Simple way to create inertia Dyadic object.

    If you don't know what a Dyadic is, just treat this like the inertia
    tensor.  Then, do the easy thing and define it in a body-fixed frame.

    Parameters
    ==========

    frame : ReferenceFrame
        The frame the inertia is defined in
    ixx : Sympifyable
        the xx element in the inertia dyadic
    iyy : Sympifyable
        the yy element in the inertia dyadic
    izz : Sympifyable
        the zz element in the inertia dyadic
    ixy : Sympifyable
        the xy element in the inertia dyadic
    iyz : Sympifyable
        the yz element in the inertia dyadic
    izx : Sympifyable
        the zx element in the inertia dyadic

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, inertia
    >>> N = ReferenceFrame('N')
    >>> inertia(N, 1, 2, 3)
    (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Need to define the inertia in a frame')
    ol = sympify(ixx) * (frame.x | frame.x)
    ol += sympify(ixy) * (frame.x | frame.y)
    ol += sympify(izx) * (frame.x | frame.z)
    ol += sympify(ixy) * (frame.y | frame.x)
    ol += sympify(iyy) * (frame.y | frame.y)
    ol += sympify(iyz) * (frame.y | frame.z)
    ol += sympify(izx) * (frame.z | frame.x)
    ol += sympify(iyz) * (frame.z | frame.y)
    ol += sympify(izz) * (frame.z | frame.z)
    return ol


def inertia_of_point_mass(mass, pos_vec, frame):
    """Inertia dyadic of a point mass realtive to point O.

    Parameters
    ==========

    mass : Sympifyable
        Mass of the point mass
    pos_vec : Vector
        Position from point O to point mass
    frame : ReferenceFrame
        Reference frame to express the dyadic in

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import ReferenceFrame, inertia_of_point_mass
    >>> N = ReferenceFrame('N')
    >>> r, m = symbols('r m')
    >>> px = r * N.x
    >>> inertia_of_point_mass(m, px, N)
    m*r**2*(N.y|N.y) + m*r**2*(N.z|N.z)

    """

    return mass * (((frame.x | frame.x) + (frame.y | frame.y) +
                   (frame.z | frame.z)) * (pos_vec & pos_vec) -
                   (pos_vec | pos_vec))


def mechanics_printing():
    """Sets up interactive printing for mechanics' derivatives.

    The main benefit of this is for printing of time derivatives;
    instead of displaying as Derivative(f(t),t), it will display f'
    This is only actually needed for when derivatives are present and are not
    in a physics.mechanics object.

    Examples
    ========

    >>> # 2 lines below are for tests to function properly
    >>> import sys
    >>> sys.displayhook = sys.__displayhook__
    >>> from sympy import Function, Symbol, diff
    >>> from sympy.physics.mechanics import mechanics_printing
    >>> f = Function('f')
    >>> t = Symbol('t')
    >>> x = Symbol('x')
    >>> diff(f(t), t)
    Derivative(f(t), t)
    >>> mechanics_printing()
    >>> diff(f(t), t)
    f'
    >>> diff(f(x), x)
    Derivative(f(x), x)
    >>> # 2 lines below are for tests to function properly
    >>> import sys
    >>> sys.displayhook = sys.__displayhook__

    """

    import sys
    sys.displayhook = mprint


def mprint(expr, **settings):
    r"""Function for printing of expressions generated in mechanics.

    Extends SymPy's StrPrinter; mprint is equivalent to:
    print sstr()
    mprint takes the same options as sstr.

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to print
    settings : args
        Same as print for SymPy

    Examples
    ========

    >>> from sympy.physics.mechanics import mprint, dynamicsymbols
    >>> u1 = dynamicsymbols('u1')
    >>> print(u1)
    u1(t)
    >>> mprint(u1)
    u1

    """

    outstr = msprint(expr, **settings)

    from sympy.core.compatibility import builtins
    if (outstr != 'None'):
        builtins._ = outstr
        print(outstr)


def msprint(expr, **settings):
    r"""Function for displaying expressions generated in mechanics.

    Returns the output of mprint() as a string.

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to print
    settings : args
        Same as print for SymPy

    Examples
    ========

    >>> from sympy.physics.mechanics import msprint, dynamicsymbols
    >>> u1, u2 = dynamicsymbols('u1 u2')
    >>> u2d = dynamicsymbols('u2', level=1)
    >>> print("%s = %s" % (u1, u2 + u2d))
    u1(t) = u2(t) + Derivative(u2(t), t)
    >>> print("%s = %s" % (msprint(u1), msprint(u2 + u2d)))
    u1 = u2 + u2'

    """

    pr = MechanicsStrPrinter(settings)
    return pr.doprint(expr)


def mpprint(expr, **settings):
    r"""Function for pretty printing of expressions generated in mechanics.

    Mainly used for expressions not inside a vector; the output of running
    scripts and generating equations of motion. Takes the same options as
    SymPy's pretty_print(); see that function for more information.

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to pretty print
    settings : args
        Same as pretty print

    Examples
    ========

    Use in the same way as pprint

    """

    mp = MechanicsPrettyPrinter(settings)
    print(mp.doprint(expr))


def mlatex(expr, **settings):
    r"""Function for printing latex representation of mechanics objects.

    For latex representation of Vectors, Dyadics, and dynamicsymbols. Takes the
    same options as SymPy's latex(); see that function for more information;

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to represent in LaTeX form
    settings : args
        Same as latex()

    Examples
    ========

    >>> from sympy.physics.mechanics import mlatex, ReferenceFrame, dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q1, q2 = dynamicsymbols('q1 q2')
    >>> q1d, q2d = dynamicsymbols('q1 q2', 1)
    >>> q1dd, q2dd = dynamicsymbols('q1 q2', 2)
    >>> mlatex(N.x + N.y)
    '\\mathbf{\\hat{n}_x} + \\mathbf{\\hat{n}_y}'
    >>> mlatex(q1 + q2)
    'q_{1} + q_{2}'
    >>> mlatex(q1d)
    '\\dot{q}_{1}'
    >>> mlatex(q1 * q2d)
    'q_{1} \\dot{q}_{2}'
    >>> mlatex(q1dd * q1 / q1d)
    '\\frac{q_{1} \\ddot{q}_{1}}{\\dot{q}_{1}}'

    """

    return MechanicsLatexPrinter(settings).doprint(expr)


def kinematic_equations(speeds, coords, rot_type, rot_order=''):
    """Gives equations relating the qdot's to u's for a rotation type.

    Supply rotation type and order as in orient. Speeds are assumed to be
    body-fixed; if we are defining the orientation of B in A using by rot_type,
    the angular velocity of B in A is assumed to be in the form: speed[0]*B.x +
    speed[1]*B.y + speed[2]*B.z

    Parameters
    ==========

    speeds : list of length 3
        The body fixed angular velocity measure numbers.
    coords : list of length 3 or 4
        The coordinates used to define the orientation of the two frames.
    rot_type : str
        The type of rotation used to create the equations. Body, Space, or
        Quaternion only
    rot_order : str
        If applicable, the order of a series of rotations.

    Examples
    ========

    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> from sympy.physics.mechanics import kinematic_equations, mprint
    >>> u1, u2, u3 = dynamicsymbols('u1 u2 u3')
    >>> q1, q2, q3 = dynamicsymbols('q1 q2 q3')
    >>> mprint(kinematic_equations([u1,u2,u3], [q1,q2,q3], 'body', '313'),
    ...     order=None)
    [-(u1*sin(q3) + u2*cos(q3))/sin(q2) + q1', -u1*cos(q3) + u2*sin(q3) + q2', (u1*sin(q3) + u2*cos(q3))*cos(q2)/sin(q2) - u3 + q3']

    """

    # Code below is checking and sanitizing input
    approved_orders = ('123', '231', '312', '132', '213', '321', '121', '131',
                       '212', '232', '313', '323', '1', '2', '3', '')
    rot_order = str(rot_order).upper()  # Now we need to make sure XYZ = 123
    rot_type = rot_type.upper()
    rot_order = [i.replace('X', '1') for i in rot_order]
    rot_order = [i.replace('Y', '2') for i in rot_order]
    rot_order = [i.replace('Z', '3') for i in rot_order]
    rot_order = ''.join(rot_order)

    if not isinstance(speeds, (list, tuple)):
        raise TypeError('Need to supply speeds in a list')
    if len(speeds) != 3:
        raise TypeError('Need to supply 3 body-fixed speeds')
    if not isinstance(coords, (list, tuple)):
        raise TypeError('Need to supply coordinates in a list')
    if rot_type.lower() in ['body', 'space']:
        if rot_order not in approved_orders:
            raise ValueError('Not an acceptable rotation order')
        if len(coords) != 3:
            raise ValueError('Need 3 coordinates for body or space')
        # Actual hard-coded kinematic differential equations
        q1, q2, q3 = coords
        q1d, q2d, q3d = [diff(i, dynamicsymbols._t) for i in coords]
        w1, w2, w3 = speeds
        s1, s2, s3 = [sin(q1), sin(q2), sin(q3)]
        c1, c2, c3 = [cos(q1), cos(q2), cos(q3)]
        if rot_type.lower() == 'body':
            if rot_order == '123':
                return [q1d - (w1 * c3 - w2 * s3) / c2, q2d - w1 * s3 - w2 *
                        c3, q3d - (-w1 * c3 + w2 * s3) * s2 / c2 - w3]
            if rot_order == '231':
                return [q1d - (w2 * c3 - w3 * s3) / c2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (- w2 * c3 + w3 * s3) * s2 / c2]
            if rot_order == '312':
                return [q1d - (-w1 * s3 + w3 * c3) / c2, q2d - w1 * c3 - w3 *
                        s3, q3d - (w1 * s3 - w3 * c3) * s2 / c2 - w2]
            if rot_order == '132':
                return [q1d - (w1 * c3 + w3 * s3) / c2, q2d + w1 * s3 - w3 *
                        c3, q3d - (w1 * c3 + w3 * s3) * s2 / c2 - w2]
            if rot_order == '213':
                return [q1d - (w1 * s3 + w2 * c3) / c2, q2d - w1 * c3 + w2 *
                        s3, q3d - (w1 * s3 + w2 * c3) * s2 / c2 - w3]
            if rot_order == '321':
                return [q1d - (w2 * s3 + w3 * c3) / c2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 - (w2 * s3 + w3 * c3) * s2 / c2]
            if rot_order == '121':
                return [q1d - (w2 * s3 + w3 * c3) / s2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 + (w2 * s3 + w3 * c3) * c2 / s2]
            if rot_order == '131':
                return [q1d - (-w2 * c3 + w3 * s3) / s2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (w2 * c3 - w3 * s3) * c2 / s2]
            if rot_order == '212':
                return [q1d - (w1 * s3 - w3 * c3) / s2, q2d - w1 * c3 - w3 *
                        s3, q3d - (-w1 * s3 + w3 * c3) * c2 / s2 - w2]
            if rot_order == '232':
                return [q1d - (w1 * c3 + w3 * s3) / s2, q2d + w1 * s3 - w3 *
                        c3, q3d + (w1 * c3 + w3 * s3) * c2 / s2 - w2]
            if rot_order == '313':
                return [q1d - (w1 * s3 + w2 * c3) / s2, q2d - w1 * c3 + w2 *
                        s3, q3d + (w1 * s3 + w2 * c3) * c2 / s2 - w3]
            if rot_order == '323':
                return [q1d - (-w1 * c3 + w2 * s3) / s2, q2d - w1 * s3 - w2 *
                        c3, q3d - (w1 * c3 - w2 * s3) * c2 / s2 - w3]
        if rot_type.lower() == 'space':
            if rot_order == '123':
                return [q1d - w1 - (w2 * s1 + w3 * c1) * s2 / c2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / c2]
            if rot_order == '231':
                return [q1d - (w1 * c1 + w3 * s1) * s2 / c2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / c2]
            if rot_order == '312':
                return [q1d - (w1 * s1 + w2 * c1) * s2 / c2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / c2]
            if rot_order == '132':
                return [q1d - w1 - (-w2 * c1 + w3 * s1) * s2 / c2, q2d - w2 *
                        s1 - w3 * c1, q3d - (w2 * c1 - w3 * s1) / c2]
            if rot_order == '213':
                return [q1d - (w1 * s1 - w3 * c1) * s2 / c2 - w2, q2d - w1 *
                        c1 - w3 * s1, q3d - (-w1 * s1 + w3 * c1) / c2]
            if rot_order == '321':
                return [q1d - (-w1 * c1 + w2 * s1) * s2 / c2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (w1 * c1 - w2 * s1) / c2]
            if rot_order == '121':
                return [q1d - w1 + (w2 * s1 + w3 * c1) * c2 / s2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / s2]
            if rot_order == '131':
                return [q1d - w1 - (w2 * c1 - w3 * s1) * c2 / s2, q2d - w2 *
                        s1 - w3 * c1, q3d - (-w2 * c1 + w3 * s1) / s2]
            if rot_order == '212':
                return [q1d - (-w1 * s1 + w3 * c1) * c2 / s2 - w2, q2d - w1 *
                        c1 - w3 * s1, q3d - (w1 * s1 - w3 * c1) / s2]
            if rot_order == '232':
                return [q1d + (w1 * c1 + w3 * s1) * c2 / s2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / s2]
            if rot_order == '313':
                return [q1d + (w1 * s1 + w2 * c1) * c2 / s2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / s2]
            if rot_order == '323':
                return [q1d - (w1 * c1 - w2 * s1) * c2 / s2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (-w1 * c1 + w2 * s1) / s2]
    elif rot_type.lower() == 'quaternion':
        if rot_order != '':
            raise ValueError('Cannot have rotation order for quaternion')
        if len(coords) != 4:
            raise ValueError('Need 4 coordinates for quaternion')
        # Actual hard-coded kinematic differential equations
        e0, e1, e2, e3 = coords
        w = Matrix(speeds + [0])
        E = Matrix([[e0, -e3, e2, e1], [e3, e0, -e1, e2], [-e2, e1, e0, e3],
            [-e1, -e2, -e3, e0]])
        edots = Matrix([diff(i, dynamicsymbols._t) for i in [e1, e2, e3, e0]])
        return list(edots.T - 0.5 * w.T * E.T)
    else:
        raise ValueError('Not an approved rotation type for this function')


def get_motion_params(frame, **kwargs):
    """
    Returns the three motion parameters - (acceleration, velocity, and
    position) as vectorial functions of time in the given frame.

    If a higher order differential function is provided, the lower order
    functions are used as boundary conditions. For example, given the
    acceleration, the velocity and position parameters are taken as
    boundary conditions.

    The values of time at which the boundary conditions are specified
    are taken from timevalue1(for position boundary condition) and
    timevalue2(for velocity boundary condition).

    If any of the boundary conditions are not provided, they are taken
    to be zero by default (zero vectors, in case of vectorial inputs). If
    the boundary conditions are also functions of time, they are converted
    to constants by substituting the time values in the dynamicsymbols._t
    time Symbol.

    This function can also be used for calculating rotational motion
    parameters. Have a look at the Parameters and Examples for more clarity.

    Parameters
    ==========

    frame : ReferenceFrame
        The frame to express the motion parameters in

    acceleration : Vector
        Acceleration of the object/frame as a function of time

    velocity : Vector
        Velocity as function of time or as boundary condition
        of velocity at time = timevalue1

    position : Vector
        Velocity as function of time or as boundary condition
        of velocity at time = timevalue1

    timevalue1 : sympyfiable
        Value of time for position boundary condition

    timevalue2 : sympyfiable
        Value of time for velocity boundary condition

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, get_motion_params, dynamicsymbols
    >>> from sympy import symbols
    >>> R = ReferenceFrame('R')
    >>> v1, v2, v3 = dynamicsymbols('v1 v2 v3')
    >>> v = v1*R.x + v2*R.y + v3*R.z
    >>> get_motion_params(R, position = v)
    (v1''*R.x + v2''*R.y + v3''*R.z, v1'*R.x + v2'*R.y + v3'*R.z, v1*R.x + v2*R.y + v3*R.z)
    >>> a, b, c = symbols('a b c')
    >>> v = a*R.x + b*R.y + c*R.z
    >>> get_motion_params(R, velocity = v)
    (0, a*R.x + b*R.y + c*R.z, a*t*R.x + b*t*R.y + c*t*R.z)
    >>> parameters = get_motion_params(R, acceleration = v)
    >>> parameters[1]
    a*t*R.x + b*t*R.y + c*t*R.z
    >>> parameters[2]
    a*t**2/2*R.x + b*t**2/2*R.y + c*t**2/2*R.z

    """

    ##Helper functions

    def _process_vector_differential(vectdiff, condition, \
                                     variable, ordinate, frame):
        """
        Helper function for get_motion methods. Finds derivative of vectdiff wrt
        variable, and its integral using the specified boundary condition at
        value of variable = ordinate.
        Returns a tuple of - (derivative, function and integral) wrt vectdiff

        """

        #Make sure boundary condition is independent of 'variable'
        if condition != 0:
            condition = express(condition, frame, variables=True)
        #Special case of vectdiff == 0
        if vectdiff == Vector(0):
            return (0, 0, condition)
        #Express vectdiff completely in condition's frame to give vectdiff1
        vectdiff1 = express(vectdiff, frame)
        #Find derivative of vectdiff
        vectdiff2 = time_derivative(vectdiff, frame)
        #Integrate and use boundary condition
        vectdiff0 = Vector(0)
        lims = (variable, ordinate, variable)
        for dim in frame:
            function1 = vectdiff1.dot(dim)
            abscissa = dim.dot(condition).subs({variable : ordinate})
            # Indefinite integral of 'function1' wrt 'variable', using
            # the given initial condition (ordinate, abscissa).
            vectdiff0 += (integrate(function1, lims) + abscissa) * dim
        #Return tuple
        return (vectdiff2, vectdiff, vectdiff0)

    ##Function body

    _check_frame(frame)
    #Decide mode of operation based on user's input
    if 'acceleration' in kwargs:
        mode = 2
    elif 'velocity' in kwargs:
        mode = 1
    else:
        mode = 0
    #All the possible parameters in kwargs
    #Not all are required for every case
    #If not specified, set to default values(may or may not be used in
    #calculations)
    conditions = ['acceleration', 'velocity', 'position',
                  'timevalue', 'timevalue1', 'timevalue2']
    for i, x in enumerate(conditions):
        if x not in kwargs:
            if i < 3:
                kwargs[x] = Vector(0)
            else:
                kwargs[x] = S(0)
        elif i < 3:
            _check_vector(kwargs[x])
        else:
            kwargs[x] = sympify(kwargs[x])
    if mode == 2:
        vel = _process_vector_differential(kwargs['acceleration'],
                                           kwargs['velocity'],
                                           dynamicsymbols._t,
                                           kwargs['timevalue2'], frame)[2]
        pos = _process_vector_differential(vel, kwargs['position'],
                                           dynamicsymbols._t,
                                           kwargs['timevalue1'], frame)[2]
        return (kwargs['acceleration'], vel, pos)
    elif mode == 1:
        return _process_vector_differential(kwargs['velocity'],
                                            kwargs['position'],
                                            dynamicsymbols._t,
                                            kwargs['timevalue1'], frame)
    else:
        vel = time_derivative(kwargs['position'], frame)
        acc = time_derivative(vel, frame)
        return (acc, vel, kwargs['position'])


def partial_velocity(vel_list, u_list, frame):
    """Returns a list of partial velocities.

    For a list of velocity or angular velocity vectors the partial derivatives
    with respect to the supplied generalized speeds are computed, in the
    specified ReferenceFrame.

    The output is a list of lists. The outer list has a number of elements
    equal to the number of supplied velocity vectors. The inner lists are, for
    each velocity vector, the partial derivatives of that velocity vector with
    respect to the generalized speeds supplied.

    Parameters
    ==========

    vel_list : list
        List of velocities of Point's and angular velocities of ReferenceFrame's
    u_list : list
        List of independent generalized speeds.
    frame : ReferenceFrame
        The ReferenceFrame the partial derivatives are going to be taken in.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, ReferenceFrame
    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> from sympy.physics.mechanics import partial_velocity
    >>> u = dynamicsymbols('u')
    >>> N = ReferenceFrame('N')
    >>> P = Point('P')
    >>> P.set_vel(N, u * N.x)
    >>> vel_list = [P.vel(N)]
    >>> u_list = [u]
    >>> partial_velocity(vel_list, u_list, N)
    [[N.x]]

    """
    if not hasattr(vel_list, '__iter__'):
        raise TypeError('Provide velocities in an iterable')
    if not hasattr(u_list, '__iter__'):
        raise TypeError('Provide speeds in an iterable')
    list_of_pvlists = []
    for i in vel_list:
        pvlist = []
        for j in u_list:
            vel = i.diff(j, frame)
            pvlist += [vel]
        list_of_pvlists += [pvlist]
    return list_of_pvlists


def linear_momentum(frame, *body):
    """Linear momentum of the system.

    This function returns the linear momentum of a system of Particle's and/or
    RigidBody's. The linear momentum of a system is equal to the vector sum of
    the linear momentum of its constituents. Consider a system, S, comprised of
    a rigid body, A, and a particle, P. The linear momentum of the system, L,
    is equal to the vector sum of the linear momentum of the particle, L1, and
    the linear momentum of the rigid body, L2, i.e-

    L = L1 + L2

    Parameters
    ==========

    frame : ReferenceFrame
        The frame in which linear momentum is desired.
    body1, body2, body3... : Particle and/or RigidBody
        The body (or bodies) whose kinetic energy is required.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, Particle, ReferenceFrame
    >>> from sympy.physics.mechanics import RigidBody, outer, linear_momentum
    >>> N = ReferenceFrame('N')
    >>> P = Point('P')
    >>> P.set_vel(N, 10 * N.x)
    >>> Pa = Particle('Pa', P, 1)
    >>> Ac = Point('Ac')
    >>> Ac.set_vel(N, 25 * N.y)
    >>> I = outer(N.x, N.x)
    >>> A = RigidBody('A', Ac, N, 20, (I, Ac))
    >>> linear_momentum(N, A, Pa)
    10*N.x + 500*N.y

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Please specify a valid ReferenceFrame')
    else:
        linear_momentum_sys = Vector(0)
        for e in body:
            if isinstance(e, (RigidBody, Particle)):
                linear_momentum_sys += e.linear_momentum(frame)
            else:
                raise TypeError('*body must have only Particle or RigidBody')
    return linear_momentum_sys


def angular_momentum(point, frame, *body):
    """Angular momentum of a system

    This function returns the angular momentum of a system of Particle's and/or
    RigidBody's. The angular momentum of such a system is equal to the vector
    sum of the angular momentum of its constituents. Consider a system, S,
    comprised of a rigid body, A, and a particle, P. The angular momentum of
    the system, H, is equal to the vector sum of the linear momentum of the
    particle, H1, and the linear momentum of the rigid body, H2, i.e-

    H = H1 + H2

    Parameters
    ==========

    point : Point
        The point about which angular momentum of the system is desired.
    frame : ReferenceFrame
        The frame in which angular momentum is desired.
    body1, body2, body3... : Particle and/or RigidBody
        The body (or bodies) whose kinetic energy is required.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, Particle, ReferenceFrame
    >>> from sympy.physics.mechanics import RigidBody, outer, angular_momentum
    >>> N = ReferenceFrame('N')
    >>> O = Point('O')
    >>> O.set_vel(N, 0 * N.x)
    >>> P = O.locatenew('P', 1 * N.x)
    >>> P.set_vel(N, 10 * N.x)
    >>> Pa = Particle('Pa', P, 1)
    >>> Ac = O.locatenew('Ac', 2 * N.y)
    >>> Ac.set_vel(N, 5 * N.y)
    >>> a = ReferenceFrame('a')
    >>> a.set_ang_vel(N, 10 * N.z)
    >>> I = outer(N.z, N.z)
    >>> A = RigidBody('A', Ac, a, 20, (I, Ac))
    >>> angular_momentum(O, N, Pa, A)
    10*N.z

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Please enter a valid ReferenceFrame')
    if not isinstance(point, Point):
        raise TypeError('Please specify a valid Point')
    else:
        angular_momentum_sys = Vector(0)
        for e in body:
            if isinstance(e, (RigidBody, Particle)):
                angular_momentum_sys += e.angular_momentum(point, frame)
            else:
                raise TypeError('*body must have only Particle or RigidBody')
    return angular_momentum_sys


def kinetic_energy(frame, *body):
    """Kinetic energy of a multibody system.

    This function returns the kinetic energy of a system of Particle's and/or
    RigidBody's. The kinetic energy of such a system is equal to the sum of
    the kinetic energies of its constituents. Consider a system, S, comprising
    a rigid body, A, and a particle, P. The kinetic energy of the system, T,
    is equal to the vector sum of the kinetic energy of the particle, T1, and
    the kinetic energy of the rigid body, T2, i.e.

    T = T1 + T2

    Kinetic energy is a scalar.

    Parameters
    ==========

    frame : ReferenceFrame
        The frame in which the velocity or angular velocity of the body is
        defined.
    body1, body2, body3... : Particle and/or RigidBody
        The body (or bodies) whose kinetic energy is required.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, Particle, ReferenceFrame
    >>> from sympy.physics.mechanics import RigidBody, outer, kinetic_energy
    >>> N = ReferenceFrame('N')
    >>> O = Point('O')
    >>> O.set_vel(N, 0 * N.x)
    >>> P = O.locatenew('P', 1 * N.x)
    >>> P.set_vel(N, 10 * N.x)
    >>> Pa = Particle('Pa', P, 1)
    >>> Ac = O.locatenew('Ac', 2 * N.y)
    >>> Ac.set_vel(N, 5 * N.y)
    >>> a = ReferenceFrame('a')
    >>> a.set_ang_vel(N, 10 * N.z)
    >>> I = outer(N.z, N.z)
    >>> A = RigidBody('A', Ac, a, 20, (I, Ac))
    >>> kinetic_energy(N, Pa, A)
    350

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Please enter a valid ReferenceFrame')
    ke_sys = S(0)
    for e in body:
        if isinstance(e, (RigidBody, Particle)):
            ke_sys += e.kinetic_energy(frame)
        else:
            raise TypeError('*body must have only Particle or RigidBody')
    return ke_sys


def potential_energy(*body):
    """Potential energy of a multibody system.

    This function returns the potential energy of a system of Particle's and/or
    RigidBody's. The potential energy of such a system is equal to the sum of
    the potential energy of its constituents. Consider a system, S, comprising
    a rigid body, A, and a particle, P. The potential energy of the system, V,
    is equal to the vector sum of the potential energy of the particle, V1, and
    the potential energy of the rigid body, V2, i.e.

    V = V1 + V2

    Potential energy is a scalar.

    Parameters
    ==========

    body1, body2, body3... : Particle and/or RigidBody
        The body (or bodies) whose potential energy is required.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, Particle, ReferenceFrame
    >>> from sympy.physics.mechanics import RigidBody, outer, potential_energy
    >>> from sympy import symbols
    >>> M, m, g, h = symbols('M m g h')
    >>> N = ReferenceFrame('N')
    >>> O = Point('O')
    >>> O.set_vel(N, 0 * N.x)
    >>> P = O.locatenew('P', 1 * N.x)
    >>> Pa = Particle('Pa', P, m)
    >>> Ac = O.locatenew('Ac', 2 * N.y)
    >>> a = ReferenceFrame('a')
    >>> I = outer(N.z, N.z)
    >>> A = RigidBody('A', Ac, a, M, (I, Ac))
    >>> Pa.set_potential_energy(m * g * h)
    >>> A.set_potential_energy(M * g * h)
    >>> potential_energy(Pa, A)
    M*g*h + g*h*m

    """

    pe_sys = S(0)
    for e in body:
        if isinstance(e, (RigidBody, Particle)):
            pe_sys += e.potential_energy
        else:
            raise TypeError('*body must have only Particle or RigidBody')
    return pe_sys


def Lagrangian(frame, *body):
    """Lagrangian of a multibody system.

    This function returns the Lagrangian of a system of Particle's and/or
    RigidBody's. The Lagrangian of such a system is equal to the difference
    between the kinetic energies and potential energies of its constituents. If
    T and V are the kinetic and potential energies of a system then it's
    Lagrangian, L, is defined as

    L = T - V

    The Lagrangian is a scalar.

    Parameters
    ==========

    frame : ReferenceFrame
        The frame in which the velocity or angular velocity of the body is
        defined to determine the kinetic energy.

    body1, body2, body3... : Particle and/or RigidBody
        The body (or bodies) whose kinetic energy is required.

    Examples
    ========

    >>> from sympy.physics.mechanics import Point, Particle, ReferenceFrame
    >>> from sympy.physics.mechanics import RigidBody, outer, Lagrangian
    >>> from sympy import symbols
    >>> M, m, g, h = symbols('M m g h')
    >>> N = ReferenceFrame('N')
    >>> O = Point('O')
    >>> O.set_vel(N, 0 * N.x)
    >>> P = O.locatenew('P', 1 * N.x)
    >>> P.set_vel(N, 10 * N.x)
    >>> Pa = Particle('Pa', P, 1)
    >>> Ac = O.locatenew('Ac', 2 * N.y)
    >>> Ac.set_vel(N, 5 * N.y)
    >>> a = ReferenceFrame('a')
    >>> a.set_ang_vel(N, 10 * N.z)
    >>> I = outer(N.z, N.z)
    >>> A = RigidBody('A', Ac, a, 20, (I, Ac))
    >>> Pa.set_potential_energy(m * g * h)
    >>> A.set_potential_energy(M * g * h)
    >>> Lagrangian(N, Pa, A)
    -M*g*h - g*h*m + 350

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Please supply a valid ReferenceFrame')
    for e in body:
        if not isinstance(e, (RigidBody, Particle)):
            raise TypeError('*body must have only Particle or RigidBody')
    return kinetic_energy(frame, *body) - potential_energy(*body)
