# This file will contain the joint library and jcalc functions required for Roy
# Featherstone's methods

# TODO figure out if I'll completely leave out generalized speeds from
# Featherstone implementation

from sympy import cos, Matrix, sin, symbols


def jcalc(jtype, q, alpha=None, jparam={}, user_func=None):
    """This function provides all the joint information needed to run
    Featherstone's different methods. Additional information on this function
    can be found in Chapter 4.4 of his book.

    In the book, Featherstone presents the following methods of using jcalc on
    page 80:

        [1] [XJ, S, vJ, cJ] = jcalc(jtype, q, alpha)
        [2] [delta, T] = jcalc(jtype, sXp)
        [3] [qdot] = jcalc(jtype, q, alpha)

    For now this version of jcalc will only perform operation [1].  Also, this
    function assumes there is no explicit time dependency of the joint. (TODO:
        Should I completely get rid of return cJ which is zero when there is no
        explicit time dependency? Could be left so that user functions can add
        explicit time dependencies and the code will know how to handle them)

    Parameters
    ==========

    jtype: string
        Indicator of the joint type. The current accepted parameters are
            'R' -> Revolute
            'P' -> Prismatic
            'H' -> Helical
            'USER' -> User defined jcalc
    q: ordered interable of functions of time
        The set of generalized coordinates for the joint. It is expected that
        the coordinates will be given in the order [w_x, w_y, w_z, x, y, z]
    alpha: ordered iterable of functions of time, optional
        The set of generarlized speeds. If not specified it will be assumed
        qdot = alpha
    jparam: dictionary, optional
        This is all of the joint parameters needed for the specific joint type
        (ex. helical joints need a pitch defined)
    user_func: function, optional
        This allows the user to define their own jcalc function and will return
        the values produced by the input function. When used all of the other
        parameters will be passed to this function.

    Returns
    =======

    TODO: XJ a 6x6 matrix? and if cJ is a motion spatial
    vector and if I'm implementing spatial vectors as their own objects. Also do
    the functions for vJ and cJ use qdot or do they really use alpha?

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, nf)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    vJ: Motion Spatial Vector
        The joint velocity, (successor body velocity minus predecessor body
        velocity, vs - vp). Defined on page 53
    cJ: Motion Spatial Vector
        The velocity product term defined on page 55

    Examples
    ========

    The first step for the examples is to form the vector of generalized
    coordinates. ::

        >>> from sympy import symbols
        >>> the_x, the_y, the_z, x, y, z = symbols('the_x the_y the_z x y z')
        >>> q = [the_x, the_y, the_z, x, y, z]

    Now use the jcalc function to recieve the necessary joint parameters. ::

        >>> from sympy.physics.mechanics.joints import jcalc
        >>> [XJ, S, vJ, cJ] = jcalc('R', q)
        >>> [XJ, S, vJ, cJ] = jcalc('P', q)
        >>> [XJ, S, vJ, cJ] = jcalc('H', q, jparam={"pitch":3})

    Notes
    =====

    nf: Degrees of relative motion freedom allowed by the joint

    References
    ==========

    .. [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
       Springer
    """

    # Pull data from the joint library based on the input joint type
    jtype.upper()
    if jtype == 'R':
        [XJ, S, T] = revolute(q[2])
    elif jtype == 'P':
        [XJ, S, T] = prismatic(q[5])
    elif jtype == 'H':
        pitch = jparam.pop('pitch', None)
        [XJ, S, T] = helical(q[2], pitch)
    elif jtype == 'USER':
        return user_func(q, alpha, jparam)
    else:
        raise ValueError("Argument %s for jtype is not a valid option" % jtype)

    # Calculate the joint velocity
    qdot = Matrix([i.diff() for i in q])
    vJ = S.transpose() * qdot

    return XJ, S, vJ, 0



################################################################################
#
# Joint Library
#
################################################################################

# In this portion of the file the parameters for different types of joints will
# be defined. These parameters are given in Table 4.1 on page 79 in the book.


def revolute(q):
    """This fucntion will return revolute joint information. This assumes
    standard link coordinates.

    Parameters
    ==========

    q: SymPy symbol
        This is the generalized angle that defines the position of the successor
        link to the predecessor link (joint angle).

    Returns
    =======

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, 1)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    T: Matrix, size(6, 5)
        The constraint force subspace matrix for the joint

    Examples
    ========

    This function is pretty straight forward as it does not need to calculate
    anything based on changing input conditions. It simply returns the matrices
    that define a revolute joint. ::

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics.joints import revolute
        >>> theta = symbols('theta')
        >>> [XJ, S, T] = revolute(theta)
    """

    XJ = Matrix([[cos(q),  sin(q), 0, 0,       0,      0],
                 [-sin(q), cos(q), 0, 0,       0,      0],
                 [0,           0,  1, 0,       0,      0],
                 [0,           0,  0, cos(q),  sin(q), 0],
                 [0,           0,  0, -sin(q), cos(q), 0],
                 [0,           0,  0, 0,       0,      1]])

    S = Matrix([0, 0, 1, 0, 0, 0])

    T = Matrix([[1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1]])

    return XJ, S, T


def prismatic(q):
    """This fucntion will return prismatic joint information. This assumes
    standard link coordinates.

    Parameters
    ==========

    q: SymPy symbol
        This is the generalized length that defines the position of the
        successor link to the predecessor link (joint offset distance).

    Returns
    =======

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, 1)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    T: Matrix, size(6, 5)
        The constraint force subspace matrix for the joint

    Examples
    ========

    This function is pretty straight forward as it does not need to calculate
    anything based on changing input conditions. It simply returns the matrices
    that define a prismatic joint. ::

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics.joints import prismatic
        >>> L = symbols('L')
        >>> [XJ, S, T] = prismatic(L)
    """

    XJ = Matrix([[1,  0, 0, 0, 0, 0],
                 [0,  1, 0, 0, 0, 0],
                 [0,  0, 1, 0, 0, 0],
                 [0,  q, 0, 1, 0, 0],
                 [-q, 0, 0, 0, 1, 0],
                 [0,  0, 0, 0, 0, 1]])

    S = Matrix([0, 0, 0, 0, 0, 1])

    T = Matrix([[1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0]])

    return XJ, S, T


def helical(q, pitch=None):
    """This fucntion will return helical joint information. This assumes
    standard link coordinates. If a pitch is given it will be substituted into
    the returns.

    Parameters
    ==========

    q: SymPy symbol
        This is the generalized angle that defines the position of the successor
        link to the predecessor link (joint angle).
    pitch: SymPy symbol or numeric, optional
        This represents the pitch of the helical joint. If not specified a sympy
        symbol of 'h' will be used.

    Returns
    =======

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, 1)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    T: Matrix, size(6, 5)
        The constraint force subspace matrix for the joint

    Examples
    ========

    This function will return the joint data for a helical joint and a specific
    number or symbol can be used for the pitch. ::

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics.joints import helical
        >>> theta, p = symbols('theta p')
        >>> [XJ, S, T] = helical(theta)
        >>> [XJ, S, T] = helical(theta, pitch=p)
        >>> [XJ, S, T] = helical(theta, pitch=5)
    """

    # Let the user have the ability to set the pitch to a specific
    # variable/value and default to the symbol 'h' otherwise
    if pitch is None:
        h = symbols('h')
    else:
        h = pitch

    rotate = Matrix([[cos(q),  sin(q), 0, 0,       0,      0],
                     [-sin(q), cos(q), 0, 0,       0,      0],
                     [0,           0,  1, 0,       0,      0],
                     [0,           0,  0, cos(q),  sin(q), 0],
                     [0,           0,  0, -sin(q), cos(q), 0],
                     [0,           0,  0, 0,       0,      1]])

    translate = Matrix([[1,    0,   0, 0, 0, 0],
                        [0,    1,   0, 0, 0, 0],
                        [0,    0,   1, 0, 0, 0],
                        [0,    h*q, 0, 1, 0, 0],
                        [-h*q, 0,   0, 0, 1, 0],
                        [0,    0,   0, 0, 0, 1]])

    XJ = rotate*translate

    S = Matrix([0, 0, 1, 0, 0, h])

    T = Matrix([[1, 0, 0, 0,  0],
                [0, 1, 0, 0,  0],
                [0, 0, 0, 0, -h],
                [0, 0, 1, 0,  0],
                [0, 0, 0, 1,  0],
                [0, 0, 0, 0,  1]])

    return XJ, S, T
