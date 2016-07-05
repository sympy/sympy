# This file will contain the joint library and jcalc functions required for Roy
# Featherstone's methods

# TODO figure out if I'll completely leave out generalized speeds from
# Featherstone implementation


def jcalc():
    """This function provides all the joint information needed to run
    Featherstone's different methods. Additional information on this function
    can be found in Chapter 4.4 of his book.

    In the book, Featherstone presents the following methods of using jcalc on
    page 80:

        [1] [XJ, S, vJ, cJ] = jcalc(jtype, q, alpha)
        [2] [delta, T] = jcalc(jtype, sXp)
        [3] [qdot] = jcalc(jtype, q, alpha)

    For now this version of jcalc will only perform operation [1].

    Parameters
    ==========

    jtype: string
        Indicator of the joint type. The current accepted parameters are
            'R' -> Revolute
            'P' -> Prismatic
            'H' -> Helical
            'USER' -> User defined jcalc
    q: ordered interable of functions of time
        The set of generalized coordinates for the joint
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
        >>> [XJ, S, vJ, cJ] = jcalc('H', q, pitch=3)

    Notes
    =====

    nf: Degrees of relative motion freedom allowed by the joint

    References
    ==========

    .. [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
       Springer
    """

    # Reminder to use the correct coordinate for rotation and translation from
    # the 6x1 vector of generalized coordinates.

################################################################################
#
# Joint Library
#
################################################################################

# In this portion of the file the parameters for different types of joints will
# be defined. These parameters are given in Table 4.1 on page 79 in the book.


def revolute():
    """This fucntion will return revolute joint information. This assumes
    standard link coordinates.

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

        >>> from sympy.physics.mechanics.joints import revolute
        >>> [Xj, S, T] = revolute()
    """


def prismatic():
    """This fucntion will return prismatic joint information. This assumes
    standard link coordinates.


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

        >>> from sympy.physics.mechanics.joints import prismatic
        >>> [Xj, S, T] = prismatic()
    """


def helical(pitch=None):
    """This fucntion will return helical joint information. This assumes
    standard link coordinates. If a pitch is given it will be substituted into
    the returns.

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
        >>> [XJ, S, T] = helical()
        >>> [XJ, S, T] = helical(pitch=symbols('p'))
        >>> [XJ, S, T] = helical(pitch=5)
    """
