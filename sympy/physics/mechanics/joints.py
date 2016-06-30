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
        The set of generalized coordinates
    alpha: ordered iterable of functions of time, optional
        The set of generarlized speeds. If not specified it will be assumed
        qdot = alpha
    jparam: dictionary, optional
        This is all of the joint parameters needed for the specific joint type
        (ex. helical joints need a pitch defined)
    user_func: function, optional
        This allows the user to define their own jcalc function and will return
        the values produced by the input function

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

    TODO: Provide an example when I figure out how I'm going to implement the
    spatial vectors

    Notes
    =====

    nf: Degrees of relative motion freedom allowed by the joint

    References
    ==========

    .. [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
       Springer
    """

################################################################################
#
# Joint Library
#
################################################################################

# In this portion of the file the parameters for different types of joints will
# be defined. These parameters are given in Table 4.1 on page 79 in the book.


def revolute():
    """This fucntion will return revolute joint information. This assumes
    standard link coordinates

    Returns
    =======

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, 1)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    T: Matrix, size(6, 5)
        The constraint force subspace matrix for the joint
    """


def prismatic():
    """This fucntion will return prismatic joint information. This assumes
    standard link coordinates


    Returns
    =======

    XJ: Matrix, size(6, 6)
        The joint transform matrix
    S: Matrix, size(6, 1)
        The motion subspace matrix describing the motion allowed by the joint
        constraint
    T: Matrix, size(6, 5)
        The constraint force subspace matrix for the joint
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
    """
