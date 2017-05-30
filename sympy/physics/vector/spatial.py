# This code will assume the spatial vectors are given such that the first three
# elements correspond to rotations and the last three elements correspond to
# translations. This is not a requirement of spatial vectors in general but
# rather a selected convention.

# References:
#     [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
#     Springer

from sympy import cos, Matrix, sin, zeros


def rx(theta):
    """This function returns the 3x3 rotation matrix for a rotation of theta
    degrees about the x axis.

    Parameters
    ==========

    theta: numeric (int, float, etc)
        This is the angle of the rotation around the x-axis. The angle is given
        in degrees

    Returns
    =======

    E: Matrix, size(3, 3)
        This is the calculated transformation matrix

    Examples
    ========

    Simply provide an angle and the function will return the rotation matrix. ::

        >>> from sympy.physics.vector.spatial import rx
        >>> theta = 35
        >>> E = rx(theta)
    """

    E = Matrix([[1,           0,          0],
                [0,  cos(theta), sin(theta)],
                [0, -sin(theta), cos(theta)]])

    return E


def ry(theta):
    """This function returns the 3x3 rotation matrix for a rotation of theta
    degrees about the y axis.

    Parameters
    ==========

    theta: numeric (int, float, etc)
        This is the angle of the rotation around the y-axis. The angle is given
        in degrees

    Returns
    =======

    E: Matrix, size(3, 3)
        This is the calculated transformation matrix

    Examples
    ========

    Simply provide an angle and the function will return the rotation matrix. ::

        >>> from sympy.physics.vector.spatial import ry
        >>> theta = 35
        >>> E = ry(theta)
    """

    E = Matrix([[cos(theta), 0, -sin(theta)],
                [0,          1,           0],
                [sin(theta), 0,  cos(theta)]])

    return E


def rz(theta):
    """This function returns the 3x3 rotation matrix for a rotation of theta
    degrees about the z axis.

    Parameters
    ==========

    theta: numeric (int, float, etc)
        This is the angle of the rotation around the z-axis. The angle is given
        in degrees

    Returns
    =======

    E: Matrix, size(3, 3)
        This is the calculated transformation matrix

    Examples
    ========

    Simply provide an angle and the function will return the rotation matrix. ::

        >>> from sympy.physics.vector.spatial import rz
        >>> theta = 35
        >>> E = rz(theta)
    """

    E = Matrix([[cos(theta),  sin(theta), 0],
                [-sin(theta), cos(theta), 0],
                [0,           0,          1]])

    return E


def rot(E):
    """This function returns the 6x6 rotation matrix for spatial vectors as
    defined in table 2.2 on page 23 of Featherstone's Rigid Body Dynamics
    Algorithms book.

    Parameters
    ==========

    E: Matrix, size(3, 3)
        This is the 3x3 rotation matrix representing the rotation in 3D
        coordinates

    Returns
    =======

    X: Matrix, size(6, 6)
        This is the 6x6 rotation matrix for spatial vectors

    Examples
    ========

    A simple example would be a single rotation of theta degrees about the x
    axis. ::

        >>> from sympy.physics.vector.spatial import rx, rot
        >>> theta = 35
        >>> E = rx(theta)
        >>> X = rot(E)

    This takes the 3x3 rotation matrix and puts it in a form to be used with
    spatial vectors.
    """

    X = E.row_join(zeros(3)).col_join(zeros(3).row_join(E))

    return X


def xlt(r):
    """This function returns the 6x6 translation matrix for spatial vectors as
    defined in table 2.2 on page 23 of Featherstones' Rigid Body Dynamics
    Algorithms book.

    Parameters
    ==========

    r: Matrix, size(3, 1)
        This is the 3D vector representing the translation from the original
        origin point to the new origin point

    Returns
    =======

    X: Matrix, size(6,6)
        This is the 6x6 translation matrix for spatial vectors

    Examples
    ========

    For a unit translation in the x direction the transformation matrix for
    spatial vectors can be formed as follows::

        >>> from sympy import Matrix
        >>> from sympy.physics.vector.spatial import xlt
        >>> r = Matrix([1, 0, 0])
        >>> X = xlt(r)
    """

    X = Matrix([[1,        0,     0, 0, 0, 0],
                [0,        1,     0, 0, 0, 0],
                [0,        0,     1, 0, 0, 0],
                [0,     r[2], -r[1], 1, 0, 0],
                [-r[2],    0,  r[0], 0, 1, 0],
                [r[1], -r[0],     0, 0, 0, 1]])

    return X


def cross_f(v):
    """This function returns the force space cross product (vx*). This operation
    is defined in Featherstone's Rigid Body Dynamics Algorithms book.

    Parameters
    ==========

    v: Matrix, size(6, 1)
        This is the input vector for the cross product vx*.

    Returns
    =======

    crossed: Matrix, size(6, 6)
        This is the result of the force cross product.

    Examples
    ========

    A simple example of the cross product of a vector v. ::

        >>> from sympy import Matrix, symbols
        >>> from sympy.physics.vector.spatial import cross_f
        >>> v1, v2, v3, v4, v5, v6 = symbols('v1 v2 v3 v4 v5 v6')
        >>> v = Matrix([v1, v2, v3, v4, v5, v6])
        >>> crossed_f = cross_f(v)

    Notes
    =====

    It is worth noting that the force cross product is the negative transpose of
    the motion cross product (vx* = -(vx)^T)
    """

    crossed = -cross_m(v).transpose()

    return crossed


def cross_m(v):
    """This function returns the motion cross product (vx). The operation will
    depend on the size of the matrix entered. If a 3x1 vector is given the 3x3
    cross product matrix will be returned. Otherwise the code will assume that a
    6x1 matrix was given and will return the 6x6 cross product matrix.

    Parameters
    ==========

    v: Matrix, size(6, 1) or size(3, 1)
        This is the input vector for the cross product vx.

    Returns
    =======

    crossed: Matrix, size(6, 6) or size(3, 3)
        This is the result of the motion cross product.

    Examples
    ========

    A simple example of the cross product of a vector v. ::

        >>> from sympy import Matrix, symbols
        >>> from sympy.physics.vector.spatial import cross_m
        >>> v1, v2, v3, v4, v5, v6 = symbols('v1 v2 v3 v4 v5 v6')
        >>> v = Matrix([v1, v2, v3, v4, v5, v6])
        >>> crossed_m = cross_m(v)
    """

    if v.shape == (3, 1):
        v1, v2, v3 = v
        crossed = Matrix([[0, -v3, v2],
                          [v3, 0, -v1],
                          [-v2, v1, 0]])
    else:
        v1, v2, v3, v4, v5, v6 = v

        crossed = Matrix([[0,  -v3,  v2,  0,   0,   0],
                          [v3,   0, -v1,  0,   0,   0],
                          [-v2, v1,   0,  0,   0,   0],
                          [0,  -v6,  v5,  0, -v3,  v2],
                          [v6,   0, -v4, v3,   0, -v1],
                          [-v5, v4,   0, -v2, v1,   0]])

    return crossed
