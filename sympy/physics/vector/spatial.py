# This file will contain spatial vector functions that will be used in the
# implementation of Featherstone's methods. The functions will be follow the
# nomenclature used in Featherstone's "Rigid Body Dynamics Algorithms" book.

# This code will assume the spatial vectors are given such that the first three
# elements correspond to rotations and the last three elements correspond to
# translations. This is not a requirement of spatial vectors in general but
# rather a selected convention.

# References:
#     [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
#     Springer

################################################################################
#
# Rotations
#
################################################################################

# Roations and nomenclature taken from table 2.2 on page 23


def rx():
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
    """


def ry():
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
    """


def rz():
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
    """


def rot():
    """This function returns the 6x6 rotation matrix for spatial vectors.

    Parameters
    ==========

    E: Matrix, size(3, 3)
        This is the 3x3 rotation matrix representing the rotation in 3D
        coordinates

    Returns
    =======

    X: Matrix, size(6, 6)
        This is the 6x6 rotation matrix for spatial vectors
    """


################################################################################
#
# Translations
#
################################################################################

# Formulas and nomeclature taken from table 2.2 on page 23

def xlt():
    """This function returns the 6x6 translation matrix for spatial vectors.

    Parameters
    ==========

    r: Matrix, size(3, 1)
        This is the 3D vector representing the translation from the original
        origin point to the new origin point

    Returns
    =======

    X: Matrix, size(6,6)
        This is the 6x6 translation matrix for spatial vectors
    """


################################################################################
#
# Spatial cross products
#
################################################################################

# The spatial cross product tables can be found in figure 2.5 on page 24

def cross_f():
    """This function returns the force cross product (vx*). If input v is a 3x1
    vector a 3x3 matrix will be returned, otherwise a 6x6 matrix will be
    returned.

    Parameters
    ==========

    v: Matrix, size(6, 1) or size(3, 1)
        This is the input vector for the cross product vx*. If v is not 6x1 it
        will be assumed that it is 3x1.

    Returns
    =======

    crossed: Matrix, size(6, 6) or size(3, 3)
        This is the result of the force cross product. It v was 6x1 the return
        is a 6x6 matrix otherwise it is a 3x3 matrix

    Notes
    =====

    It is worth noting that the force cross product is the negative transpose of
    the motion cross product (vx* = -(vx)^T)
    """


def cross_m():
    """This function returns the motion cross product (vx). If input v is a 3x1
    vector a 3x3 matrix will be returned, otherwise a 6x6 matrix will be
    returned.

    Parameters
    ==========

    v: Matrix, size(6, 1) or size(3, 1)
        This is the input vector for the cross product vx. If v is not 6x1 it
        will be assumed that it is 3x1.

    Returns
    =======

    crossed: Matrix, size(6, 6) or size(3, 3)
        This is the result of the motion cross product. It v was 6x1 the return
        is a 6x6 matrix otherwise it is a 3x3 matrix
    """


################################################################################
#
# Transformation Matrices
#
################################################################################

def X_to_Xstar():
    """This function will take a motion transformation matrix and change it to a
    force transformation matrix.

    Parameters
    ==========

    X: Matrix, size(6, 6)
        This is the transformation matrix for the motion space of spatial
        vectors

    Returns
    =======

    Xstar: Matrix, size(6, 6)
        This is the transformation matrix for the force space of spatial vectors
    """
