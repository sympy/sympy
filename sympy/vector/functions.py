from sympy.vector.scalar import BaseScalar
from sympy.vector import CoordSysCartesian
from sympy.vector.dyadic import Dyadic
from sympy import sympify


def express(expr, system, system2=None, variables=False):
    """
    Global function for 'express' functionality.

    Re-expresses a Vector, Dyadic or scalar(sympyfiable) in the given
    coordinate system.

    If 'variables' is True, then the coordinate variables (base scalars)
    of other coordinate systems present in the vector/scalar field or dyadic
    are also substituted in terms of the base scalars of the given system.

    Parameters
    ==========

    expr : Vector/Dyadic/scalar(sympyfiable)
        The expression to re-express in ReferenceFrame 'frame'

    system: CoordSysCartesian
        The coordinate system the expr is to be expressed in

    variables : boolean
        Specifies whether to substitute the coordinate variables present
        in expr, in terms of those of parameter system

    Examples
    ========

    >>> from sympy.vector import CoordSysCartesian
    >>> from sympy import Symbol
    >>> N = CoordSysCartesian('N')
    >>> q = Symbol('q')
    >>> B = N.orient_new('B', 'Axis', [q, N.k])
    >>> from sympy.vector import express
    >>> express(B.i, N)
    (cos(q))*N.i + (sin(q))*N.j
    >>> express(N.x, B, variables=True)
    B.x*cos(q) - B.y*sin(q)
    >>> d = N.i.outer(N.i)
    >>> express(d, B, N)
    cos(q)*(B.i|N.i) - sin(q)*(B.j|N.i)

    """

    from sympy.vector.vector import Vector, BaseVector
    if expr == 0 or expr == Vector.zero:
        return expr

    if not isinstance(system, CoordSysCartesian):
        raise TypeError("system should be a CoordSysCartesian \
                        instance")

    if isinstance(expr, Vector):
        if system2 is not None:
            raise ValueError("system2 should not be provided for \
                                Vectors")
        #Given expr is a Vector
        if variables:
            #If variables attribute is True, substitute
            #the coordinate variables in the Vector
            system_list = []
            for x in expr.atoms():
                if (isinstance(x, BaseScalar) or \
                    isinstance(x, BaseVector)) and \
                    x.system != system:
                    system_list.append(x.system)
            system_list = set(system_list)
            subs_dict = {}
            for f in system_list:
                subs_dict.update(f.scalar_map(system))
            expr = expr.subs(subs_dict)
        #Re-express in this frame
        outvec = Vector.zero
        parts = expr.separate()
        for x in parts:
            if x != system:
                temp = system.rotation_matrix(x) * parts[x].to_matrix(x)
                outvec += matrix_to_vector(temp, system)
            else:
                outvec += parts[x]
        return outvec

    elif isinstance(other, Dyadic):
        if system2 is None:
            system2 = system
        if not isinstance(system2, CoordSysCartesian):
            raise TypeError("system2 shoule be a CoordSysCartesian \
                            instance")
        outdyad = Dyadic.zero
        var = variables
        for k, v in expr.components.items():
            outdyad += express(v, system, variables = var) * \
                       (express(k.args[0], system, variables = var) |
                        express(k.args[1], system2, variables = var))

        return outdyad

    else:
        if system2 is not None:
            raise ValueError("system2 should not be provided for \
                                Vectors")
        if variables:
            #Given expr is a scalar field
            system_set = set([])
            expr = sympify(expr)
            #Subsitute all the coordinate variables
            for x in expr.atoms():
                if isinstance(x, BaseScalar)and x.system != system:
                    system_set.add(x.system)
            subs_dict = {}
            for f in system_set:
                subs_dict.update(f.scalar_map(system))
            return expr.subs(subs_dict)
        return expr


def matrix_to_vector(matrix, system):
    """
    Converts a vector in matrix form to a Vector instance.

    It is assumed that the elements of the Matrix represent the
    measure numbers of the components of the vector along basis
    vectors of 'system'.

    Parameters
    ==========

    matrix : SymPy Matrix, Dimensions: (1, 3)
        The matrix to be converted to a vector

    system : CoordSysCartesian
        The coordinate system the vector is to be defined in

    Examples
    ========

    >>> from sympy import ImmutableMatrix as Matrix
    >>> m = Matrix([1, 2, 3])
    >>> from sympy.vector import CoordSysCartesian, matrix_to_vector
    >>> C = CoordSysCartesian('C')
    >>> v = matrix_to_vector(m, C)
    >>> v
    C.i + 2*C.j + 3*C.k
    >>> v.to_matrix(C) == m
    True

    """

    from sympy.vector.vector import Vector
    outvec = Vector.zero
    vects = system.base_vectors()
    for i, x in enumerate(matrix):
        outvec += x * vects[i]
    return outvec


def _path(from_object, to_object):
    """
    Calculates the 'path' of objects starting from 'from_object'
    to 'to_object', along with the index of the first common
    ancestor in the tree.

    Returns (index, list) tuple.
    """

    if from_object._root != to_object._root:
        raise ValueError("No connecting path found between " +
                         str(from_object) + " and " + str(to_object))

    other_path = []
    obj = to_object
    while obj._parent is not None:
        other_path.append(obj)
        obj = obj._parent
    other_path.append(obj)
    object_set = set(other_path)
    from_path = []
    obj = from_object
    while obj not in object_set:
        from_path.append(obj)
        obj = obj._parent
    index = len(from_path)
    i = other_path.index(obj)
    while i >= 0:
        from_path.append(other_path[i])
        i -= 1
    return index, from_path
