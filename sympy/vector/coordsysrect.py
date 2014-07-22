from sympy.core.basic import Basic
from sympy.vector.scalar import BaseScalar
from sympy.vector.functions import express, _path
from sympy import (sin, cos, eye, sympify, trigsimp,
                   ImmutableMatrix as Matrix, S, Symbol, rot_axis1,
                   rot_axis2, rot_axis3)
from sympy.core.compatibility import string_types
from sympy.core.cache import cacheit


class CoordSysCartesian(Basic):
    """
    Represents a coordinate system in 3-D space.
    """

    def __new__(cls, name, vector_names=None, variable_names=None,
                location=None, rot_type=None, rot_amounts=None,
                rot_order='', parent=None):
        """
        The orientation/location parameters are necessary if this system
        is being defined at a certain orientation or location wrt another.

        For more information on the orientation parameters, please refer
        to the docs of orient_new method.

        Parameters
        ==========

        name : str
            The name of the new CoordSysCartesian instance.

        vector_names, variable_names : tuples/lists(optional)
            Tuples/Lists of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.

        location : Vector
            The position vector of the new system's origin wrt the parent
            instance.

        rot_type : str ('Axis'/'Body'/'Quaternion'/'Space')
            The type of orientation matrix that is being created.

        rot_amounts : list OR value
            The quantities that the orientation matrix will be
            defined by.

        rot_order : str (Look at the docs of orient_new for more details)
            If applicable, the order of a series of rotations.

        parent : CoordSysCartesian
            The coordinate system wrt which the orientation/location
            (or both) is being defined.

        docstring of orient_new
        =======================

        """

        if not isinstance(name, string_types):
            raise TypeError("name should be a string")

        #If orientation information has been provided, calculate
        #the DCM accordingly
        from sympy.vector.vector import BaseVector, Vector
        if rot_type is not None:
            for i, v in enumerate(rot_amounts):
                if not isinstance(v, Vector):
                    rot_amounts[i] = sympify(v)
            rot_type = rot_type.upper()
            if rot_type == 'AXIS':
                parent_orient = _orient_axis(list(rot_amounts),
                                             rot_order, parent)
            elif rot_type == 'QUATERNION':
                parent_orient = _orient_quaternion(list(rot_amounts),
                                                   rot_order)
            elif rot_type == 'BODY':
                parent_orient = _orient_body(list(rot_amounts),
                                             rot_order)
            elif rot_type == 'SPACE':
                parent_orient = _orient_space(list(rot_amounts),
                                              rot_order)
            else:
                raise NotImplementedError('Rotation not implemented')
        else:
            if not (rot_amounts == None and rot_order == ''):
                raise ValueError("No rotation type provided")
            parent_orient = Matrix(eye(3))

        #If location information is not given, adjust the default
        #location as Vector.zero
        from sympy.vector.point import Point
        if parent is not None:
            if not isinstance(parent, CoordSysCartesian):
                raise TypeError("parent should be a " +
                                "CoordSysCartesian/None")
            if location is None:
                location = Vector.zero
            origin = parent.origin.locate_new(name + '.origin',
                                              location)
            arg_parent = parent
            arg_self = Symbol('default')
        else:
            origin = Point(name + '.origin')
            arg_parent = Symbol('default')
            arg_self = Symbol(name)

        #All systems that are defined as 'roots' are unequal, unless
        #they have the same name.
        #Systems defined at same orientation/position wrt the same
        #'parent' are equal, irrespective of the name.
        #This is true even if the same orientation is provided via
        #different methods like Axis/Body/Space/Quaternion.
        #However, coincident systems may be seen as unequal if
        #positioned/oriented wrt different parents, even though
        #they may actually be 'coincident' wrt the root system.
        obj = super(CoordSysCartesian, cls).__new__(
            cls, arg_self, parent_orient, origin, arg_parent)
        obj._name = name

        #Initialize the base vectors
        if vector_names is None:
            vector_names = (name + '.i', name + '.j', name + '.k')
        obj._i = BaseVector(vector_names[0], 0, obj)
        obj._j = BaseVector(vector_names[1], 1, obj)
        obj._k = BaseVector(vector_names[2], 2, obj)

        #Initialize the base scalars
        if variable_names is None:
            variable_names = (name + '.x', name + '.y', name + '.z')
        obj._x = BaseScalar(variable_names[0], 0, obj)
        obj._y = BaseScalar(variable_names[1], 1, obj)
        obj._z = BaseScalar(variable_names[2], 2, obj)

        #Assign a Del operator instance
        from sympy.vector.deloperator import Del
        obj._del = Del(obj)

        #Assign params
        obj._parent = parent
        if obj._parent is not None:
            obj._root = obj._parent._root
        else:
            obj._root = obj

        obj._parent_rotation_matrix = parent_orient.T
        obj._origin = origin

        #Return the instance
        return obj

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__

    def __iter__(self):
        return iter([self.i, self.j, self.k])

    @property
    def origin(self):
        return self._origin

    @property
    def delop(self):
        return self._del

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j

    @property
    def k(self):
        return self._k

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    def base_vectors(self):
        return (self._i, self._j, self._k)

    def base_scalars(self):
        return (self._x, self._y, self._z)

    @cacheit
    def rotation_matrix(self, other):
        """
        Returns the direction cosine matrix(DCM), also known as the
        'rotation matrix' of this coordinate system with respect to
        another system.

        If v_a is a vector defined in system 'A' (in matrix format)
        and v_b is the same vector defined in system 'B', then
        v_a = A.rotation_matrix(B) * v_b.

        A SymPy Matrix is returned.

        Parameters
        ==========

        other : CoordSysCartesian
            The system which the DCM is generated to.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = CoordSysCartesian('N')
        >>> A = N.orient_new('A', 'Axis', [q1, N.i])
        >>> N.rotation_matrix(A)
        Matrix([
        [1,       0,        0],
        [0, cos(q1), -sin(q1)],
        [0, sin(q1),  cos(q1)]])

        """

        if not isinstance(other, CoordSysCartesian):
            raise TypeError(str(other) +
                            " is not a CoordSysCartesian")
        #Handle special cases
        if other == self:
            return eye(3)
        elif other == self._parent:
            return self._parent_rotation_matrix
        elif other._parent == self:
            return other._parent_rotation_matrix.T
        #Else, use tree to calculate position
        rootindex, path = _path(self, other)
        result = eye(3)
        i = -1
        for i in range(rootindex):
            result *= path[i]._parent_rotation_matrix
        i += 2
        while i < len(path):
            result *= path[i]._parent_rotation_matrix.T
            i += 1
        return result

    @cacheit
    def position_wrt(self, other):
        """
        Returns the position vector of the origin of this coordinate
        system with respect to another Point/CoordSysCartesian.

        Parameters
        ==========

        other : Point/CoordSysCartesian
            If other is a Point, the position of this system's origin
            wrt it is returned. If its an instance of CoordSyRect,
            the position wrt its origin is returned.

        Examples
        ========

        >>> from sympy.vector import Point, CoordSysCartesian
        >>> N = CoordSysCartesian('N')
        >>> N1 = N.locate_new('N1', 10 * N.i)
        >>> N.position_wrt(N1)
        (-10)*N.i

        """
        return self.origin.position_wrt(other)

    def scalar_map(self, other):
        """
        Returns a dictionary which expresses the coordinate variables
        (base scalars) of this frame in terms of the variables of
        otherframe.

        Parameters
        ==========

        otherframe : CoordSysCartesian
            The other system to map the variables to.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import Symbol
        >>> A = CoordSysCartesian('A')
        >>> q = Symbol('q')
        >>> B = A.orient_new('B', 'Axis', [q, A.k])
        >>> A.scalar_map(B)
        {A.x: B.x*cos(q) - B.y*sin(q), A.y: B.x*sin(q) + B.y*cos(q), A.z: B.z}

        """

        vars_matrix = (self.rotation_matrix(other) *
                       Matrix(other.base_scalars()))
        mapping = {}
        for i, x in enumerate(self.base_scalars()):
            mapping[x] = trigsimp(vars_matrix[i])
        return mapping

    def locate_new(self, name, position, vector_names=None,
                   variable_names=None):
        """
        Returns a CoordSysCartesian with its origin located at the given
        position wrt this coordinate system's origin.

        Parameters
        ==========

        name : str
            The name of the new CoordSysCartesian instance.

        position : Vector
            The position vector of the new system's origin wrt this
            one.

        vector_names, variable_names : tuples/lists(optional)
            Tuples/Lists of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> A = CoordSysCartesian('A')
        >>> B = A.locate_new('B', 10 * A.i)
        >>> B.origin.position_wrt(A.origin)
        10*A.i

        """

        return CoordSysCartesian(name, location=position,
                                 vector_names=vector_names,
                                 variable_names=variable_names,
                                 parent=self)

    def orient_new(self, name, rot_type=None, rot_amounts=None,
                   rot_order='', location=None, vector_names=None,
                   variable_names=None):
        """
        Creates a new CoordSysCartesian oriented in the user-specified way
        with respect to this system.

        Parameters
        ==========

        name : str
            The name of the new CoordSysCartesian instance.

        rot_type : str
            The type of orientation matrix that is being created.

        rot_amounts : list OR value
            The quantities that the orientation matrix will be defined
            by.

        rot_order : str
            If applicable, the order of a series of rotations.

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        vector_names, variable_names : tuples/lists(optional)
            Tuples/Lists of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = CoordSysCartesian('N')

        We have a choice of how to implement the orientation. First is
        Body. Body orientation takes this coordinate system through three
        successive simple rotations. Hence, a 'Body' fixed rotation is
        described by three angles and three body-fixed rotation axes. To
        orient a coordinate system D with respect to N, each sequential
        rotation is always about the orthogonal unit vectors fixed to D.
        For example, a '123' rotation will specify rotations about N.i,
        then D.j, then D.k. (Initially, D.i is same as N.i)
        Therefore,

        >>> D = N.orient_new('D', 'Body', [q1, q2, q3], '123')

        is same as

        >>> D = N.orient_new('D', 'Axis', [q1, N.i])
        >>> D = D.orient_new('D', 'Axis', [q2, D.j])
        >>> D = D.orient_new('D', 'Axis', [q3, D.k])

        Acceptable rotation orders are of length 3, expressed in XYZ or
        123, and cannot have a rotation about about an axis twice in a row.

        >>> B = N.orient_new('B', 'Body', [q1, q2, q3], '123')
        >>> B = N.orient_new('B', 'Body', [q1, q2, 0], 'ZXZ')
        >>> B = N.orient_new('B', 'Body', [0, 0, 0], 'XYX')

        Body fixed rotations include both Euler Angles and
        Tait-Bryan Angles, see http://en.wikipedia.org/wiki/Euler_angles.

        Next is Space. Space is like Body, but the rotations are applied
        in the opposite order. To orient a coordinate system D with
        respect to N, each sequential rotation is always about N's
        orthogonal unit vectors. For example, a '123' rotation will
        specify rotations about N.i, then N.j, then N.k.
        Therefore,

        >>> D = N.orient_new('D', 'Space', [q1, q2, q3], '312')

        is same as

        >>> B = N.orient_new('B', 'Axis', [q1, N.i])
        >>> C = B.orient_new('C', 'Axis', [q2, N.j])
        >>> D = C.orient_new('D', 'Axis', [q3, N.k])

        Next is Quaternion. This orients the new CoordSysCartesian with
        Quaternions, defined as a finite rotation about lambda, a unit
        vector, by some amount theta.
        This orientation is described by four parameters:
        q0 = cos(theta/2)
        q1 = lambda_x sin(theta/2)
        q2 = lambda_y sin(theta/2)
        q3 = lambda_z sin(theta/2)
        Quaternion does not take in a rotation order.

        >>> B = N.orient_new('B', 'Quaternion', [q0, q1, q2, q3])

        Last is Axis. This is a rotation about an arbitrary axis by
        some angle. The axis is supplied as a Vector. This is how
        simple rotations are defined.

        >>> B = N.orient_new('B', 'Axis', [q1, N.i + 2 * N.j])

        """
        return CoordSysCartesian(name, rot_type=rot_type,
                                 rot_amounts=rot_amounts,
                                 rot_order = rot_order,
                                 vector_names=vector_names,
                                 variable_names=variable_names,
                                 location = location,
                                 parent=self)
    __new__.__doc__ += orient_new.__doc__

    def orient_new_axis(self, name, angle, axis, location=None):
        """
        Axis rotation is a rotation about an arbitrary axis by
        some angle. The angle is supplied as a SymPy expr scalar, and
        the axis is supplied as a Vector.

        Parameters
        ==========

        name : string
            The name of the new coordinate system

        angle : Expr
            The angle by which the new system is to be rotated

        axis : Vector
            The axis around which the rotation has to be performed

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = CoordSysCartesian('N')
        >>> B = N.orient_new_axis('B', q1, N.i + 2 * N.j)

        """

        return self.orient_new(name, 'Axis', [angle, axis],
                               location=location)

    def orient_new_body(self, name, angle1, angle2, angle3,
                        rotation_order, location=None):
        """
        Body orientation takes this coordinate system through three
        successive simple rotations.

        Body fixed rotations include both Euler Angles and
        Tait-Bryan Angles, see http://en.wikipedia.org/wiki/Euler_angles.

        Parameters
        ==========

        name : string
            The name of the new coordinate system

        angle1, angle2, angle3 : Expr
            Three successive angles to rotate the coordinate system by

        rotation_order : string
            String defining the order of axes for rotation

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q1, q2, q3 = symbols('q1 q2 q3')
        >>> N = CoordSysCartesian('N')

        A 'Body' fixed rotation is described by three angles and
        three body-fixed rotation axes. To orient a coordinate system D
        with respect to N, each sequential rotation is always about
        the orthogonal unit vectors fixed to D. For example, a '123'
        rotation will specify rotations about N.i, then D.j, then
        D.k. (Initially, D.i is same as N.i)
        Therefore,

        >>> D = N.orient_new_body('D', q1, q2, q3, '123')

        is same as

        >>> D = N.orient_new_axis('D', q1, N.i)
        >>> D = D.orient_new_axis('D', q2, D.j)
        >>> D = D.orient_new_axis('D', q3, D.k)

        Acceptable rotation orders are of length 3, expressed in XYZ or
        123, and cannot have a rotation about about an axis twice in a row.

        >>> B = N.orient_new_body('B', q1, q2, q3, '123')
        >>> B = N.orient_new_body('B', q1, q2, 0, 'ZXZ')
        >>> B = N.orient_new_body('B', 0, 0, 0, 'XYX')

        """

        return self.orient_new(name, 'Body', [angle1, angle2, angle3],
                               rotation_order, location=location)

    def orient_new_space(self, name, angle1, angle2, angle3,
                         rotation_order, location=None):
        """
        Space rotation is similar to Body rotation, but the rotations
        are applied in the opposite order.

        Refer to the docs of orient_new_body for more information and
        examples.

        Parameters
        ==========

        name : string
            The name of the new coordinate system

        angle1, angle2, angle3 : Expr
            Three successive angles to rotate the coordinate system by

        rotation_order : string
            String defining the order of axes for rotation

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q1, q2, q3 = symbols('q1 q2 q3')
        >>> N = CoordSysCartesian('N')

        To orient a coordinate system D with respect to N, each
        sequential rotation is always about N's orthogonal unit vectors.
        For example, a '123' rotation will specify rotations about
        N.i, then N.j, then N.k.
        Therefore,

        >>> D = N.orient_new_space('D', q1, q2, q3, '312')

        is same as

        >>> B = N.orient_new_axis('B', q1, N.i)
        >>> C = B.orient_new_axis('C', q2, N.j)
        >>> D = C.orient_new_axis('D', q3, N.k)

        docs of orient_new_body
        =======================

        """

        return self.orient_new(name, 'Space', [angle1, angle2, angle3],
                               rotation_order, location=location)

    def orient_new_quaternion(self, name, q0, q1, q2, q3, location=None):
        """
        Quaternion orientation orients the new CoordSysCartesian with
        Quaternions, defined as a finite rotation about lambda, a unit
        vector, by some amount theta.
        This orientation is described by four parameters:
        q0 = cos(theta/2)
        q1 = lambda_x sin(theta/2)
        q2 = lambda_y sin(theta/2)
        q3 = lambda_z sin(theta/2)
        Quaternion does not take in a rotation order.

        Parameters
        ==========

        name : string
            The name of the new coordinate system

        q0, q1, q2, q3 : Expr
            The quaternions to rotate the coordinate system by

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = CoordSysCartesian('N')
        >>> B = N.orient_new('B', 'Quaternion', [q0, q1, q2, q3])

        """

        return self.orient_new(name, 'Quaternion', [q0, q1, q2, q3],
                               location = location)

    def __init__(self, name, vector_names=None, variable_names=None,
                location=None, rot_type=None, rot_amounts=None,
                rot_order='', parent=None):
        #Dummy initializer for setting docstring
        pass
    __init__.__doc__ = __new__.__doc__


def _rot(axis, angle):
    """DCM for simple axis 1, 2 or 3 rotations. """
    if axis == 1:
        return Matrix(rot_axis1(angle).T)
    elif axis == 2:
        return Matrix(rot_axis2(angle).T)
    elif axis == 3:
        return Matrix(rot_axis3(angle).T)


def _orient_axis(amounts, rot_order, parent):
    """
    Helper method for orientation using Axis method.
    """
    if not rot_order == '':
        raise TypeError('Axis orientation takes no' +
                        'rotation order')
    if not (isinstance(amounts, (list, tuple))
                and (len(amounts) == 2)):
        raise TypeError('Amounts should be of length 2')
    theta = amounts[0]
    axis = amounts[1]
    axis = express(axis, parent).normalize()
    axis = axis.to_matrix(parent)
    parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
            Matrix([[0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                [-axis[1], axis[0], 0]]) * sin(theta) +
                     axis * axis.T)
    return parent_orient


def _orient_quaternion(amounts, rot_order):
    """
    Helper method for orientation using Quaternion method.
    """
    if not rot_order == '':
        raise TypeError(
            'Quaternion orientation takes no rotation order')
    if not (isinstance(amounts, (list, tuple)) &
                (len(amounts) == 4)):
        raise TypeError('Amounts should be of length 4')
    q0, q1, q2, q3 = amounts
    parent_orient = (Matrix([[q0 ** 2 + q1 ** 2 - q2 ** 2 - \
                              q3 **
        2, 2 * (q1 * q2 - q0 * q3), 2 * (q0 * q2 + q1 * q3)],
        [2 * (q1 * q2 + q0 * q3), \
         q0 ** 2 - q1 ** 2 + q2 ** 2 - q3 ** 2,
        2 * (q2 * q3 - q0 * q1)], [2 * (q1 * q3 - q0 * q2), \
                                   2 * (q0 * q1 + q2 * q3), \
                                   q0 ** 2 - q1 ** 2 - \
                                   q2 ** 2 + q3 ** 2]]))
    return parent_orient


def _orient_body(amounts, rot_order):
    """
    Helper method for orientation using Body method.
    """
    approved_orders = ('123', '231', '312', '132', '213',
                       '321', '121', '131', '212', '232',
                       '313', '323', '')
    rot_order = str(
        rot_order).upper()
    if not (len(amounts) == 3 & len(rot_order) == 3):
        raise TypeError('Body orientation takes 3' +
                        'values & 3 orders')
    rot_order = [i.replace('X', '1') for i in rot_order]
    rot_order = [i.replace('Y', '2') for i in rot_order]
    rot_order = [i.replace('Z', '3') for i in rot_order]
    rot_order = ''.join(rot_order)
    if not rot_order in approved_orders:
        raise TypeError('Invalid rot_type parameter')
    a1 = int(rot_order[0])
    a2 = int(rot_order[1])
    a3 = int(rot_order[2])
    parent_orient = (_rot(a1, amounts[0]) *
                     _rot(a2, amounts[1]) *
                     _rot(a3, amounts[2]))
    return parent_orient


def _orient_space(amounts, rot_order):
    """
    Helper method for orientation using Space method.
    """
    approved_orders = ('123', '231', '312', '132', '213',
                       '321', '121', '131', '212', '232',
                       '313', '323', '')
    rot_order = str(
        rot_order).upper()
    if not (len(amounts) == 3 & len(rot_order) == 3):
        raise TypeError('Space orientation takes 3 ' +
                        'values & 3 orders')
    rot_order = [i.replace('X', '1') for i in rot_order]
    rot_order = [i.replace('Y', '2') for i in rot_order]
    rot_order = [i.replace('Z', '3') for i in rot_order]
    rot_order = ''.join(rot_order)
    if not rot_order in approved_orders:
        raise TypeError('Invalid rot_type parameter')
    a1 = int(rot_order[0])
    a2 = int(rot_order[1])
    a3 = int(rot_order[2])
    parent_orient = (_rot(a3, amounts[2]) *
                     _rot(a2, amounts[1]) *
                     _rot(a1, amounts[0]))
    return parent_orient
