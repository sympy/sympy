from sympy.core.basic import Basic
from sympy.vector.scalar import BaseScalar
from sympy.vector.functions import express, _path
from sympy import sin, cos, eye, sympify, trigsimp, \
     ImmutableMatrix as Matrix
from sympy.core.compatibility import string_types


class CoordSysRect(Basic):
    """
    Represents a coordinate system in 3-D space.
    """

    def __new__(cls, name, vector_names=None, variable_names=None, \
                location=None, rot_type=None, rot_amounts=None, \
                rot_order='', parent=None):

        #A CoordSysRect can be uniquely identified by its name
        obj = super(CoordSysRect, cls).__new__(cls, name)
        obj._name = name

        #Initialize the base vectors
        from sympy.vector.vector import BaseVector, Vector
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

        #Assign important params
        from sympy.vector.point import Point
        if not isinstance(name, string_types):
            raise TypeError("name should be a string")
        obj._parent = parent
        if obj._parent is not None:
            if not isinstance(parent, CoordSysRect):
                raise TypeError("parent should be a CoordSysRect/None")
            obj._root = obj._parent._root
            obj._origin = parent.origin.locatenew(name + '.origin', \
                                                  location)
        else:
            obj._origin = Point(name + '.origin')
            obj._root = obj

        #If orientation information has been provided, orient the
        #coordinate system accordingly
        if rot_type is not None:
            amounts = list(rot_amounts)
            for i, v in enumerate(amounts):
                if not isinstance(v, Vector):
                    amounts[i] = sympify(v)
            def _rot(axis, angle):
                """DCM for simple axis 1, 2 or 3 rotations. """
                if axis == 1:
                    return Matrix([[1, 0, 0],
                        [0, cos(angle), -sin(angle)],
                        [0, sin(angle), cos(angle)]])
                elif axis == 2:
                    return Matrix([[cos(angle), 0, sin(angle)],
                        [0, 1, 0],
                        [-sin(angle), 0, cos(angle)]])
                elif axis == 3:
                    return Matrix([[cos(angle), -sin(angle), 0],
                        [sin(angle), cos(angle), 0],
                        [0, 0, 1]])

            approved_orders = ('123', '231', '312', '132', '213', \
                               '321', '121', '131', '212', '232', \
                               '313', '323', '')
            rot_order = str(
                rot_order).upper()  # Now we need to make sure XYZ = 123
            rot_type = rot_type.upper()
            rot_order = [i.replace('X', '1') for i in rot_order]
            rot_order = [i.replace('Y', '2') for i in rot_order]
            rot_order = [i.replace('Z', '3') for i in rot_order]
            rot_order = ''.join(rot_order)
            if not rot_order in approved_orders:
                raise TypeError('Invalid rot_type parameter')
            parent_orient = []
            if rot_type == 'AXIS':
                if not rot_order == '':
                    raise TypeError('Axis orientation takes no' +  \
                                    'rotation order')
                if not (isinstance(amounts, (list, tuple)) & \
                        (len(amounts) == 2)):
                    raise TypeError('Amounts should be of length 2')
                theta = amounts[0]
                axis = amounts[1]
                axis = express(axis, parent).normalize()
                axis = axis.to_matrix(parent)
                parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                        Matrix([[0, -axis[2], axis[1]], \
                                [axis[2], 0, -axis[0]],
                            [-axis[1], axis[0], 0]]) * sin(theta) + \
                                 axis * axis.T)
            elif rot_type == 'QUATERNION':
                if not rot_order == '':
                    raise TypeError(
                        'Quaternion orientation takes no rotation order')
                if not (isinstance(amounts, (list, tuple)) & \
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
            elif rot_type == 'BODY':
                if not (len(amounts) == 3 & len(rot_order) == 3):
                    raise TypeError('Body orientation takes 3' + \
                                    'values & 3 orders')
                a1 = int(rot_order[0])
                a2 = int(rot_order[1])
                a3 = int(rot_order[2])
                parent_orient = (_rot(a1, amounts[0]) *
                                 _rot(a2, amounts[1]) *
                                 _rot(a3, amounts[2]))
            elif rot_type == 'SPACE':
                if not (len(amounts) == 3 & len(rot_order) == 3):
                    raise TypeError('Space orientation takes 3 ' + \
                                    'values & 3 orders')
                a1 = int(rot_order[0])
                a2 = int(rot_order[1])
                a3 = int(rot_order[2])
                parent_orient = (_rot(a3, amounts[2]) * \
                                 _rot(a2, amounts[1]) * \
                                 _rot(a1, amounts[0]))
            else:
                raise NotImplementedError('Rotation not implemented')
        else:
            parent_orient = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        obj._parent_dcm = parent_orient

        #If location information is given, adjust the new instance's
        #origin accordingly
        if location is None:
            location = Vector.Zero

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

    def dcm(self, other):
        """
        Returns the direction cosine matrix(DCM) of this coordinate
        system with respect to another system.

        A SymPy Matrix is returned.

        Parameters
        ==========

        other : CoordSysRect
            The system which the DCM is generated to.

        Examples
        ========

        >>> from sympy.vector import CoordSysRect
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = CoordSysRect('N')
        >>> A = N.orient_new('A', 'Axis', [q1, N.i])
        >>> N.dcm(A)
        Matrix([
        [1,       0,        0],
        [0, cos(q1), -sin(q1)],
        [0, sin(q1),  cos(q1)]])

        """

        if not isinstance(other, CoordSysRect):
            raise TypeError(str(other) + \
                            " is not a CoordSysRect")
        #Handle special cases
        if other == self:
            return eye(3)
        elif other == self._parent:
            return self._parent_dcm
        elif other._parent == self:
            return other._parent_dcm
        #Else, use tree to calculate position
        rootindex, path = self._path(other)
        result = eye(3)
        i = -1
        for i in range(rootindex):
            result *= path[i]._parent_dcm
        i += 2
        while i < len(path):
            result *= path[i]._parent_dcm.T
            i += 1
        return result

    def variable_map(self, other, simplify=False):
        """
        Returns a dictionary which expresses the coordinate variables
        (base scalars) of this frame in terms of the variables of
        otherframe.

        If 'simplify' is True, returns a simplified version of the mapped
        values. Else, returns them without simplification.

        Simplification of the expressions may take time.

        Parameters
        ==========

        otherframe : ReferenceFrame
            The other frame to map the variables to.

        simplify : bool
            Specifies whether simplification of the expressions is
            needed.

        Examples
        ========

        >>> from sympy.physics.vector import CoordSysRect
        >>> from sympy import Symbol
        >>> A = CoordSysRect('A')
        >>> q = Symbol('q')
        >>> B = A.orient_new('B', 'Axis', [q, A.z])
        >>> A.variable_map(B)
        {A.x: B.x*cos(q(t)) - B.y*sin(q(t)), A.y: B.x*sin(q(t)) + B.y*cos(q(t)), A.z: B.z}

        """

        vars_matrix = self.dcm(other) * Matrix(other.base_scalars())
        mapping = {}
        for i, x in enumerate(self.base_scalars()):
            if simplify:
                mapping[x] = trigsimp(vars_matrix[i], method='fu')
            else:
                mapping[x] = vars_matrix[i]
        return mapping

    def locate_new(self, name, position, vector_names=None, \
                   variable_names=None):
        """
        Returns a CoordSysRect with its origin located at the given
        position wrt this coordinate system's origin.

        Parameters
        ==========

        name : str
            The name of the new CoordSysRect instance

        position : Vector
            The position vector of the new system's origin wrt this
            one

        vector_names, variable_names : tuples(optional)
            Tuples of 3 strings each, with custom names for base vectors
            and base scalars of the new system respectively

        Examples
        ========

        >>> from sympy.vector import CoordSysRect
        >>> A = CoordSysRect('A')
        >>> B = A.locate_new('B', 10 * A.i)
        >>> B.origin.position_wrt(A.origin)
        10*A.i

        """

        return CoordSysRect(name, location=position, \
                            vector_names=vector_names, \
                            variable_names=variable_names, parent=self)

    def orient_new(self, vector_names=None, variable_names=None, \
                   rot_type=None, rot_amounts=None, rot_order=''):
        """
        Creates a new CoordSysRect oriented in the user-specified way
        with respect to this system.



        Parameters
        ==========

        newname : str
            The name for the new ReferenceFrame
        rot_type : str
            The type of orientation matrix that is being created.
        amounts : list OR value
            The quantities that the orientation matrix will be defined by.
        rot_order : str
            If applicable, the order of a series of rotations.


        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])

        """
        return None
