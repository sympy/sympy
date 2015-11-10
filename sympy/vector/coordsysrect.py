from sympy.core.basic import Basic
from sympy.vector.scalar import BaseScalar
from sympy import eye, trigsimp, ImmutableMatrix as Matrix, Symbol
from sympy.core.compatibility import string_types, range
from sympy.core.cache import cacheit
from sympy.vector.orienters import (Orienter, AxisOrienter, BodyOrienter,
                                    SpaceOrienter, QuaternionOrienter)
import sympy.vector


class CoordSysCartesian(Basic):
    """
    Represents a coordinate system in 3-D space.
    """

    def __new__(cls, name, location=None, rotation_matrix=None,
                parent=None, vector_names=None, variable_names=None):
        """
        The orientation/location parameters are necessary if this system
        is being defined at a certain orientation or location wrt another.

        Parameters
        ==========

        name : str
            The name of the new CoordSysCartesian instance.

        location : Vector
            The position vector of the new system's origin wrt the parent
            instance.

        rotation_matrix : SymPy ImmutableMatrix
            The rotation matrix of the new coordinate system with respect
            to the parent. In other words, the output of
            new_system.rotation_matrix(parent).

        parent : CoordSysCartesian
            The coordinate system wrt which the orientation/location
            (or both) is being defined.

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        """

        name = str(name)
        Vector = sympy.vector.Vector
        BaseVector = sympy.vector.BaseVector
        Point = sympy.vector.Point
        if not isinstance(name, string_types):
            raise TypeError("name should be a string")

        #If orientation information has been provided, store
        #the rotation matrix accordingly
        if rotation_matrix is None:
            parent_orient = Matrix(eye(3))
        else:
            if not isinstance(rotation_matrix, Matrix):
                raise TypeError("rotation_matrix should be an Immutable" +
                                "Matrix instance")
            parent_orient = rotation_matrix

        #If location information is not given, adjust the default
        #location as Vector.zero
        if parent is not None:
            if not isinstance(parent, CoordSysCartesian):
                raise TypeError("parent should be a " +
                                "CoordSysCartesian/None")
            if location is None:
                location = Vector.zero
            else:
                #Check that location does not contain base
                #scalars
                for x in location.free_symbols:
                    if isinstance(x, BaseScalar):
                        raise ValueError("location should not contain" +
                                         " BaseScalars")
            origin = parent.origin.locate_new(name + '.origin',
                                              location)
        else:
            origin = Point(name + '.origin')

        #Parent_arg to be passed to superclass
        if parent is not None:
            arg_parent = parent
        else:
            location = Vector.zero
            arg_parent = None

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
            cls, Symbol(name), location, parent_orient, arg_parent)
        obj._name = name

        #Initialize the base vectors
        if vector_names is None:
            vector_names = (name + '.i', name + '.j', name + '.k')
            latex_vects = [(r'\mathbf{\hat{i}_{%s}}' % name),
                           (r'\mathbf{\hat{j}_{%s}}' % name),
                           (r'\mathbf{\hat{k}_{%s}}' % name)]
            pretty_vects = (name + '_i', name + '_j', name + '_k')
        else:
            _check_strings('vector_names', vector_names)
            vector_names = list(vector_names)
            latex_vects = [(r'\mathbf{\hat{%s}_{%s}}' % (x, name)) for
                           x in vector_names]
            pretty_vects = [(name + '_' + x) for x in vector_names]

        obj._i = BaseVector(vector_names[0], 0, obj,
                            pretty_vects[0], latex_vects[0])
        obj._j = BaseVector(vector_names[1], 1, obj,
                            pretty_vects[1], latex_vects[1])
        obj._k = BaseVector(vector_names[2], 2, obj,
                            pretty_vects[2], latex_vects[2])

        #Initialize the base scalars
        if variable_names is None:
            variable_names = (name + '.x', name + '.y', name + '.z')
            latex_scalars = [(r"\mathbf{{x}_{%s}}" % name),
                             (r"\mathbf{{y}_{%s}}" % name),
                             (r"\mathbf{{z}_{%s}}" % name)]
            pretty_scalars = (name + '_x', name + '_y', name + '_z')
        else:
            _check_strings('variable_names', vector_names)
            variable_names = list(variable_names)
            latex_scalars = [(r"\mathbf{{%s}_{%s}}" % (x, name)) for
                           x in variable_names]
            pretty_scalars = [(name + '_' + x) for x in variable_names]

        obj._x = BaseScalar(variable_names[0], 0, obj,
                            pretty_scalars[0], latex_scalars[0])
        obj._y = BaseScalar(variable_names[1], 1, obj,
                            pretty_scalars[1], latex_scalars[1])
        obj._z = BaseScalar(variable_names[2], 2, obj,
                            pretty_scalars[2], latex_scalars[2])

        #Assign a Del operator instance
        from sympy.vector.deloperator import Del
        obj._delop = Del(obj)

        #Assign params
        obj._parent = parent
        if obj._parent is not None:
            obj._root = obj._parent._root
        else:
            obj._root = obj

        obj._parent_rotation_matrix = parent_orient
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
        return self._delop

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
        >>> A = N.orient_new_axis('A', q1, N.i)
        >>> N.rotation_matrix(A)
        Matrix([
        [1,       0,        0],
        [0, cos(q1), -sin(q1)],
        [0, sin(q1),  cos(q1)]])

        """

        from sympy.vector.functions import _path
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
        >>> B = A.orient_new_axis('B', q, A.k)
        >>> A.scalar_map(B)
        {A.x: B.x*cos(q) - B.y*sin(q), A.y: B.x*sin(q) + B.y*cos(q), A.z: B.z}

        """

        relocated_scalars = []
        origin_coords = tuple(self.position_wrt(other).to_matrix(other))
        for i, x in enumerate(other.base_scalars()):
            relocated_scalars.append(x - origin_coords[i])

        vars_matrix = (self.rotation_matrix(other) *
                       Matrix(relocated_scalars))
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

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

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

    def orient_new(self, name, orienters, location=None,
                   vector_names=None, variable_names=None):
        """
        Creates a new CoordSysCartesian oriented in the user-specified way
        with respect to this system.

        Please refer to the documentation of the orienter classes
        for more information about the orientation procedure.

        Parameters
        ==========

        name : str
            The name of the new CoordSysCartesian instance.

        orienters : iterable/Orienter
            An Orienter or an iterable of Orienters for orienting the
            new coordinate system.
            If an Orienter is provided, it is applied to get the new
            system.
            If an iterable is provided, the orienters will be applied
            in the order in which they appear in the iterable.

        location : Vector(optional)
            The location of the new coordinate system's origin wrt this
            system's origin. If not specified, the origins are taken to
            be coincident.

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = CoordSysCartesian('N')

        Using an AxisOrienter

        >>> from sympy.vector import AxisOrienter
        >>> axis_orienter = AxisOrienter(q1, N.i + 2 * N.j)
        >>> A = N.orient_new('A', (axis_orienter, ))

        Using a BodyOrienter

        >>> from sympy.vector import BodyOrienter
        >>> body_orienter = BodyOrienter(q1, q2, q3, '123')
        >>> B = N.orient_new('B', (body_orienter, ))

        Using a SpaceOrienter

        >>> from sympy.vector import SpaceOrienter
        >>> space_orienter = SpaceOrienter(q1, q2, q3, '312')
        >>> C = N.orient_new('C', (space_orienter, ))

        Using a QuaternionOrienter

        >>> from sympy.vector import QuaternionOrienter
        >>> q_orienter = QuaternionOrienter(q0, q1, q2, q3)
        >>> D = N.orient_new('D', (q_orienter, ))

        """

        if isinstance(orienters, Orienter):
            if isinstance(orienters, AxisOrienter):
                final_matrix = orienters.rotation_matrix(self)
            else:
                final_matrix = orienters.rotation_matrix()
        else:
            final_matrix = Matrix(eye(3))
            for orienter in orienters:
                if isinstance(orienter, AxisOrienter):
                    final_matrix *= orienter.rotation_matrix(self)
                else:
                    final_matrix *= orienter.rotation_matrix()

        return CoordSysCartesian(name, rotation_matrix=final_matrix,
                                 vector_names=vector_names,
                                 variable_names=variable_names,
                                 location = location,
                                 parent=self)

    def orient_new_axis(self, name, angle, axis, location=None,
                        vector_names=None, variable_names=None):
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

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = CoordSysCartesian('N')
        >>> B = N.orient_new_axis('B', q1, N.i + 2 * N.j)

        """

        orienter = AxisOrienter(angle, axis)
        return self.orient_new(name, orienter,
                               location=location,
                               vector_names=vector_names,
                               variable_names=variable_names)

    def orient_new_body(self, name, angle1, angle2, angle3,
                        rotation_order, location=None,
                        vector_names=None, variable_names=None):
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

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

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

        orienter = BodyOrienter(angle1, angle2, angle3, rotation_order)
        return self.orient_new(name, orienter,
                               location=location,
                               vector_names=vector_names,
                               variable_names=variable_names)

    def orient_new_space(self, name, angle1, angle2, angle3,
                         rotation_order, location=None,
                         vector_names=None, variable_names=None):
        """
        Space rotation is similar to Body rotation, but the rotations
        are applied in the opposite order.

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

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        See Also
        ========

        CoordSysCartesian.orient_new_body : method to orient via Euler
            angles

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

        """

        orienter = SpaceOrienter(angle1, angle2, angle3, rotation_order)
        return self.orient_new(name, orienter,
                               location=location,
                               vector_names=vector_names,
                               variable_names=variable_names)

    def orient_new_quaternion(self, name, q0, q1, q2, q3, location=None,
                              vector_names=None, variable_names=None):
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

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> from sympy import symbols
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = CoordSysCartesian('N')
        >>> B = N.orient_new_quaternion('B', q0, q1, q2, q3)

        """

        orienter = QuaternionOrienter(q0, q1, q2, q3)
        return self.orient_new(name, orienter,
                               location=location,
                               vector_names=vector_names,
                               variable_names=variable_names)

    def __init__(self, name, location=None, rotation_matrix=None,
                 parent=None, vector_names=None, variable_names=None,
                 latex_vects=None, pretty_vects=None, latex_scalars=None,
                 pretty_scalars=None):
        #Dummy initializer for setting docstring
        pass
    __init__.__doc__ = __new__.__doc__


def _check_strings(arg_name, arg):
    errorstr = arg_name + " must be an iterable of 3 string-types"
    if len(arg) != 3:
        raise ValueError(errorstr)
    try:
        for s in arg:
            if not isinstance(s, string_types):
                raise TypeError(errorstr)
    except:
        raise TypeError(errorstr)
