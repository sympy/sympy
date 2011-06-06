__all__ = ['ReferenceFrame', 'Vector']

from sympy import Matrix, Symbol, sin, cos, eye, simplify, diff, sqrt

class ReferenceFrame(object):
    """
    ReferenceFrame is a reference frame.
    It will store its basis vectors as attributes, and orientation
    information to a parent frame, or it will be at the top of a tree.
    """

    def __init__(self, name=''):
        """init for ReferenceFrame. """
        self.name = name
        self.parent = None
        self._ang_vel = None
        self._ang_vel_parent = None
        self._cur = 0
        self._x = Vector([(Matrix([1, 0, 0]), self)])
        self._y = Vector([(Matrix([0, 1, 0]), self)])
        self._z = Vector([(Matrix([0, 0, 1]), self)])

    def __iter__(self):
        """Allows iteration through the basis vectors. """
        return self

    def __str__(self):
        """Returns the name of the frame. """
        return self.name

    def __repr__(self):
        """Returns the name of the frame. """
        return self.name

    def next(self):
        """Iterator for the basis vectors. """
        if self._cur == 0:
            self._cur += 1
            return self._x
        elif self._cur == 1:
            self._cur += 1
            return self._y
        elif self._cur == 2:
            self._cur += 1
            return self._z
        else:
            self._cur = 0
            raise StopIteration

    def _common_frame(self,other):
        """This returns the first common parent between two ReferenceFrames.
        Takes in another ReferenceFrame, and returns a ReferenceFrame.
        """
        if not isinstance(other, ReferenceFrame):
            raise TypeError('You have to use a ReferenceFrame')
        leg1 = [self]
        ptr = self
        while ptr.parent != None:
            ptr = ptr.parent
            leg1.append(ptr)
        leg2 = [other]
        ptr = other
        while ptr.parent != None:
            ptr = ptr.parent
            leg2.append(ptr)
        for i,v1 in enumerate(leg1):
            for j, v2 in enumerate(leg2):
                if v1 == v2:
                    return v1
        raise ValueError('No Common Frame')

    def ang_vel_in(self, otherframe):
        """Returns the angular velocity vector of the ReferenceFrame.


        Parameters
        ==========
        otherframe : ReferenceFrame
            The ReferenceFrame which the angular velocity is returned in.

        Examples
        ========

        >>> from sympy.physics.classical.essential import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(V, N)
        >>> A.ang_vel_in(N)
        (10)*nx>

        """
        if otherframe == self._ang_vel_parent:
            return self._ang_vel
        commonframe = self._common_frame(otherframe)
        # form DCM from self to first common frame
        leg1 = 0
        ptr = self
        while ptr != commonframe:
            leg1 += ptr._ang_vel
            ptr = ptr.parent
        # form DCM from other to first common frame
        leg2 = 0
        ptr = otherframe
        while ptr != commonframe:
            leg2 -= ptr._ang_vel
            ptr = ptr.parent
        return leg1 + leg2

    def dcm(self, otherframe):
        """The direction cosine matrix between frames.

        This gives the DCM between this frame and the otherframe.
        The format is N.xyz = N.dcm(B) * B.xyz
        A SymPy Matrix is returned.

        Parameters
        ==========
        otherframe : ReferenceFrame
            The otherframe which the DCM is generated to.

        Examples
        ========

        >>> from sympy.physics.classical.essential import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 1)
        >>> N.dcm(A)
        [1,       0,        0]
        [0, cos(q1), -sin(q1)]
        [0, sin(q1),  cos(q1)]

        """

        commonframe = self._common_frame(otherframe)
        # form DCM from self to first common frame
        leg1 = eye(3)
        ptr = self
        while ptr != commonframe:
            leg1 = ptr.parent_orient * leg1
            ptr = ptr.parent
        # form DCM from other to first common frame
        leg2 = eye(3)
        ptr = otherframe
        while ptr != commonframe:
            leg2 = ptr.parent_orient * leg2
            ptr = ptr.parent
        return leg2 * leg1.T

    def orientnew(self, newname, rot_type, amounts, rot_order=''):
        """Creates a new ReferenceFrame oriented with respect to this Frame.

        See ReferenceFrame.orient() for acceptable rotation types, amounts,
        and orders. Parent is going to be self.

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

        >>> from sympy.physics.classical.essential import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 1)

        """

        newframe = ReferenceFrame(newname)
        newframe.orient(self, rot_type, amounts, rot_order)
        return newframe

    def orient(self, parent, rot_type, amounts, rot_order=''):
        """Defines the orientation of this frame relative to a parent frame.

        Supported orientation types are Simple, Body, Space, Euler, Axis.
        Examples show correct usage.

        Parameters
        ==========
        parent : ReferenceFrame
            The frame that this ReferenceFrame will have its orientation matrix
            defined in relation to.
        rot_type : str
            The type of orientation matrix that is being created.
        amounts : list OR value
            The quantities that the orientation matrix will be defined by.
        rot_order : str
            If applicable, the order of a series of rotations.

        Examples
        ========

        >>> from sympy.physics.classical.essential import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1, q2, q3, q4 = symbols('q1 q2 q3 q4')
        >>> N = ReferenceFrame('N')
        >>> B = ReferenceFrame('B')

        Now we have a choice of how to implement the orientation.  Simple is
        shown first. Simple takes in one value and one rotation axis. The
        axis can be in 123 or XYZ.  

        >>> B.orient(N, 'Simple', q1, '3')
        >>> B.orient(N, 'Simple', q1, 'X')
        >>> B.orient(N, 'Simple', 1, 'Z')

        Next is Body. Body orientation takes this reference frame through
        three successive simple rotations. Acceptable rotation orders are
        of length 3, expressed in XYZ or 123, and cannot have a rotation
        about about an axis twice in a row.

        >>> B.orient(N, 'Body', [q1, q2, q3], '123')
        >>> B.orient(N, 'Body', [q1, q2, 0], 'ZXZ')
        >>> B.orient(N, 'Body', [0, 0, 0], 'XYX')

        Next is Space. Space is like Body, but the rotations are applied in the
        opposite order.

        >>> B.orient(N, 'Space', [q1, q2, q3], '312')

        Next is Euler. This orients the new ReferenceFrame with Euler
        Parameters, defined as a finite rotation about
        :math:`\vec{\hat \lambda}`, a unit vector, by some amount
        :math:`\theta`. This orientation is described by four parameters:
        q0 = :math:`cos(\frac{\theta}{2})`
        q1 = :math:`\lambda_x' sin(\frac{\theta}{2})`
        q2 = :math:`\lambda_y' sin(\frac{\theta}{2})`
        q3 = :math:`\lambda_z' sin(\frac{\theta}{2})`
        Euler does not take in a rotation order.

        >>> B.orient(N, 'Euler', [q1, q2, q3, q4])

        Last is Axis. This is a rotation about an arbitrary, non-time-varying
        axis by some angle. The axis is supplied as a Vector.

        >>> B.orient(N, 'Axis', [q1, N.x + 2 * N.y])

        """

        def _rot(axis, angle):
            """Returns direction cosine matrix for simple axis 1,2,or 3 rotations """
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

        approved_orders = ('123', '231', '312', '132', '213', '321', '121',
                           '131', '212', '232', '313', '323', '1', '2',
                           '3', '')
        rot_order = str(rot_order).upper() # Now we need to make sure XYZ = 123
        rot_type  = rot_type.upper()
        rot_order = [i.replace('X', '1') for i in rot_order]
        rot_order = [i.replace('Y', '2') for i in rot_order]
        rot_order = [i.replace('Z', '3') for i in rot_order]
        rot_order = ''.join(rot_order)
        if not isinstance(parent, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied as a parent')
        if not rot_order in approved_orders:
            raise TypeError('The supplied order is not an approved type')

        if rot_type == 'AXIS':
            if not rot_order == '':
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            assert isinstance(axis, Vector), 'A Vector needs to be supplied'
            assert axis.dt(parent) == 0, 'Axis cannot be time-varying'
            axis = axis.express(parent).unit
            axis = axis.args[0][0]
            self.parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
        elif rot_type == 'EULER':
            if not rot_order == '':
                raise TypeError('Euler orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 4)):
                raise TypeError('Amounts are a list or tuple of length 4')
            q0 = amounts[0]
            q1 = amounts[1]
            q2 = amounts[2]
            q3 = amounts[3]
            self.parent_orient = (Matrix([[q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 **
                2, 2 * (q1 * q2 - q0 * q3), 2 * (q0 * q2 + q1 * q3)],
                [2 * (q1 * q2 + q0 * q3), q0 ** 2 - q1 ** 2 + q2 **2 - q3 ** 2,
                2 * (q2 * q3 - q0 * q1)], [2 * (q1 * q3 - q0 * q2), 2 * (q0 *
                q1 + q2 * q3), q0 ** 2 - q1 ** 2 - q2 ** 2 + q3 ** 2]]))
        elif rot_type == 'BODY':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Body orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = (_rot(a1, amounts[0]) * _rot(a2, amounts[1])
                    * _rot(a3, amounts[2]))
        elif rot_type == 'SPACE':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Space orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = (_rot(a3, amounts[2]) * _rot(a2, amounts[1])
                    * _rot(a1, amounts[0]))
        elif rot_type == 'SIMPLE':
            if ((isinstance(amounts, (list, tuple))) |
                    (isinstance(rot_order, (list, tuple)))):
                raise TypeError('Simple takes 1 value for amount and order')
            a = int(rot_order)
            self.parent_orient = _rot(a, amounts)
        else:
            raise NotImplementedError('That is not an implemented rotation')
        self.parent = parent

    def set_ang_vel(self, value, otherframe):
        """Define the angular velocity vector of the ReferenceFrame.

        Parameters
        ==========
        value : Vector
            The angular velocity vector for this frame.
        otherframe : ReferenceFrame
            The frame which the angular velocity is defined in.

        Examples
        ========

        >>> from sympy.physics.classical.essential import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(V, N)

        """

        if value != 0:
            if not isinstance(value, Vector):
                raise TypeError('Angular velocity needs to be Vector.')
        if not isinstance(otherframe, ReferenceFrame):
            raise TypeError('Need to define the angular velocity with'
                            'respect to another ReferenceFrame.')
        self._ang_vel = value
        self._ang_vel_parent = otherframe

    @property
    def x(self):
        """The basis Vector for the ReferenceFrame, in the x direction. """
        return self._x

    @property
    def y(self):
        """The basis Vector for the ReferenceFrame, in the y direction. """
        return self._y

    @property
    def z(self):
        """The basis Vector for the ReferenceFrame, in the z direction. """
        return self._z


class Vector(object):
    """
    This is the class used to define vectors.  It along with reference
    frame are the building blocks of pydy.
    Class attributes include:
    subscript_indices - a 3 character string used for printing
    """

    subscript_indices = "xyz"

    def __init__(self, inlist):
        """
        This is the constructor for the Vector class.
        It should only be used in construction of the basis vectors,
        which is part of the ReferenceFrame construction.
        It takes in a SymPy matrix and a ReferenceFrame.
        """

        self.args = []
        while len(inlist) != 0:
            added = 0
            for i, v in enumerate(self.args):
                if inlist[0][1] == self.args[i][1]:
                    self.args[i] = (self.args[i][0] +
                            inlist[0][0], inlist[0][1])
                    inlist.remove(inlist[0])
                    added = 1
                    break
            if added != 1:
                self.args.append(inlist[0])
                inlist.remove(inlist[0])
        i = 0
        # This code is to remove empty frames from the list
        while i<len(self.args):
            if ((self.args[i][0][0] == 0) & (self.args[i][0][1] == 0) &
                (self.args[i][0][2] == 0)):
                self.args.remove(self.args[i])
                i -= 1
            i += 1

    def __str__(self):
        """
        Printing method.  Uses Vector Attribute subscript_indices to choose
        how to show basis vector indices.
        """
        ar = self.args # just to shorten things
        ol = [] # output list, to be concatenated to a string
        for i, v in enumerate(ar):
            for j in 0, 1, 2:
                # if the coef of the basis vector is 1, we skip the 1
                if ar[i][0][j] == 1:
                    if len(ol) != 0:
                        ol.append(' + ')
                    ol.append( ar[i][1].name.lower() +
                              self.subscript_indices[j] + '>' )
                # if the coef of the basis vector is -1, we skip the 1
                elif ar[i][0][j] == -1:
                    if len(ol) != 0:
                        ol.append(' ')
                    ol.append( '- ' + ar[i][1].name.lower() +
                              self.subscript_indices[j] + '>' )
                elif ar[i][0][j] != 0:
                    # If the coefficient of the basis vector is not 1 or -1,
                    # we wrap it in parentheses, for readability.
                    if len(ol) != 0:
                        ol.append(' + ')
                    ol.append('(' + `ar[i][0][j]` + ')*' +
                              ar[i][1].name.lower() +
                              self.subscript_indices[j] + '>' )
        return ''.join(ol)

    def __repr__(self):
        """
        Wraps __str__
        """
        return self.__str__()

    def __add__(self, other):
        """
        The add operator for Vector.
        It checks that other is a Vector, otherwise it throws an error.
        Also works with adding a zero scalar.
        """
        if isinstance(other, int):
            if other == 0:
                return self
        assert isinstance(other, Vector), 'You can only add two Vectors'
        return Vector(self.args + other.args)

    __radd__ = __add__

    def __and__(self, other):
        """
        Dot product of two vectors.
        """
        assert isinstance(other, Vector), 'Dot product is between two vectors'
        out = 0
        for i, v1 in enumerate(self.args):
            for j, v2 in enumerate(other.args):
                out += ((v2[0].T)
                        * (v2[1].dcm(v1[1]))
                        * (v1[0]))[0]
        out2 = simplify(out) # These lines are to simplify as much as we can
        while out != out2:
            out = out2
            out2 = simplify(out2)
        return out

    def __div__(self, other):
        """
        This uses mul and inputs self and 1 divided by other.
        """
        return self.__mul__(1 / other)

    def __eq__(self, other):
        """
        Tests for equality.  If other is 0, and self is empty, returns True.
        If other is 0 and self is not empty, returns False.
        If none of the above, only accepts other as a Vector.
        """
        if isinstance(other, int):
            if other == 0:
                if self.args == []:
                    return True
                else:
                    return False
        assert isinstance(other,Vector), 'Vectors can only compare to Vectors'
        dotcheck = (self & self) == (self & other)
        crosscheck = ((self ^ other) & (self ^ other) == 0)
        return dotcheck & crosscheck

    def __mul__(self, other):
        """
        Multiplies the Vector by a scalar.
        Throws an error if another Vector is entered.
        """
        assert not(isinstance(other, Vector)), \
                'Two Vectors can\'t be multiplied'
        newlist = [v for v in self.args]
        for i, v in enumerate(newlist):
            newlist[i] = (other * newlist[i][0], newlist[i][1])
        return Vector(newlist)

    def __neg__(self):
        return self * -1

    def __rmul__(self, other):
        """
        This wraps mul.
        """
        return self.__mul__(other)

    def __sub__(self, other):
        """
        The subraction operator.
        Reuses add and multiplication operations.
        """
        return self.__add__(other * -1)

    def __xor__(self, other):
        """
        The cross product operator for two Vectors.
        Takes in two Vectors; order matters.
        Returns a Vector which is perpendicular to the two input vectors.
        This Vector is expressed in the frames of self (first vector).
        """
        if isinstance(other, int):
            if other == 0:
                return self * 0
        assert isinstance(other, Vector), 'Cross products are between Vectors'
        def _det(mat):
            """
            This is needed as a little method for to find the determinant of a
            list in python; needs to work for a 3x3 list.
            SymPy's Matrix won't take in Vector, so need a custom function.
            """
            return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                    + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] *
                    mat[2][2]) + mat[0][2] * (mat[1][0] * mat[2][1] -
                    mat[1][1] * mat[2][0]))

        outvec = Vector([])
        ar = self.args # For brevity
        for i, v in enumerate(ar):
            tempx = v[1].x
            tempy = v[1].y
            tempz = v[1].z
            tempm = ([[tempx, tempy, tempz], [Vector([ar[i]]) & tempx,
                Vector([ar[i]]) & tempy, Vector([ar[i]]) & tempz],
                [other & tempx, other & tempy, other & tempz]])
            outvec += _det(tempm)
        return outvec

    def dot(self, other):
        """The dot product method.  See dot() function. """
        return self & other

    def cross(self, other):
        """The cross product method.  See cross() function. """
        return self ^ other

    def diff(self, wrt, otherframe):
        """
        Takes in a sympifyable value, which cannot be a Vector.
        Returns the partial derivative of the self Vector with respect
        to the input value.
        """
        wrt = sympify(wrt)
        assert isinstance(otherframe, ReferenceFrame), ('Need to supply a '
            'ReferenceFrame to find the derivative in')
        outvec = 0
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(wrt), otherframe)])
            else:
                diffed = (Vector(v).express(otherframe))[0].diff(wrt)
                outvec += Vector([(diff, otherframe)]).express(v[1])

    def dt(self, otherframe):
        """
        Returns the time derivative of the self Vector in the given Reference
        Frame.
        Takes in a ReferenceFrame, and returns a Vector.
        """
        assert isinstance(otherframe, ReferenceFrame), ('Need to supply a '
            'ReferenceFrame to find the derivative in')
        outvec = 0
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(Symbol('t')), otherframe)])
            else:
                outvec += (Vector([(v[0].diff(Symbol('t')), otherframe)]) +
                    Vector([(v[1].ang_vel_in(otherframe), otherframe)]) ^
                    Vector([v]))
        return outvec

    def express(self, otherframe):
        """
        Returns the vector, expressed in the other frame.
        Uses a DCM from each basis vector triplet's current frame to the other
        frame.
        Takes in a frame.
        """
        assert isinstance(otherframe, ReferenceFrame), ('Needs a frame to '
                                                        'express in')
        outvec = Vector(self.args + [])
        for i, v in enumerate(self.args):
            if v[1] != otherframe:
                outvec += Vector([(otherframe.dcm(v[1]) * v[0], otherframe)])
                outvec -= Vector([v])
        return outvec

    @property
    def mag(self):
        """Returns the magnitude of this Vector, as a scalar. """
        return sqrt(self & self)

    @property
    def unit(self):
        """Returns this Vector, with unit length. """
        return Vector(self.args + []) / self.mag


if __name__ == "__main__":
    import doctest
    doctest.testmod()

