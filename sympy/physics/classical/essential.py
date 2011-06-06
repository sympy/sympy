__all__ = ['ReferenceFrame', 'Vector']

from sympy import Matrix, Symbol, sin, cos, eye, simplify, diff

class ReferenceFrame(object):
    """
    ReferenceFrame is a reference frame.
    It will store its basis vectors as attributes, and orientation
    information to a parent frame, or it will be at the top of a tree.
    """

    def __init__(self, name=''):
        """
        Constructor for ReferenceFrame
        """
        self.name = name
        self.parent = None
        self._ang_vel = None
        self._ang_vel_parent = None
        self._cur = 0
        self._x = Vector([(Matrix([1, 0, 0]), self)])
        self._y = Vector([(Matrix([0, 1, 0]), self)])
        self._z = Vector([(Matrix([0, 0, 1]), self)])

    def __iter__(self):
        """
        Allows iteration through the basis vectors.
        """
        return self

    def __str__(self):
        """
        Returns the name of the frame.
        """
        return self.name

    def __repr__(self):
        """
        Wraps __str__
        """
        return self.name

    def next(self):
        """
        Allows iteration through the basis vectors.
        """
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
        """
        This returns the first common parent between two ReferenceFrames.
        Takes in another ReferenceFrame, and returns a ReferenceFrame.
        """
        assert isinstance(other, ReferenceFrame), 'You have to use a \
                ReferenceFrame'
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

    def ang_vel_in(self, other):
        """
        Returns the angular velocity of the current frame relative to
        the input frame: angular velocity of self in other
        Takes in ReferenceFrame, returns Vector.
        Form is A.ang_vel_in(N) is A's angular velocity in N, or N_w_A
        """
        commonframe = self._common_frame(other)
        # form DCM from self to first common frame
        leg1 = 0
        ptr = self
        while ptr != commonframe:
            leg1 += ptr._ang_vel
            ptr = ptr.parent
        # form DCM from other to first common frame
        leg2 = 0
        ptr = other
        while ptr != commonframe:
            leg2 -= ptr._ang_vel
            ptr = ptr.parent
        return leg1 + leg2

    def dcm(self, other):
        """
        This will return the direction cosine matrix from self to other;
        other's basis = DCM * self's basis is the definition.
        All it takes in is the other frame.
        nxyz = N.dcm(B) * bxyz
        """
        commonframe = self._common_frame(other)
        # form DCM from self to first common frame
        leg1 = eye(3)
        ptr = self
        while ptr != commonframe:
            leg1 = ptr.parent_orient * leg1
            ptr = ptr.parent
        # form DCM from other to first common frame
        leg2 = eye(3)
        ptr = other
        while ptr != commonframe:
            leg2 = ptr.parent_orient * leg2
            ptr = ptr.parent
        return leg2 * leg1.T

    def orientnew(self, newname, rot_type, amounts, rot_order=''):
        newframe = ReferenceFrame(newname)
        newframe.orient(self, rot_type, amounts, rot_order)
        return newframe

    def orient(self, parent, rot_type, amounts, rot_order=''):
        """
        This function will be used to define the orientation of a
        ReferenceFrame relative to a parent.  It takes in the parent frame,
        type of rotation, amount(s) of rotation, and if applicable, order of
        rotations.
        The format for this needs to be spelled out and made explicit, with
        example.
        Simple
        Body
        Space
        Euler - 4 euler params, no order
        Axis - angle + 3 amounts for unit vector direction, no order
        """

        def _rot(axis, angle):
            """
            Returns direction cosine matrix for simple axis 1,2,or 3 rotations
            """
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

        assert isinstance(parent, ReferenceFrame), 'Parent must be a \
                ReferenceFrame'
        approved_orders = ('123', '231', '312', '132', '213', '321', '121',
                '131', '212', '232', '313', '323', '1', '2', '3', '')
        rot_order = str(rot_order).upper() # Now we need to make sure XYZ = 123
        rot_type  = rot_type.upper()
        rot_order = [i.replace('X', '1') for i in rot_order]
        rot_order = [i.replace('Y', '2') for i in rot_order]
        rot_order = [i.replace('Z', '3') for i in rot_order]
        rot_order = ''.join(rot_order)
        assert rot_order in approved_orders, 'Not approved order'

        if rot_type == 'AXIS':
            assert rot_order == '', 'Axis orientation take no rotation order'
            assert isinstance(amounts, (list, tuple)) & (len(amounts) == 4),\
                    'Amounts need to be in a list or tuple of length 4'
            theta = amounts[0]
            axis = amounts[1:]
            assert ((diff(axis[0], Symbol('t')) == 0) &
                    (diff(axis[1], Symbol('t')) == 0) &
                    (diff(axis[2], Symbol('t')) == 0)), 'Axis directions \
                    cannot be time-varying'
            axis = Matrix(axis)
            mag = simplify(simplify((axis.T * axis)[0]))
            if mag != 1:
                axis /= mag
            self.parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
        elif rot_type == 'EULER':
            assert isinstance(amounts, (list, tuple)) & len(amounts) == 4, \
                    'Amounts need to be in a list or tuple of length 4'
            assert rot_order == '', 'Euler orientation take no rotation order'
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
            assert len(amounts) == 3, 'Body orientation requires 3 values'
            assert len(rot_order) == 3, 'Body orientation requires 3 orders'
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = (_rot(a1, amounts[0]) * _rot(a2, amounts[1])
                    * _rot(a3, amounts[2]))
        elif rot_type == 'SPACE':
            assert len(amounts) == 3, 'Space orientation requires 3 values'
            assert len(rot_order) == 3, 'Space orientation requires 3 orders'
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            self.parent_orient = (_rot(a3, amounts[2]) * _rot(a2, amounts[1])
                    * _rot(a1, amounts[0]))
        elif rot_type == 'SIMPLE':
            assert not(isinstance(amounts, (list, tuple))), 'Simple takes in \
                    a single value for amount'
            assert not(isinstance(rot_order, (list, tuple))), 'Simple takes \
                    in a single value for rotation order'
            a = int(rot_order)
            self.parent_orient = _rot(a, amounts)
        else:
            raise NotImplementedError('That is not an implemented rotation')
        self.parent = parent

    def set_ang_vel(self, value, other):
        """
        Define the angular velocity vector of the current ReferenceFrame,
        with respect to the second ReferenceFrame.
        Takes in a Vector for angular velocity vector, and ReferenceFrame,
        for the frame to this to be defined in.
        """
        if value != 0:
            assert isinstance(value, Vector), 'Angular velocity needs to be a \
        Vector.'
        assert isinstance(other, ReferenceFrame), 'Need to define the \
        angular velocity with respect to another ReferenceFrame.'
        self._ang_vel = value
        self._ang_vel_parent = other

    @property
    def x(self):
        """
        The basis vector for the ReferenceFrame, in the x (or 1, or i)
        direction.  Immutable.
        Returns a Vector.
        """
        return self._x

    @property
    def y(self):
        """
        The basis vector for the ReferenceFrame, in the y (or 2, or j)
        direction.  Immutable.
        Returns a Vector.
        """
        return self._y

    @property
    def z(self):
        """
        The basis vector for the ReferenceFrame, in the z (or 3, or k)
        direction.  Immutable.
        Returns a Vector.
        """
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
        """
        Wraps around the & operator, which is the dot product operator.
        """
        return self & other

    def cross(self, other):
        """
        Wraps around the ^ operator, which is the cross product operator.
        """
        return self ^ other

    def diff(self, wrt, otherframe):
        """
        Takes in a sympifyable value, which cannot be a Vector.
        Returns the partial derivative of the self Vector with respect
        to the input value.
        """
        wrt = sympify(wrt)
        assert isinstance(otherframe, ReferenceFrame), 'Need to supply a \
                ReferenceFrame to find the derivative in'
        outvec = 0
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector(v[0].diff(wrt), otherframe)
            else:
                diffed = (Vector(v).express(otherframe))[0].diff(wrt)
                outvec += Vector(diff, otherframe).express(v[1])

    def dt(self, otherframe):
        """
        Returns the time derivative of the self Vector in the given Reference
        Frame.
        Takes in a ReferenceFrame, and returns a Vector.
        """
        assert isinstance(otherframe, ReferenceFrame), 'Need to supply a \
                ReferenceFrame to find the derivative in'
        outvec = 0
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector(v[0].diff(Symbol('t')), otherframe)
            else:
                outvec += (Vector(v[0].diff(Symbol('t')), otherframe) +
                    Vector(v[1].ang_vel_in(otherframe), otherframe) ^
                    Vector([v]))
        return outvec

    def express(self, otherframe):
        """
        Returns the vector, expressed in the other frame.
        Uses a DCM from each basis vector triplet's current frame to the other
        frame.
        Takes in a frame.
        """
        assert isinstance(otherframe, ReferenceFrame), 'Needs a frame to \
                express in'
        outvec = Vector(self.args + [])
        for i, v in enumerate(self.args):
            if v[1] != otherframe:
                outvec += Vector([(otherframe.dcm(v[1]) * v[0], otherframe)])
                outvec -= Vector([v])
        return outvec

    @property
    def mag(self):
        """
        Returns the magnitude of the Vector.
        """
        return sqrt(self & self)

    @property
    def unit(self):
        """
        Returns a vector of length one in the direction of the Vector.
        """
        return Vector(self.args + []) / self.mag
