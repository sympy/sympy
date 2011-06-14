__all__ = ['ReferenceFrame', 'Vector', 'Dyad']

from sympy import (Matrix, Symbol, sin, cos, eye, trigsimp, diff, sqrt, sympify,
                   expand, S)
from sympy.core.numbers import Zero

class Dyad(object):
    """A Dyad object.

    See:
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    A more powerful way to represent a rigid body's inertia. While it is more
    complex, by choosing Dyad components to be in body fixed basis vectors,
    the resulting matrix is equivalent to the inertia tensor.

    """

    def __init__(self, inlist):
        """
        Just like Vector's init, you shouldn't call this.
        Stores a Dyad as a list of lists; the inner list has the measure
        number and the two unit vectors; the outerlist holds each unique unit
        vector pair.

        """

        self.args = []
        while len(inlist) != 0:
            added = 0
            for i, v in enumerate(self.args):
                if ((inlist[0][1].__str__() == self.args[i][1].__str__()) and
                    (inlist[0][2].__str__() == self.args[i][2].__str__())):
                    self.args[i] = (self.args[i][0] +
                        inlist[0][0], inlist[0][1], inlist[0][2])
                    inlist.remove(inlist[0])
                    added = 1
                    break
            if added != 1:
                self.args.append(inlist[0])
                inlist.remove(inlist[0])
        i = 0
        # This code is to remove empty parts from the list
        while i < len(self.args):
            if ((self.args[i][0] == 0) | (self.args[i][1] == 0) | 
                (self.args[i][2] == 0)):
                self.args.remove(self.args[i])
                i -= 1
            i += 1

    def __str__(self):
        """Printing method. """
        ar = self.args
        ol = []
        for i, v in enumerate(ar):
            if ar[i][0] == 1:
                if len(ol) != 0:
                    ol.append(' + ')
                ol.append(ar[i][1].__str__() + ar[i][2].__str__())
            elif ar[i][0] == -1:
                if len(ol) != 0:
                    ol.append(' ')
                ol.append('- ' + ar[i][1].__str__() + ar[i][2].__str__())
            elif ar[i][0] != 0:
                if len(ol) != 0:
                    ol.append(' + ')
                ol.append('(' + `ar[i][0]` + ')*' + ar[i][1].__str__() +
                        ar[i][2].__str__())
        return ''.join(ol)

    def __add__(self, other):
        """The add operator for Dyad. """
        if isinstance(other, (int, type(Zero()))):
            if other == 0:
                return self
        self._check_dyad(other)
        return Dyad(self.args + other.args)

    def __and__(self, other):
        """The inner product operator for a Dyad and a Dyad or Vector.

        Parameters
        ==========
        other : Dyad or Vector
            The other Dyad or Vector to take the inner product with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> D1 = outer(N.x, N.y)
        >>> D2 = outer(N.y, N.y)
        >>> D1.dot(D2)
        nx>ny>
        >>> D1.dot(N.y)
        nx>

        """

        if isinstance(other, int):
            if other == 0:
                return 0
        ol = 0
        if isinstance(other, Dyad):
            for i, v in enumerate(self.args):
                for i2, v2 in enumerate(other.args):
                    ol += v[0] * v2[0] * (v[2] & v2[1]) * (v[1] | v2[2])
        elif isinstance(other, Vector):
            for i, v in enumerate(self.args):
                ol += v[0] * v[1] * (v[2] & other)
        else:
            raise TypeError('Need to supply a Vector or Dyad')
        return ol

    def __div__(self, other):
        """Divides the Dyad by a sympifyable expression. """
        return self.__mul__(1 / other)

    def __eq__(self, other):
        """Tests for equality. 

        Is currently weak; needs stronger comparison testing

        """

        if isinstance(other, (int, type(Zero()))):
            if (other == 0) & (self.args == []):
                return True
        self._check_dyad(other)
        return set(self.args) == set(other.args)

    def __mul__(self, other):
        """Multiplies the Dyad by a sympifyable expression.

        Parameters
        ==========
        other : Sympafiable
            The scalar to multiply this Dyad with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> 5 * d
        (5)*nx>nx>

        """

        newlist = [v for v in self.args]
        for i, v in enumerate(newlist):
            newlist[i] = (sympify(other) * newlist[i][0], newlist[i][1],
                    newlist[i][2])
        return Dyad(newlist)

    def __rand__(self, other):
        """The inner product operator for a Vector or Dyad, and a Dyad

        This is for: Vector dot Dyad

        Parameters
        ==========
        other : Vector
            The vector we are dotting with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, dot, outer
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> dot(N.x, d)
        nx>

        """

        if isinstance(other, int):
            if other == 0:
                return 0
        ol = 0
        if isinstance(other, Vector):
            for i, v in enumerate(self.args):
                ol += v[0] * v[2] * (v[1] & other)
        else:
            raise TypeError('Need to supply a Vector or Dyad')
        return ol

    def __rsub__(self, other):
        return (-1 * self) + other

    def __rxor__(self, other):
        """For a cross product in the form: Vector x Dyad

        Parameters
        ==========
        other : Vector
            The Vector that we are crossing this Dyad with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(N.y, d)
        - nz>nx>

        """

        self._check_vector(other)
        ol = 0
        for i, v in enumerate(self.args):
            ol += v[0] * ((other ^ v[1]) | v[2])
        return ol

    def __sub__(self, other):
        """The subtraction operator. """
        return self.__add__(other * -1)

    def __xor__(self, other):
        """For a cross product in the form: Dyad x Vector.

        Parameters
        ==========
        other : Vector
            The Vector that we are crossing this Dyad with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(d, N.y)
        nx>nz>

        """

        self._check_vector(other)
        ol = 0
        for i, v in enumerate(self.args):
            ol += v[0] * (v[1] | (v[2] ^ other))
        return ol

    __repr__ = __str__
    __radd__ = __add__
    __rmul__ = __mul__
    dot = __and__
    cross = __and__

    def _check_frame(self, other):
        if not isinstance(other, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied')

    def _check_dyad(self, other):
        if isinstance(other, (int, type(Zero()))):
            if other == 0:
                return
        if not isinstance(other, Dyad):
            raise TypeError('A Dyad must be supplied')

    def _check_vector(self, other):
        if isinstance(other, (int, type(Zero()))):
            if other == 0:
                return
        if not isinstance(other, Vector):
            raise TypeError('A Vector must be supplied')

    def express(self, frame1, frame2=None):
        """Expresses this Dyad in alternate frame(s)

        The first frame is the list side expression, the second frame is the
        right side; if Dyad is in form ax>by>, you can express it in two
        different frames. If no second frame is given, the Dyad is
        expressed in only one frame.

        Parameters
        ==========
        frame1 : ReferenceFrame
            The frame to express the left side of the Dyad in
        frame2 : ReferenceFrame
            If provided, the frame to express the right side of the Dyad in

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer,\
                DynamicSymbol
        >>> N = ReferenceFrame('N')
        >>> q = DynamicSymbol('q')
        >>> B = N.orientnew('B', 'Simple', q, 3)
        >>> d = outer(N.x, N.x)
        >>> d.express(B, N)
        (cos(q))*bx>nx> + (-sin(q))*by>nx>

        """

        if frame2 == None:
            frame2 = frame1
        self._check_frame(frame1)
        self._check_frame(frame2)
        ol = 0
        for i, v in enumerate(self.args):
            ol += v[0] * (v[1].express(frame1) | v[2].express(frame2))
        return ol

    def dt(self, frame):
        """Take the time derivative of this Dyad in a frame.

        Parameters
        ==========
        frame : ReferenceFrame
            The frame to take the time derivative in

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer,\
                DynamicSymbol
        >>> N = ReferenceFrame('N')
        >>> q = DynamicSymbol('q')
        >>> B = N.orientnew('B', 'Simple', q, 3)
        >>> d = outer(N.x, N.x)
        >>> d.dt(B)
        (-qd)*ny>nx> + (-qd)*nx>ny>
 
        """

        self._check_frame(frame)
        t = Symbol('t')
        ol = 0
        for i, v in enumerate(self.args):
            ol += (v[0].diff(t) * (v[1] | v[2]))
            ol += (v[0] * (v[1].dt(frame) | v[2]))
            ol += (v[0] * (v[1] | v[2].dt(frame)))
        return ol


class ReferenceFrame(object):
    """A reference frame in classical mechanics.

    ReferenceFrame is a class used to represent a reference frame in classical
    mechanics. It has a standard basis of three unit vectors in the frame's
    x, y, and z directions.
    It also can have a rotation relative to a parent
    frame; this rotation is defined by a direction cosine matrix relating this
    frame's basis vectors to the parent frame's basis vectors.
    It can also have an angular velocity vector, defined in another frame.

    """

    def __init__(self, name=''):
        """init for ReferenceFrame. """
        if not isinstance(name, str):
            raise TypeError('Need to supply a valid name')
        self.name = name
        self._dcm_dict = {}
        self._ang_vel_dict = {}
        self._ang_acc_dict = {}
        self._dlist = [self._dcm_dict, self._ang_vel_dict, self._ang_acc_dict]
        self._cur = 0
        self._x = Vector([(Matrix([1, 0, 0]), self)])
        self._y = Vector([(Matrix([0, 1, 0]), self)])
        self._z = Vector([(Matrix([0, 0, 1]), self)])

    def __iter__(self):
        return self

    def __str__(self):
        """Returns the name of the frame. """
        return self.name

    __repr__ = __str__

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

    def _check_frame(self, other):
        if not isinstance(other, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied')

    def _check_vector(self, other):
        if isinstance(other, (int, type(Zero()))):
            if other == 0:
                return
        if not isinstance(other, Vector):
            raise TypeError('A Vector must be supplied')

    def _w_diff_dcm(self, otherframe):
        """Angular velocity from time differentiating the DCM. """
        dcm2diff = self.dcm(otherframe)
        diffed = dcm2diff.diff(Symbol('t'))
        angvelmat = diffed * dcm2diff.T
        w1 = trigsimp(expand(angvelmat[7]), recursive=True)
        w2 = trigsimp(expand(angvelmat[2]), recursive=True)
        w3 = trigsimp(expand(angvelmat[3]), recursive=True)
        return -Vector([(Matrix([w1, w2, w3]), self)])

    def _dict_list(self, other, num):
        """Creates a list from self to other using _dcm_dict. """
        outlist = [[self]]
        oldlist = [[]]
        while outlist != oldlist:
            oldlist = outlist[:]
            for i, v in enumerate(outlist):
                templist = v[-1]._dlist[num].keys()
                for i2, v2 in enumerate(templist):
                    if not v.__contains__(v2):
                        littletemplist = v + [v2]
                        if not outlist.__contains__(littletemplist):
                            outlist.append(littletemplist)
        for i, v in enumerate(oldlist):
            if v[-1] != other:
                outlist.remove(v)
        outlist.sort(key = len)
        if len(outlist) != 0:
            return outlist[0]
        raise ValueError('No Connecting Path found between ' + self.name +
                         ' and ' + other.name)

    def ang_acc_in(self, otherframe):
        """Returns the angular acceleration Vector of the ReferenceFrame.

        Effectively returns the Vector:
        :math:`^{N} \vec{\alpha} ^{B}`
        which represent the angular acceleration of B in N, where B is self, and
        N is otherframe.

        Parameters
        ==========
        otherframe : ReferenceFrame
            The ReferenceFrame which the angular acceleration is returned in.

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_acc(N, V)
        >>> A.ang_acc_in(N)
        (10)*nx>

        """

        self._check_frame(otherframe)
        if self._ang_acc_dict.has_key(otherframe):
            return self._ang_acc_dict[otherframe]
        else:
            return self.ang_vel_in(otherframe).dt(otherframe)

    def ang_vel_in(self, otherframe):
        """Returns the angular velocity Vector of the ReferenceFrame.

        Effectively returns the Vector:
        :math:`^{N} \vec{\omega} ^{B}`
        which represent the angular velocity of B in N, where B is self, and
        N is otherframe.

        Parameters
        ==========
        otherframe : ReferenceFrame
            The ReferenceFrame which the angular velocity is returned in.

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        (10)*nx>

        """

        self._check_frame(otherframe)
        flist = self._dict_list(otherframe, 1)
        outvec = 0
        for i in range(len(flist) - 1):
            outvec += flist[i]._ang_vel_dict[flist[i + 1]]
        return outvec

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

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 1)
        >>> N.dcm(A)
        [1,       0,        0]
        [0, cos(q1), -sin(q1)]
        [0, sin(q1),  cos(q1)]

        """

        self._check_frame(otherframe)
        flist = self._dict_list(otherframe, 0)
        outdcm = eye(3)
        for i in range(len(flist) - 1):
            outdcm = outdcm * flist[i + 1]._dcm_dict[flist[i]]
        return outdcm

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

        >>> from sympy.physics.classical import ReferenceFrame, Vector
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

        Supported orientation types are Simple, Body, Space, Quaternion, Axis.
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

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q0, q1, q2, q3, q4 = symbols('q0 q1 q2 q3 q4')
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

        Next is Quaternion. This orients the new ReferenceFrame with
        Quaternions, defined as a finite rotation about
        :math:`\vec{\hat \lambda}`, a unit vector, by some amount
        :math:`\theta`. This orientation is described by four parameters:
        q0 = :math:`cos(\frac{\theta}{2})`
        q1 = :math:`\lambda_x' sin(\frac{\theta}{2})`
        q2 = :math:`\lambda_y' sin(\frac{\theta}{2})`
        q3 = :math:`\lambda_z' sin(\frac{\theta}{2})`
        Quaternion does not take in a rotation order.

        >>> B.orient(N, 'Quaternion', [q0, q1, q2, q3])

        Last is Axis. This is a rotation about an arbitrary, non-time-varying
        axis by some angle. The axis is supplied as a Vector.

        >>> B.orient(N, 'Axis', [q1, N.x + 2 * N.y])

        """

        self._check_frame(parent)
        def _rot(axis, angle):
            """DCM for simple axis 1,2,or 3 rotations. """
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
        if not rot_order in approved_orders:
            raise TypeError('The supplied order is not an approved type')
        parent_orient = []

        if rot_type == 'AXIS':
            if not rot_order == '':
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            self._check_vector(axis)
            if not axis.dt(parent) == 0:
                raise ValueError('Axis cannot be time-varying')
            axis = axis.express(parent).unit
            axis = axis.args[0][0]
            parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
        elif rot_type == 'QUATERNION':
            if not rot_order == '':
                raise TypeError('Quaternion orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 4)):
                raise TypeError('Amounts are a list or tuple of length 4')
            q0 = amounts[0]
            q1 = amounts[1]
            q2 = amounts[2]
            q3 = amounts[3]
            parent_orient = (Matrix([[q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 **
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
            parent_orient = (_rot(a1, amounts[0]) * _rot(a2, amounts[1])
                    * _rot(a3, amounts[2]))
        elif rot_type == 'SPACE':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Space orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            parent_orient = (_rot(a3, amounts[2]) * _rot(a2, amounts[1])
                    * _rot(a1, amounts[0]))
        elif rot_type == 'SIMPLE':
            if ((isinstance(amounts, (list, tuple))) |
                    (isinstance(rot_order, (list, tuple)))):
                raise TypeError('Simple takes 1 value for amount and order')
            a = int(rot_order)
            parent_orient = _rot(a, amounts)
        else:
            raise NotImplementedError('That is not an implemented rotation')
        self._dcm_dict.update({parent: parent_orient})
        parent._dcm_dict.update({self: parent_orient.T})
        # TODO double check the sign here
        if rot_type == 'QUATERNION':
            t = Symbol('t')
            q0 = amounts[0]
            q1 = amounts[1]
            q2 = amounts[2]
            q3 = amounts[3]
            q0d = diff(q0, t)
            q1d = diff(q1, t)
            q2d = diff(q2, t)
            q3d = diff(q3, t)
            w1 = 2 * (q1d * q0 + q2d * q3 - q3d * q2 - q0d * q1)
            w2 = 2 * (q2d * q0 + q3d * q1 - q1d * q3 - q0d * q2)
            w3 = 2 * (q3d * q0 + q1d * q2 - q2d * q1 - q0d * q3)
            wvec = Vector([(Matrix([w1, w2, w3]), self)])
        elif rot_type == 'AXIS':
            thetad = (amounts[0]).diff(Symbol('t'))
            wvec = thetad * amounts[1].express(parent).unit
        else:
            wvec = self._w_diff_dcm(parent)
        self._ang_vel_dict.update({parent: wvec})
        parent._ang_vel_dict.update({self: -wvec})

    def set_ang_acc(self, otherframe, value):
        """Define the angular acceleration Vector in a ReferenceFrame.

        Defines the angular acceleration of this ReferenceFrame, in another.
        Angular acceleration can be defined with respect to multiple different
        ReferenceFrames. Care must be taken to not create loops which are
        inconsistent.

        Parameters
        ==========
        otherframe : ReferenceFrame
            A ReferenceFrame to define the angular acceleration in
        value : Vector
            The Vector representing angular acceleration

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_acc(N, V)
        >>> A.ang_acc_in(N)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_frame(otherframe)
        self._ang_acc_dict.update({otherframe: value})
        otherframe._ang_acc_dict.update({self: -value})

    def set_ang_vel(self, otherframe, value):
        """Define the angular velocity vector in a ReferenceFrame.

        Defines the angular velocity of this ReferenceFrame, in another.
        Angular velocity can be defined with respect to multiple different
        ReferenceFrames. Care must be taken to not create loops which are
        inconsistent.

        Parameters
        ==========
        otherframe : ReferenceFrame
            A ReferenceFrame to define the angular velocity in
        value : Vector
            The Vector representing angular velocity

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_frame(otherframe)
        self._ang_vel_dict.update({otherframe: value})
        otherframe._ang_vel_dict.update({self: -value})

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
    """The class used to define vectors.

    It along with ReferenceFrame are the building blocks of describing a
    classical mechanics system in PyDy.

    Attributes
    ==========
    subscript_indices : str
        A 3 character string used for printing the basis vectors
        This needs to be changed as Vector.subscript_indices = "123", and not
        as SomeVectorInstance.subscript_indices = "xyz"

    """

    subscript_indices = "xyz"

    def __init__(self, inlist):
        """This is the constructor for the Vector class.
        You shouldn't be calling this, it should only be used by other
        functions. You should be treating Vectors like you would with if you
        were doing the math by hand, and getting the first 3 from the
        standard basis vectors from a ReferenceFrame.

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
        while i < len(self.args):
            if self.args[i][0] == Matrix([0, 0, 0]):
                self.args.remove(self.args[i])
                i -= 1
            i += 1

    def __str__(self):
        """Printing method. """
        ar = self.args # just to shorten things
        ol = [] # output list, to be concatenated to a string
        for i, v in enumerate(ar):
            for j in 0, 1, 2:
                # if the coef of the basis vector is 1, we skip the 1
                if ar[i][0][j] == 1:
                    if len(ol) != 0:
                        ol.append(' + ')
                    ol.append( ar[i][1].name.lower() +
                              Vector.subscript_indices[j] + '>' )
                # if the coef of the basis vector is -1, we skip the 1
                elif ar[i][0][j] == -1:
                    if len(ol) != 0:
                        ol.append(' ')
                    ol.append('- ' + ar[i][1].name.lower() +
                              Vector.subscript_indices[j] + '>' )
                elif ar[i][0][j] != 0:
                    # If the coefficient of the basis vector is not 1 or -1,
                    # we wrap it in parentheses, for readability.
                    if len(ol) != 0:
                        ol.append(' + ')
                    ol.append('(' + `ar[i][0][j]` + ')*' +
                              ar[i][1].name.lower() +
                              Vector.subscript_indices[j] + '>' )
        return ''.join(ol)

    def __add__(self, other):
        """The add operator for Vector. """
        if isinstance(other, int):
            if other == 0:
                return self
        self._check_vector(other)
        return Vector(self.args + other.args)

    def __and__(self, other):
        """Dot product of two vectors.

        Returns a scalar, the dot product of the two Vectors

        Parameters
        ==========
        other : Vector
            The Vector which we are dotting with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector, dot
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> dot(N.x, N.x)
        1
        >>> dot(N.x, N.y)
        0
        >>> A = N.orientnew('A', 'Simple', q1, 1)
        >>> dot(N.y, A.y)
        cos(q1)

        """

        if isinstance(other, Dyad):
            return NotImplemented
        self._check_vector(other)
        out = 0
        for i, v1 in enumerate(self.args):
            for j, v2 in enumerate(other.args):
                out += ((v2[0].T)
                        * (v2[1].dcm(v1[1]))
                        * (v1[0]))[0]
        return trigsimp(sympify(out), recursive=True)

    def __div__(self, other):
        """This uses mul and inputs self and 1 divided by other. """
        return self.__mul__(1 / other)

    def __eq__(self, other):
        """Tests for equality.

        It is very import to note that this is only as good as the SymPy
        equality test; False does not always mean they are not equivalent
        Vectors.
        If other is 0, and self is empty, returns True.
        If other is 0 and self is not empty, returns False.
        If none of the above, only accepts other as a Vector.

        """

        if isinstance(other, int):
            if other == 0:
                if self.args == []:
                    return True
                else:
                    return False
        check = True
        frame = self.args[0][1]
        for i, v in enumerate(frame):
            if expand((self - other) & v) != 0:
                return False
        return True

    def __mul__(self, other):
        """Multiplies the Vector by a sympifyable expression.

        Parameters
        ==========
        other : Sympifyable
            The scalar to multiply this Vector with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> from sympy import Symbol
        >>> N = ReferenceFrame('N')
        >>> t = Symbol('t')
        >>> V = 10 * t * N.x
        >>> print V
        (10*t)*nx>

        """

        newlist = [v for v in self.args]
        for i, v in enumerate(newlist):
            newlist[i] = (sympify(other) * newlist[i][0], newlist[i][1])
        return Vector(newlist)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __neg__(self):
        return self * -1

    def __or__(self, other):
        """Outer product between two Vectors.

        A rank increasing operation, which returns a Dyad from two Vectors

        Parameters
        ==========
        other : Vector
            The Vector to take the outer product with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> outer(N.x, N.x)
        nx>nx>
        
        """

        if other == 0:
            return 0
        ol = 0
        for i, v in enumerate(self.args):
            for i2, v2 in enumerate(other.args):
                # it looks this way because if we are in the same frame and
                # use the enumerate function on the same frame in a nested
                # fashion, then bad things happen
                ol += Dyad([(v[0][0] * v2[0][0], v[1].x, v2[1].x)])
                ol += Dyad([(v[0][0] * v2[0][1], v[1].x, v2[1].y)])
                ol += Dyad([(v[0][0] * v2[0][2], v[1].x, v2[1].z)])
                ol += Dyad([(v[0][1] * v2[0][0], v[1].y, v2[1].x)])
                ol += Dyad([(v[0][1] * v2[0][1], v[1].y, v2[1].y)])
                ol += Dyad([(v[0][1] * v2[0][2], v[1].y, v2[1].z)])
                ol += Dyad([(v[0][2] * v2[0][0], v[1].z, v2[1].x)])
                ol += Dyad([(v[0][2] * v2[0][1], v[1].z, v2[1].y)])
                ol += Dyad([(v[0][2] * v2[0][2], v[1].z, v2[1].z)])
        return ol

    def __ror__(self, other):
        """Outer product between two Vectors.

        A rank increasing operation, which returns a Dyad from two Vectors

        Parameters
        ==========
        other : Vector
            The Vector to take the outer product with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> outer(N.x, N.x)
        nx>nx>
        
        """

        if other == 0:
            return 0
        ol = 0
        for i, v in enumerate(other.args):
            for i2, v2 in enumerate(self.args):
                # it looks this way because if we are in the same frame and
                # use the enumerate function on the same frame in a nested
                # fashion, then bad things happen
                ol += Dyad([(v[0][0] * v2[0][0], v[1].x, v2[1].x)])
                ol += Dyad([(v[0][0] * v2[0][1], v[1].x, v2[1].y)])
                ol += Dyad([(v[0][0] * v2[0][2], v[1].x, v2[1].z)])
                ol += Dyad([(v[0][1] * v2[0][0], v[1].y, v2[1].x)])
                ol += Dyad([(v[0][1] * v2[0][1], v[1].y, v2[1].y)])
                ol += Dyad([(v[0][1] * v2[0][2], v[1].y, v2[1].z)])
                ol += Dyad([(v[0][2] * v2[0][0], v[1].z, v2[1].x)])
                ol += Dyad([(v[0][2] * v2[0][1], v[1].z, v2[1].y)])
                ol += Dyad([(v[0][2] * v2[0][2], v[1].z, v2[1].z)])
        return ol

    def __rsub__(self, other):
        return (-1 * self) + other

    def __sub__(self, other):
        """The subraction operator. """
        return self.__add__(other * -1)

    def __xor__(self, other):
        """The cross product operator for two Vectors.

        Returns a Vector, expressed in the same ReferenceFrames as self.

        Parameters
        ==========
        other : Vector
            The Vector which we are crossing with

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> N.x ^ N.y
        nz>
        >>> A = N.orientnew('A', 'Simple', q1, 1)
        >>> A.x ^ N.y
        (sin(q1))*ay> + (cos(q1))*az>
        >>> N.y ^ A.x
        - nz>

        """

        if isinstance(other, Dyad):
            return NotImplemented
        if isinstance(other, int):
            if other == 0:
                return self * 0
        self._check_vector(other)

        def _det(mat):
            """This is needed as a little method for to find the determinant
            of a list in python; needs to work for a 3x3 list.
            SymPy's Matrix won't take in Vector, so need a custom function.
            You shouldn't be calling this.

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

    def _check_frame(self, other):
        if not isinstance(other, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied')

    def _check_vector(self, other):
        if isinstance(other, (int, type(Zero()))):
            if other == 0:
                return
        if not isinstance(other, Vector):
            raise TypeError('A Vector must be supplied')

    __repr__ = __str__
    __radd__ = __add__
    __rmul__ = __mul__

    def dot(self, other):
        return self & other
    dot.__doc__ = __and__.__doc__

    def cross(self, other):
        return self ^ other
    cross.__doc__ = __xor__.__doc__

    def outer(self, other):
        return self | other
    outer.__doc__ = __or__.__doc__

    def diff(self, wrt, otherframe):
        """Takes the partial derivative, with respect to a value, in a frame.

        Returns a Vector.

        Parameters
        ==========
        wrt : Symbol
            What the partial derivative is taken with respect to.
        otherframe : ReferenceFrame
            The ReferenceFrame that the partial derivative is taken in.

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector,\
                DynamicSymbol
        >>> from sympy import Symbol
        >>> t = Symbol('t')
        >>> q1 = DynamicSymbol('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 2)
        >>> A.x.diff(t, N)
        (-q1d)*az>

        """

        wrt = sympify(wrt)
        self._check_frame(otherframe)
        outvec = 0
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(wrt), otherframe)])
            else:
                diffed = (Vector([v]).express(otherframe)).args[0][0].diff(wrt)
                outvec += Vector([(diffed, otherframe)]).express(v[1])
        return outvec

    def dt(self, otherframe):
        """Returns the time derivative of the Vector in a ReferenceFrame.

        Returns a Vector which is the time derivative of the self Vector, taken
        in frame otherframe.

        Parameters
        ==========
        otherframe : ReferenceFrame
            The ReferenceFrame that the partial derivative is taken in.

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector,\
                DynamicSymbol
        >>> from sympy import Symbol
        >>> q1 = Symbol('q1')
        >>> u1 = DynamicSymbol('u1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 1)
        >>> v = u1 * N.x
        >>> A.set_ang_vel(N, 10*A.x)
        >>> A.x.dt(N) == 0
        True
        >>> v.dt(N)
        (u1d)*nx>

        """

        outvec = 0
        self._check_frame(otherframe)
        for i,v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(Symbol('t')), otherframe)])
            else:
                outvec += (Vector([v]).dt(v[1]) -
                    (Vector([v]) ^ v[1].ang_vel_in(otherframe)))
        return outvec

    def express(self, otherframe):
        """Returns a vector, expressed in the other frame.

        A new Vector is returned, equalivalent to this Vector, but its
        components are all defined in only the otherframe.

        Parameters
        ==========
        otherframe : ReferenceFrame
            The frame for this Vector to be described in

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Vector,\
                DynamicSymbol
        >>> q1 = DynamicSymbol('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Simple', q1, 2)
        >>> A.x.express(N)
        (cos(q1))*nx> + (-sin(q1))*nz>

        """

        self._check_frame(otherframe)
        outvec = Vector(self.args + [])
        for i, v in enumerate(self.args):
            if v[1] != otherframe:
                temp = otherframe.dcm(v[1]) * v[0]
                temp.simplify()
                for i2, v2 in enumerate(temp):
                    temp[i2] = trigsimp(v2, recursive=True)
                outvec += Vector([(temp, otherframe)])
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

