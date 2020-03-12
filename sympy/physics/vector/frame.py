from sympy.core.backend import (diff, expand, sin, cos, sympify,
                   eye, symbols, ImmutableMatrix as Matrix, MatrixBase)
from sympy import (trigsimp, solve, Symbol, Dummy)
from sympy.physics.vector.vector import Vector, _check_vector
from sympy.utilities.misc import translate

__all__ = ['CoordinateSym', 'ReferenceFrame']


class CoordinateSym(Symbol):
    """
    A coordinate symbol/base scalar associated wrt a Reference Frame.

    Ideally, users should not instantiate this class. Instances of
    this class must only be accessed through the corresponding frame
    as 'frame[index]'.

    CoordinateSyms having the same frame and index parameters are equal
    (even though they may be instantiated separately).

    Parameters
    ==========

    name : string
        The display name of the CoordinateSym

    frame : ReferenceFrame
        The reference frame this base scalar belongs to

    index : 0, 1 or 2
        The index of the dimension denoted by this coordinate variable

    Examples
    ========

    >>> from sympy.physics.vector import ReferenceFrame, CoordinateSym
    >>> A = ReferenceFrame('A')
    >>> A[1]
    A_y
    >>> type(A[0])
    <class 'sympy.physics.vector.frame.CoordinateSym'>
    >>> a_y = CoordinateSym('a_y', A, 1)
    >>> a_y == A[1]
    True

    """

    def __new__(cls, name, frame, index):
        # We can't use the cached Symbol.__new__ because this class depends on
        # frame and index, which are not passed to Symbol.__xnew__.
        assumptions = {}
        super(CoordinateSym, cls)._sanitize(assumptions, cls)
        obj = super(CoordinateSym, cls).__xnew__(cls, name, **assumptions)
        _check_frame(frame)
        if index not in range(0, 3):
            raise ValueError("Invalid index specified")
        obj._id = (frame, index)
        return obj

    @property
    def frame(self):
        return self._id[0]

    def __eq__(self, other):
        #Check if the other object is a CoordinateSym of the same frame
        #and same index
        if isinstance(other, CoordinateSym):
            if other._id == self._id:
                return True
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return tuple((self._id[0].__hash__(), self._id[1])).__hash__()


class ReferenceFrame(object):
    """A reference frame in classical mechanics.

    ReferenceFrame is a class used to represent a reference frame in classical
    mechanics. It has a standard basis of three unit vectors in the frame's
    x, y, and z directions.

    It also can have a rotation relative to a parent frame; this rotation is
    defined by a direction cosine matrix relating this frame's basis vectors to
    the parent frame's basis vectors.  It can also have an angular velocity
    vector, defined in another frame.

    """
    _count = 0

    def __init__(self, name, indices=None, latexs=None, variables=None):
        """ReferenceFrame initialization method.

        A ReferenceFrame has a set of orthonormal basis vectors, along with
        orientations relative to other ReferenceFrames and angular velocities
        relative to other ReferenceFrames.

        Parameters
        ==========

        indices : tuple of str
            Enables the reference frame's basis unit vectors to be accessed by
            Python's square bracket indexing notation using the provided three
            indice strings and alters the printing of the unit vectors to
            reflect this choice.
        latexs : tuple of str
            Alters the LaTeX printing of the reference frame's basis unit
            vectors to the provided three valid LaTeX strings.

        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, vlatex
        >>> N = ReferenceFrame('N')
        >>> N.x
        N.x
        >>> O = ReferenceFrame('O', indices=('1', '2', '3'))
        >>> O.x
        O['1']
        >>> O['1']
        O['1']
        >>> P = ReferenceFrame('P', latexs=('A1', 'A2', 'A3'))
        >>> vlatex(P.x)
        'A1'

        symbols() can be used to create multiple Reference Frames in one step, for example:

        >>> from sympy.physics.vector import ReferenceFrame
        >>> from sympy import symbols
        >>> A, B, C = symbols('A B C', cls=ReferenceFrame)
        >>> D, E = symbols('D E', cls=ReferenceFrame, indices=('1', '2', '3'))
        >>> A[0]
        A_x
        >>> D.x
        D['1']
        >>> E.y
        E['2']
        >>> type(A) == type(D)
        True

        """

        if not isinstance(name, str):
            raise TypeError('Need to supply a valid name')
        # The if statements below are for custom printing of basis-vectors for
        # each frame.
        # First case, when custom indices are supplied
        if indices is not None:
            if not isinstance(indices, (tuple, list)):
                raise TypeError('Supply the indices as a list')
            if len(indices) != 3:
                raise ValueError('Supply 3 indices')
            for i in indices:
                if not isinstance(i, str):
                    raise TypeError('Indices must be strings')
            self.str_vecs = [(name + '[\'' + indices[0] + '\']'),
                             (name + '[\'' + indices[1] + '\']'),
                             (name + '[\'' + indices[2] + '\']')]
            self.pretty_vecs = [(name.lower() + u"_" + indices[0]),
                                (name.lower() + u"_" + indices[1]),
                                (name.lower() + u"_" + indices[2])]
            self.latex_vecs = [(r"\mathbf{\hat{%s}_{%s}}" % (name.lower(),
                               indices[0])), (r"\mathbf{\hat{%s}_{%s}}" %
                               (name.lower(), indices[1])),
                               (r"\mathbf{\hat{%s}_{%s}}" % (name.lower(),
                               indices[2]))]
            self.indices = indices
        # Second case, when no custom indices are supplied
        else:
            self.str_vecs = [(name + '.x'), (name + '.y'), (name + '.z')]
            self.pretty_vecs = [name.lower() + u"_x",
                                name.lower() + u"_y",
                                name.lower() + u"_z"]
            self.latex_vecs = [(r"\mathbf{\hat{%s}_x}" % name.lower()),
                               (r"\mathbf{\hat{%s}_y}" % name.lower()),
                               (r"\mathbf{\hat{%s}_z}" % name.lower())]
            self.indices = ['x', 'y', 'z']
        # Different step, for custom latex basis vectors
        if latexs is not None:
            if not isinstance(latexs, (tuple, list)):
                raise TypeError('Supply the indices as a list')
            if len(latexs) != 3:
                raise ValueError('Supply 3 indices')
            for i in latexs:
                if not isinstance(i, str):
                    raise TypeError('Latex entries must be strings')
            self.latex_vecs = latexs
        self.name = name
        self._var_dict = {}
        #The _dcm_dict dictionary will only store the dcms of parent-child
        #relationships. The _dcm_cache dictionary will work as the dcm
        #cache.
        self._dcm_dict = {}
        self._dcm_cache = {}
        self._ang_vel_dict = {}
        self._ang_acc_dict = {}
        self._dlist = [self._dcm_dict, self._ang_vel_dict, self._ang_acc_dict]
        self._cur = 0
        self._x = Vector([(Matrix([1, 0, 0]), self)])
        self._y = Vector([(Matrix([0, 1, 0]), self)])
        self._z = Vector([(Matrix([0, 0, 1]), self)])
        #Associate coordinate symbols wrt this frame
        if variables is not None:
            if not isinstance(variables, (tuple, list)):
                raise TypeError('Supply the variable names as a list/tuple')
            if len(variables) != 3:
                raise ValueError('Supply 3 variable names')
            for i in variables:
                if not isinstance(i, str):
                    raise TypeError('Variable names must be strings')
        else:
            variables = [name + '_x', name + '_y', name + '_z']
        self.varlist = (CoordinateSym(variables[0], self, 0), \
                        CoordinateSym(variables[1], self, 1), \
                        CoordinateSym(variables[2], self, 2))
        ReferenceFrame._count += 1
        self.index = ReferenceFrame._count

    def __getitem__(self, ind):
        """
        Returns basis vector for the provided index, if the index is a string.

        If the index is a number, returns the coordinate variable correspon-
        -ding to that index.
        """
        if not isinstance(ind, str):
            if ind < 3:
                return self.varlist[ind]
            else:
                raise ValueError("Invalid index provided")
        if self.indices[0] == ind:
            return self.x
        if self.indices[1] == ind:
            return self.y
        if self.indices[2] == ind:
            return self.z
        else:
            raise ValueError('Not a defined index')

    def __iter__(self):
        return iter([self.x, self.y, self.z])

    def __str__(self):
        """Returns the name of the frame. """
        return self.name

    __repr__ = __str__

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
        outlist.sort(key=len)
        if len(outlist) != 0:
            return outlist[0]
        raise ValueError('No Connecting Path found between ' + self.name +
                         ' and ' + other.name)

    def _w_diff_dcm(self, otherframe):
        """Angular velocity from time differentiating the DCM. """
        from sympy.physics.vector.functions import dynamicsymbols
        dcm2diff = otherframe.dcm(self)
        diffed = dcm2diff.diff(dynamicsymbols._t)
        angvelmat = diffed * dcm2diff.T
        w1 = trigsimp(expand(angvelmat[7]), recursive=True)
        w2 = trigsimp(expand(angvelmat[2]), recursive=True)
        w3 = trigsimp(expand(angvelmat[3]), recursive=True)
        return Vector([(Matrix([w1, w2, w3]), otherframe)])

    def variable_map(self, otherframe):
        """
        Returns a dictionary which expresses the coordinate variables
        of this frame in terms of the variables of otherframe.

        If Vector.simp is True, returns a simplified version of the mapped
        values. Else, returns them without simplification.

        Simplification of the expressions may take time.

        Parameters
        ==========

        otherframe : ReferenceFrame
            The other frame to map the variables to

        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, dynamicsymbols
        >>> A = ReferenceFrame('A')
        >>> q = dynamicsymbols('q')
        >>> B = A.orientnew('B', 'Axis', [q, A.z])
        >>> A.variable_map(B)
        {A_x: B_x*cos(q(t)) - B_y*sin(q(t)), A_y: B_x*sin(q(t)) + B_y*cos(q(t)), A_z: B_z}

        """

        _check_frame(otherframe)
        if (otherframe, Vector.simp) in self._var_dict:
            return self._var_dict[(otherframe, Vector.simp)]
        else:
            vars_matrix = self.dcm(otherframe) * Matrix(otherframe.varlist)
            mapping = {}
            for i, x in enumerate(self):
                if Vector.simp:
                    mapping[self.varlist[i]] = trigsimp(vars_matrix[i], method='fu')
                else:
                    mapping[self.varlist[i]] = vars_matrix[i]
            self._var_dict[(otherframe, Vector.simp)] = mapping
            return mapping

    def ang_acc_in(self, otherframe):
        """Returns the angular acceleration Vector of the ReferenceFrame.

        Effectively returns the Vector:
        ^N alpha ^B
        which represent the angular acceleration of B in N, where B is self, and
        N is otherframe.

        Parameters
        ==========

        otherframe : ReferenceFrame
            The ReferenceFrame which the angular acceleration is returned in.

        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_acc(N, V)
        >>> A.ang_acc_in(N)
        10*N.x

        """

        _check_frame(otherframe)
        if otherframe in self._ang_acc_dict:
            return self._ang_acc_dict[otherframe]
        else:
            return self.ang_vel_in(otherframe).dt(otherframe)

    def ang_vel_in(self, otherframe):
        """Returns the angular velocity Vector of the ReferenceFrame.

        Effectively returns the Vector:
        ^N omega ^B
        which represent the angular velocity of B in N, where B is self, and
        N is otherframe.

        Parameters
        ==========

        otherframe : ReferenceFrame
            The ReferenceFrame which the angular velocity is returned in.

        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        10*N.x

        """

        _check_frame(otherframe)
        flist = self._dict_list(otherframe, 1)
        outvec = Vector(0)
        for i in range(len(flist) - 1):
            outvec += flist[i]._ang_vel_dict[flist[i + 1]]
        return outvec

    def dcm(self, otherframe):
        r"""Returns the direction cosine matrix relative to the provided
        reference frame.

        The returned matrix can be used to express the orthogonal unit vectors
        of this frame in terms of the orthogonal unit vectors of
        ``otherframe``.

        Parameters
        ==========

        otherframe : ReferenceFrame
            The reference frame which the direction cosine matrix of this frame
            is formed relative to.

        Examples
        ========

        The following example rotates the reference frame A relative to N by a
        simple rotation and then calculates the direction cosine matrix of N
        relative to A.

        >>> from sympy import symbols, sin, cos
        >>> from sympy.physics.vector import ReferenceFrame
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', (q1, N.x))
        >>> N.dcm(A)
        Matrix([
        [1,       0,        0],
        [0, cos(q1), -sin(q1)],
        [0, sin(q1),  cos(q1)]])

        The second row of the above direction cosine matrix represents the
        ``N.y`` unit vector in N expressed in A. Like so:

        >>> Ny = 0*A.x + cos(q1)*A.y - sin(q1)*A.z

        Thus, expressing ``N.y`` in A should return the same result:

        >>> N.y.express(A)
        cos(q1)*A.y - sin(q1)*A.z

        Notes
        =====

        It is import to know what form of the direction cosine matrix is
        returned. If ``B.dcm(A)`` is called, it means the "direction cosine
        matrix of B relative to A". This is the matrix :math:`{}^A\mathbf{R}^B`
        shown in the following relationship:

        .. math::

           \begin{bmatrix}
             \hat{\mathbf{b}}_1 \\
             \hat{\mathbf{b}}_2 \\
             \hat{\mathbf{b}}_3
           \end{bmatrix}
           =
           {}^A\mathbf{R}^B
           \begin{bmatrix}
             \hat{\mathbf{a}}_1 \\
             \hat{\mathbf{a}}_2 \\
             \hat{\mathbf{a}}_3
           \end{bmatrix}.

        :math:`^{}A\mathbf{R}^B` is the matrix that expresses the B unit
        vectors in terms of the A unit vectors.

        """

        _check_frame(otherframe)
        # Check if the dcm wrt that frame has already been calculated
        if otherframe in self._dcm_cache:
            return self._dcm_cache[otherframe]
        flist = self._dict_list(otherframe, 0)
        outdcm = eye(3)
        for i in range(len(flist) - 1):
            outdcm = outdcm * flist[i]._dcm_dict[flist[i + 1]]
        # After calculation, store the dcm in dcm cache for faster future
        # retrieval
        self._dcm_cache[otherframe] = outdcm
        otherframe._dcm_cache[self] = outdcm.T
        return outdcm

    def orient(self, parent, rot_type, amounts, rot_order=''):
        """Sets the orientation of this reference frame relative to another
        (parent) reference frame.

        Parameters
        ==========

        parent : ReferenceFrame
            Reference frame that this reference frame will be rotated relative
            to.
        rot_type : str
            The method used to generate the direction cosine matrix. Supported
            methods are:

            - ``'Axis'``: simple rotations about a single common axis
            - ``'DCM'``: for setting the direction cosine matrix directly
            - ``'Body'``: three successive rotations about new intermediate
              axes, also called "Euler and Tait-Bryan angles"
            - ``'Space'``: three successive rotations about the parent
              frames' unit vectors
            - ``'Quaternion'``: rotations defined by four parameters which
              result in a singularity free direction cosine matrix

        amounts :
            Expressions defining the rotation angles or direction cosine
            matrix. These must match the ``rot_type``. See examples below for
            details. The input types are:

            - ``'Axis'``: 2-tuple (expr/sym/func, Vector)
            - ``'DCM'``: Matrix, shape(3,3)
            - ``'Body'``: 3-tuple of expressions, symbols, or functions
            - ``'Space'``: 3-tuple of expressions, symbols, or functions
            - ``'Quaternion'``: 4-tuple of expressions, symbols, or
              functions

        rot_order : str or int, optional
            If applicable, the order of the successive of rotations. The string
            ``'123'`` and integer ``123`` are equivalent, for example. Required
            for ``'Body'`` and ``'Space'``.

        Examples
        ========

        Setup variables for the examples:

        >>> from sympy import symbols
        >>> from sympy.physics.vector import ReferenceFrame
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = ReferenceFrame('N')
        >>> B = ReferenceFrame('B')
        >>> B1 = ReferenceFrame('B')
        >>> B2 = ReferenceFrame('B2')

        Axis
        ----

        ``rot_type='Axis'`` creates a direction cosine matrix defined by a
        simple rotation about a single axis fixed in both reference frames.
        This is a rotation about an arbitrary, non-time-varying
        axis by some angle. The axis is supplied as a Vector. This is how
        simple rotations are defined.

        >>> B.orient(N, 'Axis', (q1, N.x))

        The ``orient()`` method generates a direction cosine matrix and its
        transpose which defines the orientation of B relative to N and vice
        versa. Once orient is called, ``dcm()`` outputs the appropriate
        direction cosine matrix.

        >>> B.dcm(N)
        Matrix([
        [1,       0,      0],
        [0,  cos(q1), sin(q1)],
        [0, -sin(q1), cos(q1)]])

        The following two lines show how the sense of the rotation can be
        defined. Both lines produce the same result.

        >>> B.orient(N, 'Axis', (q1, -N.x))
        >>> B.orient(N, 'Axis', (-q1, N.x))

        The axis does not have to be defined by a unit vector, it can be any
        vector in the parent frame.

        >>> B.orient(N, 'Axis', (q1, N.x + 2 * N.y))

        DCM
        ---

        The direction cosine matrix can be set directly. The orientation of a
        frame A can be set to be the same as the frame B above like so:

        >>> B.orient(N, 'Axis', (q1, N.x))
        >>> A = ReferenceFrame('A')
        >>> A.orient(N, 'DCM', N.dcm(B))
        >>> A.dcm(N)
        Matrix([
        [1,       0,      0],
        [0,  cos(q1), sin(q1)],
        [0, -sin(q1), cos(q1)]])

        **Note carefully that** ``N.dcm(B)`` **was passed into** ``orient()``
        **for** ``A.dcm(N)`` **to match** ``B.dcm(N)``.

        Body
        ----

        ``rot_type='Body'`` rotates this reference frame relative to the
        provided reference frame by rotating through three successive simple
        rotations.  Each subsequent axis of rotation is about the "body fixed"
        unit vectors of the new intermediate reference frame. This type of
        rotation is also referred to rotating through the `Euler and Tait-Bryan
        Angles <https://en.wikipedia.org/wiki/Euler_angles>`_.

        For example, the classic Euler Angle rotation can be done by:

        >>> B.orient(N, 'Body', (q1, q2, q3), 'XYX')
        >>> B.dcm(N)
        Matrix([
        [        cos(q2),                            sin(q1)*sin(q2),                           -sin(q2)*cos(q1)],
        [sin(q2)*sin(q3), -sin(q1)*sin(q3)*cos(q2) + cos(q1)*cos(q3),  sin(q1)*cos(q3) + sin(q3)*cos(q1)*cos(q2)],
        [sin(q2)*cos(q3), -sin(q1)*cos(q2)*cos(q3) - sin(q3)*cos(q1), -sin(q1)*sin(q3) + cos(q1)*cos(q2)*cos(q3)]])

        This rotates B relative to N through ``q1`` about ``N.x``, then rotates
        B again through q2 about B.y, and finally through q3 about B.x. It is
        equivalent to:

        >>> B1.orient(N, 'Axis', (q1, N.x))
        >>> B2.orient(B1, 'Axis', (q2, B1.y))
        >>> B.orient(B2, 'Axis', (q3, B2.x))
        >>> B.dcm(N)
        Matrix([
        [        cos(q2),                            sin(q1)*sin(q2),                           -sin(q2)*cos(q1)],
        [sin(q2)*sin(q3), -sin(q1)*sin(q3)*cos(q2) + cos(q1)*cos(q3),  sin(q1)*cos(q3) + sin(q3)*cos(q1)*cos(q2)],
        [sin(q2)*cos(q3), -sin(q1)*cos(q2)*cos(q3) - sin(q3)*cos(q1), -sin(q1)*sin(q3) + cos(q1)*cos(q2)*cos(q3)]])

        Acceptable rotation orders are of length 3, expressed in as a string
        ``'XYZ'`` or ``'123'`` or integer ``123``. Rotations about an axis
        twice in a row are prohibited.

        >>> B.orient(N, 'Body', (q1, q2, 0), 'ZXZ')
        >>> B.orient(N, 'Body', (q1, q2, 0), '121')
        >>> B.orient(N, 'Body', (q1, q2, q3), 123)

        Space
        -----

        ``rot_type='Space'`` also rotates the reference frame in three
        successive simple rotations but the axes of rotation are the
        "Space-fixed" axes. For example:

        >>> B.orient(N, 'Space', (q1, q2, q3), '312')
        >>> B.dcm(N)
        Matrix([
        [ sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3), sin(q1)*cos(q2), sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1)],
        [-sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1), cos(q1)*cos(q2), sin(q1)*sin(q3) + sin(q2)*cos(q1)*cos(q3)],
        [                           sin(q3)*cos(q2),        -sin(q2),                           cos(q2)*cos(q3)]])

        is equivalent to:

        >>> B1.orient(N, 'Axis', (q1, N.z))
        >>> B2.orient(B1, 'Axis', (q2, N.x))
        >>> B.orient(B2, 'Axis', (q3, N.y))
        >>> B.dcm(N).simplify()  # doctest: +SKIP
        Matrix([
        [ sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3), sin(q1)*cos(q2), sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1)],
        [-sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1), cos(q1)*cos(q2), sin(q1)*sin(q3) + sin(q2)*cos(q1)*cos(q3)],
        [                           sin(q3)*cos(q2),        -sin(q2),                           cos(q2)*cos(q3)]])

        It is worth noting that space-fixed and body-fixed rotations are
        related by the order of the rotations, i.e. the reverse order of body
        fixed will give space fixed and vice versa.

        >>> B.orient(N, 'Space', (q1, q2, q3), '231')
        >>> B.dcm(N)
        Matrix([
        [cos(q1)*cos(q2), sin(q1)*sin(q3) + sin(q2)*cos(q1)*cos(q3), -sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1)],
        [       -sin(q2),                           cos(q2)*cos(q3),                            sin(q3)*cos(q2)],
        [sin(q1)*cos(q2), sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1),  sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3)]])

        >>> B.orient(N, 'Body', (q3, q2, q1), '132')
        >>> B.dcm(N)
        Matrix([
        [cos(q1)*cos(q2), sin(q1)*sin(q3) + sin(q2)*cos(q1)*cos(q3), -sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1)],
        [       -sin(q2),                           cos(q2)*cos(q3),                            sin(q3)*cos(q2)],
        [sin(q1)*cos(q2), sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1),  sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3)]])

        Quaternion
        ----------

        ``rot_type='Quaternion'`` orients the reference frame using
        quaternions. Quaternion rotation is defined as a finite rotation about
        lambda, a unit vector, by an amount theta. This orientation is
        described by four parameters:

        - ``q0 = cos(theta/2)``
        - ``q1 = lambda_x sin(theta/2)``
        - ``q2 = lambda_y sin(theta/2)``
        - ``q3 = lambda_z sin(theta/2)``

        This type does not need a ``rot_order``.

        >>> B.orient(N, 'Quaternion', (q0, q1, q2, q3))
        >>> B.dcm(N)
        Matrix([
        [q0**2 + q1**2 - q2**2 - q3**2,             2*q0*q3 + 2*q1*q2,            -2*q0*q2 + 2*q1*q3],
        [           -2*q0*q3 + 2*q1*q2, q0**2 - q1**2 + q2**2 - q3**2,             2*q0*q1 + 2*q2*q3],
        [            2*q0*q2 + 2*q1*q3,            -2*q0*q1 + 2*q2*q3, q0**2 - q1**2 - q2**2 + q3**2]])

        """

        from sympy.physics.vector.functions import dynamicsymbols
        _check_frame(parent)

        # Allow passing a rotation matrix manually.
        if rot_type == 'DCM':
            # When rot_type == 'DCM', then amounts must be a Matrix type object
            # (e.g. sympy.matrices.dense.MutableDenseMatrix).
            if not isinstance(amounts, MatrixBase):
                raise TypeError("Amounts must be a sympy Matrix type object.")
        else:
            amounts = list(amounts)
            for i, v in enumerate(amounts):
                if not isinstance(v, Vector):
                    amounts[i] = sympify(v)

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
                           '131', '212', '232', '313', '323', '')
        # make sure XYZ => 123 and rot_type is in upper case
        rot_order = translate(str(rot_order), 'XYZxyz', '123123')
        rot_type = rot_type.upper()
        if rot_order not in approved_orders:
            raise TypeError('The supplied order is not an approved type')
        parent_orient = []
        if rot_type == 'AXIS':
            if not rot_order == '':
                raise TypeError('Axis orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 2)):
                raise TypeError('Amounts are a list or tuple of length 2')
            theta = amounts[0]
            axis = amounts[1]
            axis = _check_vector(axis)
            if not axis.dt(parent) == 0:
                raise ValueError('Axis cannot be time-varying')
            axis = axis.express(parent).normalize()
            axis = axis.args[0][0]
            parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                             Matrix([[0, -axis[2], axis[1]],
                                     [axis[2], 0, -axis[0]],
                                     [-axis[1], axis[0], 0]]) *
                             sin(theta) + axis * axis.T)
        elif rot_type == 'QUATERNION':
            if not rot_order == '':
                raise TypeError(
                    'Quaternion orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 4)):
                raise TypeError('Amounts are a list or tuple of length 4')
            q0, q1, q2, q3 = amounts
            parent_orient = (Matrix([[q0**2 + q1**2 - q2**2 - q3**2,
                                      2 * (q1 * q2 - q0 * q3),
                                      2 * (q0 * q2 + q1 * q3)],
                                     [2 * (q1 * q2 + q0 * q3),
                                      q0**2 - q1**2 + q2**2 - q3**2,
                                      2 * (q2 * q3 - q0 * q1)],
                                     [2 * (q1 * q3 - q0 * q2),
                                      2 * (q0 * q1 + q2 * q3),
                                      q0**2 - q1**2 - q2**2 + q3**2]]))
        elif rot_type == 'BODY':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Body orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            parent_orient = (_rot(a1, amounts[0]) * _rot(a2, amounts[1]) *
                             _rot(a3, amounts[2]))
        elif rot_type == 'SPACE':
            if not (len(amounts) == 3 & len(rot_order) == 3):
                raise TypeError('Space orientation takes 3 values & 3 orders')
            a1 = int(rot_order[0])
            a2 = int(rot_order[1])
            a3 = int(rot_order[2])
            parent_orient = (_rot(a3, amounts[2]) * _rot(a2, amounts[1]) *
                             _rot(a1, amounts[0]))
        elif rot_type == 'DCM':
            parent_orient = amounts
        else:
            raise NotImplementedError('That is not an implemented rotation')
        # Reset the _dcm_cache of this frame, and remove it from the
        # _dcm_caches of the frames it is linked to. Also remove it from the
        # _dcm_dict of its parent
        frames = self._dcm_cache.keys()
        dcm_dict_del = []
        dcm_cache_del = []
        for frame in frames:
            if frame in self._dcm_dict:
                dcm_dict_del += [frame]
            dcm_cache_del += [frame]
        for frame in dcm_dict_del:
            del frame._dcm_dict[self]
        for frame in dcm_cache_del:
            del frame._dcm_cache[self]
        # Add the dcm relationship to _dcm_dict
        self._dcm_dict = self._dlist[0] = {}
        self._dcm_dict.update({parent: parent_orient.T})
        parent._dcm_dict.update({self: parent_orient})
        # Also update the dcm cache after resetting it
        self._dcm_cache = {}
        self._dcm_cache.update({parent: parent_orient.T})
        parent._dcm_cache.update({self: parent_orient})
        if rot_type == 'QUATERNION':
            t = dynamicsymbols._t
            q0, q1, q2, q3 = amounts
            q0d = diff(q0, t)
            q1d = diff(q1, t)
            q2d = diff(q2, t)
            q3d = diff(q3, t)
            w1 = 2 * (q1d * q0 + q2d * q3 - q3d * q2 - q0d * q1)
            w2 = 2 * (q2d * q0 + q3d * q1 - q1d * q3 - q0d * q2)
            w3 = 2 * (q3d * q0 + q1d * q2 - q2d * q1 - q0d * q3)
            wvec = Vector([(Matrix([w1, w2, w3]), self)])
        elif rot_type == 'AXIS':
            thetad = (amounts[0]).diff(dynamicsymbols._t)
            wvec = thetad * amounts[1].express(parent).normalize()
        elif rot_type == 'DCM':
            wvec = self._w_diff_dcm(parent)
        else:
            try:
                from sympy.polys.polyerrors import CoercionFailed
                from sympy.physics.vector.functions import kinematic_equations
                q1, q2, q3 = amounts
                u1, u2, u3 = symbols('u1, u2, u3', cls=Dummy)
                templist = kinematic_equations([u1, u2, u3], [q1, q2, q3],
                                               rot_type, rot_order)
                templist = [expand(i) for i in templist]
                td = solve(templist, [u1, u2, u3])
                u1 = expand(td[u1])
                u2 = expand(td[u2])
                u3 = expand(td[u3])
                wvec = u1 * self.x + u2 * self.y + u3 * self.z
            except (CoercionFailed, AssertionError):
                wvec = self._w_diff_dcm(parent)
        self._ang_vel_dict.update({parent: wvec})
        parent._ang_vel_dict.update({self: -wvec})
        self._var_dict = {}

    def orientnew(self, newname, rot_type, amounts, rot_order='',
                  variables=None, indices=None, latexs=None):
        r"""Returns a new reference frame oriented with respect to this
        reference frame.

        See ``ReferenceFrame.orient()`` for detailed examples of how to orient
        reference frames.

        Parameters
        ==========

        newname : str
            Name for the new reference frame.
        rot_type : str
            The method used to generate the direction cosine matrix. Supported
            methods are:

            - ``'Axis'``: simple rotations about a single common axis
            - ``'DCM'``: for setting the direction cosine matrix directly
            - ``'Body'``: three successive rotations about new intermediate
              axes, also called "Euler and Tait-Bryan angles"
            - ``'Space'``: three successive rotations about the parent
              frames' unit vectors
            - ``'Quaternion'``: rotations defined by four parameters which
              result in a singularity free direction cosine matrix

        amounts :
            Expressions defining the rotation angles or direction cosine
            matrix. These must match the ``rot_type``. See examples below for
            details. The input types are:

            - ``'Axis'``: 2-tuple (expr/sym/func, Vector)
            - ``'DCM'``: Matrix, shape(3,3)
            - ``'Body'``: 3-tuple of expressions, symbols, or functions
            - ``'Space'``: 3-tuple of expressions, symbols, or functions
            - ``'Quaternion'``: 4-tuple of expressions, symbols, or
              functions

        rot_order : str or int, optional
            If applicable, the order of the successive of rotations. The string
            ``'123'`` and integer ``123`` are equivalent, for example. Required
            for ``'Body'`` and ``'Space'``.
        indices : tuple of str
            Enables the reference frame's basis unit vectors to be accessed by
            Python's square bracket indexing notation using the provided three
            indice strings and alters the printing of the unit vectors to
            reflect this choice.
        latexs : tuple of str
            Alters the LaTeX printing of the reference frame's basis unit
            vectors to the provided three valid LaTeX strings.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.vector import ReferenceFrame, vlatex
        >>> q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
        >>> N = ReferenceFrame('N')

        Create a new reference frame A rotated relative to N through a simple
        rotation.

        >>> A = N.orientnew('A', 'Axis', (q0, N.x))

        Create a new reference frame B rotated relative to N through body-fixed
        rotations.

        >>> B = N.orientnew('B', 'Body', (q1, q2, q3), '123')

        Create a new reference frame C rotated relative to N through a simple
        rotation with unique indices and LaTeX printing.

        >>> C = N.orientnew('C', 'Axis', (q0, N.x), indices=('1', '2', '3'),
        ... latexs=(r'\hat{\mathbf{c}}_1',r'\hat{\mathbf{c}}_2',
        ... r'\hat{\mathbf{c}}_3'))
        >>> C['1']
        C['1']
        >>> print(vlatex(C['1']))
        \hat{\mathbf{c}}_1

        """

        newframe = self.__class__(newname, variables=variables,
                                  indices=indices, latexs=latexs)
        newframe.orient(self, rot_type, amounts, rot_order)
        return newframe

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

        >>> from sympy.physics.vector import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_acc(N, V)
        >>> A.ang_acc_in(N)
        10*N.x

        """

        if value == 0:
            value = Vector(0)
        value = _check_vector(value)
        _check_frame(otherframe)
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

        >>> from sympy.physics.vector import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        10*N.x

        """

        if value == 0:
            value = Vector(0)
        value = _check_vector(value)
        _check_frame(otherframe)
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

    def partial_velocity(self, frame, *gen_speeds):
        """Returns the partial angular velocities of this frame in the given
        frame with respect to one or more provided generalized speeds.

        Parameters
        ==========
        frame : ReferenceFrame
            The frame with which the angular velocity is defined in.
        gen_speeds : functions of time
            The generalized speeds.

        Returns
        =======
        partial_velocities : tuple of Vector
            The partial angular velocity vectors corresponding to the provided
            generalized speeds.

        Examples
        ========

        >>> from sympy.physics.vector import ReferenceFrame, dynamicsymbols
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> u1, u2 = dynamicsymbols('u1, u2')
        >>> A.set_ang_vel(N, u1 * A.x + u2 * N.y)
        >>> A.partial_velocity(N, u1)
        A.x
        >>> A.partial_velocity(N, u1, u2)
        (A.x, N.y)

        """

        partials = [self.ang_vel_in(frame).diff(speed, frame, var_in_dcm=False)
                    for speed in gen_speeds]

        if len(partials) == 1:
            return partials[0]
        else:
            return tuple(partials)


def _check_frame(other):
    from .vector import VectorTypeError
    if not isinstance(other, ReferenceFrame):
        raise VectorTypeError(other, ReferenceFrame('A'))
