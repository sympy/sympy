__all__ = ['ReferenceFrame', 'Vector', 'Dyadic', 'dynamicsymbols',
           'MechanicsStrPrinter', 'MechanicsPrettyPrinter',
           'MechanicsLatexPrinter']

from sympy import (
    Symbol, sin, cos, eye, trigsimp, diff, sqrt, sympify,
    expand, zeros, Derivative, Function, symbols, Add,
    solve, S, ImmutableMatrix as Matrix)
from sympy.core import C
from sympy.core.compatibility import reduce, u, string_types
from sympy.core.function import UndefinedFunction
from sympy.printing.conventions import split_super_sub
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.printing.str import StrPrinter
from sympy.utilities import group


class Dyadic(object):
    """A Dyadic object.

    See:
    http://en.wikipedia.org/wiki/Dyadic_tensor
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    A more powerful way to represent a rigid body's inertia. While it is more
    complex, by choosing Dyadic components to be in body fixed basis vectors,
    the resulting matrix is equivalent to the inertia tensor.

    """

    def __init__(self, inlist):
        """
        Just like Vector's init, you shouldn't call this.

        Stores a Dyadic as a list of lists; the inner list has the measure number
        and the two unit vectors; the outerlist holds each unique unit vector
        pair.

        """

        self.args = []
        while len(inlist) != 0:
            added = 0
            for i, v in enumerate(self.args):
                if ((str(inlist[0][1]) == str(self.args[i][1])) and
                        (str(inlist[0][2]) == str(self.args[i][2]))):
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

    def __add__(self, other):
        """The add operator for Dyadic. """
        other = _check_dyadic(other)
        return Dyadic(self.args + other.args)

    def __and__(self, other):
        """The inner product operator for a Dyadic and a Dyadic or Vector.

        Parameters
        ==========

        other : Dyadic or Vector
            The other Dyadic or Vector to take the inner product with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> D1 = outer(N.x, N.y)
        >>> D2 = outer(N.y, N.y)
        >>> D1.dot(D2)
        (N.x|N.y)
        >>> D1.dot(N.y)
        N.x

        """

        if isinstance(other, Dyadic):
            other = _check_dyadic(other)
            ol = Dyadic([])
            for i, v in enumerate(self.args):
                for i2, v2 in enumerate(other.args):
                    ol += v[0] * v2[0] * (v[2] & v2[1]) * (v[1] | v2[2])
        else:
            other = _check_vector(other)
            ol = Vector([])
            for i, v in enumerate(self.args):
                ol += v[0] * v[1] * (v[2] & other)
        return ol

    def __div__(self, other):
        """Divides the Dyadic by a sympifyable expression. """
        return self.__mul__(1 / other)

    __truediv__ = __div__

    def __eq__(self, other):
        """Tests for equality.

        Is currently weak; needs stronger comparison testing

        """

        other = _check_dyadic(other)
        if (self.args == []) and (other.args == []):
            return True
        elif (self.args == []) or (other.args == []):
            return False
        return set(self.args) == set(other.args)

    def __mul__(self, other):
        """Multiplies the Dyadic by a sympifyable expression.

        Parameters
        ==========

        other : Sympafiable
            The scalar to multiply this Dyadic with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> 5 * d
        5*(N.x|N.x)

        """

        newlist = [v for v in self.args]
        for i, v in enumerate(newlist):
            newlist[i] = (sympify(other) * newlist[i][0], newlist[i][1],
                    newlist[i][2])
        return Dyadic(newlist)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __neg__(self):
        return self * -1

    def _latex(self, printer=None):
        ar = self.args  # just to shorten things
        if len(ar) == 0:
            return str(0)
        ol = []  # output list, to be concatenated to a string
        mlp = MechanicsLatexPrinter()
        for i, v in enumerate(ar):
            # if the coef of the dyadic is 1, we skip the 1
            if ar[i][0] == 1:
                ol.append(' + ' + mlp.doprint(ar[i][1]) + r"\otimes " +
                        mlp.doprint(ar[i][2]))
            # if the coef of the dyadic is -1, we skip the 1
            elif ar[i][0] == -1:
                ol.append(' - ' +
                          mlp.doprint(ar[i][1]) +
                          r"\otimes " +
                          mlp.doprint(ar[i][2]))
            # If the coefficient of the dyadic is not 1 or -1,
            # we might wrap it in parentheses, for readability.
            elif ar[i][0] != 0:
                arg_str = mlp.doprint(ar[i][0])
                if isinstance(ar[i][0], Add):
                    arg_str = '(%s)' % arg_str
                if arg_str.startswith('-'):
                    arg_str = arg_str[1:]
                    str_start = ' - '
                else:
                    str_start = ' + '
                ol.append(str_start + arg_str + r" " +
                          mlp.doprint(ar[i][1]) +
                          r"\otimes " +
                          mlp.doprint(ar[i][2]))
        outstr = ''.join(ol)
        if outstr.startswith(' + '):
            outstr = outstr[3:]
        elif outstr.startswith(' '):
            outstr = outstr[1:]
        return outstr

    def _pretty(self, printer=None):
        e = self

        class Fake(object):
            baseline = 0

            def render(self, *args, **kwargs):
                self = e
                ar = self.args  # just to shorten things
                mpp = MechanicsPrettyPrinter()
                if len(ar) == 0:
                    return unicode(0)
                ol = []  # output list, to be concatenated to a string
                for i, v in enumerate(ar):
                    # if the coef of the dyadic is 1, we skip the 1
                    if ar[i][0] == 1:
                        ol.append(u(" + ") +
                                  mpp.doprint(ar[i][1]) +
                                  u("\u2a02 ") +
                                  mpp.doprint(ar[i][2]))
                    # if the coef of the dyadic is -1, we skip the 1
                    elif ar[i][0] == -1:
                        ol.append(u(" - ") +
                                  mpp.doprint(ar[i][1]) +
                                  u("\u2a02 ") +
                                  mpp.doprint(ar[i][2]))
                    # If the coefficient of the dyadic is not 1 or -1,
                    # we might wrap it in parentheses, for readability.
                    elif ar[i][0] != 0:
                        arg_str = mpp.doprint(ar[i][0])
                        if isinstance(ar[i][0], Add):
                            arg_str = u("(%s)") % arg_str
                        if arg_str.startswith(u("-")):
                            arg_str = arg_str[1:]
                            str_start = u(" - ")
                        else:
                            str_start = u(" + ")
                        ol.append(str_start + arg_str + u(" ") +
                                  mpp.doprint(ar[i][1]) +
                                  u("\u2a02 ") +
                                  mpp.doprint(ar[i][2]))
                outstr = u("").join(ol)
                if outstr.startswith(u(" + ")):
                    outstr = outstr[3:]
                elif outstr.startswith(" "):
                    outstr = outstr[1:]
                return outstr
        return Fake()

    def __rand__(self, other):
        """The inner product operator for a Vector or Dyadic, and a Dyadic

        This is for: Vector dot Dyadic

        Parameters
        ==========

        other : Vector
            The vector we are dotting with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, dot, outer
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> dot(N.x, d)
        N.x

        """

        other = _check_vector(other)
        ol = Vector([])
        for i, v in enumerate(self.args):
            ol += v[0] * v[2] * (v[1] & other)
        return ol

    def __rsub__(self, other):
        return (-1 * self) + other

    def __rxor__(self, other):
        """For a cross product in the form: Vector x Dyadic

        Parameters
        ==========

        other : Vector
            The Vector that we are crossing this Dyadic with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(N.y, d)
        - (N.z|N.x)

        """

        other = _check_vector(other)
        ol = Dyadic([])
        for i, v in enumerate(self.args):
            ol += v[0] * ((other ^ v[1]) | v[2])
        return ol

    def __str__(self, printer=None):
        """Printing method. """
        ar = self.args  # just to shorten things
        if len(ar) == 0:
            return str(0)
        ol = []  # output list, to be concatenated to a string
        for i, v in enumerate(ar):
            # if the coef of the dyadic is 1, we skip the 1
            if ar[i][0] == 1:
                ol.append(' + (' + str(ar[i][1]) + '|' + str(ar[i][2]) + ')')
            # if the coef of the dyadic is -1, we skip the 1
            elif ar[i][0] == -1:
                ol.append(' - (' + str(ar[i][1]) + '|' + str(ar[i][2]) + ')')
            # If the coefficient of the dyadic is not 1 or -1,
            # we might wrap it in parentheses, for readability.
            elif ar[i][0] != 0:
                arg_str = MechanicsStrPrinter().doprint(ar[i][0])
                if isinstance(ar[i][0], Add):
                    arg_str = "(%s)" % arg_str
                if arg_str[0] == '-':
                    arg_str = arg_str[1:]
                    str_start = ' - '
                else:
                    str_start = ' + '
                ol.append(str_start + arg_str + '*(' + str(ar[i][1]) +
                          '|' + str(ar[i][2]) + ')')
        outstr = ''.join(ol)
        if outstr.startswith(' + '):
            outstr = outstr[3:]
        elif outstr.startswith(' '):
            outstr = outstr[1:]
        return outstr

    def __sub__(self, other):
        """The subtraction operator. """
        return self.__add__(other * -1)

    def __xor__(self, other):
        """For a cross product in the form: Dyadic x Vector.

        Parameters
        ==========

        other : Vector
            The Vector that we are crossing this Dyadic with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(d, N.y)
        (N.x|N.z)

        """

        other = _check_vector(other)
        ol = Dyadic([])
        for i, v in enumerate(self.args):
            ol += v[0] * (v[1] | (v[2] ^ other))
        return ol

    _sympystr = __str__
    _sympyrepr = _sympystr
    __repr__ = __str__
    __radd__ = __add__
    __rmul__ = __mul__

    def express(self, frame1, frame2=None):
        """Expresses this Dyadic in alternate frame(s)

        The first frame is the list side expression, the second frame is the
        right side; if Dyadic is in form A.x|B.y, you can express it in two
        different frames. If no second frame is given, the Dyadic is expressed in
        only one frame.

        Parameters
        ==========

        frame1 : ReferenceFrame
            The frame to express the left side of the Dyadic in
        frame2 : ReferenceFrame
            If provided, the frame to express the right side of the Dyadic in

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, dynamicsymbols
        >>> N = ReferenceFrame('N')
        >>> q = dynamicsymbols('q')
        >>> B = N.orientnew('B', 'Axis', [q, N.z])
        >>> d = outer(N.x, N.x)
        >>> d.express(B, N)
        cos(q)*(B.x|N.x) - sin(q)*(B.y|N.x)

        """

        if frame2 is None:
            frame2 = frame1
        _check_frame(frame1)
        _check_frame(frame2)
        out = Dyadic([])
        for i, v in enumerate(self.args):
            out += v[0] * (v[1].express(frame1) | v[2].express(frame2))
        return out

    def doit(self, **hints):
        """Calls .doit() on each term in the Dyadic"""
        return sum([Dyadic( [ (v[0].doit(**hints), v[1], v[2]) ]) for
                    v in self.args])

    def dt(self, frame):
        """Take the time derivative of this Dyadic in a frame.

        Parameters
        ==========

        frame : ReferenceFrame
            The frame to take the time derivative in

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, dynamicsymbols
        >>> N = ReferenceFrame('N')
        >>> q = dynamicsymbols('q')
        >>> B = N.orientnew('B', 'Axis', [q, N.z])
        >>> d = outer(N.x, N.x)
        >>> d.dt(B)
        - q'*(N.y|N.x) - q'*(N.x|N.y)

        """

        _check_frame(frame)
        t = dynamicsymbols._t
        ol = S(0)
        for i, v in enumerate(self.args):
            ol += (v[0].diff(t) * (v[1] | v[2]))
            ol += (v[0] * (v[1].dt(frame) | v[2]))
            ol += (v[0] * (v[1] | v[2].dt(frame)))
        return ol

    def simplify(self):
        """Returns a simplified Dyadic."""
        out = Dyadic([])
        for v in self.args:
            out += Dyadic([(v[0].simplify(), v[1], v[2])])
        return out

    def subs(self, *args, **kwargs):
        """Substituion on the Dyadic.

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame
        >>> from sympy import Symbol
        >>> N = ReferenceFrame('N')
        >>> s = Symbol('s')
        >>> a = s * (N.x|N.x)
        >>> a.subs({s: 2})
        2*(N.x|N.x)

        """

        return sum([ Dyadic([(v[0].subs(*args, **kwargs), v[1], v[2])])
                     for v in self.args])

    dot = __and__
    cross = __xor__


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

    def __init__(self, name, indices=None, latexs=None):
        """ReferenceFrame initialization method.

        A ReferenceFrame has a set of orthonormal basis vectors, along with
        orientations relative to other ReferenceFrames and angular velocities
        relative to other ReferenceFrames.

        Parameters
        ==========

        indices : list (of strings)
            If custom indices are desired for console, pretty, and LaTeX
            printing, supply three as a list. The basis vectors can then be
            accessed with the get_item method.
        latexs : list (of strings)
            If custom names are desired for LaTeX printing of each basis
            vector, supply the names here in a list.

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, mlatex
        >>> N = ReferenceFrame('N')
        >>> N.x
        N.x
        >>> O = ReferenceFrame('O', ('1', '2', '3'))
        >>> O.x
        O['1']
        >>> O['1']
        O['1']
        >>> P = ReferenceFrame('P', latexs=('A1', 'A2', 'A3'))
        >>> mlatex(P.x)
        'A1'

        """

        if not isinstance(name, string_types):
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
                if not isinstance(i, string_types):
                    raise TypeError('Indices must be strings')
            self.str_vecs = [(name + '[\'' + indices[0] + '\']'),
                             (name + '[\'' + indices[1] + '\']'),
                             (name + '[\'' + indices[2] + '\']')]
            self.pretty_vecs = [(u("\033[94m\033[1m") + name.lower() + u("_") +
                                indices[0] + u("\033[0;0m\x1b[0;0m")),
                                (u("\033[94m\033[1m") + name.lower() + u("_") +
                                indices[1] + u("\033[0;0m\x1b[0;0m")),
                                (u("\033[94m\033[1m") + name.lower() + u("_") +
                                indices[2] + u("\033[0;0m\x1b[0;0m"))]
            self.latex_vecs = [(r"\mathbf{\hat{%s}_{%s}}" % (name.lower(),
                               indices[0])), (r"\mathbf{\hat{%s}_{%s}}" %
                               (name.lower(), indices[1])),
                               (r"\mathbf{\hat{%s}_{%s}}" % (name.lower(),
                               indices[2]))]
            self.indices = indices
        # Second case, when no custom indices are supplied
        else:
            self.str_vecs = [(name + '.x'), (name + '.y'), (name + '.z')]
            self.pretty_vecs = [(u("\033[94m\033[1m") + name.lower() +
                                u("_x\033[0;0m\x1b[0;0m")),
                                (u("\033[94m\033[1m") + name.lower() +
                                u("_y\033[0;0m\x1b[0;0m")),
                                (u("\033[94m\033[1m") + name.lower() +
                                u("_z\033[0;0m\x1b[0;0m"))]
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
                if not isinstance(i, string_types):
                    raise TypeError('Latex entries must be strings')
            self.latex_vecs = latexs
        self.name = name
        self._dcm_dict = {}
        self._ang_vel_dict = {}
        self._ang_acc_dict = {}
        self._dlist = [self._dcm_dict, self._ang_vel_dict, self._ang_acc_dict]
        self._cur = 0
        self._x = Vector([(Matrix([1, 0, 0]), self)])
        self._y = Vector([(Matrix([0, 1, 0]), self)])
        self._z = Vector([(Matrix([0, 0, 1]), self)])

    def __getitem__(self, ind):
        """Returns basis vector for the provided index (index being an str)"""
        if not isinstance(ind, string_types):
            raise TypeError('Supply a valid str for the index')
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
        dcm2diff = self.dcm(otherframe)
        diffed = dcm2diff.diff(dynamicsymbols._t)
        angvelmat = diffed * dcm2diff.T
        w1 = trigsimp(expand(angvelmat[7]), recursive=True)
        w2 = trigsimp(expand(angvelmat[2]), recursive=True)
        w3 = trigsimp(expand(angvelmat[3]), recursive=True)
        return -Vector([(Matrix([w1, w2, w3]), self)])

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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        10*N.x

        """

        _check_frame(otherframe)
        flist = self._dict_list(otherframe, 1)
        outvec = Vector([])
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])
        >>> N.dcm(A)
        Matrix([
        [1,       0,        0],
        [0, cos(q1), -sin(q1)],
        [0, sin(q1),  cos(q1)]])

        """

        _check_frame(otherframe)
        flist = self._dict_list(otherframe, 0)
        outdcm = eye(3)
        for i in range(len(flist) - 1):
            outdcm = outdcm * flist[i + 1]._dcm_dict[flist[i]]
        return outdcm

    def orient(self, parent, rot_type, amounts, rot_order=''):
        """Defines the orientation of this frame relative to a parent frame.

        Parameters
        ==========

        parent : ReferenceFrame
            The frame that this ReferenceFrame will have its orientation matrix
            defined in relation to.
        rot_type : str
            The type of orientation matrix that is being created. Supported
            types are 'Body', 'Space', 'Quaternion', and 'Axis'. See examples
            for correct usage.
        amounts : list OR value
            The quantities that the orientation matrix will be defined by.
        rot_order : str
            If applicable, the order of a series of rotations.

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q0, q1, q2, q3, q4 = symbols('q0 q1 q2 q3 q4')
        >>> N = ReferenceFrame('N')
        >>> B = ReferenceFrame('B')

        Now we have a choice of how to implement the orientation. First is
        Body. Body orientation takes this reference frame through three
        successive simple rotations. Acceptable rotation orders are of length
        3, expressed in XYZ or 123, and cannot have a rotation about about an
        axis twice in a row.

        >>> B.orient(N, 'Body', [q1, q2, q3], '123')
        >>> B.orient(N, 'Body', [q1, q2, 0], 'ZXZ')
        >>> B.orient(N, 'Body', [0, 0, 0], 'XYX')

        Next is Space. Space is like Body, but the rotations are applied in the
        opposite order.

        >>> B.orient(N, 'Space', [q1, q2, q3], '312')

        Next is Quaternion. This orients the new ReferenceFrame with
        Quaternions, defined as a finite rotation about lambda, a unit vector,
        by some amount theta.
        This orientation is described by four parameters:
        q0 = cos(theta/2)
        q1 = lambda_x sin(theta/2)
        q2 = lambda_y sin(theta/2)
        q3 = lambda_z sin(theta/2)
        Quaternion does not take in a rotation order.

        >>> B.orient(N, 'Quaternion', [q0, q1, q2, q3])

        Last is Axis. This is a rotation about an arbitrary, non-time-varying
        axis by some angle. The axis is supplied as a Vector. This is how
        simple rotations are defined.

        >>> B.orient(N, 'Axis', [q1, N.x + 2 * N.y])

        """

        _check_frame(parent)
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
        rot_order = str(
            rot_order).upper()  # Now we need to make sure XYZ = 123
        rot_type = rot_type.upper()
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
            axis = _check_vector(axis)
            if not axis.dt(parent) == 0:
                raise ValueError('Axis cannot be time-varying')
            axis = axis.express(parent).normalize()
            axis = axis.args[0][0]
            parent_orient = ((eye(3) - axis * axis.T) * cos(theta) +
                    Matrix([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]]) * sin(theta) + axis * axis.T)
        elif rot_type == 'QUATERNION':
            if not rot_order == '':
                raise TypeError(
                    'Quaternion orientation takes no rotation order')
            if not (isinstance(amounts, (list, tuple)) & (len(amounts) == 4)):
                raise TypeError('Amounts are a list or tuple of length 4')
            q0, q1, q2, q3 = amounts
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
        else:
            raise NotImplementedError('That is not an implemented rotation')
        self._dcm_dict.update({parent: parent_orient})
        parent._dcm_dict.update({self: parent_orient.T})
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
        else:
            try:
                from sympy.polys.polyerrors import CoercionFailed
                from sympy.physics.mechanics.functions import kinematic_equations
                q1, q2, q3 = amounts
                u1, u2, u3 = dynamicsymbols('u1, u2, u3')
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

    def orientnew(self, newname, rot_type, amounts, rot_order='', indices=None,
            latexs=None):
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])


        .orient() documentation:\n
        ========================

        """

        newframe = ReferenceFrame(newname, indices, latexs)
        newframe.orient(self, rot_type, amounts, rot_order)
        return newframe

    orientnew.__doc__ += orient.__doc__

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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_acc(N, V)
        >>> A.ang_acc_in(N)
        10*N.x

        """

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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> V = 10 * N.x
        >>> A.set_ang_vel(N, V)
        >>> A.ang_vel_in(N)
        10*N.x

        """

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


class Vector(object):
    """The class used to define vectors.

    It along with ReferenceFrame are the building blocks of describing a
    classical mechanics system in PyDy.

    Attributes
    ==========

    simp : Boolean
        Let certain methods use trigsimp on their outputs

    """

    simp = False

    def __init__(self, inlist):
        """This is the constructor for the Vector class.  You shouldn't be
        calling this, it should only be used by other functions. You should be
        treating Vectors like you would with if you were doing the math by
        hand, and getting the first 3 from the standard basis vectors from a
        ReferenceFrame.

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

    def __hash__(self):
        return hash(tuple(self.args))

    def __add__(self, other):
        """The add operator for Vector. """
        other = _check_vector(other)
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector, dot
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> dot(N.x, N.x)
        1
        >>> dot(N.x, N.y)
        0
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])
        >>> dot(N.y, A.y)
        cos(q1)

        """

        if isinstance(other, Dyadic):
            return NotImplemented
        other = _check_vector(other)
        out = S(0)
        for i, v1 in enumerate(self.args):
            for j, v2 in enumerate(other.args):
                out += ((v2[0].T)
                        * (v2[1].dcm(v1[1]))
                        * (v1[0]))[0]
        if Vector.simp is True:
            return trigsimp(sympify(out), recursive=True)
        else:
            return sympify(out)

    def __div__(self, other):
        """This uses mul and inputs self and 1 divided by other. """
        return self.__mul__(1 / other)

    __truediv__ = __div__

    def __eq__(self, other):
        """Tests for equality.

        It is very import to note that this is only as good as the SymPy
        equality test; False does not always mean they are not equivalent
        Vectors.
        If other is 0, and self is empty, returns True.
        If other is 0 and self is not empty, returns False.
        If none of the above, only accepts other as a Vector.

        """

        other = _check_vector(other)
        if (self.args == []) and (other.args == []):
            return True
        elif (self.args == []) or (other.args == []):
            return False

        frame = self.args[0][1]
        for v in frame:
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> from sympy import Symbol
        >>> N = ReferenceFrame('N')
        >>> b = Symbol('b')
        >>> V = 10 * b * N.x
        >>> print(V)
        10*b*N.x

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

        A rank increasing operation, which returns a Dyadic from two Vectors

        Parameters
        ==========

        other : Vector
            The Vector to take the outer product with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> outer(N.x, N.x)
        (N.x|N.x)

        """

        other = _check_vector(other)
        ol = Dyadic([])
        for i, v in enumerate(self.args):
            for i2, v2 in enumerate(other.args):
                # it looks this way because if we are in the same frame and
                # use the enumerate function on the same frame in a nested
                # fashion, then bad things happen
                ol += Dyadic([(v[0][0] * v2[0][0], v[1].x, v2[1].x)])
                ol += Dyadic([(v[0][0] * v2[0][1], v[1].x, v2[1].y)])
                ol += Dyadic([(v[0][0] * v2[0][2], v[1].x, v2[1].z)])
                ol += Dyadic([(v[0][1] * v2[0][0], v[1].y, v2[1].x)])
                ol += Dyadic([(v[0][1] * v2[0][1], v[1].y, v2[1].y)])
                ol += Dyadic([(v[0][1] * v2[0][2], v[1].y, v2[1].z)])
                ol += Dyadic([(v[0][2] * v2[0][0], v[1].z, v2[1].x)])
                ol += Dyadic([(v[0][2] * v2[0][1], v[1].z, v2[1].y)])
                ol += Dyadic([(v[0][2] * v2[0][2], v[1].z, v2[1].z)])
        return ol

    def _latex(self, printer=None):
        """Latex Printing method. """
        ar = self.args  # just to shorten things
        if len(ar) == 0:
            return str(0)
        ol = []  # output list, to be concatenated to a string
        for i, v in enumerate(ar):
            for j in 0, 1, 2:
                # if the coef of the basis vector is 1, we skip the 1
                if ar[i][0][j] == 1:
                    ol.append(' + ' + ar[i][1].latex_vecs[j])
                # if the coef of the basis vector is -1, we skip the 1
                elif ar[i][0][j] == -1:
                    ol.append(' - ' + ar[i][1].latex_vecs[j])
                elif ar[i][0][j] != 0:
                    # If the coefficient of the basis vector is not 1 or -1;
                    # also, we might wrap it in parentheses, for readability.
                    arg_str = MechanicsStrPrinter().doprint(ar[i][0][j])
                    if isinstance(ar[i][0][j], Add):
                        arg_str = "(%s)" % arg_str
                    if arg_str[0] == '-':
                        arg_str = arg_str[1:]
                        str_start = ' - '
                    else:
                        str_start = ' + '
                    ol.append(str_start + arg_str + '*' +
                              ar[i][1].latex_vecs[j])
        outstr = ''.join(ol)
        if outstr.startswith(' + '):
            outstr = outstr[3:]
        elif outstr.startswith(' '):
            outstr = outstr[1:]
        return outstr

    def _pretty(self, printer=None):
        """Pretty Printing method. """
        e = self

        class Fake(object):
            baseline = 0

            def render(self, *args, **kwargs):
                self = e
                ar = self.args  # just to shorten things
                if len(ar) == 0:
                    return unicode(0)
                ol = []  # output list, to be concatenated to a string
                for i, v in enumerate(ar):
                    for j in 0, 1, 2:
                        # if the coef of the basis vector is 1, we skip the 1
                        if ar[i][0][j] == 1:
                            ol.append(u(" + ") + ar[i][1].pretty_vecs[j])
                        # if the coef of the basis vector is -1, we skip the 1
                        elif ar[i][0][j] == -1:
                            ol.append(u(" - ") + ar[i][1].pretty_vecs[j])
                        elif ar[i][0][j] != 0:
                            # If the basis vector coeff is not 1 or -1,
                            # we might wrap it in parentheses, for readability.
                            arg_str = (MechanicsPrettyPrinter().doprint(
                                ar[i][0][j]))
                            if isinstance(ar[i][0][j], Add):
                                arg_str = u("(%s)") % arg_str
                            if arg_str[0] == u("-"):
                                arg_str = arg_str[1:]
                                str_start = u(" - ")
                            else:
                                str_start = u(" + ")
                            ol.append(str_start + arg_str + '*' +
                                      ar[i][1].pretty_vecs[j])
                outstr = u("").join(ol)
                if outstr.startswith(u(" + ")):
                    outstr = outstr[3:]
                elif outstr.startswith(" "):
                    outstr = outstr[1:]
                return outstr
        return Fake()

    def __ror__(self, other):
        """Outer product between two Vectors.

        A rank increasing operation, which returns a Dyadic from two Vectors

        Parameters
        ==========

        other : Vector
            The Vector to take the outer product with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> outer(N.x, N.x)
        (N.x|N.x)

        """

        other = _check_vector(other)
        ol = Dyadic([])
        for i, v in enumerate(other.args):
            for i2, v2 in enumerate(self.args):
                # it looks this way because if we are in the same frame and
                # use the enumerate function on the same frame in a nested
                # fashion, then bad things happen
                ol += Dyadic([(v[0][0] * v2[0][0], v[1].x, v2[1].x)])
                ol += Dyadic([(v[0][0] * v2[0][1], v[1].x, v2[1].y)])
                ol += Dyadic([(v[0][0] * v2[0][2], v[1].x, v2[1].z)])
                ol += Dyadic([(v[0][1] * v2[0][0], v[1].y, v2[1].x)])
                ol += Dyadic([(v[0][1] * v2[0][1], v[1].y, v2[1].y)])
                ol += Dyadic([(v[0][1] * v2[0][2], v[1].y, v2[1].z)])
                ol += Dyadic([(v[0][2] * v2[0][0], v[1].z, v2[1].x)])
                ol += Dyadic([(v[0][2] * v2[0][1], v[1].z, v2[1].y)])
                ol += Dyadic([(v[0][2] * v2[0][2], v[1].z, v2[1].z)])
        return ol

    def __rsub__(self, other):
        return (-1 * self) + other

    def __str__(self, printer=None):
        """Printing method. """
        ar = self.args  # just to shorten things
        if len(ar) == 0:
            return str(0)
        ol = []  # output list, to be concatenated to a string
        for i, v in enumerate(ar):
            for j in 0, 1, 2:
                # if the coef of the basis vector is 1, we skip the 1
                if ar[i][0][j] == 1:
                    ol.append(' + ' + ar[i][1].str_vecs[j])
                # if the coef of the basis vector is -1, we skip the 1
                elif ar[i][0][j] == -1:
                    ol.append(' - ' + ar[i][1].str_vecs[j])
                elif ar[i][0][j] != 0:
                    # If the coefficient of the basis vector is not 1 or -1;
                    # also, we might wrap it in parentheses, for readability.
                    arg_str = MechanicsStrPrinter().doprint(ar[i][0][j])
                    if isinstance(ar[i][0][j], Add):
                        arg_str = "(%s)" % arg_str
                    if arg_str[0] == '-':
                        arg_str = arg_str[1:]
                        str_start = ' - '
                    else:
                        str_start = ' + '
                    ol.append(str_start + arg_str + '*' + ar[i][1].str_vecs[j])
        outstr = ''.join(ol)
        if outstr.startswith(' + '):
            outstr = outstr[3:]
        elif outstr.startswith(' '):
            outstr = outstr[1:]
        return outstr

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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector
        >>> from sympy import symbols
        >>> q1 = symbols('q1')
        >>> N = ReferenceFrame('N')
        >>> N.x ^ N.y
        N.z
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])
        >>> A.x ^ N.y
        N.z
        >>> N.y ^ A.x
        - sin(q1)*A.y - cos(q1)*A.z

        """

        if isinstance(other, Dyadic):
            return NotImplemented
        other = _check_vector(other)
        if other.args == []:
            return self * S(0)

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
        ar = other.args  # For brevity
        for i, v in enumerate(ar):
            tempx = v[1].x
            tempy = v[1].y
            tempz = v[1].z
            tempm = ([[tempx, tempy, tempz], [self & tempx, self & tempy,
                self & tempz], [Vector([ar[i]]) & tempx,
                Vector([ar[i]]) & tempy, Vector([ar[i]]) & tempz]])
            outvec += _det(tempm)
        return outvec

    _sympystr = __str__
    _sympyrepr = _sympystr
    __repr__ = __str__
    __radd__ = __add__
    __rand__ = __and__
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector, dynamicsymbols
        >>> from sympy import Symbol
        >>> Vector.simp = True
        >>> t = Symbol('t')
        >>> q1 = dynamicsymbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.y])
        >>> A.x.diff(t, N)
        - q1'*A.z

        """

        wrt = sympify(wrt)
        _check_frame(otherframe)
        outvec = S(0)
        for i, v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(wrt), otherframe)])
            else:
                if otherframe.dcm(v[1]).diff(wrt) == zeros(3, 3):
                    d = v[0].diff(wrt)
                    outvec += Vector([(d, v[1])])
                else:
                    d = (Vector([v]).express(otherframe)).args[0][0].diff(wrt)
                    outvec += Vector([(d, otherframe)]).express(v[1])
        return outvec

    def doit(self, **hints):
        """Calls .doit() on each term in the Vector"""
        ov = S(0)
        for i, v in enumerate(self.args):
            ov += Vector([(v[0].applyfunc(lambda x: x.doit(**hints)), v[1])])
        return ov

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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector, dynamicsymbols
        >>> from sympy import Symbol
        >>> q1 = Symbol('q1')
        >>> u1 = dynamicsymbols('u1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.x])
        >>> v = u1 * N.x
        >>> A.set_ang_vel(N, 10*A.x)
        >>> A.x.dt(N) == 0
        True
        >>> v.dt(N)
        u1'*N.x

        """

        outvec = S(0)
        _check_frame(otherframe)
        for i, v in enumerate(self.args):
            if v[1] == otherframe:
                outvec += Vector([(v[0].diff(dynamicsymbols._t), otherframe)])
            else:
                outvec += (Vector([v]).dt(v[1]) +
                    (v[1].ang_vel_in(otherframe) ^ Vector([v])))
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

        >>> from sympy.physics.mechanics import ReferenceFrame, Vector, dynamicsymbols
        >>> q1 = dynamicsymbols('q1')
        >>> N = ReferenceFrame('N')
        >>> A = N.orientnew('A', 'Axis', [q1, N.y])
        >>> A.x.express(N)
        cos(q1)*N.x - sin(q1)*N.z

        """

        _check_frame(otherframe)
        outvec = Vector([])
        for i, v in enumerate(self.args):
            if v[1] != otherframe:
                temp = otherframe.dcm(v[1]) * v[0]
                if Vector.simp is True:
                    temp = temp.applyfunc(lambda x: trigsimp(x, method='fu'))
                outvec += Vector([(temp, otherframe)])
            else:
                outvec += Vector([v])
        return outvec

    def simplify(self):
        """Returns a simplified Vector."""
        outvec = Vector([])
        for i in self.args:
            outvec += Vector([(i[0].simplify(), i[1])])
        return outvec

    def subs(self, *args, **kwargs):
        """Substituion on the Vector.

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame
        >>> from sympy import Symbol
        >>> N = ReferenceFrame('N')
        >>> s = Symbol('s')
        >>> a = N.x * s
        >>> a.subs({s: 2})
        2*N.x

        """

        ov = S(0)
        for i, v in enumerate(self.args):
            ov += Vector([(v[0].subs(*args, **kwargs), v[1])])
        return ov

    def magnitude(self):
        """Returns the magnitude (Euclidean norm) of self."""
        return sqrt(self & self)

    def normalize(self):
        """Returns a Vector of magnitude 1, codirectional with self."""
        return Vector(self.args + []) / self.magnitude()


class MechanicsStrPrinter(StrPrinter):
    """String Printer for mechanics. """

    def _print_Derivative(self, e):
        t = dynamicsymbols._t
        if (bool(sum([i == t for i in e.variables])) &
                isinstance(type(e.args[0]), UndefinedFunction)):
            ol = str(e.args[0].func)
            for i, v in enumerate(e.variables):
                ol += dynamicsymbols._str
            return ol
        else:
            return StrPrinter().doprint(e)

    def _print_Function(self, e):
        t = dynamicsymbols._t
        if isinstance(type(e), UndefinedFunction):
            return StrPrinter().doprint(e).replace("(%s)" % t, '')
        return e.func.__name__ + "(%s)" % self.stringify(e.args, ", ")


class MechanicsLatexPrinter(LatexPrinter):
    """Latex Printer for mechanics. """

    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__
        t = dynamicsymbols._t

        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        elif isinstance(type(expr), UndefinedFunction) and (expr.args == (t,)):
            name, sup, sub = split_super_sub(func)
            if len(sup) != 0:
                sup = r"^{%s}" % "".join(sup)
            else:
                sup = r""
            if len(sub) != 0:
                sub = r"_{%s}" % "".join(sub)
            else:
                sub = r""
            return r"%s" % (name + sup + sub)
        else:
            args = [ str(self._print(arg)) for arg in expr.args ]
            # How inverse trig functions should be displayed, formats are:
            # abbreviated: asin, full: arcsin, power: sin^-1
            inv_trig_style = self._settings['inv_trig_style']
            # If we are dealing with a power-style inverse trig function
            inv_trig_power_case = False
            # If it is applicable to fold the argument brackets
            can_fold_brackets = self._settings['fold_func_brackets'] and \
                len(args) == 1 and \
                not self._needs_function_brackets(expr.args[0])

            inv_trig_table = ["asin", "acos", "atan", "acot"]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    func = func
                elif inv_trig_style == "full":
                    func = "arc" + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                name = r"\operatorname{%s}^{%s}" % (func, exp)
            else:
                name = r"\operatorname{%s}" % func

            if can_fold_brackets:
                name += r"%s"
            else:
                name += r"\left(%s\right)"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            return name % ",".join(args)

    def _print_Derivative(self, der_expr):
        # make sure it is an the right form
        der_expr = der_expr.doit()
        if not isinstance(der_expr, Derivative):
            return self.doprint(der_expr)

        # check if expr is a dynamicsymbol
        from sympy.core.function import AppliedUndef
        t = dynamicsymbols._t
        expr = der_expr.expr
        red = expr.atoms(AppliedUndef)
        syms = der_expr.variables
        test1 = not all([True for i in red if i.atoms() == set([t])])
        test2 = not all([(t == i) for i in syms])
        if test1 or test2:
            return LatexPrinter().doprint(der_expr)

        # done checking
        dots = len(syms)
        base = self._print_Function(expr)
        base_split = base.split('_', 1)
        base = base_split[0]
        if dots == 1:
            base = r"\dot{%s}" % base
        elif dots == 2:
            base = r"\ddot{%s}" % base
        elif dots == 3:
            base = r"\dddot{%s}" % base
        if len(base_split) is not 1:
            base += '_' + base_split[1]
        return base


class MechanicsPrettyPrinter(PrettyPrinter):
    """Pretty Printer for mechanics. """

    def _print_Derivative(self, deriv):
        # XXX use U('PARTIAL DIFFERENTIAL') here ?
        t = dynamicsymbols._t
        dots = 0
        can_break = True
        syms = list(reversed(deriv.variables))
        x = None

        while len(syms) > 0:
            if syms[-1] == t:
                syms.pop()
                dots += 1
            else:
                break

        f = prettyForm(binding=prettyForm.FUNC, *self._print(deriv.expr))
        if not (isinstance(type(deriv.expr), UndefinedFunction)
                and (deriv.expr.args == (t,))):
            dots = 0
            can_break = False
            f = prettyForm(binding=prettyForm.FUNC,
                    *self._print(deriv.expr).parens())

        if dots == 0:
            dots = u("")
        elif dots == 1:
            dots = u("\u0307")
        elif dots == 2:
            dots = u("\u0308")
        elif dots == 3:
            dots = u("\u20db")
        elif dots == 4:
            dots = u("\u20dc")

        uni_subs = [u("\u2080"), u("\u2081"), u("\u2082"), u("\u2083"), u("\u2084"),
                    u("\u2085"), u("\u2086"), u("\u2087"), u("\u2088"), u("\u2089"),
                    u("\u208a"), u("\u208b"), u("\u208c"), u("\u208d"), u("\u208e"),
                    u("\u208f"), u("\u2090"), u("\u2091"), u("\u2092"), u("\u2093"),
                    u("\u2094"), u("\u2095"), u("\u2096"), u("\u2097"), u("\u2098"),
                    u("\u2099"), u("\u209a"), u("\u209b"), u("\u209c"), u("\u209d"),
                    u("\u209e"), u("\u209f")]

        fpic = f.__dict__['picture']
        funi = f.__dict__['unicode']
        ind = len(funi)
        val = ""

        for i in uni_subs:
            cur_ind = funi.find(i)
            if (cur_ind != -1) and (cur_ind < ind):
                ind = cur_ind
                val = i
        if ind == len(funi):
            funi += dots
        else:
            funi = funi.replace(val, dots + val)

        if f.__dict__['picture'] == [f.__dict__['unicode']]:
            fpic = [funi]
        f.__dict__['picture'] = fpic
        f.__dict__['unicode'] = funi

        if (len(syms)) == 0 and can_break:
            return f

        for sym, num in group(syms, multiple=False):
            s = self._print(sym)
            ds = prettyForm(*s.left('d'))

            if num > 1:
                ds = ds**prettyForm(str(num))

            if x is None:
                x = ds
            else:
                x = prettyForm(*x.right(' '))
                x = prettyForm(*x.right(ds))
        pform = prettyForm('d')
        if len(syms) > 1:
            pform = pform**prettyForm(str(len(syms)))
        pform = prettyForm(*pform.below(stringPict.LINE, x))
        pform.baseline = pform.baseline + 1
        pform = prettyForm(*stringPict.next(pform, f))
        return pform

    def _print_Function(self, e):
        t = dynamicsymbols._t
        # XXX works only for applied functions
        func = e.func
        args = e.args
        func_name = func.__name__
        prettyFunc = self._print(C.Symbol(func_name))
        prettyArgs = prettyForm(*self._print_seq(args).parens())
        # If this function is an Undefined function of t, it is probably a
        # dynamic symbol, so we'll skip the (t). The rest of the code is
        # identical to the normal PrettyPrinter code
        if isinstance(func, UndefinedFunction) and (args == (t,)):
            pform = prettyForm(binding=prettyForm.FUNC,
                       *stringPict.next(prettyFunc))
        else:
            pform = prettyForm(binding=prettyForm.FUNC,
                       *stringPict.next(prettyFunc, prettyArgs))
        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs
        return pform


def _check_dyadic(other):
    if not isinstance(other, Dyadic):
        other = sympify(other)
        if other != S(0):
            raise TypeError('A Dyadic must be supplied')
        else:
            other = Dyadic([])
    return other


def _check_frame(other):
    if not isinstance(other, ReferenceFrame):
        raise TypeError('A ReferenceFrame must be supplied')


def _check_vector(other):
    if not isinstance(other, Vector):
        other = sympify(other)
        if other != S(0):
            raise TypeError('A Vector must be supplied')
        else:
            other = Vector([])
    return other


def dynamicsymbols(names, level=0):
    """Uses symbols and Function for functions of time.

    Creates a SymPy UndefinedFunction, which is then initialized as a function
    of a variable, the default being Symbol('t').

    Parameters
    ==========

    names : str
        Names of the dynamic symbols you want to create; works the same way as
        inputs to symbols
    level : int
        Level of differentiation of the returned function; d/dt once of t,
        twice of t, etc.

    Examples
    =======

    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> from sympy import diff, Symbol
    >>> q1 = dynamicsymbols('q1')
    >>> q1
    q1(t)
    >>> diff(q1, Symbol('t'))
    Derivative(q1(t), t)

    """

    esses = symbols(names, cls=Function)
    t = dynamicsymbols._t
    if hasattr(esses, '__iter__'):
        esses = [reduce(diff, [t]*level, e(t)) for e in esses]
        return esses
    else:
        return reduce(diff, [t]*level, esses(t))

dynamicsymbols._t = Symbol('t')
dynamicsymbols._str = '\''
