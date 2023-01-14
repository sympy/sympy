import warnings
from sympy.utilities.iterables import iterable
from sympy.matrices.dense import MutableDenseMatrix as Matrix
from sympy.core.numbers import pi
from sympy.algebras.quaternion import Quaternion
from sympy.core.sympify import sympify, _sympify
from sympy.core.expr import Expr

from mpmath.libmp.libmpf import prec_to_dps


def _check_sequence(seq):
    """validate seq and return info"""
    if type(seq) != str:
        raise ValueError('Expected seq to be a string.')
    if len(seq) != 3:
        raise ValueError("Expected 3 axes, got `{}`.".format(seq))

    intrinsic = seq.isupper()
    extrinsic = seq.islower()
    if not (intrinsic or extrinsic):
        raise ValueError("seq must either be fully uppercase (for extrinsic "
                         "rotations), or fully lowercase, for intrinsic "
                         "rotations).")

    i, j, k = seq.lower()
    if (i == j) or (j == k):
        raise ValueError("Consecutive axes must be different")

    bad = set(seq) - set('xyzXYZ')
    if bad:
        raise ValueError("Expected axes from `seq` to be from "
                         "['x', 'y', 'z'] or ['X', 'Y', 'Z'], "
                         "got {}".format(''.join(bad)))

    # get elementary basis vectors
    ei = [1 if n == i else 0 for n in 'xyz']
    ej = [1 if n == j else 0 for n in 'xyz']
    ek = [1 if n == k else 0 for n in 'xyz']

    # get indices
    i = 'xyz'.index(i) + 1
    j = 'xyz'.index(j) + 1
    k = 'xyz'.index(k) + 1

    info = {
            'seq': seq,
            'extrinsic': extrinsic,
            'intrinsic': intrinsic,
            'symmetric': i == k,
            'i': i,
            'j': j,
            'k': k,
            'ei': ei,
            'ej': ej,
            'ek': ek
            }

    return info


def _wrap(angle):
    if angle.is_number:
        return (angle + pi) % (2 * pi) - pi
    else:
        return angle


class Euler(Expr):
    """Provides Euler angles operations.

    Euler objects can be instantiated as Euler(alpha, beta, gamma, seq).

    Parameters
    ==========

    seq : string of length 3
        Represents the sequence of rotations.
        For intrinsic rotations, seq must be all lowercase and its elements
        must be from the set `{'x', 'y', 'z'}`
        For extrinsic rotations, seq must be all uppercase and its elements
        must be from the set `{'X', 'Y', 'Z'}`

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Euler_angles

    """
    _op_priority = 11.0

    is_commutative = False

    def __new__(cls, alpha=0, beta=0, gamma=0, seq='ZYX'):
        info = _check_sequence(seq)
        alpha, beta, gamma = map(sympify, (alpha, beta, gamma))
        alpha, beta, gamma = map(_wrap, (alpha, beta, gamma))

        if any(i.is_commutative is False for i in [alpha, beta, gamma]):
            raise ValueError("arguments have to be commutative")
        else:
            obj = Expr.__new__(cls, alpha, beta, gamma)
            obj._alpha = alpha
            obj._beta = beta
            obj._gamma = gamma
            obj._info = info
            obj._q = None
            obj._R = None
            return obj

    @property
    def alpha(self):
        return self._alpha

    @property
    def beta(self):
        return self._beta

    @property
    def gamma(self):
        return self._gamma

    @property
    def angles(self):
        return self.args

    @property
    def info(self):
        return self._info

    @property
    def seq(self):
        return self._info['seq']

    @property
    def extrinsic(self):
        return self._info['extrinsic']

    @property
    def intrinsic(self):
        return self._info['intrinsic']

    @property
    def symmetric(self):
        return self._info['symmetric']

    @property
    def i(self):
        return self._info['i']

    @property
    def j(self):
        return self._info['j']

    @property
    def k(self):
        return self._info['k']

    @property
    def ei(self):
        return self._info['ei']

    @property
    def ej(self):
        return self._info['ej']

    @property
    def ek(self):
        return self._info['ek']

    def __getitem__(self, key):
        return self.args[key]

    def to_Matrix(self):
        return Matrix(self.args)

    def to_quaternion(self):
        """Returns the equivalent full rotation quaternion equivalent to Euler
        angles.

        Returns
        =======

        Quaternion
            The equivalent rotation quaternion

        Examples
        ========

        >>> from sympy import Euler
        >>> from sympy import pi
        >>> angles = Euler(pi/2, 0, 0, seq='ZYX')
        >>> angles
        Euler(pi/2,  0,  0, ZYX)

        >>> angles.to_rotation_matrix()
        Matrix([
        [0, -1, 0],
        [1,  0, 0],
        [0,  0, 1]])

        Generates a 4x4 transformation matrix (used for rotation about a point
        other than the origin) if the point(v) is passed as an argument.

        Examples
        ========

        >>> from sympy import Euler
        >>> from sympy import pi
        >>> angles = Euler(pi/2, 0, 0, seq='ZYX')
        >>> angles.to_quaternion()
        sqrt(2)/2 + 0*i + 0*j + sqrt(2)/2*k
        """
        if self._q is None:
            self._q = Quaternion.from_euler(self.args, self.seq)
        return self._q

    @classmethod
    def from_quaternion(self, q, seq,
                        angle_addition=True, avoid_square_root=False):
        r"""Returns Euler angles representing same rotation as the quaternion,
        in the sequence given by `seq`. This implements the method described
        in [1]_.

        Parameters
        ==========

        q : Quaternion or list
            Rotation quaternion elements.

        seq : string of length 3
            Represents the sequence of rotations.
            For intrinsic rotations, seq must be all lowercase and its elements
            must be from the set `{'x', 'y', 'z'}`
            For extrinsic rotations, seq must be all uppercase and its elements
            must be from the set `{'X', 'Y', 'Z'}`

        angle_addition : bool
            Default : True
            When True, first and third angles are given as an addition and
            subtraction of two simpler `atan2` expressions. When False, the
            first and third angles are each given by a single more complicated
            `atan2` expression. This equivalent is given by:

            --math::

                \operatorname{atan_2} (b,a) \pm \operatorname{atan_2} (d,c) =
                \operatorname{atan_2} (bc\pm ad, ac\mp bd)

        avoid_square_root : bool
            Default : False
            When True, the second angle is calculated with an expression based
            on acos`, which is slightly more complicated but avoids a square
            root. When False, second angle is calculated with `atan2`, which
            is simpler and can be better for numerical reasons (some
            numerical implementations of `acos` have problems near zero).

        Returns
        =======

        Tuple
            The Euler angles calculated from the quaternion

        Examples
        ========

        >>> from sympy import Quaternion
        >>> from sympy.abc import a, b, c, d
        >>> euler = Quaternion(a, b, c, d).to_euler('zyz')
        >>> euler
        (-atan2(-b, c) + atan2(d, a),
         2*atan2(sqrt(b**2 + c**2), sqrt(a**2 + d**2)),
         atan2(-b, c) + atan2(d, a))

        References
        ==========

        .. [1] https://doi.org/10.1371/journal.pone.0276302

        """
        if not isinstance(q, Quaternion):
            if not iterable(q):
                q = [q]
            q = Quaternion(*q)

        angles = q.to_euler(seq,
                            angle_addition=angle_addition,
                            avoid_square_root=avoid_square_root)
        return Euler(*angles, seq=seq)

    def add(self, other):
        """Adds Euler angles elements element-wise, treating it as a Matrix.

        This does not have much mathematical meaning, are you sure this is what
        you want?

        Parameters
        ==========

        other : Euler
            The Euler angles to add.

        Returns
        =======

        Quaternion
            The Euler angles resulting from the addition.

        Examples
        ========

        >>> from sympy import Euler, pi
        >>> from sympy.abc import alpha, beta, gamma
        >>> euler = Euler(alpha, beta, gamma)
        >>> euler + Euler(pi, 0, pi)
        Euler((alpha - pi),  beta,  (gamma - pi), ZYX)

        """
        warnings.warn('Adding Euler as if it were a Matrix, '
                      'is this what you want?')

        if isinstance(other, Euler):
            if other.seq != self.seq:
                warnings.warn('Adding Euler elements of different sequence, '
                              'assuming first sequence.')
            other = other.to_Matrix()

        euler = self.to_Matrix() + other
        return Euler(*euler, seq=self.seq)

    def to_rotation_matrix(self, v=None, homogeneous=True):
        """Returns the equivalent full rotation transformation matrix of Euler
        angles.

        Parameters
        ==========

        v : tuple or None
            Default value: None
        homogeneous : bool
            When True, gives an expression that may be more efficient for
            symbolic calculations but less so for direct evaluation. Both
            formulas are mathematically equivalent.
            Default value: True

        Returns
        =======

        tuple
            Returns the equivalent rotation transformation matrix of the angles
            which represents rotation about the origin if v is not passed.

        Examples
        ========

        >>> from sympy import Euler
        >>> from sympy import pi
        >>> angles = Euler(pi/2, 0, 0, seq='ZYX')
        >>> angles
        Euler(pi/2,  0,  0, ZYX)

        >>> angles.to_rotation_matrix()
        Matrix([
        [0, -1, 0],
        [1,  0, 0],
        [0,  0, 1]])

        Generates a 4x4 transformation matrix (used for rotation about a point
        other than the origin) if the point(v) is passed as an argument.

        Examples
        ========

        >>> from sympy import Euler
        >>> from sympy import pi
        >>> angles = Euler(pi/2, 0, 0, seq='ZYX')
        >>> angles.to_rotation_matrix((1, 1, 1))
        Matrix([
        [0, -1, 0, 2],
        [1,  0, 0, 0],
        [0,  0, 1, 0],
        [0,  0, 0, 1]])

        """
        if self._R is None:
            self._R = self.to_quaternion().to_rotation_matrix(v, homogeneous)
        return self._R

    @staticmethod
    def _generic_mul(lhs, rhs):
        """Generic multiplication.

        Parameters
        ==========

        lhs : Euler or symbol
        rhs : Euler or symbol

        It is important to note that if neither lhs nor rhs is an Euler
        instance, this function simply returns lhs * rhs.

        Returns
        =======

        Euler
            The resultant Euler after multiplying lhs and rhs.

        Examples
        ========

        Multiplication with a number multiplies every angle by this number:

        >>> from sympy import Euler, pi
        >>> euler = Euler(pi / 4, 0, 0, 'ZYX')
        >>> Euler._generic_mul(2, euler)
        Euler(pi/2,  0,  0, ZYX)

        Wrapping always happens, and the angles are between [-pi and pi],
        so the result might seem unexpected:

        >>> Euler._generic_mul(4, euler)
        Euler((-pi),  0,  0, ZYX)

        Possible to compose Euler instances. For simple cases it might work
        well:
        >>> euler = Euler((-pi),  0,  0, 'ZYX')
        >>> Euler._generic_mul(euler, euler)
        Euler(0,  0,  0, ZYX)

        Most of the time, it is not very efficient:

        >>> euler = Euler(pi / 4, 0, 0, 'ZYX')
        >>> Euler.mul(euler, euler)
        Euler(2*atan(2*sqrt(2)*sqrt(1/2 - sqrt(2)/4)*sqrt(sqrt(2)/4 + 1/2)),  0,  0, ZYX)

        By evaluating it, we can see it is equivalent to the expected answer:

        >>> Euler.mul(euler, euler).evalf()
        Euler(1.5707963267949,  0,  0, ZYX)
        """
        lhs = sympify(lhs)
        rhs = sympify(rhs)

        # None is a Euler:
        if not isinstance(lhs, Euler) and not isinstance(rhs, Euler):
            return lhs * rhs

        # If lhs is a number or a SymPy expression instead of a Euler
        if not isinstance(lhs, Euler):
            if lhs.is_commutative:
                return Euler(*(lhs * i for i in rhs.args), seq=rhs.seq)
            else:
                raise ValueError('Only commutative expressions can be '
                                 'multiplied with an Euler.')

        # If rhs is a number or a SymPy expression instead of a Euler
        if not isinstance(rhs, Euler):
            if rhs.is_commutative:
                return Euler(*(rhs * i for i in lhs.args), seq=lhs.seq)
            else:
                raise ValueError('Only commutative expressions can be '
                                 'multiplied with an Euler.')

        if lhs.seq != rhs.seq:
            warnings.warn('Different sequences detected, '
                          'returning sequence of lhs.')
        else:
            warnings.warn('Euler angles are bad for composing, '
                          'consider using quaternions instead.')

        q = lhs.to_quaternion() * rhs.to_quaternion()
        return Euler.from_quaternion(q, lhs.seq)

    def set_sequence(self, new_sequence):
        """Gets the equivalent Euler angles in the new sequence."""
        if self.seq == new_sequence:
            return self

        q = self.to_quaternion()
        return Euler.from_quaternion(q, new_sequence)

    def switch_set(self):
        """Gets a different set of Euler angles int the same sequence
        that correspond to the same rotation."""
        lamb = 0 if self.symmetric else pi
        alpha = self.alpha + pi
        beta = -self.beta + lamb
        gamma = self.gamma + pi
        return Euler(alpha, beta, gamma, self.seq)

    def mul(self, other):
        """Multiplies Euler angles.

        Euler angles are NOT suited for angle composition.
        When multiplying two Euler angles elements, both are converted to
        quaternions, multiplied, then converted back to Euler angles.

        Parameters
        ==========

        other : Euler or symbol
            The Euler to multiply to current (self) Euler.

        Returns
        =======

        Euler
            The resultant Euler angles after multiplying self with other

        Examples
        ========

        Multiplication with a number multiplies every angle by this number:

        >>> from sympy import Euler, pi
        >>> euler = Euler(pi / 4, 0, 0, 'ZYX')
        >>> euler.mul(2)
        Euler(pi/2,  0,  0, ZYX)

        Wrapping always happens, and the angles are between [-pi and pi],
        so the result might seem unexpected:

        >>> euler.mul(4)
        Euler((-pi),  0,  0, ZYX)

        Possible to compose Euler instances. For simple cases it might work
        well:
        >>> euler = Euler((-pi),  0,  0, 'ZYX')
        >>> euler.mul(euler)
        Euler(0,  0,  0, ZYX)

        Most of the time, it is not very efficient:

        >>> euler = Euler(pi / 4, 0, 0, 'ZYX')
        >>> euler.mul(euler)
        Euler(2*atan(2*sqrt(2)*sqrt(1/2 - sqrt(2)/4)*sqrt(sqrt(2)/4 + 1/2)),  0,  0, ZYX)

        By evaluating it, we can see it is equivalent to the expected answer:

        >>> euler.mul(euler).evalf()
        Euler(1.5707963267949,  0,  0, ZYX)
        """
        return self._generic_mul(self, _sympify(other))

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.add(-other)

    def __mul__(self, other):
        return self._generic_mul(self, _sympify(other))

    def __rmul__(self, other):
        return self._generic_mul(_sympify(other), self)

    def __pow__(self, p):
        return self.pow(p)

    def __neg__(self):
        return Euler(-self.alpha, -self.beta, -self.gamma, seq=self.seq)

    def __truediv__(self, other):
        return self * sympify(other)**-1

    def __rtruediv__(self, other):
        raise ValueError('Division by an Euler is not well defined.')

    def inverse(self):
        """Returns inverse Euler angles

        In this implementation, the sequence is also reversed.
        """
        seq = self.seq[::-1]
        alpha, beta, gamma = [-i for i in self.args[::-1]]
        return Euler(alpha, beta, gamma, seq)

    def pow(self, p):
        """Finds the pth power of the Euler angles.

        Parameters
        ==========

        p : int
            Power to be applied on Euler.

        Returns
        =======

        Euler
            Returns the p-th power of the current Euler.
            Returns the inverse if p = -1.

        Examples
        ========

        >>> from sympy import Euler, pi
        >>> euler = Euler(0, 0, pi / 2, 'ZYX')
        >>> euler.pow(2)
        Euler(0,  0,  (-pi), ZYX)

        For negative powers, returns the inverse:

        >>> euler.pow(-1)
        Euler((-pi/2),  0,  0, XYZ)

        """
        p = sympify(p)
        if p == -1:
            return self.inverse()

        q = self.to_quaternion() ** p
        return Euler.from_quaternion(q, seq=self.seq)

    @property
    def is_zero_rotation(self):
        """Returns true if the angles are all zero"""
        return all(i.is_zero for i in self.args)

    def _eval_evalf(self, prec):
        """Returns the floating point approximations (decimal numbers) of the
        angles.

        Returns
        =======
        Euler
            Floating point approximations of Euler(self)
        """
        nprec = prec_to_dps(prec)
        return Euler(*[arg.evalf(n=nprec) for arg in self.args], seq=self.seq)
