from sympy.matrices.dense import MutableDenseMatrix as Matrix
from sympy.algebras.quaternion import Quaternion
from sympy.core.sympify import sympify
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
    Quaternion objects can be instantiated as Euler(alpha, beta, gamma, seq).
    as in (a + b*i + c*j + d*k).

    Parameters
    ==========

    seq : string of length 3
        Represents the sequence of rotations.
        For intrinsic rotations, seq must be all lowercase and its elements
        must be from the set `{'x', 'y', 'z'}`
        For extrinsic rotations, seq must be all uppercase and its elements
        must be from the set `{'X', 'Y', 'Z'}`

    Examples
    ========


    References
    ==========

    .. [2] https://en.wikipedia.org/wiki/Quaternion

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
        Euler(pi/2,  0,  0)^{ZYX}

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
        Euler(pi/2,  0,  0)^{ZYX}

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

    def inverse(self):
        """Returns inverse Euler angles"""
        seq = self.seq[::-1]
        alpha, beta, gamma = [-i for i in self.args[::-1]]
        return Euler(alpha, beta, gamma, seq)

    @property
    def is_zero_rotation(self):
        """Returns true if the angles are all zero"""
        return all(i.is_zero for i in self.args)

    def _eval_evalf(self, prec):
        """Returns the floating point approximations (decimal numbers) of the angles.

        Returns
        =======

        Euler
            Floating point approximations of Euler(self)

        """
        nprec = prec_to_dps(prec)
        return Euler(*[arg.evalf(n=nprec) for arg in self.args], seq=self.seq)
