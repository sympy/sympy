# References :
# http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/
# https://en.wikipedia.org/wiki/Quaternion

from __future__ import print_function

from sympy.core.expr import Expr
from sympy import Rational
from sympy import re, im, conjugate
from sympy import sqrt, sin, cos, acos, asin, exp, ln
from sympy import trigsimp
from sympy import diff, integrate
from sympy import Matrix, Add, Mul
from sympy import symbols, sympify
from sympy.printing.latex import latex
from sympy.printing import StrPrinter


class Quaternion(Expr):
    """Provides basic quaternion operations.
    Quaternion objects can be instantiated as Quaternion(a, b, c, d)
    as in (a + b*i + c*j + d*k).

    Example
    ========

    >>> from sympy.algebras.quaternion import Quaternion
    >>> q = Quaternion(1, 2, 3, 4)
    >>> q
    1 + 2*i + 3*j + 4*k
    """
    _op_priority = 11.0

    is_commutative = False

    def __new__(cls, a=0, b=0, c=0, d=0):
        a = sympify(a)
        b = sympify(b)
        c = sympify(c)
        d = sympify(d)

        if any(i.is_real is False for i in [a, b, c, d]):
            raise ValueError("arguments have to be real")
        else:
            obj = Expr.__new__(cls, a, b, c, d)

            obj._a = a
            obj._b = b
            obj._c = c
            obj._d = d
            return obj

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    @classmethod
    def from_axis_angle(cls, vector, angle):
        """Returns a rotation quaternion given the axis and the angle of rotation.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion.from_axis_angle((2, 4, 6), 2)
        >>> q
        cos(1) + sqrt(14)*sin(1)/14*i + sqrt(14)*sin(1)/7*j + 3*sqrt(14)*sin(1)/14*k
        """
        (x, y, z) = vector
        norm = sqrt(x**2 + y**2 + z**2)
        (x, y, z) = (x / norm, y / norm, z / norm)
        s = sin(angle * Rational(1, 2))
        a = cos(angle * Rational(1, 2))
        b = x * s
        c = y * s
        d = z * s

        return cls(a, b, c, d)

    @classmethod
    def from_matrix(cls, M):
        """Returns the equivalent quaternion of a matrix. The quaternion will be normalized
        only if the matrix is special orthogonal (orthogonal and det(M) = 1).

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> from sympy import Matrix
        >>> M = Matrix([[4, 6, 2], [1, 8, 6], [5, 2, 7]])
        >>> q = Quaternion.from_matrix(M)
        >>> q
        sqrt(238**(1/3) + 19)/2
        + 0*i
        + (-sqrt(-3 + 238**(1/3))/2)*j
        + (-sqrt(-5 + 238**(1/3))/2)*k
        """

        absQ = M.det()**Rational(1, 3)
        # The max( 0, ... ) is just a safeguard against rounding error.

        a = sqrt(max(0, absQ + M[0, 0] + M[1, 1] + M[2, 2])) / 2
        b = sqrt(max(0, absQ + M[0, 0] - M[1, 1] - M[2, 2])) / 2
        c = sqrt(max(0, absQ - M[0, 0] + M[1, 1] - M[2, 2])) / 2
        d = sqrt(max(0, absQ - M[0, 0] - M[1, 1] + M[2, 2])) / 2

        b = Quaternion.__copysign(b, M[2, 1] - M[1, 2])
        c = Quaternion.__copysign(c, M[0, 2] - M[2, 0])
        d = Quaternion.__copysign(d, M[1, 0] - M[0, 1])

        return Quaternion(a, b, c, d)

    @staticmethod
    def __copysign(x, y):

        # Takes the sign from the second term and sets the sign of the first
        # without altering the magnitude.

        if y == 0:
            return 0
        return x if x*y > 0 else -x

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.add(other*-1)

    def __mul__(self, other):
        return self._generic_mul(self, other)

    def __rmul__(self, other):
        return self._generic_mul(other, self)

    def __pow__(self, p):
        return self.pow(p)

    def __neg__(self):
        return Quaternion(-self._a, -self._b, -self._c, -self.d)

    def _eval_Integral(self, *args):
        return self.integrate(*args)

    def _eval_diff(self, *symbols, **kwargs):
        return self.diff(*symbols)

    def add(self, other):
        """Adds quaternions.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> from sympy import symbols
        >>> q1 = Quaternion(1, 2, 3, 4)
        >>> q2 = Quaternion(5, 6, 7, 8)
        >>> q1.add(q2)
        6 + 8*i + 10*j + 12*k
        >>> q1 + 5
        6 + 2*i + 3*j + 4*k
        >>> x = symbols('x', real = True)
        >>> q1.add(x)
        (x + 1) + 2*i + 3*j + 4*k
        """
        q1 = self
        q2 = sympify(other)

        # If q2 is a number or a sympy expression instead of a quaternion
        if not isinstance(q2, Quaternion):
            if q2.is_complex:
                return Quaternion(re(q2) + q1.a, im(q2) + q1.b, q1.c, q1.d)
            else:
                # q2 is something strange, do not evaluate:
                return Add(q1, q2)

        return Quaternion(q1.a + q2.a, q1.b + q2.b, q1.c + q2.c, q1.d
                          + q2.d)

    def mul(self, other):
        """Multiplies quaternions.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> from sympy import symbols
        >>> q1 = Quaternion(1, 2, 3, 4)
        >>> q2 = Quaternion(5, 6, 7, 8)
        >>> q1.mul(q2)
        (-60) + 12*i + 30*j + 24*k
        >>> q1.mul(2)
        2 + 4*i + 6*j + 8*k
        >>> x = symbols('x', real = True)
        >>> q1.mul(x)
        x + 2*x*i + 3*x*j + 4*x*k
        """
        return self._generic_mul(self, other)

    @staticmethod
    def _generic_mul(q1, q2):

        q1 = sympify(q1)
        q2 = sympify(q2)

        # None is a Quaternion:
        if not isinstance(q1, Quaternion) and not isinstance(q2, Quaternion):
            return q1 * q2

        # If q1 is a number or a sympy expression instead of a quaternion
        if not isinstance(q1, Quaternion):
            if q1.is_complex:
                return Quaternion(-im(q1)*q2.b + re(q1)*q2.a,
                                  im(q1)* q2.a + re(q1)*q2.b,
                                  -im(q1)*q2.d + re(q1)*q2.c,
                                  im(q1) * q2.c + re(q1)*q2.d)
            else:
                return Mul(q1, q2)

        # If q2 is a number or a sympy expression instead of a quaternion
        if not isinstance(q2, Quaternion):
            if q2.is_complex:
                return Quaternion(-q1.b*im(q2) + q1.a*re(q2),
                                  q1.b*re(q2) + q1.a*im(q2),
                                  q1.c*re(q2)+ q1.d*im(q2),
                                  -q1.c*im(q2) + q1.d*re(q2))
            else:
                return Mul(q1, q2)

        return Quaternion(-q1.b*q2.b - q1.c*q2.c - q1.d*q2.d + q1.a*q2.a,
                          q1.b*q2.a + q1.c*q2.d - q1.d*q2.c + q1.a*q2.b,
                          -q1.b*q2.d + q1.c*q2.a + q1.d*q2.b + q1.a*q2.c,
                          q1.b*q2.c - q1.c*q2.b + q1.d*q2.a + q1.a * q2.d)

    def _eval_conjugate(self):
        """Returns the conjugate of the quaternion."""
        q = self
        return Quaternion(q.a, -q.b, -q.c, -q.d)

    def norm(self):
        """Returns the norm of the quaternion."""
        q = self
        # trigsimp is used to simplify sin(x)^2 + cos(x)^2 (these terms
        # arise when from_axis_angle is used).
        return sqrt(trigsimp(q.a**2 + q.b**2 + q.c**2 + q.d**2))

    def normalize(self):
        """Returns the normalized form of the quaternion."""
        q = self
        return q * (1/q.norm())

    def inverse(self):
        """Returns the inverse of the quaternion."""
        q = self
        return conjugate(q) * (1/q.norm()**2)

    def pow(self, p):
        """Finds the pth power of the quaternion.
        Returns the inverse if p = -1.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.pow(4)
        668 + (-224)*i + (-336)*j + (-448)*k
        """
        q = self
        if p == -1:
            return q.inverse()
        res = 1
        while p > 0:
            if p & 1:
                res = q * res

            p = p >> 1
            q = q * q

        return res

    def exp(self):
        """Returns the exponential of q (e^q).

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.exp()
        E*cos(sqrt(29))
        + 2*sqrt(29)*E*sin(sqrt(29))/29*i
        + 3*sqrt(29)*E*sin(sqrt(29))/29*j
        + 4*sqrt(29)*E*sin(sqrt(29))/29*k
        """
        # exp(q) = e^a(cos||v|| + v/||v||*sin||v||)
        q = self
        vector_norm = sqrt(q.b**2 + q.c**2 + q.d**2)
        a = exp(q.a) * cos(vector_norm)
        b = exp(q.a) * sin(vector_norm) * q.b / vector_norm
        c = exp(q.a) * sin(vector_norm) * q.c / vector_norm
        d = exp(q.a) * sin(vector_norm) * q.d / vector_norm

        return Quaternion(a, b, c, d)

    def ln(self):
        """Returns the natural logarithm of the quaternion (ln(q)).

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.ln()
        log(sqrt(30))
        + 2*sqrt(29)*acos(sqrt(30)/30)/29*i
        + 3*sqrt(29)*acos(sqrt(30)/30)/29*j
        + 4*sqrt(29)*acos(sqrt(30)/30)/29*k
        """
        # ln(q) = ln||q|| + v/||v||*arccos(a/||q||)
        q = self
        vector_norm = sqrt(q.b**2 + q.c**2 + q.d**2)
        q_norm = q.norm()
        a = ln(q_norm)
        b = q.b * acos(q.a / q_norm) / vector_norm
        c = q.c * acos(q.a / q_norm) / vector_norm
        d = q.d * acos(q.a / q_norm) / vector_norm

        return Quaternion(a, b, c, d)

    def pow_cos_sin(self, p):
        """Computes the pth power in the cos-sin form.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.pow_cos_sin(4)
        900*cos(4*acos(sqrt(30)/30))
        + 1800*sqrt(29)*sin(4*acos(sqrt(30)/30))/29*i
        + 2700*sqrt(29)*sin(4*acos(sqrt(30)/30))/29*j
        + 3600*sqrt(29)*sin(4*acos(sqrt(30)/30))/29*k
        """
        # q = ||q||*(cos(a) + u*sin(a))
        # q^p = ||q||^p * (cos(p*a) + u*sin(p*a))

        q = self
        (v, angle) = q.to_axis_angle()
        q2 = Quaternion.from_axis_angle(v, p * angle)
        return q2 * (q.norm()**p)

    def diff(self, *args):
        return Quaternion(diff(self.a, *args), diff(self.b, *args),
                          diff(self.c, *args), diff(self.d, *args))

    def integrate(self, *args):
        # TODO: is this expression correct?
        return Quaternion(integrate(self.a, *args), integrate(self.b, *args),
                          integrate(self.c, *args), integrate(self.d, *args))

    @staticmethod
    def rotate(pin, r):
        """Returns the coordinates of the point pin(a 3 tuple) after rotation.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> Quaternion.rotate((1, 2, 4), q)
        (38/15, 8/3, 41/15)
        >>> (axis, angle) = q.to_axis_angle()
        >>> Quaternion.rotate((1, 2, 4), (axis, angle))
        (38/15, 8/3, 41/15)
        """
        if isinstance(r, tuple):
            # if r is of the form (vector, angle)
            q = Quaternion.from_axis_angle(r[0], r[1])
        else:
            # if r is a quaternion
            q = r.normalize()
        pout = q * Quaternion(0, pin[0], pin[1], pin[2]) * conjugate(q)
        return (pout.b, pout.c, pout.d)

    def to_axis_angle(self):
        """Returns the axis and angle of rotation of a quaternion

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> (axis, angle) = q.to_axis_angle()
        >>> axis
        (2*sqrt(29)/29, 3*sqrt(29)/29, 4*sqrt(29)/29)
        >>> angle
        2*acos(sqrt(30)/30)
        """
        q = self
        try:
            # Skips it if it doesn't know whether q.a is negative
            if q.a < 0:
                # avoid error with acos
                # axis and angle of rotation of q and q*-1 will be the same
                q = q * -1
        except BaseException:
            pass

        q = q.normalize()
        angle = trigsimp(2 * acos(q.a))

        # Since quaternion is normalised, q.a is less than 1.
        s = sqrt(1 - q.a*q.a)

        x = trigsimp(q.b / s)
        y = trigsimp(q.c / s)
        z = trigsimp(q.d / s)

        v = (x, y, z)
        t = (v, angle)

        return t

    def to_matrix(self):
        """ Returns the equivalent rotation transformation matrix of the quaternion

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.to_matrix()
        Matrix([
        [-2/5,  2/5, 11/5],
        [   2, -1/5,    2],
        [   1, 14/5,    0]])
        """
        q = self
        sqa = sqrt(q.a**2)
        sqb = sqrt(q.b**2)
        sqc = sqrt(q.c**2)
        sqd = sqrt(q.d**2)

        invs = 1 / (sqa + sqb + sqc + sqd)
        m00 = (sqb - sqc - sqd + sqa) * invs
        m11 = (-sqb + sqc - sqd + sqa) * invs
        m22 = (-sqb - sqc + sqd + sqa) * invs

        tmp1 = q.b * q.c
        tmp2 = q.d * q.a
        m10 = 2 * (tmp1 + tmp2) * invs
        m01 = 2 * (tmp1 - tmp2) * invs

        tmp1 = q.b * q.d
        tmp2 = q.c * q.a
        m20 = 2 * (tmp1 - tmp2) * invs
        m02 = 2 * (tmp1 + tmp2) * invs

        tmp1 = q.c * q.d
        tmp2 = q.b * q.a
        m21 = 2 * (tmp1 + tmp2) * invs
        m12 = 2 * (tmp1 - tmp2) * invs

        return Matrix([[m00, m01, m02], [m10, m11, m12], [m20, m21,
                                                          m22]])

    def to_matrix_4x4(self, v):
        """ Generates a 4x4 transformation matrix (used for rotation about a point
        other than the origin) from a quaternion and the position vector of the point.

        Example
        ========

        >>> from sympy.algebras.quaternion import Quaternion
        >>> q = Quaternion(1, 2, 3, 4)
        >>> q.to_matrix_4x4((1, 4, 8))
        Matrix([
        [-2/5,  2/5, 11/5, -89/5],
        [   2, -1/5,    2, -66/5],
        [   1, 14/5,    0, -21/5],
        [   0,    0,    0,     1]])
        """

        q = self
        sqa = sqrt(q.a**2)
        sqb = sqrt(q.b**2)
        sqc = sqrt(q.c**2)
        sqd = sqrt(q.d**2)

        invs = 1 / (sqa + sqb + sqc + sqd)
        m00 = (sqb - sqc - sqd + sqa) * invs
        m11 = (-sqb + sqc - sqd + sqa) * invs
        m22 = (-sqb - sqc + sqd + sqa) * invs

        tmp1 = q.b * q.c
        tmp2 = q.d * q.a
        m10 = 2 * (tmp1 + tmp2) * invs
        m01 = 2 * (tmp1 - tmp2) * invs

        tmp1 = q.b * q.d
        tmp2 = q.c * q.a
        m20 = 2 * (tmp1 - tmp2) * invs
        m02 = 2 * (tmp1 + tmp2) * invs

        tmp1 = q.c * q.d
        tmp2 = q.b * q.a
        m21 = 2 * (tmp1 + tmp2) * invs
        m12 = 2 * (tmp1 - tmp2) * invs

        (x, y, z) = v

        m03 = x - x*m00 - y*m01 - z*m02
        m13 = y - x*m10 - y*m11 - z*m12
        m23 = z - x*m20 - y*m21 - z*m22
        m30 = m31 = m32 = 0
        m33 = 1

        return Matrix([[m00, m01, m02, m03], [m10, m11, m12, m13],
                       [m20, m21, m22, m23], [m30, m31, m32, m33]])
