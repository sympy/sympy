# References :
# https://en.wikipedia.org/wiki/Dual_number
from sympy import conjugate, exp, ln, sympify
from sympy.core.expr import Expr
from sympy.algebras.quaternion import Quaternion


class Dual(Expr):
    """Provides basic dual number operations. Dual objects can be instantiated
    as Dual(a, b), as in (a + b*d), where d^2 = 0.

    Dual numbers can be initiated with Quaternions to create Dual Quaternions,
    the algebra of all rigid body kinematics in 3D. In other words, translations
    and rotations. Another important use of dual numbers is automatic
    differentiation.

    Examples
    ========

    Dual Quaternions are created by initiating a dual number from quaternions:

    >>> a, b, c, d, e, f, g, h = symbols('a:h')
    >>> p = Quaternion(a, b, c, d)
    >>> q = Quaternion(e, f, g, h)
    >>> x = Dual(p, q)
    >>> x
    a + b*i + c*j + d*k + (e + f*i + g*j + h*k)*ð
    >>> x*x.inverse()
    1 + 0*i + 0*j + 0*k + (0 + 0*i + 0*j + 0*k)*ð

    Take note that operations like conjugate only perform the dual conjugate,
    not the quaternionic conjugate.

    >>> conjugate(x)
    a + b*i + c*j + d*k + ((-e) + (-f)*i + (-g)*j + (-h)*k)*ð

    As an example of automatic differentiation consider the simple function x^2.
    Evaluating the function with a dual number automatically performs the chain
    rulewhile the function itself is being evaluated.

    >>> from sympy import symbols
    >>> from sympy.algebras.dual import Dual
    >>> a = symbols('a')
    >>> d = Dual(a, 1)
    >>> d
    a + 1*d
    >>> f = lambda x: x**2
    >>> f(d)
    a**2 + 2*a*d

    We recognize the term proportional to d as the derivative of the function.
    To obtain only the derivative we select only the dual component:

    >>> f(d).b
    2*a
    """
    _op_priority = 11.0

    is_commutative = True

    def __new__(cls, a=0, b=0):
        a = sympify(a)
        b = sympify(b)

        obj = Expr.__new__(cls, a, b)
        obj._a = a
        obj._b = b
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            # Dual quaterions are not a commutative algebra
            obj.is_commutative = False
        return obj

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

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
        return Dual(-self._a, -self._b)

    def __truediv__(self, other):
        return self * sympify(other)**-1

    __div__ = __truediv__

    def __rtruediv__(self, other):
        return sympify(other) * self**-1

    __rdiv__ = __rtruediv__

    def diff(self, *symbols, **kwargs):
        kwargs.setdefault('evaluate', True)
        return self.func(*[a.diff(*symbols, **kwargs) for a in self.args])

    def add(self, other):
        """Adds dual numbers.

        Parameters
        ==========

        other : Dual
            The dual number to add to current (self) dual number.

        Returns
        =======

        Dual
            The resultant dual number after adding self to other

        Examples
        ========

        >>> from sympy.algebras.Dual import Dual
        >>> d1 = Dual(1, 2)
        >>> d2 = Dual(5, 6)
        >>> d1.add(d2)
        6 + 8*d
        >>> d1 + 5
        6 + 2*d

        """
        d1 = self
        d2 = sympify(other)

        if not isinstance(d2, Dual):
            return Dual(d1.a + d2, d1.b)

        return Dual(d1.a + d2.a, d1.b + d2.b)

    def mul(self, other):
        """Multiplies dual numbers.

        Parameters
        ==========

        other : Dual or symbol
            The dual number to multiply to current (self) dual number.

        Returns
        =======

        Dual
            The resultant dual number after multiplying self with other

        Examples
        ========

        >>> from sympy.algebras.dual import Dual
        >>> d1 = Dual(1, 2)
        >>> d2 = Dual(5, 6)
        >>> d1.mul(d2)
        5 + 16*d
        >>> d1.mul(2)
        2 + 4*d

        """
        return self._generic_mul(self, other)

    @staticmethod
    def _generic_mul(d1, d2):
        """Generic multiplication.

        Parameters
        ==========

        d1 : Dual or symbol
        d2 : Dual or symbol

        It's important to note that if neither d1 nor d2 is a Dual,
        this function simply returns d1 * d2.

        Returns
        =======

        Dual
            The resultant dual number after multiplying d1 and d2

        Examples
        ========

        >>> from sympy.algebras.dual import Dual
        >>> d1 = Dual(1, 2)
        >>> d2 = Dual(5, 6)
        >>> d1.mul(d2)
        5 + 16*d
        >>> d1.mul(2)
        2 + 4*d

        """
        d1 = sympify(d1)
        d2 = sympify(d2)

        # None is a Dual:
        if not isinstance(d1, Dual) and not isinstance(d2, Dual):
            return d1 * d2

        # If d1 is a number or a sympy expression instead of a dual number
        if not isinstance(d1, Dual):
            return Dual(d1 * d2.a, d1 * d2.b)

        # If d2 is a number or a sympy expression instead of a dual number
        if not isinstance(d2, Dual):
            return Dual(d1.a * d2, d1.b * d2)

        return Dual(d1.a * d2.a, d1.a * d2.b + d1.b * d2.a)

    def _eval_conjugate(self):
        """Returns the conjugate of the dual number."""
        q = self
        return Dual(q.a, - q.b)

    def norm(self):
        """Returns the norm of the dual number."""
        q = self
        return q.a

    def normalize(self):
        """Returns the normalized form of the dual number."""
        q = self
        return q * q.norm()**-1

    def inverse(self):
        """Returns the inverse of the dual number."""
        q = self
        if not q.norm():
            raise ValueError("Cannot compute inverse for a dual number with zero norm")
        qainv = q.a ** -1
        return Dual(qainv, -qainv * q.b * qainv)

    def pow(self, p):
        """Finds the pth power of the dual number.

        Parameters
        ==========

        p : int
            Power to be applied on dual number.

        Returns
        =======

        Dual
            Returns the p-th power of the current dual number.
            Returns the inverse if p = -1.

        Examples
        ========

        >>> from sympy.algebras.dual import Dual
        >>> q = Dual(1, 2)
        >>> q.pow(4)
        1 + 8*d

        """
        p = sympify(p)
        q = self
        if p == -1:
            return q.inverse()


        if isinstance(q.a, Quaternion) and isinstance(q.b, Quaternion):
            res = 1

            if not p.is_Integer:
                return NotImplemented

            if p < 0:
                q, p = q.inverse(), -p

            while p > 0:
                res *= q
                p -= 1

            return res
        else:
            if p.is_Integer:
                return Dual(q.a**p, q.b * p * q.a ** (p-1))
            else:
                return (p*q._ln()).exp()

    def exp(self):
        """Returns the exponential of q (e^q).

        Returns
        =======

        Dual
            Exponential of q (e^q).

        Examples
        ========

        >>> from sympy.algebras.dual import Dual
        >>> q = Dual(3, 4)
        >>> q.exp()
        exp(3) + 4*exp(3)*d

        """
        q = self
        exp_qa = exp(q.a) if not isinstance(q.a, Dual) else q.a.exp()
        return Dual(exp_qa, q.b * exp_qa)

    def _ln(self):
        """Returns the natural logarithm of the dual number (_ln(q)).

        Examples
        ========

        >>> from sympy.algebras.dual import Dual
        >>> q = Dual(3, 4)
        >>> q._ln()
        log(3) + 4/3*d

        """
        q = self
        ln_qa = ln(q.a) if not isinstance(q.a, Dual) else q.a._ln()
        return Dual(ln_qa, q.b/q.a)
