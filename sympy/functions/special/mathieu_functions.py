""" This module contains the Mathieu functions.
"""

from __future__ import print_function, division

from sympy.core import S
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import sin, cos


class mathieus(Function):
    r"""
    The Mathieu Sine function. This function is one solution
    of the Mathieu differential equation:

    .. math ::
        y(x)^{\prime\prime} + (a - 2 q \cos(2 x) y(x) = 0

    The other solution is the Mathieu Cosine function.

    Examples
    ========

    >>> from sympy import I, oo, mathieus
    >>> from sympy.abc import a, q, z

    >>> mathieus(a, q, z)
    mathieus(a, q, z)

    >>> mathieus(a, 0, z)
    sin(sqrt(a)*z)

    See Also
    ========

    mathieuc: Mathieu cosine function
    mathieusprime: Derivative of Mathieu sine function
    mathieucprime: Derivative of Mathieu cosine function

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathieu_function
    .. [2] http://dlmf.nist.gov/28
    .. [3] http://mathworld.wolfram.com/MathieuFunction.html
    .. [4] http://functions.wolfram.com/MathieuandSpheroidalFunctions/MathieuS/
    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 3:
            a, q, z = self.args
            return mathieusprime(a, q, z)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, a, q, z):
        if q.is_Number and q is S.Zero:
            return sin(sqrt(a)*z)
        # Try to pull out factors of -1
        if z.could_extract_minus_sign():
            return -cls(a, q, -z)

    def _eval_conjugate(self):
        a, q, z = self.args
        return self.func(a.conjugate(), q.conjugate(), z.conjugate())


class mathieuc(Function):
    r"""
    The Mathieu Cosine function. This function is one solution
    of the Mathieu differential equation:

    .. math ::
        y(x)^{\prime\prime} + (a - 2 q \cos(2 x) y(x) = 0

    The other solution is the Mathieu Sine function.

    Examples
    ========

    >>> from sympy import I, oo, mathieus
    >>> from sympy.abc import a, q, z

    >>> mathieuc(a, q, z)
    mathieuc(a, q, z)

    >>> mathieuc(a, 0, z)
    cos(sqrt(a)*z)

    See Also
    ========

    mathieus: Mathieu sine function
    mathieusprime: Derivative of Mathieu sine function
    mathieucprime: Derivative of Mathieu cosine function

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathieu_function
    .. [2] http://dlmf.nist.gov/28
    .. [3] http://mathworld.wolfram.com/MathieuFunction.html
    .. [4] http://functions.wolfram.com/MathieuandSpheroidalFunctions/MathieuC/
    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 3:
            a, q, z = self.args
            return mathieucprime(a, q, z)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, a, q, z):
        if q.is_Number and q is S.Zero:
            return cos(sqrt(a)*z)
        # Try to pull out factors of -1
        if z.could_extract_minus_sign():
            return cls(a, q, -z)

    def _eval_conjugate(self):
        a, q, z = self.args
        return self.func(a.conjugate(), q.conjugate(), z.conjugate())
