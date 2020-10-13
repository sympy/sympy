from sympy.core.sympify import _sympify
from sympy.matrices.expressions import MatrixExpr
from sympy import S, I, sqrt, exp


class DFT(MatrixExpr):
    r"""
    Returns a discrete Fourier transform matrix.

    Parameters
    ==========

    n : integer or Symbol
        Size of the transform.
    unitary : bool
        If True (default), the matrix is scaled with :math:`\frac{1}{\sqrt{N}}`.
    twiddle : Symbol or None
        Use a symbol instead of :math:`\exp(-j2\pi/n)` as twiddle factor.
        Default value: None

    Examples
    ========

    >>> from sympy.abc import n
    >>> from sympy.matrices.expressions.fourier import DFT
    >>> DFT(3)
    DFT(3)
    >>> DFT(3).as_explicit()
    Matrix([
    [sqrt(3)/3,                sqrt(3)/3,                sqrt(3)/3],
    [sqrt(3)/3, sqrt(3)*exp(-2*I*pi/3)/3,  sqrt(3)*exp(2*I*pi/3)/3],
    [sqrt(3)/3,  sqrt(3)*exp(2*I*pi/3)/3, sqrt(3)*exp(-2*I*pi/3)/3]])
    >>> DFT(n).shape
    (n, n)
    >>> DFT(3, unitary=False).as_explicit()
    Matrix([
    [1,              1,              1],
    [1, exp(-2*I*pi/3),  exp(2*I*pi/3)],
    [1,  exp(2*I*pi/3), exp(-2*I*pi/3)]])

    It is possible to use a symbolic twiddle factor

    >>> from sympy.abc import omega
    >>> DFT(4, unitary=False, twiddle=omega).as_explicit()
    Matrix([
    [1,        1,        1,        1],
    [1,    omega, omega**2, omega**3],
    [1, omega**2, omega**4, omega**6],
    [1, omega**3, omega**6, omega**9]])

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/DFT_matrix

    """

    def __new__(cls, n, unitary=True, twiddle=None):
        n = _sympify(n)
        cls._check_dim(n)

        obj = super().__new__(cls, n)
        cls._unitary = unitary
        cls._twiddle = twiddle
        return obj

    n = property(lambda self: self.args[0])  # type: ignore
    shape = property(lambda self: (self.n, self.n))  # type: ignore

    def _entry(self, i, j, **kwargs):
        if self._twiddle is None:
            w = exp(-2*S.Pi*I/self.n)
        else:
            w = self._twiddle
        if self._unitary:
            return w**(i*j) / sqrt(self.n)
        else:
            return w**(i*j)

    def _eval_inverse(self):
        return IDFT(self.n, unitary=self._unitary, twiddle=self._twiddle)


class IDFT(DFT):
    r"""
    Returns an inverse discrete Fourier transform matrix.

    Parameters
    ==========

    n : integer or Symbol
        Size of the transform
    unitary : bool
        If True (default), the matrix is scaled with :math:`\frac{1}{\sqrt{N}}`. If
        False, the matrix is scaled with :math:`\frac{1}{N}`.
    twiddle : Symbol or None
        Use a symbol instead of :math:`\exp(-j2\pi/n)` as twiddle factor.
        Default value: None

    Examples
    ========

    >>> from sympy.matrices.expressions.fourier import DFT, IDFT
    >>> IDFT(3)
    IDFT(3)
    >>> IDFT(4)*DFT(4)
    I

    See Also
    ========

    sympy.matrices.expressions.fourier.DFT

    """
    def _entry(self, i, j, **kwargs):
        if self._twiddle is None:
            w = exp(-2*S.Pi*I/self.n)
        else:
            w = self._twiddle
        if self._unitary:
            return w**(-i*j) / sqrt(self.n)
        else:
            return w**(-i*j) / self.n

    def _eval_inverse(self):
        return DFT(self.n, unitary=self._unitary, twiddle=self._twiddle)
