from typing import Type
from sympy import Interval, numer, Rational, solveset
from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.evalf import EvalfMixin
from sympy.core.expr import Expr
from sympy.core.function import expand
from sympy.core.mul import Mul
from sympy.core.numbers import I, pi, oo
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import Dummy, Symbol
from sympy.functions import Abs
from sympy.core.sympify import sympify, _sympify
from sympy.matrices import (Matrix, ImmutableMatrix, ImmutableDenseMatrix, eye,
                            ShapeError, zeros)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.matrices.expressions import MatMul, MatAdd
from sympy.polys import Poly, rootof
from sympy.polys.polyroots import roots
from sympy.polys.polytools import (cancel, degree)
from sympy.polys.domains import EXRAW
from sympy.series import limit
from sympy.utilities.misc import filldedent
from sympy.solvers.ode.systems import linodesolve
from sympy.solvers.solveset import linsolve, linear_eq_to_matrix
from sympy.logic.boolalg import false, true, Boolean
from sympy.solvers.inequalities import reduce_inequalities
from abc import ABC, abstractmethod

from mpmath.libmp.libmpf import prec_to_dps

__all__ = ['TransferFunction', 'DiscreteTransferFunction',
           'create_transfer_function','PIDController', 'Series',
           'MIMOSeries', 'Parallel', 'MIMOParallel', 'Feedback', 'MIMOFeedback',
           'TransferFunctionMatrix', 'StateSpace', 'DiscreteStateSpace',
           'create_state_space', 'gbt', 'bilinear', 'forward_diff',
           'backward_diff', 'phase_margin', 'gain_margin']

def _roots(poly, var):
    """ like roots, but works on higher-order polynomials. """
    r = roots(poly, var, multiple=True)
    n = degree(poly)
    if len(r) != n:
        r = [rootof(poly, var, k) for k in range(n)]
    return r

def gbt(tf, sample_per, alpha):
    r"""
    Returns falling coefficients of H(z) from numerator and denominator.

    Parameters
    ==========

    tf : TransferFunction
        The continuous transfer function H(s) to be discretized.
    sample_per : Symbol, Number
        Time interval between two consecutive sampling instants.
    alpha: Symbol, Number
        The parameter for the generalised bilinear transformation method.

    Explanation
    ===========

    Where H(z) is the corresponding discretized transfer function,
    discretized with the generalised bilinear transformation method.
    H(z) is obtained from the continuous transfer function H(s)
    by substituting $s(z) = \frac{z-1}{T(\alpha z + (1-\alpha))}$ into H(s),
    where T is the sample period.
    Coefficients are falling, i.e. $H(z) = \frac{az+b}{cz+d}$ is returned
    as [a, b], [c, d].

    Examples
    ========

    >>> from sympy.physics.control.lti import TransferFunction, gbt
    >>> from sympy.abc import s, L, R, T

    >>> tf = TransferFunction(1, s*L + R, s)
    >>> numZ, denZ = gbt(tf, T, 0.5)
    >>> numZ
    [T/(2*(L + R*T/2)), T/(2*(L + R*T/2))]
    >>> denZ
    [1, (-L + R*T/2)/(L + R*T/2)]

    >>> numZ, denZ = gbt(tf, T, 0)
    >>> numZ
    [T/L]
    >>> denZ
    [1, (-L + R*T)/L]

    >>> numZ, denZ = gbt(tf, T, 1)
    >>> numZ
    [T/(L + R*T), 0]
    >>> denZ
    [1, -L/(L + R*T)]

    >>> numZ, denZ = gbt(tf, T, 0.3)
    >>> numZ
    [3*T/(10*(L + 3*R*T/10)), 7*T/(10*(L + 3*R*T/10))]
    >>> denZ
    [1, (-L + 7*R*T/10)/(L + 3*R*T/10)]

    References
    ==========

    .. [1] https://www.polyu.edu.hk/ama/profile/gfzhang/Research/ZCC09_IJC.pdf
    """
    if not tf.is_SISO:
        raise NotImplementedError("Not implemented for MIMO systems.")

    T = sample_per  # and sample period T
    s = tf.var
    z =  s         # dummy discrete variable z

    np = tf.num.as_poly(s).all_coeffs()
    dp = tf.den.as_poly(s).all_coeffs()
    alpha = Rational(alpha).limit_denominator(1000)

    # The next line results from multiplying H(z) with z^N/z^N
    N = max(len(np), len(dp)) - 1
    num = Add(*[ T**(N-i) * c * (z-1)**i * (alpha * z + 1 - alpha)**(N-i) for c, i in zip(np[::-1], range(len(np))) ])
    den = Add(*[ T**(N-i) * c * (z-1)**i * (alpha * z + 1 - alpha)**(N-i) for c, i in zip(dp[::-1], range(len(dp))) ])

    num_coefs = num.as_poly(z).all_coeffs()
    den_coefs = den.as_poly(z).all_coeffs()

    para = den_coefs[0]
    num_coefs = [coef/para for coef in num_coefs]
    den_coefs = [coef/para for coef in den_coefs]

    return num_coefs, den_coefs

def bilinear(tf, sample_per):
    r"""
    Returns falling coefficients of H(z) from numerator and denominator.

    Parameters
    ==========

    tf : TransferFunction
        The continuous transfer function H(s) to be discretized.
    sampling_time : Symbol, Number
        Time interval between two consecutive sampling instants.

    Explanation
    ===========

    Where H(z) is the corresponding discretized transfer function,
    discretized with the bilinear transform method.
    H(z) is obtained from the continuous transfer function H(s)
    by substituting $s(z) = \frac{2}{T}\frac{z-1}{z+1}$ into H(s), where T is
    the sample period.
    Coefficients are falling, i.e. $H(z) = \frac{az+b}{cz+d}$ is returned
    as [a, b], [c, d].

    Examples
    ========

    >>> from sympy.physics.control.lti import TransferFunction, bilinear
    >>> from sympy.abc import s, L, R, T

    >>> tf = TransferFunction(1, s*L + R, s)
    >>> numZ, denZ = bilinear(tf, T)
    >>> numZ
    [T/(2*(L + R*T/2)), T/(2*(L + R*T/2))]
    >>> denZ
    [1, (-L + R*T/2)/(L + R*T/2)]
    """
    return gbt(tf, sample_per, S.Half)

def forward_diff(tf, sample_per):
    r"""
    Returns falling coefficients of H(z) from numerator and denominator.

    Parameters
    ==========

    tf : TransferFunction
        The continuous transfer function H(s) to be discretized.
    sampling_time : Symbol, Number
        Time interval between two consecutive sampling instants.

    Explanation
    ===========

    Where H(z) is the corresponding discretized transfer function,
    discretized with the forward difference transform method.
    H(z) is obtained from the continuous transfer function H(s)
    by substituting $s(z) = \frac{z-1}{T}$ into H(s), where T is the
    sample period.
    Coefficients are falling, i.e. $H(z) = \frac{az+b}{cz+d}$ is returned
    as [a, b], [c, d].

    Examples
    ========

    >>> from sympy.physics.control.lti import TransferFunction, forward_diff
    >>> from sympy.abc import s, L, R, T

    >>> tf = TransferFunction(1, s*L + R, s)
    >>> numZ, denZ = forward_diff(tf, T)
    >>> numZ
    [T/L]
    >>> denZ
    [1, (-L + R*T)/L]
    """
    return gbt(tf, sample_per, S.Zero)

def backward_diff(tf, sample_per):
    r"""
    Returns falling coefficients of H(z) from numerator and denominator.

    Parameters
    ==========

    tf : TransferFunction
        The continuous transfer function H(s) to be discretized.
    sampling_time : Symbol, Number
        Time interval between two consecutive sampling instants.

    Explanation
    ===========

    Where H(z) is the corresponding discretized transfer function,
    discretized with the backward difference transform method.
    H(z) is obtained from the continuous transfer function H(s)
    by substituting $s(z) =  \frac{z-1}{Tz}$ into H(s), where T is the
    sample period.
    Coefficients are falling, i.e. $H(z) = \frac{az+b}{cz+d}$ is returned
    as [a, b], [c, d].

    Examples
    ========

    >>> from sympy.physics.control.lti import TransferFunction, backward_diff
    >>> from sympy.abc import s, L, R, T

    >>> tf = TransferFunction(1, s*L + R, s)
    >>> numZ, denZ = backward_diff(tf, T)
    >>> numZ
    [T/(L + R*T), 0]
    >>> denZ
    [1, -L/(L + R*T)]
    """
    return gbt(tf, sample_per, S.One)

def phase_margin(system):
    r"""
    Returns the phase margin of a continuous time system.
    Only applicable to Transfer Functions which can generate valid bode plots.

    Raises
    ======

    NotImplementedError
        When time delay terms are present in the system.

    ValueError
        When a SISO LTI system is not passed.

        When more than one free symbol is present in the system.
        The only variable in the transfer function should be
        the variable of the Laplace transform.

    Examples
    ========

    >>> from sympy.physics.control import TransferFunction, phase_margin
    >>> from sympy.abc import s

    >>> tf = TransferFunction(1, s**3 + 2*s**2 + s, s)
    >>> phase_margin(tf)
    180*(-pi + atan((-1 + (-2*18**(1/3)/(9 + sqrt(93))**(1/3) + 12**(1/3)*(9 + sqrt(93))**(1/3))**2/36)/(-12**(1/3)*(9 + sqrt(93))**(1/3)/3 + 2*18**(1/3)/(3*(9 + sqrt(93))**(1/3)))))/pi + 180
    >>> phase_margin(tf).n()
    21.3863897518751

    >>> tf1 = TransferFunction(s**3, s**2 + 5*s, s)
    >>> phase_margin(tf1)
    -180 + 180*(atan(sqrt(2)*(-51/10 - sqrt(101)/10)*sqrt(1 + sqrt(101))/(2*(sqrt(101)/2 + 51/2))) + pi)/pi
    >>> phase_margin(tf1).n()
    -25.1783920627277

    >>> tf2 = TransferFunction(1, s + 1, s)
    >>> phase_margin(tf2)
    -180

    See Also
    ========

    gain_margin

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Phase_margin

    """
    from sympy.functions import arg

    if not isinstance(system, SISOLinearTimeInvariant):
        raise ValueError("Margins are only applicable for SISO LTI systems.")

    _w = Dummy("w", real=True)
    repl = I*_w
    expr = system.to_expr()
    len_free_symbols = len(expr.free_symbols)
    if expr.has(exp):
        raise NotImplementedError("Margins for systems with Time delay terms are not supported.")
    elif len_free_symbols > 1:
        raise ValueError("Extra degree of freedom found. Make sure"
            " that there are no free symbols in the dynamical system other"
            " than the variable of Laplace transform.")

    w_expr = expr.subs({system.var: repl})

    mag = 20*log(Abs(w_expr), 10)
    mag_sol = list(solveset(mag, _w, Interval(0, oo, left_open=True)))

    if (len(mag_sol) == 0):
        pm = S(-180)
    else:
        wcp = mag_sol[0]
        pm = ((arg(w_expr)*S(180)/pi).subs({_w:wcp}) + S(180)) % 360

    if(pm >= 180):
        pm = pm - 360

    return pm

def gain_margin(system):
    r"""
    Returns the gain margin of a continuous time system.
    Only applicable to Transfer Functions which can generate valid bode plots.

    Raises
    ======

    NotImplementedError
        When time delay terms are present in the system.

    ValueError
        When a SISO LTI system is not passed.

        When more than one free symbol is present in the system.
        The only variable in the transfer function should be
        the variable of the Laplace transform.

    Examples
    ========

    >>> from sympy.physics.control import TransferFunction, gain_margin
    >>> from sympy.abc import s

    >>> tf = TransferFunction(1, s**3 + 2*s**2 + s, s)
    >>> gain_margin(tf)
    20*log(2)/log(10)
    >>> gain_margin(tf).n()
    6.02059991327962

    >>> tf1 = TransferFunction(s**3, s**2 + 5*s, s)
    >>> gain_margin(tf1)
    oo

    See Also
    ========

    phase_margin

    References
    ==========

    https://en.wikipedia.org/wiki/Bode_plot

    """
    if not isinstance(system, SISOLinearTimeInvariant):
        raise ValueError("Margins are only applicable for SISO LTI systems.")

    _w = Dummy("w", real=True)
    repl = I*_w
    expr = system.to_expr()
    len_free_symbols = len(expr.free_symbols)
    if expr.has(exp):
        raise NotImplementedError("Margins for systems with Time delay terms are not supported.")
    elif len_free_symbols > 1:
        raise ValueError("Extra degree of freedom found. Make sure"
            " that there are no free symbols in the dynamical system other"
            " than the variable of Laplace transform.")

    w_expr = expr.subs({system.var: repl})

    mag = 20*log(Abs(w_expr), 10)
    phase = w_expr
    phase_sol = list(solveset(numer(phase.as_real_imag()[1].cancel()),_w, Interval(0, oo, left_open = True)))

    if (len(phase_sol) == 0):
        gm = oo
    else:
        wcg = phase_sol[0]
        gm = -mag.subs({_w:wcg})

    return gm

class LinearTimeInvariant(Basic, EvalfMixin, ABC):
    """A common class for all the Linear Time-Invariant Dynamical Systems."""

    _clstype: Type

    # Users should not directly interact with this class.
    def __new__(cls, *system, **kwargs):
        if cls is LinearTimeInvariant:
            raise NotImplementedError('The LTICommon class is not meant to be used directly.')
        return super(LinearTimeInvariant, cls).__new__(cls, *system, **kwargs)

    @classmethod
    def _check_args(cls, args):
        """
        Check if the arguments passed to the class are valid.
        Every argument must be of the same type as the class (_clstype).
        All arguments must have the same complex variable of the Laplace
        transform.

        """
        if not args:
            raise ValueError("At least 1 argument must be passed.")
        if not all(isinstance(arg, cls._clstype) for arg in args):
            raise TypeError(f"All arguments must be of type {cls._clstype}.")
        var_set = {arg.var for arg in args}
        if len(var_set) != 1:
            raise ValueError(filldedent(f"""
                All transfer functions should use the same complex variable
                of the Laplace transform or z-transform. {len(var_set)} different
                values found."""))

    @property
    def is_continuous(self):
        """Returns ``True`` if the passed LTI system is continuous,
        else returns ``False``."""
        return self._is_continuous

    @property
    @abstractmethod
    def sampling_time(self):
        """Returns the sampling time of the passed LTI system."""
        pass

    @property
    def is_SISO(self):
        """Returns ``True`` if the passed LTI system is SISO else returns False."""
        return self._is_SISO


class SISOLinearTimeInvariant(LinearTimeInvariant):
    """A common class for all the SISO Linear Time-Invariant Dynamical Systems."""
    # Users should not directly interact with this class.

    @property
    def num_inputs(self):
        """Return the number of inputs for SISOLinearTimeInvariant."""
        return 1

    @property
    def num_outputs(self):
        """Return the number of outputs for SISOLinearTimeInvariant."""
        return 1

    _is_SISO = True


class MIMOLinearTimeInvariant(LinearTimeInvariant):
    """A common class for all the MIMO Linear Time-Invariant Dynamical Systems."""
    # Users should not directly interact with this class.
    _is_SISO = False


SISOLinearTimeInvariant._clstype = SISOLinearTimeInvariant
MIMOLinearTimeInvariant._clstype = MIMOLinearTimeInvariant


def _check_other_SISO(func):
    def wrapper(*args, **kwargs):
        if not isinstance(args[-1], SISOLinearTimeInvariant):
            return NotImplemented
        else:
            return func(*args, **kwargs)
    return wrapper


def _check_other_MIMO(func):
    def wrapper(*args, **kwargs):
        if not isinstance(args[-1], MIMOLinearTimeInvariant):
            return NotImplemented
        else:
            return func(*args, **kwargs)
    return wrapper


def _check_time_compatibility(systems):
    """
    Checks compatibility between different systems
    Every system must be either continuous-time or discrete-time and
    must have the same sampling time if discrete-time.

    """
    if not all(isinstance(system, LinearTimeInvariant) \
               for system in systems):
        return

    continuous = systems[0].is_continuous
    sampling_time = systems[0].sampling_time
    for system in systems[1:]:
        if system.is_continuous != continuous:
            raise TypeError("""
                Incompatible systems detected.
                All systems should be either continuous-time or discrete-time.
                Found: {} and {}""".format(systems[0], system,))
        if system.sampling_time != sampling_time:
            raise TypeError("""
                Incompatible sampling times detected.
                All systems should have the same sampling time.
                Found sampling time: {} and {}""".\
                    format(systems[0].sampling_time, system.sampling_time))

def _compatibility_decorator(func):
    """
    Decorator to check compatibility of systems before performing operations.
    """
    def wrapper(self, *args, **kwargs):
        _check_time_compatibility([self] + list(args))
        return func(self, *args, **kwargs)
    return wrapper

def create_transfer_function(num, den, var, sampling_time=0):
    """
    Creates a new transfer function object.
    sampling_time == 0 means continuous time transfer function.
    sampling_time > 0 means discrete time transfer function.

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Complex variable of the Laplace or z transform used by the
        polynomials of the transfer function.
    sampling_time : Symbol, Number, optional
        Default is 0.
        Time interval between two consecutive sampling instants.
        If sampling_time == 0, it is a continuous time transfer function,
        else it is a discrete time transfer function.

    Examples
    ========

    >>> from sympy.abc import s, z
    >>> from sympy.physics.control.lti import create_transfer_function
    >>> num = s + 5
    >>> den = 3*s**2 + 2*s + 1
    >>> tf = create_transfer_function(num, den, s)
    >>> tf
    TransferFunction(s + 5, 3*s**2 + 2*s + 1, s)
    >>> num = z
    >>> den = z + 1
    >>> dtf = create_transfer_function(num, den, z, 0.1)
    >>> dtf
    DiscreteTransferFunction(z, z + 1, z, 0.1)

    See Also
    ========

    TransferFunction, DiscreteTransferFunction

    """
    if sampling_time == 0:
        return TransferFunction(num, den, var)
    return DiscreteTransferFunction(num, den, var, sampling_time)

class TransferFunctionBase(SISOLinearTimeInvariant, ABC):
    r"""
    Base class for transfer tunction objects.
    This class is not meant to be used directly.

    Explanation
    ===========

    LTI systems can be described by linear differential equations (for
    continuous-time systems) or linear difference equations (for discrete-time
    systems).
    Mathematical transforms are employed to convert these equations into simpler
    algebraic forms in a complex variable domain.
    For continuous-time systems, the Laplace transform (using complex variable
    :math:`s`) is used. For discrete-time systems, the z-transform (using
    complex variable :math:`z`) is used.

    We will call the generic transformation used as :math:`\mathcal{T\{\cdot\}}`
    and the variable of the transform as :math:`p`.

    Consider a continuous-time LTI system described by a linear ordinary
    differential equation:

    .. math::
        b_{m}y^{\left(m\right)}+b_{m-1}y^{\left(m-1\right)}+\dots+b_{1}
        y^{\left(1\right)}+b_{0}y=
        a_{n}x^{\left(n\right)}+a_{n-1}x^{\left(n-1\right)}+\dots+a_{1}
        x^{\left(1\right)}+a_{0}x

    or a discrete-time LTI system described by a linear difference equation:

    .. math::
        b_{m}y[k-m]+b_{m-1}y[k-m+1]+\dots+b_{1}y[k-1]+b_{0}y[k]=
        a_{n}x[k-n]+a_{n-1}x[k-n+1]+\dots+a_{1}x[k-1]+a_{0}x[k]

    Here, :math:`x` is the input signal and :math:`y` is the output signal.

    It is not feasible to analyse the properties of such systems in their native
    form therefore, we use mathematical tools like Laplace transform or
    z-transform to get a better perspective.
    Taking the transform :math:`\mathcal{T}\{\cdot\}` of both the sides in the
    equation (at zero initial conditions), we get:

    .. math::
        \mathcal{T}[b_{m}y^{\left(m\right)}+b_{m-1}y^{\left(m-1\right)}+\dots+
        b_{1}y^{\left(1\right)}+b_{0}y]=
        \mathcal{T}[a_{n}x^{\left(n\right)}+a_{n-1}x^{\left(n-1\right)}+\dots+
        a_{1}x^{\left(1\right)}+a_{0}x]

    and using its linearity and differentiation properties
    (:math:`\mathcal{L}\{f^{(k)}(t)\} = s^k F(s)` under zero initial conditions),
    we get:

    .. math::
        b_{m}p^{m}\mathcal{T}[y]+\dots+b_{1}p\mathcal{T}[y]+b_{0}\mathcal{T}[y]=
        a_{n}p^{n}\mathcal{T}[x]+\dots+a_{1}p\mathcal{T}[x]+a_{0}\mathcal{T}[x]

    Note that the zero initial conditions assumption, mentioned above,
    is very important and cannot be ignored otherwise the dynamical system
    cannot be considered time-independent and the simplified equation above
    cannot be reached.

    The numerator of the transfer function is, therefore, the transform of the
    output signal and similarly, the denominator of the transfer function is
    the transform of the input signal.
    It is also a convention to denote the input and output signal's transform
    with capital alphabets like shown below:

    .. math::
        H(p) = \frac{Y(p)}{X(p)} = \frac{ \mathcal{T}\left\{y(t)\right\}}
        {\mathcal{T}\left\{x(t)\right\}}

    Transfer functions are sometimes also referred to as the transform of the
    system's impulse response. Transfer function, :math:`H`, is represented as a
    rational function like,

    .. math::
        H(p) =\ \frac{a_{n}p^{n}+a_{n-1}p^{n-1}+\dots+a_{1}p+a_{0}}{b_{m}p^{m}+
        b_{m-1}p^{m-1}+\dots+b_{1}p+b_{0}}

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Complex variable of the Laplace transform used by the
        polynomials of the transfer function.
    *args, **kwargs:
        Additional arguments and keyword arguments that are passed to the
        parent class such as sampling time for discrete-time systems.

    Raises
    ======

    TypeError
        When ``var`` is not a Symbol or when ``num`` or ``den`` is not a
        number or a polynomial.
    ValueError
        When ``den`` is zero.

    See Also
    ========

    TransferFunction, DiscreteTransferFunction, Feedback, Series, Parallel

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Transfer_function
    .. [2] https://en.wikipedia.org/wiki/Laplace_transform
    .. [3] https://en.wikipedia.org/wiki/Z-transform

    """
    def __new__(cls, num, den, var, *args, **kwargs):
        if cls is TransferFunctionBase:
            raise NotImplementedError(
                """
                The TransferFunctionBase class is not meant to be used directly.
                """)
        num, den = _sympify(num), _sympify(den)

        if not isinstance(var, Symbol):
            raise TypeError("Variable input must be a Symbol.")

        if den == 0:
            raise ValueError("TransferFunction cannot have a zero denominator.")

        accepted_istances = (Expr, TransferFunctionBase, Series, Parallel)
        num_accepted = ((isinstance(num, accepted_istances) and num.has(Symbol))
                    or num.is_number)
        den_accepted = ((isinstance(den, accepted_istances) and den.has(Symbol))
                    or den.is_number)

        if not num_accepted or not den_accepted:
            raise TypeError("""Unsupported type for numerator or denominator of
                            TransferFunction.""")

        obj = super(TransferFunctionBase, cls).__new__(cls, num, den, var,
                                                        *args, **kwargs)
        obj.is_StateSpace_object = False
        return obj


    @classmethod
    def from_rational_expression(cls, expr, var=None, *args, **kwargs):
        r"""
        Creates a new transfer function efficiently from a rational
        expression.

        Parameters
        ==========

        expr : Expr, Number
            The rational expression representing the transfer function.
        var : Symbol, optional
            Complex variable used by the polynomials of the transfer function.
        *args, **kwargs :
            Additional positional and keyword arguments to be passed to the
            constructor of the :class:`~.TransferFunctionBase`, such as
            sampling_time.

        Raises
        ======

        ValueError
            When ``expr`` is of type ``Number`` and optional parameter ``var``
            is not passed.

            When ``expr`` has more than one variables and an optional parameter
            ``var`` is not passed.
        ZeroDivisionError
            When denominator of ``expr`` is zero or it has ``ComplexInfinity``
            in its numerator.

        Examples
        ========

        >>> from sympy.abc import s, p, a, z
        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> expr1 = (s + 5)/(3*s**2 + 2*s + 1)
        >>> tf1 = TransferFunction.from_rational_expression(expr1)
        >>> tf1
        TransferFunction(s + 5, 3*s**2 + 2*s + 1, s)
        >>> expr2 = (a*p**3 - a*p**2 + s*p)/(p + a**2)  # Expr with more than one variables
        >>> tf2 = TransferFunction.from_rational_expression(expr2, p)
        >>> tf2
        TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)
        >>> expr3 = (z + 1)/(z**2 + 2*z + 1)  # Discrete time transfer function
        >>> dtf = DiscreteTransferFunction.from_rational_expression(expr3, z, sampling_time=0.1)
        >>> dtf
        DiscreteTransferFunction(z + 1, z**2 + 2*z + 1, z, 0.1)

        In case of conflict between two or more variables in a expression, SymPy will
        raise a ``ValueError``, if ``var`` is not passed by the user.

        >>> tf = TransferFunction.from_rational_expression((a + a*s)/(s**2 + s + 1))
        Traceback (most recent call last):
        ...
        ValueError: Conflicting values found for positional argument `var` ({a, s}). Specify it manually.

        This can be corrected by specifying the ``var`` parameter manually.

        >>> tf = TransferFunction.from_rational_expression((a + a*s)/(s**2 + s + 1), s)
        >>> tf
        TransferFunction(a*s + a, s**2 + s + 1, s)

        ``var`` also need to be specified when ``expr`` is a ``Number``

        >>> tf3 = TransferFunction.from_rational_expression(10, s)
        >>> tf3
        TransferFunction(10, 1, s)

        """
        expr = _sympify(expr)
        if var is None:
            _free_symbols = expr.free_symbols
            _len_free_symbols = len(_free_symbols)
            if _len_free_symbols == 1:
                var = list(_free_symbols)[0]
            elif _len_free_symbols == 0:
                raise ValueError(filldedent("""
                    Positional argument `var` not found in the
                    TransferFunction defined. Specify it manually."""))
            else:
                raise ValueError(filldedent("""
                    Conflicting values found for positional argument `var` ({}).
                    Specify it manually.""".format(_free_symbols)))

        _num, _den = expr.as_numer_denom()
        if _den == 0 or _num.has(S.ComplexInfinity):
            raise ZeroDivisionError("TransferFunction cannot have a zero denominator.")
        return cls(_num, _den, var, *args, **kwargs)

    @classmethod
    def from_coeff_lists(cls, num_list, den_list, var, *args, **kwargs):
        r"""
        Creates a new transfer function efficiently from a list of coefficients.

        Parameters
        ==========

        num_list : Sequence
            Sequence comprising of numerator coefficients.
        den_list : Sequence
            Sequence comprising of denominator coefficients.
        var : Symbol
            Complex variable used by the polynomials of the transfer function.
        *args, **kwargs :
            Additional positional and keyword arguments to be passed to the
            constructor of the :class:`~.TransferFunctionBase`, such as
            sampling_time.

        Raises
        ======

        ZeroDivisionError
            When the constructed denominator is zero.

        Examples
        ========

        >>> from sympy.abc import s, p, z
        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> num = [1, 0, 2]
        >>> den = [3, 2, 2, 1]
        >>> tf = TransferFunction.from_coeff_lists(num, den, s)
        >>> tf
        TransferFunction(s**2 + 2, 3*s**3 + 2*s**2 + 2*s + 1, s)
        >>> #Create a Transfer Function with more than one variable
        >>> tf1 = TransferFunction.from_coeff_lists([p, 1], [2*p, 0, 4], s)
        >>> tf1
        TransferFunction(p*s + 1, 2*p*s**2 + 4, s)
        >>> dtf = DiscreteTransferFunction.from_coeff_lists([2, 1, -3], [1, 4, 3, 2], z, sampling_time=0.1)
        >>> dtf
        DiscreteTransferFunction(2*z**2 + z - 3, z**3 + 4*z**2 + 3*z + 2, z, 0.1)

        """
        num_list = num_list[::-1]
        den_list = den_list[::-1]
        num_var_powers = [var**i for i in range(len(num_list))]
        den_var_powers = [var**i for i in range(len(den_list))]

        _num = sum(coeff * var_power for coeff, var_power in zip(num_list, num_var_powers))
        _den = sum(coeff * var_power for coeff, var_power in zip(den_list, den_var_powers))

        if _den == 0:
            raise ZeroDivisionError("TransferFunction cannot have a zero denominator.")

        return cls(_num, _den, var, *args, **kwargs)

    @classmethod
    def from_zpk(cls, zeros, poles, gain, var, *args, **kwargs):
        r"""
        Creates a new transfer function from given zeros, poles and gain.

        Parameters
        ==========

        zeros : Sequence
            Sequence comprising of zeros of transfer function.
        poles : Sequence
            Sequence comprising of poles of transfer function.
        gain : Number, Symbol, Expression
            A scalar value specifying gain of the model.
        var : Symbol
            Complex variable used by the polynomials of the transfer function.
        *args, **kwargs :
            Additional positional and keyword arguments to be passed to the
            constructor of the :class:`~.TransferFunctionBase`, such as
            sampling_time.

        Examples
        ========

        >>> from sympy.abc import s, p, k
        >>> from sympy.physics.control.lti import TransferFunction
        >>> zeros = [1, 2, 3]
        >>> poles = [6, 5, 4]
        >>> gain = 7
        >>> tf = TransferFunction.from_zpk(zeros, poles, gain, s)
        >>> tf
        TransferFunction(7*(s - 3)*(s - 2)*(s - 1), (s - 6)*(s - 5)*(s - 4), s)
        >>> #Create a Transfer Function with variable poles and zeros
        >>> tf1 = TransferFunction.from_zpk([p, k], [p + k, p - k], 2, s)
        >>> tf1
        TransferFunction(2*(-k + s)*(-p + s), (-k - p + s)*(k - p + s), s)
        >>> #Complex poles or zeros are acceptable
        >>> tf2 = TransferFunction.from_zpk([0], [1-1j, 1+1j, 2], -2, s)
        >>> tf2
        TransferFunction(-2*s, (s - 2)*(s - 1.0 - 1.0*I)*(s - 1.0 + 1.0*I), s)

        """
        num_poly = 1
        den_poly = 1
        for zero in zeros:
            num_poly *= var - zero
        for pole in poles:
            den_poly *= var - pole

        return cls(gain*num_poly, den_poly, var, *args, **kwargs)

    @property
    def num(self):
        """
        Returns the numerator polynomial of the transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(s**2 + p*s + 3, s - 4, s)
        >>> G1.num
        p*s + s**2 + 3
        >>> G2 = TransferFunction((p + 5)*(p - 3), (p - 3)*(p + 1), p)
        >>> G2.num
        (p - 3)*(p + 5)

        """
        return self.args[0]

    @property
    def den(self):
        """
        Returns the denominator polynomial of the transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(s + 4, p**3 - 2*p + 4, s)
        >>> G1.den
        p**3 - 2*p + 4
        >>> G2 = TransferFunction(3, 4, s)
        >>> G2.den
        4

        """
        return self.args[1]

    @property
    def var(self):
        """
        Returns the complex variable used by the polynomials of the transfer
        function.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G1.var
        p
        >>> G2 = TransferFunction(0, s - 5, s)
        >>> G2.var
        s

        """
        return self.args[2]

    def _eval_subs(self, old, new):
        if old == self.var:
            return self

        arg_num = self.num.subs(old, new)
        arg_den = self.den.subs(old, new)

        argnew = create_transfer_function(arg_num, arg_den, self.var, self.sampling_time)

        return argnew

    def _eval_evalf(self, prec):
        return create_transfer_function(
            self.num._eval_evalf(prec),
            self.den._eval_evalf(prec),
            self.var, self.sampling_time)

    def _eval_simplify(self, **kwargs):
        tf = cancel(Mul(self.num, 1/self.den, evaluate=False),
                    expand=False).as_numer_denom()
        num_, den_ = tf[0], tf[1]

        return create_transfer_function(num_, den_, self.var, self.sampling_time)

    def expand(self):
        """
        Returns the transfer function with numerator and denominator
        in expanded form.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> G1 = TransferFunction((a - s)**2, (s**2 + a)**2, s)
        >>> G1.expand()
        TransferFunction(a**2 - 2*a*s + s**2, a**2 + 2*a*s**2 + s**4, s)
        >>> G2 = DiscreteTransferFunction((p + 3*b)*(p - b), (p - b)*(p + 2*b), p, 12)
        >>> G2.expand()
        DiscreteTransferFunction(-3*b**2 + 2*b*p + p**2, -2*b**2 + b*p + p**2, p, 12)

        """
        return create_transfer_function(expand(self.num), expand(self.den), self.var,
                              self.sampling_time)

    @abstractmethod
    def dc_gain(self):
        """
        Computes the gain of the response as the frequency approaches zero.

        The DC gain is infinite for systems with pure integrators.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b, z
        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> tf1 = TransferFunction(s + 3, s**2 - 9, s)
        >>> tf1.dc_gain()
        -1/3
        >>> tf2 = TransferFunction(p**2, p - 3 + p**3, p)
        >>> tf2.dc_gain()
        0
        >>> tf3 = TransferFunction(a*p**2 - b, s + b, s)
        >>> tf3.dc_gain()
        (a*p**2 - b)/b
        >>> tf4 = TransferFunction(1, s, s)
        >>> tf4.dc_gain()
        oo
        >>> dtf1 = DiscreteTransferFunction(z, z - 1, z, 0.1)
        >>> dtf1.dc_gain()
        oo

        """
        pass

    def poles(self):
        """
        Returns the poles of a transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
        >>> tf1.poles()
        [-5, 1]
        >>> tf2 = TransferFunction((1 - s)**2, (s**2 + 1)**2, s)
        >>> tf2.poles()
        [I, I, -I, -I]
        >>> tf3 = TransferFunction(s**2, a*s + p, s)
        >>> tf3.poles()
        [-p/a]

        """
        return _roots(Poly(self.den, self.var), self.var)

    def zeros(self):
        """
        Returns the zeros of a transfer function.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
        >>> tf1.zeros()
        [-3, 1]
        >>> tf2 = TransferFunction((1 - s)**2, (s**2 + 1)**2, s)
        >>> tf2.zeros()
        [1, 1]
        >>> tf3 = TransferFunction(s**2, a*s + p, s)
        >>> tf3.zeros()
        [0, 0]

        """
        return _roots(Poly(self.num, self.var), self.var)

    def eval_frequency(self, other):
        """
        Returns the system response at any point in the real or complex plane.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy.physics.control.lti import TransferFunction
        >>> from sympy import I
        >>> tf1 = TransferFunction(1, s**2 + 2*s + 1, s)
        >>> omega = 0.1
        >>> tf1.eval_frequency(I*omega)
        0.970493088912852 - 0.196059209881384*I
        >>> tf2 = TransferFunction(s**2, a*s + p, s)
        >>> tf2.eval_frequency(2)
        4/(2*a + p)
        >>> tf2.eval_frequency(I*2)
        -4/(2*I*a + p)
        """
        arg_num = self.num.subs(self.var, other)
        arg_den = self.den.subs(self.var, other)
        return Mul(arg_num, S.One / arg_den).expand()

    def is_stable(self, cancel_poles_zeros=False):
        """
        Returns True if the transfer function is asymptotically stable;
        else False.

        This would not check the marginal or conditional stability
        of the system.

        Note: Also with cancel_poles_zeros = True, there could be unaccounted
        pole-zero cancellations.

        Parameters
        ==========

        cancel_poles_zeros : Boolean
            If True, cancels common factors between numerator and denominator
            before checking stability.

        Examples
        ========

        >>> from sympy.abc import s, p, a
        >>> from sympy import symbols
        >>> from sympy.physics.control.lti import TransferFunction
        >>> q, r = symbols('q, r', negative=True)
        >>> tf1 = TransferFunction((1 - s)**2, (s + 1)**2, s)
        >>> tf1.is_stable()
        True
        >>> tf2 = TransferFunction((1 - p)**2, (s**2 + 1)**2, s)
        >>> tf2.is_stable()
        False
        >>> tf3 = TransferFunction(4, q*s - r, s)
        >>> tf3.is_stable()
        False
        >>> tf4 = TransferFunction(p + 1, a*p - s**2, p)
        >>> tf4.is_stable() is None   # Not enough info about the symbols to determine stability
        True
        >>> tf5 = TransferFunction((s+1)*(s-1), (s-1)*(s+2)*(s+4), s)
        >>> tf5.is_stable()
        False
        >>> tf5.is_stable(cancel_poles_zeros = True)
        True

        """
        tf = self.to_standard_form(cancel_poles_zeros)

        conditions = tf.get_asymptotic_stability_conditions(
            cancel_poles_zeros = False
        )

        try:
            output = reduce_inequalities(conditions)
        except NotImplementedError:
            # If there are more than one symbols,
            # reduce_inequalities could fail
            return None

        if output in (true, false):
            return bool(output)

        return None

    def to_standard_form(self, cancel_poles_zeros=False):
        r"""
        Return the transfer function in its standard form.

        Standard form:

        .. math::
            \frac{a_n s^n + a_{n-1} s^{n-1} + \cdots + a_1 s + a_0}
            {b_m s^m + b_{m-1} s^{m-1} + \cdots + b_1 s + b_0}

        Note: Also with cancel_poles_zeros = True, there could be unaccounted
        pole-zero cancellations.

        Examples
        ========
        >>> from sympy import symbols
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> s,k = symbols('s k')
        >>> tf1 = TransferFunction(-2*s + 12, s**3 + 2*s**2 + 100*s, s)
        >>> tf2 = TransferFunction(1, k, s)
        >>> feedback = Feedback(tf1, tf2).doit()
        >>> feedback
        TransferFunction(k*(12 - 2*s)*(s**3 + 2*s**2 + 100*s), (s**3 + 2*s**2 + 100*s)*(k*(s**3 + 2*s**2 + 100*s) - 2*s + 12), s)
        >>> feedback.to_standard_form(cancel_poles_zeros=True)
        TransferFunction(-2*k*s + 12*k, k*s**3 + 2*k*s**2 + s*(100*k - 2) + 12, s)

        """
        tf = self.expand()

        num = tf.num
        den = tf.den
        if cancel_poles_zeros:
            num, den = cancel(tf.num / tf.den).as_numer_denom()

        return create_transfer_function(num.collect(self.var), den.collect(self.var),
                      self.var, self.sampling_time)

    @abstractmethod
    def get_asymptotic_stability_conditions(
        self, cancel_poles_zeros=False, fast=False
    ) -> list[Boolean]:
        """
        Returns the asymptotic stability conditions for the transfer function.

        This is a convenient shorthand, based on the system type (continuous or
        discrete), for:

        - ``[c > 0 for c in Poly(sys.den, sys.var).hurwitz_conditions()]``
          which gives conditions for stability such that the poles of the
          transfer function are in the left half of the complex plane.
        - ``[c > 0 for c in Poly(sys.den, sys.var).schur_conditions()]``
          which gives conditions for stability such that the poles of the
          transfer function lie inside the unit circle.

        For some systems that are larger or have more complicated
        expressions it might be useful to set ``fast`` to ``True``.
        The algorithm will use the EXRAW domain which will quickly generate
        large unsimplified expressions that are mostly only suitable for use
        with lambdify.

        Notes
        =====

        - Also with ``cancel_poles_zeros = True``, the set generated by the
          inequalities may be a subset of the real asymptotic stability
          conditions set due to potentially unaccounted pole-zero cancellations.
        - This method assumes that the leading coefficient is non-zero.
          In the opposite case, additional verification is required.
        - In the discrete case, using ``fast = True`` may lead to significant
          precision issues.

        Parameters
        ==========

        cancel_poles_zeros : Boolean
            If True, cancels common factors between numerator and denominator
            before checking stability.
        fast : Boolean
            If True, uses the EXRAW domain to quickly generate large
            unsimplified expressions. If False, uses the default domain
            which is suitable for symbolic manipulation.

        Examples
        ========

        >>> from sympy import symbols, solve, reduce_inequalities
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, Feedback

        >>> b1, b2, b3, b4 = symbols('b_1 b_2 b_3 b_4')
        >>> p1 = b1*s**3 + b2*s**2 + b3*s + b4
        >>> tf1 = TransferFunction(1, p1, s)
        >>> tf1.get_asymptotic_stability_conditions()
        [b_1*b_2 > 0, -b_1*b_4 + b_2*b_3 > 0, b_1*b_4 > 0]

        >>> p2 = s**4 + 3*s**3 + 6*s**2 + 12*s + 8
        >>> solve(p2)
        [-2, -1, -2*I, 2*I]
        >>> tf2 = TransferFunction(1, p2, s)
        >>> tf2.get_asymptotic_stability_conditions()
        [False]

        >>> p3 = s**4 + 17*s**3 + 137/2*s**2 + 213/2*s + 54
        >>> solve(p3)
        [-12.0, -1.0, -2.0 - 0.707106781186548*I, -2.0 + 0.707106781186548*I]
        >>> tf3 = TransferFunction(1, p3, s)
        >>> tf3.get_asymptotic_stability_conditions()
        [True, True, True, True]

        >>> k = symbols('k')
        >>> tf4 = TransferFunction(-20*s + 20, s**3 + 2*s**2 + 100*s, s)
        >>> tf5 = TransferFunction(1, k, s)
        >>> feedback = Feedback(tf4, tf5).doit()
        >>> ineq = feedback.get_asymptotic_stability_conditions(cancel_poles_zeros = True)
        >>> ineq
        [2*k**2 > 0, 200*k**2 - 60*k > 0, 20*k > 0]
        >>> reduce_inequalities(ineq)
        (3/10 < k) & (k < oo)

        """
        pass

    @_compatibility_decorator
    def __add__(self, other):
        if hasattr(other, "is_StateSpace_object") and other.is_StateSpace_object:
            return Parallel(self, other)
        if isinstance(other, (self.__class__, Series, Feedback)):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace domain or z-domain."""))
            return Parallel(self, other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace domain or z-domain."""))
            arg_list = list(other.args)
            return Parallel(self, *arg_list)
        else:
            raise ValueError("{} cannot be added with {}.".
                format(type(self), type(other)))

    @_compatibility_decorator
    def __radd__(self, other):
        return self + other

    @_compatibility_decorator
    def __sub__(self, other):
        if hasattr(other, "is_StateSpace_object") and other.is_StateSpace_object:
            return Parallel(self, -other)
        elif isinstance(other, (self.__class__, Series)):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace domain or z-domain."""))
            return Parallel(self, -other)
        elif isinstance(other, Parallel):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace transform."""))
            arg_list = [-i for i in list(other.args)]
            return Parallel(self, *arg_list)
        else:
            raise ValueError("{} cannot be subtracted from a {}."
                .format(type(self), type(other)))

    @_compatibility_decorator
    def __rsub__(self, other):
        return -self + other

    @_compatibility_decorator
    def __mul__(self, other):
        if hasattr(other, "is_StateSpace_object") and other.is_StateSpace_object:
            return Series(self, other)
        elif isinstance(other, (self.__class__, Parallel, Feedback)):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex variable
                    of the Laplace domain or z-domain."""))
            return Series(self, other)
        elif isinstance(other, Series):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex variable
                    of the Laplace domain or z-domain."""))
            arg_list = list(other.args)
            return Series(self, *arg_list)
        else:
            raise ValueError("{} cannot be multiplied with {}."
                .format(type(self), type(other)))

    __rmul__ = __mul__

    @_compatibility_decorator
    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace domain or z-domain."""))
            return Series(self, create_transfer_function(other.den, other.num,
                                       self.var, self.sampling_time))
        elif (isinstance(other, Parallel) and len(other.args
                ) == 2 and isinstance(other.args[0], self.__class__)
            and isinstance(other.args[1], (Series, self.__class__))):

            if not self.var == other.var:
                raise ValueError(filldedent("""
                    Both TransferFunction and Parallel should use the
                    same complex variable of the Laplace domain or z-domain."""))
            if other.args[1] == self:
                # plant and controller with unit feedback.
                return Feedback(self, other.args[0])
            other_arg_list = list(other.args[1].args) if isinstance(
                other.args[1], Series) else other.args[1]
            if other_arg_list == other.args[1]:
                return Feedback(self, other_arg_list)
            elif self in other_arg_list:
                other_arg_list.remove(self)
            else:
                return Feedback(self, Series(*other_arg_list))

            if len(other_arg_list) == 1:
                return Feedback(self, *other_arg_list)
            else:
                return Feedback(self, Series(*other_arg_list))
        else:
            raise ValueError("{} cannot be divided by {}.".
                format(type(self), type(other)))

    __rtruediv__ = __truediv__

    def __pow__(self, p):
        p = sympify(p)
        if not p.is_Integer:
            raise ValueError("Exponent must be an integer.")
        if p is S.Zero:
            return create_transfer_function(1, 1, self.var, self.sampling_time)
        elif p > 0:
            num_, den_ = self.num**p, self.den**p
        else:
            p = abs(p)
            num_, den_ = self.den**p, self.num**p

        return create_transfer_function(num_, den_, self.var, self.sampling_time)

    def __neg__(self):
        return create_transfer_function(-self.num, self.den, self.var, self.sampling_time)

    def _StateSpace_matrices_equivalent(self):
        """
        Returns the state-space matrices equivalent
        to the transfer function.
        The state-space matrices are returned in the form of (A, B, C, D) in
        the controllable canonical form.
        """
        if not self.is_proper:
            raise ValueError("Transfer Function must be proper.")

        num_poly = Poly(self.num, self.var)
        den_poly = Poly(self.den, self.var)
        n = den_poly.degree()

        num_coeffs = num_poly.all_coeffs()
        den_coeffs = den_poly.all_coeffs()
        diff = n - num_poly.degree()
        num_coeffs = [0]*diff + num_coeffs

        a = den_coeffs[1:]
        a_mat = Matrix([[(-1)*coefficient/den_coeffs[0] for coefficient in reversed(a)]])
        vert = zeros(n-1, 1)
        mat = eye(n-1)
        A = vert.row_join(mat)
        A = A.col_join(a_mat)

        B = zeros(n, 1)
        B[n-1] = 1

        i = n
        C = []
        while(i > 0):
            C.append(num_coeffs[i] - den_coeffs[i]*num_coeffs[0])
            i -= 1
        C = Matrix([C])

        D = Matrix([num_coeffs[0]])

        return A, B, C, D

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial is less than
        or equal to degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf1.is_proper
        False
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*p + 2, p)
        >>> tf2.is_proper
        True

        """
        return degree(self.num, self.var) <= degree(self.den, self.var)

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial is strictly less
        than degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf1.is_strictly_proper
        False
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf2.is_strictly_proper
        True

        """
        return degree(self.num, self.var) < degree(self.den, self.var)

    @property
    def is_biproper(self):
        """
        Returns True if degree of the numerator polynomial is equal to
        degree of the denominator polynomial, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf1.is_biproper
        True
        >>> tf2 = TransferFunction(p**2, p + a, p)
        >>> tf2.is_biproper
        False

        """
        return degree(self.num, self.var) == degree(self.den, self.var)

    def to_expr(self):
        """
        Converts a :obj:`~.TransferFunction` object to SymPy
        :obj:`sympy.core.expr.Expr`.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction
        >>> from sympy import Expr
        >>> tf1 = TransferFunction(s, a*s**2 + 1, s)
        >>> tf1.to_expr()
        s/(a*s**2 + 1)
        >>> isinstance(_, Expr)
        True
        >>> tf2 = TransferFunction(1, (p + 3*b)*(b - p), p)
        >>> tf2.to_expr()
        1/((b - p)*(3*b + p))
        >>> tf3 = TransferFunction((s - 2)*(s - 3), (s - 1)*(s - 2)*(s - 3), s)
        >>> tf3.to_expr()
        ((s - 3)*(s - 2))/(((s - 3)*(s - 2)*(s - 1)))

        """

        if self.num != 1:
            return Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)
        else:
            return Pow(self.den, -1, evaluate=False)

class TransferFunction(TransferFunctionBase):
    r"""
    A class for representing LTI (Linear, time-invariant) systems that can be
    strictly described by ratio of polynomials in the Laplace transform complex variable. The arguments
    are ``num``, ``den``, and ``var``, where ``num`` and ``den`` are numerator and
    denominator polynomials of the ``TransferFunction`` respectively, and the third argument is
    a complex variable of the Laplace transform used by these polynomials of the transfer function.
    ``num`` and ``den`` can be either polynomials or numbers, whereas ``var``
    has to be a :py:class:`~.Symbol`.

    See :class:`TransferFunctionBase` for more information.

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Complex variable of the Laplace transform used by the
        polynomials of the transfer function.

    Raises
    ======

    TypeError
        When ``var`` is not a Symbol or when ``num`` or ``den`` is not a
        number or a polynomial.
    ValueError
        When ``den`` is zero.

    Examples
    ========

    >>> from sympy.abc import s, p, a
    >>> from sympy.physics.control.lti import TransferFunction
    >>> tf1 = TransferFunction(s + a, s**2 + s + 1, s)
    >>> tf1
    TransferFunction(a + s, s**2 + s + 1, s)
    >>> tf1.num
    a + s
    >>> tf1.den
    s**2 + s + 1
    >>> tf1.var
    s
    >>> tf1.args
    (a + s, s**2 + s + 1, s)

    Any complex variable can be used for ``var``.

    >>> tf2 = TransferFunction(a*p**3 - a*p**2 + s*p, p + a**2, p)
    >>> tf2
    TransferFunction(a*p**3 - a*p**2 + p*s, a**2 + p, p)
    >>> tf3 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    >>> tf3
    TransferFunction((p - 1)*(p + 3), (p - 1)*(p + 5), p)

    To negate a transfer function the ``-`` operator can be prepended:

    >>> tf4 = TransferFunction(-a + s, p**2 + s, p)
    >>> -tf4
    TransferFunction(a - s, p**2 + s, p)
    >>> tf5 = TransferFunction(s**4 - 2*s**3 + 5*s + 4, s + 4, s)
    >>> -tf5
    TransferFunction(-s**4 + 2*s**3 - 5*s - 4, s + 4, s)

    You can use a float or an integer (or other constants) as numerator and denominator:

    >>> tf6 = TransferFunction(1/2, 4, s)
    >>> tf6.num
    0.500000000000000
    >>> tf6.den
    4
    >>> tf6.var
    s
    >>> tf6.args
    (0.5, 4, s)

    You can take the integer power of a transfer function using the ``**`` operator:

    >>> tf7 = TransferFunction(s + a, s - a, s)
    >>> tf7**3
    TransferFunction((a + s)**3, (-a + s)**3, s)
    >>> tf7**0
    TransferFunction(1, 1, s)
    >>> tf8 = TransferFunction(p + 4, p - 3, p)
    >>> tf8**-1
    TransferFunction(p - 3, p + 4, p)

    Addition, subtraction, and multiplication of transfer functions can form
    unevaluated ``Series`` or ``Parallel`` objects.

    >>> tf9 = TransferFunction(s + 1, s**2 + s + 1, s)
    >>> tf10 = TransferFunction(s - p, s + 3, s)
    >>> tf11 = TransferFunction(4*s**2 + 2*s - 4, s - 1, s)
    >>> tf12 = TransferFunction(1 - s, s**2 + 4, s)
    >>> tf9 + tf10
    Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf10 - tf11
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-4*s**2 - 2*s + 4, s - 1, s))
    >>> tf9 * tf10
    Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-p + s, s + 3, s))
    >>> tf10 - (tf9 + tf12)
    Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(-s - 1, s**2 + s + 1, s), TransferFunction(s - 1, s**2 + 4, s))
    >>> tf10 - (tf9 * tf12)
    Parallel(TransferFunction(-p + s, s + 3, s), Series(TransferFunction(-1, 1, s), TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> tf11 * tf10 * tf9
    Series(TransferFunction(4*s**2 + 2*s - 4, s - 1, s), TransferFunction(-p + s, s + 3, s), TransferFunction(s + 1, s**2 + s + 1, s))
    >>> tf9 * tf11 + tf10 * tf12
    Parallel(Series(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)), Series(TransferFunction(-p + s, s + 3, s), TransferFunction(1 - s, s**2 + 4, s)))
    >>> (tf9 + tf12) * (tf10 + tf11)
    Series(Parallel(TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(1 - s, s**2 + 4, s)), Parallel(TransferFunction(-p + s, s + 3, s), TransferFunction(4*s**2 + 2*s - 4, s - 1, s)))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function using ``.doit()`` method or by ``.rewrite(TransferFunction)``.

    >>> ((tf9 + tf10) * tf12).doit()
    TransferFunction((1 - s)*((-p + s)*(s**2 + s + 1) + (s + 1)*(s + 3)), (s + 3)*(s**2 + 4)*(s**2 + s + 1), s)
    >>> (tf9 * tf10 - tf11 * tf12).rewrite(TransferFunction)
    TransferFunction(-(1 - s)*(s + 3)*(s**2 + s + 1)*(4*s**2 + 2*s - 4) + (-p + s)*(s - 1)*(s + 1)*(s**2 + 4), (s - 1)*(s + 3)*(s**2 + 4)*(s**2 + s + 1), s)

    See Also
    ========

    TransferFunctionBase, DiscreteTransferFunction, Feedback, Series, Parallel

    """
    def __new__(cls, num, den, var):
        return super(TransferFunction, cls).__new__(cls, num, den, var)

    @classmethod
    def from_rational_expression(cls, expr, var=None):
        r"""
        See :func:`TransferFunctionBase.from_rational_expression`.

        """
        return super().from_rational_expression(expr, var)

    @classmethod
    def from_coeff_lists(cls, num_list, den_list, var):
        r"""
        See :func:`TransferFunctionBase.from_coeff_lists`.

        """
        return super().from_coeff_lists(num_list, den_list, var)

    @classmethod
    def from_zpk(cls, zeros, poles, gain, var):
        r"""
        See :func:`TransferFunctionBase.from_zpk`.

        """
        return super().from_zpk(zeros, poles, gain, var)

    def dc_gain(self):
        r"""
        See :func:`TransferFunctionBase.dc_gain`.

        """
        m = Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)
        return limit(m, self.var, 0)

    def get_asymptotic_stability_conditions(
        self, cancel_poles_zeros=False, fast=False
    ) -> list[Boolean]:
        r"""
        See :func:`TransferFunctionBase.get_asymptotic_stability_conditions`.

        """
        standard_form = self.to_standard_form(cancel_poles_zeros)

        domain = EXRAW if fast else None

        p = Poly(standard_form.den, self.var, domain = domain)

        return [c > 0 for c in p.hurwitz_conditions()]

    def _eval_rewrite_as_StateSpace(self, *args):
        """
        Returns the equivalent space model of the transfer function model.
        The state space model will be returned in the controllable canonical
        form.

        Unlike the space state to transfer function model conversion, the
        transfer function to state space model conversion is not unique.
        There can be multiple state space representations of a given transfer
        function model.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control import TransferFunction, StateSpace
        >>> tf = TransferFunction(s**2 + 1, s**3 + 2*s + 10, s)
        >>> tf.rewrite(StateSpace)
        StateSpace(Matrix([
        [  0,  1, 0],
        [  0,  0, 1],
        [-10, -2, 0]]), Matrix([
        [0],
        [0],
        [1]]), Matrix([[1, 0, 1]]), Matrix([[0]]))

        """
        A, B, C, D = self._StateSpace_matrices_equivalent()
        return StateSpace(A, B, C, D)

    def _eval_rewrite_as_DiscreteStateSpace(self, *args):
        raise TypeError("""
            The continuous transfer function model cannot be rewritten as a
            discrete-time state space model.
            """)

    @property
    def sampling_time(self):
        """The sampling time of the transfer function is zero."""
        return S.Zero

    _is_continuous = True

class DiscreteTransferFunction(TransferFunctionBase):
    r"""
    A class for representing LTI (Linear, time-invariant) systems that can be
    strictly described by ratio of polynomials in the z-transform complex
    variable.
    The arguments are ``num``, ``den``, ``var``, and ``sampling_time``,
    where ``num`` and ``den`` are numerator and denominator polynomials of the
    ``DiscreteTransferFunction`` respectively, the third argument is
    a complex variable of the z-transform used by these polynomials of the
    transfer function, and the fourth represents the time interval between two
    consecutive sampling instants. If not specified, it defaults to 1 and, if
    it's set to 0, an instance of :class:`TransferFunction` is returned.
    ``num`` and ``den`` can be either polynomials or numbers, ``var``
    has to be a :py:class:`~.Symbol` and ``sampling_time`` can be both
    :py:class:`~.Symbol` or numbers.

    See :class:`TransferFunctionBase` for more information.

    Parameters
    ==========

    num : Expr, Number
        The numerator polynomial of the transfer function.
    den : Expr, Number
        The denominator polynomial of the transfer function.
    var : Symbol
        Complex variable of the z-transform used by the
        polynomials of the transfer function.
    sampling_time : Symbol, Number
        Time interval between two consecutive sampling instants
        Defaults to 1.
        If sampling_time == 0, an instance of :class:`TransferFunction` is
        returned.

    Raises
    ======

    TypeError
        When ``var`` is not a Symbol or when ``num`` or ``den`` is not a
        number or a polynomial.
    ValueError
        When ``den`` is zero.

    Examples
    ========

    >>> from sympy.abc import z, p, a, k
    >>> from sympy.physics.control.lti import DiscreteTransferFunction
    >>> dtf1 = DiscreteTransferFunction(z**2 + a*z, z**2 + z + 1, z, 0.1)
    >>> dtf1
    DiscreteTransferFunction(a*z + z**2, z**2 + z + 1, z, 0.1)
    >>> dtf1.num
    a*z + z**2
    >>> dtf1.den
    z**2 + z + 1
    >>> dtf1.var
    z
    >>> dtf1.sampling_time
    0.1
    >>> dtf1.args
    (a*z + z**2, z**2 + z + 1, z, 0.1)

    Any complex variable can be used for ``var``.

    >>> dtf2 = DiscreteTransferFunction(a*p**3 - a*p**2 + z*p, p + a**2, p, k)
    >>> dtf2
    DiscreteTransferFunction(a*p**3 - a*p**2 + p*z, a**2 + p, p, k)
    >>> dtf3 = DiscreteTransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p, k)
    >>> dtf3
    DiscreteTransferFunction((p - 1)*(p + 3), (p - 1)*(p + 5), p, k)

    To negate a transfer function the ``-`` operator can be prepended:

    >>> dtf4 = DiscreteTransferFunction(-a + z, p**2 + z, p, 0.2)
    >>> -dtf4
    DiscreteTransferFunction(a - z, p**2 + z, p, 0.2)
    >>> dtf5 = DiscreteTransferFunction(z**4 - 2*z**3 + 5*z + 4, z + 4, z, 12)
    >>> -dtf5
    DiscreteTransferFunction(-z**4 + 2*z**3 - 5*z - 4, z + 4, z, 12)

    You can use a float or an integer (or other constants) as numerator and
    denominator:

    >>> dtf6 = DiscreteTransferFunction(1/2, 4, z, k)
    >>> dtf6.num
    0.500000000000000
    >>> dtf6.den
    4
    >>> dtf6.var
    z
    >>> dtf6.sampling_time
    k
    >>> dtf6.args
    (0.5, 4, z, k)

    You can take the integer power of a transfer function using the ``**``
    operator:

    >>> dtf7 = DiscreteTransferFunction(z + a, z - a, z, 0.3)
    >>> dtf7**3
    DiscreteTransferFunction((a + z)**3, (-a + z)**3, z, 0.3)
    >>> dtf7**0
    DiscreteTransferFunction(1, 1, z, 0.3)
    >>> dtf8 = DiscreteTransferFunction(p + 4, p - 3, p, k)
    >>> dtf8**-1
    DiscreteTransferFunction(p - 3, p + 4, p, k)

    Addition, subtraction, and multiplication of transfer functions can form
    unevaluated ``Series`` or ``Parallel`` objects.

    >>> dtf9 = DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k)
    >>> dtf10 = DiscreteTransferFunction(z - p, z + 3, z, k)
    >>> dtf11 = DiscreteTransferFunction(4*z**2 + 2*z - 4, z - 1, z, k)
    >>> dtf12 = DiscreteTransferFunction(1 - z, z**2 + 4, z, k)
    >>> dtf9 + dtf10
    Parallel(DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k), DiscreteTransferFunction(-p + z, z + 3, z, k))
    >>> dtf10 - dtf11
    Parallel(DiscreteTransferFunction(-p + z, z + 3, z, k), DiscreteTransferFunction(-4*z**2 - 2*z + 4, z - 1, z, k))
    >>> dtf9 * dtf10
    Series(DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k), DiscreteTransferFunction(-p + z, z + 3, z, k))
    >>> dtf10 - (dtf9 + dtf12)
    Parallel(DiscreteTransferFunction(-p + z, z + 3, z, k), DiscreteTransferFunction(-z - 1, z**2 + z + 1, z, k), DiscreteTransferFunction(z - 1, z**2 + 4, z, k))
    >>> dtf10 - (dtf9 * dtf12)
    Parallel(DiscreteTransferFunction(-p + z, z + 3, z, k), Series(DiscreteTransferFunction(-1, 1, z, k), DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k), DiscreteTransferFunction(1 - z, z**2 + 4, z, k)))
    >>> dtf11 * dtf10 * dtf9
    Series(DiscreteTransferFunction(4*z**2 + 2*z - 4, z - 1, z, k), DiscreteTransferFunction(-p + z, z + 3, z, k), DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k))
    >>> dtf9 * dtf11 + dtf10 * dtf12
    Parallel(Series(DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k), DiscreteTransferFunction(4*z**2 + 2*z - 4, z - 1, z, k)), Series(DiscreteTransferFunction(-p + z, z + 3, z, k), DiscreteTransferFunction(1 - z, z**2 + 4, z, k)))
    >>> (dtf9 + dtf12) * (dtf10 + dtf11)
    Series(Parallel(DiscreteTransferFunction(z + 1, z**2 + z + 1, z, k), DiscreteTransferFunction(1 - z, z**2 + 4, z, k)), Parallel(DiscreteTransferFunction(-p + z, z + 3, z, k), DiscreteTransferFunction(4*z**2 + 2*z - 4, z - 1, z, k)))

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function using ``.doit()`` method or by ``.rewrite(DiscreteTransferFunction)``.

    >>> ((dtf9 + dtf10) * dtf12).doit()
    DiscreteTransferFunction((1 - z)*((-p + z)*(z**2 + z + 1) + (z + 1)*(z + 3)), (z + 3)*(z**2 + 4)*(z**2 + z + 1), z, k)
    >>> (dtf9 * dtf10 - dtf11 * dtf12).rewrite(DiscreteTransferFunction)
    DiscreteTransferFunction(-(1 - z)*(z + 3)*(z**2 + z + 1)*(4*z**2 + 2*z - 4) + (-p + z)*(z - 1)*(z + 1)*(z**2 + 4), (z - 1)*(z + 3)*(z**2 + 4)*(z**2 + z + 1), z, k)

    See Also
    ========

    TransferFunctionBase, TransferFunction, Feedback, Series, Parallel

    """
    def __new__(cls, num, den, var, sampling_time=1):
        if sampling_time == 0:
            raise ValueError(filldedent("""
                The sampling time cannot be zero.
                If you want to create a continuous transfer function,
                use the TransferFunction class instead."""))

        sampling_time = sympify(sampling_time)
        obj = super(DiscreteTransferFunction, cls).__new__(cls, num, den, var,
                                                      sampling_time)
        obj._sampling_time = sampling_time
        return obj

    @classmethod
    def from_rational_expression(cls, expr, var=None, sampling_time=1):
        r"""
        See :func:`TransferFunctionBase.from_rational_expression`.

        """
        return super().from_rational_expression(expr, var,
                                                sampling_time=sampling_time)

    @classmethod
    def from_coeff_lists(cls, num_list, den_list, var, sampling_time=1):
        r"""
        See :func:`TransferFunctionBase.from_coeff_lists`.

        """
        return super().from_coeff_lists(num_list, den_list, var,
                                        sampling_time=sampling_time)

    @classmethod
    def from_zpk(cls, zeros, poles, gain, var, sampling_time=1):
        r"""
        See :func:`TransferFunctionBase.from_zpk`.
        """
        return super().from_zpk(zeros, poles, gain,
                                var, sampling_time = sampling_time)

    @classmethod
    def from_gbt(cls, cont_tf, sampling_time, alpha, var):
        r"""
        Returns the discretized transfer function H(z) from a continuous
        transfer function H(s).

        Parameters
        ==========

        cont_tf : TransferFunction
            The continuous transfer function H(s) to be discretized.
        sampling_time : Symbol, Number
            Time interval between two consecutive sampling instants.
        alpha: Symbol, Number
            The parameter for the generalised bilinear transformation method.
        var: Symbol
            The complex variable of the z-transform used by the polynomials
            of the transfer function.

        Explanation
        ===========

        Where H(z) is the corresponding discretized transfer function,
        discretized with the generalised bilinear transformation method.
        H(z) is obtained from the continuous transfer function H(s)
        by substituting $s(z) = \frac{z-1}{T(\alpha z + (1-\alpha))}$ into H(s),
        where T is the sample time.

        Examples
        ========

        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> from sympy.abc import s, z, L, R, T

        >>> tf = TransferFunction(1, s*L + R, s)
        >>> dttf1 = DiscreteTransferFunction.from_gbt(tf, T, 0.5, z)
        >>> dttf1.num
        T*z/(2*(L + R*T/2)) + T/(2*(L + R*T/2))
        >>> dttf1.den
        z + (-L + R*T/2)/(L + R*T/2)
        >>> dttf1.sampling_time
        T

        >>> dttf2 = DiscreteTransferFunction.from_gbt(tf, T, 0, z)
        >>> dttf2.num
        T/L
        >>> dttf2.den
        z + (-L + R*T)/L
        >>> dttf2.sampling_time
        T

        >>> dttf3 = DiscreteTransferFunction.from_gbt(tf, T, 1, z)
        >>> dttf3.num
        T*z/(L + R*T)
        >>> dttf3.den
        -L/(L + R*T) + z
        >>> dttf3.sampling_time
        T

        >>> dttf4 = DiscreteTransferFunction.from_gbt(tf, T, 0.3, z)
        >>> dttf4.num
        3*T*z/(10*(L + 3*R*T/10)) + 7*T/(10*(L + 3*R*T/10))
        >>> dttf4.den
        z + (-L + 7*R*T/10)/(L + 3*R*T/10)
        >>> dttf4.sampling_time
        T

        References
        ==========

        .. [1] https://www.polyu.edu.hk/ama/profile/gfzhang/Research/ZCC09_IJC.pdf
        """
        num, den = gbt(cont_tf, sampling_time, alpha)
        return cls.from_coeff_lists(num, den, var, sampling_time)

    @classmethod
    def from_bilinear(cls, cont_tf, sampling_time, var):
        r"""
        Returns the discretized transfer function H(z) from a continuous
        transfer function H(s).

        Parameters
        ==========

        cont_tf : TransferFunction
            The continuous transfer function H(s) to be discretized.
        sampling_time : Symbol, Number
            Time interval between two consecutive sampling instants.
        var: Symbol
            The complex variable of the z-transform used by the polynomials
            of the transfer function.

        Explanation
        ===========

        Where H(z) is the corresponding discretized transfer function,
        discretized with the bilinear transform method.
        H(z) is obtained from the continuous transfer function H(s)
        by substituting $s(z) = \frac{2}{T}\frac{z-1}{z+1}$ into H(s), where T
        is the sample time.

        Examples
        ========

        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> from sympy.abc import s, z, L, R, T

        >>> tf = TransferFunction(1, s*L + R, s)
        >>> dttf = DiscreteTransferFunction.from_bilinear(tf, T, z)
        >>> dttf.num
        T*z/(2*(L + R*T/2)) + T/(2*(L + R*T/2))
        >>> dttf.den
        z + (-L + R*T/2)/(L + R*T/2)
        >>> dttf.sampling_time
        T

        """
        num, den = bilinear(cont_tf, sampling_time)
        return cls.from_coeff_lists(num, den, var, sampling_time)

    @classmethod
    def from_forward_diff(cls, cont_tf, sampling_time, var):
        r"""
        Returns the discretized transfer function H(z) from a continuous
        transfer function H(s).

        Parameters
        ==========

        cont_tf : TransferFunction
            The continuous transfer function H(s) to be discretized.
        sampling_time : Symbol, Number
            Time interval between two consecutive sampling instants.
        var: Symbol
            The complex variable of the z-transform used by the polynomials
            of the transfer function.

        Explanation
        ===========

        Where H(z) is the corresponding discretized transfer function,
        discretized with the forward difference transform method.
        H(z) is obtained from the continuous transfer function H(s)
        by substituting $s(z) = \frac{z-1}{T}$ into H(s), where T is the
        sample time.

        Examples
        ========

        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> from sympy.abc import s, z, L, R, T

        >>> tf = TransferFunction(1, s*L + R, s)
        >>> dttf = DiscreteTransferFunction.from_forward_diff(tf, T, z)
        >>> dttf.num
        T/L
        >>> dttf.den
        z + (-L + R*T)/L
        >>> dttf.sampling_time
        T

        """
        num, den = forward_diff(cont_tf, sampling_time)
        return cls.from_coeff_lists(num, den, var, sampling_time)

    @classmethod
    def from_backward_diff(cls, cont_tf, sampling_time, var):
        r"""
        Returns the discretized transfer function H(z) from a continuous
        transfer function H(s).

        Parameters
        ==========

        cont_tf : TransferFunction
            The continuous transfer function H(s) to be discretized.
        sampling_time : Symbol, Number
            Time interval between two consecutive sampling instants.
        var: Symbol
            The complex variable of the z-transform used by the polynomials
            of the transfer function.

        Explanation
        ===========

        Where H(z) is the corresponding discretized transfer function,
        discretized with the backward difference transform method.
        H(z) is obtained from the continuous transfer function H(s)
        by substituting $s(z) =  \frac{z-1}{Tz}$ into H(s), where T is the
        sample time.

        Examples
        ========

        >>> from sympy.physics.control.lti import TransferFunction, DiscreteTransferFunction
        >>> from sympy.abc import s, z, L, R, T

        >>> tf = TransferFunction(1, s*L + R, s)
        >>> dttf = DiscreteTransferFunction.from_backward_diff(tf, T, z)
        >>> dttf.num
        T*z/(L + R*T)
        >>> dttf.den
        -L/(L + R*T) + z
        >>> dttf.sampling_time
        T

        """
        num, den = backward_diff(cont_tf, sampling_time)
        return cls.from_coeff_lists(num, den, var, sampling_time)

    def dc_gain(self):
        r"""
        See :func:`TransferFunctionBase.dc_gain`.

        """
        m = Mul(self.num, Pow(self.den, -1, evaluate=False), evaluate=False)
        return limit(m, self.var, 1)

    def get_asymptotic_stability_conditions(
        self, cancel_poles_zeros=False, fast=False
    ) -> list[Boolean]:
        r"""
        See :func:`TransferFunctionBase.get_asymptotic_stability_conditions`.

        """
        standard_form = self.to_standard_form(cancel_poles_zeros)

        domain = EXRAW if fast else None

        p = Poly(standard_form.den, self.var, domain = domain)

        return [c > 0 for c in p.schur_conditions()]

    def _eval_rewrite_as_DiscreteStateSpace(self, *args):
        """
        Returns the equivalent space model of the transfer function model.
        The state space model will be returned in the controllable canonical
        form.

        Unlike the space state to transfer function model conversion, the
        transfer function to state space model conversion is not unique.
        There can be multiple state space representations of a given transfer function model.

        Examples
        ========

        >>> from sympy.abc import z
        >>> from sympy.physics.control import DiscreteTransferFunction, DiscreteStateSpace
        >>> dtf = DiscreteTransferFunction(z**2 + 1, z**3 + z*2 + 10, z, 0.1)
        >>> dtf.rewrite(DiscreteStateSpace)
        DiscreteStateSpace(Matrix([
        [  0,  1, 0],
        [  0,  0, 1],
        [-10, -2, 0]]), Matrix([
        [0],
        [0],
        [1]]), Matrix([[1, 0, 1]]), Matrix([[0]]), 0.1)

        """
        A, B, C, D = self._StateSpace_matrices_equivalent()
        return DiscreteStateSpace(A, B, C, D, self.sampling_time)

    def _eval_rewrite_as_StateSpace(self, *args):
        raise TypeError("""
            The discrete transfer function model cannot be rewritten as a
            continuous-time state space model.
            """)

    @property
    def sampling_time(self):
        """Returns the sampling time of the DiscreteTransferFunction."""
        return self._sampling_time

    _is_continuous = False

class PIDController(TransferFunction):
    r"""
    A class for representing PID (Proportional-Integral-Derivative)
    controllers in control systems. The PIDController class is a subclass
    of TransferFunction, representing the controller's transfer function
    in the Laplace domain. The arguments are ``kp``, ``ki``, ``kd``,
    ``tf``, and ``var``, where ``kp``, ``ki``, and ``kd`` are the
    proportional, integral, and derivative gains respectively.``tf``
    is the derivative filter time constant, which can be used to
    filter out the noise and ``var`` is the complex variable used in
    the transfer function.

    Parameters
    ==========

    kp : Expr, Number
        Proportional gain. Defaults to ``Symbol('kp')`` if not specified.
    ki : Expr, Number
        Integral gain. Defaults to ``Symbol('ki')`` if not specified.
    kd : Expr, Number
        Derivative gain. Defaults to ``Symbol('kd')`` if not specified.
    tf : Expr, Number
        Derivative filter time constant.  Defaults to ``0`` if not specified.
    var : Symbol
        The complex frequency variable.  Defaults to ``s`` if not specified.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.control.lti import PIDController
    >>> kp, ki, kd = symbols('kp ki kd')
    >>> p1 = PIDController(kp, ki, kd)
    >>> print(p1)
    PIDController(kp, ki, kd, 0, s)
    >>> p1.doit()
    TransferFunction(kd*s**2 + ki + kp*s, s, s)
    >>> p1.kp
    kp
    >>> p1.ki
    ki
    >>> p1.kd
    kd
    >>> p1.tf
    0
    >>> p1.var
    s
    >>> p1.to_expr()
    (kd*s**2 + ki + kp*s)/s

    See Also
    ========

    TransferFunction

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/PID_controller
    .. [2] https://in.mathworks.com/help/control/ug/proportional-integral-derivative-pid-controllers.html

    """
    def __new__(cls, kp=Symbol('kp'), ki=Symbol('ki'), kd=Symbol('kd'), tf=0, var=Symbol('s')):
        kp, ki, kd, tf = _sympify(kp), _sympify(ki), _sympify(kd), _sympify(tf)
        num = kp*tf*var**2 + kp*var + ki*tf*var + ki + kd*var**2
        den = tf*var**2 + var
        obj = TransferFunction.__new__(cls, num, den, var)
        return obj

    def __init__(self, kp=Symbol('kp'), ki=Symbol('ki'), kd=Symbol('kd'), tf=0, var=Symbol('s')):
        self._kp = kp
        self._ki = ki
        self._kd = kd
        self._tf = tf

    def __repr__(self):
        return f"PIDController({self.kp}, {self.ki}, {self.kd}, {self.tf}, {self.var})"

    __str__ = __repr__

    @property
    def kp(self):
        """
        Returns the Proportional gain (kp) of the PIDController.
        """
        return self._kp

    @property
    def ki(self):
        """
        Returns the Integral gain (ki) of the PIDController.
        """
        return self._ki

    @property
    def kd(self):
        """
        Returns the Derivative gain (kd) of the PIDController.
        """
        return self._kd

    @property
    def tf(self):
        """
        Returns the Derivative filter time constant (tf) of the PIDController.
        """
        return self._tf

    def doit(self):
        """
        Convert the PIDController into TransferFunction.
        """
        return TransferFunction(self.num, self.den, self.var)

    _is_continuous = True

def _flatten_args(args, _cls):
    temp_args = []
    for arg in args:
        if isinstance(arg, _cls):
            temp_args.extend(arg.args)
        else:
            temp_args.append(arg)
    return tuple(temp_args)


def _dummify_args(_arg, var):
    dummy_dict = {}
    dummy_arg_list = []

    for arg in _arg:
        _s = Dummy()
        dummy_dict[_s] = var
        dummy_arg = arg.subs({var: _s})
        dummy_arg_list.append(dummy_arg)

    return dummy_arg_list, dummy_dict

def _any_state_space_systems(arg_list: tuple) -> bool:
    """
    Check if there are any state space objects in the argument list.
    If there are, return True, else return False.

    Systems are considered to be state space systems if they are instances of
    StateSpaceBase or have an attribute `is_StateSpace_object` that is True.

    """
    return any(isinstance(arg, StateSpaceBase) or
               (hasattr(arg, 'is_StateSpace_object') and
                arg.is_StateSpace_object) for arg in arg_list)

def _are_input_output_compatible(args):
    """
    Check if the input and output of the systems are compatible for series
    connection. The input of the second system should match the output of the
    first system.

    """
    for i in range(1, len(args)):
        if args[i].num_inputs != args[i-1].num_outputs:
            return False
    return True

class Series(SISOLinearTimeInvariant):
    r"""
    A class for representing a series configuration of SISO systems.

    Parameters
    ==========

    args : SISOLinearTimeInvariant
        SISO systems in a series configuration.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``Series(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, SISO in this case.

    Examples
    ========

    >>> from sympy.abc import s, p, a, b
    >>> from sympy import Matrix
    >>> from sympy.physics.control.lti import TransferFunction, Series, Parallel, StateSpace
    >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
    >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
    >>> tf3 = TransferFunction(p**2, p + s, s)
    >>> S1 = Series(tf1, tf2)
    >>> S1
    Series(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s))
    >>> S1.var
    s
    >>> S2 = Series(tf2, Parallel(tf3, -tf1))
    >>> S2
    Series(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), Parallel(TransferFunction(p**2, p + s, s), TransferFunction(-a*p**2 - b*s, -p + s, s)))
    >>> S2.var
    s
    >>> S3 = Series(Parallel(tf1, tf2), Parallel(tf2, tf3))
    >>> S3
    Series(Parallel(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)), Parallel(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), TransferFunction(p**2, p + s, s)))
    >>> S3.var
    s

    You can get the resultant transfer function by using ``.doit()`` method:

    >>> S3 = Series(tf1, tf2, -tf3)
    >>> S3.doit()
    TransferFunction(-p**2*(s**3 - 2)*(a*p**2 + b*s), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)
    >>> S4 = Series(tf2, Parallel(tf1, -tf3))
    >>> S4.doit()
    TransferFunction((s**3 - 2)*(-p**2*(-p + s) + (p + s)*(a*p**2 + b*s)), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)

    You can also connect StateSpace which results in SISO

    >>> A1 = Matrix([[-1]])
    >>> B1 = Matrix([[1]])
    >>> C1 = Matrix([[-1]])
    >>> D1 = Matrix([1])
    >>> A2 = Matrix([[0]])
    >>> B2 = Matrix([[1]])
    >>> C2 = Matrix([[1]])
    >>> D2 = Matrix([[0]])
    >>> ss1 = StateSpace(A1, B1, C1, D1)
    >>> ss2 = StateSpace(A2, B2, C2, D2)
    >>> S5 = Series(ss1, ss2)
    >>> S5
    Series(StateSpace(Matrix([[-1]]), Matrix([[1]]), Matrix([[-1]]), Matrix([[1]])), StateSpace(Matrix([[0]]), Matrix([[1]]), Matrix([[1]]), Matrix([[0]])))
    >>> S5.doit()
    StateSpace(Matrix([
    [-1,  0],
    [-1, 0]]), Matrix([
    [1],
    [1]]), Matrix([[0, 1]]), Matrix([[0]]))

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    MIMOSeries, Parallel, TransferFunction, Feedback

    """
    def __new__(cls, *args, evaluate=False):
        args = _flatten_args(args, Series)
        _check_time_compatibility(args)
        obj = super().__new__(cls, *args)

        # For StateSpace series connection
        if args and _any_state_space_systems(args):
            # Check for SISO
            if any(not arg.is_SISO for arg in args):
                raise ValueError(
                    filldedent("""To use Series connection for MIMO systems use
                               MIMOSeries instead."""))
            # Check the interconnection
            if not _are_input_output_compatible(args):
                raise ValueError(
                    filldedent("""
                        Systems with incompatible inputs and outputs
                        cannot be connected in Series.
                        """))
            obj._is_series_StateSpace = True
        else:
            obj._is_series_StateSpace = False
            obj._check_args(args)

        obj._is_continuous = args[0].is_continuous

        return obj.doit() if evaluate else obj

    def __repr__(self):
        systems_repr = ', '.join(repr(system) for system in self.args)
        return f"Series({systems_repr})"

    __str__ = __repr__

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Series, Parallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> Series(G1, G2).var
        p
        >>> Series(-G3, Parallel(G1, G2)).var
        p

        """
        return self.args[0].var

    def doit(self, **hints):
        """
        Returns the resultant transfer function or state space obtained after
        evaluating the series interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> Series(tf2, tf1).doit()
        TransferFunction((s**3 - 2)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Series(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-a*p**2 - b*s), (-p + s)*(s**4 + 5*s + 6), s)

        Notes
        =====

        If a series connection contains only TransferFunctionBase components,
        the equivalent system returned will be a transfer function. However, if
        a StateSpaceBase object is used in any of the arguments,
        the output will be a state space object.

        """
        # Check if the system is a StateSpace
        if self._is_series_StateSpace:
            # Return the equivalent StateSpace model
            ss_class = StateSpace if self.is_continuous else DiscreteStateSpace

            res = self.args[0]
            if not isinstance(res, ss_class):
                res = res.doit().rewrite(ss_class)
            for arg in self.args[1:]:
                if not isinstance(arg, ss_class):
                    arg = arg.doit().rewrite(ss_class)
                res = arg * res
            return res

        _num_arg = (arg.doit().num for arg in self.args)
        _den_arg = (arg.doit().den for arg in self.args)
        res_num = Mul(*_num_arg, evaluate=True)
        res_den = Mul(*_den_arg, evaluate=True)

        sampling_time = self.args[0].sampling_time
        return create_transfer_function(res_num, res_den, self.var, sampling_time)

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        if not self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a discrete-time Series object as a
                    TransferFunction.""")
        if self._is_series_StateSpace:
            return self.doit().rewrite(TransferFunction)[0][0]
        return self.doit()

    def _eval_rewrite_as_DiscreteTransferFunction(self, *args, **kwargs):
        if self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a continuous-time Series object as a
                    DiscreteTransferFunction.""")
        if self._is_series_StateSpace:
            return self.doit().rewrite(DiscreteTransferFunction)[0][0]
        return self.doit()

    @_compatibility_decorator
    @_check_other_SISO
    def __add__(self, other):

        if isinstance(other, Parallel):
            arg_list = list(other.args)
            return Parallel(self, *arg_list)

        return Parallel(self, other)

    __radd__ = __add__

    @_compatibility_decorator
    @_check_other_SISO
    def __sub__(self, other):
        return self + (-other)

    @_compatibility_decorator
    def __rsub__(self, other):
        return -self + other

    @_compatibility_decorator
    @_check_other_SISO
    def __mul__(self, other):

        arg_list = list(self.args)
        return Series(*arg_list, other)

    @_compatibility_decorator
    def __truediv__(self, other):
        tf_class = TransferFunction if self.is_continuous\
                            else DiscreteTransferFunction

        if isinstance(other, TransferFunctionBase):
            return Series(*self.args, create_transfer_function(other.den, other.num,
                                             other.var, other.sampling_time))
        elif isinstance(other, Series):
            tf_self = self.rewrite(tf_class)
            tf_other = other.rewrite(tf_class)
            return tf_self / tf_other
        elif (isinstance(other, Parallel) and len(other.args) == 2
            and isinstance(other.args[0], tf_class) and \
                isinstance(other.args[1], Series)):

            if not self.var == other.var:
                raise ValueError(filldedent("""
                    All the transfer functions should use the same complex
                    variable of the Laplace transform or z-transform."""))
            self_arg_list = set(self.args)
            other_arg_list = set(other.args[1].args)
            res = list(self_arg_list ^ other_arg_list)
            if len(res) == 0:
                return Feedback(self, other.args[0])
            elif len(res) == 1:
                return Feedback(self, *res)
            else:
                return Feedback(self, Series(*res))
        else:
            raise ValueError("This transfer function expression is invalid.")

    def __neg__(self):
        return Series(create_transfer_function(-1, 1, self.var, self.sampling_time), self)

    def to_expr(self):
        """Returns the equivalent ``Expr`` object."""
        return Mul(*(arg.to_expr() for arg in self.args), evaluate=False)

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is less than or equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*s + 2, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> S1 = Series(-tf2, tf1)
        >>> S1.is_proper
        False
        >>> S2 = Series(tf1, tf2, tf3)
        >>> S2.is_proper
        True

        """
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is strictly less than degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**2 + 5*s + 6, s)
        >>> tf3 = TransferFunction(1, s**2 + s + 1, s)
        >>> S1 = Series(tf1, tf2)
        >>> S1.is_strictly_proper
        False
        >>> S2 = Series(tf1, tf2, tf3)
        >>> S2.is_strictly_proper
        True

        """
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        r"""
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Series
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(p, s**2, s)
        >>> tf3 = TransferFunction(s**2, 1, s)
        >>> S1 = Series(tf1, -tf2)
        >>> S1.is_biproper
        False
        >>> S2 = Series(tf2, tf3)
        >>> S2.is_biproper
        True

        """
        return self.doit().is_biproper

    @property
    def is_StateSpace_object(self):
        return self._is_series_StateSpace

    @property
    def sampling_time(self):
        return self.args[0].sampling_time

class MIMOSeries(MIMOLinearTimeInvariant):
    r"""
    A class for representing a series configuration of MIMO systems.

    Parameters
    ==========

    args : MIMOLinearTimeInvariant
        MIMO systems in a series configuration.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``MIMOSeries(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.

        ``num_outputs`` of the MIMO system is not equal to the
        ``num_inputs`` of its adjacent MIMO system. (Matrix
        multiplication constraint, basically)
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, MIMO in this case.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import MIMOSeries, TransferFunctionMatrix, StateSpace
    >>> from sympy import Matrix, pprint
    >>> mat_a = Matrix([[5*s], [5]])  # 2 Outputs 1 Input
    >>> mat_b = Matrix([[5, 1/(6*s**2)]])  # 1 Output 2 Inputs
    >>> mat_c = Matrix([[1, s], [5/s, 1]])  # 2 Outputs 2 Inputs
    >>> tfm_a = TransferFunctionMatrix.from_Matrix(mat_a, s)
    >>> tfm_b = TransferFunctionMatrix.from_Matrix(mat_b, s)
    >>> tfm_c = TransferFunctionMatrix.from_Matrix(mat_c, s)
    >>> MIMOSeries(tfm_c, tfm_b, tfm_a)
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, 1, s), TransferFunction(s, 1, s)), (TransferFunction(5, s, s), TransferFunction(1, 1, s)))), TransferFunctionMatrix(((TransferFunction(5, 1, s), TransferFunction(1, 6*s**2, s)),)), TransferFunctionMatrix(((TransferFunction(5*s, 1, s),), (TransferFunction(5, 1, s),))))
    >>> pprint(_, use_unicode=False)  #  For Better Visualization
    [5*s]                 [1  s]
    [---]    [5   1  ]    [-  -]
    [ 1 ]    [-  ----]    [1  1]
    [   ]   *[1     2]   *[    ]
    [ 5 ]    [   6*s ]{t} [5  1]
    [ - ]                 [-  -]
    [ 1 ]{t}              [s  1]{t}
    >>> MIMOSeries(tfm_c, tfm_b, tfm_a).doit()
    TransferFunctionMatrix(((TransferFunction(150*s**4 + 25*s, 6*s**3, s), TransferFunction(150*s**4 + 5*s, 6*s**2, s)), (TransferFunction(150*s**3 + 25, 6*s**3, s), TransferFunction(150*s**3 + 5, 6*s**2, s))))
    >>> pprint(_, use_unicode=False)  # (2 Inputs -A-> 2 Outputs) -> (2 Inputs -B-> 1 Output) -> (1 Input -C-> 2 Outputs) is equivalent to (2 Inputs -Series Equivalent-> 2 Outputs).
    [     4              4      ]
    [150*s  + 25*s  150*s  + 5*s]
    [-------------  ------------]
    [        3             2    ]
    [     6*s           6*s     ]
    [                           ]
    [      3              3     ]
    [ 150*s  + 25    150*s  + 5 ]
    [ -----------    ---------- ]
    [        3             2    ]
    [     6*s           6*s     ]{t}
    >>> a1 = Matrix([[4, 1], [2, -3]])
    >>> b1 = Matrix([[5, 2], [-3, -3]])
    >>> c1 = Matrix([[2, -4], [0, 1]])
    >>> d1 = Matrix([[3, 2], [1, -1]])
    >>> a2 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    >>> b2 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    >>> c2 = Matrix([[4, 2, -3], [1, 4, 3]])
    >>> d2 = Matrix([[-2, 4], [0, 1]])
    >>> ss1 = StateSpace(a1, b1, c1, d1) #2 inputs, 2 outputs
    >>> ss2 = StateSpace(a2, b2, c2, d2) #2 inputs, 2 outputs
    >>> S1 = MIMOSeries(ss1, ss2) #(2 inputs, 2 outputs) -> (2 inputs, 2 outputs)
    >>> S1
    MIMOSeries(StateSpace(Matrix([
            [4,  1],
            [2, -3]]), Matrix([
            [ 5,  2],
            [-3, -3]]), Matrix([
            [2, -4],
            [0,  1]]), Matrix([
            [3,  2],
            [1, -1]])), StateSpace(Matrix([
            [-3,  4, 2],
            [-1, -3, 0],
            [ 2,  5, 3]]), Matrix([
            [ 1,  4],
            [-3, -3],
            [-2,  1]]), Matrix([
            [4, 2, -3],
            [1, 4,  3]]), Matrix([
            [-2, 4],
            [ 0, 1]])))
    >>> S1.doit()
    StateSpace(Matrix([
    [ 4,  1,  0,  0, 0],
    [ 2, -3,  0,  0, 0],
    [ 2,  0, -3,  4, 2],
    [-6,  9, -1, -3, 0],
    [-4,  9,  2,  5, 3]]), Matrix([
    [  5,  2],
    [ -3, -3],
    [  7, -2],
    [-12, -3],
    [ -5, -5]]), Matrix([
    [-4, 12, 4, 2, -3],
    [ 0,  1, 1, 4,  3]]), Matrix([
    [-2, -8],
    [ 1, -1]]))

    Notes
    =====

    All the transfer function matrices should use the same complex variable ``var`` of the Laplace transform.

    ``MIMOSeries(A, B)`` is not equivalent to ``A*B``. It is always in the reverse order, that is ``B*A``.

    See Also
    ========

    Series, MIMOParallel

    """
    def __new__(cls, *args, evaluate=False):
        args = _flatten_args(args, MIMOSeries)
        obj = super().__new__(cls, *args)
        if args and _any_state_space_systems(args):
            # Check compatibility
            if not _are_input_output_compatible(args):
                raise ValueError(filldedent("""
                    Systems with incompatible inputs and outputs
                    cannot be connected in MIMOSeries."""))
            obj._is_series_StateSpace = True

        else:
            cls._check_args(args)
            obj._is_series_StateSpace = False

            if not _are_input_output_compatible(args):
                raise ValueError(filldedent("""
                    Number of input signals do not match the number
                    of output signals of adjacent systems for some args."""))

        _check_time_compatibility(args)
        obj._is_continuous = args[0].is_continuous

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, MIMOSeries, TransferFunctionMatrix
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> tfm_1 = TransferFunctionMatrix([[G1, G2, G3]])
        >>> tfm_2 = TransferFunctionMatrix([[G1], [G2], [G3]])
        >>> MIMOSeries(tfm_2, tfm_1).var
        p

        """
        return self.args[0].var

    @property
    def num_inputs(self):
        """Returns the number of input signals of the series system."""
        return self.args[0].num_inputs

    @property
    def num_outputs(self):
        """Returns the number of output signals of the series system."""
        return self.args[-1].num_outputs

    @property
    def shape(self):
        """Returns the shape of the equivalent MIMO system."""
        return self.num_outputs, self.num_inputs

    @property
    def is_StateSpace_object(self):
        return self._is_series_StateSpace

    def doit(self, cancel=False, **kwargs):
        """
        Returns the resultant obtained after evaluating the MIMO systems arranged
        in a series configuration.
        For transfer function systems it returns a TransferFunctionMatrix
        and for state space systems it returns the resultant state space system.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, MIMOSeries, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tfm1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf2]])
        >>> tfm2 = TransferFunctionMatrix([[tf2, tf1], [tf1, tf1]])
        >>> MIMOSeries(tfm2, tfm1).doit()
        TransferFunctionMatrix(((TransferFunction(2*(-p + s)*(s**3 - 2)*(a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)**2*(s**4 + 5*s + 6)**2, s), TransferFunction((-p + s)**2*(s**3 - 2)*(a*p**2 + b*s) + (-p + s)*(a*p**2 + b*s)**2*(s**4 + 5*s + 6), (-p + s)**3*(s**4 + 5*s + 6), s)), (TransferFunction((-p + s)*(s**3 - 2)**2*(s**4 + 5*s + 6) + (s**3 - 2)*(a*p**2 + b*s)*(s**4 + 5*s + 6)**2, (-p + s)*(s**4 + 5*s + 6)**3, s), TransferFunction(2*(s**3 - 2)*(a*p**2 + b*s), (-p + s)*(s**4 + 5*s + 6), s))))

        """
        if self._is_series_StateSpace:
            # Return the equivalent StateSpace model
            ss_class = StateSpace if self.is_continuous else DiscreteStateSpace
            res = self.args[0]
            if not isinstance(res, ss_class):
                res = res.doit().rewrite(ss_class)
            for arg in self.args[1:]:
                if not isinstance(arg, ss_class):
                    arg = arg.doit().rewrite(ss_class)
                res = arg * res
            return res

        _arg = (arg.doit()._expr_mat for arg in reversed(self.args))

        if cancel:
            res = MatMul(*_arg, evaluate=True)
            return TransferFunctionMatrix.from_Matrix(res, self.var)

        _dummy_args, _dummy_dict = _dummify_args(_arg, self.var)
        res = MatMul(*_dummy_args, evaluate=True)
        temp_tfm = TransferFunctionMatrix.from_Matrix(res, self.var,
                                                      self.sampling_time)
        return temp_tfm.subs(_dummy_dict)

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        tf_class = TransferFunction if self.is_continuous \
            else DiscreteTransferFunction
        if self._is_series_StateSpace:
            return self.doit().rewrite(tf_class)
        return self.doit()

    @_check_other_MIMO
    def __add__(self, other):

        if isinstance(other, MIMOParallel):
            arg_list = list(other.args)
            return MIMOParallel(self, *arg_list)

        return MIMOParallel(self, other)

    __radd__ = __add__

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_MIMO
    def __mul__(self, other):

        if isinstance(other, MIMOSeries):
            self_arg_list = list(self.args)
            other_arg_list = list(other.args)
            return MIMOSeries(*other_arg_list, *self_arg_list)  # A*B = MIMOSeries(B, A)

        arg_list = list(self.args)
        return MIMOSeries(other, *arg_list)

    def __neg__(self):
        arg_list = list(self.args)
        arg_list[0] = -arg_list[0]
        return MIMOSeries(*arg_list)

    @property
    def sampling_time(self):
        return self.args[0].sampling_time

class Parallel(SISOLinearTimeInvariant):
    r"""
    A class for representing a parallel configuration of SISO systems.

    Parameters
    ==========

    args : SISOLinearTimeInvariant
        SISO systems in a parallel arrangement.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``Parallel(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.abc import s, p, a, b
    >>> from sympy.physics.control.lti import TransferFunction, Parallel, Series, StateSpace
    >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
    >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
    >>> tf3 = TransferFunction(p**2, p + s, s)
    >>> P1 = Parallel(tf1, tf2)
    >>> P1
    Parallel(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s))
    >>> P1.var
    s
    >>> P2 = Parallel(tf2, Series(tf3, -tf1))
    >>> P2
    Parallel(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), Series(TransferFunction(p**2, p + s, s), TransferFunction(-a*p**2 - b*s, -p + s, s)))
    >>> P2.var
    s
    >>> P3 = Parallel(Series(tf1, tf2), Series(tf2, tf3))
    >>> P3
    Parallel(Series(TransferFunction(a*p**2 + b*s, -p + s, s), TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)), Series(TransferFunction(s**3 - 2, s**4 + 5*s + 6, s), TransferFunction(p**2, p + s, s)))
    >>> P3.var
    s

    You can get the resultant transfer function by using ``.doit()`` method:

    >>> Parallel(tf1, tf2, -tf3).doit()
    TransferFunction(-p**2*(-p + s)*(s**4 + 5*s + 6) + (-p + s)*(p + s)*(s**3 - 2) + (p + s)*(a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)
    >>> Parallel(tf2, Series(tf1, -tf3)).doit()
    TransferFunction(-p**2*(a*p**2 + b*s)*(s**4 + 5*s + 6) + (-p + s)*(p + s)*(s**3 - 2), (-p + s)*(p + s)*(s**4 + 5*s + 6), s)

    Parallel can be used to connect SISO ``StateSpace`` systems together.

    >>> A1 = Matrix([[-1]])
    >>> B1 = Matrix([[1]])
    >>> C1 = Matrix([[-1]])
    >>> D1 = Matrix([1])
    >>> A2 = Matrix([[0]])
    >>> B2 = Matrix([[1]])
    >>> C2 = Matrix([[1]])
    >>> D2 = Matrix([[0]])
    >>> ss1 = StateSpace(A1, B1, C1, D1)
    >>> ss2 = StateSpace(A2, B2, C2, D2)
    >>> P4 = Parallel(ss1, ss2)
    >>> P4
    Parallel(StateSpace(Matrix([[-1]]), Matrix([[1]]), Matrix([[-1]]), Matrix([[1]])), StateSpace(Matrix([[0]]), Matrix([[1]]), Matrix([[1]]), Matrix([[0]])))

    ``doit()`` can be used to find ``StateSpace`` equivalent for the system containing ``StateSpace`` objects.

    >>> P4.doit()
    StateSpace(Matrix([
    [-1, 0],
    [ 0, 0]]), Matrix([
    [1],
    [1]]), Matrix([[-1, 1]]), Matrix([[1]]))
    >>> P4.rewrite(TransferFunction)
    TransferFunction(s*(s + 1) + 1, s*(s + 1), s)

    Notes
    =====

    All the transfer functions should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Series, TransferFunction, Feedback

    """

    def __new__(cls, *args, evaluate=False):
        args = _flatten_args(args, Parallel)
        obj = super().__new__(cls, *args)
        # For StateSpace parallel connection
        if args and _any_state_space_systems(args):
            # Check for SISO
            if any(not arg.is_SISO for arg in args):
                raise ValueError(filldedent("""
                    To use Parallel connection for MIMO systems use
                    MIMOParallel instead."""))
            obj._is_parallel_StateSpace = True
        else:
            obj._is_parallel_StateSpace = False
            obj._check_args(args)

        _check_time_compatibility(args)
        obj._is_continuous = args[0].is_continuous

        return obj.doit() if evaluate else obj

    def __repr__(self):
        systems_repr = ', '.join(repr(system) for system in self.args)
        return f"Parallel({systems_repr})"

    __str__ = __repr__

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Parallel, Series
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> Parallel(G1, G2).var
        p
        >>> Parallel(-G3, Series(G1, G2)).var
        p

        """
        return self.args[0].var

    def doit(self, **hints):
        """
        Returns the resultant transfer function or state space obtained by
        parallel connection of transfer functions or state space objects.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> Parallel(tf2, tf1).doit()
        TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)
        >>> Parallel(-tf1, -tf2).doit()
        TransferFunction((2 - s**3)*(-p + s) + (-a*p**2 - b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)

        """
        if self._is_parallel_StateSpace:
            # Return the equivalent StateSpace model
            ss_class = StateSpace if self.is_continuous else DiscreteStateSpace

            res = self.args[0].doit()
            if not isinstance(res, ss_class):
                res = res.rewrite(ss_class)
            for arg in self.args[1:]:
                if not isinstance(arg, ss_class):
                    arg = arg.doit().rewrite(ss_class)
                res += arg
            return res

        _arg = (arg.doit().to_expr() for arg in self.args)
        res = Add(*_arg).as_numer_denom()

        sampling_time = self.args[0].sampling_time
        return create_transfer_function(*res, self.var, sampling_time)

    def _eval_rewrite_as_TransferFunction(self, *args, **kwargs):
        if not self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a discrete-time Parallel object as a
                    TransferFunction.""")
        if self._is_parallel_StateSpace:
            return self.doit().rewrite(TransferFunction)[0][0]
        return self.doit()

    def _eval_rewrite_as_DiscreteTransferFunction(self, *args, **kwargs):
        if self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a continuous-time Parallel object as a
                    DiscreteTransferFunction.""")
        if self._is_parallel_StateSpace:
            return self.doit().rewrite(DiscreteTransferFunction)[0][0]
        return self.doit()

    @_compatibility_decorator
    @_check_other_SISO
    def __add__(self, other):

        self_arg_list = list(self.args)
        return Parallel(*self_arg_list, other)

    __radd__ = __add__

    @_compatibility_decorator
    @_check_other_SISO
    def __sub__(self, other):
        return self + (-other)

    @_compatibility_decorator
    def __rsub__(self, other):
        return -self + other

    @_compatibility_decorator
    @_check_other_SISO
    def __mul__(self, other):

        if isinstance(other, Series):
            arg_list = list(other.args)
            return Series(self, *arg_list)

        return Series(self, other)

    def __neg__(self):
        return Series(create_transfer_function(-1, 1, self.var, self.sampling_time), self)

    def to_expr(self):
        """Returns the equivalent ``Expr`` object."""
        return Add(*(arg.to_expr() for arg in self.args), evaluate=False)

    @property
    def is_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is less than or equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
        >>> tf2 = TransferFunction(p**2 - 4*p, p**3 + 3*s + 2, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(-tf2, tf1)
        >>> P1.is_proper
        False
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_proper
        True

        """
        return self.doit().is_proper

    @property
    def is_strictly_proper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is strictly less than degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(tf1, tf2)
        >>> P1.is_strictly_proper
        False
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_strictly_proper
        True

        """
        return self.doit().is_strictly_proper

    @property
    def is_biproper(self):
        """
        Returns True if degree of the numerator polynomial of the resultant transfer
        function is equal to degree of the denominator polynomial of
        the same, else False.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Parallel
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(p**2, p + s, s)
        >>> tf3 = TransferFunction(s, s**2 + s + 1, s)
        >>> P1 = Parallel(tf1, -tf2)
        >>> P1.is_biproper
        True
        >>> P2 = Parallel(tf2, tf3)
        >>> P2.is_biproper
        False

        """
        return self.doit().is_biproper

    @property
    def is_StateSpace_object(self):
        return self._is_parallel_StateSpace

    @property
    def sampling_time(self):
        return self.args[0].sampling_time

class MIMOParallel(MIMOLinearTimeInvariant):
    r"""
    A class for representing a parallel configuration of MIMO systems.

    Parameters
    ==========

    args : MIMOLinearTimeInvariant
        MIMO Systems in a parallel arrangement.
    evaluate : Boolean, Keyword
        When passed ``True``, returns the equivalent
        ``MIMOParallel(*args).doit()``. Set to ``False`` by default.

    Raises
    ======

    ValueError
        When no argument is passed.

        ``var`` attribute is not same for every system.

        All MIMO systems passed do not have same shape.
    TypeError
        Any of the passed ``*args`` has unsupported type

        A combination of SISO and MIMO systems is
        passed. There should be homogeneity in the
        type of systems passed, MIMO in this case.

    Examples
    ========

    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import TransferFunctionMatrix, MIMOParallel, StateSpace
    >>> from sympy import Matrix, pprint
    >>> expr_1 = 1/s
    >>> expr_2 = s/(s**2-1)
    >>> expr_3 = (2 + s)/(s**2 - 1)
    >>> expr_4 = 5
    >>> tfm_a = TransferFunctionMatrix.from_Matrix(Matrix([[expr_1, expr_2], [expr_3, expr_4]]), s)
    >>> tfm_b = TransferFunctionMatrix.from_Matrix(Matrix([[expr_2, expr_1], [expr_4, expr_3]]), s)
    >>> tfm_c = TransferFunctionMatrix.from_Matrix(Matrix([[expr_3, expr_4], [expr_1, expr_2]]), s)
    >>> MIMOParallel(tfm_a, tfm_b, tfm_c)
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)), (TransferFunction(s + 2, s**2 - 1, s), TransferFunction(5, 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s**2 - 1, s), TransferFunction(1, s, s)), (TransferFunction(5, 1, s), TransferFunction(s + 2, s**2 - 1, s)))), TransferFunctionMatrix(((TransferFunction(s + 2, s**2 - 1, s), TransferFunction(5, 1, s)), (TransferFunction(1, s, s), TransferFunction(s, s**2 - 1, s)))))
    >>> pprint(_, use_unicode=False)  #  For Better Visualization
    [  1       s   ]      [  s       1   ]      [s + 2     5   ]
    [  -     ------]      [------    -   ]      [------    -   ]
    [  s      2    ]      [ 2        s   ]      [ 2        1   ]
    [        s  - 1]      [s  - 1        ]      [s  - 1        ]
    [              ]    + [              ]    + [              ]
    [s + 2     5   ]      [  5     s + 2 ]      [  1       s   ]
    [------    -   ]      [  -     ------]      [  -     ------]
    [ 2        1   ]      [  1      2    ]      [  s      2    ]
    [s  - 1        ]{t}   [        s  - 1]{t}   [        s  - 1]{t}
    >>> MIMOParallel(tfm_a, tfm_b, tfm_c).doit()
    TransferFunctionMatrix(((TransferFunction(s**2 + s*(2*s + 2) - 1, s*(s**2 - 1), s), TransferFunction(2*s**2 + 5*s*(s**2 - 1) - 1, s*(s**2 - 1), s)), (TransferFunction(s**2 + s*(s + 2) + 5*s*(s**2 - 1) - 1, s*(s**2 - 1), s), TransferFunction(5*s**2 + 2*s - 3, s**2 - 1, s))))
    >>> pprint(_, use_unicode=False)
    [       2                              2       / 2    \    ]
    [      s  + s*(2*s + 2) - 1         2*s  + 5*s*\s  - 1/ - 1]
    [      --------------------         -----------------------]
    [             / 2    \                       / 2    \      ]
    [           s*\s  - 1/                     s*\s  - 1/      ]
    [                                                          ]
    [ 2                   / 2    \             2               ]
    [s  + s*(s + 2) + 5*s*\s  - 1/ - 1      5*s  + 2*s - 3     ]
    [---------------------------------      --------------     ]
    [              / 2    \                      2             ]
    [            s*\s  - 1/                     s  - 1         ]{t}

    ``MIMOParallel`` can also be used to connect MIMO ``StateSpace`` systems.

    >>> A1 = Matrix([[4, 1], [2, -3]])
    >>> B1 = Matrix([[5, 2], [-3, -3]])
    >>> C1 = Matrix([[2, -4], [0, 1]])
    >>> D1 = Matrix([[3, 2], [1, -1]])
    >>> A2 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    >>> B2 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    >>> C2 = Matrix([[4, 2, -3], [1, 4, 3]])
    >>> D2 = Matrix([[-2, 4], [0, 1]])
    >>> ss1 = StateSpace(A1, B1, C1, D1)
    >>> ss2 = StateSpace(A2, B2, C2, D2)
    >>> p1 = MIMOParallel(ss1, ss2)
    >>> p1
    MIMOParallel(StateSpace(Matrix([
    [4,  1],
    [2, -3]]), Matrix([
    [ 5,  2],
    [-3, -3]]), Matrix([
    [2, -4],
    [0,  1]]), Matrix([
    [3,  2],
    [1, -1]])), StateSpace(Matrix([
    [-3,  4, 2],
    [-1, -3, 0],
    [ 2,  5, 3]]), Matrix([
    [ 1,  4],
    [-3, -3],
    [-2,  1]]), Matrix([
    [4, 2, -3],
    [1, 4,  3]]), Matrix([
    [-2, 4],
    [ 0, 1]])))

    ``doit()`` can be used to find ``StateSpace`` equivalent for the system containing ``StateSpace`` objects.

    >>> p1.doit()
    StateSpace(Matrix([
    [4,  1,  0,  0, 0],
    [2, -3,  0,  0, 0],
    [0,  0, -3,  4, 2],
    [0,  0, -1, -3, 0],
    [0,  0,  2,  5, 3]]), Matrix([
    [ 5,  2],
    [-3, -3],
    [ 1,  4],
    [-3, -3],
    [-2,  1]]), Matrix([
    [2, -4, 4, 2, -3],
    [0,  1, 1, 4,  3]]), Matrix([
    [1, 6],
    [1, 0]]))

    Notes
    =====

    All the transfer function matrices should use the same complex variable
    ``var`` of the Laplace transform.

    See Also
    ========

    Parallel, MIMOSeries

    """

    def __new__(cls, *args, evaluate=False):
        args = _flatten_args(args, MIMOParallel)
        obj = super().__new__(cls, *args)
        # For StateSpace Parallel connection
        if args and _any_state_space_systems(args):
            if not _are_input_output_compatible(args):
                raise ShapeError("Systems with incompatible inputs and outputs"
                                 "cannot be connected in MIMOParallel.")
            obj._is_parallel_StateSpace = True
        else:
            obj._check_args(args)
            if any(arg.shape != args[0].shape for arg in args):
                raise TypeError("Shape of all the args is not equal.")
            obj._is_parallel_StateSpace = False

        _check_time_compatibility(args)
        obj._is_continuous = args[0].is_continuous

        return obj.doit() if evaluate else obj

    @property
    def var(self):
        """
        Returns the complex variable used by all the systems.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOParallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> G4 = TransferFunction(p**2, p**2 - 1, p)
        >>> tfm_a = TransferFunctionMatrix([[G1, G2], [G3, G4]])
        >>> tfm_b = TransferFunctionMatrix([[G2, G1], [G4, G3]])
        >>> MIMOParallel(tfm_a, tfm_b).var
        p

        """
        return self.args[0].var

    @property
    def num_inputs(self):
        """Returns the number of input signals of the parallel system."""
        return self.args[0].num_inputs

    @property
    def num_outputs(self):
        """Returns the number of output signals of the parallel system."""
        return self.args[0].num_outputs

    @property
    def shape(self):
        """Returns the shape of the equivalent MIMO system."""
        return self.num_outputs, self.num_inputs

    @property
    def is_StateSpace_object(self):
        return self._is_parallel_StateSpace

    def doit(self, **hints):
        """
        Returns the resultant transfer function matrix or StateSpace obtained
        after evaluating the MIMO systems arranged in a parallel configuration.

        Examples
        ========

        >>> from sympy.abc import s, p, a, b
        >>> from sympy.physics.control.lti import TransferFunction, MIMOParallel, TransferFunctionMatrix
        >>> tf1 = TransferFunction(a*p**2 + b*s, s - p, s)
        >>> tf2 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> tfm_2 = TransferFunctionMatrix([[tf2, tf1], [tf1, tf2]])
        >>> MIMOParallel(tfm_1, tfm_2).doit()
        TransferFunctionMatrix(((TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s)), (TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s), TransferFunction((-p + s)*(s**3 - 2) + (a*p**2 + b*s)*(s**4 + 5*s + 6), (-p + s)*(s**4 + 5*s + 6), s))))

        """
        if self._is_parallel_StateSpace:
            # Return the equivalent StateSpace model.
            res = self.args[0]
            if not isinstance(res, StateSpace):
                res = res.doit().rewrite(StateSpace)
            for arg in self.args[1:]:
                if not isinstance(arg, StateSpace):
                    arg = arg.doit().rewrite(StateSpace)
                else:
                    arg = arg.doit()
                res += arg
            return res
        _arg = (arg.doit()._expr_mat for arg in self.args)
        res = MatAdd(*_arg, evaluate=True)
        return TransferFunctionMatrix.from_Matrix(res, self.var,
                                                  self.sampling_time)

    def _eval_rewrite_as_TransferFunctionMatrix(self, *args, **kwargs):
        tf_class = TransferFunction if self.is_continuous \
                                    else DiscreteTransferFunction
        if self._is_parallel_StateSpace:
            return self.doit().rewrite(tf_class)
        return self.doit()

    @_check_other_MIMO
    def __add__(self, other):

        self_arg_list = list(self.args)
        return MIMOParallel(*self_arg_list, other)

    __radd__ = __add__

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    @_check_other_MIMO
    def __mul__(self, other):

        if isinstance(other, MIMOSeries):
            arg_list = list(other.args)
            return MIMOSeries(*arg_list, self)

        return MIMOSeries(other, self)

    def __neg__(self):
        arg_list = [-arg for arg in list(self.args)]
        return MIMOParallel(*arg_list)

    @property
    def sampling_time(self):
        return self.args[0].sampling_time

class Feedback(SISOLinearTimeInvariant):
    r"""
    A class for representing closed-loop feedback interconnection between two
    SISO input/output systems.

    The first argument, ``sys1``, is the feedforward part of the closed-loop
    system or in simple words, the dynamical model representing the process
    to be controlled. The second argument, ``sys2``, is the feedback system
    and controls the fed back signal to ``sys1``. Both ``sys1`` and ``sys2``
    can either be ``Series``, state space or transfer function objects.

    Parameters
    ==========

    sys1 : Series, StateSpaceBase, TransferFunctionBase
        The feedforward path system.
    sys2 : Series, StateSpaceBase, TransferFunctionBase, optional
        The feedback path system (often a feedback controller).
        It is the model sitting on the feedback path.

        If not specified explicitly, the sys2 is
        assumed to be unit (1.0) transfer function.
    sign : int, optional
        The sign of feedback. Can either be ``1``
        (for positive feedback) or ``-1`` (for negative feedback).
        Default value is `-1`.

    Raises
    ======

    ValueError
        When ``sys1`` and ``sys2`` are not using the
        same complex variable of the Laplace transform or z-transform.

        When a combination of ``sys1`` and ``sys2`` yields
        zero denominator.

    TypeError
        When either ``sys1`` or ``sys2`` is not a ``Series``, ``StateSpaceBase``
        or ``TransferFunctionBase`` object.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import StateSpace, TransferFunction, Feedback
    >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> controller = TransferFunction(5*s - 10, s + 7, s)
    >>> F1 = Feedback(plant, controller)
    >>> F1
    Feedback(TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s), -1)
    >>> F1.var
    s
    >>> F1.args
    (TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s), TransferFunction(5*s - 10, s + 7, s), -1)

    You can get the feedforward and feedback path systems by using ``.sys1`` and ``.sys2`` respectively.

    >>> F1.sys1
    TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
    >>> F1.sys2
    TransferFunction(5*s - 10, s + 7, s)

    You can get the resultant closed loop transfer function obtained by negative feedback
    interconnection using ``.doit()`` method.

    >>> F1.doit()
    TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
    >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
    >>> C = TransferFunction(5*s + 10, s + 10, s)
    >>> F2 = Feedback(G*C, TransferFunction(1, 1, s))
    >>> F2.doit()
    TransferFunction((s + 10)*(5*s + 10)*(s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s + 10)*((s + 10)*(s**2 + 2*s + 3) + (5*s + 10)*(2*s**2 + 5*s + 1))*(s**2 + 2*s + 3), s)

    To negate a ``Feedback`` object, the ``-`` operator can be prepended:

    >>> -F1
    Feedback(TransferFunction(-3*s**2 - 7*s + 3, s**2 - 4*s + 2, s), TransferFunction(10 - 5*s, s + 7, s), -1)
    >>> -F2
    Feedback(Series(TransferFunction(-1, 1, s), TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s), TransferFunction(5*s + 10, s + 10, s)), TransferFunction(-1, 1, s), -1)

    ``Feedback`` can also be used to connect SISO ``StateSpace`` systems together.

    >>> A1 = Matrix([[-1]])
    >>> B1 = Matrix([[1]])
    >>> C1 = Matrix([[-1]])
    >>> D1 = Matrix([1])
    >>> A2 = Matrix([[0]])
    >>> B2 = Matrix([[1]])
    >>> C2 = Matrix([[1]])
    >>> D2 = Matrix([[0]])
    >>> ss1 = StateSpace(A1, B1, C1, D1)
    >>> ss2 = StateSpace(A2, B2, C2, D2)
    >>> F3 = Feedback(ss1, ss2)
    >>> F3
    Feedback(StateSpace(Matrix([[-1]]), Matrix([[1]]), Matrix([[-1]]), Matrix([[1]])), StateSpace(Matrix([[0]]), Matrix([[1]]), Matrix([[1]]), Matrix([[0]])), -1)

    ``doit()`` can be used to find ``StateSpace`` equivalent for the system containing ``StateSpace`` objects.

    >>> F3.doit()
    StateSpace(Matrix([
    [-1, -1],
    [-1, -1]]), Matrix([
    [1],
    [1]]), Matrix([[-1, -1]]), Matrix([[1]]))

    We can also find the equivalent ``TransferFunction`` by using ``rewrite(TransferFunction)`` method.

    >>> F3.rewrite(TransferFunction)
    TransferFunction(s, s + 2, s)

    See Also
    ========

    MIMOFeedback, Series, Parallel

    """
    def __new__(cls, sys1, sys2=None, sign=-1):
        if not sys2:
            sys2 = create_transfer_function(1, 1, sys1.var, sys1.sampling_time)

        if not isinstance(sys1, (TransferFunctionBase, Series, StateSpaceBase,
                                 Feedback)):
            raise TypeError("Unsupported type for `sys1` in Feedback.")

        if not isinstance(sys2, (TransferFunctionBase, Series, StateSpaceBase,
                                 Feedback)):
            raise TypeError("Unsupported type for `sys2` in Feedback.")

        if not (sys1.num_inputs == sys1.num_outputs == sys2.num_inputs ==
                sys2.num_outputs == 1):
            raise ValueError(filldedent("""To use Feedback connection for MIMO systems
                            use MIMOFeedback instead."""))

        if sign not in [-1, 1]:
            raise ValueError(filldedent("""
                Unsupported type for feedback. `sign` arg should
                either be 1 (positive feedback loop) or -1
                (negative feedback loop)."""))

        obj = super(SISOLinearTimeInvariant, cls).__new__(cls, sys1, sys2,
                                                          _sympify(sign))

        if sys1.is_StateSpace_object or sys2.is_StateSpace_object:
            obj.is_StateSpace_object = True
        else:
            if Mul(sys1.to_expr(), sys2.to_expr()).simplify() == sign:
                raise ValueError(filldedent("""The equivalent system will have zero
                                 denominator."""))
            if sys1.var != sys2.var:
                raise ValueError(filldedent("""Both `sys1` and `sys2` should be
                    using the same complex variable."""))
            obj.is_StateSpace_object = False

        _check_time_compatibility([sys1, sys2])
        obj._is_continuous = sys1.is_continuous

        return obj

    def __repr__(self):
        return f"Feedback({self.sys1}, {self.sys2}, {self.sign})"

    __str__ = __repr__

    @property
    def sys1(self):
        """
        Returns the feedforward system of the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.sys1
        TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.sys1
        TransferFunction(1, 1, p)

        """
        return self.args[0]

    @property
    def sys2(self):
        """
        Returns the feedback controller of the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.sys2
        TransferFunction(5*s - 10, s + 7, s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.sys2
        Series(TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p), TransferFunction(5*p + 10, p + 10, p), TransferFunction(1 - s, p + 2, p))

        """
        return self.args[1]

    @property
    def var(self):
        """
        Returns the complex variable of the Laplace transform used by all
        the transfer functions involved in the feedback interconnection.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.var
        s
        >>> G = TransferFunction(2*s**2 + 5*s + 1, p**2 + 2*p + 3, p)
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - s, p + 2, p)
        >>> F2 = Feedback(TransferFunction(1, 1, p), G*C*P)
        >>> F2.var
        p

        """
        return self.sys1.var

    @property
    def sign(self):
        """
        Returns the type of MIMO Feedback model. ``1``
        for Positive and ``-1`` for Negative.
        """
        return self.args[2]

    @property
    def num(self):
        """
        Returns the numerator of the closed loop feedback system.
        """
        return self.sys1

    @property
    def den(self):
        """
        Returns the denominator of the closed loop feedback model.
        """
        unit = create_transfer_function(1, 1, self.var, self.args[0].sampling_time)
        arg_list = list(self.sys1.args) if isinstance(self.sys1, Series) else [self.sys1]
        if self.sign == 1:
            return Parallel(unit, -Series(self.sys2, *arg_list))
        return Parallel(unit, Series(self.sys2, *arg_list))

    @property
    def sensitivity(self):
        """
        Returns the sensitivity function of the feedback loop.

        Sensitivity of a Feedback system is the ratio
        of change in the open loop gain to the change in
        the closed loop gain.

        .. note::
            This method would not return the complementary
            sensitivity function.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> C = TransferFunction(5*p + 10, p + 10, p)
        >>> P = TransferFunction(1 - p, p + 2, p)
        >>> F_1 = Feedback(P, C)
        >>> F_1.sensitivity
        1/((1 - p)*(5*p + 10)/((p + 2)*(p + 10)) + 1)

        """

        return 1/(1 - self.sign*self.sys1.to_expr()*self.sys2.to_expr())

    def doit(self, cancel=False, expand=False, **hints):
        """
        Returns the resultant transfer function or state space obtained by
        feedback connection of transfer functions or state space objects.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy import Matrix
        >>> from sympy.physics.control.lti import TransferFunction, Feedback, StateSpace
        >>> plant = TransferFunction(3*s**2 + 7*s - 3, s**2 - 4*s + 2, s)
        >>> controller = TransferFunction(5*s - 10, s + 7, s)
        >>> F1 = Feedback(plant, controller)
        >>> F1.doit()
        TransferFunction((s + 7)*(s**2 - 4*s + 2)*(3*s**2 + 7*s - 3), ((s + 7)*(s**2 - 4*s + 2) + (5*s - 10)*(3*s**2 + 7*s - 3))*(s**2 - 4*s + 2), s)
        >>> G = TransferFunction(2*s**2 + 5*s + 1, s**2 + 2*s + 3, s)
        >>> F2 = Feedback(G, TransferFunction(1, 1, s))
        >>> F2.doit()
        TransferFunction((s**2 + 2*s + 3)*(2*s**2 + 5*s + 1), (s**2 + 2*s + 3)*(3*s**2 + 7*s + 4), s)

        Use kwarg ``expand=True`` to expand the resultant transfer function.
        Use ``cancel=True`` to cancel out the common terms in numerator and
        denominator.

        >>> F2.doit(cancel=True, expand=True)
        TransferFunction(2*s**2 + 5*s + 1, 3*s**2 + 7*s + 4, s)
        >>> F2.doit(expand=True)
        TransferFunction(2*s**4 + 9*s**3 + 17*s**2 + 17*s + 3, 3*s**4 + 13*s**3 + 27*s**2 + 29*s + 12, s)

        If the connection contain any ``StateSpace`` object then ``doit()``
        will return the equivalent ``StateSpace`` object.

        >>> A1 = Matrix([[-1.5, -2], [1, 0]])
        >>> B1 = Matrix([0.5, 0])
        >>> C1 = Matrix([[0, 1]])
        >>> A2 = Matrix([[0, 1], [-5, -2]])
        >>> B2 = Matrix([0, 3])
        >>> C2 = Matrix([[0, 1]])
        >>> ss1 = StateSpace(A1, B1, C1)
        >>> ss2 = StateSpace(A2, B2, C2)
        >>> F3 = Feedback(ss1, ss2)
        >>> F3.doit()
        StateSpace(Matrix([
        [-1.5, -2,  0, -0.5],
        [   1,  0,  0,    0],
        [   0,  0,  0,    1],
        [   0,  3, -5,   -2]]), Matrix([
        [0.5],
        [  0],
        [  0],
        [  0]]), Matrix([[0, 1, 0, 0]]), Matrix([[0]]))

        """
        if self.is_StateSpace_object:
            ss_class = StateSpace if self.is_continuous else DiscreteStateSpace

            sys1_ss = self.sys1.doit().rewrite(ss_class)
            sys2_ss = self.sys2.doit().rewrite(ss_class)
            A1, B1, C1, D1 = sys1_ss.A, sys1_ss.B, sys1_ss.C, sys1_ss.D
            A2, B2, C2, D2 = sys2_ss.A, sys2_ss.B, sys2_ss.C, sys2_ss.D

            # Create identity matrices
            I_inputs = eye(self.num_inputs)
            I_outputs = eye(self.num_outputs)

            # Compute F and its inverse
            F = I_inputs - self.sign * D2 * D1
            E = F.inv()

            # Compute intermediate matrices
            E_D2 = E * D2
            E_C2 = E * C2
            T1 = I_outputs + self.sign * D1 * E_D2
            T2 = I_inputs + self.sign * E_D2 * D1
            A = Matrix.vstack(
            Matrix.hstack(A1 + self.sign * B1 * E_D2 * C1, self.sign * B1 * E_C2),
            Matrix.hstack(B2 * T1 * C1, A2 + self.sign * B2 * D1 * E_C2)
            )
            B = Matrix.vstack(B1 * T2, B2 * D1 * T2)
            C = Matrix.hstack(T1 * C1, self.sign * D1 * E_C2)
            D = D1 * T2
            return create_state_space(A, B, C, D, self.sampling_time)

        arg_list = list(self.sys1.args) if isinstance(self.sys1, Series) else [self.sys1]
        # F_n and F_d are resultant TFs of num and den of Feedback.
        F_n = self.sys1.doit()
        unit = create_transfer_function(1, 1, self.sys1.var, self.sys1.sampling_time)

        if self.sign == -1:
            F_d = Parallel(unit, Series(self.sys2, *arg_list)).doit()
        else:
            F_d = Parallel(unit, -Series(self.sys2, *arg_list)).doit()

        _resultant_tf = create_transfer_function(F_n.num * F_d.den, F_n.den * F_d.num,
                               F_n.var, self.sys1.sampling_time)

        if cancel:
            _resultant_tf = _resultant_tf.simplify()

        if expand:
            _resultant_tf = _resultant_tf.expand()

        return _resultant_tf

    def _eval_rewrite_as_TransferFunction(self, num, den, sign, **kwargs):
        if not self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a discrete-time Feedback object as a
                    TransferFunction.""")
        if self.is_StateSpace_object:
            return self.doit().rewrite(TransferFunction)[0][0]
        return self.doit()

    def _eval_rewrite_as_DiscreteTransferFunction(self, *args, **kwargs):
        if self.is_continuous:
            raise TypeError("""
                    Cannot rewrite a continuous-time Feedback object as a
                    DiscreteTransferFunction.""")
        if self.is_StateSpace_object:
            return self.doit().rewrite(DiscreteTransferFunction)[0][0]
        return self.doit()

    def to_expr(self):
        """
        Converts a ``Feedback`` object to SymPy Expr.

        Examples
        ========

        >>> from sympy.abc import s, a, b
        >>> from sympy.physics.control.lti import TransferFunction, Feedback
        >>> from sympy import Expr
        >>> tf1 = TransferFunction(a+s, 1, s)
        >>> tf2 = TransferFunction(b+s, 1, s)
        >>> fd1 = Feedback(tf1, tf2)
        >>> fd1.to_expr()
        (a + s)/((a + s)*(b + s) + 1)
        >>> isinstance(_, Expr)
        True
        """

        return self.doit().to_expr()

    def __neg__(self):
        return Feedback(-self.sys1, -self.sys2, self.sign)

    @property
    def sampling_time(self):
        return self.sys1.sampling_time


def _is_invertible(a, b, sign):
    """
    Checks whether a given pair of MIMO
    systems passed is invertible or not.
    """
    _mat = eye(a.num_outputs) - sign*(a.doit()._expr_mat)*(b.doit()._expr_mat)
    _det = _mat.det()

    return _det != 0


class MIMOFeedback(MIMOLinearTimeInvariant):
    r"""
    A class for representing closed-loop feedback interconnection between two
    MIMO input/output systems.

    Parameters
    ==========

    sys1 : MIMOSeries, TransferFunctionMatrix, StateSpaceBase
        The MIMO system placed on the feedforward path.
    sys2 : MIMOSeries, TransferFunctionMatrix, StateSpaceBase
        The system placed on the feedback path
        (often a feedback controller).
    sign : int, optional
        The sign of feedback. Can either be ``1``
        (for positive feedback) or ``-1`` (for negative feedback).
        Default value is `-1`.

    Raises
    ======

    ValueError
        When ``sys1`` and ``sys2`` are not using the
        same complex variable of the Laplace transform or z-transform.

        Forward path model should have an equal number of inputs/outputs
        to the feedback path outputs/inputs.

        When product of ``sys1`` and ``sys2`` is not a square matrix.

        When the equivalent MIMO system is not invertible.

    TypeError
        When either ``sys1`` or ``sys2`` is not a ``MIMOSeries``,
        ``TransferFunctionMatrix`` or a ``StateSpaceBase`` object.

    Examples
    ========

    >>> from sympy import Matrix, pprint
    >>> from sympy.abc import s
    >>> from sympy.physics.control.lti import StateSpace, TransferFunctionMatrix, MIMOFeedback
    >>> plant_mat = Matrix([[1, 1/s], [0, 1]])
    >>> controller_mat = Matrix([[10, 0], [0, 10]])  # Constant Gain
    >>> plant = TransferFunctionMatrix.from_Matrix(plant_mat, s)
    >>> controller = TransferFunctionMatrix.from_Matrix(controller_mat, s)
    >>> feedback = MIMOFeedback(plant, controller)  # Negative Feedback (default)
    >>> pprint(feedback, use_unicode=False)
    /    [1  1]    [10  0 ]   \-1   [1  1]
    |    [-  -]    [--  - ]   |     [-  -]
    |    [1  s]    [1   1 ]   |     [1  s]
    |I + [    ]   *[      ]   |   * [    ]
    |    [0  1]    [0   10]   |     [0  1]
    |    [-  -]    [-   --]   |     [-  -]
    \    [1  1]{t} [1   1 ]{t}/     [1  1]{t}

    To get the equivalent system matrix, use either ``doit`` or ``rewrite`` method.

    >>> pprint(feedback.doit(), use_unicode=False)
    [1     1  ]
    [--  -----]
    [11  121*s]
    [         ]
    [0    1   ]
    [-    --  ]
    [1    11  ]{t}

    To negate the ``MIMOFeedback`` object, use ``-`` operator.

    >>> neg_feedback = -feedback
    >>> pprint(neg_feedback.doit(), use_unicode=False)
    [-1    -1  ]
    [---  -----]
    [11   121*s]
    [          ]
    [ 0    -1  ]
    [ -    --- ]
    [ 1    11  ]{t}

    ``MIMOFeedback`` can also be used to connect MIMO ``StateSpace`` systems.

    >>> A1 = Matrix([[4, 1], [2, -3]])
    >>> B1 = Matrix([[5, 2], [-3, -3]])
    >>> C1 = Matrix([[2, -4], [0, 1]])
    >>> D1 = Matrix([[3, 2], [1, -1]])
    >>> A2 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    >>> B2 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    >>> C2 = Matrix([[4, 2, -3], [1, 4, 3]])
    >>> D2 = Matrix([[-2, 4], [0, 1]])
    >>> ss1 = StateSpace(A1, B1, C1, D1)
    >>> ss2 = StateSpace(A2, B2, C2, D2)
    >>> F1 = MIMOFeedback(ss1, ss2)
    >>> F1
    MIMOFeedback(StateSpace(Matrix([
    [4,  1],
    [2, -3]]), Matrix([
    [ 5,  2],
    [-3, -3]]), Matrix([
    [2, -4],
    [0,  1]]), Matrix([
    [3,  2],
    [1, -1]])), StateSpace(Matrix([
    [-3,  4, 2],
    [-1, -3, 0],
    [ 2,  5, 3]]), Matrix([
    [ 1,  4],
    [-3, -3],
    [-2,  1]]), Matrix([
    [4, 2, -3],
    [1, 4,  3]]), Matrix([
    [-2, 4],
    [ 0, 1]])), -1)

    ``doit()`` can be used to find ``StateSpace`` equivalent for the system containing ``StateSpace`` objects.

    >>> F1.doit()
    StateSpace(Matrix([
    [   3,  -3/4, -15/4, -37/2, -15],
    [ 7/2, -39/8,   9/8,  39/4,   9],
    [   3, -41/4, -45/4, -51/2, -19],
    [-9/2, 129/8,  73/8, 171/4,  36],
    [-3/2,  47/8,  31/8,  85/4,  18]]), Matrix([
    [-1/4,  19/4],
    [ 3/8, -21/8],
    [ 1/4,  29/4],
    [ 3/8, -93/8],
    [ 5/8, -35/8]]), Matrix([
    [  1, -15/4,  -7/4, -21/2, -9],
    [1/2, -13/8, -13/8, -19/4, -3]]), Matrix([
    [-1/4, 11/4],
    [ 1/8,  9/8]]))

    See Also
    ========

    Feedback, MIMOSeries, MIMOParallel

    """
    def __new__(cls, sys1, sys2, sign=-1):
        if not isinstance(sys1,
                          (TransferFunctionMatrix, MIMOSeries, StateSpaceBase)):
            raise TypeError("Unsupported type for `sys1` in MIMO Feedback.")

        if not isinstance(sys2,
                          (TransferFunctionMatrix, MIMOSeries, StateSpaceBase)):
            raise TypeError("Unsupported type for `sys2` in MIMO Feedback.")

        if sys1.num_inputs != sys2.num_outputs or \
            sys1.num_outputs != sys2.num_inputs:
            raise ValueError(filldedent("""
                Product of `sys1` and `sys2` must
                yield a square matrix."""))

        if sign not in (-1, 1):
            raise ValueError(filldedent("""
                Unsupported type for feedback. `sign` arg should
                either be 1 (positive feedback loop) or -1
                (negative feedback loop)."""))

        obj = super().__new__(cls, sys1, sys2, _sympify(sign))

        if sys1.is_StateSpace_object or sys2.is_StateSpace_object:
            obj.is_StateSpace_object = True
        else:
            if not _is_invertible(sys1, sys2, sign):
                raise ValueError("Non-Invertible system inputted.")
            obj.is_StateSpace_object = False

        if not obj.is_StateSpace_object and sys1.var != sys2.var:
            raise ValueError(filldedent("""
                Both `sys1` and `sys2` should be using the
                same complex variable."""))

        _check_time_compatibility([sys1, sys2])
        obj._is_continuous = sys1.is_continuous

        return obj

    @property
    def sys1(self):
        r"""
        Returns the system placed on the feedforward path of the MIMO feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s**2 + s + 1, s**2 - s + 1, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(1, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf3, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)
        >>> F_1.sys1
        TransferFunctionMatrix(((TransferFunction(s**2 + s + 1, s**2 - s + 1, s), TransferFunction(1, s, s)), (TransferFunction(1, s, s), TransferFunction(s**2 + s + 1, s**2 - s + 1, s))))
        >>> pprint(_, use_unicode=False)
        [ 2                    ]
        [s  + s + 1      1     ]
        [----------      -     ]
        [ 2              s     ]
        [s  - s + 1            ]
        [                      ]
        [             2        ]
        [    1       s  + s + 1]
        [    -       ----------]
        [    s        2        ]
        [            s  - s + 1]{t}

        """
        return self.args[0]

    @property
    def sys2(self):
        r"""
        Returns the feedback controller of the MIMO feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s**2, s**3 - s + 1, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(1, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2)
        >>> F_1.sys2
        TransferFunctionMatrix(((TransferFunction(s**2, s**3 - s + 1, s), TransferFunction(1, 1, s)), (TransferFunction(1, 1, s), TransferFunction(1, s, s))))
        >>> pprint(_, use_unicode=False)
        [     2       ]
        [    s       1]
        [----------  -]
        [ 3          1]
        [s  - s + 1   ]
        [             ]
        [    1       1]
        [    -       -]
        [    1       s]{t}

        """
        return self.args[1]

    @property
    def var(self):
        r"""
        Returns the complex variable of the Laplace transform used by all
        the transfer functions involved in the MIMO feedback loop.

        Examples
        ========

        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(p, 1 - p, p)
        >>> tf2 = TransferFunction(1, p, p)
        >>> tf3 = TransferFunction(1, 1, p)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)  # Positive feedback
        >>> F_1.var
        p

        """
        return self.sys1.var

    @property
    def sign(self):
        r"""
        Returns the type of feedback interconnection of two models. ``1``
        for Positive and ``-1`` for Negative.
        """
        return self.args[2]

    @property
    def sensitivity(self):
        r"""
        Returns the sensitivity function matrix of the feedback loop.

        Sensitivity of a closed-loop system is the ratio of change
        in the open loop gain to the change in the closed loop gain.

        .. note::
            This method would not return the complementary
            sensitivity function.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(p, 1 - p, p)
        >>> tf2 = TransferFunction(1, p, p)
        >>> tf3 = TransferFunction(1, 1, p)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf2, tf1]])
        >>> sys2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf2]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)  # Positive feedback
        >>> F_2 = MIMOFeedback(sys1, sys2)  # Negative feedback
        >>> pprint(F_1.sensitivity, use_unicode=False)
        [   4      3      2               5      4      2           ]
        [- p  + 3*p  - 4*p  + 3*p - 1    p  - 2*p  + 3*p  - 3*p + 1 ]
        [----------------------------  -----------------------------]
        [  4      3      2              5      4      3      2      ]
        [ p  + 3*p  - 8*p  + 8*p - 3   p  + 3*p  - 8*p  + 8*p  - 3*p]
        [                                                           ]
        [       4    3    2                  3      2               ]
        [      p  - p  - p  + p           3*p  - 6*p  + 4*p - 1     ]
        [ --------------------------    --------------------------  ]
        [  4      3      2               4      3      2            ]
        [ p  + 3*p  - 8*p  + 8*p - 3    p  + 3*p  - 8*p  + 8*p - 3  ]
        >>> pprint(F_2.sensitivity, use_unicode=False)
        [ 4      3      2           5      4      2          ]
        [p  - 3*p  + 2*p  + p - 1  p  - 2*p  + 3*p  - 3*p + 1]
        [------------------------  --------------------------]
        [   4      3                   5      4      2       ]
        [  p  - 3*p  + 2*p - 1        p  - 3*p  + 2*p  - p   ]
        [                                                    ]
        [     4    3    2               4      3             ]
        [    p  - p  - p  + p        2*p  - 3*p  + 2*p - 1   ]
        [  -------------------       ---------------------   ]
        [   4      3                   4      3              ]
        [  p  - 3*p  + 2*p - 1        p  - 3*p  + 2*p - 1    ]

        """
        _sys1_mat = self.sys1.doit()._expr_mat
        _sys2_mat = self.sys2.doit()._expr_mat

        return (eye(self.sys1.num_inputs) - \
            self.sign*_sys1_mat*_sys2_mat).inv()

    @property
    def num_inputs(self):
        """Returns the number of inputs of the system."""
        return self.sys1.num_inputs

    @property
    def num_outputs(self):
        """Returns the number of outputs of the system."""
        return self.sys1.num_outputs

    def doit(self, cancel=True, expand=False, **hints):
        r"""
        Returns the resultant transfer function matrix obtained by the
        feedback interconnection.

        Examples
        ========

        >>> from sympy import pprint
        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, MIMOFeedback
        >>> tf1 = TransferFunction(s, 1 - s, s)
        >>> tf2 = TransferFunction(1, s, s)
        >>> tf3 = TransferFunction(5, 1, s)
        >>> tf4 = TransferFunction(s - 1, s, s)
        >>> tf5 = TransferFunction(0, 1, s)
        >>> sys1 = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
        >>> sys2 = TransferFunctionMatrix([[tf3, tf5], [tf5, tf5]])
        >>> F_1 = MIMOFeedback(sys1, sys2, 1)
        >>> pprint(F_1, use_unicode=False)
        /    [  s      1  ]    [5  0]   \-1   [  s      1  ]
        |    [-----    -  ]    [-  -]   |     [-----    -  ]
        |    [1 - s    s  ]    [1  1]   |     [1 - s    s  ]
        |I - [            ]   *[    ]   |   * [            ]
        |    [  5    s - 1]    [0  0]   |     [  5    s - 1]
        |    [  -    -----]    [-  -]   |     [  -    -----]
        \    [  1      s  ]{t} [1  1]{t}/     [  1      s  ]{t}
        >>> pprint(F_1.doit(), use_unicode=False)
        [  -s           s - 1       ]
        [-------     -----------    ]
        [6*s - 1     s*(6*s - 1)    ]
        [                           ]
        [5*s - 5  (s - 1)*(6*s + 24)]
        [-------  ------------------]
        [6*s - 1     s*(6*s - 1)    ]{t}

        If the user wants the resultant ``TransferFunctionMatrix`` object without
        canceling the common factors then the ``cancel`` kwarg should be passed ``False``.

        >>> pprint(F_1.doit(cancel=False), use_unicode=False)
        [             s*(s - 1)                              s - 1               ]
        [         -----------------                       -----------            ]
        [         (1 - s)*(6*s - 1)                       s*(6*s - 1)            ]
        [                                                                        ]
        [s*(25*s - 25) + 5*(1 - s)*(6*s - 1)  s*(s - 1)*(6*s - 1) + s*(25*s - 25)]
        [-----------------------------------  -----------------------------------]
        [         (1 - s)*(6*s - 1)                        2                     ]
        [                                                 s *(6*s - 1)           ]{t}

        If the user wants the expanded form of the resultant transfer function matrix,
        the ``expand`` kwarg should be passed as ``True``.

        >>> pprint(F_1.doit(expand=True), use_unicode=False)
        [  -s          s - 1      ]
        [-------      --------    ]
        [6*s - 1         2        ]
        [             6*s  - s    ]
        [                         ]
        [            2            ]
        [5*s - 5  6*s  + 18*s - 24]
        [-------  ----------------]
        [6*s - 1         2        ]
        [             6*s  - s    ]{t}

        """
        if self.is_StateSpace_object:
            ss_class = StateSpace if self.is_continuous else DiscreteStateSpace
            sys1_ss = self.sys1.doit().rewrite(ss_class)
            sys2_ss = self.sys2.doit().rewrite(ss_class)
            A1, B1, C1, D1 = sys1_ss.A, sys1_ss.B, sys1_ss.C, sys1_ss.D
            A2, B2, C2, D2 = sys2_ss.A, sys2_ss.B, sys2_ss.C, sys2_ss.D

            # Create identity matrices
            I_inputs = eye(self.num_inputs)
            I_outputs = eye(self.num_outputs)

            # Compute F and its inverse
            F = I_inputs - self.sign * D2 * D1
            E = F.inv()

            # Compute intermediate matrices
            E_D2 = E * D2
            E_C2 = E * C2
            T1 = I_outputs + self.sign * D1 * E_D2
            T2 = I_inputs + self.sign * E_D2 * D1
            A = Matrix.vstack(
            Matrix.hstack(A1 + self.sign * B1 * E_D2 * C1, self.sign * B1 * E_C2),
            Matrix.hstack(B2 * T1 * C1, A2 + self.sign * B2 * D1 * E_C2)
            )
            B = Matrix.vstack(B1 * T2, B2 * D1 * T2)
            C = Matrix.hstack(T1 * C1, self.sign * D1 * E_C2)
            D = D1 * T2
            return create_state_space(A, B, C, D, self.sampling_time)

        _mat = self.sensitivity * self.sys1.doit()._expr_mat

        _resultant_tfm = _to_TFM(_mat, self.var, self.sampling_time)

        if cancel:
            _resultant_tfm = _resultant_tfm.simplify()

        if expand:
            _resultant_tfm = _resultant_tfm.expand()

        return _resultant_tfm

    def _eval_rewrite_as_TransferFunctionMatrix(self, sys1, sys2, sign, **kwargs):
        return self.doit()

    def __neg__(self):
        return MIMOFeedback(-self.sys1, -self.sys2, self.sign)

    @property
    def sampling_time(self):
        return self.sys1.sampling_time


def _to_TFM(mat, var, sampling_time):
    """Private method to convert ImmutableMatrix to TransferFunctionMatrix
    efficiently"""
    if sampling_time == 0:
        to_tf = lambda expr: \
            TransferFunction.from_rational_expression(expr,var)
    else:
        to_tf = lambda expr: \
            DiscreteTransferFunction.from_rational_expression(expr,var, sampling_time)
    arg = [[to_tf(expr) for expr in row] for row in mat.tolist()]
    return TransferFunctionMatrix(arg)


class TransferFunctionMatrix(MIMOLinearTimeInvariant):
    r"""
    A class for representing the MIMO (multiple-input and multiple-output)
    generalization of the SISO (single-input and single-output) transfer function.

    It is a matrix of transfer functions (``TransferFunction``, SISO-``Series`` or SISO-``Parallel``).
    There is only one argument, ``arg`` which is also the compulsory argument.
    ``arg`` is expected to be strictly of the type list of lists
    which holds the transfer functions or reducible to transfer functions.

    Parameters
    ==========

    arg : Nested ``List`` (strictly).
        Users are expected to input a nested list of ``TransferFunction``, ``Series``
        and/or ``Parallel`` objects.

    Examples
    ========

    .. note::
        ``pprint()`` can be used for better visualization of ``TransferFunctionMatrix`` objects.

    >>> from sympy.abc import s, p, a
    >>> from sympy import pprint
    >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, Series, Parallel
    >>> tf_1 = TransferFunction(s + a, s**2 + s + 1, s)
    >>> tf_2 = TransferFunction(p**4 - 3*p + 2, s + p, s)
    >>> tf_3 = TransferFunction(3, s + 2, s)
    >>> tf_4 = TransferFunction(-a + p, 9*s - 9, s)
    >>> tfm_1 = TransferFunctionMatrix([[tf_1], [tf_2], [tf_3]])
    >>> tfm_1
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(3, s + 2, s),)))
    >>> tfm_1.var
    s
    >>> tfm_1.num_inputs
    1
    >>> tfm_1.num_outputs
    3
    >>> tfm_1.shape
    (3, 1)
    >>> tfm_1.args
    (((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(3, s + 2, s),)),)
    >>> tfm_2 = TransferFunctionMatrix([[tf_1, -tf_3], [tf_2, -tf_1], [tf_3, -tf_2]])
    >>> tfm_2
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(-a - s, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-p**4 + 3*p - 2, p + s, s))))
    >>> pprint(tfm_2, use_unicode=False)  # pretty-printing for better visualization
    [   a + s           -3       ]
    [ ----------       -----     ]
    [  2               s + 2     ]
    [ s  + s + 1                 ]
    [                            ]
    [ 4                          ]
    [p  - 3*p + 2      -a - s    ]
    [------------    ----------  ]
    [   p + s         2          ]
    [                s  + s + 1  ]
    [                            ]
    [                 4          ]
    [     3        - p  + 3*p - 2]
    [   -----      --------------]
    [   s + 2          p + s     ]{t}

    TransferFunctionMatrix can be transposed, if user wants to switch the input and output transfer functions

    >>> tfm_2.transpose()
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(p**4 - 3*p + 2, p + s, s), TransferFunction(3, s + 2, s)), (TransferFunction(-3, s + 2, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s))))
    >>> pprint(_, use_unicode=False)
    [             4                          ]
    [  a + s     p  - 3*p + 2        3       ]
    [----------  ------------      -----     ]
    [ 2             p + s          s + 2     ]
    [s  + s + 1                              ]
    [                                        ]
    [                             4          ]
    [   -3          -a - s     - p  + 3*p - 2]
    [  -----      ----------   --------------]
    [  s + 2       2               p + s     ]
    [             s  + s + 1                 ]{t}

    >>> tf_5 = TransferFunction(5, s, s)
    >>> tf_6 = TransferFunction(5*s, (2 + s**2), s)
    >>> tf_7 = TransferFunction(5, (s*(2 + s**2)), s)
    >>> tf_8 = TransferFunction(5, 1, s)
    >>> tfm_3 = TransferFunctionMatrix([[tf_5, tf_6], [tf_7, tf_8]])
    >>> tfm_3
    TransferFunctionMatrix(((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)), (TransferFunction(5, s*(s**2 + 2), s), TransferFunction(5, 1, s))))
    >>> pprint(tfm_3, use_unicode=False)
    [    5        5*s  ]
    [    -       ------]
    [    s        2    ]
    [            s  + 2]
    [                  ]
    [    5         5   ]
    [----------    -   ]
    [  / 2    \    1   ]
    [s*\s  + 2/        ]{t}
    >>> tfm_3.var
    s
    >>> tfm_3.shape
    (2, 2)
    >>> tfm_3.num_outputs
    2
    >>> tfm_3.num_inputs
    2
    >>> tfm_3.args
    (((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)), (TransferFunction(5, s*(s**2 + 2), s), TransferFunction(5, 1, s))),)

    To access the ``TransferFunction`` at any index in the ``TransferFunctionMatrix``, use the index notation.

    >>> tfm_3[1, 0]  # gives the TransferFunction present at 2nd Row and 1st Col. Similar to that in Matrix classes
    TransferFunction(5, s*(s**2 + 2), s)
    >>> tfm_3[0, 0]  # gives the TransferFunction present at 1st Row and 1st Col.
    TransferFunction(5, s, s)
    >>> tfm_3[:, 0]  # gives the first column
    TransferFunctionMatrix(((TransferFunction(5, s, s),), (TransferFunction(5, s*(s**2 + 2), s),)))
    >>> pprint(_, use_unicode=False)
    [    5     ]
    [    -     ]
    [    s     ]
    [          ]
    [    5     ]
    [----------]
    [  / 2    \]
    [s*\s  + 2/]{t}
    >>> tfm_3[0, :]  # gives the first row
    TransferFunctionMatrix(((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)),))
    >>> pprint(_, use_unicode=False)
    [5   5*s  ]
    [-  ------]
    [s   2    ]
    [   s  + 2]{t}

    To negate a transfer function matrix, ``-`` operator can be prepended:

    >>> tfm_4 = TransferFunctionMatrix([[tf_2], [-tf_1], [tf_3]])
    >>> -tfm_4
    TransferFunctionMatrix(((TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(-3, s + 2, s),)))
    >>> tfm_5 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, -tf_1]])
    >>> -tfm_5
    TransferFunctionMatrix(((TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(-p**4 + 3*p - 2, p + s, s)), (TransferFunction(-3, s + 2, s), TransferFunction(a + s, s**2 + s + 1, s))))

    ``subs()`` returns the ``TransferFunctionMatrix`` object with the value substituted in the expression. This will not
    mutate your original ``TransferFunctionMatrix``.

    >>> tfm_2.subs(p, 2)  #  substituting p everywhere in tfm_2 with 2.
    TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(12, s + 2, s), TransferFunction(-a - s, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-12, s + 2, s))))
    >>> pprint(_, use_unicode=False)
    [  a + s        -3     ]
    [----------    -----   ]
    [ 2            s + 2   ]
    [s  + s + 1            ]
    [                      ]
    [    12        -a - s  ]
    [  -----     ----------]
    [  s + 2      2        ]
    [            s  + s + 1]
    [                      ]
    [    3          -12    ]
    [  -----       -----   ]
    [  s + 2       s + 2   ]{t}
    >>> pprint(tfm_2, use_unicode=False) # State of tfm_2 is unchanged after substitution
    [   a + s           -3       ]
    [ ----------       -----     ]
    [  2               s + 2     ]
    [ s  + s + 1                 ]
    [                            ]
    [ 4                          ]
    [p  - 3*p + 2      -a - s    ]
    [------------    ----------  ]
    [   p + s         2          ]
    [                s  + s + 1  ]
    [                            ]
    [                 4          ]
    [     3        - p  + 3*p - 2]
    [   -----      --------------]
    [   s + 2          p + s     ]{t}

    ``subs()`` also supports multiple substitutions.

    >>> tfm_2.subs({p: 2, a: 1})  # substituting p with 2 and a with 1
    TransferFunctionMatrix(((TransferFunction(s + 1, s**2 + s + 1, s), TransferFunction(-3, s + 2, s)), (TransferFunction(12, s + 2, s), TransferFunction(-s - 1, s**2 + s + 1, s)), (TransferFunction(3, s + 2, s), TransferFunction(-12, s + 2, s))))
    >>> pprint(_, use_unicode=False)
    [  s + 1        -3     ]
    [----------    -----   ]
    [ 2            s + 2   ]
    [s  + s + 1            ]
    [                      ]
    [    12        -s - 1  ]
    [  -----     ----------]
    [  s + 2      2        ]
    [            s  + s + 1]
    [                      ]
    [    3          -12    ]
    [  -----       -----   ]
    [  s + 2       s + 2   ]{t}

    Users can reduce the ``Series`` and ``Parallel`` elements of the matrix to ``TransferFunction`` by using
    ``doit()``.

    >>> tfm_6 = TransferFunctionMatrix([[Series(tf_3, tf_4), Parallel(tf_3, tf_4)]])
    >>> tfm_6
    TransferFunctionMatrix(((Series(TransferFunction(3, s + 2, s), TransferFunction(-a + p, 9*s - 9, s)), Parallel(TransferFunction(3, s + 2, s), TransferFunction(-a + p, 9*s - 9, s))),))
    >>> pprint(tfm_6, use_unicode=False)
    [-a + p    3    -a + p      3  ]
    [-------*-----  ------- + -----]
    [9*s - 9 s + 2  9*s - 9   s + 2]{t}
    >>> tfm_6.doit()
    TransferFunctionMatrix(((TransferFunction(-3*a + 3*p, (s + 2)*(9*s - 9), s), TransferFunction(27*s + (-a + p)*(s + 2) - 27, (s + 2)*(9*s - 9), s)),))
    >>> pprint(_, use_unicode=False)
    [    -3*a + 3*p     27*s + (-a + p)*(s + 2) - 27]
    [-----------------  ----------------------------]
    [(s + 2)*(9*s - 9)       (s + 2)*(9*s - 9)      ]{t}
    >>> tf_9 = TransferFunction(1, s, s)
    >>> tf_10 = TransferFunction(1, s**2, s)
    >>> tfm_7 = TransferFunctionMatrix([[Series(tf_9, tf_10), tf_9], [tf_10, Parallel(tf_9, tf_10)]])
    >>> tfm_7
    TransferFunctionMatrix(((Series(TransferFunction(1, s, s), TransferFunction(1, s**2, s)), TransferFunction(1, s, s)), (TransferFunction(1, s**2, s), Parallel(TransferFunction(1, s, s), TransferFunction(1, s**2, s)))))
    >>> pprint(tfm_7, use_unicode=False)
    [ 1      1   ]
    [----    -   ]
    [   2    s   ]
    [s*s         ]
    [            ]
    [ 1    1    1]
    [ --   -- + -]
    [  2    2   s]
    [ s    s     ]{t}
    >>> tfm_7.doit()
    TransferFunctionMatrix(((TransferFunction(1, s**3, s), TransferFunction(1, s, s)), (TransferFunction(1, s**2, s), TransferFunction(s**2 + s, s**3, s))))
    >>> pprint(_, use_unicode=False)
    [1     1   ]
    [--    -   ]
    [ 3    s   ]
    [s         ]
    [          ]
    [     2    ]
    [1   s  + s]
    [--  ------]
    [ 2     3  ]
    [s     s   ]{t}

    Addition, subtraction, and multiplication of transfer function matrices can form
    unevaluated ``Series`` or ``Parallel`` objects.

    - For addition and subtraction:
      All the transfer function matrices must have the same shape.

    - For multiplication (C = A * B):
      The number of inputs of the first transfer function matrix (A) must be equal to the
      number of outputs of the second transfer function matrix (B).

    Also, use pretty-printing (``pprint``) to analyse better.

    >>> tfm_8 = TransferFunctionMatrix([[tf_3], [tf_2], [-tf_1]])
    >>> tfm_9 = TransferFunctionMatrix([[-tf_3]])
    >>> tfm_10 = TransferFunctionMatrix([[tf_1], [tf_2], [tf_4]])
    >>> tfm_11 = TransferFunctionMatrix([[tf_4], [-tf_1]])
    >>> tfm_12 = TransferFunctionMatrix([[tf_4, -tf_1, tf_3], [-tf_2, -tf_4, -tf_3]])
    >>> tfm_8 + tfm_10
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a + p, 9*s - 9, s),))))
    >>> pprint(_, use_unicode=False)
    [     3      ]      [   a + s    ]
    [   -----    ]      [ ---------- ]
    [   s + 2    ]      [  2         ]
    [            ]      [ s  + s + 1 ]
    [ 4          ]      [            ]
    [p  - 3*p + 2]      [ 4          ]
    [------------]    + [p  - 3*p + 2]
    [   p + s    ]      [------------]
    [            ]      [   p + s    ]
    [   -a - s   ]      [            ]
    [ ---------- ]      [   -a + p   ]
    [  2         ]      [  -------   ]
    [ s  + s + 1 ]{t}   [  9*s - 9   ]{t}
    >>> -tfm_10 - tfm_8
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(-a - s, s**2 + s + 1, s),), (TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a - p, 9*s - 9, s),))), TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),), (TransferFunction(-p**4 + 3*p - 2, p + s, s),), (TransferFunction(a + s, s**2 + s + 1, s),))))
    >>> pprint(_, use_unicode=False)
    [    -a - s    ]      [     -3       ]
    [  ----------  ]      [    -----     ]
    [   2          ]      [    s + 2     ]
    [  s  + s + 1  ]      [              ]
    [              ]      [   4          ]
    [   4          ]      [- p  + 3*p - 2]
    [- p  + 3*p - 2]    + [--------------]
    [--------------]      [    p + s     ]
    [    p + s     ]      [              ]
    [              ]      [    a + s     ]
    [    a - p     ]      [  ----------  ]
    [   -------    ]      [   2          ]
    [   9*s - 9    ]{t}   [  s  + s + 1  ]{t}
    >>> tfm_12 * tfm_8
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)), (TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)))))
    >>> pprint(_, use_unicode=False)
                                           [     3      ]
                                           [   -----    ]
    [    -a + p        -a - s      3  ]    [   s + 2    ]
    [   -------      ----------  -----]    [            ]
    [   9*s - 9       2          s + 2]    [ 4          ]
    [                s  + s + 1       ]    [p  - 3*p + 2]
    [                                 ]   *[------------]
    [   4                             ]    [   p + s    ]
    [- p  + 3*p - 2    a - p      -3  ]    [            ]
    [--------------   -------    -----]    [   -a - s   ]
    [    p + s        9*s - 9    s + 2]{t} [ ---------- ]
                                           [  2         ]
                                           [ s  + s + 1 ]{t}
    >>> tfm_12 * tfm_8 * tfm_9
    MIMOSeries(TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),),)), TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),))), TransferFunctionMatrix(((TransferFunction(-a + p, 9*s - 9, s), TransferFunction(-a - s, s**2 + s + 1, s), TransferFunction(3, s + 2, s)), (TransferFunction(-p**4 + 3*p - 2, p + s, s), TransferFunction(a - p, 9*s - 9, s), TransferFunction(-3, s + 2, s)))))
    >>> pprint(_, use_unicode=False)
                                           [     3      ]
                                           [   -----    ]
    [    -a + p        -a - s      3  ]    [   s + 2    ]
    [   -------      ----------  -----]    [            ]
    [   9*s - 9       2          s + 2]    [ 4          ]
    [                s  + s + 1       ]    [p  - 3*p + 2]    [ -3  ]
    [                                 ]   *[------------]   *[-----]
    [   4                             ]    [   p + s    ]    [s + 2]{t}
    [- p  + 3*p - 2    a - p      -3  ]    [            ]
    [--------------   -------    -----]    [   -a - s   ]
    [    p + s        9*s - 9    s + 2]{t} [ ---------- ]
                                           [  2         ]
                                           [ s  + s + 1 ]{t}
    >>> tfm_10 + tfm_8*tfm_9
    MIMOParallel(TransferFunctionMatrix(((TransferFunction(a + s, s**2 + s + 1, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a + p, 9*s - 9, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(-3, s + 2, s),),)), TransferFunctionMatrix(((TransferFunction(3, s + 2, s),), (TransferFunction(p**4 - 3*p + 2, p + s, s),), (TransferFunction(-a - s, s**2 + s + 1, s),)))))
    >>> pprint(_, use_unicode=False)
    [   a + s    ]      [     3      ]
    [ ---------- ]      [   -----    ]
    [  2         ]      [   s + 2    ]
    [ s  + s + 1 ]      [            ]
    [            ]      [ 4          ]
    [ 4          ]      [p  - 3*p + 2]    [ -3  ]
    [p  - 3*p + 2]    + [------------]   *[-----]
    [------------]      [   p + s    ]    [s + 2]{t}
    [   p + s    ]      [            ]
    [            ]      [   -a - s   ]
    [   -a + p   ]      [ ---------- ]
    [  -------   ]      [  2         ]
    [  9*s - 9   ]{t}   [ s  + s + 1 ]{t}

    These unevaluated ``Series`` or ``Parallel`` objects can convert into the
    resultant transfer function matrix using ``.doit()`` method or by
    ``.rewrite(TransferFunctionMatrix)``.

    >>> (-tfm_8 + tfm_10 + tfm_8*tfm_9).doit()
    TransferFunctionMatrix(((TransferFunction((a + s)*(s + 2)**3 - 3*(s + 2)**2*(s**2 + s + 1) - 9*(s + 2)*(s**2 + s + 1), (s + 2)**3*(s**2 + s + 1), s),), (TransferFunction((p + s)*(-3*p**4 + 9*p - 6), (p + s)**2*(s + 2), s),), (TransferFunction((-a + p)*(s + 2)*(s**2 + s + 1)**2 + (a + s)*(s + 2)*(9*s - 9)*(s**2 + s + 1) + (3*a + 3*s)*(9*s - 9)*(s**2 + s + 1), (s + 2)*(9*s - 9)*(s**2 + s + 1)**2, s),)))
    >>> (-tfm_12 * -tfm_8 * -tfm_9).rewrite(TransferFunctionMatrix)
    TransferFunctionMatrix(((TransferFunction(3*(-3*a + 3*p)*(p + s)*(s + 2)*(s**2 + s + 1)**2 + 3*(-3*a - 3*s)*(p + s)*(s + 2)*(9*s - 9)*(s**2 + s + 1) + 3*(a + s)*(s + 2)**2*(9*s - 9)*(-p**4 + 3*p - 2)*(s**2 + s + 1), (p + s)*(s + 2)**3*(9*s - 9)*(s**2 + s + 1)**2, s),), (TransferFunction(3*(-a + p)*(p + s)*(s + 2)**2*(-p**4 + 3*p - 2)*(s**2 + s + 1) + 3*(3*a + 3*s)*(p + s)**2*(s + 2)*(9*s - 9) + 3*(p + s)*(s + 2)*(9*s - 9)*(-3*p**4 + 9*p - 6)*(s**2 + s + 1), (p + s)**2*(s + 2)**3*(9*s - 9)*(s**2 + s + 1), s),)))

    See Also
    ========

    DiscreteTransferFunction, TransferFunction, MIMOSeries, MIMOParallel, Feedback

    """
    def __new__(cls, arg):

        expr_mat_arg = []
        try:
            var = arg[0][0].var
        except TypeError:
            raise ValueError(filldedent("""
                `arg` param in TransferFunctionMatrix should
                strictly be a nested list containing TransferFunctionBase
                objects."""))
        for row in arg:
            temp = []
            for element in row:
                if not isinstance(element, SISOLinearTimeInvariant):
                    raise TypeError(filldedent("""
                        Each element is expected to be of
                        type `SISOLinearTimeInvariant`."""))

                if var != element.var:
                    raise ValueError(filldedent("""
                        Conflicting value(s) found for `var`.
                        All TransferFunction instances in TransferFunctionMatrix
                        should use the same complex variable in Laplace domain
                        or z-domain."""))

                temp.append(element.to_expr())
            expr_mat_arg.append(temp)

        _check_time_compatibility([sys for row in arg for sys in row])

        if isinstance(arg, (tuple, list, Tuple)):
            # Making nested Tuple (sympy.core.containers.Tuple) from nested list or nested Python tuple
            arg = Tuple(*(Tuple(*r, sympify=False) for r in arg), sympify=False)

        obj = super(TransferFunctionMatrix, cls).__new__(cls, arg)
        obj._expr_mat = ImmutableMatrix(expr_mat_arg)
        obj.is_StateSpace_object = False
        obj._is_continuous = arg[0][0].is_continuous

        return obj

    @classmethod
    def from_Matrix(cls, matrix, var, sampling_time=0):
        """
        Creates a new ``TransferFunctionMatrix`` efficiently from a SymPy Matrix
        of ``Expr`` objects.

        Parameters
        ==========

        matrix : ``ImmutableMatrix`` having ``Expr``/``Number`` elements.
        var : Symbol
            Complex variable of the Laplace transform or z-transform which will
            be used by the all the transfer function objects in the
            ``TransferFunctionMatrix``.
        sampling_time : Number, Symbol, optional
            Sampling time for the discrete-time transfer function matrix.
            Default is 0, which means that the transfer function matrix will be
            treated as a continuous-time transfer function matrix.

        Examples
        ========

        >>> from sympy.abc import s, z
        >>> from sympy.physics.control.lti import TransferFunctionMatrix
        >>> from sympy import Matrix, pprint
        >>> M = Matrix([[s, 1/s], [1/(s+1), s]])
        >>> M_tf = TransferFunctionMatrix.from_Matrix(M, s)
        >>> pprint(M_tf, use_unicode=False)
        [  s    1]
        [  -    -]
        [  1    s]
        [        ]
        [  1    s]
        [-----  -]
        [s + 1  1]{t}
        >>> M_tf.elem_poles()
        [[[], [0]], [[-1], []]]
        >>> M_tf.elem_zeros()
        [[[0], []], [[], [0]]]
        >>> M_2 = Matrix([[z/(z-1), z/(z-8)], [z**2/(z**2-2+1), z]])
        >>> M2_tf = TransferFunctionMatrix.from_Matrix(M_2, z, 0.1)
        >>> pprint(M2_tf, use_unicode=False)
          [  z       z  ]
          [-----   -----]
          [z - 1   z - 8]
          [             ]
          [   2         ]
          [  z       z  ]
          [------    -  ]
          [ 2        1  ]
          [z  - 1       ]{k}
        [st: 0.100000000000000]



        """
        return _to_TFM(matrix, var, sampling_time)

    @property
    def var(self):
        """
        Returns the complex variable used by all the transfer functions or
        ``Series``/``Parallel`` objects in a transfer function matrix.

        Examples
        ========

        >>> from sympy.abc import p, s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix, Series, Parallel
        >>> G1 = TransferFunction(p**2 + 2*p + 4, p - 6, p)
        >>> G2 = TransferFunction(p, 4 - p, p)
        >>> G3 = TransferFunction(0, p**4 - 1, p)
        >>> G4 = TransferFunction(s + 1, s**2 + s + 1, s)
        >>> S1 = Series(G1, G2)
        >>> S2 = Series(-G3, Parallel(G2, -G1))
        >>> tfm1 = TransferFunctionMatrix([[G1], [G2], [G3]])
        >>> tfm1.var
        p
        >>> tfm2 = TransferFunctionMatrix([[-S1, -S2], [S1, S2]])
        >>> tfm2.var
        p
        >>> tfm3 = TransferFunctionMatrix([[G4]])
        >>> tfm3.var
        s

        """
        return self.args[0][0][0].var

    @property
    def num_inputs(self):
        """
        Returns the number of inputs of the system.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> G1 = TransferFunction(s + 3, s**2 - 3, s)
        >>> G2 = TransferFunction(4, s**2, s)
        >>> G3 = TransferFunction(p**2 + s**2, p - 3, s)
        >>> tfm_1 = TransferFunctionMatrix([[G2, -G1, G3], [-G2, -G1, -G3]])
        >>> tfm_1.num_inputs
        3

        See Also
        ========

        num_outputs

        """
        return self._expr_mat.shape[1]

    @property
    def num_outputs(self):
        """
        Returns the number of outputs of the system.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunctionMatrix
        >>> from sympy import Matrix
        >>> M_1 = Matrix([[s], [1/s]])
        >>> TFM = TransferFunctionMatrix.from_Matrix(M_1, s)
        >>> print(TFM)
        TransferFunctionMatrix(((TransferFunction(s, 1, s),), (TransferFunction(1, s, s),)))
        >>> TFM.num_outputs
        2

        See Also
        ========

        num_inputs

        """
        return self._expr_mat.shape[0]

    @property
    def shape(self):
        """
        Returns the shape of the transfer function matrix, that is,
        ``(# of outputs, # of inputs)``.

        Examples
        ========

        >>> from sympy.abc import s, p
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf1 = TransferFunction(p**2 - 1, s**4 + s**3 - p, p)
        >>> tf2 = TransferFunction(1 - p, p**2 - 3*p + 7, p)
        >>> tf3 = TransferFunction(3, 4, p)
        >>> tfm1 = TransferFunctionMatrix([[tf1, -tf2]])
        >>> tfm1.shape
        (1, 2)
        >>> tfm2 = TransferFunctionMatrix([[-tf2, tf3], [tf1, -tf1]])
        >>> tfm2.shape
        (2, 2)

        """
        return self._expr_mat.shape

    def __neg__(self):
        neg = -self._expr_mat
        return _to_TFM(neg, self.var, self.sampling_time)

    @_check_other_MIMO
    def __add__(self, other):

        if not isinstance(other, MIMOParallel):
            return MIMOParallel(self, other)
        other_arg_list = list(other.args)
        return MIMOParallel(self, *other_arg_list)

    @_check_other_MIMO
    def __sub__(self, other):
        return self + (-other)

    @_check_other_MIMO
    def __mul__(self, other):

        if not isinstance(other, MIMOSeries):
            return MIMOSeries(other, self)
        other_arg_list = list(other.args)
        return MIMOSeries(*other_arg_list, self)

    def __getitem__(self, key):
        trunc = self._expr_mat.__getitem__(key)
        if isinstance(trunc, ImmutableMatrix):
            return _to_TFM(trunc, self.var, self.sampling_time)

        if self.sampling_time == 0:
            to_tf = lambda expr: \
                TransferFunction.from_rational_expression(expr, self.var)
        else:
            to_tf = lambda expr: \
                DiscreteTransferFunction.from_rational_expression(expr, self.var,
                                                            self.sampling_time)
        return to_tf(trunc)

    def transpose(self):
        """
        Returns the transpose of the ``TransferFunctionMatrix``
        (switched input and output layers).

        """
        transposed_mat = self._expr_mat.transpose()
        return _to_TFM(transposed_mat, self.var, self.sampling_time)

    def elem_poles(self):
        """
        Returns the poles of each element of the ``TransferFunctionMatrix``.

        .. note::
            Actual poles of a MIMO system are NOT the poles of individual
            elements.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf_1 = TransferFunction(3, (s + 1), s)
        >>> tf_2 = TransferFunction(s + 6, (s + 1)*(s + 2), s)
        >>> tf_3 = TransferFunction(s + 3, s**2 + 3*s + 2, s)
        >>> tf_4 = TransferFunction(s + 2, s**2 + 5*s - 10, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, tf_4]])
        >>> tfm_1
        TransferFunctionMatrix(((TransferFunction(3, s + 1, s), TransferFunction(s + 6, (s + 1)*(s + 2), s)), (TransferFunction(s + 3, s**2 + 3*s + 2, s), TransferFunction(s + 2, s**2 + 5*s - 10, s))))
        >>> tfm_1.elem_poles()
        [[[-1], [-2, -1]], [[-2, -1], [-5/2 + sqrt(65)/2, -sqrt(65)/2 - 5/2]]]

        See Also
        ========

        elem_zeros

        """
        return [[element.poles() for element in row] for row in \
                                                        self.doit().args[0]]

    def elem_zeros(self):
        """
        Returns the zeros of each element of the ``TransferFunctionMatrix``.

        .. note::
            Actual zeros of a MIMO system are NOT the zeros of individual
            elements.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> tf_1 = TransferFunction(3, (s + 1), s)
        >>> tf_2 = TransferFunction(s + 6, (s + 1)*(s + 2), s)
        >>> tf_3 = TransferFunction(s + 3, s**2 + 3*s + 2, s)
        >>> tf_4 = TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, tf_4]])
        >>> tfm_1
        TransferFunctionMatrix(((TransferFunction(3, s + 1, s), TransferFunction(s + 6, (s + 1)*(s + 2), s)), (TransferFunction(s + 3, s**2 + 3*s + 2, s), TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s))))
        >>> tfm_1.elem_zeros()
        [[[], [-6]], [[-3], [4, 5]]]

        See Also
        ========

        elem_poles

        """
        return [[element.zeros() for element in row] for row in \
                                                        self.doit().args[0]]

    def eval_frequency(self, other):
        """
        Evaluates system response of each transfer function in the
        ``TransferFunctionMatrix`` at any point in the real or complex plane.

        Examples
        ========

        >>> from sympy.abc import s
        >>> from sympy.physics.control.lti import TransferFunction, TransferFunctionMatrix
        >>> from sympy import I
        >>> tf_1 = TransferFunction(3, (s + 1), s)
        >>> tf_2 = TransferFunction(s + 6, (s + 1)*(s + 2), s)
        >>> tf_3 = TransferFunction(s + 3, s**2 + 3*s + 2, s)
        >>> tf_4 = TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s)
        >>> tfm_1 = TransferFunctionMatrix([[tf_1, tf_2], [tf_3, tf_4]])
        >>> tfm_1
        TransferFunctionMatrix(((TransferFunction(3, s + 1, s), TransferFunction(s + 6, (s + 1)*(s + 2), s)), (TransferFunction(s + 3, s**2 + 3*s + 2, s), TransferFunction(s**2 - 9*s + 20, s**2 + 5*s - 10, s))))
        >>> tfm_1.eval_frequency(2)
        Matrix([
        [   1, 2/3],
        [5/12, 3/2]])
        >>> tfm_1.eval_frequency(I*2)
        Matrix([
        [   3/5 - 6*I/5,                -I],
        [3/20 - 11*I/20, -101/74 + 23*I/74]])
        """
        mat = self._expr_mat.subs(self.var, other)
        return mat.expand()

    def _flat(self):
        """Returns flattened list of args in TransferFunctionMatrix"""
        return [elem for tup in self.args[0] for elem in tup]

    def _eval_evalf(self, prec):
        """
        Calls evalf() on each transfer function in the transfer function
        matrix

        """
        dps = prec_to_dps(prec)
        mat = self._expr_mat.applyfunc(lambda a: a.evalf(n=dps))
        return _to_TFM(mat, self.var, self.sampling_time)

    def _eval_simplify(self, **kwargs):
        """Simplifies the transfer function matrix"""
        simp_mat = self._expr_mat.applyfunc(lambda a: cancel(a, expand=False))
        return _to_TFM(simp_mat, self.var, self.sampling_time)

    def expand(self, **hints):
        """Expands the transfer function matrix"""
        expand_mat = self._expr_mat.expand(**hints)
        return _to_TFM(expand_mat, self.var, self.sampling_time)

    @property
    def sampling_time(self):
        return self.args[0][0][0].sampling_time

def create_state_space(A, B, C, D, sampling_time=0):
    """
    Creates a new state space object.
    sampling_time == 0 means continuous time state space.
    sampling_time > 0 means discrete time state space.

    Parameters
    ==========


    sampling_time : Symbol, Number, optional
        Default is 0.
        Time interval between two consecutive sampling instants.
        If sampling_time == 0, it is a continuous time state space,
        else it is a discrete time state space.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.abc import t
    >>> from sympy.physics.control.lti import create_state_space
    >>> A = Matrix([[1,0],[0,1]])
    >>> B = Matrix([1,0])
    >>> C = Matrix([1,0]).T
    >>> D = Matrix([0])
    >>> create_state_space(A, B, C, D)
    StateSpace(Matrix([
    [1, 0],
    [0, 1]]), Matrix([
    [1],
    [0]]), Matrix([[1, 0]]), Matrix([[0]]))
    >>> create_state_space(A, B, C, D, t)
    DiscreteStateSpace(Matrix([
    [1, 0],
    [0, 1]]), Matrix([
    [1],
    [0]]), Matrix([[1, 0]]), Matrix([[0]]), t)

    See Also
    ========

    StateSpace, DiscreteStateSpace

    """
    if sampling_time == 0:
        return StateSpace(A, B, C, D)

    return DiscreteStateSpace(A, B, C, D, sampling_time)

class StateSpaceBase(LinearTimeInvariant, ABC):
    r"""
    Base class for state space objects.
    This class is not meant to be used directly.

    Explanation
    ===========

    State space model (ssm) of a linear, time invariant control system.

    Represents the standard state-space model with A, B, C, D as state-space
    matrices. This makes the linear control system:

    For continuous-time systems:
        (1) x'(t) = A * x(t) + B * u(t);    x in R^n , u in R^k
        (2) y(t)  = C * x(t) + D * u(t);    y in R^m

    For discrete-time systems:
        (1) x[k+1] = A * x[k] + B * u[k];    x in R^n , u in R^k
        (2) y[k]   = C * x[k] + D * u[k];    y in R^m

    where u(t) or u[k] is any input signal, y(t) or y[k] the corresponding
    output, and x(t) or x[k] the system's state.

    Parameters
    ==========

    A : Matrix, optional
        The State matrix of the state space model.
    B : Matrix, optional
        The Input-to-State matrix of the state space model.
    C : Matrix, optional
        The State-to-Output matrix of the state space model.
    D : Matrix, optional
        The Feedthrough matrix of the state space model.
    *args, **kwargs:
        Additional arguments and keyword arguments that are passed to the
        parent class such as sampling time for discrete-time systems.

    See Also
    ========

    StateSpace, DiscreteStateSpace, TransferFunction, DiscreteTransferFunction,
    TransferFunctionMatrix

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/State-space_representation
    .. [2] https://in.mathworks.com/help/control/ref/ss.html

    """
    def __new__(cls, A=None, B=None, C=None, D=None, *args, **kwargs):
        if cls is StateSpaceBase:
            raise NotImplementedError(
                """
                The StateSpaceBase class is not meant to be used directly.
                """)
        if A is None:
            A = zeros(1)
        if B is None:
            B = zeros(A.rows, 1)
        if C is None:
            C = zeros(1, A.cols)
        if D is None:
            D = zeros(C.rows, B.cols)

        A = _sympify(A)
        B = _sympify(B)
        C = _sympify(C)
        D = _sympify(D)

        if (isinstance(A, ImmutableDenseMatrix) and isinstance(B, ImmutableDenseMatrix) and
            isinstance(C, ImmutableDenseMatrix) and isinstance(D, ImmutableDenseMatrix)):
            # Check State Matrix is square
            if A.rows != A.cols:
                raise ShapeError("Matrix A must be a square matrix.")

            # Check State and Input matrices have same rows
            if A.rows != B.rows:
                raise ShapeError("Matrices A and B must have the same number of rows.")

            # Check Output and Feedthrough matrices have same rows
            if C.rows != D.rows:
                raise ShapeError("Matrices C and D must have the same number of rows.")

            # Check State and Output matrices have same columns
            if A.cols != C.cols:
                raise ShapeError("Matrices A and C must have the same number of columns.")

            # Check Input and Feedthrough matrices have same columns
            if B.cols != D.cols:
                raise ShapeError("Matrices B and D must have the same number of columns.")

            obj = super(StateSpaceBase, cls).__new__(cls, A, B, C, D,
                                                     *args, **kwargs)
            obj._A = A
            obj._B = B
            obj._C = C
            obj._D = D

            # Determine if the system is SISO or MIMO
            num_outputs = D.rows
            num_inputs = D.cols
            if num_inputs == 1 and num_outputs == 1:
                obj._is_SISO = True
                obj._clstype = SISOLinearTimeInvariant
            else:
                obj._is_SISO = False
                obj._clstype = MIMOLinearTimeInvariant
            obj.is_StateSpace_object = True
            return obj

        else:
            raise TypeError("A, B, C and D inputs must all be sympy Matrices.")

    @property
    def state_matrix(self):
        """
        Returns the state matrix of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.state_matrix
        Matrix([
        [1, 2],
        [1, 0]])

        """
        return self._A

    @property
    def input_matrix(self):
        """
        Returns the input matrix of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.input_matrix
        Matrix([
        [1],
        [1]])

        """
        return self._B

    @property
    def output_matrix(self):
        """
        Returns the output matrix of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.output_matrix
        Matrix([[0, 1]])

        """
        return self._C

    @property
    def feedforward_matrix(self):
        """
        Returns the feedforward matrix of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.feedforward_matrix
        Matrix([[0]])

        """
        return self._D

    A = state_matrix
    B = input_matrix
    C = output_matrix
    D = feedforward_matrix

    @property
    def num_states(self):
        """
        Returns the number of states of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.num_states
        2

        """
        return self.A.rows

    @property
    def num_inputs(self):
        """
        Returns the number of inputs of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.num_inputs
        1

        """
        return self.D.cols

    @property
    def num_outputs(self):
        """
        Returns the number of outputs of the model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 2], [1, 0]])
        >>> B = Matrix([1, 1])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.num_outputs
        1

        """
        return self.D.rows


    @property
    def shape(self):
        """Returns the shape of the equivalent StateSpace system."""
        return self.num_outputs, self.num_inputs

    def _eval_evalf(self, prec):
        """
        Returns state space model where numerical expressions are evaluated into floating point numbers.
        """
        dps = prec_to_dps(prec)
        return create_state_space(
            self._A.evalf(n = dps),
            self._B.evalf(n = dps),
            self._C.evalf(n = dps),
            self._D.evalf(n = dps),
            self.sampling_time)

    def __add__(self, other):
        """
        Add two State Space systems (parallel connection).

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A1 = Matrix([[1]])
        >>> B1 = Matrix([[2]])
        >>> C1 = Matrix([[-1]])
        >>> D1 = Matrix([[-2]])
        >>> A2 = Matrix([[-1]])
        >>> B2 = Matrix([[-2]])
        >>> C2 = Matrix([[1]])
        >>> D2 = Matrix([[2]])
        >>> ss1 = StateSpace(A1, B1, C1, D1)
        >>> ss2 = StateSpace(A2, B2, C2, D2)
        >>> ss1 + ss2
        StateSpace(Matrix([
        [1,  0],
        [0, -1]]), Matrix([
        [ 2],
        [-2]]), Matrix([[-1, 1]]), Matrix([[0]]))

        """
        # Check for scalars
        if isinstance(other, (int, float, complex, Symbol)):
            A = self.A
            B = self.B
            C = self.C
            D = self.D.applyfunc(lambda element: element + other)
            return create_state_space(A, B, C, D, self.sampling_time)

        # Check nature of system
        if not isinstance(other, StateSpaceBase):
            raise ValueError("Addition is only supported for 2 State Space models.")
        # Check dimensions of system
        elif ((self.num_inputs != other.num_inputs) or (self.num_outputs != other.num_outputs)):
            raise ShapeError("Systems with incompatible inputs and outputs cannot be added.")

        _check_time_compatibility([self, other])

        m1 = (self.A).row_join(zeros(self.A.shape[0], other.A.shape[-1]))
        m2 = zeros(other.A.shape[0], self.A.shape[-1]).row_join(other.A)

        A = m1.col_join(m2)
        B = self.B.col_join(other.B)
        C = self.C.row_join(other.C)
        D = self.D + other.D

        return create_state_space(A, B, C, D, self.sampling_time)

    def __radd__(self, other):
        """
        Right add two State Space systems.

        Examples
        ========

        >>> from sympy.physics.control import StateSpace
        >>> s = StateSpace()
        >>> 5 + s
        StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[5]]))

        """
        return self + other

    def __sub__(self, other):
        """
        Subtract two State Space systems.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A1 = Matrix([[1]])
        >>> B1 = Matrix([[2]])
        >>> C1 = Matrix([[-1]])
        >>> D1 = Matrix([[-2]])
        >>> A2 = Matrix([[-1]])
        >>> B2 = Matrix([[-2]])
        >>> C2 = Matrix([[1]])
        >>> D2 = Matrix([[2]])
        >>> ss1 = StateSpace(A1, B1, C1, D1)
        >>> ss2 = StateSpace(A2, B2, C2, D2)
        >>> ss1 - ss2
        StateSpace(Matrix([
        [1,  0],
        [0, -1]]), Matrix([
        [ 2],
        [-2]]), Matrix([[-1, -1]]), Matrix([[-4]]))

        """
        return self + (-other)

    def __rsub__(self, other):
        """
        Right subtract two tate Space systems.

        Examples
        ========

        >>> from sympy.physics.control import StateSpace
        >>> s = StateSpace()
        >>> 5 - s
        StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[5]]))

        """
        return other + (-self)

    def __neg__(self):
        """
        Returns the negation of the state space model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-5, -1], [3, -1]])
        >>> B = Matrix([2, 5])
        >>> C = Matrix([[1, 2]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> -ss
        StateSpace(Matrix([
        [-5, -1],
        [ 3, -1]]), Matrix([
        [2],
        [5]]), Matrix([[-1, -2]]), Matrix([[0]]))

        """
        return create_state_space(self.A, self.B, -self.C, -self.D,
                               self.sampling_time)

    def __mul__(self, other):
        """
        Multiplication of two State Space systems (serial connection).

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-5, -1], [3, -1]])
        >>> B = Matrix([2, 5])
        >>> C = Matrix([[1, 2]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss*5
        StateSpace(Matrix([
        [-5, -1],
        [ 3, -1]]), Matrix([
        [2],
        [5]]), Matrix([[5, 10]]), Matrix([[0]]))

        """
        # Check for scalars
        if isinstance(other, (int, float, complex, Symbol)):
            A = self.A
            B = self.B
            C = self.C.applyfunc(lambda element: element*other)
            D = self.D.applyfunc(lambda element: element*other)
            return create_state_space(A, B, C, D, self.sampling_time)

        # Check nature of system
        if not isinstance(other, StateSpaceBase):
            raise ValueError("Multiplication is only supported for 2 State Space models.")
        # Check dimensions of system
        elif self.num_inputs != other.num_outputs:
            raise ShapeError("Systems with incompatible inputs and outputs cannot be multiplied.")

        _check_time_compatibility([self, other])

        m1 = (other.A).row_join(zeros(other.A.shape[0], self.A.shape[1]))
        m2 = (self.B * other.C).row_join(self.A)

        A = m1.col_join(m2)
        B = (other.B).col_join(self.B * other.D)
        C = (self.D * other.C).row_join(self.C)
        D = self.D * other.D

        return create_state_space(A, B, C, D, self.sampling_time)

    def __rmul__(self, other):
        """
        Right multiply two tate Space systems.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-5, -1], [3, -1]])
        >>> B = Matrix([2, 5])
        >>> C = Matrix([[1, 2]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> 5*ss
        StateSpace(Matrix([
        [-5, -1],
        [ 3, -1]]), Matrix([
        [10],
        [25]]), Matrix([[1, 2]]), Matrix([[0]]))

        """
        if isinstance(other, (int, float, complex, Symbol)):
            A = self.A
            C = self.C
            B = self.B.applyfunc(lambda element: element*other)
            D = self.D.applyfunc(lambda element: element*other)
            return create_state_space(A, B, C, D, self.sampling_time)
        else:
            return self*other

    def __repr__(self):
        A_str = self.A.__repr__()
        B_str = self.B.__repr__()
        C_str = self.C.__repr__()
        D_str = self.D.__repr__()

        return f"StateSpaceBase(\n{A_str},\n\n{B_str},\n\n{C_str},\n\n{D_str})"

    @_compatibility_decorator
    def append(self, other):
        """
        Returns the first model appended with the second model. The order is preserved.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A1 = Matrix([[1]])
        >>> B1 = Matrix([[2]])
        >>> C1 = Matrix([[-1]])
        >>> D1 = Matrix([[-2]])
        >>> A2 = Matrix([[-1]])
        >>> B2 = Matrix([[-2]])
        >>> C2 = Matrix([[1]])
        >>> D2 = Matrix([[2]])
        >>> ss1 = StateSpace(A1, B1, C1, D1)
        >>> ss2 = StateSpace(A2, B2, C2, D2)
        >>> ss1.append(ss2)
        StateSpace(Matrix([
        [1,  0],
        [0, -1]]), Matrix([
        [2,  0],
        [0, -2]]), Matrix([
        [-1, 0],
        [ 0, 1]]), Matrix([
        [-2, 0],
        [ 0, 2]]))

        """
        n = self.num_states + other.num_states
        m = self.num_inputs + other.num_inputs
        p = self.num_outputs + other.num_outputs

        A = zeros(n, n)
        B = zeros(n, m)
        C = zeros(p, n)
        D = zeros(p, m)

        A[:self.num_states, :self.num_states] = self.A
        A[self.num_states:, self.num_states:] = other.A
        B[:self.num_states, :self.num_inputs] = self.B
        B[self.num_states:, self.num_inputs:] = other.B
        C[:self.num_outputs, :self.num_states] = self.C
        C[self.num_outputs:, self.num_states:] = other.C
        D[:self.num_outputs, :self.num_inputs] = self.D
        D[self.num_outputs:, self.num_inputs:] = other.D
        return create_state_space(A, B, C, D, self.sampling_time)

    def _calc_orthogonal_complement(self, M, dim):
        """
        Returns a basis of the orthogonal complement of a subspace
        represented by the matrix M.
        The orthogonal complement is computed as the null space of the
        transpose of M.

        """
        if M.shape[0] == 0:
            return eye(dim).columnspace()

        return M.T.nullspace()

    def apply_similarity_transform(self, transform_matrix):
        r"""
        Returns an algebrically equivalent state space model, based on the
        transformation matrix `T` such that:

        .. math::
            \begin{cases}
            \bar A=T^{-1}AT\\
            \bar B=T^{-1}B\\
            \bar C=CT\\
            \bar D=D\end{cases}

        Parameters
        ==========

        transform_matrix : Matrix
            The transformation matrix `T` to be applied to the state space
            model.
            The transformation matrix must be invertible and have the same
            dimensions as the state matrix `A`.

        Returns
        =======

        StateSpace
            The transformed state space model.

        Examples
        ========

        >>> from sympy import Matrix, Rational
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[5, -2, 2], [0, 1, -4], [0, 0, -3]])
        >>> B = Matrix([1,0,0])
        >>> C = Matrix([1, 0, 0]).T
        >>> ss = StateSpace(A, B, C)
        >>> T = Matrix([[0, Rational(1, 2), 1], [1, 1, 0], [1, 0, 0]])
        >>> ss.apply_similarity_transform(T).A
            Matrix([
            [-3, 0, 0],
            [ 0, 1, 0],
            [ 0, 0, 5]])

        """
        T_inv = transform_matrix.inv()
        A_decomp = T_inv * self._A * transform_matrix
        B_decomp = T_inv * self._B
        C_decomp = self._C * transform_matrix

        return StateSpace(A_decomp, B_decomp, C_decomp, self._D)

    def to_observable_form(self):
        r"""
        Returns an equivalent state space model decomposed in observable and
        unobservable parts.
        The returned system is algebraically similar to the original but with
        the A and C matrices in block triangular form showing the observable
        and unobservable subsystems.

        .. math::
            \begin{bmatrix}
            A_{O} & 0\\ A_{O\bar O} & A_{\bar O}
            \end{bmatrix}

        .. math::
            \begin{bmatrix}
            B_{O} \\ B_{\bar O}
            \end{bmatrix}

        .. math::
            \begin{bmatrix}
            C_{O} & 0
            \end{bmatrix}

        Examples
        ========

        >>> from sympy import Matrix, Rational
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 0, 1], [0, 0, 0],[0, 0, -2]])
        >>> B = Matrix([1, 1, 0])
        >>> C = Matrix([1, 1, Rational(1,3)]).T
        >>> ss = StateSpace(A, B, C)
        >>> ss = ss.to_observable_form()
        >>> ss.A
        Matrix([
        [0,  0,  0],
        [0,  1,  0],
        [0, -3, -2]])
        >>> ss.B
        Matrix([
        [    1],
        [ 3/10],
        [-3/10]])
        >>> ss.C
        Matrix([[1, 10/3, 0]])

        """
        obs_subsp = Matrix.hstack(*self.observable_subspace())
        unobs_subsp = Matrix.hstack(*self.unobservable_subspace())

        if unobs_subsp.shape[1] == 0:
            # fully observable system
            return self

        if obs_subsp.shape[1] == 0:
            # fully unobservable system
            return self

        T = Matrix.hstack(obs_subsp, unobs_subsp)
        return self.apply_similarity_transform(T)

    def to_controllable_form(self):
        r"""
        Returns an equivalent state space model decomposed in controllable and
        uncontrollable parts.
        The returned system is algebraically similar to the original but with
        the A and B matrices in block triangular form showing the controllable
        and uncontrollable subsystems.

        .. math::
            \begin{bmatrix}
            A_{R} & A_{R\bar R}\\0  & A_{\bar R}
            \end{bmatrix}

        .. math::
            \begin{bmatrix}
            B_{R} \\ 0
            \end{bmatrix}

        .. math::
            \begin{bmatrix}
            C_{R} & C_{\bar R}
            \end{bmatrix}

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[1, 0, 1], [0, 0, 0],[0, 0, -2]])
        >>> B = Matrix([1, 1, 0])
        >>> C = Matrix([1, 1, 0]).T
        >>> ss = StateSpace(A, B, C)
        >>> ss = ss.to_controllable_form()
        >>> ss.A
        Matrix([
        [0, 0,  0],
        [1, 1,  1],
        [0, 0, -2]])
        >>> ss.B
        Matrix([
        [1],
        [0],
        [0]])
        >>> ss.C
        Matrix([[2, 1, 0]])

        """
        contr_subsp = Matrix.hstack(*self.controllable_subspace())
        uncontr_subsp = Matrix.hstack(*self.uncontrollable_subspace())

        if uncontr_subsp.shape[1] == 0:
            # fully controllable system
            return self

        if contr_subsp.shape[1] == 0:
            # fully uncontrollable system
            return self

        T = Matrix.hstack(contr_subsp, uncontr_subsp)
        return self.apply_similarity_transform(T)

    def observability_matrix(self):
        """
        Returns the observability matrix of the state space model:
            [C, C * A^1, C * A^2, .. , C * A^(n-1)]; A in R^(n x n),
            C in R^(m x k)

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> ob = ss.observability_matrix()
        >>> ob
        Matrix([
        [0, 1],
        [1, 0]])

        References
        ==========
        .. [1] https://in.mathworks.com/help/control/ref/statespacemodel.obsv.html

        """
        n = self.num_states
        ob = self.C
        for i in range(1,n):
            ob = ob.col_join(self.C * self.A**i)

        return Matrix(ob)

    def unobservable_subspace(self):
        """
        Returns the unobservable subspace of the state space model.

        Note: this method raises an error for matrices with symbolic entries
        because it is not possible determine the unobservable subspace
        in general.

        Raises
        ======
        NotImplementedError
            If the state space model has symbolic entries, the unobservable
            subspace cannot be determined in general.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> ob_subspace = ss.observable_subspace()
        >>> ob_subspace
        [Matrix([
        [1],
        [0]]), Matrix([
        [0],
        [1]])]
        """
        M = self.observability_matrix()
        if len(M.free_symbols) > 0:
            raise NotImplementedError(
                "Unbservable subspace cannot be determined for " \
                "symbolic matrices."
            )
        return M.nullspace()

    def observable_subspace(self):
        """
        Returns the observable subspace of the state space model.

        Note: this method raises an error for matrices with symbolic entries
        because it is not possible determine the observable subspace
        in general.

        Raises
        ======
        NotImplementedError
            If the state space model has symbolic entries, the observable
            subspace cannot be determined in general.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> ob_subspace = ss.observable_subspace()
        >>> ob_subspace
        [Matrix([
        [1],
        [0]]), Matrix([
        [0],
        [1]])]

        """
        M = Matrix.hstack(*self.unobservable_subspace())
        if len(M.free_symbols) > 0:
                raise NotImplementedError(
                    "Observable subspace cannot be determined for " \
                    "symbolic matrices."
                )
        return self._calc_orthogonal_complement(M, self._A.shape[0])

    def is_observable(self):
        """
        Returns conditions for the state space model to be observable.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A1 = Matrix([[-1.5, -2], [1, 0]])
        >>> B1 = Matrix([0.5, 0])
        >>> C1 = Matrix([[0, 1]])
        >>> D1 = Matrix([1])
        >>> ss1 = StateSpace(A1, B1, C1, D1)
        >>> ss1.is_observable()
        True

        """
        return self.observability_matrix().rank() == self.num_states

    def controllability_matrix(self):
        """
        Returns the controllability matrix of the system:
            [B, A * B, A^2 * B, .. , A^(n-1) * B]; A in R^(n x n), B in R^(n x m)

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.controllability_matrix()
        Matrix([
        [0.5, -0.75],
        [  0,   0.5]])

        References
        ==========
        .. [1] https://in.mathworks.com/help/control/ref/statespacemodel.ctrb.html

        """
        co = self.B
        n = self.A.shape[0]
        for i in range(1, n):
            co = co.row_join(((self.A)**i) * self.B)

        return Matrix(co)

    def uncontrollable_subspace(self):
        """
        Returns the uncontrollable subspace of the state space model.

        Note: this method raises an error for matrices with symbolic entries
        because it is not possible determine the uncontrollable subspace
        in general.

        Raises
        ======
        NotImplementedError
            If the state space model has symbolic entries, the uncontrollable
            subspace cannot be determined in general.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> co_subspace = ss.controllable_subspace()
        >>> co_subspace
        [Matrix([
        [0.5],
        [  0]]), Matrix([
        [-0.75],
        [  0.5]])]

        """
        M = Matrix.hstack(*self.controllable_subspace())
        if len(M.free_symbols) > 0:
            raise NotImplementedError(
                "Uncontrollable subspace cannot be determined for " \
                "symbolic matrices."
            )
        return self._calc_orthogonal_complement(M, self._A.shape[0])

    def controllable_subspace(self):
        """
        Returns the controllable subspace of the state space model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-1.5, -2], [1, 0]])
        >>> B = Matrix([0.5, 0])
        >>> C = Matrix([[0, 1]])
        >>> D = Matrix([1])
        >>> ss = StateSpace(A, B, C, D)
        >>> co_subspace = ss.controllable_subspace()
        >>> co_subspace
        [Matrix([
        [0.5],
        [  0]]), Matrix([
        [-0.75],
        [  0.5]])]

        """
        return self.controllability_matrix().columnspace()

    def is_controllable(self):
        """
        Returns conditions for the state space model to be controllable.

        Examples
        ========

        >>> from sympy import symbols, Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A1 = Matrix([[-1.5, -2], [1, 0]])
        >>> B1 = Matrix([0.5, 0])
        >>> C1 = Matrix([[0, 1]])
        >>> D1 = Matrix([1])
        >>> ss1 = StateSpace(A1, B1, C1, D1)
        >>> ss1.is_controllable()
        True

        >>> a = symbols("a")
        >>> A2 = Matrix([[1, 0, 1], [1, 1, 0], [0, 0, a]])
        >>> B2 = Matrix([0.5, 1, 1])
        >>> C2 = Matrix([[0, 1, 0]])
        >>> ss2 = StateSpace(A2, B2, C2)

        """
        return self.controllability_matrix().rank() == self.num_states

    @abstractmethod
    def get_asymptotic_stability_conditions(self, fast=False) -> list[Boolean]:
        """
        Returns the asymptotic stability conditions for the state space.

        This is a convenient shorthand, based on the system type (continuous or
        discrete), for:

        - ``[c > 0 for c in sys.A.charpoly().hurwitz_conditions()]``
          which gives conditions for stability such that the eigenvalues of
          the A matrix are in the left half of the complex plane.
        - ``[c > 0 for c in sys.A.charpoly().schur_conditions()]``
          which gives conditions for stability such that the eigenvalues of
          the A matrix lie inside the unit circle.

        For some systems that are larger or have more complicated
        expressions it might be useful to set ``fast`` to ``True``.
        The algorithm will use the EXRAW domain which will quickly generate
        large unsimplified expressions that are mostly only suitable for use
        with lambdify.

        Notes
        =====
        - In the discrete case, using ``fast = True`` may lead to significant
          precision issues.

        Parameters
        ==========

        fast : Boolean
            If True, uses the EXRAW domain to quickly generate large
            unsimplified expressions. If False, uses the default domain
            which is suitable for symbolic manipulation.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> from sympy import symbols, reduce_inequalities
        >>> k = symbols('k')
        >>> A = Matrix([[0,1,0],[0,0,1], [k-1, -2*k, -1]])
        >>> B = Matrix([1, 0, 0])
        >>> C = Matrix([[0, 1, 0]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ineq = ss.get_asymptotic_stability_conditions()
        >>> ineq
        [True, 3*k - 1 > 0, 1 - k > 0]
        >>> reduce_inequalities(ineq)
        (1/3 < k) & (k < 1)

        """
        pass

class StateSpace(StateSpaceBase):
    """
    State space model of a linear, continuous-time, time invariant control
    system.

    See py:class:`~.StateSpaceBase` for more details.

    Parameters
    ==========

    A : Matrix, optional
        The State matrix of the state space model.
    B : Matrix, optional
        The Input-to-State matrix of the state space model.
    C : Matrix, optional
        The State-to-Output matrix of the state space model.
    D : Matrix, optional
        The Feedthrough matrix of the state space model.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.physics.control import StateSpace

    The easiest way to create a StateSpaceModel is via four matrices:

    >>> A = Matrix([[1, 2], [1, 0]])
    >>> B = Matrix([1, 1])
    >>> C = Matrix([[0, 1]])
    >>> D = Matrix([0])
    >>> StateSpace(A, B, C, D)
    StateSpace(Matrix([
    [1, 2],
    [1, 0]]), Matrix([
    [1],
    [1]]), Matrix([[0, 1]]), Matrix([[0]]))

    One can use less matrices. The rest will be filled with a minimum of zeros:

    >>> StateSpace(A, B)
    StateSpace(Matrix([
    [1, 2],
    [1, 0]]), Matrix([
    [1],
    [1]]), Matrix([[0, 0]]), Matrix([[0]]))

    See Also
    ========

    StateSpaceBase, DiscreteStateSpace, TransferFunction,
    DiscreteTransferFunction

    """
    def __new__(cls, A=None, B=None, C=None, D=None):
        return super(StateSpace, cls).__new__(cls, A, B, C, D)

    def __repr__(self):
        A_str = self.A.__repr__()
        B_str = self.B.__repr__()
        C_str = self.C.__repr__()
        D_str = self.D.__repr__()

        return f"StateSpace(\n{A_str},\n\n{B_str},\n\n{C_str},\n\n{D_str})"

    def dsolve(self, initial_conditions=None, input_vector=None, var=Symbol('t')):
        r"""
        Returns `y(t)` or output of StateSpace given by the solution of equations:

        .. math::
            \begin{aligned}
            \dot{x}(t) &= Ax(t) + Bu(t) \\
            y(t) &= Cx(t) + Du(t)
            \end{aligned}

        Parameters
        ============

        initial_conditions : Matrix
            The initial conditions of `x` state vector. If not provided, it defaults to a zero vector.
        input_vector : Matrix
            The input vector for state space. If not provided, it defaults to a zero vector.
        var : Symbol
            The symbol representing time. If not provided, it defaults to `t`.

        Examples
        ==========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import StateSpace
        >>> A = Matrix([[-2, 0], [1, -1]])
        >>> B = Matrix([[1], [0]])
        >>> C = Matrix([[2, 1]])
        >>> ip = Matrix([5])
        >>> i = Matrix([0, 0])
        >>> ss = StateSpace(A, B, C)
        >>> ss.dsolve(input_vector=ip, initial_conditions=i).simplify()
        Matrix([[15/2 - 5*exp(-t) - 5*exp(-2*t)/2]])

        If no input is provided it defaults to solving the system with zero initial conditions and zero input.

        >>> ss.dsolve()
        Matrix([[0]])

        References
        ==========
        .. [1] https://web.mit.edu/2.14/www/Handouts/StateSpaceResponse.pdf
        .. [2] https://docs.sympy.org/latest/modules/solvers/ode.html#sympy.solvers.ode.systems.linodesolve

        """

        if not isinstance(var, Symbol):
            raise ValueError("Variable for representing time must be a Symbol.")
        if not initial_conditions:
            initial_conditions = zeros(self._A.shape[0], 1)
        elif initial_conditions.shape != (self._A.shape[0], 1):
            raise ShapeError("Initial condition vector should have the same number of "
                             "rows as the state matrix.")
        if not input_vector:
            input_vector = zeros(self._B.shape[1], 1)
        elif input_vector.shape != (self._B.shape[1], 1):
            raise ShapeError("Input vector should have the same number of "
                             "columns as the input matrix.")
        sol = linodesolve(A=self._A, t=var, b=self._B*input_vector, type='type2', doit=True)
        mat1 = Matrix(sol)
        mat2 = mat1.replace(var, 0)
        free1 = self._A.free_symbols | self._B.free_symbols | input_vector.free_symbols
        free2 = mat2.free_symbols
        # Get all the free symbols form the matrix
        dummy_symbols = list(free2-free1)
        # Convert the matrix to a Coefficient matrix
        r1, r2 = linear_eq_to_matrix(mat2, dummy_symbols)
        s = linsolve((r1, initial_conditions+r2))
        res_tuple = next(iter(s))
        for ind, v in enumerate(res_tuple):
            mat1 = mat1.replace(dummy_symbols[ind], v)
        res = self._C*mat1 + self._D*input_vector
        return res

    def _eval_rewrite_as_TransferFunction(self, *args):
        """
        Returns the equivalent :class:`~.TransferFunction` of the state space
        model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import TransferFunction, StateSpace
        >>> A = Matrix([[-5, -1], [3, -1]])
        >>> B = Matrix([2, 5])
        >>> C = Matrix([[1, 2]])
        >>> D = Matrix([0])
        >>> ss = StateSpace(A, B, C, D)
        >>> ss.rewrite(TransferFunction)
        [[TransferFunction(12*s + 59, s**2 + 6*s + 8, s)]]

        """
        s = Symbol('s')
        n = self.A.shape[0]
        I = eye(n)
        G = self.C*(s*I - self.A).solve(self.B) + self.D
        G = G.simplify()
        to_tf = lambda expr: TransferFunction.from_rational_expression(expr, s)
        tf_mat = [[to_tf(expr) for expr in sublist] for sublist in G.tolist()]
        return tf_mat

    def _eval_rewrite_as_DiscreteTransferFunction(self, *args):
        raise TypeError("""
            The continuous state space model cannot be rewritten as a
            discrete-time transfer function model.
            """)

    def get_asymptotic_stability_conditions(self, fast=False) -> list[Boolean]:
        r"""
        See :func:`StateSpaceBase.get_asymptotic_stability_conditions`.

        """
        s = Symbol('s')

        domain = EXRAW if fast else None
        # if domain is None, to_DM will find the domain automatically
        _A = self.A.to_DM(domain = domain) # type: ignore
        charpoly = _A.charpoly()
        charpoly = Poly(charpoly, s, domain = _A.domain)

        return [c > 0 for c in charpoly.hurwitz_conditions()]

    @property
    def sampling_time(self):
        return S.Zero

    _is_continuous = True

class DiscreteStateSpace(StateSpaceBase):
    """
    State space model of a linear, discrete-time, time invariant control
    system.

    See py:class:`~.StateSpaceBase` for more details.

    Parameters
    ==========

    A : Matrix, optional
        The State matrix of the state space model.
    B : Matrix, optional
        The Input-to-State matrix of the state space model.
    C : Matrix, optional
        The State-to-Output matrix of the state space model.
    D : Matrix, optional
        The Feedthrough matrix of the state space model.
    sampling_time : Symbol, Number, optional
        Time interval between two consecutive sampling instants
        Defaults to 1.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.physics.control import DiscreteStateSpace

    The easiest way to create a StateSpaceModel is via four matrices:

    >>> A = Matrix([[1, 2], [1, 0]])
    >>> B = Matrix([1, 1])
    >>> C = Matrix([[0, 1]])
    >>> D = Matrix([0])
    >>> DiscreteStateSpace(A, B, C, D)
    DiscreteStateSpace(Matrix([
    [1, 2],
    [1, 0]]), Matrix([
    [1],
    [1]]), Matrix([[0, 1]]), Matrix([[0]]), 1)

    One can use less matrices. The rest will be filled with a minimum of zeros:

    >>> DiscreteStateSpace(A, B, sampling_time=0.2)
    DiscreteStateSpace(Matrix([
    [1, 2],
    [1, 0]]), Matrix([
    [1],
    [1]]), Matrix([[0, 0]]), Matrix([[0]]), 0.2)

    See Also
    ========

    StateSpaceBase, StateSpace, TransferFunction,
    DiscreteTransferFunction

    """
    def __new__(cls, A=None, B=None, C=None, D=None, sampling_time=1):
        if sampling_time == 0:
            raise ValueError(filldedent("""
                The sampling time cannot be zero.
                If you want to create a continuous state space,
                use the StateSpace class instead."""))

        sampling_time = sympify(sampling_time)
        obj = super(DiscreteStateSpace, cls).__new__(cls, A, B, C, D, sampling_time)
        obj._sampling_time = sampling_time

        return obj

    def __repr__(self):
        A_str = self.A.__repr__()
        B_str = self.B.__repr__()
        C_str = self.C.__repr__()
        D_str = self.D.__repr__()

        return f"""DiscreteStateSpace(\n{A_str},
        \n{B_str},
        \n{C_str},
        \n{D_str},
        \nst: {self.sampling_time})"""

    def _eval_rewrite_as_DiscreteTransferFunction(self, *args):
        """
        Returns the equivalent :class:`~.DiscreteTransferFunction` of the state
        space model.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.physics.control import DiscreteTransferFunction, DiscreteStateSpace
        >>> A = Matrix([[-5, -1], [3, -1]])
        >>> B = Matrix([2, 5])
        >>> C = Matrix([[1, 2]])
        >>> D = Matrix([0])
        >>> ss = DiscreteStateSpace(A, B, C, D)
        >>> ss.rewrite(DiscreteTransferFunction)
        [[DiscreteTransferFunction(12*z + 59, z**2 + 6*z + 8, z, 1)]]

        """
        z = Symbol('z')
        n = self.A.shape[0]
        I = eye(n)
        G = self.C*(z*I - self.A).solve(self.B) + self.D
        G = G.simplify()
        to_tf = lambda expr: DiscreteTransferFunction.\
            from_rational_expression(expr, z, self.sampling_time)
        tf_mat = [[to_tf(expr) for expr in sublist] for sublist in G.tolist()]
        return tf_mat

    def _eval_rewrite_as_TransferFunction(self, *args):
        raise TypeError("""
            The discrete state space model cannot be rewritten as a
            continuous-time transfer function model.
            """)

    def get_asymptotic_stability_conditions(
        self, fast=False
    ) -> list[Boolean]:
        r"""
        See :func:`StateSpaceBase.get_asymptotic_stability_conditions`.

        """
        s = Symbol('s')

        domain = EXRAW if fast else None
        # if domain is None, to_DM will find the domain automatically
        _A = self.A.to_DM(domain = domain) # type: ignore
        charpoly = _A.charpoly()
        charpoly = Poly(charpoly, s, domain = _A.domain)

        return [c > 0 for c in charpoly.schur_conditions()]

    @property
    def sampling_time(self):
        return self._sampling_time

    _is_continuous = False
