from sympy import (I, log, sqrt, symbols, apart, Wild,
    RootSum, Lambda, together, exp, gamma)
from sympy.core.symbol import Dummy
from sympy.external import import_module
from sympy.functions import arg
from sympy.physics.control.lti import SISOLinearTimeInvariant
from sympy.plotting.plot import LineOver1DRangeSeries
from sympy.polys.polytools import Poly
from sympy.printing.latex import latex

matplotlib = import_module(
        'matplotlib', import_kwargs={'fromlist': ['pyplot']},
        catch=(RuntimeError,))

numpy = import_module('numpy')

if matplotlib:
    plt = matplotlib.pyplot

if numpy:
    np = numpy  # Matplotlib already has numpy as a compulsory dependency. No need to install it separately.


def _check_system(system):
    """Function to check whether the dynamical system passed for plots is
    compatible or not."""
    if not isinstance(system, SISOLinearTimeInvariant):
        raise NotImplementedError("Only SISO LTI systems are currently supported.")
    sys = system.to_expr()
    len_free_symbols = len(sys.free_symbols)
    if len_free_symbols > 1:
        raise ValueError("Extra degree of freedom found. Make sure"
            " that there are no free symbols in the dynamical system other"
            " than the variable of Laplace transform.")
    if sys.has(exp):
        raise NotImplementedError("Time delay terms are not supported.")

    return


def pole_zero_numerical_data(system):
    """Returns the numerical data of poles and zeros of the system."""
    _check_system(system)
    system = system.doit()  # Get the equivalent TransferFunction object.

    num_poly = Poly(system.num, system.var).all_coeffs()
    den_poly = Poly(system.den, system.var).all_coeffs()

    num_poly = np.array(num_poly, dtype=np.float64)
    den_poly = np.array(den_poly, dtype=np.float64)

    zeros = np.roots(num_poly)
    poles = np.roots(den_poly)

    return zeros, poles


def pole_zero_plot(system, pole_color='blue', pole_markersize=10,
    zero_color='orange', zero_markersize=7, grid=True, show_axes=True,
    show=True, **kwargs):
    r"""
    Returns the Pole-Zero plot (also known as PZ Plot or PZ Map) of a system.
    """
    zeros, poles = pole_zero_numerical_data(system)

    zero_real = np.real(zeros)
    zero_imag = np.imag(zeros)

    pole_real = np.real(poles)
    pole_imag = np.imag(poles)

    plt.plot(pole_real, pole_imag, 'x', mfc='none',
        markersize=pole_markersize, color=pole_color)
    plt.plot(zero_real, zero_imag, 'o', markersize=zero_markersize,
        color=zero_color)
    plt.xlabel('Real Axis')
    plt.ylabel('Imaginary Axis')
    plt.title(f'Poles and Zeros of ${latex(system)}$', pad=20)

    if grid:
        plt.grid()
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def step_response_numerical_data(system, prec=8, lower_limit=0,
    upper_limit=10, **kwargs):
    if lower_limit < 0:
        raise ValueError("Lower limit of time must be greater "
            "than or equal to zero.")
    _check_system(system)
    _x = Dummy("x")
    expr = system.to_expr()/(system.var)
    expr = apart(expr, system.var, full=True)
    _y = _fast_inverse_laplace(expr, system.var, _x).evalf(prec)
    return LineOver1DRangeSeries(_y, (_x, lower_limit, upper_limit),
        **kwargs).get_points()


def step_response_plot(system, color='b', show=True, lower_limit=0,
    upper_limit=10, prec=8, show_axes=False, grid=True, **kwargs):
    r"""
    Returns the unit step response of a continuous-time system. It is
    the response of the system when the input is a step input.
    """
    x, y = step_response_numerical_data(system, prec=prec,
        lower_limit=lower_limit, upper_limit=upper_limit, **kwargs)
    plt.plot(x, y, color=color)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Unit Step Response of ${latex(system)}$', pad=20)

    if grid:
        plt.grid()
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def impulse_response_numerical_data(system, prec=8, lower_limit=0,
    upper_limit=10, **kwargs):
    if lower_limit < 0:
        raise ValueError("Lower limit of time must be greater "
            "than or equal to zero.")
    _check_system(system)
    _x = Dummy("x")
    expr = system.to_expr()
    expr = apart(expr, system.var, full=True)
    _y = _fast_inverse_laplace(expr, system.var, _x).evalf(prec)
    return LineOver1DRangeSeries(_y, (_x, lower_limit, upper_limit),
        **kwargs).get_points()


def impulse_response_plot(system, color='b', show=True, lower_limit=0,
    upper_limit=10, prec=8, show_axes=False, grid=True, **kwargs):
    r"""
    Returns the unit impulse response (Input is the Dirac-Delta Function) of a
    continuous-time system.
    """
    x, y = impulse_response_numerical_data(system, prec=prec,
        lower_limit=lower_limit, upper_limit=upper_limit, **kwargs)
    plt.plot(x, y, color=color)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Unit Step Response of ${latex(system)}$', pad=20)

    if grid:
        plt.grid()
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def ramp_response_numerical_data(system, slope=1, prec=8,
    lower_limit=0, upper_limit=10, **kwargs):
    if slope < 0:
        raise ValueError("Slope must be greater than or equal"
            " to zero.")
    if lower_limit < 0:
        raise ValueError("Lower limit of time must be greater "
            "than or equal to zero.")
    _check_system(system)
    _x = Dummy("x")
    expr = (slope*system.to_expr())/((system.var)**2)
    expr = apart(expr, system.var, full=True)
    _y = _fast_inverse_laplace(expr, system.var, _x).evalf(prec)
    return LineOver1DRangeSeries(_y, (_x, lower_limit, upper_limit),
        **kwargs).get_points()


def ramp_response_plot(system, slope=1, color='b', show=True, lower_limit=0,
    upper_limit=10, prec=8, show_axes=False, grid=True, **kwargs):
    r"""
    Returns the ramp response of a continuous-time system. Ramp
    function is defined as the the straight line passing through
    origin ($$f(x) = mx$$). The slope of the ramp function can be
    varied by the user and the default value is 1.
    """
    x, y = ramp_response_numerical_data(system, slope=slope, prec=prec,
        lower_limit=lower_limit, upper_limit=upper_limit, **kwargs)
    plt.plot(x, y, color=color)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Ramp Response of ${latex(system)}$ [Slope = {slope}]', pad=20)

    if grid:
        plt.grid()
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def bode_magnitude_numerical_data(system, initial_exp=-5, final_exp=5):
    _check_system(system)
    expr = system.to_expr()
    _w = Dummy("w", real=True)
    w_expr = expr.subs({system.var: I*_w})
    real, imag = w_expr.as_real_imag()

    mag = 20*log(sqrt(real**2 + imag**2), 10)

    return LineOver1DRangeSeries(mag,
        (_w, 10**initial_exp, 10**final_exp), xscale='log').get_points()


def bode_magnitude_plot(system, initial_exp=-5, final_exp=5,
    color='b', show=True, show_axes=False, grid=True, **kwargs):
    r"""
    Returns the Bode magnitude plot of a continuous-time system.
    """
    x, y = bode_magnitude_numerical_data(system, initial_exp=initial_exp,
        final_exp=final_exp)
    plt.plot(x, y, color=color, **kwargs)
    plt.xscale('log')

    plt.xlabel('Frequency (Hz) [Log Scale]')
    plt.ylabel('Magnitude (dB)')
    plt.title(f'Bode Plot (Magnitude) of ${latex(system)}$', pad=20)

    if grid:
        plt.grid(True)
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def bode_phase_numerical_data(system, initial_exp=-5, final_exp=5):
    _check_system(system)
    expr = system.to_expr()
    _w = Dummy("w", real=True)
    w_expr = expr.subs({system.var: I*_w})

    phase = arg(w_expr)

    return LineOver1DRangeSeries(phase,
        (_w, 10**initial_exp, 10**final_exp), xscale='log').get_points()


def bode_phase_plot(system, initial_exp=-5, final_exp=5,
    color='b', show=True, show_axes=False, grid=True, **kwargs):
    r"""
    Returns the Bode phase plot of a continuous-time system.
    """
    x, y = bode_phase_numerical_data(system, initial_exp=initial_exp,
        final_exp=final_exp)
    plt.plot(x, y, color=color, **kwargs)
    plt.xscale('log')

    plt.xlabel('Frequency (Hz) [Log Scale]')
    plt.ylabel('Phase (rad)')
    plt.title(f'Bode Plot (Phase) of ${latex(system)}$', pad=20)

    if grid:
        plt.grid(True)
    if show_axes:
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')
    if show:
        plt.show()
        return

    return plt


def bode_plot(system, initial_exp=-5, final_exp=5, show=True, grid=True, show_axes=False, **kwargs):
    r"""
    Returns the Bode plot of a continuous-time system.
    """
    plt.subplot(211)
    bode_magnitude_plot(system, initial_exp=initial_exp, final_exp=final_exp,
        show=False, grid=grid, show_axes=show_axes,
        **kwargs).title(f'Bode Plot of ${latex(system)}$', pad=20)
    plt.subplot(212)
    bode_phase_plot(system, initial_exp=initial_exp, final_exp=final_exp,
        show=False, grid=grid, show_axes=show_axes, **kwargs).title(None)

    if show:
        plt.show()
        return

    return plt


def _fast_inverse_laplace(e, s, t):
    """Fast inverse Laplace transform of rational function including RootSum"""
    a, b, n = symbols('a, b, n', cls=Wild, exclude=[s])

    def _ilt(e):
        if not e.has(s):
            return e
        elif e.is_Add:
            return _ilt_add(e)
        elif e.is_Mul:
            return _ilt_mul(e)
        elif e.is_Pow:
            return _ilt_pow(e)
        elif isinstance(e, RootSum):
            return _ilt_rootsum(e)
        else:
            raise NotImplementedError

    def _ilt_add(e):
        return e.func(*map(_ilt, e.args))

    def _ilt_mul(e):
        coeff, expr = e.as_independent(s)
        if expr.is_Mul:
            raise NotImplementedError
        return coeff * _ilt(expr)

    def _ilt_pow(e):
        match = e.match((a*s + b)**n)
        if match is not None:
            nm, am, bm = match[n], match[a], match[b]
            if nm.is_Integer and nm < 0:
                if nm == 1:
                    return exp(-(bm/am)*t) / am
                else:
                    return t**(-nm-1)*exp(-(bm/am)*t)/(am**-nm*gamma(-nm))
        raise NotImplementedError

    def _ilt_rootsum(e):
        expr = e.fun.expr
        [variable] = e.fun.variables
        return RootSum(e.poly, Lambda(variable, together(_ilt(expr))))

    return _ilt(e)
