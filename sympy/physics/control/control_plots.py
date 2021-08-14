from sympy import (I, log, sqrt, symbols, apart, Wild,
    RootSum, Lambda, together, exp, gamma)
from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.functions import arg
from sympy.integrals.transforms import inverse_laplace_transform
from sympy.physics.control.lti import SISOLinearTimeInvariant
from sympy.plotting import PlotGrid, plot
from sympy.plotting.plot import LineOver1DRangeSeries
from sympy.polys.polytools import Poly

matplotlib = import_module(
        'matplotlib', import_kwargs={'fromlist': ['pyplot']},
        catch=(RuntimeError,))

numpy = import_module('numpy')

if matplotlib:
    plt = matplotlib.pyplot

if numpy:
    np = numpy  # Matplotlib already has numpy as a compulsory dependency. No need to install it separately.


def pole_zero_numerical_data(system):
    """Returns the numerical data of poles and zeros of the system."""
    system = system.doit()  # Get the equivalent TransferFunction object.
    if not isinstance(system, SISOLinearTimeInvariant):
        
    sys = system.to_expr()
    len_free_symbols = len(sys.free_symbols)
    if len_free_symbols > 1:
        raise ValueError("Extra degree of freedom found. Make sure"
            "There are no free symbols in the dynamical system other"
            " than the variable of Laplace transform.")

    num_poly = Poly(system.num, system.var).all_coeffs()
    den_poly = Poly(system.den, system.var).all_coeffs()

    num_poly = np.array(num_poly, dtype=np.float64)
    den_poly = np.array(den_poly, dtype=np.float64)

    zeros = np.roots(num_poly)
    poles = np.roots(den_poly)

    return zeros, poles


def pole_zero_plot(system, pole_colour='r', zero_colour='b', grid=True, show=True, **kwargs):
    r"""
    Returns the Pole-Zero plot (also known as PZ Plot or PZ Map) of a system.
    """
    zeros, poles = pole_zero_numerical_data(system)

    zero_real = np.real(zeros)
    zero_imag = np.imag(zeros)

    pole_real = np.real(poles)
    pole_imag = np.imag(poles)

    plt.plot(pole_real, pole_imag, 'x', mfc='none', markersize=10)
    plt.plot(zero_real, zero_imag, 'o', markersize=7)
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Poles and Zeros')
    if grid:
        plt.grid()
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    if show:
        plt.show()
        return
    return plt


def step_response_numerical_data(system):
    pass


def step_response_plot(system, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit step response of a continuous-time system.
    """
    x = Symbol("x")
    expr = system.to_expr()/(system.var)
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, upper_limit), show=show, title="Unit Step Response", \
        xlabel="Time (s)", ylabel="Amplitude")


def impulse_response_numerical_data(system):
    pass


def impulse_response_plot(system, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit impulse response (Input is the Dirac-Delta Function) of a
    continuous-time system.
    """
    x = Symbol("x")
    expr = system.to_expr()
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, 60), show=True, title="Impulse Response", \
        xlabel="Time (s)", ylabel="Amplitude")


def ramp_response_numerical_data(system):
    pass


def ramp_response_plot(system, slope=1, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit impulse response (Input is the Dirac-Delta Function) of a
    continuous-time system.
    """
    x = Symbol("x")
    expr = (slope*system.to_expr())/((system.var)**2)
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, upper_limit), show=show, title="Ramp Response", \
        xlabel="Time (s)", ylabel="Amplitude")


def bode_numerical_data(system):
    pass


def bode_plot(system, initial_exp=-5, final_exp=5, show=True, **kwargs):
    r"""
    Plot a Bode plot of a continuous-time system.
    """
    expr = system.to_expr()
    w = Symbol("w", real=True)
    w_expr = expr.subs({system.var: I*w})
    real, imag = w_expr.as_real_imag()

    mag = 20*log(sqrt(real**2 + imag**2), 10)

    mag_plot = plot(mag, (w, 10**initial_exp, 10**final_exp), xscale='log', title="Bode Plot (Magnitude)", \
        xlabel="Frequency (Hz) [Log Scale]", ylabel="Magnitude (dB)", show=False)

    phase = arg(w_expr)

    phase_plot = plot(phase, (w, 10**initial_exp, 10**final_exp), xscale='log', title="Bode Plot (Phase)", \
        xlabel="Frequency (Hz) [Log Scale]", ylabel="Phase (rad)", show=False)

    return PlotGrid(2, 1, mag_plot, phase_plot, show=show)


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
