from sympy import I, log, sqrt
from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.functions import arg
from sympy.integrals.transforms import inverse_laplace_transform
from sympy.plotting import PlotGrid, plot

matplotlib = import_module(
        'matplotlib', import_kwargs={'fromlist': ['pyplot']},
        catch=(RuntimeError,))

if matplotlib:
    plt = matplotlib.pyplot


def pole_zero(transfer_function, pole_colour='r', zero_colour='b', grid=True, **kwargs):
    r"""
    Returns the Pole-Zero plot (also known as PZ Plot or PZ Map) of a system.
    """
    poles = transfer_function.poles()
    zeros = transfer_function.zeros()

    pole_points = [pole.as_real_imag() for pole in poles]
    zero_points = [zero.as_real_imag() for zero in zeros]

    x_poles = list(map(lambda x: x[0], pole_points))
    y_poles = list(map(lambda x: x[1], pole_points))

    x_zeros = list(map(lambda x: x[0], zero_points))
    y_zeros = list(map(lambda x: x[1], zero_points))

    plt.plot(x_poles, y_poles, 'x', mfc='none', markersize=10)
    plt.plot(x_zeros, y_zeros, 'o', markersize=7)
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Poles and Zeros')
    if grid:
        plt.grid()
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.show()
    return


def step_response(system, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit step response of a continuous-time system.
    """
    x = Symbol("x")
    expr = system.to_expr()/(system.var)
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, upper_limit), show=show, title="Unit Step Response", \
        xlabel="Time (s)", ylabel="Amplitude")


def impulse_response(system, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit impulse response (Input is the Dirac-Delta Function) of a
    continuous-time system.
    """
    x = Symbol("x")
    expr = system.to_expr()
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, 6), show=True, title="Impulse Response", \
        xlabel="Time (s)", ylabel="Amplitude")


def ramp_response(system, slope=1, show_input=False, colour='r', show=True, grid=True, upper_limit=6, **kwargs):
    r"""
    Return the unit impulse response (Input is the Dirac-Delta Function) of a
    continuous-time system.
    """
    x = Symbol("x")
    expr = (slope*system.to_expr())/((system.var)**2)
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, upper_limit), show=show, title="Ramp Response", \
        xlabel="Time (s)", ylabel="Amplitude")


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
