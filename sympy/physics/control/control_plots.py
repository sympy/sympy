from sympy import I, log, sqrt
from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.functions import arg
from sympy.integrals.transforms import inverse_laplace_transform
from sympy.plotting import PlotGrid, plot
from sympy.polys.polytools import Poly

matplotlib = import_module(
        'matplotlib', import_kwargs={'fromlist': ['pyplot']},
        catch=(RuntimeError,))

numpy = import_module('numpy')

if matplotlib:
    plt = matplotlib.pyplot
else:
    raise ImportError("Matplotlib is required for plotting (as an external dependency).")

if numpy:
    np = numpy  # Matplotlib already has numpy as a compulsory dependency. No need to install it separately.
else:
    raise ImportError("NumPy is required for this function.")


def pole_zero(system, pole_colour='r', zero_colour='b', grid=True, show=True, **kwargs):
    r"""
    Returns the Pole-Zero plot (also known as PZ Plot or PZ Map) of a system.
    """
    system = system.doit()  # Get the equivalent TransferFunction object.

    num_poly = Poly(system.num, system.var).all_coeffs()
    den_poly = Poly(system.den, system.var).all_coeffs()

    num_poly = np.array(num_poly, dtype=np.float64)
    den_poly = np.array(den_poly, dtype=np.float64)

    zeros, poles = _find_roots(num_poly, den_poly)

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


def _find_roots(num, den, k=None):
    """Private method to find roots efficiently."""
    if k is not None:
        num_roots = max(len(np.roots(den)), len(np.roots(num)))
        roots = []
        for gain in k:
            individual_roots = np.roots(np.polyadd(den, (gain*num)))
            if len(individual_roots) == num_roots:
                individual_roots.sort()
                roots.append(individual_roots)

        return np.vstack(roots)

    num_roots = np.roots(num)
    den_roots = np.roots(den)

    return num_roots, den_roots


def root_locus(system, k_max=400.0, num=10000, show=True, grid=False, **kwargs):
    """
    Generates a root locus plot of a continuous-time system.
    """
    system = system.doit()

    k = np.linspace(0.0, float(k_max), num=num)

    num_poly = Poly(system.num, system.var).all_coeffs()
    den_poly = Poly(system.den, system.var).all_coeffs()

    num_poly = np.array(num_poly, dtype=np.float64)
    den_poly = np.array(den_poly, dtype=np.float64)

    roots = _find_roots(num_poly, den_poly, k)

    real = np.real(roots)
    imag = np.imag(roots)

    pz = pole_zero(system, show=False, grid=grid)

    pz.title('Root Loci')
    pz.plot(real, imag)

    if show:
        pz.show()
        return
    return pz
