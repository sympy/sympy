from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.integrals.transforms import inverse_laplace_transform
from sympy.plotting import plot

matplotlib = import_module(
        'matplotlib', import_kwargs={'fromlist': ['pyplot']},
        catch=(RuntimeError,))

if matplotlib:
    plt = matplotlib.pyplot

def pole_zero(transfer_function, **kwargs):
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

    plt.plot(x_poles, y_poles, 'x', mfc='none', markersize=15)
    plt.plot(x_zeros, y_zeros, 'o', markersize=10)
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Poles and Zeros')
    plt.grid()
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.show()
    return

def step_response(system, **kwargs):
    r"""
    Return the unit step response of a continuous-time system.
    """
    x = Symbol("x")
    expr = system.to_expr()/(system.var)
    y = inverse_laplace_transform(expr, system.var, x)
    return plot(y, (x, 0, 6), show=True, title="Unit Step Response", \
        xlabel="Time (Seconds)", ylabel="Amplitude")
