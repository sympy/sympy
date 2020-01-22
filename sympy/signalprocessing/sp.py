from sympy.simplify.radsimp import fraction
from sympy.solvers import solve
from sympy.core.numbers import oo, I
from sympy.core.symbol import Symbol
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.complexes import Abs, arg, re, im
from sympy.utilities.lambdify import lambdify
from sympy.core.singleton import S
from sympy.core.compatibility import ordered
import math
import matplotlib.pyplot as plt
import numpy as np


# Function for obtaining the poles and zeros of a transfer function
# Credits scmichr
def pole_zero(expression, symbol):
    num, den = expression.as_numer_denom()
    poles = set(solve(den, symbol))
    zeros = set(solve(num, symbol))
  # keep unique elements
    common = poles & zeros
    poles = list(ordered(poles - common))
    zeros = list(ordered(zeros - common))
    poles_extra = len(poles) - len(zeros)
  # make them the same length
    poles += [S.Zero]*-poles_extra
    zeros += [S.Infinity]*poles_extra
    return poles, zeros

# Function for obtaining the poles and zeros of a transfer function
def pole_zero_plot(expression, symbol):
    poles, zeros = np.asarray(pole_zero(expression, symbol)[0]), np.asarray(pole_zero(expression, symbol)[1])
    circle_points = np.linspace(0.1, 200, 200)
    circle_points = math.pi / 100 * circle_points
    plt.figure(figsize=(5, 5))
    plt.plot(np.cos(circle_points), np.sin(circle_points))
    x_poles, y_poles, x_zeros, y_zeros = [], [], [], []
    for i in range(len(poles)):
        x_poles.append(re(poles[i]))
        y_poles.append(im(poles[i]))
        x_zeros.append(re(zeros[i]))
        y_zeros.append(im(zeros[i]))
    plt.plot(x_zeros, y_zeros, 'o', label = 'Zeros')
    plt.plot(x_poles, y_poles, 'x', label = 'Poles')
    plt.legend()
    plt.show()


# Function for obtaining the Magnitude and Phase Response of a transfer function
def magnitude_and_phase_response(expression, symbol):
    w = Symbol('w')
    expr = expression.subs(symbol, exp(I*w))
    magnitude_resp = Abs(expr)
    phase_resp = arg(expr)
    return magnitude_resp, phase_resp, w


# Function for plotting the Magnitude Response
def plot_mag(expression,  symbol):
    magnitude_resp, _, w = magnitude_and_phase_response(expression, symbol)
    x_val = np.linspace(0.1, 100, 100)
    x_val = math.pi/100 * x_val
    f = lambdify(w, magnitude_resp, "numpy")
    y_val = f(x_val)
    plt.plot(x_val, y_val)
    plt.title('Magnitude Response')
    plt.ylabel('Magnitude')
    plt.xlabel('Freq (in Rad/s)')
    plt.show()


# Function for plotting the Phase Response
def plot_phase(expression,  symbol):
    _, phase_resp, w = magnitude_and_phase_response(expression, symbol)
    x_val = np.linspace(0.1, 100, 100)
    x_val = math.pi/100 * x_val
    f = lambdify(w, phase_resp, "numpy")
    y_val = f(x_val)
    plt.plot(x_val, y_val)
    plt.title('Phase Response')
    plt.ylabel('Phase')
    plt.xlabel('Freq (in Rad/s)')
    plt.show()
