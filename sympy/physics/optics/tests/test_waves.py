from sympy import symbols, Symbol, pi, sqrt, cos
from sympy.physics.optics import TWave


def test_twave():
    A1, phi1, A2, phi2, f = symbols('A1, phi1, A2, phi2, f')
    c = Symbol('c')  # Speed of light in vacuum
    n = Symbol('n')  # Refractive index
    t = Symbol('t')  # Time
    w1 = TWave(A1, f, phi1)
    w2 = TWave(A2, f, phi2)
    assert w1.amplitude == A1
    assert w1.frequency == f
    assert w1.phase == phi1
    assert w1.wavelength == c/(f*n)
    assert w1.time_period == 1/f
    w3 = w1 + w2
    assert w3.amplitude == sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2)
    assert w3.frequency == f
    assert w3.wavelength == c/(f*n)
    assert w3.time_period == 1/f
    assert w3.angular_velocity == 2*pi*f
    assert w3.equation() == sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2)*cos(2*pi*f*t + phi1 + phi2)
