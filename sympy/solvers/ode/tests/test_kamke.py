"""
File to test Kamke ODEs. To run test suite, run

    python test_kamke.py

Information about failing ODEs is logged in test_report.txt.
"""


from contextlib import contextmanager
import threading
import _thread
import time
import os

from sympy import (symbols, Function, Eq, dsolve, checkodesol,
    classify_ode, latex, Abs, sqrt, cos, exp, log, sin, tan, Derivative,
    Ci, Ei, Piecewise, Si, cosh, cot, coth, sinh, tanh, Integral, I, Subs,
    diff, Ne, Sum, IndexedBase, Pow)
from sympy.core.function import AppliedUndef
from sympy.solvers.deutils import ode_order

class TimeOutError(Exception):
    def __init__(self, msg=''):
        self.msg = msg


@contextmanager
def time_limit(seconds, msg=''):
    timer = threading.Timer(seconds, lambda: _thread.interrupt_main())
    timer.start()
    try:
        yield
    except KeyboardInterrupt:
        raise TimeOutError("Timed out for operation: {}".format(msg))
    finally:
        # Cancel timer if process finishes in time
        timer.cancel()

A, A0, A1, A2, B, C, DD, EE, a, a0, a1, a11, a12, a2, a21, a22, a3, a4, al1, al2, al3, alpha, b, b0, b1, b11, b12, b2, b21, b22, b3, b4, bbeta, bl1, bl2, bl3, c, c0, c1, c11, c12, c2, c21, c22, c3, c4, d, d1, ddf, ddg, delta, df, dg, e, e1, e2, fs, f1, f2, f3, gs, g2, g3, ggamma, hs, k, l, lambda_, m, mu, n, nu, omega, ps, qs, r, rho, s, sigma, t, v, x, x4, xp, z = symbols('A, A0, A1, A2, B, C, DD, EE, a, a0, a1, a11, a12, a2, a21, a22, a3, a4, al1, al2, al3, alpha, b, b0, b1, b11, b12, b2, b21, b22, b3, b4, bbeta, bl1, bl2, bl3, c, c0, c1, c11, c12, c2, c21, c22, c3, c4, d, d1, ddf, ddg, delta, df, dg, e, e1, e2, f, f1, f2, f3, g, g2, g3, ggamma, h, k, l, lambda_, m, mu, n, nu, omega, p, q, r, rho, s, sigma, t, v, x, x4, xp, z')
Af, Cf, F, F00, F01, F02, F10, F11, F12, F20, F21, F22, JacobiCN, JacobiDN, JacobiSN, LegendreP, LegendreQ, P, R1, R2, WeierstrassP, WeierstrassPPrime, _F1, arctan, b_0, c_0, d0, e0, f, f0, f_1, f_2, f_3, f4, f5, f_nu, g, g0, g1, g_nu, h, j, kf, p, phi, q, sint, tg, xf, x1, x2, x3, xf4, x5, y, zf = symbols('A, C, F, F00, F01, F02, F10, F11, F12, F20, F21, F22, JacobiCN, JacobiDN, JacobiSN, LegendreP, LegendreQ, P, R1, R2, WeierstrassP, WeierstrassPPrime, _F1, arctan, b0, c0, d0, e0, f, f0, f1, f2, f3, f4, f5, f_nu, g, g0, g1, g_nu, h, j, k, p, phi, q, sint, tg, x, x1, x2, x3, x4, x5, y, z', cls=Function)
AA = IndexedBase('AA')
AAA = 6*(a - 1)**2 - 2*c**2*(mu**2 + nu**2) + 1
BBB = 3*c - 2*a + 1
CCC = 2*c**2*(mu**2 + nu**2) - 2*a*(a - 1) - 1
DDD = (a - c)*(a - 2*c)
EEE = (mu*c + nu*c + a)*(mu*c + nu*c - a)*(mu*c - nu*c + a)*(mu*c - nu*c - a)

class Kamke:

    chapter_1 = {
        "kamke_1.1": Derivative(y(x), x) - 1/sqrt(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4),
        "kamke_1.2": a*y(x) - c*exp(b*x) + Derivative(y(x), x),
        "kamke_1.3": a*y(x) - b*sin(c*x) + Derivative(y(x), x),
        "kamke_1.4": 2*x*y(x) - x*exp(-x**2) + Derivative(y(x), x),
        "kamke_1.5": y(x)*cos(x) - exp(2*x) + Derivative(y(x), x),
        "kamke_1.6": y(x)*cos(x) - sin(2*x)/2 + Derivative(y(x), x),
        "kamke_1.7": y(x)*cos(x) + Derivative(y(x), x) - exp(-sin(x)),
        "kamke_1.8": y(x)*tan(x) - sin(2*x) + Derivative(y(x), x),
        "kamke_1.9": -(a + sin(log(x)) + cos(log(x)))*y(x) + Derivative(y(x), x),
        "kamke_1.10": -f(x)*Derivative(f(x), x) + y(x)*Derivative(f(x), x) + Derivative(y(x), x),
        "kamke_1.11": f(x)*y(x) - g(x) + Derivative(y(x), x),
        "kamke_1.12": y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.13": -a*x - b + y(x)**2 + Derivative(y(x), x),
        "kamke_1.15": x**4 - 2*x**2*y(x) - 2*x + y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.16": (x*y(x) - 1)*f(x) + y(x)**2 + Derivative(y(x), x),
        "kamke_1.17": -y(x)**2 - 3*y(x) + Derivative(y(x), x) + 4,
        "kamke_1.18": -x*y(x) - x - y(x)**2 + Derivative(y(x), x) + 1,
        "kamke_1.19": -(x + y(x))**2 + Derivative(y(x), x),
        "kamke_1.20": -2*x + (x**2 + 1)*y(x) - y(x)**2 + Derivative(y(x), x),
        "kamke_1.21": -y(x)**2 + y(x)*sin(x) - cos(x) + Derivative(y(x), x),
        "kamke_1.22": -y(x)**2 - y(x)*sin(2*x) - cos(2*x) + Derivative(y(x), x),
        "kamke_1.23": a*y(x)**2 - b + Derivative(y(x), x),
        "kamke_1.25": a*y(x)**2 - b*x**(2*nu) - c*x**(nu - 1) + Derivative(y(x), x),
        "kamke_1.26": -(A*y(x) - a)*(B*y(x) - b) + Derivative(y(x), x),
        "kamke_1.27": a*(-x + y(x))*y(x) + Derivative(y(x), x) - 1,
        "kamke_1.28": -x**3*y(x) + x*y(x)**2 - 2*x + Derivative(y(x), x),
        "kamke_1.29": -x*y(x)**2 - 3*x*y(x) + Derivative(y(x), x),
        "kamke_1.30": -x**a + x**(-a - 1)*y(x)**2 + Derivative(y(x), x),
        "kamke_1.32": y(x)**2*sin(x) - 2*sin(x)/cos(x)**2 + Derivative(y(x), x),
        "kamke_1.33": Derivative(y(x), x) - y(x)**2*Derivative(f(x), x)/g(x) + Derivative(g(x), x)/f(x),
        "kamke_1.34": f(x)*y(x)**2 + g(x)*y(x) + Derivative(y(x), x),
        "kamke_1.35": (2*a*y(x) + b + y(x)**2)*f(x) + Derivative(y(x), x),
        "kamke_1.36": a*x*y(x)**2 + y(x)**3 + Derivative(y(x), x),
        "kamke_1.37": -a*y(x)**2*exp(x) - y(x)**3 + Derivative(y(x), x),
        "kamke_1.38": -a*y(x)**3 - b/x**(3/2) + Derivative(y(x), x),
        "kamke_1.39": -a0 - a1*y(x) - a2*y(x)**2 - a3*y(x)**3 + Derivative(y(x), x),
        "kamke_1.40": 6*a*x*y(x)**2 + 3*a*y(x)**3 + Derivative(y(x), x),
        "kamke_1.41": a*x*y(x)**3 + b*y(x)**2 + Derivative(y(x), x),
        "kamke_1.42": -x*(x + 2)*y(x)**3 - (x + 3)*y(x)**2 + Derivative(y(x), x),
        "kamke_1.43": 3*x*y(x)**2 + (4*a**2*x + 3*a*x**2 + b)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.44": 2*a*x**3*y(x)**3 + 2*x*y(x) + Derivative(y(x), x),
        "kamke_1.45": 3*b*y(x)**2 + (2*a**2*x**3 - 2*b**2*x)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.46": a*x**(-a - 1) - x**a*y(x)**3 + 3*y(x)**2 + Derivative(y(x), x) - y(x)/x**a - 1/x**(2*a),
        "kamke_1.49": 6*a*phi(x)*y(x)**2 + a*y(x)**3*Derivative(phi(x), x) + 2*a + (2*a + 1)*y(x)*Derivative(phi(x), (x, 2))/Derivative(phi(x), x) + Derivative(y(x), x) + 2,
        "kamke_1.50": -f0(x) - f_1(x)*y(x) - f_2(x)*y(x)**2 - f_3(x)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.51": -(-f(x) + y(x))*(-g(x) + y(x))*(y(x) - (a*f(x) + b*g(x))/(a + b))*h(x) + Derivative(y(x), x) - (-g(x) + y(x))*Derivative(f(x), x)/(f(x) - g(x)) - (-f(x) + y(x))*Derivative(g(x), x)/(-f(x) + g(x)),
        "kamke_1.52": -a*y(x)**n - b*x**(n/(1 - n)) + Derivative(y(x), x),
        "kamke_1.53": -f(x)*Derivative(g(x), x) + Derivative(y(x), x) - y(x)*Derivative(f(x), x)/f(x) - f(x)**(1 - n)*y(x)**n*Derivative(g(x), x)/(a*g(x) + b)**n,
        "kamke_1.54": -a**n*f(x)**(1 - n)*y(x)**n*Derivative(g(x), x) - f(x)*Derivative(g(x), x) + Derivative(y(x), x) - y(x)*Derivative(f(x), x)/f(x),
        "kamke_1.55": -f(x)*y(x)**n - g(x)*y(x) - h(x) + Derivative(y(x), x),
        "kamke_1.56": -f(x)*y(x)**a - g(x)*y(x)**b + Derivative(y(x), x),
        "kamke_1.57": -sqrt(Abs(y(x))) + Derivative(y(x), x),
        "kamke_1.58": -a*sqrt(y(x)) - b*x + Derivative(y(x), x),
        "kamke_1.59": -a*sqrt(y(x)**2 + 1) - b + Derivative(y(x), x),
        "kamke_1.60": Derivative(y(x), x) - sqrt(y(x)**2 - 1)/sqrt(x**2 - 1),
        "kamke_1.61": -sqrt(x**2 - 1)/sqrt(y(x)**2 - 1) + Derivative(y(x), x),
        "kamke_1.62": -(-x**2*sqrt(x**2 - y(x)**2) + y(x))/(x*sqrt(x**2 - y(x)**2)*y(x) + x) + Derivative(y(x), x),
        "kamke_1.63": Derivative(y(x), x) - (y(x)**2 + 1)/((x + 1)**(3/2)*Abs(sqrt(y(x) + 1) + y(x))),
        "kamke_1.64": -sqrt((a*y(x)**2 + b*y(x) + c)/(a*x**2 + b*x + c)) + Derivative(y(x), x),
        "kamke_1.65": -sqrt((y(x)**3 + 1)/(x**3 + 1)) + Derivative(y(x), x),
        "kamke_1.66": Derivative(y(x), x) - sqrt(Abs((a*y(x) - 1)*(y(x) - 1)*y(x)))/sqrt(Abs(x*(x - 1)*(a*x - 1))),
        "kamke_1.67": Derivative(y(x), x) - sqrt(1 - y(x)**4)/sqrt(1 - x**4),
        "kamke_1.68": -sqrt((a*y(x)**4 + b*y(x)**2 + 1)/(a*x**4 + b*x**2 + 1)) + Derivative(y(x), x),
        "kamke_1.69": -sqrt((a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)*(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.70": -sqrt((a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)/(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.71": -sqrt((b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)/(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)) + Derivative(y(x), x),
        "kamke_1.72": -R1(x, sqrt(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4))*R2(y(x), sqrt(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.73": -((a0 + a1*x + a2*x**2 + a3*x**3)/(a0 + a1*y(x) + a2*y(x)**2 + a3*y(x)**3))**(2/3) + Derivative(y(x), x),
        "kamke_1.74": -sqrt((-a + y(x))*(-b + y(x)))*(-g(x) + y(x))*f(x) + Derivative(y(x), x),
        "kamke_1.75": exp(x) - exp(x - y(x)) + Derivative(y(x), x),
        "kamke_1.76": -a*cos(y(x)) + b + Derivative(y(x), x),
        "kamke_1.77": -cos(a*y(x) + b*x) + Derivative(y(x), x),
        "kamke_1.78": a*sin(alpha*y(x) + bbeta*x) + b + Derivative(y(x), x),
        "kamke_1.79": f(x)*cos(a*y(x)) + g(x)*sin(a*y(x)) + h(x) + Derivative(y(x), x),
        "kamke_1.80": (1 - Derivative(f(x), x))*cos(y(x)) + f(x)*sin(y(x)) - Derivative(f(x), x) + Derivative(y(x), x) - 1,
        "kamke_1.81": 2*tan(x)*tan(y(x)) + Derivative(y(x), x) - 1,
        "kamke_1.82": -a*(tan(y(x))**2 + 1) + tan(x)*tan(y(x)) + Derivative(y(x), x),
        "kamke_1.83": -tan(x*y(x)) + Derivative(y(x), x),
        "kamke_1.84": -f(a*x + b*y(x)) + Derivative(y(x), x),
        "kamke_1.85": -x**(a - 1)*f(y(x)**b/b + x**a/a)*y(x)**(1 - b) + Derivative(y(x), x),
        "kamke_1.86": -(-x*f(a*y(x)**2 + x**2) + y(x))/(a*f(a*y(x)**2 + x**2)*y(x) + x) + Derivative(y(x), x),
        "kamke_1.87": -(a*f(x**c*y(x))*y(x) + c*x**a*y(x)**b)/(b*x*f(x**c*y(x)) - x**a*y(x)**b) + Derivative(y(x), x),
        "kamke_1.88": -4*a*y(x) - b - c*exp(-2*a*x) - 3*y(x)**2 + 2*Derivative(y(x), x),
        "kamke_1.89": x*Derivative(y(x), x) - sqrt(a**2 - x**2),
        "kamke_1.90": -x*sin(x) + x*Derivative(y(x), x) + y(x),
        "kamke_1.91": x*Derivative(y(x), x) - x/log(x) - y(x),
        "kamke_1.92": -x**2*sin(x) + x*Derivative(y(x), x) - y(x),
        "kamke_1.93": x*Derivative(y(x), x) - x*cos(log(log(x)))/log(x) - y(x),
        "kamke_1.95": x**2 + x*Derivative(y(x), x) + y(x)**2,
        "kamke_1.96": x*Derivative(y(x), x) - y(x)**2 + 1,
        "kamke_1.97": a*y(x)**2 + b*x**2 + x*Derivative(y(x), x) - y(x),
        "kamke_1.98": a*y(x)**2 - b*y(x) + c*x**(2*b) + x*Derivative(y(x), x),
        "kamke_1.99": a*y(x)**2 - b*y(x) - c*x**bbeta + x*Derivative(y(x), x),
        "kamke_1.100": a + x*y(x)**2 + x*Derivative(y(x), x),
        "kamke_1.101": x*y(x)**2 + x*Derivative(y(x), x) - y(x),
        "kamke_1.102": -a*x**3 + x*y(x)**2 + x*Derivative(y(x), x) - y(x),
        "kamke_1.103": -x**3 + x*y(x)**2 + x*Derivative(y(x), x) - (2*x**2 + 1)*y(x),
        "kamke_1.104": a*x*y(x)**2 + b*x + x*Derivative(y(x), x) + 2*y(x),
        "kamke_1.105": a*x*y(x)**2 + b*y(x) + c*x + d + x*Derivative(y(x), x),
        "kamke_1.106": x*Derivative(y(x), x) + x**a*y(x)**2 + x**b + (a - b)*y(x)/2,
        "kamke_1.107": a*x**alpha*y(x)**2 + b*y(x) - c*x**bbeta + x*Derivative(y(x), x),
        "kamke_1.108": x*Derivative(y(x), x) - y(x)**2*log(x) + y(x),
        "kamke_1.109": x*Derivative(y(x), x) - (2*y(x)*log(x) - 1)*y(x),
        "kamke_1.110": x*Derivative(y(x), x) + (-x**2 + y(x)**2)*f(x),
        "kamke_1.111": 3*x*y(x)**2 + x*Derivative(y(x), x) + y(x)**3,
        "kamke_1.112": x*Derivative(y(x), x) - sqrt(x**2 + y(x)**2) - y(x),
        "kamke_1.114": -x*sqrt(x**2 + y(x)**2) + x*Derivative(y(x), x) - y(x),
        "kamke_1.115": -x*(-x + y(x))*sqrt(x**2 + y(x)**2) + x*Derivative(y(x), x) - y(x),
        "kamke_1.116": -x*sqrt((-4*x**2 + y(x)**2)*(-x**2 + y(x)**2)) + x*Derivative(y(x), x) - y(x),
        "kamke_1.117": -x*exp(y(x)/x) + x*Derivative(y(x), x) - x - y(x),
        "kamke_1.118": x*Derivative(y(x), x) - y(x)*log(y(x)),
        "kamke_1.119": x*Derivative(y(x), x) - (log(x*y(x)) - 1)*y(x),
        "kamke_1.120": x*Derivative(y(x), x) - (x*log(x**2/y(x)) + 2)*y(x),
        "kamke_1.121": x*Derivative(y(x), x) - sin(x - y(x)),
        "kamke_1.122": x*Derivative(y(x), x) + (-3*x**2*cos(y(x)) + sin(y(x)))*cos(y(x)),
        "kamke_1.123": -x*sin(y(x)/x) + x*Derivative(y(x), x) - y(x),
        "kamke_1.124": x*cos(y(x)/x) + x*Derivative(y(x), x) + x - y(x),
        "kamke_1.125": x*tan(y(x)/x) + x*Derivative(y(x), x) - y(x),
        "kamke_1.126": x*Derivative(y(x), x) - f(x*y(x))*y(x),
        "kamke_1.127": x*Derivative(y(x), x) - f(x**a*y(x)**b)*y(x),
        "kamke_1.128": a*y(x) + x*Derivative(y(x), x) - f(x)*g(x**a*y(x)),
        "kamke_1.129": (-x + y(x))*y(x) + (x + 1)*Derivative(y(x), x),
        "kamke_1.130": -2*x**3 + 2*x*Derivative(y(x), x) - y(x),
        "kamke_1.131": (2*x + 1)*Derivative(y(x), x) + 2 - 4*exp(-y(x)),
        "kamke_1.132": -3*x*y(x)**4*log(x) + 3*x*Derivative(y(x), x) - y(x),
        "kamke_1.133": x**2*Derivative(y(x), x) - x + y(x),
        "kamke_1.134": x**2*exp(x - 1/x) + x**2*Derivative(y(x), x) - y(x),
        "kamke_1.135": x**2*Derivative(y(x), x) - (x - 1)*y(x),
        "kamke_1.136": x**2*Derivative(y(x), x) + x**2 + x*y(x) + y(x)**2,
        "kamke_1.137": x**2*Derivative(y(x), x) - x*y(x) - y(x)**2,
        "kamke_1.138": x**2*Derivative(y(x), x) - x**2 - x*y(x) - y(x)**2,
        "kamke_1.139": a*x**k - b*(b - 1) + x**2*(y(x)**2 + Derivative(y(x), x)),
        "kamke_1.140": x**2*(y(x)**2 + Derivative(y(x), x)) + 4*x*y(x) + 2,
        "kamke_1.141": a*x*y(x) + b + x**2*(y(x)**2 + Derivative(y(x), x)),
        "kamke_1.142": -a*x**2*y(x) + a*x + x**2*(-y(x)**2 + Derivative(y(x), x)) + 2,
        "kamke_1.143": -b + x**2*(a*y(x)**2 + Derivative(y(x), x)),
        "kamke_1.144": b*x**alpha + c + x**2*(a*y(x)**2 + Derivative(y(x), x)),
        "kamke_1.145": -a*x**2*y(x)**2 + a*y(x)**3 + x**2*Derivative(y(x), x),
        "kamke_1.146": a*y(x)**2 + x**2*Derivative(y(x), x) + x*y(x)**3,
        "kamke_1.147": a*x**2*y(x)**3 + b*y(x)**2 + x**2*Derivative(y(x), x),
        "kamke_1.148": x*y(x) + (x**2 + 1)*Derivative(y(x), x) - 1,
        "kamke_1.149": -x*(x**2 + 1) + x*y(x) + (x**2 + 1)*Derivative(y(x), x),
        "kamke_1.150": -2*x**2 + 2*x*y(x) + (x**2 + 1)*Derivative(y(x), x),
        "kamke_1.151": (x**2 + 1)*Derivative(y(x), x) + (2*x*y(x) - 1)*(y(x)**2 + 1),
        "kamke_1.152": -x*(x**2 + 1)*cos(y(x))**2 + x*sin(y(x))*cos(y(x)) + (x**2 + 1)*Derivative(y(x), x),
        "kamke_1.153": a - x*y(x) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_1.154": 2*x*y(x) + (x**2 - 1)*Derivative(y(x), x) - cos(x),
        "kamke_1.155": -2*x*y(x) + (x**2 - 1)*Derivative(y(x), x) + y(x)**2 + 1,
        "kamke_1.156": -(-x + y(x))*y(x) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_1.157": a*(-2*x*y(x) + y(x)**2 + 1) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_1.158": a*x*y(x)**2 + x*y(x) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_1.159": -2*x*y(x)*log(y(x)) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_1.160": (x + 2)*y(x)**2 + (x**2 - 4)*Derivative(y(x), x) - 4*y(x),
        "kamke_1.161": x**2 + 3*x*y(x) + (x**2 - 5*x + 6)*Derivative(y(x), x) - 8*y(x),
        "kamke_1.162": k*(-a + x + y(x))*(-b + x + y(x)) + (-a + x)*(-b + x)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.163": 2*a**2*x + 2*x**2*Derivative(y(x), x) - x*y(x) - 2*y(x)**2,
        "kamke_1.164": 2*a**2*x + 2*x**2*Derivative(y(x), x) - 3*x*y(x) - 2*y(x)**2,
        "kamke_1.165": x*(2*x - 1)*Derivative(y(x), x) + 4*x - (4*x + 1)*y(x) + y(x)**2,
        "kamke_1.166": 2*x*(x - 1)*Derivative(y(x), x) - x + (x - 1)*y(x)**2,
        "kamke_1.167": 3*x**2*Derivative(y(x), x) - x**2 - 3*x*y(x) - 7*y(x)**2,
        "kamke_1.168": -x*y(x) + (3*x**2 - 12)*Derivative(y(x), x) + y(x)**2 - 3,
        "kamke_1.169": c*y(x)**2 + (a*x + b)**2*Derivative(y(x), x) + (a*x + b)*y(x)**3,
        "kamke_1.170": -x**4 + x**3*Derivative(y(x), x) - y(x)**2,
        "kamke_1.171": x**3*Derivative(y(x), x) - x**2*y(x) - y(x)**2,
        "kamke_1.172": -x**4*y(x)**2 + x**3*Derivative(y(x), x) + x**2*y(x) + 20,
        "kamke_1.173": -x**6*y(x)**2 + x**3*Derivative(y(x), x) - x**2*(2*x - 3)*y(x) + 3,
        "kamke_1.174": x**2*y(x) + x*(x**2 + 1)*Derivative(y(x), x),
        "kamke_1.175": a*x**3 + x*(x**2 - 1)*Derivative(y(x), x) - (2*x**2 - 1)*y(x),
        "kamke_1.176": -x**2 + x*(x**2 - 1)*Derivative(y(x), x) + (x**2 - 1)*y(x)**2,
        "kamke_1.177": x**2*(x - 1)*Derivative(y(x), x) - x*(x - 2)*y(x) - y(x)**2,
        "kamke_1.178": x**2 + 2*x*(x**2 - 1)*Derivative(y(x), x) + (2*x**2 - 2)*y(x)**2 - (3*x**2 - 5)*y(x) - 3,
        "kamke_1.179": 3*x*(x**2 - 1)*Derivative(y(x), x) + x*y(x)**2 - 3*x - (x**2 + 1)*y(x),
        "kamke_1.180": x**2 + (x*Derivative(y(x), x) - y(x))*(a*x**2 + b*x + c) - y(x)**2,
        "kamke_1.181": a + x**4*(y(x)**2 + Derivative(y(x), x)),
        "kamke_1.182": x**2 + x*(x**3 - 1)*Derivative(y(x), x) - 2*x*y(x)**2 + y(x),
        "kamke_1.183": -(2*x**3 - 2)*y(x) + (2*x**4 - x)*Derivative(y(x), x),
        "kamke_1.184": A + (y(x)**2 + Derivative(y(x), x))*(a*x**2 + b*x + c)**2,
        "kamke_1.185": x**7*Derivative(y(x), x) + 5*x**3*y(x)**2 + (2*x**2 + 2)*y(x)**3,
        "kamke_1.186": x**n*Derivative(y(x), x) - x**(n - 1)*(n - 1)*y(x) + x**(2*n - 2) + y(x)**2,
        "kamke_1.187": -a*y(x)**2 - b*x**(2*n - 2) + x**n*Derivative(y(x), x),
        "kamke_1.188": -a*y(x)**3 - b*x**(3*n) + x**(2*n + 1)*Derivative(y(x), x),
        "kamke_1.189": -a*y(x)**n - b*x**(n*(m + 1)) + x**(m*(n - 1) + n)*Derivative(y(x), x),
        "kamke_1.190": sqrt(x**2 - 1)*Derivative(y(x), x) - sqrt(y(x)**2 - 1),
        "kamke_1.191": sqrt(1 - x**2)*Derivative(y(x), x) - sqrt(y(x)**2 - 1)*y(x),
        "kamke_1.192": x + sqrt(a**2 + x**2)*Derivative(y(x), x) - sqrt(a**2 + x**2) + y(x),
        "kamke_1.193": -a*x*(log(x) + 1) + x*log(x)*Derivative(y(x), x) + y(x),
        "kamke_1.194": x*log(x)*Derivative(y(x), x) - (2*log(x)**2 + 1)*y(x) - y(x)**2*log(x) - log(x)**3,
        "kamke_1.195": (-3*sin(x) + cos(x))*y(x) - y(x)**2*sin(x)**2 + sin(x)*Derivative(y(x), x) + 4,
        "kamke_1.196": (sin(x) + 1)*cos(x) + y(x) + cos(x)*Derivative(y(x), x),
        "kamke_1.197": -y(x)**4 - y(x)*sin(x) + cos(x)*Derivative(y(x), x),
        "kamke_1.198": -y(x) - sin(x)**3 + sin(x)*cos(x)*Derivative(y(x), x),
        "kamke_1.199": sin(2*x)*Derivative(y(x), x) + sin(2*y(x)),
        "kamke_1.200": A*x*(a*sin(x)**2 + c) + a*y(x)*sin(2*x) + (a*sin(x)**2 + b)*Derivative(y(x), x),
        "kamke_1.201": -2*f(x)**2 + 2*f(x)*y(x)**2 + 2*f(x)*Derivative(y(x), x) - y(x)*Derivative(f(x), x),
        "kamke_1.202": f(x)*Derivative(y(x), x) + g(x)*tg(y(x)) + h(x),
        "kamke_1.203": x**3 + y(x)*Derivative(y(x), x) + y(x),
        "kamke_1.204": a*y(x) + x + y(x)*Derivative(y(x), x),
        "kamke_1.206": a*y(x) - 2*a + b*exp(x) + y(x)*Derivative(y(x), x),
        "kamke_1.207": 4*x*(x + 1) + y(x)**2 + y(x)*Derivative(y(x), x),
        "kamke_1.208": a*y(x)**2 - b*cos(c + x) + y(x)*Derivative(y(x), x),
        "kamke_1.209": -sqrt(a*y(x)**2 + b) + y(x)*Derivative(y(x), x),
        "kamke_1.210": x*y(x)**2 - 4*x + y(x)*Derivative(y(x), x),
        "kamke_1.211": -x*exp(x/y(x)) + y(x)*Derivative(y(x), x),
        "kamke_1.212": x + f(x**2 + y(x)**2)*g(x) + y(x)*Derivative(y(x), x),
        "kamke_1.213": -x + (y(x) + 1)*Derivative(y(x), x) - y(x),
        "kamke_1.214": 2*x + (x + y(x) - 1)*Derivative(y(x), x) - y(x) + 3,
        "kamke_1.215": x + (2*x + y(x) - 2)*Derivative(y(x), x) - y(x) + 1,
        "kamke_1.216": x + (-2*x + y(x) + 1)*Derivative(y(x), x) + y(x),
        "kamke_1.217": -x + (-x**2 + y(x))*Derivative(y(x), x),
        "kamke_1.218": 4*x*y(x) + (-x**2 + y(x))*Derivative(y(x), x),
        "kamke_1.219": (g(x) + y(x))*Derivative(y(x), x) - f0(x) - f_1(x)*y(x) - f_2(x)*y(x)**2,
        "kamke_1.220": -x**3 - x*y(x)**2 + 2*y(x)*Derivative(y(x), x),
        "kamke_1.221": -x + (x + 2*y(x) + 1)*Derivative(y(x), x) - 2*y(x) + 1,
        "kamke_1.222": 2*x + (x + 2*y(x) + 7)*Derivative(y(x), x) - y(x) + 4,
        "kamke_1.223": -2*x + (-x + 2*y(x))*Derivative(y(x), x) - y(x),
        "kamke_1.224": 3*x + (-6*x + 2*y(x))*Derivative(y(x), x) - y(x) + 2,
        "kamke_1.225": -x + (2*x + 4*y(x) + 3)*Derivative(y(x), x) - 2*y(x) - 1,
        "kamke_1.226": -x + (-2*x + 4*y(x) - 3)*Derivative(y(x), x) + 2*y(x) - 1,
        "kamke_1.227": 7*x + (-3*x + 4*y(x) - 5)*Derivative(y(x), x) - 3*y(x) + 2,
        "kamke_1.228": -8*x + (11*x + 4*y(x) - 11)*Derivative(y(x), x) - 25*y(x) + 62,
        "kamke_1.229": 2*x + (-5*x + 12*y(x) - 8)*Derivative(y(x), x) - 5*y(x) + 3,
        "kamke_1.230": a*y(x)*Derivative(y(x), x) + b*y(x)**2 + f(x),
        "kamke_1.231": alpha*y(x) + bbeta*x + ggamma + (a*y(x) + b*x + c)*Derivative(y(x), x),
        "kamke_1.232": x**2 + x*y(x)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.233": a*x**3*cos(x) + x*y(x)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.234": x**3 - 2*x**2 + x*y(x)*Derivative(y(x), x) + x*y(x) - y(x)**2,
        "kamke_1.235": b*y(x) + (a + x*y(x))*Derivative(y(x), x),
        "kamke_1.236": x*(y(x) + 4)*Derivative(y(x), x) - 2*x - y(x)**2 - 2*y(x),
        "kamke_1.237": b*y(x) + c*x + x*(a + y(x))*Derivative(y(x), x),
        "kamke_1.238": -b + (a + x*(x + y(x)))*Derivative(y(x), x) - (x + y(x))*y(x),
        "kamke_1.239": -2*x**2 - 3*x*y(x) + (-x**2 + x*y(x))*Derivative(y(x), x) + y(x)**2,
        "kamke_1.240": a*x + 2*x*y(x)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.241": a*x**2 + 2*x*y(x)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.242": 2*x*y(x)*Derivative(y(x), x) + 2*y(x)**2 + 1,
        "kamke_1.243": x*(x + 2*y(x) - 1)*Derivative(y(x), x) - (2*x + y(x) + 1)*y(x),
        "kamke_1.244": x*(-x + 2*y(x) - 1)*Derivative(y(x), x) + (2*x - y(x) - 1)*y(x),
        "kamke_1.245": 112*x**2*y(x) + (4*x**3 + 2*x*y(x))*Derivative(y(x), x) + y(x)**2,
        "kamke_1.246": x*(2*x + 3*y(x))*Derivative(y(x), x) + 3*(x + y(x))**2,
        "kamke_1.247": -7*x**2 + x*y(x) - 9*x + (3*x + 2)*(-2*x + y(x) - 1)*Derivative(y(x), x) - y(x)**2 - 3,
        "kamke_1.248": 2*x*y(x) + 2*x + (x**2 + 6*x*y(x) + 3)*Derivative(y(x), x) + 3*y(x)**2,
        "kamke_1.249": alpha*y(x)**3 + bbeta*y(x)**2 + (a*x*y(x) + b*x**n)*Derivative(y(x), x),
        "kamke_1.250": A*x*y(x) - B*g(x)**2 + alpha*x + bbeta*y(x) + ggamma + (A*x**2 + B*x*y(x) + a*x + b*y(x) + c)*Derivative(y(x), x),
        "kamke_1.251": x*y(x)**2 + (x**2*y(x) - 1)*Derivative(y(x), x) - 1,
        "kamke_1.252": -x*y(x)**2 + (x**2*y(x) - 1)*Derivative(y(x), x) + 1,
        "kamke_1.253": 8*x*y(x)**2 + (x**2*y(x) - 1)*Derivative(y(x), x) - 8,
        "kamke_1.254": x**2*y(x)**3 + x*(x*y(x) - 2)*Derivative(y(x), x) + x*y(x)**2 - 2*y(x),
        "kamke_1.255": x*(x*y(x) - 3)*Derivative(y(x), x) + x*y(x)**2 - y(x),
        "kamke_1.256": x**2*(y(x) - 1)*Derivative(y(x), x) + (x - 1)*y(x),
        "kamke_1.257": x*(x**4 + x*y(x) - 1)*Derivative(y(x), x) - (-x**4 + x*y(x) - 1)*y(x),
        "kamke_1.258": -2*x**3 + 2*x**2*y(x)*Derivative(y(x), x) - x**2 + y(x)**2,
        "kamke_1.259": 2*x**2*y(x)*Derivative(y(x), x) - x**2*exp(x - 1/x) - y(x)**2,
        "kamke_1.260": -x**2*y(x)**3 + 2*x*y(x)**2 + (2*x**2*y(x) + x)*Derivative(y(x), x) + y(x),
        "kamke_1.261": -2*x*y(x)**2 + (2*x**2*y(x) - x)*Derivative(y(x), x) - y(x),
        "kamke_1.262": 2*x**3 - 4*x*y(x)**2 + (-x**3 + 2*x**2*y(x))*Derivative(y(x), x) + y(x)**3,
        "kamke_1.263": 2*x**3 + 3*x**2*y(x)**2 + y(x)*Derivative(y(x), x) + 7,
        "kamke_1.264": 2*x*(x**3*y(x) + 1)*Derivative(y(x), x) + (3*x**3*y(x) - 1)*y(x),
        "kamke_1.265": 2*x**(n - 1)*(n + 1)**2*(x**(n**2)*y(x)**2 - 1) + (x**(n*(n + 1))*y(x) - 1)*Derivative(y(x), x),
        "kamke_1.266": -a*sqrt((y(x)**2 + 1)**3) + (-x + y(x))*sqrt(x**2 + 1)*Derivative(y(x), x),
        "kamke_1.267": y(x)**2*sin(x)*cos(x) + y(x)*sin(x)**2*Derivative(y(x), x) - 1,
        "kamke_1.268": f(x)*y(x)*Derivative(y(x), x) + g(x)*y(x)**2 + h(x),
        "kamke_1.269": (g0(x) + g1(x)*y(x))*Derivative(y(x), x) - f0(x) - f_1(x)*y(x) - f_2(x)*y(x)**2 - f_3(x)*y(x)**3,
        "kamke_1.270": x**2 + (-x + y(x)**2)*Derivative(y(x), x) - y(x),
        "kamke_1.271": 2*x*(2*x + y(x)) + (x**2 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.272": (x**2 + y(x)**2)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.273": 2*x*y(x) + (a + x**2 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.274": b + x**2 + 2*x*y(x) + (a + x**2 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.275": (x**2 + x + y(x)**2)*Derivative(y(x), x) - y(x),
        "kamke_1.276": 2*x*y(x) + (-x**2 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.277": -4*x**3*y(x) + (x**4 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.278": (y(x)**2 + 4*sin(x))*Derivative(y(x), x) - cos(x),
        "kamke_1.279": (x + y(x))**2*y(x)**2 + (y(x) + 1)*y(x) + (x + y(x)**2 + 2*y(x))*Derivative(y(x), x),
        "kamke_1.280": -a**2 + (x + y(x))**2*Derivative(y(x), x),
        "kamke_1.281": x**2 + 2*x*y(x) + (-x**2 + 2*x*y(x) + y(x)**2)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.282": -(2*y(x) - 1)*(6*x + 4*y(x) - 3) + (3*x + y(x) - 1)**2*Derivative(y(x), x),
        "kamke_1.283": -6*x*(x + 1)*y(x) + (-3*x**2 + 3*y(x)**2)*Derivative(y(x), x) + 2*y(x)**3 - 3*exp(x),
        "kamke_1.284": -x*y(x) + (x**2 + 4*y(x)**2)*Derivative(y(x), x),
        "kamke_1.285": 2*x**2 + 6*x*y(x) + (3*x**2 + 2*x*y(x) + 4*y(x)**2)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.286": (-3*x + 2*y(x) + 1)**2*Derivative(y(x), x) - (-2*x + 3*y(x) - 4)**2,
        "kamke_1.287": -(-2*x + y(x))**2 + (-4*x + 2*y(x) + 1)**2*Derivative(y(x), x),
        "kamke_1.288": -3*x*y(x)**2 + x + (-3*x**2*y(x) + 6*y(x)**2 + 1)*Derivative(y(x), x),
        "kamke_1.289": a + 2*x*y(x) + (-x + 6*y(x))**2*Derivative(y(x), x) - 6*y(x)**2,
        "kamke_1.290": b*y(x)**2 + 2*c*x*y(x) + d*x**2 + (a*y(x)**2 + 2*b*x*y(x) + c*x**2)*Derivative(y(x), x),
        "kamke_1.291": a*(alpha*x + bbeta*y(x))**2 - alpha*(a*x + b*y(x)) + (b*(alpha*x + bbeta*y(x))**2 - bbeta*(a*x + b*y(x)))*Derivative(y(x), x),
        "kamke_1.292": (a*y(x) + b*x + c)**2*Derivative(y(x), x) + (alpha*y(x) + bbeta*x + ggamma)**2,
        "kamke_1.293": x*(-3*x + y(x)**2)*Derivative(y(x), x) - 5*x*y(x) + 2*y(x)**3,
        "kamke_1.294": x*(-a + x**2 + y(x)**2)*Derivative(y(x), x) - (a + x**2 + y(x)**2)*y(x),
        "kamke_1.295": x**2*y(x) + x*(-x**2 + x*y(x) + y(x)**2)*Derivative(y(x), x) + x*y(x)**2 - y(x)**3,
        "kamke_1.296": x**4 - 2*x**2*y(x)**2 + x*(x**2*y(x) + x**2 + y(x)**2)*Derivative(y(x), x) - 2*y(x)**3,
        "kamke_1.297": -x**2*y(x) + 2*x*(5*x**2 + y(x)**2)*Derivative(y(x), x) + y(x)**3,
        "kamke_1.298": 3*x*y(x)**2*Derivative(y(x), x) - 2*x + y(x)**3,
        "kamke_1.299": -2*x*y(x) + (-x**2 + 3*x*y(x)**2)*Derivative(y(x), x) + y(x)**3,
        "kamke_1.300": 6*x*y(x)**2*Derivative(y(x), x) + x + 2*y(x)**3,
        "kamke_1.301": -(-x + 3*y(x)**2)*y(x) + (x**2 + 6*x*y(x)**2)*Derivative(y(x), x),
        "kamke_1.302": (x**2*y(x)**2 + x)*Derivative(y(x), x) + y(x),
        "kamke_1.303": x*(x*y(x) - 1)**2*Derivative(y(x), x) + (x**2*y(x)**2 + 1)*y(x),
        "kamke_1.304": 5*x**2*y(x)**3 + x*y(x)**2 + (10*x**3*y(x)**2 + x**2*y(x) + 2*x)*Derivative(y(x), x),
        "kamke_1.305": x**2 + (-3*x + y(x)**3)*Derivative(y(x), x) - 3*y(x),
        "kamke_1.306": -x**2*y(x) + (-x**3 + y(x)**3)*Derivative(y(x), x),
        "kamke_1.307": x*(-a + x**2 + y(x)**2) + (a + x**2 + y(x)**2)*y(x)*Derivative(y(x), x),
        "kamke_1.308": x*y(x)**2 + 2*y(x)**3*Derivative(y(x), x),
        "kamke_1.309": -2*x**3 - x + (2*y(x)**3 + y(x))*Derivative(y(x), x),
        "kamke_1.310": x**3 + 5*x*y(x)**2 + (5*x**2*y(x) + 2*y(x)**3)*Derivative(y(x), x),
        "kamke_1.311": 4*x**3 + 9*x**2*y(x) + 6*x*y(x)**2 + (3*x**3 + 6*x**2*y(x) - 3*x*y(x)**2 + 20*y(x)**3)*Derivative(y(x), x) - y(x)**3,
        "kamke_1.312": (a - b)*(-x + y(x)*Derivative(y(x), x))/(a + b) + (x + y(x)*Derivative(y(x), x))*(y(x)**2/b + x**2/a),
        "kamke_1.313": -a*y(x)**3 + 2*b*x**3 + 3*b*x**2*y(x) + c*y(x)**2 + (3*a*x*y(x)**2 + 2*a*y(x)**3 - b*x**3 + c*x**2)*Derivative(y(x), x),
        "kamke_1.314": x*y(x)**3*Derivative(y(x), x) - x*sin(x) + y(x)**4,
        "kamke_1.315": 2*x**3*y(x) + (-x**4 + 2*x*y(x)**3)*Derivative(y(x), x) - y(x)**4,
        "kamke_1.316": (2*x*y(x)**3 + y(x))*Derivative(y(x), x) + 2*y(x)**2,
        "kamke_1.317": -x*y(x) + (x**2 + 2*x*y(x)**3 + x*y(x))*Derivative(y(x), x) + y(x)**2,
        "kamke_1.318": (y(x)**2 - 2)*y(x)**2 + (3*x*y(x)**3 - 4*x*y(x) + y(x))*Derivative(y(x), x),
        "kamke_1.319": (7*x*y(x)**3 - 5*x + y(x))*Derivative(y(x), x) + y(x)**4 - 5*y(x),
        "kamke_1.320": (x**2*y(x)**3 + x*y(x))*Derivative(y(x), x) - 1,
        "kamke_1.321": (2*x**2*y(x)**3 + x**2*y(x)**2 - 2*x)*Derivative(y(x), x) - 2*y(x) - 1,
        "kamke_1.322": 5*x*y(x)**4 + x + (10*x**2*y(x)**3 - 3*y(x)**2 - 2)*Derivative(y(x), x),
        "kamke_1.323": x*(a*x*y(x)**3 + c)*Derivative(y(x), x) + (b*x**3*y(x) + c)*y(x),
        "kamke_1.324": 2*x**3*y(x)**3 + (2*x**3*y(x)**3 - x)*Derivative(y(x), x) - y(x),
        "kamke_1.325": x*(-x**3 + 2*y(x)**3) + (-2*x**3 + y(x)**3)*y(x)*Derivative(y(x), x),
        "kamke_1.326": x*(a*y(x)**3 + (a*y(x) + b*x)**3) + (b*x**3 + (a*y(x) + b*x)**3)*y(x)*Derivative(y(x), x),
        "kamke_1.327": (2*x**2*y(x)**3 + x*y(x)**4 + x + 2*y(x))*Derivative(y(x), x) + y(x)**5 + y(x),
        "kamke_1.328": a*x**2*y(x)**n*Derivative(y(x), x) - 2*x*Derivative(y(x), x) + y(x),
        "kamke_1.329": alpha*x*Derivative(y(x), x) + bbeta*y(x) + x**n*(a*x*Derivative(y(x), x) + b*y(x))*y(x)**m,
        "kamke_1.330": (f(x + y(x)) + 1)*Derivative(y(x), x) + f(x + y(x)),
        "kamke_1.331": (-y(x) + y(x)**(ps + 1))*f_nu(x)*Derivative(y(x), x)/(y(x) - 1) - (-y(x) + y(x)**(qs + 1))*g_nu(x)/(y(x) - 1),
        "kamke_1.332": x*(sqrt(x*y(x)) - 1)*Derivative(y(x), x) - (sqrt(x*y(x)) + 1)*y(x),
        "kamke_1.333": -x**(3/2)*y(x)**(5/2) + x*y(x)**2 + (2*x**(5/2)*y(x)**(3/2) + x**2*y(x) - x)*Derivative(y(x), x) - y(x),
        "kamke_1.334": (sqrt(x + y(x)) + 1)*Derivative(y(x), x) + 1,
        "kamke_1.335": -sqrt(x**2 - 1) + sqrt(y(x)**2 - 1)*Derivative(y(x), x),
        "kamke_1.336": a*y(x) + sqrt(x**2 + 1) + (a*x + sqrt(y(x)**2 + 1))*Derivative(y(x), x),
        "kamke_1.337": (x + sqrt(x**2 + y(x)**2))*Derivative(y(x), x) - y(x),
        "kamke_1.338": x*sqrt(x**2 + y(x)**2) + 2*x*y(x)*sin(alpha) + (-x**2 + y(x)**2)*cos(alpha) + (-2*x*y(x)*cos(alpha) + (-x**2 + y(x)**2)*sin(alpha) + sqrt(x**2 + y(x)**2)*y(x))*Derivative(y(x), x),
        "kamke_1.339": -x*(x**2 + y(x)**2) + (x*sqrt(x**2 + y(x)**2 + 1) - (x**2 + y(x)**2)*y(x))*Derivative(y(x), x) - sqrt(x**2 + y(x)**2 + 1)*y(x),
        "kamke_1.340": -(e1/((a + x)**2 + y(x)**2)**(3/2) + e2/((-a + x)**2 + y(x)**2)**(3/2))*y(x) + (e1*(a + x)/((a + x)**2 + y(x)**2)**(3/2) + e2*(-a + x)/((-a + x)**2 + y(x)**2)**(3/2))*Derivative(y(x), x),
        "kamke_1.341": (x*exp(y(x)) + exp(x))*Derivative(y(x), x) + y(x)*exp(x) + exp(y(x)),
        "kamke_1.342": x*(x*Derivative(y(x), x) + y(x))*(3*exp(x*y(x)) + 2*exp(-x*y(x))) + 1,
        "kamke_1.343": (x + log(y(x)))*Derivative(y(x), x) - 1,
        "kamke_1.344": (2*x + log(y(x)) - 1)*Derivative(y(x), x) - 2*y(x),
        "kamke_1.345": x*(2*x**2*y(x)*log(y(x)) + 1)*Derivative(y(x), x) - 2*y(x),
        "kamke_1.346": x*(-a*x + y(x)*log(x*y(x)) + y(x))*Derivative(y(x), x) - (a*x*log(x*y(x)) + a*x - y(x))*y(x),
        "kamke_1.347": (sin(x) + 1)*sin(y(x))*Derivative(y(x), x) + (cos(y(x)) - 1)*cos(x),
        "kamke_1.348": (x*cos(y(x)) + sin(x))*Derivative(y(x), x) + y(x)*cos(x) + sin(y(x)),
        "kamke_1.349": 2*x*sin(y(x)/x) + x*cot(y(x)/x)*Derivative(y(x), x) - y(x)*cot(y(x)/x),
        "kamke_1.350": -sin(y(x))**2*cos(x) - sin(y(x)) + cos(y(x))*Derivative(y(x), x),
        "kamke_1.351": x*sin(y(x))*cos(y(x))**2 - sin(y(x))**3 + cos(y(x))*Derivative(y(x), x),
        "kamke_1.352": (-sin(alpha)*sin(x) + cos(y(x)))*cos(y(x))*Derivative(y(x), x) + (-sin(alpha)*sin(y(x)) + cos(x))*cos(x),
        "kamke_1.353": x*cos(y(x))*Derivative(y(x), x) + sin(y(x)),
        "kamke_1.354": (x*sin(y(x)) - 1)*Derivative(y(x), x) + cos(y(x)),
        "kamke_1.355": (x*cos(y(x)) + cos(x))*Derivative(y(x), x) - y(x)*sin(x) + sin(y(x)),
        "kamke_1.356": 2*x*sin(y(x)) + (x**2*cos(y(x)) + 2*y(x)*sin(x))*Derivative(y(x), x) + y(x)**2*cos(x),
        "kamke_1.357": x*log(x)*sin(y(x))*Derivative(y(x), x) + (-x*cos(y(x)) + 1)*cos(y(x)),
        "kamke_1.358": sin(x)*cos(y(x)) + sin(y(x))*cos(x)*Derivative(y(x), x),
        "kamke_1.359": 5*y(x)*cos(x)**4 + 3*sin(x)*sin(y(x))*Derivative(y(x), x),
        "kamke_1.360": -b*(-c*cos(a*y(x)) + 1)*sqrt(c*cos(a*y(x)) + cos(a*y(x))**2 - 1) + cos(a*y(x))*Derivative(y(x), x),
        "kamke_1.361": (x*sin(x*y(x)) - sin(y(x)) + cos(x + y(x)))*Derivative(y(x), x) + y(x)*sin(x*y(x)) + cos(x) + cos(x + y(x)),
        "kamke_1.362": x*y(x)**2*sin(x*y(x)) + (x**2*y(x)*sin(x*y(x)) - 4*x)*Derivative(y(x), x) - y(x),
        "kamke_1.363": x + (x*Derivative(y(x), x) - y(x))*cos(y(x)/x)**2,
        "kamke_1.364": x*(-x*cos(y(x)/x) + y(x)*sin(y(x)/x))*Derivative(y(x), x) - (x*cos(y(x)/x) + y(x)*sin(y(x)/x))*y(x),
        "kamke_1.365": x*f(x**2 + y(x)**2) + (-x + f(x**2 + y(x)**2)*y(x))*Derivative(y(x), x) + y(x),
        "kamke_1.366": -x*Derivative(y(x), x) + (a*y(x)*Derivative(y(x), x) + x)*f(a*y(x)**2 + x**2) - y(x),
        "kamke_1.367": -x**a*(c*y(x) + x*Derivative(y(x), x))*y(x)**b + (-a + b*x*Derivative(y(x), x))*f(x**c*y(x)),
        "kamke_1.368": a*y(x) + b*x**2 + Derivative(y(x), x)**2,
        "kamke_1.369": -a**2 + y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.370": -f(x)**2 + y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.371": -y(x)**3 + y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.372": a*y(x) + b - 4*y(x)**3 + Derivative(y(x), x)**2,
        "kamke_1.373": a**2*(log(y(x))**2 - 1)*y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.374": -y(x)**2 + Derivative(y(x), x)**2 - 2*Derivative(y(x), x),
        "kamke_1.375": a*Derivative(y(x), x) + b*x + Derivative(y(x), x)**2,
        "kamke_1.376": a*Derivative(y(x), x) + b*y(x) + Derivative(y(x), x)**2,
        "kamke_1.377": (x - 2)*Derivative(y(x), x) - y(x) + Derivative(y(x), x)**2 + 1,
        "kamke_1.378": (a + x)*Derivative(y(x), x) - y(x) + Derivative(y(x), x)**2,
        "kamke_1.379": -(x + 1)*Derivative(y(x), x) + y(x) + Derivative(y(x), x)**2,
        "kamke_1.380": 2*x*Derivative(y(x), x) - y(x) + Derivative(y(x), x)**2,
        "kamke_1.381": -2*x*Derivative(y(x), x) + y(x) + Derivative(y(x), x)**2,
        "kamke_1.382": a*x*Derivative(y(x), x) - b*x**2 - c + Derivative(y(x), x)**2,
        "kamke_1.383": a*x*Derivative(y(x), x) + b*y(x) + c*x**2 + Derivative(y(x), x)**2,
        "kamke_1.384": -a*y(x) + c + (a*x + b)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.385": -2*x**2*Derivative(y(x), x) + 2*x*y(x) + Derivative(y(x), x)**2,
        "kamke_1.386": a*x**3*Derivative(y(x), x) - 2*a*x**2*y(x) + Derivative(y(x), x)**2,
        "kamke_1.387": (-y(x) + Derivative(y(x), x))*exp(x) + Derivative(y(x), x)**2,
        "kamke_1.388": -2*x - 2*y(x)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.389": (4*y(x) + 1)*y(x) - (4*y(x) + 1)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.390": a*y(x)*Derivative(y(x), x) - b*x - c + Derivative(y(x), x)**2,
        "kamke_1.391": a*b*x*y(x) + (a*y(x) + b*x)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.392": -x*y(x)*Derivative(y(x), x) + y(x)**2*log(a*y(x)) + Derivative(y(x), x)**2,
        "kamke_1.393": -y(x)**2 + 2*y(x)*cot(x)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.394": -(-f(x)**2 + g(x))*exp(-2*Integral(f(xp), (xp, a, x))) + 2*f(x)*y(x)*Derivative(y(x), x) + g(x)*y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.395": 2*f(x)*y(x)*Derivative(y(x), x) + g(x)*y(x)**2 + h(x) + Derivative(y(x), x)**2,
        "kamke_1.396": -x*y(x)**3 + (-x + y(x))*y(x)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.397": -2*x**3*y(x)**2*Derivative(y(x), x) - 4*x**2*y(x)**3 + Derivative(y(x), x)**2,
        "kamke_1.398": -3*x*y(x)**(2/3)*Derivative(y(x), x) + 9*y(x)**(5/3) + Derivative(y(x), x)**2,
        "kamke_1.399": (x - 1)*Derivative(y(x), x) - y(x) + 2*Derivative(y(x), x)**2,
        "kamke_1.400": -2*x**2*Derivative(y(x), x) + 3*x*y(x) + 2*Derivative(y(x), x)**2,
        "kamke_1.401": -2*x*Derivative(y(x), x) + y(x) + 3*Derivative(y(x), x)**2,
        "kamke_1.402": x**2 + 4*x*Derivative(y(x), x) - y(x) + 3*Derivative(y(x), x)**2,
        "kamke_1.403": a*Derivative(y(x), x)**2 + b*Derivative(y(x), x) - y(x),
        "kamke_1.404": a*Derivative(y(x), x)**2 + b*x**2*Derivative(y(x), x) + c*x*y(x),
        "kamke_1.405": a*Derivative(y(x), x)**2 - x + y(x)*Derivative(y(x), x),
        "kamke_1.406": a*Derivative(y(x), x)**2 - x - y(x)*Derivative(y(x), x),
        "kamke_1.407": x*Derivative(y(x), x)**2 - y(x),
        "kamke_1.408": x*Derivative(y(x), x)**2 + x - 2*y(x),
        "kamke_1.409": x*Derivative(y(x), x)**2 - y(x) - 2*Derivative(y(x), x),
        "kamke_1.410": x*Derivative(y(x), x)**2 - 2*y(x) + 4*Derivative(y(x), x),
        "kamke_1.411": x*Derivative(y(x), x)**2 + x*Derivative(y(x), x) - y(x),
        "kamke_1.412": a + x*Derivative(y(x), x)**2 + y(x)*Derivative(y(x), x),
        "kamke_1.413": -x**2 + x*Derivative(y(x), x)**2 + y(x)*Derivative(y(x), x),
        "kamke_1.414": x**3 + x*Derivative(y(x), x)**2 + y(x)*Derivative(y(x), x),
        "kamke_1.415": x*Derivative(y(x), x)**2 - y(x)**4 + y(x)*Derivative(y(x), x),
        "kamke_1.416": x*Derivative(y(x), x)**2 + (-3*x + y(x))*Derivative(y(x), x) + y(x),
        "kamke_1.417": a + x*Derivative(y(x), x)**2 - y(x)*Derivative(y(x), x),
        "kamke_1.418": a*y(x) + x*Derivative(y(x), x)**2 - y(x)*Derivative(y(x), x),
        "kamke_1.419": x*Derivative(y(x), x)**2 - x + 2*y(x)*Derivative(y(x), x),
        "kamke_1.420": a + x*Derivative(y(x), x)**2 - 2*y(x)*Derivative(y(x), x),
        "kamke_1.421": x*Derivative(y(x), x)**2 - x - 2*y(x)*Derivative(y(x), x),
        "kamke_1.422": x*Derivative(y(x), x)**2 + 4*x - 2*y(x)*Derivative(y(x), x),
        "kamke_1.423": x*Derivative(y(x), x)**2 + x - 2*y(x)*Derivative(y(x), x) + 2*y(x),
        "kamke_1.424": a*y(x)*Derivative(y(x), x) + b*x + x*Derivative(y(x), x)**2,
        "kamke_1.425": (x + 1)*Derivative(y(x), x)**2 - (x + y(x))*Derivative(y(x), x) + y(x),
        "kamke_1.426": (3*x + 1)*Derivative(y(x), x)**2 - (3*y(x) + 6)*Derivative(y(x), x) + 9,
        "kamke_1.427": -(x + 3*y(x))*Derivative(y(x), x) + (3*x + 5)*Derivative(y(x), x)**2 + y(x),
        "kamke_1.428": a*x*Derivative(y(x), x)**2 - b*y(x) + (-a*y(x) + b*x + c)*Derivative(y(x), x),
        "kamke_1.429": a*x*Derivative(y(x), x)**2 + b*y(x) - (a*y(x) - a + b*x - b)*Derivative(y(x), x),
        "kamke_1.430": a0*x + b0*y(x) + c0 + (a2*x + c2)*Derivative(y(x), x)**2 + (a1*x + b1*y(x) + c1)*Derivative(y(x), x),
        "kamke_1.431": x**2*Derivative(y(x), x)**2 - y(x)**4 + y(x)**2,
        "kamke_1.432": -2*a*y(x) + x**2 + (a + x*Derivative(y(x), x))**2,
        "kamke_1.433": -4*a - 4*x**2 - 4*x*y(x) + (x*Derivative(y(x), x) + 2*x + y(x))**2,
        "kamke_1.434": Derivative(y(x), x) - 1,
        "kamke_1.435": x**2*Derivative(y(x), x)**2 - 2*x*y(x)*Derivative(y(x), x) - x + (y(x) + 1)*y(x),
        "kamke_1.436": -x**4 + x**2*Derivative(y(x), x)**2 - 2*x*y(x)*Derivative(y(x), x) + (1 - x**2)*y(x)**2,
        "kamke_1.437": x**2*Derivative(y(x), x)**2 - (a + 2*x*y(x))*Derivative(y(x), x) + y(x)**2,
        "kamke_1.438": x**2*Derivative(y(x), x)**2 + 3*x*y(x)*Derivative(y(x), x) + 2*y(x)**2,
        "kamke_1.439": x**2*Derivative(y(x), x)**2 + 3*x*y(x)*Derivative(y(x), x) + 3*y(x)**2,
        "kamke_1.440": x**2*Derivative(y(x), x)**2 + 4*x*y(x)*Derivative(y(x), x) - 5*y(x)**2,
        "kamke_1.441": x**2*Derivative(y(x), x)**2 - 4*x*(y(x) + 2)*Derivative(y(x), x) + 4*(y(x) + 2)*y(x),
        "kamke_1.442": x**2*Derivative(y(x), x)**2 + (1 - x)*(-x**2*y(x) + y(x)**2) + (x**3 + x**2*y(x) - 2*x*y(x))*Derivative(y(x), x),
        "kamke_1.443": x*(x*Derivative(y(x), x) - y(x))**2 - Derivative(y(x), x),
        "kamke_1.444": x**2*Derivative(y(x), x)**2 - (-2*x + y(x))*y(x)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.445": a*b*y(x)**3 + x**2*Derivative(y(x), x)**2 + (a*x**2*y(x)**3 + b)*Derivative(y(x), x),
        "kamke_1.446": -2*x*y(x)*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), x)**2 + y(x)**2 - 1,
        "kamke_1.447": (x**2 - 1)*Derivative(y(x), x)**2 - 1,
        "kamke_1.448": (x**2 - 1)*Derivative(y(x), x)**2 - y(x)**2 + 1,
        "kamke_1.449": 2*x*y(x)*Derivative(y(x), x) + (-a**2 + x**2)*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.450": -x**2 - 2*x*y(x)*Derivative(y(x), x) + (-a**2 + x**2)*Derivative(y(x), x)**2,
        "kamke_1.451": b - 2*x*y(x)*Derivative(y(x), x) + (a + x**2)*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.452": (2*x**2 + 1)*Derivative(y(x), x)**2 + (x**2 + 2*x*y(x) + y(x)**2 + 2)*Derivative(y(x), x) + 2*y(x)**2 + 1,
        "kamke_1.453": a**2*x**2 + x**2*(a**2 - 1)*Derivative(y(x), x)**2 + 2*x*y(x)*Derivative(y(x), x) - y(x)**2,
        "kamke_1.454": -a*x**2*(a - 1) + a*x**2*Derivative(y(x), x)**2 - 2*a*x*y(x)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.455": a + x**3*Derivative(y(x), x)**2 + x**2*y(x)*Derivative(y(x), x),
        "kamke_1.456": x*(x**2 - 1)*Derivative(y(x), x)**2 + x*y(x)**2 - x + (2 - 2*x**2)*y(x)*Derivative(y(x), x),
        "kamke_1.457": x**4*Derivative(y(x), x)**2 - x*Derivative(y(x), x) - y(x),
        "kamke_1.458": x**2*(-a**2 + x**2)*Derivative(y(x), x)**2 - 1,
        "kamke_1.459": -(Derivative(y(x), x) - 1)**2 + exp(-2*y(x)) + exp(-2*x)*Derivative(y(x), x)**2,
        "kamke_1.460": -a**2 + (y(x)**2 + Derivative(y(x), x)**2)*cos(x)**4,
        "kamke_1.461": 2*b_0(x)*y(x)*Derivative(y(x), x) + c_0(x)*y(x)**2 + d0(x)*Derivative(y(x), x)**2 + 2*d0(x)*Derivative(y(x), x) + 2*e0(x)*y(x) + f0(x),
        "kamke_1.462": y(x)*Derivative(y(x), x)**2 - 1,
        "kamke_1.463": y(x)*Derivative(y(x), x)**2 - exp(2*x),
        "kamke_1.464": 2*x*Derivative(y(x), x) + y(x)*Derivative(y(x), x)**2 - y(x),
        "kamke_1.465": 2*x*Derivative(y(x), x) + y(x)*Derivative(y(x), x)**2 - 9*y(x),
        "kamke_1.466": -2*x*Derivative(y(x), x) + y(x)*Derivative(y(x), x)**2 + y(x),
        "kamke_1.467": -4*x*Derivative(y(x), x) + y(x)*Derivative(y(x), x)**2 + y(x),
        "kamke_1.468": -4*a**2*x*Derivative(y(x), x) + a**2*y(x) + y(x)*Derivative(y(x), x)**2,
        "kamke_1.469": a*x*Derivative(y(x), x) + b*y(x) + y(x)*Derivative(y(x), x)**2,
        "kamke_1.470": x**3*Derivative(y(x), x) - x**2*y(x) + y(x)*Derivative(y(x), x)**2,
        "kamke_1.471": -x - (-x + y(x))*Derivative(y(x), x) + y(x)*Derivative(y(x), x)**2,
        "kamke_1.472": 2*x*Derivative(y(x), x) + (x + y(x))*Derivative(y(x), x)**2 - y(x),
        "kamke_1.473": (-2*x + y(x))*Derivative(y(x), x)**2 - (2*x - 2)*Derivative(y(x), x) + y(x) - 2,
        "kamke_1.474": -(4*x - 5)*Derivative(y(x), x) + 2*y(x)*Derivative(y(x), x)**2 + 2*y(x),
        "kamke_1.475": 2*x*Derivative(y(x), x) + 4*y(x)*Derivative(y(x), x)**2 - y(x),
        "kamke_1.476": 4*x**3*Derivative(y(x), x) - 4*x**2*y(x) + 9*y(x)*Derivative(y(x), x)**2,
        "kamke_1.477": a*y(x)*Derivative(y(x), x)**2 + (-b + 2*x)*Derivative(y(x), x) - y(x),
        "kamke_1.478": -c + (a*y(x) + b)*(Derivative(y(x), x)**2 + 1),
        "kamke_1.479": a0*x + b0*y(x) + c0 + (a1*x + b1*y(x) + c1)*Derivative(y(x), x) + (a2*x + b2*y(x) + c2)*Derivative(y(x), x)**2,
        "kamke_1.480": 2*x*y(x)*Derivative(y(x), x)**2 + (a*y(x) - x**2)*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.481": x*y(x)*Derivative(y(x), x)**2 + x*y(x) + (x**2 + y(x)**2)*Derivative(y(x), x),
        "kamke_1.482": x*y(x)*Derivative(y(x), x)**2 - x*y(x) + (a + x**22 - y(x)**2)*Derivative(y(x), x),
        "kamke_1.483": 2*x*y(x)*Derivative(y(x), x) + 2*x*y(x) + (-x**2 + 2*x*y(x))*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.484": -6*x*y(x)*Derivative(y(x), x) + 2*x*y(x) + (-x**2 + 2*x*y(x))*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.485": a*x*y(x)*Derivative(y(x), x)**2 + b*x*y(x) - (a*y(x)**2 + b*x**2 + c)*Derivative(y(x), x),
        "kamke_1.486": -a**2 + y(x)**2*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.487": -6*x**3*Derivative(y(x), x) + 4*x**2*y(x) + y(x)**2*Derivative(y(x), x)**2,
        "kamke_1.488": 4*a**2 - 4*a*x - 4*a*y(x)*Derivative(y(x), x) + y(x)**2*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.489": a*y(x)**2 + b*x + c + 2*x*y(x)*Derivative(y(x), x) + y(x)**2*Derivative(y(x), x)**2,
        "kamke_1.490": a - x**2 - 2*x*y(x)*Derivative(y(x), x) + y(x)**2*Derivative(y(x), x)**2 + 2*y(x)**2,
        "kamke_1.491": a*x**2 + 2*a*x*y(x)*Derivative(y(x), x) + b*(a - 1) + (1 - a)*y(x)**2 + y(x)**2*Derivative(y(x), x)**2,
        "kamke_1.492": (-a**2 + y(x)**2)*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.493": 2*a*y(x)*Derivative(y(x), x) + (a**2 - 2*a*x + y(x)**2)*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.494": x**2*(1 - a**2) + 2*x*y(x)*Derivative(y(x), x) + (-a**2*x**2 + y(x)**2)*Derivative(y(x), x)**2,
        "kamke_1.495": 2*a*x*y(x)*Derivative(y(x), x) + x**2 + (1 - a)*y(x)**2 + (x**2*(1 - a) + y(x)**2)*Derivative(y(x), x)**2,
        "kamke_1.496": -a**2*(Derivative(y(x), x) + 1)**2 + (-x + y(x))**2*(Derivative(y(x), x)**2 + 1),
        "kamke_1.497": -x**2 - 2*x*y(x)*Derivative(y(x), x) + 3*y(x)**2*Derivative(y(x), x)**2 + 4*y(x)**2,
        "kamke_1.498": (3*y(x) - 2)*Derivative(y(x), x)**2 + 4*y(x) - 4,
        "kamke_1.499": -a**2*x**2 - 2*a**2*x*y(x)*Derivative(y(x), x) + (1 - a**2)*y(x)**2*Derivative(y(x), x)**2 + y(x)**2,
        "kamke_1.500": -a*b + a*y(x)**2 - b*x**2 - 2*b*x*y(x)*Derivative(y(x), x) + (a - b)*y(x)**2*Derivative(y(x), x)**2,
        "kamke_1.501": -b*y(x)*Derivative(y(x), x) + d*y(x)**2 + (a*y(x)**2 + b*x + c)*Derivative(y(x), x)**2,
        "kamke_1.502": -c**2*(a*Derivative(y(x), x) + b)**2 + (a*y(x) - b*x)**2*(a**2*Derivative(y(x), x)**2 + b**2),
        "kamke_1.503": a0 + b0*y(x) + c0 + (a1*x + b1*y(x) + c1)*Derivative(y(x), x) + (a2*x + b2*y(x) + c2)**2*Derivative(y(x), x)**2,
        "kamke_1.504": x**2*y(x) + x*y(x)**2*Derivative(y(x), x)**2 - (-a + x**3 + y(x)**3)*Derivative(y(x), x),
        "kamke_1.505": -x**3 + x*y(x)**2*Derivative(y(x), x)**2 + 2*x*y(x)**2 - 2*y(x)**3*Derivative(y(x), x),
        "kamke_1.506": 2*x**2*(-x + y(x))*y(x)**2*Derivative(y(x), x) + x**2*(x*y(x)**2 - 1)*Derivative(y(x), x)**2 - (x**2*y(x) - 1)*y(x)**2,
        "kamke_1.507": 2*a**2*x*y(x)*Derivative(y(x), x) + (-a**2 + y(x)**2)*y(x)**2 + (-a**2*x**2 + y(x)**4)*Derivative(y(x), x)**2,
        "kamke_1.508": 2*x*y(x)*Derivative(y(x), x) + (x**2*y(x)**2 - x**2 + y(x)**4)*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.509": -4*x**2 - 6*x*y(x)**5*Derivative(y(x), x) + 9*(x**2 - 1)*y(x)**4*Derivative(y(x), x)**2,
        "kamke_1.510": 2*x**3*(-x**2 + y(x)**2)*y(x)**3*Derivative(y(x), x) + x**2*(x**2*y(x)**4 - 1)*Derivative(y(x), x)**2 - (x**4*y(x)**2 - 1)*y(x)**2,
        "kamke_1.511": a**2*sqrt(x**2 + y(x)**2) + 2*x*y(x)*Derivative(y(x), x) + (a**2*sqrt(x**2 + y(x)**2) - x**2)*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.512": a*(x**2 + y(x)**2)**(3/2) + 2*x*y(x)*Derivative(y(x), x) + (a*(x**2 + y(x)**2)**(3/2) - x**2)*Derivative(y(x), x)**2 - y(x)**2,
        "kamke_1.513": 2*x*cos(y(x))**3*Derivative(y(x), x) - sin(y(x))*cos(y(x))**4 + sin(y(x))*Derivative(y(x), x)**2,
        "kamke_1.514": -c*cos(y(x)) + d + (a*cos(y(x)) + b)*Derivative(y(x), x)**2,
        "kamke_1.515": -(x*Derivative(y(x), x) - y(x))**2 + (Derivative(y(x), x)**2 + 1)*f(x**2 + y(x)**2),
        "kamke_1.516": (x**2 + y(x)**2)*(Derivative(y(x), x)**2 + 1)*f(x/sqrt(x**2 + y(x)**2)) - (x*Derivative(y(x), x) - y(x))**2,
        "kamke_1.517": (x**2 + y(x)**2)*(Derivative(y(x), x)**2 + 1)*f(y(x)/sqrt(x**2 + y(x)**2)) - (x*Derivative(y(x), x) - y(x))**2,
        "kamke_1.518": -(-a + y(x))**2*(-b + y(x))**2 + Derivative(y(x), x)**3,
        "kamke_1.519": -(a*y(x)**2 + b*y(x) + c)**2*f(x) + Derivative(y(x), x)**3,
        "kamke_1.520": -y(x) + Derivative(y(x), x)**3 + Derivative(y(x), x),
        "kamke_1.521": x*Derivative(y(x), x) - y(x) + Derivative(y(x), x)**3,
        "kamke_1.522": -(x + 5)*Derivative(y(x), x) + y(x) + Derivative(y(x), x)**3,
        "kamke_1.523": -a*x*Derivative(y(x), x) + x**3 + Derivative(y(x), x)**3,
        "kamke_1.524": y(x)**2 - 2*y(x)*Derivative(y(x), x) + Derivative(y(x), x)**3,
        "kamke_1.525": -a*x*y(x)*Derivative(y(x), x) + 2*a*y(x)**2 + Derivative(y(x), x)**2,
        "kamke_1.526": -x**3*y(x)**3 - (x**2 + x*y(x) + y(x)**2)*Derivative(y(x), x)**2 + (x**3*y(x) + x**2*y(x)**2 + x*y(x)**3)*Derivative(y(x), x) + Derivative(y(x), x)**3,
        "kamke_1.527": -x*y(x)**4*Derivative(y(x), x) - y(x)**5 + Derivative(y(x), x)**3,
        "kamke_1.528": a*b*x + a*Derivative(y(x), x)**2 + b*y(x) + Derivative(y(x), x)**3,
        "kamke_1.529": x*Derivative(y(x), x)**2 - y(x) + Derivative(y(x), x)**3,
        "kamke_1.530": y(x)**2 - y(x)*Derivative(y(x), x)**2 + Derivative(y(x), x)**3,
        "kamke_1.531": -x**3*y(x)**6 - (x**2 + x*y(x)**2 + y(x)**4)*Derivative(y(x), x)**2 + (x**3*y(x)**2 + x**2*y(x)**4 + x*y(x)**6)*Derivative(y(x), x) + Derivative(y(x), x)**2,
        "kamke_1.532": a*Derivative(y(x), x)**3 + b*Derivative(y(x), x)**2 + c*Derivative(y(x), x) - d - y(x),
        "kamke_1.533": a + x*Derivative(y(x), x)**3 - y(x)*Derivative(y(x), x)**2,
        "kamke_1.534": 4*x*Derivative(y(x), x)**3 - x - 6*y(x)*Derivative(y(x), x)**2 + 3*y(x),
        "kamke_1.535": 8*x*Derivative(y(x), x)**3 - 12*y(x)*Derivative(y(x), x)**2 + 9*y(x),
        "kamke_1.536": b*x*(-a**2 + x**2)*Derivative(y(x), x)**2 + b*x + (-a**2 + x**2)*Derivative(y(x), x)**3 + Derivative(y(x), x),
        "kamke_1.537": -2*x**5*y(x) + x**3*Derivative(y(x), x)**3 - 3*x**2*y(x)*Derivative(y(x), x)**2 + (x**6 + 3*x*y(x)**2)*Derivative(y(x), x) - y(x)**3,
        "kamke_1.538": 2*(x*Derivative(y(x), x) + y(x))**3 - y(x)*Derivative(y(x), x),
        "kamke_1.539": -(y(x)*sin(x) - cos(x)**2)*Derivative(y(x), x)**2 - (y(x)*cos(x)**2 + sin(x))*Derivative(y(x), x) + y(x)*sin(x) + sin(x)*Derivative(y(x), x)**3,
        "kamke_1.540": 2*x*Derivative(y(x), x) - x + 2*y(x)*Derivative(y(x), x)**3 - y(x)*Derivative(y(x), x)**2,
        "kamke_1.541": 2*x*Derivative(y(x), x) + y(x)**2*Derivative(y(x), x)**3 - y(x),
        "kamke_1.542": 2*x*Derivative(y(x), x) + 16*y(x)**2*Derivative(y(x), x)**3 - y(x),
        "kamke_1.543": -x**2*y(x) + x*(x**2 + 1)*Derivative(y(x), x) + x*y(x)**2*Derivative(y(x), x)**3 - y(x)**3*Derivative(y(x), x)**2,
        "kamke_1.544": x**7*y(x)**2*Derivative(y(x), x)**3 + 3*x**5*y(x)**4*Derivative(y(x), x) - x**4*y(x)**5 - (3*x**6*y(x)**3 - 1)*Derivative(y(x), x)**2,
        "kamke_1.545": -(-a + y(x))**3*(-b + y(x))**2 + Derivative(y(x), x)**4,
        "kamke_1.546": 3*x + (3*x - 3)*Derivative(y(x), x)**2 - (6*y(x) - 3)*Derivative(y(x), x) + Derivative(y(x), x)**4,
        "kamke_1.547": -4*(x*Derivative(y(x), x) - 2*y(x))**2*y(x) + Derivative(y(x), x)**4,
        "kamke_1.548": -(-a + y(x))**4*(-b + y(x))**3 + Derivative(y(x), x)**6,
        "kamke_1.549": -a**2 + x**2*(Derivative(y(x), x)**2 + 1)**3,
        "kamke_1.550": -a*y(x)**s - b*x**(r*s/(r - s)) + Derivative(y(x), x)**r,
        "kamke_1.551": -(-a + y(x))**(n + 1)*(-b + y(x))**(n - 1)*f(x)**n + Derivative(y(x), x)**n,
        "kamke_1.552": -f(x)*g(y(x)) + Derivative(y(x), x)**n,
        "kamke_1.553": a*Derivative(y(x), x)**m + b*Derivative(y(x), x)**n - y(x),
        "kamke_1.554": -n*x*Derivative(y(x), x) + x**(n - 1)*Derivative(y(x), x)**n + y(x),
        "kamke_1.555": x*Derivative(y(x), x) + sqrt(Derivative(y(x), x)**2 + 1) - y(x),
        "kamke_1.556": x*Derivative(y(x), x)**2 + sqrt(Derivative(y(x), x)**2 + 1) + y(x),
        "kamke_1.557": x*(sqrt(Derivative(y(x), x)**2 + 1) + Derivative(y(x), x)) - y(x),
        "kamke_1.558": a*x*sqrt(Derivative(y(x), x)**2 + 1) + x*Derivative(y(x), x) - y(x),
        "kamke_1.559": -a*x - a*y(x)*Derivative(y(x), x) + sqrt(Derivative(y(x), x)**2 + 1)*y(x),
        "kamke_1.560": a*sqrt(Derivative(y(x), x)**2 + 1)*y(x) - x**2 - 2*x*y(x)*Derivative(y(x), x) + y(x)**2,
        "kamke_1.561": -x*Derivative(y(x), x) + sqrt(Derivative(y(x), x)**2 + 1)*f(x**2 + y(x)**2) + y(x),
        "kamke_1.562": a*(Derivative(y(x), x)**3 + 1)**(1/3) + b*x*Derivative(y(x), x) - y(x),
        "kamke_1.563": a*y(x) + b + x*Derivative(y(x), x) + log(Derivative(y(x), x)),
        "kamke_1.564": a*(x*Derivative(y(x), x) - y(x)) + log(Derivative(y(x), x)),
        "kamke_1.565": -x*y(x) - y(x)*log(y(x)) + y(x)*log(Derivative(y(x), x)) + Derivative(y(x), x),
        "kamke_1.566": -x + sin(Derivative(y(x), x)) + Derivative(y(x), x),
        "kamke_1.567": a*cos(Derivative(y(x), x)) + b*Derivative(y(x), x) + x,
        "kamke_1.568": -y(x) + sin(Derivative(y(x), x))*Derivative(y(x), x)**2,
        "kamke_1.569": (Derivative(y(x), x)**2 + 1)*sin(x*Derivative(y(x), x) - y(x))**2 - 1,
        "kamke_1.570": (a*x + arctan(Derivative(y(x), x)))*(Derivative(y(x), x)**2 + 1) + Derivative(y(x), x),
        "kamke_1.571": a*x**n*f(Derivative(y(x), x)) + x*Derivative(y(x), x) - y(x),
        "kamke_1.572": x*h(Derivative(y(x), x)) + (x*Derivative(y(x), x) - y(x))**n*f(Derivative(y(x), x)) + g(Derivative(y(x), x))*y(x),
        "kamke_1.573": 2*x*Derivative(y(x), x) + f(x*Derivative(y(x), x)**2) - y(x),
        "kamke_1.574": f(x - 3*Derivative(y(x), x)**2/2) - y(x) + Derivative(y(x), x)**3,
        "kamke_1.575": -x**2*Derivative(y(x), x) + x*y(x) + f(x*y(x)*Derivative(y(x), x) - y(x)**2)*Derivative(y(x), x),
        "kamke_1.576": phi(f(x, y(x), Derivative(y(x), x)), g(x, y(x), Derivative(y(x), x))),
        "kamke_1.577": -F(y(x)/(a + x)) + Derivative(y(x), x),
        "kamke_1.578": -2*x - F(-x**2 + y(x)) + Derivative(y(x), x),
        "kamke_1.579": a*x/2 - F(a*x**2/4 + b*x/2 + y(x)) + Derivative(y(x), x),
        "kamke_1.580": -F(y(x)*exp(-b*x))*exp(b*x) + Derivative(y(x), x),
        "kamke_1.581": Derivative(y(x), x) - (x*F((x**2*y(x) + 1/4)/x**2) + 1/2)/x**3,
        "kamke_1.582": Derivative(y(x), x) - (a*x**2*F((a*x*y(x) + 1)/(a*x)) + 1)/(a*x**2),
        "kamke_1.583": -x*(-a*x**2/2 + F(a*x**4/8 + y(x))) + Derivative(y(x), x),
        "kamke_1.584": -2*a/(2*a*F(-4*a*x + y(x)**2) + y(x)) + Derivative(y(x), x),
        "kamke_1.585": -F(-log(x) + log(log(y(x))))*y(x) + Derivative(y(x), x),
        "kamke_1.586": -x*F(y(x)/sqrt(x**2 + 1))/sqrt(x**2 + 1) + Derivative(y(x), x),
        "kamke_1.587": -sqrt(x)*(x**(3/2)/2 + F(-x**3/6 + y(x))) + Derivative(y(x), x),
        "kamke_1.588": -(x + F((-x + y(x))*(x + y(x))))/y(x) + Derivative(y(x), x),
        "kamke_1.589": Derivative(y(x), x) - F((-y(x)*log(x) + 1)/y(x))*y(x)**2/x,
        "kamke_1.590": -x/(F(x**2 + y(x)**2) - y(x)) + Derivative(y(x), x),
        "kamke_1.591": Derivative(y(x), x) - x*F((a*y(x)**2 + b*x**2)/a)/(sqrt(a)*y(x)),
        "kamke_1.592": Derivative(y(x), x) - (sqrt(x) + 6*x**3/5 + F(-2*sqrt(x) - 2*x**3/5 + y(x)))/x,
        "kamke_1.593": -F(y(x)**(3/2) - 3*exp(x)/2)*exp(x)/sqrt(y(x)) + Derivative(y(x), x),
        "kamke_1.594": -x*F((-b + y(x)**2)/x**2)/y(x) + Derivative(y(x), x),
        "kamke_1.595": Derivative(y(x), x) - F((x*y(x)**2 + 1)/x)/(x**2*y(x)),
        "kamke_1.596": Derivative(y(x), x) - (-2*x**2 + x + F(x**2 - x + y(x)))/x,
        "kamke_1.597": -2*a/(x**2*(2*a*F((-4*a + x*y(x)**2)/x) - y(x))) + Derivative(y(x), x),
        "kamke_1.598": Derivative(y(x), x) - (F(y(x)/x) + y(x))/(x - 1),
        "kamke_1.599": -(-x + F(x**2 + y(x)**2))/y(x) + Derivative(y(x), x),
        "kamke_1.600": Derivative(y(x), x) - F((-2*y(x)*log(x) + 1)/y(x))*y(x)**2/x,
        "kamke_1.601": -x*F((-x + y(x))*(x + y(x)))/y(x) + Derivative(y(x), x),
        "kamke_1.602": Derivative(y(x), x) - (x**2*F((x**2 - y(x))/(x**2*y(x))) + 2)*y(x)**2/x**3,
        "kamke_1.603": Derivative(y(x), x) - (2*x*F(y(x) + log(2*x + 1)) + F(y(x) + log(2*x + 1)) - 2)/(2*x + 1),
        "kamke_1.604": Derivative(y(x), x) - 2*y(x)**3/(2*F((4*x*y(x)**2 + 1)/y(x)**2)*y(x) + 1),
        "kamke_1.605": Derivative(y(x), x) + (2*x - F((-x*y(x)/2 + 1)/y(x)))*y(x)**2/(4*x),
        "kamke_1.606": -x*(-x**2*exp(-x**2) + F(-x**2*exp(-x**2)/2 + y(x)) + exp(-x**2)) + Derivative(y(x), x),
        "kamke_1.607": Derivative(y(x), x) - (x**3*F(y(x)/x**2) + 2*y(x))/x,
        "kamke_1.608": Derivative(y(x), x) - sqrt(y(x))/(F((x - y(x))/sqrt(y(x))) + sqrt(y(x))),
        "kamke_1.609": Derivative(y(x), x) - (-3*x**2*y(x) + F(x**3*y(x)))/x**3,
        "kamke_1.610": Derivative(y(x), x) - (x**2*F(y(x)/x) + y(x))/x,
        "kamke_1.611": Derivative(y(x), x) - (-2*x + F(x*(x + y(x))) - y(x))/x,
        "kamke_1.612": -(x*y(x)*exp(-x**2/4)/2 + F(y(x)*exp(-x**2/4)))*exp(x**2/4) + Derivative(y(x), x),
        "kamke_1.613": Derivative(y(x), x) - (x**2*F((-x*log(x) + y(x))/x) + x + y(x))/x,
        "kamke_1.614": -x*(a - 1)*(a + 1)/(a**2*F(-a**2*x**2/2 + x**2/2 + y(x)**2/2) - F(-a**2*x**2/2 + x**2/2 + y(x)**2/2) + y(x)) + Derivative(y(x), x),
        "kamke_1.615": Derivative(y(x), x) - y(x)/(x*(F(x*y(x))*y(x) - 1)),
        "kamke_1.616": Derivative(y(x), x) + (2*x**3*y(x) - x**2 - F(x*(x*y(x) - 1)))/x**4,
        "kamke_1.617": -x*F((y(x)/3 + 1)*exp(3*x**2/2)/y(x))*y(x)**2*exp(-3*x**2/2)/9 + Derivative(y(x), x),
        "kamke_1.618": Derivative(y(x), x) - (x*(y(x) - log(x) - log(y(x) + 1)) + 1)*(y(x) + 1)/(x*y(x)),
        "kamke_1.619": Derivative(y(x), x) - 6*y(x)/(-F(x - y(x)**4/3 - y(x)**3/2 - y(x)**2 - y(x)) + 8*y(x)**4 + 9*y(x)**3 + 12*y(x)**2 + 6*y(x)),
        "kamke_1.620": Derivative(y(x), x) - (x**2 + 2*x*y(x) + y(x)**2 + exp(2*F((-x + y(x))*(x + y(x)))))/(x**2 + 2*x*y(x) + y(x)**2 - exp(2*F((-x + y(x))*(x + y(x))))),
        "kamke_1.621": Derivative(y(x), x) - 1/(sqrt(x) + y(x)),
        "kamke_1.622": Derivative(y(x), x) - 1/(sqrt(3*x + 1) + y(x) + 2),
        "kamke_1.623": -x**2/(x**(3/2) + y(x)) + Derivative(y(x), x),
        "kamke_1.624": -x**(5/3)/(x**(4/3) + y(x)) + Derivative(y(x), x),
        "kamke_1.625": -I*x**2*(-2*sqrt(-x**3 + 6*y(x)) + I)/2 + Derivative(y(x), x),
        "kamke_1.626": -x/(sqrt(x**2 + 1) + y(x)) + Derivative(y(x), x),
        "kamke_1.627": Derivative(y(x), x) - (y(x)*log(x) - 1)**2/x,
        "kamke_1.628": -x*(3*sqrt(x**2 + 3*y(x)) - 2)/3 + Derivative(y(x), x),
        "kamke_1.629": Derivative(y(x), x) - (2*y(x)*log(x) - 1)**2/x,
        "kamke_1.630": Derivative(y(x), x) - exp(b*x)/(y(x)*exp(-b*x) + 1),
        "kamke_1.631": -x**2*(2*sqrt(x**3 - 6*y(x)) + 1)/2 + Derivative(y(x), x),
        "kamke_1.632": Derivative(y(x), x) - exp(x)/(y(x)*exp(-x) + 1),
        "kamke_1.633": Derivative(y(x), x) - exp(2*x/3)/(y(x)*exp(-2*x/3) + 1),
        "kamke_1.634": Derivative(y(x), x) - (x**5*sqrt(4*x**2*y(x) + 1) + 1/2)/x**3,
        "kamke_1.635": -x*(x + 2*sqrt(x**3 - 6*y(x)))/2 + Derivative(y(x), x),
        "kamke_1.636": -(x**2 - log(y(x)))*y(x) + Derivative(y(x), x),
        "kamke_1.637": -x*exp(-x**2)/(y(x)*exp(x**2) + 1) + Derivative(y(x), x),
        "kamke_1.638": -(-log(x) + log(log(y(x))))*y(x) + Derivative(y(x), x),
        "kamke_1.639": -(log(x) - log(log(y(x))))**2*y(x) + Derivative(y(x), x),
        "kamke_1.640": Derivative(y(x), x) - y(x)/(-log(x) + log(log(y(x))) + 1),
        "kamke_1.641": Derivative(y(x), x) - (x**4*sqrt(4*x**2*y(x) + 1) + 1/2)/x**3,
        "kamke_1.642": -(4*a*x - y(x)**2)**2/y(x) + Derivative(y(x), x),
        "kamke_1.643": -x*(3*x*sqrt(x**2 + 3*y(x)) - 2)/3 + Derivative(y(x), x),
        "kamke_1.644": x**2*(a*x - 2*sqrt(a*(a*x**4 + 8*y(x))))/2 + Derivative(y(x), x),
        "kamke_1.645": -(x - log(y(x)))*y(x) + Derivative(y(x), x),
        "kamke_1.646": Derivative(y(x), x) - (x**3/2 + x**2/2 + sqrt(x**3 - 6*y(x)))/(x + 1),
        "kamke_1.647": Derivative(y(x), x) - x*(a*y(x)**2 + b*x**2)**2/(a**(5/2)*y(x)),
        "kamke_1.648": sqrt(a)*x**3*(sqrt(a)*x + sqrt(a) - 2*sqrt(a*x**4 + 8*y(x)))/(2*(x + 1)) + Derivative(y(x), x),
        "kamke_1.649": -x*sqrt(x**2 - 2*x + 8*y(x) + 1) + x/4 + Derivative(y(x), x) - 1/4,
        "kamke_1.650": a/2 - x*sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) + x/2 + Derivative(y(x), x),
        "kamke_1.651": Derivative(y(x), x) - (x**2 + log(y(x)))*y(x)/x,
        "kamke_1.652": -(2*a + x*sqrt(4*a*x - y(x)**2))/y(x) + Derivative(y(x), x),
        "kamke_1.653": -x*sqrt(x**2 - 4*x + 4*y(x)) + x/2 + Derivative(y(x), x) - 1,
        "kamke_1.654": Derivative(y(x), x) - (-2*x**2/3 - 2*x/3 + sqrt(x**2 + 3*y(x)))/(x + 1),
        "kamke_1.655": Derivative(y(x), x) - y(x)**3*exp(-4*x/3)/(y(x)*exp(-2*x/3) + 1),
        "kamke_1.656": Derivative(y(x), x) - (x**3 + log(y(x)))*y(x)/x,
        "kamke_1.657": -x**2*sqrt(x**2 - 2*x + 8*y(x) + 1) + x/4 + Derivative(y(x), x) - 1/4,
        "kamke_1.658": Derivative(y(x), x) - (-x**2/4 + sqrt(x**2 - 2*x + 8*y(x) + 1) + 1/4)/(x + 1),
        "kamke_1.659": a*x/2 + b/2 - x*sqrt(a**2*x**2 + 2*a*b*x + 4*a*y(x) + b**2 - 4*c) + Derivative(y(x), x),
        "kamke_1.660": a/2 - x**2*sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) + x/2 + Derivative(y(x), x),
        "kamke_1.661": a*x/2 + b/2 - x**2*sqrt(a**2*x**2 + 2*a*b*x + 4*a*y(x) + b**2 - 4*c) + Derivative(y(x), x),
        "kamke_1.662": -x**2*sqrt(x**2 + 2*x - 4*y(x) + 1) - x/2 + Derivative(y(x), x) - 1/2,
        "kamke_1.663": -(2*a + x**2*sqrt(4*a*x - y(x)**2))/y(x) + Derivative(y(x), x),
        "kamke_1.664": -x**2*sqrt(x**2 - 4*x + 4*y(x)) + x/2 + Derivative(y(x), x) - 1,
        "kamke_1.665": -sqrt(a)*(-sqrt(a)*x**4/2 - sqrt(a)*x**3/2 + sqrt(a*x**4 + 8*y(x)))/(x + 1) + Derivative(y(x), x),
        "kamke_1.666": -(x**3 + x**2 - log(y(x)) + 1)*y(x) + Derivative(y(x), x),
        "kamke_1.667": Derivative(y(x), x) - y(x)**3*exp(-2*b*x)/(y(x)*exp(-b*x) + 1),
        "kamke_1.668": Derivative(y(x), x) - y(x)**3*exp(-2*x)/(y(x)*exp(-x) + 1),
        "kamke_1.669": -(-2*y(x)**(3/2) + 3*exp(x))**2*exp(x)/(4*sqrt(y(x))) + Derivative(y(x), x),
        "kamke_1.670": -I*x*(-2*sqrt(-x**2 + 4*log(a) + 4*log(y(x))) + I)*y(x)/2 + Derivative(y(x), x),
        "kamke_1.671": Derivative(y(x), x) - (x*y(x)**2 + 1)**2/(x**4*y(x)),
        "kamke_1.672": -x**2*(3*x + sqrt(-9*x**4 + 4*y(x)**3))/y(x)**2 + Derivative(y(x), x),
        "kamke_1.673": Derivative(y(x), x) - (x**2*cos(2*y(x))/2 + x**2/2 - sin(2*y(x))/2)/x,
        "kamke_1.674": Derivative(y(x), x) - (-x**2/2 + x/2 + sqrt(x**2 - 4*x + 4*y(x)) + 1)/(x + 1),
        "kamke_1.675": Derivative(y(x), x) - (a*x**4 + a*x**3*exp(x) + a*x**3 - x**2*y(x)**2 - x*y(x)**2*exp(x) - x*y(x)**2 + y(x))/x,
        "kamke_1.676": Derivative(y(x), x) - (x**6*sqrt(4*x**2*y(x) + 1) + x/2 + 1/2)/(x**3*(x + 1)),
        "kamke_1.677": Derivative(y(x), x) - (a*x**4 + a*x**3*log(x + 1) + a*x**3 - x**2*y(x)**2 - x*y(x)**2*log(x + 1) - x*y(x)**2 + y(x))/x,
        "kamke_1.678": -x**2*(2*x*sqrt(x**3 - 6*y(x)) + x + 1)/(2*(x + 1)) + Derivative(y(x), x),
        "kamke_1.679": Derivative(y(x), x) - (x**4 + x**3*log(x) + x**3 + 7*x**2*y(x)**2 + 7*x*y(x)**2*log(x) + 7*x*y(x)**2 + y(x))/x,
        "kamke_1.680": Derivative(y(x), x) - (x**2/2 + x + sqrt(x**2 + 2*x - 4*y(x) + 1) + 1/2)/(x + 1),
        "kamke_1.681": Derivative(y(x), x) - (a*x**2*y(x)**2 + a*x*y(x)**2*log(1/x) + a*x*y(x)**2 + b*x**4 + b*x**3*log(1/x) + b*x**3 + y(x))/x,
        "kamke_1.682": -2*a/(x*(-8*a**2 + 2*a*x*y(x)**2 - x*y(x))) + Derivative(y(x), x),
        "kamke_1.683": Derivative(y(x), x) - (x**4*y(x)*log(x*(x + 1)) - x**3*log(x*(x + 1)) - 1)*y(x)/x,
        "kamke_1.684": Derivative(y(x), x) - (x**2*sqrt(x**2 + y(x)**2) + y(x))/x,
        "kamke_1.685": Derivative(y(x), x) - (x**3*log((x - 1)*(x + 1)) + 7*x*y(x)**2*log((x - 1)*(x + 1)) + y(x))/x,
        "kamke_1.686": -x*y(x)**3*exp(2*x**2)/(y(x)*exp(x**2) + 1) + Derivative(y(x), x),
        "kamke_1.687": Derivative(y(x), x) - (-x**3*log((x + 1)/(x - 1)) + x*y(x)**2*log((x + 1)/(x - 1)) + y(x))/x,
        "kamke_1.688": Derivative(y(x), x) - (x**3*exp((x + 1)/(x - 1)) + x*y(x)**2*exp((x + 1)/(x - 1)) + y(x))/x,
        "kamke_1.689": Derivative(y(x), x) - (-x**3*exp(x + 1) + x*y(x)**2*exp(x + 1) + x*y(x) - y(x))/(x*(x - 1)),
        "kamke_1.690": Derivative(y(x), x) - (x**3*sqrt(x**2 - 2*x + 8*y(x) + 1) - x**2/4 + 1/4)/(x + 1),
        "kamke_1.691": Derivative(y(x), x) - (x**3*cos(2*y(x))/2 + x**3/2 - sin(2*y(x))/2)/x,
        "kamke_1.692": Derivative(y(x), x) - (x**3*sqrt(x**2 + y(x)**2) + y(x))/x,
        "kamke_1.693": -(y(x)**3*exp(-3*b*x) + y(x)**2*exp(-2*b*x) + 1)*exp(b*x) + Derivative(y(x), x),
        "kamke_1.694": Derivative(y(x), x) - (x**3*sqrt(4*x**2*y(x) + 1) + x/2 + 1/2)/(x**3*(x + 1)),
        "kamke_1.695": Derivative(y(x), x) - (x**4 + x**3 + x**2*y(x)**2 + x*y(x)**2 + y(x)*log(x - 1))/(x*log(x - 1)),
        "kamke_1.696": Derivative(y(x), x) - (x**3*exp(x + 1) + 7*x*y(x)**2*exp(x + 1) + y(x)*log(x - 1))/(x*log(x - 1)),
        "kamke_1.697": -(y(x)**3*exp(-2*x) + y(x)**2*exp(-4*x/3) + 1)*exp(2*x/3) + Derivative(y(x), x),
        "kamke_1.698": -(y(x)**3*exp(-3*x) + y(x)**2*exp(-2*x) + 1)*exp(x) + Derivative(y(x), x),
        "kamke_1.699": -x*(3*x**2*sqrt(x**2 + 3*y(x)) - 2*x - 2)/(3*(x + 1)) + Derivative(y(x), x),
        "kamke_1.700": Derivative(y(x), x) - 1/(x*(x*y(x)**2 + x + 1)*y(x)),
        "kamke_1.701": Derivative(y(x), x) - (x**4*log(x) + x**4 - 2*x**2*y(x)*log(x) - 2*x**2*y(x) + 2*x*exp(x) - 2*x + y(x)**2*log(x) + y(x)**2 - log(x) - 1)/(exp(x) - 1),
        "kamke_1.702": Derivative(y(x), x) - (-x**3*log(x) - x**3 - x*y(x)**2*log(x) - x*y(x)**2 + x*y(x) - y(x)*exp(x))/(x*(x - exp(x))),
        "kamke_1.703": Derivative(y(x), x) - (x**3*y(x) + x**2*y(x)*log(x) - x**2 - x*log(x) - x + 1)*y(x)/(x*(x - 1)),
        "kamke_1.704": Derivative(y(x), x) - (2*a*x**3*y(x)**2 + 2*b*x**5 + x*y(x)*log(x) - y(x))/(x*(x*log(x) - 1)),
        "kamke_1.705": Derivative(y(x), x) - (x**4 + x**3 + x + log(y(x)))*y(x)/x,
        "kamke_1.706": -x*(y(x) + 1)**2*(-log(x)/4 + log(y(x) - 1)/8 - log(y(x) + 1)/8) + Derivative(y(x), x),
        "kamke_1.707": -x*(y(x) + 1)**2*(2*log(x) - log(y(x) - 1) + log(y(x) + 1))**2/16 + Derivative(y(x), x),
        "kamke_1.708": -(4*a*x - y(x)**2)**3/((4*a*x - y(x)**2 - 1)*y(x)) + Derivative(y(x), x),
        "kamke_1.709": Derivative(y(x), x) - (2*a*x + 2*a + x**3*sqrt(4*a*x - y(x)**2))/((x + 1)*y(x)),
        "kamke_1.710": Derivative(y(x), x) - (2*x**3 + 4*x**2*y(x) + 2*x*y(x)**2 + 2*x + exp(1/x) - log(x))/(-exp(1/x) + log(x)),
        "kamke_1.711": Derivative(y(x), x) - (-x*log(y(x)) - log(y(x)) + 1)*y(x)/(x + 1),
        "kamke_1.712": Derivative(y(x), x) - (x**3*sqrt(x**2 + 2*x - 4*y(x) + 1) + x**2/2 + x + 1/2)/(x + 1),
        "kamke_1.713": Derivative(y(x), x) - (-a**2 - a*b*sqrt(x) - a*b*y(x) + a*b + b**2*x + b**2)/(a*(-a*sqrt(x) - a*y(x) + a + b*x + b)),
        "kamke_1.714": Derivative(y(x), x) + (x**3*y(x) + x**2*y(x)*log(x) - x**2 - x*log(x) + exp(x) - log(1/x))*y(x)/(x*(exp(x) - log(1/x))),
        "kamke_1.715": Derivative(y(x), x) - (x**3*sqrt(x**2 - 4*x + 4*y(x)) - x**2/2 + x/2 + 1)/(x + 1),
        "kamke_1.716": Derivative(y(x), x) - (3*x**4 + 3*x**3 + sqrt(9*x**4 - 4*y(x)**3))/((x + 1)*y(x)**2),
        "kamke_1.717": Derivative(y(x), x) - (-a*x/2 - a/2 - x**2/2 - x/2 + sqrt(a**2 + 2*a*x + x**2 + 4*y(x)))/(x + 1),
        "kamke_1.718": -x*(y(x)**3*exp(3*x**2) + y(x)**2*exp(2*x**2) + 1)*exp(-x**2) + Derivative(y(x), x),
        "kamke_1.719": Derivative(y(x), x) - (x**2*y(x)*log(2*x) - x*log(2*x) - exp(x))*y(x)*exp(-x)/x,
        "kamke_1.720": -x**3*(3*x + sqrt(9*x**4 - 4*y(x)**3) + 3)/((x + 1)*y(x)**2) + Derivative(y(x), x),
        "kamke_1.721": -sqrt(x)*(x**(3/2)/2 + x**6/36 - x**3*y(x)/3 + y(x)**2) + Derivative(y(x), x),
        "kamke_1.722": Derivative(y(x), x) + y(x)**3/(x*(2*y(x)*log(x) - y(x) - 1)),
        "kamke_1.723": -2*a/(32*a**3*x**2 - 16*a**2*x*y(x)**2 + 2*a*y(x)**4 + y(x)) + Derivative(y(x), x),
        "kamke_1.724": Derivative(y(x), x) + y(x)**3/(x*(y(x)*log(x) - y(x) - 1)),
        "kamke_1.725": -(x**2*log(2*x) + 2*x*y(x)*log(2*x) + y(x)**2*log(2*x) - log(x) + log(2*x))/log(x) + Derivative(y(x), x),
        "kamke_1.726": Derivative(y(x), x) - (a**2 - a*b*sqrt(x) - a*b*y(x) - b**2*x + b*c)/(a*(a*sqrt(x) + a*y(x) + b*x - c)),
        "kamke_1.727": Derivative(y(x), x) - (2*x + y(x) + 2)*y(x)/((x + 1)*(2*x + log(y(x)) - 1)),
        "kamke_1.728": Derivative(y(x), x) - (x**3 + 3*y(x)**2)*y(x)/(x*(x + 6*y(x)**2)),
        "kamke_1.729": Derivative(y(x), x) - (x - y(x))*y(x)/(x*(x - y(x)**3)),
        "kamke_1.730": -(2*y(x)**(3/2) - 3*exp(x))**3*exp(x)/(4*(2*y(x)**(3/2) - 3*exp(x) + 2)*sqrt(y(x))) + Derivative(y(x), x),
        "kamke_1.731": Derivative(y(x), x) - (2*y(x) + 1)/(x*(2*x*y(x)**3 + x*y(x)**2 - 2)),
        "kamke_1.732": Derivative(y(x), x) - (-a*x/2 - a/2 + x**3*sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) - x**2/2 - x/2)/(x + 1),
        "kamke_1.733": -(x**4*log(2*x) - 2*x**2*y(x)*log(2*x) + 2*x*sin(x) + y(x)**2*log(2*x) - log(2*x))/sin(x) + Derivative(y(x), x),
        "kamke_1.734": Derivative(y(x), x) - (x**3 - x*log(y(x)) - log(y(x)))*y(x)/(x + 1),
        "kamke_1.735": Derivative(y(x), x) - (2*y(x)*log(x) - 1)**3/(x*(2*y(x)*log(x) - y(x) - 1)),
        "kamke_1.736": Derivative(y(x), x) - (x**4 - 2*x**2*y(x) + 2*x**2 + 2*x + y(x)**2 - 1)/(x + 1),
        "kamke_1.737": -x*(2*x**3 - 2*x*y(x) + x - 1)/(x**2 - y(x)) + Derivative(y(x), x),
        "kamke_1.738": -2*a/(32*a**3 - 16*a**2*x*y(x)**2 + 2*a*x**2*y(x)**4 - x**2*y(x)) + Derivative(y(x), x),
        "kamke_1.739": Derivative(y(x), x) - (2*y(x) + 1)/(x*(2*x*y(x)**2 + x*y(x) - 2)),
        "kamke_1.740": -(x**4 - 2*x**2*y(x)**2 + x + y(x)**4)/y(x) + Derivative(y(x), x),
        "kamke_1.741": Derivative(y(x), x) - x*(a*y(x)**2 + b*x**2)**3/(a**(5/2)*(a*y(x)**2 + a + b*x**2)*y(x)),
        "kamke_1.742": Derivative(y(x), x) + (x - cos(y(x)) + 1)*cos(y(x))/((x + 1)*(x*sin(y(x)) - 1)),
        "kamke_1.743": I*(x**4 + 8*x**2*y(x)**2 + 8*I*x + 16*y(x)**4)/(32*y(x)) + Derivative(y(x), x),
        "kamke_1.744": -x/(x**4 + 2*x**2*y(x)**2 + y(x)**4 - y(x)) + Derivative(y(x), x),
        "kamke_1.745": Derivative(y(x), x) - (y(x)*log(x) - 1)**3/(x*(y(x)*log(x) - y(x) - 1)),
        "kamke_1.746": I*(x**4 + 2*x**2*y(x)**2 + I*x + y(x)**4)/y(x) + Derivative(y(x), x),
        "kamke_1.747": Derivative(y(x), x) + (-x**2*y(x)*log(2*x) + x*log(2*x) + tan(x))*y(x)/(x*tan(x)),
        "kamke_1.748": Derivative(y(x), x) - (x + y(x))*y(x)/(x*(x + y(x)**3)),
        "kamke_1.749": -x*(x - y(x))**2*(x + y(x))**2/y(x) + Derivative(y(x), x),
        "kamke_1.750": Derivative(y(x), x) - (x**2 + 3*y(x)**2)*y(x)/(x*(x + 6*y(x)**2)),
        "kamke_1.751": Derivative(y(x), x) - (x**4 + x*log(y(x)) + log(y(x)))*y(x)/(x*(x + 1)),
        "kamke_1.752": Derivative(y(x), x) - (x**3*cos(y(x)) - x - 1)*cos(y(x))/((x + 1)*(x*sin(y(x)) - 1)),
        "kamke_1.753": Derivative(y(x), x) - (x**4*log(y(x)) + x + 1)*y(x)*log(y(x))/(x*(x + 1)),
        "kamke_1.754": Derivative(y(x), x) - (x**3 + x*y(x)**2 + x*y(x) + y(x)**3)/x**2,
        "kamke_1.755": Derivative(y(x), x) - y(x)**(3/2)/(x**2 - 2*x*y(x) + y(x)**(3/2) + y(x)**2),
        "kamke_1.756": Derivative(y(x), x) - (x**6 + 2*x**3*y(x) + x**2*y(x)**2 + y(x)**3)/x**4,
        "kamke_1.757": Derivative(y(x), x) - (x**3 + 2*x**2 - 4*x*y(x) - 4*x - 8)/(2*x**2 + 4*x - 8*y(x) - 8),
        "kamke_1.758": Derivative(y(x), x) - (x**3*y(x) + 2*x + 2)*y(x)/((x + 1)*(2*x + log(y(x)) - 1)),
        "kamke_1.759": I*x*(x**8 + 18*x**4*y(x)**2 + 54*I*x**2 + 81*y(x)**4)/(243*y(x)) + Derivative(y(x), x),
        "kamke_1.760": Derivative(y(x), x) - (x*y(x)**2 + 1)**3/(x**4*(x*y(x)**2 + x + 1)*y(x)),
        "kamke_1.761": Derivative(y(x), x) - (-x**3 + 4*x**2 - 4*x*y(x) - 4*x + 8)/(2*x**2 - 8*x + 8*y(x) + 8),
        "kamke_1.762": Derivative(y(x), x) - (-x*log(y(x)) + x - log(y(x)))*y(x)/(x*(x + 1)),
        "kamke_1.763": Derivative(y(x), x) - (x*log(y(x)) + x + log(y(x)))*y(x)/(x*(x + 1)),
        "kamke_1.764": Derivative(y(x), x) - (x**4 - x*log(y(x)) - log(y(x)))*y(x)/(x*(x + 1)),
        "kamke_1.765": Derivative(y(x), x) - (x*y(x)*log((x - 1)*(x + 1)/x) - log((x - 1)*(x + 1)/x) - 1)*y(x)/x,
        "kamke_1.766": Derivative(y(x), x) - (x**2*y(x)*log((x - 1)*(x + 1)/x) - x*log((x - 1)*(x + 1)/x) - log(x))*y(x)/(x*log(x)),
        "kamke_1.767": Derivative(y(x), x) - (-x**3 + 2*x**2 - 8*x*y(x) - 8*x + 32)/(4*x**2 - 8*x + 32*y(x) + 32),
        "kamke_1.768": Derivative(y(x), x) - (y(x) + 1)*y(x)/(x*(x*y(x) - y(x) - 1)),
        "kamke_1.769": I*x*(x**8 + 8*x**4*y(x)**2 + 16*I*x**2 + 16*y(x)**4)/(32*y(x)) + Derivative(y(x), x),
        "kamke_1.770": Derivative(y(x), x) - 2*y(x)**6/(32*x**2*y(x)**4 + 16*x*y(x)**2 + y(x)**3 + 2),
        "kamke_1.771": Derivative(y(x), x) - (-a**2*x**3 - 2*a*b*x**2 - 4*a*x*y(x) - 4*a*x + 8)/(2*a*x**2 + 4*b*x + 8*y(x) + 8),
        "kamke_1.772": Derivative(y(x), x) - (x*log(y(x)) + x + 1)*y(x)*log(y(x))/(x*(x + 1)),
        "kamke_1.773": Derivative(y(x), x) - (x*y(x) + x + y(x)**2)/((x - 1)*(x + y(x))),
        "kamke_1.774": Derivative(y(x), x) - (-2*a*x**2 - x**3 - 4*x*y(x) - 4*x + 8)/(4*a*x + 2*x**2 + 8*y(x) + 8),
        "kamke_1.775": -(x + sqrt(y(x)) - y(x))/(x + sqrt(y(x)) - y(x) + 1) + Derivative(y(x), x),
        "kamke_1.776": Derivative(y(x), x) - (x**2*y(x)*log((x**2 + 1)/x) - x*log((x**2 + 1)/x) - log(1/x))*y(x)/(x*log(1/x)),
        "kamke_1.777": Derivative(y(x), x) - (y(x) + 1)*y(x)/(x*(x*y(x)**4 - y(x) - 1)),
        "kamke_1.778": Derivative(y(x), x) - (x**9*y(x)**3 + x**6*y(x)**2 - 3*x**2*y(x) + 1)/x**3,
        "kamke_1.779": Derivative(y(x), x) - (x**3*y(x) + x**3 + x*y(x)**2 + y(x)**3)/(x**3*(x - 1)),
        "kamke_1.780": Derivative(y(x), x) - (x*sqrt(x**2 + y(x)**2) + x*y(x) + y(x))/(x*(x + 1)),
        "kamke_1.781": Derivative(y(x), x) - (x**4 + x**3 + x + 3*y(x)**2)*y(x)/(x*(x + 6*y(x)**2)),
        "kamke_1.782": Derivative(y(x), x) - (x**2*y(x)*log((x**2 + 1)/x) - x*log((x**2 + 1)/x) - tanh(1/x))*y(x)/(x*tanh(1/x)),
        "kamke_1.783": Derivative(y(x), x) + (-x**2*y(x)*log(2*x) + x*log(2*x) + tanh(x))*y(x)/(x*tanh(x)),
        "kamke_1.784": -(x**2*log(x) + 2*x*y(x)*log(x) + y(x)**2*log(x) + log(x) - sinh(x))/sinh(x) + Derivative(y(x), x),
        "kamke_1.785": -(x**2*sinh(x) + 2*x*y(x)*sinh(x) + y(x)**2*sinh(x) - log(x) + sinh(x))/log(x) + Derivative(y(x), x),
        "kamke_1.786": Derivative(y(x), x) - (a*x*y(x)**2*cosh(x) + b*x**3*cosh(x) + y(x)*log(x))/(x*log(x)),
        "kamke_1.787": -x*(2*x**4 - 2*x**2*y(x) + x**2 - x - 1)/((x + 1)*(x**2 - y(x))) + Derivative(y(x), x),
        "kamke_1.788": Derivative(y(x), x) + (-x**2*y(x)*coth(x + 1) + x*coth(x + 1) + log(x - 1))*y(x)/(x*log(x - 1)),
        "kamke_1.789": -(x**2*coth(x + 1) + 2*x*y(x)*coth(x + 1) + y(x)**2*coth(x + 1) - log(x - 1) + coth(x + 1))/log(x - 1) + Derivative(y(x), x),
        "kamke_1.790": -(x**4*coth((x + 1)/(x - 1)) - 2*x**2*y(x)*coth((x + 1)/(x - 1)) + 2*x*log(1/(x - 1)) + y(x)**2*coth((x + 1)/(x - 1)) - coth((x + 1)/(x - 1)))/log(1/(x - 1)) + Derivative(y(x), x),
        "kamke_1.791": Derivative(y(x), x) - (x**5 + x**4 - 2*x**3*y(x) - 2*x**2*y(x) + 2*x**2*cosh(1/(x - 1)) + x*y(x)**2 - 2*x*cosh(1/(x - 1)) - x + y(x)**2 - 1)/((x - 1)*cosh(1/(x - 1))),
        "kamke_1.792": Derivative(y(x), x) - (x**3*y(x) + x**2*y(x) - x**2 - x*cosh(1/(x + 1)) - x + cosh(1/(x + 1)))*y(x)/(x*(x - 1)*cosh(1/(x + 1))),
        "kamke_1.793": Derivative(y(x), x) + (x*y(x) + 1)*y(x)/(x*(x*y(x) - y(x) + 1)),
        "kamke_1.794": Derivative(y(x), x) - y(x)/(x*(x**3*y(x)**4 + x**2*y(x)**3 + y(x) - 1)),
        "kamke_1.795": Derivative(y(x), x) - (a**3 + 3*a**2*x + 3*a*x**2 + a*y(x)**2 + x**3 + x*y(x)**2 + y(x)**3)/(a + x)**3,
        "kamke_1.796": -x*y(x)**3*exp(-3*x**2/2)/(3*(y(x)*exp(3*x**2/2) + 3*y(x) + 3*exp(3*x**2/2))) + Derivative(y(x), x),
        "kamke_1.797": Derivative(y(x), x) - (x**3*y(x)*cosh((x + 1)/(x - 1)) + x**2*y(x)*cosh((x + 1)/(x - 1)) - x**2*cosh((x + 1)/(x - 1)) - x*cosh((x + 1)/(x - 1)) - 1)*y(x)/x,
        "kamke_1.798": Derivative(y(x), x) - (x + y(x) + 1)*y(x)/((x + 1)*(x + 2*y(x)**3 + y(x))),
        "kamke_1.799": Derivative(y(x), x) - (x**3*y(x)*exp((x + 1)/(x - 1)) + x**2*y(x)*exp((x + 1)/(x - 1)) - x**2*exp((x + 1)/(x - 1)) - x*exp((x + 1)/(x - 1)) - 1)*y(x)/x,
        "kamke_1.800": Derivative(y(x), x) - (-b**3 + 6*b**2*x - 12*b*x**2 - 4*b*y(x)**2 + 8*x**3 + 8*x*y(x)**2 + 8*y(x)**3)/(-b + 2*x)**3,
        "kamke_1.801": -(x*y(x)*exp(-x**2/4)/2 + y(x)**3*exp(-3*x**2/4) + y(x)**2*exp(-x**2/2) + 1)*exp(x**2/4) + Derivative(y(x), x),
        "kamke_1.802": Derivative(y(x), x) - (_F1(y(x) + 1/x) + 1/x)/x,
        "kamke_1.803": Derivative(y(x), x) - _F1(y(x)**2 - 2*log(x))/(x*sqrt(y(x)**2)),
        "kamke_1.804": Derivative(y(x), x) - (x**4*cos(2*y(x))/2 + x**4/2 - x*sin(2*y(x))/2 - sin(2*y(x))/2)/(x*(x + 1)),
        "kamke_1.805": Derivative(y(x), x) - (x**4*sqrt(x**2 + y(x)**2) + x*y(x) + y(x))/(x*(x + 1)),
        "kamke_1.806": Derivative(y(x), x) - (-x*sin(2*y(x))/2 + x*cos(2*y(x))/2 + x/2 - sin(2*y(x))/2)/(x*(x + 1)),
        "kamke_1.807": Derivative(y(x), x) + 1/(-x - _F1(y(x) - log(x))*y(x)*exp(y(x))),
        "kamke_1.808": Derivative(y(x), x) - (y(x) + 1)*(2*y(x) + 1)/(x*(2*x*y(x) + x - 2*y(x) - 2)),
        "kamke_1.809": Derivative(y(x), x) - (64*x**3 - 240*x**2 + 64*x*y(x)**2 + 300*x + 64*y(x)**3 - 80*y(x)**2 - 125)/(4*x - 5)**3,
        "kamke_1.810": Derivative(y(x), x) - (x**2*log(x)**2 - 2*x*y(x)*log(x) + x + y(x)**2 + y(x))/x,
        "kamke_1.811": Derivative(y(x), x) - (x**4 + x**3*exp(y(x)) + x*y(x) - x*log(x + exp(y(x))) + x + y(x)*exp(y(x)) - exp(y(x))*log(x + exp(y(x))))/x**2,
        "kamke_1.812": -x**3*sqrt(x**3 - 6*y(x)) - x**2*sqrt(x**3 - 6*y(x)) - x**2/2 - sqrt(x**3 - 6*y(x)) + Derivative(y(x), x),
        "kamke_1.813": -sqrt(a)*(-sqrt(a)*x**3/2 + x**3*sqrt(a*x**4 + 8*y(x)) + x**2*sqrt(a*x**4 + 8*y(x)) + sqrt(a*x**4 + 8*y(x))) + Derivative(y(x), x),
        "kamke_1.814": Derivative(y(x), x) - (x**7*y(x)**2 - 3*x**3*y(x) - 3)*y(x)/(x*(x**3*y(x) + 1)),
        "kamke_1.815": -x*(y(x) + 3)**3*exp(3*x**2)/(81*(y(x)*exp(3*x**2/2) + 3*y(x) + 3*exp(3*x**2/2))) + Derivative(y(x), x),
        "kamke_1.816": -x*(x - y(x))**3*(x + y(x))**3/((x**2 - y(x)**2 - 1)*y(x)) + Derivative(y(x), x),
        "kamke_1.817": Derivative(y(x), x) - (x**3*log(x)*cos(2*y(x))/2 + x**3*log(x)/2 - cos(y(x)))/(x*log(x)*sin(y(x))),
        "kamke_1.818": Derivative(y(x), x) - y(x)/(x*(x*y(x)**4 + x*y(x)**3 + x*y(x) - 1)),
        "kamke_1.819": -x**3*sqrt(x**2 + 3*y(x)) - x**2*sqrt(x**2 + 3*y(x)) + 2*x/3 - sqrt(x**2 + 3*y(x)) + Derivative(y(x), x),
        "kamke_1.820": Derivative(y(x), x) - (x**2*log(x)*cos(2*y(x))/2 + x**2*log(x)/2 - cos(y(x)))/(x*log(x)*sin(y(x))),
        "kamke_1.821": Derivative(y(x), x) - (x*y(x) + 1)*y(x)/(x*(x**3*y(x)**4 - x*y(x) - 1)),
        "kamke_1.822": -x*(x**4*exp(-2*x**2)/4 - x**2*y(x)*exp(-x**2) - x**2*exp(-x**2) + y(x)**2 + exp(-x**2)) + Derivative(y(x), x),
        "kamke_1.823": Derivative(y(x), x) - (x + y(x))*y(x)/(x*(x + y(x)**4 + y(x)**3 + y(x))),
        "kamke_1.824": Derivative(y(x), x) - (x**3 + x**2*y(x) + y(x)**2)*y(x)/(x**2*(x - 1)*(x + y(x))),
        "kamke_1.825": -x*(x**2*(x**2 + 1)**(3/2) + x**2*y(x)**3 + (x**2 + 1)**(3/2)*y(x)**2 + (x**2 + 1)**(3/2) + y(x)**3)/(x**2 + 1)**3 + Derivative(y(x), x),
        "kamke_1.826": Derivative(y(x), x) - (3*x*y(x)**2 + x + 3*y(x)**2)*y(x)/(x*(x + 1)*(x + 6*y(x)**2)),
        "kamke_1.827": Derivative(y(x), x) - (-x**3*sqrt(x**2 + y(x)**2) + x**2*sqrt(x**2 + y(x)**2)*y(x) + y(x))/x,
        "kamke_1.828": Derivative(y(x), x) - (y(x) + 1)*(2*y(x) + 1)/(x*(2*x*y(x)**4 + x*y(x)**3 - 2*y(x) - 2)),
        "kamke_1.829": Derivative(y(x), x) - (x**6*sqrt(4*x**2*y(x) + 1) + x**5*sqrt(4*x**2*y(x) + 1) + x**3*sqrt(4*x**2*y(x) + 1) + 1/2)/x**3,
        "kamke_1.830": Derivative(y(x), x) - (x - y(x))*y(x)/(x*(x - y(x)**4 - y(x)**3 - y(x))),
        "kamke_1.831": -(2*a + x**3*sqrt(4*a*x - y(x)**2) + x**2*sqrt(4*a*x - y(x)**2) + sqrt(4*a*x - y(x)**2))/y(x) + Derivative(y(x), x),
        "kamke_1.832": Derivative(y(x), x) - (x + y(x) + 1)*y(x)/((x + 1)*(x + y(x)**4 + y(x)**3 + y(x)**2)),
        "kamke_1.833": Derivative(y(x), x) - (-x**4*sqrt(x**2 + y(x)**2) + x**3*sqrt(x**2 + y(x)**2)*y(x) + y(x))/x,
        "kamke_1.834": Derivative(y(x), x) - (x**4 + 3*x*y(x)**2 + 3*y(x)**2)*y(x)/(x*(x + 1)*(x + 6*y(x)**2)),
        "kamke_1.835": Derivative(y(x), x) + 1/(-x*(y(x)**3)**(2/3) - x*(y(x)**3)**(1/3)*_F1(y(x)**3 - 3*log(x))),
        "kamke_1.836": Derivative(y(x), x) - (x - y(x))*(y(x) + 1)*y(x)/(x*(x*y(x) + x - y(x))),
        "kamke_1.837": Derivative(y(x), x) + 1/(-(y(x)**3)**(2/3)*log(x) - (y(x)**3)**(1/3)*_F1(y(x)**3 + 3*Ei(-log(x)))*log(x)),
        "kamke_1.838": Derivative(y(x), x) - (8*x**(7/2)/5 - 4*sqrt(x)*y(x) + sqrt(x) + 4*x**6/25 - 4*x**3*y(x)/5 + 6*x**3/5 + 4*x + y(x)**2)/x,
        "kamke_1.839": Derivative(y(x), x) - (x**2 + x*exp(-y(x)/x) + y(x)*exp(-y(x)/x))*exp(y(x)/x)/x,
        "kamke_1.840": Derivative(y(x), x) - (x**3 + x*exp(-y(x)/x) + y(x)*exp(-y(x)/x))*exp(y(x)/x)/x,
        "kamke_1.841": Derivative(y(x), x) - (a**(5/2)*y(x)**4 - 2*a**(3/2)*b*x**2*y(x)**2 + 2*a**(3/2)*c*y(x)**2 + sqrt(a)*b**2*x**4 - 2*sqrt(a)*b*c*x**2 + sqrt(a)*c**2 + b*x**3)/(a*x**2*y(x)),
        "kamke_1.842": Derivative(y(x), x) - (x**2*y(x)**2*log(x) + 2*x**2*y(x)*log(x)**2 + x**2*log(x)**3 + y(x))/(x*log(x)),
        "kamke_1.843": Derivative(y(x), x) - (x**3*y(x)**2*log(x) + 2*x**3*y(x)*log(x)**2 + x**3*log(x)**3 + y(x))/(x*log(x)),
        "kamke_1.844": Derivative(y(x), x) - (x + y(x))*(y(x) + 1)*y(x)/(x*(x*y(x) + x + y(x))),
        "kamke_1.845": -(x**3*sqrt(-9*x**4 + 4*y(x)**3) + 3*x**3 + x**2*sqrt(-9*x**4 + 4*y(x)**3) + sqrt(-9*x**4 + 4*y(x)**3))/y(x)**2 + Derivative(y(x), x),
        "kamke_1.846": Derivative(y(x), x) - 1/(-x**2*(1 + 1/y(x))*_F1(x*(1 + 1/y(x))) + x**2*_F1(x*(1 + 1/y(x))) + x*(1 + 1/y(x)) - x),
        "kamke_1.847": -x**3*sqrt(x**2 + 2*x - 4*y(x) + 1) - x**2*sqrt(x**2 + 2*x - 4*y(x) + 1) - x/2 - sqrt(x**2 + 2*x - 4*y(x) + 1) + Derivative(y(x), x) - 1/2,
        "kamke_1.848": -_F1(y(x) - log(sinh(x))) + Derivative(y(x), x) - cosh(x)/sinh(x),
        "kamke_1.849": -x**3*sqrt(x**2 - 4*x + 4*y(x)) - x**2*sqrt(x**2 - 4*x + 4*y(x)) + x/2 - sqrt(x**2 - 4*x + 4*y(x)) + Derivative(y(x), x) - 1,
        "kamke_1.850": -_F1(y(x) + log(cos(x) + 1) - log(sin(x))) + Derivative(y(x), x) - 1/sin(x),
        "kamke_1.851": Derivative(y(x), x) - (a**3*x**3 + 3*a**2*b*x**2*y(x) + a**2*b*x**2 + 3*a*b**2*x*y(x)**2 + 2*a*b**2*x*y(x) + b**3*y(x)**3 + b**3*y(x)**2 + b**3)/b**3,
        "kamke_1.852": Derivative(y(x), x) - (alpha**3*y(x)**3 + alpha**3*y(x)**2 + alpha**3 + 3*alpha**2*bbeta*x*y(x)**2 + 2*alpha**2*bbeta*x*y(x) + 3*alpha*bbeta**2*x**2*y(x) + alpha*bbeta**2*x**2 + bbeta**3*x**3)/alpha**3,
        "kamke_1.853": Derivative(y(x), x) - (x**3*y(x)**3 + 6*x**2*y(x)**2 + 14*x*y(x) + 2*x + 12)/(x**2*(x*y(x) + x + 2)),
        "kamke_1.854": Derivative(y(x), x) - (x**2*log(x)**2 + 2*x**2*log(x)*log(y(x)) + x**2*log(y(x))**2 + log(x) + log(y(x)) - 1)*y(x)/x,
        "kamke_1.855": Derivative(y(x), x) - (x**3*log(x)**2 + 2*x**3*log(x)*log(y(x)) + x**3*log(y(x))**2 + log(x) + log(y(x)) - 1)*y(x)/x,
        "kamke_1.856": -x*(_F1(-2*x + y(x)**2) + 1/x)/sqrt(y(x)**2) + Derivative(y(x), x),
        "kamke_1.857": -x**3*sqrt(x**2 - 2*x + 8*y(x) + 1) - x**2*sqrt(x**2 - 2*x + 8*y(x) + 1) + x/4 - sqrt(x**2 - 2*x + 8*y(x) + 1) + Derivative(y(x), x) - 1/4,
        "kamke_1.858": Derivative(y(x), x) - (a**3*y(x)**3 + a**3*y(x)**2 + a**3 + 3*a**2*b*x*y(x)**2 + 2*a**2*b*x*y(x) + 3*a*b**2*x**2*y(x) + a*b**2*x**2 + b**3*x**3)/a**3,
        "kamke_1.859": Derivative(y(x), x) - (x + _F1(-2*x + y(x)**2))/(x*sqrt(y(x)**2)),
        "kamke_1.860": Derivative(y(x), x) - (x**4*cos(2*y(x))/2 + x**4/2 + x**3*cos(2*y(x))/2 + x**3/2 + x*cos(2*y(x))/2 + x/2 - sin(2*y(x))/2)/x,
        "kamke_1.861": Derivative(y(x), x) - (_F1(y(x)*exp(1/x)) + y(x)*exp(1/x)/x)*exp(-1/x)/x,
        "kamke_1.862": -(_F1(x) - Ei(-log(y(x) - 1))/x)*log(y(x) - 1) + Derivative(y(x), x),
        "kamke_1.863": Derivative(y(x), x) - (x**4*sqrt(x**2 + y(x)**2) + x**3*sqrt(x**2 + y(x)**2) + x*sqrt(x**2 + y(x)**2) + y(x))/x,
        "kamke_1.864": Derivative(y(x), x) - (x*y(x)*exp(-x**2/2) + x*exp(-x**2/4) + 2*y(x)**2*exp(-3*x**2/4))*y(x)*exp(x**2/4)/(2*y(x)*exp(-x**2/4) + 2),
        "kamke_1.865": -(1 - y(x))*(-f(x) + y(x)*log(y(x) - 1)/(x*(1 - y(x))*log(x)) - log(y(x) - 1)/(x*(1 - y(x))*log(x))) + Derivative(y(x), x),
        "kamke_1.866": a/2 - x**3*sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) - x**2*sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) + x/2 - sqrt(a**2 + 2*a*x + x**2 + 4*y(x)) + Derivative(y(x), x),
        "kamke_1.867": -x**6/27 - x**4*y(x)/3 - x**4/9 - x**2*y(x)**2 - 2*x**2*y(x)/3 + 2*x/3 - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.868": x**6 - 3*x**4*y(x) - x**4 + 3*x**2*y(x)**2 + 2*x**2*y(x) - 2*x - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.869": Derivative(y(x), x) - (2*x**5 + 2*x**4 - 2*x**3*y(x) + x**3 - 2*x**2*y(x) + 3*x**2 - x - 2*y(x) + 1)/(x**2 - y(x)),
        "kamke_1.870": Derivative(y(x), x) - (x**4 + x**3 + x + x*exp(-y(x)/x) + y(x)*exp(-y(x)/x))*exp(y(x)/x)/x,
        "kamke_1.871": Derivative(y(x), x) - (2*x*y(x)**2 + 4*x*y(x)*log(2*x + 1) + 2*x*log(2*x + 1)**2 + y(x)**2 + 2*y(x)*log(2*x + 1) + log(2*x + 1)**2 - 2)/(2*x + 1),
        "kamke_1.872": Derivative(y(x), x) - (14*x**(7/2) - 5*sqrt(x)*y(x) - 5*sqrt(x) + 12*x**6/5 - 6*x**3*y(x) - 6*x**3 + 10*x - 5)/(x*(10*sqrt(x) + 2*x**3 - 5*y(x) - 5)),
        "kamke_1.873": Derivative(y(x), x) - (2*y(x) + 1)/(x*(2*x*y(x)**4 + 3*x*y(x)**3 + x*y(x)**2 + 2*x*y(x) + x - 2)),
        "kamke_1.874": -x*(a**3*x**12/512 + 3*a**2*x**8*y(x)/64 + a**2*x**8/64 + 3*a*x**4*y(x)**2/8 + a*x**4*y(x)/4 - a*x**2/2 + y(x)**3 + y(x)**2 + 1) + Derivative(y(x), x),
        "kamke_1.875": Derivative(y(x), x) - (-x**5*sqrt(x**2 + y(x)**2) + x**4*sqrt(x**2 + y(x)**2)*y(x) + x*y(x) + y(x))/(x*(x + 1)),
        "kamke_1.876": Derivative(y(x), x) + (x**2*y(x) - 2*x*y(x) - 2*x + y(x))*y(x)**2/(2*x*(x*y(x) - 2*y(x) - 2)),
        "kamke_1.877": Derivative(y(x), x) - (x**6 - 3*x**4*y(x) + 2*x**3 + 3*x**2*y(x)**2 - 2*x*y(x) - 2*x - y(x)**3)/(x**2 - y(x) - 1),
        "kamke_1.878": -(-64*a**3*x**3 + 48*a**2*x**2*y(x)**2 + 16*a**2*x**2 - 12*a*x*y(x)**4 - 8*a*x*y(x)**2 + y(x)**6 + y(x)**4 + 1)/y(x) + Derivative(y(x), x),
        "kamke_1.879": Derivative(y(x), x) - (-x**2*sqrt(x**2 + y(x)**2) + x*sqrt(x**2 + y(x)**2)*y(x) + x*y(x) + y(x))/(x*(x + 1)),
        "kamke_1.880": 2*a/(128*a**4*x**3 - 96*a**3*x**2*y(x)**2 - 32*a**3*x**2 + 24*a**2*x*y(x)**4 + 16*a**2*x*y(x)**2 - 2*a*y(x)**6 - 2*a*y(x)**4 - 2*a - y(x)) + Derivative(y(x), x),
        "kamke_1.881": Derivative(y(x), x) - (x**6 + 9*x**4*y(x) - 6*x**3 + 27*x**2*y(x)**2 - 18*x*y(x) - 18*x + 27*y(x)**3)/(9*x**2 + 27*y(x) + 27),
        "kamke_1.882": -sqrt(x)*(x**(3/2)/2 - x**9/216 + x**6*y(x)/12 + x**6/36 - x**3*y(x)**2/2 - x**3*y(x)/3 + y(x)**3 + y(x)**2 + 1) + Derivative(y(x), x),
        "kamke_1.883": Derivative(y(x), x) - x*(a**3*y(x)**6 + a**3*y(x)**4 + a**3 + 3*a**2*b*x**2*y(x)**4 + 2*a**2*b*x**2*y(x)**2 + 3*a*b**2*x**4*y(x)**2 + a*b**2*x**4 + b**3*x**6)/(a**(7/2)*y(x)),
        "kamke_1.884": -x*(-x**6 + 3*x**4*y(x)**2 + x**4 - 3*x**2*y(x)**4 - 2*x**2*y(x)**2 + y(x)**6 + y(x)**4 + 1)/y(x) + Derivative(y(x), x),
        "kamke_1.885": I*(x**6 + 12*x**4*y(x)**2 + 4*x**4 + 48*x**2*y(x)**4 + 32*x**2*y(x)**2 + 32*I*x + 64*y(x)**6 + 64*y(x)**4 + 64)/(128*y(x)) + Derivative(y(x), x),
        "kamke_1.886": Derivative(y(x), x) - (x**6*y(x)**3 - 3*x**5*y(x)**2 + x**4*y(x)**2 + 3*x**4*y(x) - 4*x**3*y(x) - x**3 + 2*x**2 + 1)/x**4,
        "kamke_1.887": Derivative(y(x), x) - (a**3*x**3*y(x)**3 + 3*a**2*x**2*y(x)**2 + a**2*x*y(x) + a**2*x + 3*a*x*y(x) + a + 1)/(a**2*x**2*(a*x*y(x) + a*x + 1)),
        "kamke_1.888": Derivative(y(x), x) - (x**4*y(x)**3 - 5*x**3*y(x)**2 + 6*x**2*y(x) - 2*x*y(x) - 2*x + 1)/(x**2*(x**2*y(x) - x + 1)),
        "kamke_1.889": -(y(x)**(9/2) + 27*y(x)**(3/2)*exp(2*x)/4 - 3*y(x)**(3/2)*exp(x) - 9*y(x)**3*exp(x)/2 + y(x)**3 - 27*exp(3*x)/8 + 9*exp(2*x)/4 + 1)*exp(x)/sqrt(y(x)) + Derivative(y(x), x),
        "kamke_1.890": -x/(x**6 + 3*x**4*y(x)**2 + x**4 + 3*x**2*y(x)**4 + 2*x**2*y(x)**2 + y(x)**6 + y(x)**4 - y(x) + 1) + Derivative(y(x), x),
        "kamke_1.891": Derivative(y(x), x) - (x**4*y(x) + 2*x**2*y(x) + 2*x**2 - 2*y(x))*y(x)**2/(x**3*(x**2*y(x) + x**2 - y(x))),
        "kamke_1.892": Derivative(y(x), x) - (x**2 + 2*x*y(x) + y(x)**2 + exp(-2/(x**2 - y(x)**2 - 1)))/(x**2 + 2*x*y(x) + y(x)**2 - exp(-2/(x**2 - y(x)**2 - 1))),
        "kamke_1.893": Derivative(y(x), x) - (x**3*y(x)**3 + x**3*y(x)**2 + x**3 + 6*x**2*y(x)**2 + 4*x**2*y(x) + 12*x*y(x) + 6*x + 8)/x**3,
        "kamke_1.894": I*(x**6 + 3*x**4*y(x)**2 + x**4 + 3*x**2*y(x)**4 + 2*x**2*y(x)**2 + I*x + y(x)**6 + y(x)**4 + 1)/y(x) + Derivative(y(x), x),
        "kamke_1.895": -x*(a**3*x**12 + 24*a**2*x**8*y(x) - 32*a**2*x**6 + 192*a*x**4*y(x)**2 - 256*a*x**2*y(x) - 256*a*x**2 + 512*y(x)**3)/(64*a*x**4 + 512*y(x) + 512) + Derivative(y(x), x),
        "kamke_1.896": -(-x**6 + 3*x**4*y(x)**2 + x**4 - 3*x**2*y(x)**4 - 2*x**2*y(x)**2 + x + y(x)**6 + y(x)**4 + 1)/y(x) + Derivative(y(x), x),
        "kamke_1.897": -sqrt(x)*(18*x**(9/2) - 108*x**(3/2)*y(x) - 108*x**(3/2) + x**9 - 18*x**6*y(x) + 108*x**3*y(x)**2 - 216*y(x)**3)/(36*x**3 - 216*y(x) - 216) + Derivative(y(x), x),
        "kamke_1.898": Derivative(y(x), x) - (64*x**6*y(x)**3 + 32*x**5*y(x) + 32*x**5 + 48*x**4*y(x)**2 + 8*x**3 + 12*x**2*y(x) + 1)/(16*x**6*(4*x**2*y(x) + 4*x**2 + 1)),
        "kamke_1.899": Derivative(y(x), x) - (x**6*y(x)**3 + x**6*y(x)**2 + x**6 + x**5/2 + 3*x**4*y(x)**2/4 + x**4*y(x)/2 + 3*x**2*y(x)/16 + x**2/16 + 1/64)/x**8,
        "kamke_1.900": -2*a*(4*a*x - y(x)**2 - 1)/(128*a**4*x**3 - 96*a**3*x**2*y(x)**2 + 24*a**2*x*y(x)**4 + 4*a*x*y(x) - 2*a*y(x)**6 - y(x)**3 - y(x)) + Derivative(y(x), x),
        "kamke_1.901": Derivative(y(x), x) - (-a*x*log(y(x)) + x**2 + y(x))*y(x)/(x*(a*x - y(x)*log(x) - y(x)*log(y(x)) - y(x))),
        "kamke_1.902": Derivative(y(x), x) - (x**6 - 3*x**4*y(x)**2 + x**3 + 3*x**2*y(x)**4 - x*y(x)**2 - x - y(x)**6)/((x**2 - y(x)**2 - 1)*y(x)),
        "kamke_1.903": Derivative(y(x), x) - (2*x**2*sin(y(x)/(2*x))*cos(y(x)/(2*x)) + y(x))*sin(y(x)/x)/(2*x*sin(y(x)/(2*x))*cos(y(x)/(2*x))),
        "kamke_1.904": Derivative(y(x), x) - (2*x**3*sin(y(x)/(2*x))*cos(y(x)/(2*x)) + y(x))*sin(y(x)/x)/(2*x*sin(y(x)/(2*x))*cos(y(x)/(2*x))),
        "kamke_1.905": Derivative(y(x), x) - (a**3*x**3*y(x)**3 + a**3*x**3*y(x)**2 + a**3*x**3 + 3*a**2*x**2*y(x)**2 + 2*a**2*x**2*y(x) + a**2*x + 3*a*x*y(x) + a*x + 1)/(a**3*x**3),
        "kamke_1.906": -x*(x**2 + y(x)**2 + 1)/(x**6 + 3*x**4*y(x)**2 + 3*x**2*y(x)**4 - x**2*y(x) + y(x)**6 - y(x)**3 - y(x)) + Derivative(y(x), x),
        "kamke_1.907": Derivative(y(x), x) - (x**2*sin(x) - 2*x**2*cos(x) + x**2*cos(2*x)/2 + 3*x**2/2 + 2*x*y(x)*cos(x) - 2*x*y(x) - x*cos(x) + x + y(x)**2)/x,
        "kamke_1.908": -4*x*(a - 1)*(a + 1)/(a**6*x**4 - 3*a**4*x**4 - 2*a**4*x**2*y(x)**2 + 3*a**2*x**4 + 4*a**2*x**2*y(x)**2 + a**2*y(x)**4 - x**4 - 2*x**2*y(x)**2 - y(x)**4 + 4*y(x)) + Derivative(y(x), x),
        "kamke_1.909": Derivative(y(x), x) - (x**3*y(x)**6 + x**3*y(x)**4 + x**3 + 3*x**2*y(x)**4 + 2*x**2*y(x)**2 + 3*x*y(x)**2 + x + 1)/(x**5*y(x)),
        "kamke_1.910": Derivative(y(x), x) - (x**6 + 3*x**5*y(x) + 3*x**4*y(x)**2 + x**4 + x**3*y(x)**3 + 2*x**3*y(x) + x**2*y(x)**2 - 2*x - y(x) + 1)/x,
        "kamke_1.911": -(_F1(x) - log(y(x))*cos(x)/sin(x) + log(y(x))/x)*y(x) + Derivative(y(x), x),
        "kamke_1.912": -2*a*x/(-128*a**4 + 96*a**3*x*y(x)**2 + 32*a**3*x - 24*a**2*x**2*y(x)**4 - 16*a**2*x**2*y(x)**2 + 2*a*x**3*y(x)**6 + 2*a*x**3*y(x)**4 + 2*a*x**3 - x**3*y(x)) + Derivative(y(x), x),
        "kamke_1.913": Derivative(y(x), x) - (-y(x)**3*log(x)**3 + y(x)**3*log(x)**2 + y(x)**3 + 3*y(x)**2*log(x)**2 - 2*y(x)**2*log(x) - 3*y(x)*log(x) + y(x) + 1)/(x*y(x)),
        "kamke_1.914": -2*a*(-4*a + x*y(x)**2 + x)/(-128*a**4 + 96*a**3*x*y(x)**2 - 24*a**2*x**2*y(x)**4 + 2*a*x**3*y(x)**6 + 4*a*x**2*y(x) - x**3*y(x)**3 - x**3*y(x)) + Derivative(y(x), x),
        "kamke_1.915": Derivative(y(x), x) - (-8*y(x)**3*log(x)**3 + 4*y(x)**3*log(x)**2 + y(x)**3 + 12*y(x)**2*log(x)**2 - 4*y(x)**2*log(x) - 6*y(x)*log(x) + y(x) + 1)/(x*y(x)),
        "kamke_1.916": Derivative(y(x), x) - (x**4*log(x)**2 + 2*x**4*log(x)*log(y(x)) + x**4*log(y(x))**2 + x*log(x) + x*log(y(x)) - x + log(x) + log(y(x)) - 1)*y(x)/(x*(x + 1)),
        "kamke_1.917": Derivative(y(x), x) - (x*log(x)**2 + 2*x*log(x)*log(y(x)) + x*log(x) + x*log(y(x))**2 + x*log(y(x)) - x + log(x) + log(y(x)) - 1)*y(x)/(x*(x + 1)),
        "kamke_1.918": Derivative(y(x), x) - 2*y(x)**8/(128*x**3*y(x)**6 + 32*x**2*y(x)**6 + 96*x**2*y(x)**4 + 16*x*y(x)**4 + 24*x*y(x)**2 + 2*y(x)**6 + y(x)**5 + 2*y(x)**2 + 2),
        "kamke_1.919": -(x + sqrt(y(x)) - y(x))*y(x)**(3/2)/(x**3 - 3*x**2*y(x) + x*y(x)**(3/2) + 3*x*y(x)**2 - y(x)**(5/2) - y(x)**3 + y(x)**2) + Derivative(y(x), x),
        "kamke_1.920": -2*(4*x*y(x)**2 + y(x)**2 + 1)*y(x)**6/(128*x**3*y(x)**6 + 96*x**2*y(x)**4 + 4*x*y(x)**5 + 24*x*y(x)**2 + y(x)**5 + y(x)**3 + 2) + Derivative(y(x), x),
        "kamke_1.921": -(_F1(x) + log(y(x))/x - log(y(x))/(x*log(x)))*y(x) + Derivative(y(x), x),
        "kamke_1.922": Derivative(y(x), x) - y(x)**2/(x**3 + x**2*sqrt(y(x)) - 3*x**2*y(x) - 2*x*y(x)**(3/2) + 3*x*y(x)**2 + y(x)**(5/2) + y(x)**(3/2) - y(x)**3 + y(x)**2),
        "kamke_1.923": Derivative(y(x), x) - (x**2 + 2*x*y(x) + y(x)**2 + exp((-2*x + 2*y(x))*(x + y(x))))/(x**2 + 2*x*y(x) + y(x)**2 - exp((-2*x + 2*y(x))*(x + y(x)))),
        "kamke_1.924": -(_F1(x) + log(y(x))**2/(2*x))*y(x)/log(y(x)) + Derivative(y(x), x),
        "kamke_1.925": Derivative(y(x), x) - (x**2 + 2*x*y(x) + y(x)**2 + exp(2*(x - y(x))**2*(x + y(x))**2))/(x**2 + 2*x*y(x) + y(x)**2 - exp(2*(x - y(x))**2*(x + y(x))**2)),
        "kamke_1.926": Derivative(y(x), x) - (x**3*y(x)**3/16 - x**2*y(x)**3/2 - 3*x**2*y(x)**2/8 + x*y(x)**3 + x*y(x)**2 + 3*x*y(x)/4 - 1/2)/(x*(x*y(x) - 2*y(x) - 2)),
        "kamke_1.927": -x*(-x**6*exp(-3*x**2)/8 + 3*x**4*y(x)*exp(-2*x**2)/4 + x**4*exp(-2*x**2)/4 - 3*x**2*y(x)**2*exp(-x**2)/2 - x**2*y(x)*exp(-x**2) - x**2*exp(-x**2) + y(x)**3 + y(x)**2 + 1 + exp(-x**2)) + Derivative(y(x), x),
        "kamke_1.928": Derivative(y(x), x) - (x**2*exp(-y(x)/x) + x*y(x)*exp(-y(x)/x) + x + x*exp(-y(x)/x) + y(x)*exp(-y(x)/x))*exp(y(x)/x)/(x*(x + 1)),
        "kamke_1.929": Derivative(y(x), x) + (x**3*y(x)**3 - 2*x**2*y(x)**3 - 6*x**2*y(x)**2 + 16*x*y(x)**3 + 8*x*y(x)**2 + 12*x*y(x) - 8*y(x)**3 - 8*y(x) - 8)/(32*x*y(x)),
        "kamke_1.930": Derivative(y(x), x) - (x**4 + x**2*exp(-y(x)/x) + x*y(x)*exp(-y(x)/x) + x*exp(-y(x)/x) + y(x)*exp(-y(x)/x))*exp(y(x)/x)/(x*(x + 1)),
        "kamke_1.931": Derivative(y(x), x) - (x**6 + 3*x**5*y(x) + 3*x**4*y(x)**2 + x**3*y(x)**3 - 2*x**3 - 3*x**2*y(x) - x*y(x)**2 - 2*x - y(x))/(x*(x**2 + x*y(x) + 1)),
        "kamke_1.932": -x*(y(x)**3*exp(9*x**2/2)/243 + y(x)**3*exp(3*x**2)/81 + y(x)**3/9 + y(x)**2*exp(9*x**2/2)/27 + 2*y(x)**2*exp(3*x**2)/27 + y(x)*exp(9*x**2/2)/9 + y(x)*exp(3*x**2)/9 + exp(9*x**2/2)/9)*exp(-3*x**2/2)/y(x) + Derivative(y(x), x),
        "kamke_1.933": Derivative(y(x), x) - (-x**3*log(x)**3 + x**3*log(x)**2 + x**3 + 3*x**2*y(x)*log(x)**2 - 2*x**2*y(x)*log(x) + x**2 - 3*x*y(x)**2*log(x) + x*y(x)**2 + x*y(x) + y(x)**3)/x**2,
        "kamke_1.934": x**6/64 + 3*x**5/32 - 3*x**4*y(x)/16 + x**4/8 - 3*x**3*y(x)/4 - x**3/8 + 3*x**2*y(x)**2/4 - x**2*y(x)/4 - x**2/4 + 3*x*y(x)**2/2 + x*y(x) - x/2 - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.935": -x**6/64 + 3*x**5/16 - 3*x**4*y(x)/16 - 13*x**4/16 + 3*x**3*y(x)/2 + 3*x**3/2 - 3*x**2*y(x)**2/4 - 7*x**2*y(x)/2 - x**2 + 3*x*y(x)**2 + 2*x*y(x) + x/2 - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.936": -x**6/512 + 3*x**5/256 - 3*x**4*y(x)/64 - 5*x**4/128 + 3*x**3*y(x)/16 + 5*x**3/64 - 3*x**2*y(x)**2/8 - 7*x**2*y(x)/16 - x**2/16 + 3*x*y(x)**2/4 + x*y(x)/2 + x/4 - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.937": Derivative(y(x), x) - (2*x*y(x)**3 + 6*x*y(x)**2*log(2*x + 1) + 6*x*y(x)*log(2*x + 1)**2 + 2*x*log(2*x + 1)**3 + y(x)**3 + 3*y(x)**2*log(2*x + 1) + 3*y(x)*log(2*x + 1)**2 - 2*y(x) + log(2*x + 1)**3 - 2*log(2*x + 1) - 2)/((2*x + 1)*(y(x) + log(2*x + 1) + 1)),
        "kamke_1.938": Derivative(y(x), x) - (x**6 - 3*x**5 + 3*x**4*y(x) + 4*x**4 - 6*x**3*y(x) - 3*x**3 + 3*x**2*y(x)**2 + 5*x**2*y(x) - x**2 - 3*x*y(x)**2 - 2*x*y(x) + x + y(x)**3 + y(x)**2 + 1)/x,
        "kamke_1.939": Derivative(y(x), x) - (x**6 + 6*x**5 - 12*x**4*y(x) + 12*x**4 - 48*x**3*y(x) + 16*x**3 + 48*x**2*y(x)**2 - 48*x**2*y(x) + 16*x**2 + 96*x*y(x)**2 - 32*x*y(x) - 32*x - 64*y(x)**3)/(16*x**2 + 32*x - 64*y(x) - 64),
        "kamke_1.940": Derivative(y(x), x) - (x**3*log(x)**3 - 3*x**2*y(x)*log(x)**2 + x**2*log(x) - x**2 + 3*x*y(x)**2*log(x) + x*y(x)*log(x) - 2*x*y(x) - y(x)**3 - y(x)**2)/(x*(x*log(x) - x - y(x))),
        "kamke_1.941": Derivative(y(x), x) - (x**6 - 12*x**5 + 12*x**4*y(x) + 48*x**4 - 96*x**3*y(x) - 72*x**3 + 48*x**2*y(x)**2 + 192*x**2*y(x) + 32*x**2 - 192*x*y(x)**2 - 32*x*y(x) - 32*x + 64*y(x)**3)/(16*x**2 - 64*x + 64*y(x) + 64),
        "kamke_1.942": -(-x**2 - 2*x*y(x) - y(x)**2 - exp(2*(x - y(x))**3*(x + y(x))**3/(x**2 - y(x)**2 - 1)))/(-x**2 - 2*x*y(x) - y(x)**2 + exp(2*(x - y(x))**3*(x + y(x))**3/(x**2 - y(x)**2 - 1))) + Derivative(y(x), x),
        "kamke_1.943": Derivative(y(x), x) - (x**6 - 6*x**5 + 24*x**4*y(x) + 12*x**4 - 96*x**3*y(x) - 24*x**3 + 192*x**2*y(x)**2 + 96*x**2*y(x) + 32*x**2 - 384*x*y(x)**2 - 128*x*y(x) - 128*x + 512*y(x)**3)/(64*x**2 - 128*x + 512*y(x) + 512),
        "kamke_1.944": Derivative(y(x), x) - (a**3*x**6 + 6*a**2*b*x**5 + 12*a**2*x**4*y(x) - 8*a**2*x**3 + 12*a*b**2*x**4 + 48*a*b*x**3*y(x) - 16*a*b*x**2 + 48*a*x**2*y(x)**2 - 32*a*x*y(x) - 32*a*x + 8*b**3*x**3 + 48*b**2*x**2*y(x) + 96*b*x*y(x)**2 + 64*y(x)**3)/(16*a*x**2 + 32*b*x + 64*y(x) + 64),
        "kamke_1.945": Derivative(y(x), x) - (8*a**3*x**3 + 12*a**2*x**4 + 48*a**2*x**2*y(x) + 6*a*x**5 + 48*a*x**3*y(x) - 16*a*x**2 + 96*a*x*y(x)**2 + x**6 + 12*x**4*y(x) - 8*x**3 + 48*x**2*y(x)**2 - 32*x*y(x) - 32*x + 64*y(x)**3)/(32*a*x + 16*x**2 + 64*y(x) + 64),
        "kamke_1.946": -x*(x**6*exp(-3*x**2) - 6*x**4*y(x)*exp(-2*x**2) - 4*x**4*exp(-2*x**2) + 12*x**2*y(x)**2*exp(-x**2) + 8*x**2*y(x)*exp(-x**2) + 8*x**2*exp(-x**2) + 4*x**2*exp(-2*x**2) - 8*y(x)**3 - 8*y(x)*exp(-x**2) - 8*exp(-x**2))/(4*x**2*exp(-x**2) - 8*y(x) - 8) + Derivative(y(x), x),
        "kamke_1.947": Derivative(y(x), x) - (x**3*sin(x) + x**2*y(x)**2 + 2*x**2*y(x)*cos(x) + x**2*cos(x) + x**2*cos(2*x)/2 + x**2/2 - 2*x*y(x)*sin(x) + 2*x*y(x) - x*sin(x) - x*sin(2*x) + 2*x*cos(x) + x - 2*sin(x) - cos(2*x)/2 + 3/2)/x**3,
        "kamke_1.948": Derivative(y(x), x) + 216*y(x)/(36*x**2 - 24*x*y(x)**4 - 36*x*y(x)**3 - 72*x*y(x)**2 - 72*x*y(x) + 4*y(x)**8 + 12*y(x)**7 + 33*y(x)**6 + 60*y(x)**5 - 216*y(x)**4 - 252*y(x)**3 - 396*y(x)**2 - 216*y(x)),
        "kamke_1.949": Derivative(y(x), x) - (x**6 - 3*x**5 + 3*x**4*y(x) + x**4 - 6*x**3*y(x) + 2*x**3 + 3*x**2*y(x)**2 + x**2*y(x) - 3*x**2 - 3*x*y(x)**2 + x*y(x) + x + y(x)**3)/(x*(x**2 - x + y(x) + 1)),
        "kamke_1.950": -a**3*x**6/64 - 3*a**2*b*x**5/32 - 3*a**2*x**4*y(x)/16 - a**2*x**4/16 - 3*a*b**2*x**4/16 - 3*a*b*x**3*y(x)/4 - a*b*x**3/4 - 3*a*x**2*y(x)**2/4 - a*x**2*y(x)/2 + a*x/2 - b**3*x**3/8 - 3*b**2*x**2*y(x)/4 - b**2*x**2/4 - 3*b*x*y(x)**2/2 - b*x*y(x) - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.951": -a**3*x**3/8 - 3*a**2*x**4/16 - 3*a**2*x**2*y(x)/4 - a**2*x**2/4 - 3*a*x**5/32 - 3*a*x**3*y(x)/4 - a*x**3/4 - 3*a*x*y(x)**2/2 - a*x*y(x) - x**6/64 - 3*x**4*y(x)/16 - x**4/16 - 3*x**2*y(x)**2/4 - x**2*y(x)/2 + x/2 - y(x)**3 - y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.952": Derivative(y(x), x) - (-x**5*sqrt(x**2 + y(x)**2) + x**4*sqrt(x**2 + y(x)**2)*y(x) - x**4*sqrt(x**2 + y(x)**2) + x**3*sqrt(x**2 + y(x)**2)*y(x) - x**2*sqrt(x**2 + y(x)**2) + x*sqrt(x**2 + y(x)**2)*y(x) + y(x))/x,
        "kamke_1.953": Derivative(y(x), x) - (x**4*log(x)**2 + 2*x**4*log(x)*log(y(x)) + x**4*log(y(x))**2 + x**3*log(x)**2 + 2*x**3*log(x)*log(y(x)) + x**3*log(y(x))**2 + x*log(x)**2 + 2*x*log(x)*log(y(x)) + x*log(y(x))**2 + log(x) + log(y(x)) - 1)*y(x)/x,
        "kamke_1.954": Derivative(y(x), x) - (-24*x**(13/2)/25 + 24*x**(7/2)*y(x)/5 + 8*x**(7/2)/5 - 8*x**(3/2) - 6*sqrt(x)*y(x)**2 - 4*sqrt(x)*y(x) + sqrt(x) - 8*x**9/125 + 12*x**6*y(x)/25 + 4*x**6/25 - 24*x**4/5 - 6*x**3*y(x)**2/5 - 4*x**3*y(x)/5 + 6*x**3/5 + 12*x*y(x) + 4*x + y(x)**3 + y(x)**2 + 1)/x,
        "kamke_1.955": Derivative(y(x), x) - (24*x**(13/2)/5 - 24*x**(7/2)*y(x) + 14*x**(7/2) + 40*x**(3/2) + 30*sqrt(x)*y(x)**2 - 5*sqrt(x)*y(x) - 5*sqrt(x) + 8*x**9/25 - 12*x**6*y(x)/5 + 12*x**6/5 + 24*x**4 + 6*x**3*y(x)**2 - 6*x**3*y(x) - 6*x**3 - 60*x*y(x) + 10*x - 5*y(x)**3)/(x*(10*sqrt(x) + 2*x**3 - 5*y(x) - 5)),
        "kamke_1.956": Derivative(y(x), x) - (x**2*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1))*log(x)**2 + 2*x**2*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1))*log(x) + x**2*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1)) - x**2*x**(2/(log(x) + 1))*exp(2*log(x)**2/(log(x) + 1))*log(x) - x**2*x**(2/(log(x) + 1))*exp(2*log(x)**2/(log(x) + 1)) - 1)*y(x)/(x*(log(x) + 1)),
        "kamke_1.957": Derivative(y(x), x) - (x**3*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1))*log(x)**2 + 2*x**3*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1))*log(x) + x**3*x**(2/(log(x) + 1))*y(x)*exp(2*log(x)**2/(log(x) + 1)) - x**3*x**(2/(log(x) + 1))*exp(2*log(x)**2/(log(x) + 1))*log(x) - x**3*x**(2/(log(x) + 1))*exp(2*log(x)**2/(log(x) + 1)) - 1)*y(x)/(x*(log(x) + 1)),
        "kamke_1.958": Derivative(y(x), x) - (2*x*y(x)**3 + 6*x*y(x)**2*log(2*x + 1) + 2*x*y(x)**2 + 6*x*y(x)*log(2*x + 1)**2 + 4*x*y(x)*log(2*x + 1) + 2*x*log(2*x + 1)**3 + 2*x*log(2*x + 1)**2 + 2*x + y(x)**3 + 3*y(x)**2*log(2*x + 1) + y(x)**2 + 3*y(x)*log(2*x + 1)**2 + 2*y(x)*log(2*x + 1) + log(2*x + 1)**3 + log(2*x + 1)**2 - 1)/(2*x + 1),
        "kamke_1.959": Derivative(y(x), x) - (x**3*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - y(x)*sin(y(x)/x)/2 + y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2)/(x*sin(y(x)/(2*x))*cos(y(x)/(2*x))*cos(y(x)/x)),
        "kamke_1.960": Derivative(y(x), x) - (x**2*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - y(x)*sin(y(x)/x)/2 + y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2)/(x*sin(y(x)/(2*x))*cos(y(x)/(2*x))*cos(y(x)/x)),
        "kamke_1.961": Derivative(y(x), x) - (x**2 + 2*x*y(x) + y(x)**2 + exp(-2*x**6 + 6*x**4*y(x)**2 + 2*x**4 - 6*x**2*y(x)**4 - 4*x**2*y(x)**2 + 2*y(x)**6 + 2*y(x)**4 + 2))/(x**2 + 2*x*y(x) + y(x)**2 - exp(-2*x**6 + 6*x**4*y(x)**2 + 2*x**4 - 6*x**2*y(x)**4 - 4*x**2*y(x)**2 + 2*y(x)**6 + 2*y(x)**4 + 2)),
        "kamke_1.962": -4*x*(a - 1)*(a + 1)*(a**2*x**2 - x**2 - y(x)**2 - 2)/(a**8*x**6 - 4*a**6*x**6 - 3*a**6*x**4*y(x)**2 + 6*a**4*x**6 + 9*a**4*x**4*y(x)**2 + 3*a**4*x**2*y(x)**4 - 4*a**2*x**6 - 9*a**2*x**4*y(x)**2 - 6*a**2*x**2*y(x)**4 + 4*a**2*x**2*y(x) - a**2*y(x)**6 + x**6 + 3*x**4*y(x)**2 + 3*x**2*y(x)**4 - 4*x**2*y(x) + y(x)**6 - 4*y(x)**3 - 8*y(x)) + Derivative(y(x), x),
        "kamke_1.963": Derivative(y(x), x) - (15*x**3*cos(x)/4 - 3*x**3*cos(2*x)/2 + x**3*cos(3*x)/4 - 5*x**3/2 - 6*x**2*y(x)*cos(x) + 3*x**2*y(x)*cos(2*x)/2 + 9*x**2*y(x)/2 + x**2*sin(x) - 2*x**2*cos(x) + x**2*cos(2*x)/2 + 3*x**2/2 + 3*x*y(x)**2*cos(x) - 3*x*y(x)**2 + 2*x*y(x)*cos(x) - 2*x*y(x) - x*cos(x) + x + y(x)**3 + y(x)**2 + 1)/x,
        "kamke_1.964": 8*x*(a - 1)*(a + 1)/(a**8*x**6 - 4*a**6*x**6 - 3*a**6*x**4*y(x)**2 - 2*a**6*x**4 + 6*a**4*x**6 + 9*a**4*x**4*y(x)**2 + 6*a**4*x**4 + 3*a**4*x**2*y(x)**4 + 4*a**4*x**2*y(x)**2 - 4*a**2*x**6 - 9*a**2*x**4*y(x)**2 - 6*a**2*x**4 - 6*a**2*x**2*y(x)**4 - 8*a**2*x**2*y(x)**2 - a**2*y(x)**6 - 2*a**2*y(x)**4 - 8*a**2 + x**6 + 3*x**4*y(x)**2 + 2*x**4 + 3*x**2*y(x)**4 + 4*x**2*y(x)**2 + y(x)**6 + 2*y(x)**4 - 8*y(x) + 8) + Derivative(y(x), x),
        "kamke_1.965": Derivative(y(x), x) - (x**4*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + x**3*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + x*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - y(x)*sin(y(x)/x)/2 + y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2)/(x*sin(y(x)/(2*x))*cos(y(x)/(2*x))*cos(y(x)/x)),
        "kamke_1.966": Derivative(y(x), x) + 1296*y(x)/(216*x**3 - 216*x**2*y(x)**4 - 324*x**2*y(x)**3 - 648*x**2*y(x)**2 - 648*x**2*y(x) + 216*x**2 + 72*x*y(x)**8 + 216*x*y(x)**7 + 594*x*y(x)**6 + 1080*x*y(x)**5 + 1152*x*y(x)**4 + 1080*x*y(x)**3 + 216*x*y(x)**2 - 432*x*y(x) - 8*y(x)**12 - 36*y(x)**11 - 126*y(x)**10 - 315*y(x)**9 - 570*y(x)**8 - 846*y(x)**7 - 882*y(x)**6 - 612*y(x)**5 - 1944*y(x)**4 - 1728*y(x)**3 - 2376*y(x)**2 - 1296*y(x) + 216),
        "kamke_1.967": x*(64*x**9 - 288*x**8*y(x) - 96*x**8 + 432*x**7*y(x)**2 + 288*x**7*y(x) - 144*x**7 - 216*x**6*y(x)**3 - 216*x**6*y(x)**2 - 288*x**6*y(x) - 456*x**6 + 864*x**5*y(x)**2 + 1008*x**5*y(x) - 576*x**5 - 648*x**4*y(x)**3 - 972*x**4*y(x)**2 - 216*x**4*y(x) - 864*x**4 + 432*x**3*y(x)**2 + 720*x**3*y(x) - 756*x**3 - 648*x**2*y(x)**3 - 1296*x**2*y(x)**2 - 594*x**2*y(x) - 1134*x**2 - 432*x - 216*y(x)**3 - 540*y(x)**2 - 378*y(x) - 513)/(216*(x**2 + 1)**4) + Derivative(y(x), x),
        "kamke_1.968": Derivative(y(x), x) - (x**4*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + x*y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - x*y(x)*sin(y(x)/x)/2 + x*y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2 + y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - y(x)*sin(y(x)/x)/2 + y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2)/(x*(x + 1)*sin(y(x)/(2*x))*cos(y(x)/(2*x))*cos(y(x)/x)),
        "kamke_1.969": Derivative(y(x), x) - (x*y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - x*y(x)*sin(y(x)/x)/2 + x*y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2 + x*sin(y(x)/(2*x))*sin(y(x)/x)*cos(y(x)/(2*x)) + y(x)*sin(y(x)/(2*x))*cos(y(x)/(2*x))/2 - y(x)*sin(y(x)/x)/2 + y(x)*sin(3*y(x)/(2*x))*cos(y(x)/(2*x))/2)/(x*(x + 1)*sin(y(x)/(2*x))*cos(y(x)/(2*x))*cos(y(x)/x)),
        "kamke_1.970": 216*(6*x - 2*y(x)**4 - 3*y(x)**3 - 6*y(x)**2 - 6*y(x) + 6)*y(x)/(216*x**3 - 216*x**2*y(x)**4 - 324*x**2*y(x)**3 - 648*x**2*y(x)**2 - 648*x**2*y(x) + 72*x*y(x)**8 + 216*x*y(x)**7 + 594*x*y(x)**6 + 1080*x*y(x)**5 - 432*x*y(x)**4 - 648*x*y(x)**3 - 1944*x*y(x)**2 - 1296*x*y(x) - 8*y(x)**12 - 36*y(x)**11 - 126*y(x)**10 - 315*y(x)**9 - 18*y(x)**8 + 594*y(x)**7 + 2484*y(x)**6 + 4428*y(x)**5 + 2808*y(x)**4 + 1728*y(x)**3 - 1296*y(x)**2 - 1296*y(x)) + Derivative(y(x), x),
        "kamke_1.971": Derivative(y(x), x) - (x*y(x) + 1)**3/x**5,
        "kamke_1.972": -x*(-2*x**4 + 2*x**2*y(x) - x**2 + 1)/(-x**2 + y(x)) + Derivative(y(x), x),
        "kamke_1.973": -(y(x)**2 + y(x)*exp(b*x) + exp(2*b*x))*y(x)*exp(-2*b*x) + Derivative(y(x), x),
        "kamke_1.974": x**6 - 3*x**4*y(x) + 3*x**2*y(x)**2 - 2*x - y(x)**3 + Derivative(y(x), x),
        "kamke_1.975": -x**6/27 - x**4*y(x)/3 - x**2*y(x)**2 + 2*x/3 - y(x)**3 + Derivative(y(x), x),
        "kamke_1.976": Derivative(y(x), x) - (x**7*y(x)**2 + x**4*y(x) + x - 3)*y(x)/x,
        "kamke_1.977": -x*(y(x)**2 + y(x)*exp(-x**2) + exp(-2*x**2))*y(x)*exp(2*x**2) + Derivative(y(x), x),
        "kamke_1.978": Derivative(y(x), x) - (x**2 + x*y(x) + x + y(x)**2)*y(x)/x**2,
        "kamke_1.979": Derivative(y(x), x) - (-x**3 + 3*x**2*y(x) - 3*x*y(x)**2 + x + y(x)**3)/x,
        "kamke_1.980": Derivative(y(x), x) - (x**3*y(x)**3 + 6*x**2*y(x)**2 + 12*x*y(x) + 2*x + 8)/x**3,
        "kamke_1.981": Derivative(y(x), x) - (a**3*x**3*y(x)**3 + 3*a**2*x**2*y(x)**2 + a**2*x + 3*a*x*y(x) + 1)/(a**3*x**3),
        "kamke_1.982": -(x*exp(x**2/2) + 2*y(x)**2 + 2*y(x)*exp(x**2/4) + 2*exp(x**2/2))*y(x)*exp(-x**2/2)/2 + Derivative(y(x), x),
        "kamke_1.983": Derivative(y(x), x) - (-x**3 + 3*x**2*y(x) + x**2 - 3*x*y(x)**2 + y(x)**3)/((x - 1)*(x + 1)),
        "kamke_1.984": Derivative(y(x), x) - (x - 1)*(x**2*y(x)**2 + x*y(x)*exp(x) + exp(2*x))*y(x)*exp(-2*x)/x,
        "kamke_1.985": Derivative(y(x), x) - (x*y(x) + 1)*(x**2*y(x)**2 + x**2*y(x) + x**2 + 2*x*y(x) + x + 1)/x**5,
        "kamke_1.986": Derivative(y(x), x) - (-x**3*log(x)**3 + 3*x**2*y(x)*log(x)**2 + x**2 - 3*x*y(x)**2*log(x) + x*y(x) + y(x)**3)/x**2,
        "kamke_1.987": (-a*x**2 + y(x)**2)*F(x) + Derivative(y(x), x) - y(x)/x,
        "kamke_1.988": (-x**2 - 2*x*y(x) + y(x)**2)*F(x) + Derivative(y(x), x) - y(x)/x,
        "kamke_1.989": (-a*y(x)**2 - b*x**2)*F(x) + Derivative(y(x), x) - y(x)/x,
        "kamke_1.990": -2*x + (-x**4 + 2*x**2*y(x) - y(x)**2 + 1)*F(x) + Derivative(y(x), x),
        "kamke_1.991": (x**2 + 2*x*y(x) - y(x)**2)*F(x) + Derivative(y(x), x) - y(x)/x,
        "kamke_1.992": (-x**3 - 7*x*y(x)**2)*F(x) + Derivative(y(x), x) - y(x)/x,
        "kamke_1.993": (-y(x)**2 - 2*y(x)*log(x) - log(x)**2)*F(x) + Derivative(y(x), x) - y(x)/(x*log(x)),
        "kamke_1.994": x**3*(-y(x)**2 - 2*y(x)*log(x) - log(x)**2) + Derivative(y(x), x) - y(x)/(x*log(x)),
        "kamke_1.995": -(y(x) - exp(x))**2 - exp(x) + Derivative(y(x), x),
        "kamke_1.996": Derivative(y(x), x) - ((y(x) - Si(x))**2 + sin(x))/x,
        "kamke_1.997": -(y(x) + cos(x))**2 - sin(x) + Derivative(y(x), x),
        "kamke_1.998": Derivative(y(x), x) - ((y(x) - log(x) - Ci(x))**2 + cos(x))/x,
        "kamke_1.999": Derivative(y(x), x) - (x + (-x + y(x) + log(x + 1))**2)/(x + 1),
        "kamke_1.1000": Derivative(y(x), x) - (x**3 + 2*x**2*y(x) + x*y(x)*log(x) - x*y(x) - y(x)**2)/(x**2*(x + log(x)))
    }

    chapter_2 = {
        "kamke_2.1": Derivative(y(x), (x, 2)),
        "kamke_2.2": y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.3": y(x) - sin(n*x) + Derivative(y(x), (x, 2)),
        "kamke_2.4": -a*cos(b*x) + y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.5": y(x) - sin(a*x)*sin(b*x) + Derivative(y(x), (x, 2)),
        "kamke_2.6": -y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.7": -4*x**2*exp(x**2) - 2*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.8": a**2*y(x) - cot(a*x) + Derivative(y(x), (x, 2)),
        "kamke_2.9": l*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.10": (a*x + b)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.11": -(x**2 + 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.12": -(a + x**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.13": -(a**2*x**2 + a)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.14": -c*x**a*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.15": -(a**2*x**(2*n) - 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.16": (a*x**(2*c) + b*x**(c - 1))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.17": (-v**2 + exp(2*x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.18": a*y(x)*exp(b*x) + Derivative(y(x), (x, 2)),
        "kamke_2.19": -(4*a**2*b**2*x**2*exp(2*b*x**2) - 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.20": (a*exp(2*x) + b*exp(x) + c)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.21": (a*cosh(x)**2 + b)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.22": (a*cos(2*x) + b)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.23": (a*cos(x)**2 + b)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.24": -(2*tan(x)**2 + 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.25": -(a + m*(m - 1)/cos(x)**2 + n*(n - 1)/sin(x)**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.26": -(B + n*(n + 1)*WeierstrassP(x, g2, g3))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.27": -(b + k**2*n*(n + 1)*JacobiSN(x, k)**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.28": -(a*p(x) + b + 7*Derivative(p(x), (x, 2))/3 + Derivative(p(x), (x, 4))/30)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.29": -(f(x)**2 + Derivative(f(x), x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.30": (l + P(x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.31": -f(x)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.32": ((1/4 - v**2)*Derivative(g(x), x)**2/g(x) + Derivative(g(x), x)**2 + Derivative(g(x), (x, 3))/(2*Derivative(g(x), x)) - 3*Derivative(g(x), (x, 2))**2/(4*Derivative(g(x), x)**2))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.33": a*y(x)*exp(-2*x) + Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.34": y(x)*exp(2*x) - Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.35": a*Derivative(y(x), x) + b*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.36": a*Derivative(y(x), x) + b*y(x) - f(x) + Derivative(y(x), (x, 2)),
        "kamke_2.37": a*Derivative(y(x), x) - (b**2*x**2 + c)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.38": 2*a*Derivative(y(x), x) + f(x)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.39": x*Derivative(y(x), x) + y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.40": x*Derivative(y(x), x) - y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.41": x*Derivative(y(x), x) + (n + 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.42": -n*y(x) + x*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.43": -x*Derivative(y(x), x) + 2*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.44": -a*y(x) - x*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.45": -x*Derivative(y(x), x) + (x - 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.46": a*y(x) - 2*x*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.47": 4*x*Derivative(y(x), x) + (4*x**2 + 2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.48": -4*x*Derivative(y(x), x) + (2*n + 3*x**2 - 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.49": -4*x*Derivative(y(x), x) + (4*x**2 - 1)*y(x) - exp(x) + Derivative(y(x), (x, 2)),
        "kamke_2.50": -4*x*Derivative(y(x), x) + (4*x**2 - 2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.51": -4*x*Derivative(y(x), x) + (4*x**2 - 3)*y(x) - exp(x**2) + Derivative(y(x), (x, 2)),
        "kamke_2.52": a*x*Derivative(y(x), x) + b*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.53": a**2*x**2*y(x) + 2*a*x*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.54": (a*x + b)*Derivative(y(x), x) + (c*x + d)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.55": (a*x + b)*Derivative(y(x), x) + (a1*x**2 + b1*x + c1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.56": -x**2*Derivative(y(x), x) + x*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.57": -x**2*Derivative(y(x), x) - (x + 1)**2*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.58": -x**2*(x + 1)*Derivative(y(x), x) + x*(x**4 - 2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.59": x**4*Derivative(y(x), x) - x**3*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.60": a*x**(qs - 1)*Derivative(y(x), x) + b*x**(qs - 2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.61": sqrt(x)*Derivative(y(x), x) - x*exp(-x**(3/2)/3) + (x/4 - 9 + 1/(4*sqrt(x)))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.62": Derivative(y(x), (x, 2)) + (sqrt(x) + x - 8)*y(x)/(4*x**2) - Derivative(y(x), x)/sqrt(x),
        "kamke_2.63": -(2*exp(x) + 1)*Derivative(y(x), x) + y(x)*exp(2*x) - exp(3*x) + Derivative(y(x), (x, 2)),
        "kamke_2.64": a*Derivative(y(x), x) + b*y(x) + tan(x) + Derivative(y(x), (x, 2)),
        "kamke_2.65": 2*n*cot(x)*Derivative(y(x), x) + (-a**2 + n**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.66": y(x)*cos(x)**2 + tan(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.67": -y(x)*cos(x)**2 + tan(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.68": v*(v + 1)*y(x) + cot(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.69": y(x)*sin(x)**2 - cot(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.70": a*tan(x)*Derivative(y(x), x) + b*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.71": 2*a*cot(a*x)*Derivative(y(x), x) + (-a**2 + b**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.72": a*Derivative(p(x), (x, 2))*Derivative(y(x), x) + (-4*a*n*p(x)**2 + a + b*p(x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.73": Subs(Derivative(y(x), (x, 2)) + (-p(x)*Derivative(p(x), x) - Derivative(p(x), (x, 2)) + Derivative(p(x), (x, 3)))*Derivative(y(x), x)/(p(x)**2 + Derivative(p(x), x)) + (-p(x)**2*Derivative(p(x), x) - p(x)*Derivative(p(x), (x, 2)) + Derivative(p(x), x)**2)*y(x)/(p(x)**2 + Derivative(p(x), x)), p(x), WeierstrassP(x, a, b)),
        "kamke_2.74": k**2*JacobiCN(x, k)*JacobiSN(x, k)*Derivative(y(x), x)/JacobiDN(x, k) + n**2*JacobiDN(x, k)**2*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.75": f(x)*Derivative(y(x), x) + g(x)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.76": (a + Derivative(f(x), x))*y(x) + f(x)*Derivative(y(x), x) - g(x) + Derivative(y(x), (x, 2)),
        "kamke_2.77": (a*f(x) + b)*Derivative(y(x), x) + (c*f(x) + d)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.78": (a + f(x)**2/4 + Derivative(f(x), x)/2)*y(x) + f(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_2.79": -a*Derivative(f(x), x)*Derivative(y(x), x)/f(x) + b*f(x)**(2*a)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.80": -(2*a + Derivative(f(x), x)/f(x))*Derivative(y(x), x) + (a**2 + a*Derivative(f(x), x)/f(x) - b**2*f(x)**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.81": -a**2*y(x)*Derivative(f(x), x)**2/(b**2 + f(x)**2) + Derivative(y(x), (x, 2)) + f(x)*Derivative(f(x), (x, 3))*Derivative(y(x), x)/(b**2 + f(x)**2),
        "kamke_2.82": -((2*m - 1)*Derivative(g(x), x)/g(x) + Derivative(g(x), (x, 2))/Derivative(g(x), x))*Derivative(y(x), x) + ((m**2 - v**2)*Derivative(g(x), x)**2/g(x) + Derivative(g(x), x)**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.83": ((1/4 - v**2)*Derivative(g(x), x)**2/g(x)**2 + Derivative(g(x), x)**2 + Derivative(g(x), (x, 3))/(2*Derivative(g(x), x)) - 3*Derivative(g(x), (x, 2))**2/(4*Derivative(g(x), x)**2) - Derivative(f(x), (x, 2))/(2*f(x)) + 3*Derivative(f(x), x)**2/(4*f(x)**2))*y(x) + Derivative(y(x), (x, 2)) - Derivative(f(x), x)*Derivative(y(x), x)/f(x),
        "kamke_2.84": -(Derivative(g(x), (x, 2))/Derivative(g(x), x) - Derivative(g(x), x)/g(x) + 2*Derivative(f(x), x)/f(x))*Derivative(y(x), x) + (-v**2*Derivative(g(x), x)**2/g(x)**2 + (Derivative(g(x), (x, 2))/Derivative(g(x), x) - Derivative(g(x), x)/g(x) + 2*Derivative(f(x), x)/f(x))*Derivative(f(x), x)/f(x) + Derivative(g(x), x)**2 - Derivative(f(x), (x, 2))/f(x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.85": -((2*v - 1)*Derivative(g(x), x)/g(x) + Derivative(g(x), (x, 2))/Derivative(g(x), x) + 2*Derivative(h(x), x)/h(x))*Derivative(y(x), x) + (((2*v - 1)*Derivative(g(x), x)/g(x) + Derivative(g(x), (x, 2))/Derivative(g(x), x) + 2*Derivative(h(x), x)/h(x))*Derivative(h(x), x)/h(x) + Derivative(g(x), x)**2 - Derivative(h(x), (x, 2))/h(x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.86": 9*x*y(x) + 4*Derivative(y(x), (x, 2)),
        "kamke_2.87": -(a + x**2)*y(x) + 4*Derivative(y(x), (x, 2)),
        "kamke_2.88": -(5*tan(x)**2 + 2)*y(x) + 4*tan(x)*Derivative(y(x), x) + 4*Derivative(y(x), (x, 2)),
        "kamke_2.89": a*Derivative(y(x), (x, 2)) + (b*(c + x) + d)*y(x) - (a*b + c + x)*Derivative(y(x), x),
        "kamke_2.90": a**2*Derivative(y(x), (x, 2)) + a*(a**2 - 2*b*exp(-a*x))*Derivative(y(x), x) + b**2*y(x)*exp(-2*a*x),
        "kamke_2.91": x*(y(x) + Derivative(y(x), (x, 2))) - cos(x),
        "kamke_2.92": x*Derivative(y(x), (x, 2)) + (a + x)*y(x),
        "kamke_2.93": x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_2.94": a*y(x) + x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_2.95": l*x*y(x) + x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_2.96": x*Derivative(y(x), (x, 2)) + (a + x)*y(x) + Derivative(y(x), x),
        "kamke_2.97": a*y(x) + x*Derivative(y(x), (x, 2)) - Derivative(y(x), x),
        "kamke_2.98": -a*x**3*y(x) + x*Derivative(y(x), (x, 2)) - Derivative(y(x), x),
        "kamke_2.99": x**3*(-v**2 + exp(x**2))*y(x) + x*Derivative(y(x), (x, 2)) - Derivative(y(x), x),
        "kamke_2.100": -x*y(x) + x*Derivative(y(x), (x, 2)) - exp(x) + 2*Derivative(y(x), x),
        "kamke_2.101": a*x*y(x) + x*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x),
        "kamke_2.102": a*x**2*y(x) + x*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x),
        "kamke_2.103": a*y(x) + x*Derivative(y(x), (x, 2)) - 2*Derivative(y(x), x),
        "kamke_2.104": a*y(x) + v*Derivative(y(x), x) + x*Derivative(y(x), (x, 2)),
        "kamke_2.105": a*Derivative(y(x), x) + b*x*y(x) + x*Derivative(y(x), (x, 2)),
        "kamke_2.106": a*Derivative(y(x), x) + b*x**a1*y(x) + x*Derivative(y(x), (x, 2)),
        "kamke_2.107": a*y(x) + x*Derivative(y(x), (x, 2)) + (b + x)*Derivative(y(x), x),
        "kamke_2.108": a*y(x) + x*Derivative(y(x), (x, 2)) + (a + b + x)*Derivative(y(x), x),
        "kamke_2.109": -x*(x + 1)*exp(x) - x*Derivative(y(x), x) + x*Derivative(y(x), (x, 2)) - y(x),
        "kamke_2.110": -a*y(x) - x*Derivative(y(x), x) + x*Derivative(y(x), (x, 2)),
        "kamke_2.111": x*Derivative(y(x), (x, 2)) - (x + 1)*Derivative(y(x), x) + y(x),
        "kamke_2.112": x*Derivative(y(x), (x, 2)) - (x + 1)*Derivative(y(x), x) - (2*x - 2)*y(x),
        "kamke_2.113": -a*y(x) + x*Derivative(y(x), (x, 2)) + (b - x)*Derivative(y(x), x),
        "kamke_2.114": x*Derivative(y(x), (x, 2)) - (2*x - 2)*Derivative(y(x), x) - y(x),
        "kamke_2.115": x*Derivative(y(x), (x, 2)) - (2*x - 3)*y(x) - (3*x - 2)*Derivative(y(x), x),
        "kamke_2.116": a*n*y(x) + x*Derivative(y(x), (x, 2)) + (a*x + b + n)*Derivative(y(x), x),
        "kamke_2.117": a*b*x*y(x) + x*Derivative(y(x), (x, 2)) - (a + b)*(x + 1)*Derivative(y(x), x),
        "kamke_2.118": x*Derivative(y(x), (x, 2)) + (m + n + x*(a + b))*Derivative(y(x), x) + (a*b*x + a*n + b*m)*y(x),
        "kamke_2.119": x*Derivative(y(x), (x, 2)) - (2*a*x + 2*b)*Derivative(y(x), x) + (a**2*x + 2*a*b)*y(x),
        "kamke_2.120": x*Derivative(y(x), (x, 2)) + (a*x + b)*Derivative(y(x), x) + (c*x + d)*y(x),
        "kamke_2.121": x*Derivative(y(x), (x, 2)) + (x - 1)*y(x) - (x**2 - x)*Derivative(y(x), x),
        "kamke_2.122": -x*(x + 3)*y(x) + x*Derivative(y(x), (x, 2)) - (x**2 - x - 2)*Derivative(y(x), x),
        "kamke_2.123": b*x**3*y(x) + x*Derivative(y(x), (x, 2)) - (2*a*x**2 + 1)*Derivative(y(x), x),
        "kamke_2.124": 2*n*x*y(x) + x*Derivative(y(x), (x, 2)) - (-2*a + 2*x**2)*Derivative(y(x), x),
        "kamke_2.125": -4*x**5 - 4*x**3*y(x) + x*Derivative(y(x), (x, 2)) + (4*x**2 - 1)*Derivative(y(x), x),
        "kamke_2.126": x**2*(a**2*x**3 + a)*y(x) + x*Derivative(y(x), (x, 2)) + (2*a*x**3 - 1)*Derivative(y(x), x),
        "kamke_2.127": x*Derivative(y(x), (x, 2)) + (2*a*x*log(x) + 1)*Derivative(y(x), x) + (a**2*x*log(x)**2 + a*log(x) + a)*y(x),
        "kamke_2.128": x*Derivative(y(x), (x, 2)) + (x*f(x) + 2)*Derivative(y(x), x) + f(x)*y(x),
        "kamke_2.129": (x - 3)*Derivative(y(x), (x, 2)) + (3*x - 6)*y(x) - (4*x - 9)*Derivative(y(x), x),
        "kamke_2.130": a*y(x) + 2*x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_2.131": a*y(x) + 2*x*Derivative(y(x), (x, 2)) - (x - 1)*Derivative(y(x), x),
        "kamke_2.132": a*y(x) + 2*x*Derivative(y(x), (x, 2)) - (2*x - 1)*Derivative(y(x), x),
        "kamke_2.133": (x - 3)*y(x) + (2*x - 1)*Derivative(y(x), (x, 2)) - (3*x - 4)*Derivative(y(x), x),
        "kamke_2.134": 4*x*Derivative(y(x), (x, 2)) - (a + x)*y(x),
        "kamke_2.135": 4*x*Derivative(y(x), (x, 2)) - y(x) + 2*Derivative(y(x), x),
        "kamke_2.136": 4*x*Derivative(y(x), (x, 2)) - (x + 2)*y(x) + 4*Derivative(y(x), x),
        "kamke_2.137": l*y(x) + 4*x*Derivative(y(x), (x, 2)) - (x + 2)*y(x) + 4*y(x),
        "kamke_2.138": 4*m*Derivative(y(x), x) + 4*x*Derivative(y(x), (x, 2)) - (-2*m - 4*n + x)*y(x),
        "kamke_2.139": 16*x*Derivative(y(x), (x, 2)) - (a + x)*y(x) + 8*Derivative(y(x), x),
        "kamke_2.140": a*x*Derivative(y(x), (x, 2)) + b*Derivative(y(x), x) + c*y(x),
        "kamke_2.141": a*x*Derivative(y(x), (x, 2)) + 3*b*y(x) + (3*a + b*x)*Derivative(y(x), x),
        "kamke_2.142": 8*a*Derivative(y(x), x) + c*(a*x + b)**(1/5)*y(x) + (5*a*x + 5*b)*Derivative(y(x), (x, 2)),
        "kamke_2.143": 2*a*x*Derivative(y(x), (x, 2)) + c*y(x) + (a + b*x)*Derivative(y(x), x),
        "kamke_2.144": 2*a*x*Derivative(y(x), (x, 2)) + c*y(x) + (3*a + b*x)*Derivative(y(x), x),
        "kamke_2.145": (a0*x + b0)*y(x) + (a1*x + b1)*Derivative(y(x), x) + (a2*x + b2)*Derivative(y(x), (x, 2)),
        "kamke_2.146": x**2*Derivative(y(x), (x, 2)) - 6*y(x),
        "kamke_2.147": x**2*Derivative(y(x), (x, 2)) - 12*y(x),
        "kamke_2.148": a*y(x) + x**2*Derivative(y(x), (x, 2)),
        "kamke_2.149": x**2*Derivative(y(x), (x, 2)) + (a*x + b)*y(x),
        "kamke_2.150": x**2*Derivative(y(x), (x, 2)) + (x**2 - 2)*y(x),
        "kamke_2.151": x**2*Derivative(y(x), (x, 2)) - (a*x**2 + 2)*y(x),
        "kamke_2.152": x**2*Derivative(y(x), (x, 2)) + (a**2*x**2 - 6)*y(x),
        "kamke_2.153": x**2*Derivative(y(x), (x, 2)) + (a*x**2 - v*(v - 1))*y(x),
        "kamke_2.154": x**2*Derivative(y(x), (x, 2)) + (a*x**2 + b*x + c)*y(x),
        "kamke_2.155": x**2*Derivative(y(x), (x, 2)) + (a*x**k - b*(b - 1))*y(x),
        "kamke_2.156": x**2*Derivative(y(x), (x, 2)) - x*(x*log(x) + 2)*exp(x) + y(x)/log(x),
        "kamke_2.157": a*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) - x*y(x),
        "kamke_2.158": a*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) - (a*b + b**2*x**2)*y(x),
        "kamke_2.159": -a*x**2 + x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) - y(x),
        "kamke_2.160": a*y(x) + x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x),
        "kamke_2.161": x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) - (a + x)*y(x),
        "kamke_2.162": x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) + (-v**2 + x**2)*y(x),
        "kamke_2.163": x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) + (-v**2 + x**2)*y(x) - f(x),
        "kamke_2.164": x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) + (l*x**2 - v**2)*y(x),
        "kamke_2.165": x**2*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) + (-4*v**2 + 4*x**4)*y(x),
        "kamke_2.166": -3*x**3 + x**2*Derivative(y(x), (x, 2)) - x*Derivative(y(x), x) + y(x),
        "kamke_2.167": x**2*Derivative(y(x), (x, 2)) - x*Derivative(y(x), x) + (a*x**m + b)*y(x),
        "kamke_2.168": x**2*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x),
        "kamke_2.169": x**2*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x) + (a*x - b**2)*y(x),
        "kamke_2.170": x**2*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x) + (a*x**2 + b)*y(x),
        "kamke_2.171": x**2*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x) + (a*x + l*x**2 - n*(n + 1))*y(x),
        "kamke_2.172": a*y(x) + x**2*Derivative(y(x), (x, 2)) + (2*x - 2)*Derivative(y(x), x),
        "kamke_2.173": -b*(b - 1)*y(x) + x**2*Derivative(y(x), (x, 2)) + (2*a + 2*x)*Derivative(y(x), x),
        "kamke_2.174": -x**5*log(x) + x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + 2*y(x),
        "kamke_2.175": x**2*Derivative(y(x), (x, 2)) - x*sin(x) - 2*x*Derivative(y(x), x) - (a*x**2 + 12*a + 4)*cos(x) - 4*y(x),
        "kamke_2.176": x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + (x**2 + 2)*y(x),
        "kamke_2.177": x**2*Derivative(y(x), (x, 2)) - x**2/cos(x) - 2*x*Derivative(y(x), x) + (x**2 + 2)*y(x),
        "kamke_2.178": -x**3/cos(x) + x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + (x**2 + 2)*y(x),
        "kamke_2.179": x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + (a**2*x**2 + 2)*y(x),
        "kamke_2.180": x**2*Derivative(y(x), (x, 2)) + 3*x*Derivative(y(x), x) + (-v**2 + x**2 + 1)*y(x) - f(x),
        "kamke_2.181": x**2*Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x) + y(x),
        "kamke_2.182": x**2*Derivative(y(x), (x, 2)) - 3*x*Derivative(y(x), x) - 5*x + 4*y(x),
        "kamke_2.183": -x**2*log(x) + x**2*Derivative(y(x), (x, 2)) - 3*x*Derivative(y(x), x) - 5*y(x),
        "kamke_2.184": -x**4 + x**2*Derivative(y(x), (x, 2)) + x**2 - 4*x*Derivative(y(x), x) + 6*y(x),
        "kamke_2.185": x**2*Derivative(y(x), (x, 2)) + 5*x*Derivative(y(x), x) - (2*x**3 - 4)*y(x),
        "kamke_2.186": -x**3*sin(x) + x**2*Derivative(y(x), (x, 2)) - 5*x*Derivative(y(x), x) + 8*y(x),
        "kamke_2.187": a*x*Derivative(y(x), x) + b*y(x) + x**2*Derivative(y(x), (x, 2)),
        "kamke_2.188": c*y(x) + x**2*Derivative(y(x), (x, 2)) + (a*x + b)*Derivative(y(x), x),
        "kamke_2.189": a*x*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (b*x**m + c)*y(x),
        "kamke_2.190": x**2*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (a*x + b)*y(x),
        "kamke_2.191": x**2*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) - 2*y(x),
        "kamke_2.192": x**2*Derivative(y(x), (x, 2)) + (x**2 - 1)*Derivative(y(x), x) - y(x),
        "kamke_2.193": x**2*Derivative(y(x), (x, 2)) + x*(x + 1)*Derivative(y(x), x) + (x - 9)*y(x),
        "kamke_2.194": x**2*Derivative(y(x), (x, 2)) + x*(x + 1)*Derivative(y(x), x) + (3*x - 1)*y(x),
        "kamke_2.195": x**2*Derivative(y(x), (x, 2)) + x*(x + 3)*Derivative(y(x), x) - y(x),
        "kamke_2.196": x**2*Derivative(y(x), (x, 2)) - x*(x - 1)*Derivative(y(x), x) + (x - 1)*y(x),
        "kamke_2.197": x**2*Derivative(y(x), (x, 2)) - (a + x)*y(x) - (x**2 - 2*x)*Derivative(y(x), x),
        "kamke_2.198": x**2*Derivative(y(x), (x, 2)) - (3*x + 2)*y(x) - (x**2 - 2*x)*Derivative(y(x), x),
        "kamke_2.199": x**2*Derivative(y(x), (x, 2)) - x*(x + 4)*Derivative(y(x), x) + 4*y(x),
        "kamke_2.200": -v*(v - 1)*y(x) + 2*x**2*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)),
        "kamke_2.201": x**2*Derivative(y(x), (x, 2)) + x*(2*x + 1)*Derivative(y(x), x) - 4*y(x),
        "kamke_2.202": x**2*Derivative(y(x), (x, 2)) - 2*x*(x + 1)*Derivative(y(x), x) + (2*x + 2)*y(x),
        "kamke_2.203": a*x**2*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) - 2*y(x),
        "kamke_2.204": x**2*(a + 2*b)*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (b*x**2*(a + b) - 2)*y(x),
        "kamke_2.205": a*x**2*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + f(x)*y(x),
        "kamke_2.206": x**2*Derivative(y(x), (x, 2)) + x*(2*a*x + b)*Derivative(y(x), x) + (a*b*x + c*x**2 + d)*y(x),
        "kamke_2.207": x**2*Derivative(y(x), (x, 2)) + x*(a*x + b)*Derivative(y(x), x) + (a1*x**2 + b1*x + c1)*y(x),
        "kamke_2.208": x**3*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (x**2 - 2)*y(x),
        "kamke_2.209": x**2*Derivative(y(x), (x, 2)) + x*(x**2 + 2)*Derivative(y(x), x) + (x**2 - 2)*y(x),
        "kamke_2.210": x**2*Derivative(y(x), (x, 2)) - 2*x*(-a + x**2)*Derivative(y(x), x) + (a*((-1)**n - 1) + 2*n*x**2)*y(x),
        "kamke_2.211": 4*x**3*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (4*x**4 + 2*x**2 + 1)*y(x),
        "kamke_2.212": x**2*Derivative(y(x), (x, 2)) + x*(a*x**2 + b)*Derivative(y(x), x) + f(x)*y(x),
        "kamke_2.213": x**2*Derivative(y(x), (x, 2)) + x*(x**3 + 1)*Derivative(y(x), x) - y(x),
        "kamke_2.214": x**2*Derivative(y(x), (x, 2)) + ((-1)**n*a - a**2 - x**4 + x**2*(2*a + 2*n + 1))*y(x),
        "kamke_2.215": x**2*Derivative(y(x), (x, 2)) + x*(a*x**n + b)*Derivative(y(x), x) + (a1*x**(2*n) + b1*x**n + c1)*y(x),
        "kamke_2.216": x**2*Derivative(y(x), (x, 2)) + x*(a*x**a1 + b)*Derivative(y(x), x) + (A*x**(2*a1) + B*x**a1 + C*x**b1 + DD)*y(x),
        "kamke_2.217": x**2*Derivative(y(x), (x, 2)) - (a + x*tan(x))*y(x) - (2*x**2*tan(x) - x)*Derivative(y(x), x),
        "kamke_2.218": x**2*Derivative(y(x), (x, 2)) + (a + x*cot(x))*y(x) + (2*x**2*cot(x) + x)*Derivative(y(x), x),
        "kamke_2.219": x**2*Derivative(y(x), (x, 2)) + 2*x*f(x)*Derivative(y(x), x) + (a*x**2 + b*x + c + x*Derivative(f(x), x) + f(x)**2 - f(x))*y(x),
        "kamke_2.220": 2*x**2*f(x)*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 2)) + (-v*(v - 1) + x**2*(a + f(x)**2 + Derivative(f(x), x)))*y(x),
        "kamke_2.221": x**2*Derivative(y(x), (x, 2)) + (-2*x**2*f(x) + x)*Derivative(y(x), x) + (-v**2 + x**2*(f(x)**2 - Derivative(f(x), x) + 1) - x*f(x))*y(x),
        "kamke_2.222": x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)) + 2*y(x),
        "kamke_2.223": x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)) - 9*y(x),
        "kamke_2.224": a*y(x) + x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.225": -x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)) + y(x),
        "kamke_2.226": -v*(v - 1)*y(x) + 2*x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.227": -2*x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)) + 2*y(x),
        "kamke_2.228": a*y(x) + 3*x*Derivative(y(x), x) + (x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.229": 4*x*Derivative(y(x), x) + 2*x + (x**2 + 1)*Derivative(y(x), (x, 2)) + 2*y(x) - 2*cos(x),
        "kamke_2.230": a*x*Derivative(y(x), x) + (a - 2)*y(x) + (x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.231": -v*(v + 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.232": -n*(n + 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + Derivative(LegendreP(n, x), x),
        "kamke_2.233": -n*(n + 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + Derivative(LegendreQ(n, x), x),
        "kamke_2.234": x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + 2,
        "kamke_2.235": a*y(x) + x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.236": x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + f(x)*y(x),
        "kamke_2.237": 2*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.238": -a + 2*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.239": -l*y(x) + 2*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.240": -v*(v + 1)*y(x) + 2*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.241": -2*x*Derivative(y(x), x) - (v - 1)*(v + 2)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.242": -(3*x + 1)*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)) - (x**2 - x)*y(x),
        "kamke_2.243": 4*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + (x**2 + 1)*y(x),
        "kamke_2.244": x*(2*n + 2)*Derivative(y(x), x) - (-n + v)*(n + v + 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.245": -x*(2*n - 2)*Derivative(y(x), x) - (n + v)*(-n + v + 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.246": -2*v*y(x) - x*(2*v - 2)*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.247": 2*a*x*Derivative(y(x), x) + a*(a - 1)*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.248": a*x*Derivative(y(x), x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + (b*x**2 + c*x + d)*y(x),
        "kamke_2.249": c*y(x) + (x**2 - 1)*Derivative(y(x), (x, 2)) + (a*x + b)*Derivative(y(x), x),
        "kamke_2.250": 8*x*Derivative(y(x), x) + (-a**2 + x**2)*Derivative(y(x), (x, 2)) + 12*y(x),
        "kamke_2.251": x*(x + 1)*Derivative(y(x), (x, 2)) - (x - 1)*Derivative(y(x), x) + y(x),
        "kamke_2.252": c*y(x) + x*(x + 1)*Derivative(y(x), (x, 2)) + (a*x + b)*Derivative(y(x), x),
        "kamke_2.253": x*(x + 1)*Derivative(y(x), (x, 2)) + (3*x + 2)*Derivative(y(x), x) + y(x),
        "kamke_2.254": (x**2 - x)*Derivative(y(x), x) - (6*x**2 + 7*x)*y(x) + (x**2 + x - 2)*Derivative(y(x), (x, 2)),
        "kamke_2.255": a*Derivative(y(x), x) + x*(x - 1)*Derivative(y(x), (x, 2)) - 2*y(x),
        "kamke_2.256": -v*(v + 1)*y(x) + x*(x - 1)*Derivative(y(x), (x, 2)) + (2*x - 1)*Derivative(y(x), x),
        "kamke_2.257": x*(x - 1)*Derivative(y(x), (x, 2)) + (b + x*(a + 1))*Derivative(y(x), x),
        "kamke_2.258": c*y(x) + x*(x - 1)*Derivative(y(x), (x, 2)) + (a*x + b)*Derivative(y(x), x),
        "kamke_2.259": -l*y(x) + x*(x - 1)*Derivative(y(x), (x, 2)) + (b + x*(a + 1))*Derivative(y(x), x),
        "kamke_2.260": a1*b1*d1 + x*(x - 1)*Derivative(y(x), (x, 2)) + (-d1 + x*(a1 + b1 + 1))*Derivative(y(x), x),
        "kamke_2.261": x*(x + 2)*Derivative(y(x), (x, 2)) + (2*l*ps + 2*l*x*(-n + ps - 1) + m)*y(x) + (-2*l*x**2 + 2*n + 2*x*(-2*l + n + 1) + 2)*Derivative(y(x), x),
        "kamke_2.262": (x + 1)**2*Derivative(y(x), (x, 2)) - (x + 2)*y(x) + (x**2 + x - 1)*Derivative(y(x), x),
        "kamke_2.263": x*(x + 3)*Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x) - (20*x + 30)*(x**2 + 3*x)**(7/3) + y(x),
        "kamke_2.264": -(2*x + 3)*y(x) + (x**2 + x + 1)*Derivative(y(x), x) + (x**2 + 3*x + 4)*Derivative(y(x), (x, 2)),
        "kamke_2.265": (x - 2)*(x - 1)*Derivative(y(x), (x, 2)) - (2*x - 3)*Derivative(y(x), x) + y(x),
        "kamke_2.266": (x - 2)**2*Derivative(y(x), (x, 2)) - (x - 2)*Derivative(y(x), x) - 3*y(x),
        "kamke_2.267": 2*x**2*Derivative(y(x), (x, 2)) - (4*x - 1)*y(x) - (l + 2*x**2 - 5*x)*Derivative(y(x), x),
        "kamke_2.268": 2*x*(x - 1)*Derivative(y(x), (x, 2)) + (2*x - 1)*Derivative(y(x), x) + (a*x + b)*y(x),
        "kamke_2.269": 2*x*(x - 1)*Derivative(y(x), (x, 2)) + (v + 1)*y(x) + (-2*v + x*(2*v + 5) - 3)*Derivative(y(x), x),
        "kamke_2.270": (2*x**2 + 6*x + 4)*Derivative(y(x), (x, 2)) + (10*x**2 + 21*x + 8)*Derivative(y(x), x) + (12*x**2 + 17*x + 8)*y(x),
        "kamke_2.271": 4*x**2*Derivative(y(x), (x, 2)) + y(x),
        "kamke_2.272": 4*x**2*Derivative(y(x), (x, 2)) + (4*a**2*x**2 + 1)*y(x),
        "kamke_2.273": 4*x**2*Derivative(y(x), (x, 2)) - (-4*k*x + 4*m**2 + x**2 - 1)*y(x),
        "kamke_2.274": 4*x**2*Derivative(y(x), (x, 2)) + 4*x*Derivative(y(x), x) + (-v**2 + x)*y(x),
        "kamke_2.275": 4*x**2*Derivative(y(x), (x, 2)) + 4*x*Derivative(y(x), x) + (-m**2 - x**2 + x*(4*l - 2*m + 2) + 1)*y(x),
        "kamke_2.276": 4*x**2*Derivative(y(x), (x, 2)) + 4*x*Derivative(y(x), x) - (4*x**2 + 1)*y(x) - 4*sqrt(x**3)*exp(x),
        "kamke_2.277": 4*x**2*Derivative(y(x), (x, 2)) + 4*x*Derivative(y(x), x) - (a*x**2 + 1)*y(x),
        "kamke_2.278": 4*x**2*Derivative(y(x), (x, 2)) + 4*x*Derivative(y(x), x) + f(x)*y(x),
        "kamke_2.279": 4*x**2*Derivative(y(x), (x, 2)) + 5*x*Derivative(y(x), x) - y(x) - log(x),
        "kamke_2.280": 4*x**2*Derivative(y(x), (x, 2)) + 8*x*Derivative(y(x), x) - (4*x**2 + 12*x + 3)*y(x),
        "kamke_2.281": 4*x**2*Derivative(y(x), (x, 2)) - 4*x*(2*x - 1)*Derivative(y(x), x) + (4*x**2 - 4*x - 1)*y(x),
        "kamke_2.282": 4*x**3*Derivative(y(x), x) + 4*x**2*Derivative(y(x), (x, 2)) + (x**2 - 4)*(x**2 + 6)*y(x),
        "kamke_2.283": -4*x**2*sqrt(exp(x)/x**x) + 4*x**2*log(x)*Derivative(y(x), x) + 4*x**2*Derivative(y(x), (x, 2)) + (x**2*log(x)**2 + 2*x - 8)*y(x),
        "kamke_2.284": -3*x + (2*x + 1)**2*Derivative(y(x), (x, 2)) - (4*x + 2)*Derivative(y(x), x) - 12*y(x) - 1,
        "kamke_2.285": a*(a - 1)*y(x) + x*(4*x - 1)*Derivative(y(x), (x, 2)) + (-a + x*(4*a + 2))*Derivative(y(x), x),
        "kamke_2.286": (3*x - 1)**2*Derivative(y(x), (x, 2)) + (9*x - 3)*Derivative(y(x), x) - 9*y(x) - log(3*x - 1)**2,
        "kamke_2.287": 9*x*(x - 1)*Derivative(y(x), (x, 2)) + (6*x - 3)*Derivative(y(x), x) - 20*y(x),
        "kamke_2.288": 16*x**2*Derivative(y(x), (x, 2)) + (4*x + 3)*y(x),
        "kamke_2.289": 16*x**2*Derivative(y(x), (x, 2)) + 32*x*Derivative(y(x), x) - (4*x + 5)*y(x),
        "kamke_2.290": 27*x*Derivative(y(x), x) + (27*x**2 + 4)*Derivative(y(x), (x, 2)) - 3*y(x),
        "kamke_2.291": 48*x*(x - 1)*Derivative(y(x), (x, 2)) + (152*x - 40)*Derivative(y(x), x) + 53*y(x),
        "kamke_2.292": 50*x*(x - 1)*Derivative(y(x), (x, 2)) + (50*x - 25)*Derivative(y(x), x) - 2*y(x),
        "kamke_2.293": 144*x*(x - 1)*Derivative(y(x), (x, 2)) + (120*x - 48)*Derivative(y(x), x) + y(x),
        "kamke_2.294": 144*x*(x - 1)*Derivative(y(x), (x, 2)) + (168*x - 96)*Derivative(y(x), x) + y(x),
        "kamke_2.295": a*x**2*Derivative(y(x), (x, 2)) + b*x*Derivative(y(x), x) + (c*x**2 + d*x + fs)*y(x),
        "kamke_2.296": a2*x**2*Derivative(y(x), (x, 2)) + (a1*x**2 + b1*x)*Derivative(y(x), x) + (a0*x**2 + b0*x + c0)*y(x),
        "kamke_2.297": a*x*Derivative(y(x), x) + b*y(x) + (a*x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.298": b*x*Derivative(y(x), x) + c*y(x) + (a*x**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_2.299": 2*a**2*x*Derivative(y(x), x) + (a**2*x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.300": 2*a**2*x*Derivative(y(x), x) - 2*a**2*y(x) + (a**2*x**2 - 1)*Derivative(y(x), (x, 2)),
        "kamke_2.301": -2*a*y(x) + 2*b*Derivative(y(x), x) + (a*x**2 + b*x)*Derivative(y(x), (x, 2)),
        "kamke_2.302": A0*(a*x + b)*y(x) + A1*(a*x + b)*Derivative(y(x), x) + A2*(a*x + b)**2*Derivative(y(x), (x, 2)),
        "kamke_2.303": gs*y(x) + (d*x + fs)*Derivative(y(x), x) + (a*x**2 + b*x + c)*Derivative(y(x), (x, 2)),
        "kamke_2.304": x**3*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x) - (2*x + 3)*y(x),
        "kamke_2.305": x**3*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x) - y(x),
        "kamke_2.306": x**3*Derivative(y(x), (x, 2)) + x**2*Derivative(y(x), x) + (a*x**2 + a + b*x)*y(x),
        "kamke_2.307": x**3*Derivative(y(x), (x, 2)) + x*(x + 1)*Derivative(y(x), x) - 2*y(x),
        "kamke_2.308": x**3*Derivative(y(x), (x, 2)) - x**2*Derivative(y(x), x) + x*y(x) - log(x)**3,
        "kamke_2.309": x**3*Derivative(y(x), (x, 2)) + x*y(x) - (x**2 - 1)*Derivative(y(x), x),
        "kamke_2.310": x**3*Derivative(y(x), (x, 2)) + 3*x**2*Derivative(y(x), x) + x*y(x) - 1,
        "kamke_2.311": -v*x*(v + 1)*y(x) + x*(x**2 + 1)*Derivative(y(x), (x, 2)) + (2*x**2 + 1)*Derivative(y(x), x),
        "kamke_2.312": x*(x**2 + 1)*Derivative(y(x), (x, 2)) - 2*x*y(x) + (2*x**2 - 2)*Derivative(y(x), x),
        "kamke_2.313": -x*(-n + v)*(n + v + 1)*y(x) + x*(x**2 + 1)*Derivative(y(x), (x, 2)) + (2*n + x**2*(2*n + 2) + 1)*Derivative(y(x), x),
        "kamke_2.314": x*(n + v)*(n - v - 1)*y(x) + x*(x**2 + 1)*Derivative(y(x), (x, 2)) - (2*n + x**2*(2*n - 2) - 1)*Derivative(y(x), x),
        "kamke_2.315": a*x**3*y(x) + x*(x**2 - 1)*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_2.316": x*(x**2 - 1)*Derivative(y(x), (x, 2)) - x*y(x) + (x**2 - 1)*Derivative(y(x), x),
        "kamke_2.317": x*(x**2 - 1)*Derivative(y(x), (x, 2)) + x*y(x) + (3*x**2 - 1)*Derivative(y(x), x),
        "kamke_2.318": c*x*y(x) + x*(x**2 - 1)*Derivative(y(x), (x, 2)) + (a*x**2 + b)*Derivative(y(x), x),
        "kamke_2.319": x*(x**2 + 2)*Derivative(y(x), (x, 2)) - 6*x*y(x) - Derivative(y(x), x),
        "kamke_2.320": x*(x**2 - 2)*Derivative(y(x), (x, 2)) + (x**2 + 4*x + 2)*y(x) - (x**3 + 3*x**2 - 2*x - 2)*Derivative(y(x), x),
        "kamke_2.321": x**2*(x + 1)*Derivative(y(x), (x, 2)) - x*(2*x + 1)*Derivative(y(x), x) + (2*x + 1)*y(x),
        "kamke_2.322": x**2*(x + 1)*Derivative(y(x), (x, 2)) + 2*x*(3*x + 2)*Derivative(y(x), x),
        "kamke_2.323": Derivative(y(x), (x, 2)) + 2*(x - 2)*Derivative(y(x), x)/(x*(x - 1)) - 2*(x + 1)*y(x)/(x**2*(x - 1)),
        "kamke_2.324": Derivative(y(x), (x, 2)) - (5*x - 4)*Derivative(y(x), x)/(x*(x - 1)) + (9*x - 6)*y(x)/(x**2*(x - 1)),
        "kamke_2.325": Derivative(y(x), (x, 2)) - (-alpha - bbeta - x*(a + b + 1) + 1)*Derivative(y(x), x)/(x*(x - 1)) + (a*b*x - alpha*bbeta)*y(x)/(x**2*(x - 1)),
        "kamke_2.326": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/(x + 1) + y(x)/(x*(x + 1)**2),
        "kamke_2.327": Derivative(y(x), (x, 2)) - 2*Derivative(y(x), x)/(x*(x - 2)) + y(x)/(x**2*(x - 2)),
        "kamke_2.328": Derivative(y(x), (x, 2)) - 2*y(x)/(x*(x - 1)**2),
        "kamke_2.329": Derivative(y(x), (x, 2)) + (alpha*bbeta*x - qs)*y(x)/(x*(-a + x)*(x - 1)) - (-a*ggamma - x**2*(alpha + bbeta + 1) + x*(a*(delta + ggamma) + alpha + bbeta - delta + 1))*Derivative(y(x), x)/(x*(-a + x)*(x - 1)),
        "kamke_2.330": Derivative(y(x), (x, 2)) + (DD*x + EE)*y(x)/((-a + x)*(-b + x)*(-c + x)) - (-A*x**2 - B*x - C)*Derivative(y(x), x)/((-a + x)*(-b + x)*(-c + x)),
        "kamke_2.331": Derivative(y(x), (x, 2)) - (x - 4)*Derivative(y(x), x)/(2*x*(x - 2)) + (x/2 - 3/2)*y(x)/(x**2*(x - 2)),
        "kamke_2.332": Derivative(y(x), (x, 2)) - Derivative(y(x), x)/(x + 1) + (3*x/4 + 1/4)*y(x)/(x**2*(x + 1)),
        "kamke_2.333": -v*(v + 1)*y(x)/(4*x**2) + Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x)/(2*x*(x - 1)),
        "kamke_2.334": Derivative(y(x), (x, 2)) - (-x*(a + 1) + 1)*Derivative(y(x), x)/(x*(x - 1)) + (c**2/4 + x*(a**2 - b**2)/4)*y(x)/(x**2*(x - 1)),
        "kamke_2.335": Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x)/(2*x*(x - 1)) + (a*x/4 + b/4)*y(x)/(x*(x - 1)**2),
        "kamke_2.336": Derivative(y(x), (x, 2)) - (3*x - 1)*y(x)/((x - 1)*(2*x - 1)**2),
        "kamke_2.337": (a/4 - b/4)*y(x)/((a + x)**2*(b + x)) + Derivative(y(x), (x, 2)) + (a + 2*b + 3*x)*Derivative(y(x), x)/(2*(a + x)*(b + x)),
        "kamke_2.338": Derivative(y(x), (x, 2)) - (6*x - 1)*Derivative(y(x), x)/(3*x*(x - 2)) - y(x)/(3*x**2*(x - 2)),
        "kamke_2.339": Derivative(y(x), (x, 2)) + (a*b*x - c*d)*y(x)/(x**2*(a*x + 1)) - (-a*x**2*(b + 2) - x*(c - d + 1))*Derivative(y(x), x)/(x**2*(a*x + 1)),
        "kamke_2.340": Derivative(y(x), (x, 2)) - 2*(a*x + 2*b)*Derivative(y(x), x)/(x*(a*x + b)) + (2*a*x + 6*b)*y(x)/(x**2*(a*x + b)),
        "kamke_2.341": -A*x + Derivative(y(x), (x, 2)) + (2*a*x + b)*Derivative(y(x), x)/(x*(a*x + b)) + (a*v*x - b)*y(x)/(x**2*(a*x + b)),
        "kamke_2.342": a*y(x)/x**4 + Derivative(y(x), (x, 2)),
        "kamke_2.343": Derivative(y(x), (x, 2)) - (-a*x**2*(1 - a) + b*(b + x))*y(x)/x**4,
        "kamke_2.344": Derivative(y(x), (x, 2)) - (v**2 - exp(2/x))*y(x)/x**4,
        "kamke_2.345": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x**3 - 2*y(x)/x**4,
        "kamke_2.346": Derivative(y(x), (x, 2)) - (a + b)*Derivative(y(x), x)/x**2 + (a*b + x*(a + b))*y(x)/x**4,
        "kamke_2.347": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x + y(x)/x**4,
        "kamke_2.348": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x + (a*(x**4 + 1) + b*x**2)*y(x)/x**4,
        "kamke_2.349": Derivative(y(x), (x, 2)) - (-x**2 - 1)*Derivative(y(x), x)/x**3 + y(x)/x**4,
        "kamke_2.350": a**2*y(x)/x**4 + Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x)/x,
        "kamke_2.351": Derivative(y(x), (x, 2)) - (-2*x**2 - 1)*Derivative(y(x), x)/x**3 - y(x)/x**4,
        "kamke_2.352": b*y(x)/x**4 + Derivative(y(x), (x, 2)) + 2*(a + x)*Derivative(y(x), x)/x**2,
        "kamke_2.353": Derivative(y(x), (x, 2)) - (2*x**2 - 1)*Derivative(y(x), x)/x**3 + y(x)/x**4,
        "kamke_2.354": Derivative(y(x), (x, 2)) - (2*x**2 - 1)*Derivative(y(x), x)/x**3 + 2*y(x)/x**4,
        "kamke_2.355": -x*y(x)/(x**3 + 1) + Derivative(y(x), (x, 2)) - (1 - x**3)*Derivative(y(x), x)/(x*(x**3 + 1)),
        "kamke_2.356": Derivative(y(x), (x, 2)) - (-2*x**2 - 1)*Derivative(y(x), x)/(x*(x**2 + 1)) + (-n**2 - v*x**2*(v + 1))*y(x)/(x**2*(x**2 + 1)),
        "kamke_2.357": Derivative(y(x), (x, 2)) + (a*x**2 + a - 1)*Derivative(y(x), x)/(x*(x**2 + 1)) + (b*x**2 + c)*y(x)/(x**2*(x**2 + 1)),
        "kamke_2.358": Derivative(y(x), (x, 2)) - (x**2 - 2)*Derivative(y(x), x)/(x*(x**2 - 1)) + (x**2 - 2)*y(x)/(x**2*(x**2 - 1)),
        "kamke_2.359": v*(v + 1)*y(x)/(x**2*(x**2 - 1)) + 2*x*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)),
        "kamke_2.360": -v*(v + 1)*y(x)/x**2 + 2*x*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)),
        "kamke_2.361": -2*x*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)) + (2*a*x**2 + n*(n + 1)*(x**2 - 1) + x**2*(a - n)*(x**2 - 1)*(a + n + 1))*y(x)/(x**2*(x**2 - 1)),
        "kamke_2.362": -2*x**3*Derivative(y(x), x) + x**2*(x**2 - 1)*Derivative(y(x), (x, 2)) - (2*a*x**2 + n*(n + 1)*(x**2 - 1) + x**2*(a - n)*(x**2 - 1)*(a + n + 1))*y(x),
        "kamke_2.363": b*y(x)/x**2 + Derivative(y(x), (x, 2)) + (a*x**2 + a - 2)*Derivative(y(x), x)/(x*(x**2 - 1)),
        "kamke_2.364": Derivative(y(x), (x, 2)) - (-2*a + 2*b*c*x**c*(x**2 - 1) + x**2*(2*a - 2))*Derivative(y(x), x)/(x*(x**2 - 1)) + (-a*(a + 1) + b**2*c**2*x**(2*c)*(x**2 - 1) - b*c*x**c*(2*a - c + 1) + b*c*x**(c + 2)*(2*a - c - 1) + x**2*(a*(a - 1) - v*(v + 1)))*y(x)/(x**2*(x**2 - 1)),
        "kamke_2.365": a*y(x)/(x**2 + 1)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.366": 2*x*Derivative(y(x), x)/(x**2 + 1) + Derivative(y(x), (x, 2)) + y(x)/(x**2 + 1)**2,
        "kamke_2.367": 2*x*Derivative(y(x), x)/(x**2 + 1) + Derivative(y(x), (x, 2)) + (a**2*(x**2 + 1)**2 + m**2 - n*(n + 1)*(x**2 + 1))*y(x)/(x**2 + 1)**2,
        "kamke_2.368": a*x*Derivative(y(x), x)/(x**2 + 1) + b*y(x)/(x**2 + 1)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.369": a*y(x)/(x**2 - 1)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.370": -a**2*y(x)/(x**2 - 1)**2 + 2*x*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)),
        "kamke_2.371": 2*x*Derivative(y(x), x)/(x**2 - 1) + (-a**2 - lambda_*(x**2 - 1))*y(x)/(x**2 - 1)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.372": 2*x*Derivative(y(x), x)/(x**2 - 1) + (-k**2 + (x**2 - 1)*(a*x**2 + b*x + c))*y(x)/(x**2 - 1)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.373": 2*x*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)) + (-a**2*(x**2 - 1)**2 - m**2 - n*(n + 1)*(x**2 - 1))*y(x)/(x**2 - 1)**2,
        "kamke_2.374": -2*x*(2*a - 1)*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)) + (2*a + v*(v + 1) + x**2*(2*a*(2*a - 1) - v*(v + 1)))*y(x)/(x**2 - 1)**2,
        "kamke_2.375": 2*x*(-2*a + n + 1)*Derivative(y(x), x)/(x**2 - 1) + Derivative(y(x), (x, 2)) + (4*a*x**2*(a - n) - (2*a + (-n + v)*(n + v + 1))*(x**2 - 1))*y(x)/(x**2 - 1)**2,
        "kamke_2.376": b*y(x)/(x**2*(a + x**2)) + Derivative(y(x), (x, 2)) + (a + 2*x**2)*Derivative(y(x), x)/(x*(a + x**2)),
        "kamke_2.377": b**2*y(x)/(a**2 + x**2)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.378": Derivative(y(x), (x, 2)) + 2*(x**2 - 1)*Derivative(y(x), x)/(x*(x - 1)**2) + (-2*x**2 + 2*x + 2)*y(x)/(x**2*(x - 1)**2),
        "kamke_2.379": Derivative(y(x), (x, 2)) - 12*y(x)/((x + 1)**2*(x**2 + 2*x + 3)),
        "kamke_2.380": b*y(x)/(x**2*(-a + x)**2) + Derivative(y(x), (x, 2)),
        "kamke_2.381": b*y(x)/(x**2*(-a + x)**2) - c + Derivative(y(x), (x, 2)),
        "kamke_2.382": -c*y(x)/((-a + x)**2*(-b + x)**2) + Derivative(y(x), (x, 2)),
        "kamke_2.383": alpha*bbeta*(a - b)**2*y(x)/((-a + x)**2*(-b + x)**2) + Derivative(y(x), (x, 2)) - (-(-a + x)**2*(-b + x)*(alpha + bbeta + 1) - (-a + x)*(-b + x)**2*(-alpha - bbeta + 1))*Derivative(y(x), x)/((-a + x)**2*(-b + x)**2),
        "kamke_2.384": Derivative(y(x), (x, 2)) - (b**2/4 - b*x*(2*a + 6)/4 + x**2*(a**2 - 1)/4)*y(x)/x**2,
        "kamke_2.385": Derivative(y(x), (x, 2)) - (-a*x**2/4 - a/4 + 3/4)*y(x)/(x**2 + 1)**2,
        "kamke_2.386": Derivative(y(x), (x, 2)) - 18*y(x)/((2*x + 1)**2*(x**2 + x + 1)),
        "kamke_2.387": Derivative(y(x), (x, 2)) - 3*y(x)/(4*(x**2 + x + 1)**2),
        "kamke_2.388": Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x)/(2*x*(x - 1)) + (-a**2*x/4 + v*(v + 1)*(x - 1)/4)*y(x)/(x**2*(x - 1)**2),
        "kamke_2.389": Derivative(y(x), (x, 2)) + (3*x - 1)*Derivative(y(x), x)/(2*x*(x - 1)) + (-n**2*x - v*(v + 1)*(x - 1)**2/4)*y(x)/(x**2*(x - 1)**2),
        "kamke_2.390": Derivative(y(x), (x, 2)) + 3*y(x)/(16*x**2*(x - 1)**2),
        "kamke_2.391": Derivative(y(x), (x, 2)) - (7*a*x**2 + 5)*Derivative(y(x), x)/(x*(a*x**2 + 1)) + (15*a*x**2 + 5)*y(x)/(x**2*(a*x**2 + 1)),
        "kamke_2.392": Derivative(y(x), (x, 2)) + b*x*Derivative(y(x), x)/(a*(x**2 - 1)) + (c*x**2 + d*x + e)*y(x)/(a*(x**2 - 1)**2),
        "kamke_2.393": Derivative(y(x), (x, 2)) - (-b*x**2 - c*x - d)*y(x)/(a*x**2*(x - 1)**2),
        "kamke_2.394": c*y(x)/(x**2*(a*x + b)**2) + Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x)/x,
        "kamke_2.395": Derivative(y(x), (x, 2)) + y(x)/(a*x + b)**4,
        "kamke_2.396": A*y(x)/(a*x**2 + b*x + c)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.397": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x**4 - y(x)/x**5,
        "kamke_2.398": Derivative(y(x), (x, 2)) + (x**2 - (2*v + 1)**2 - 1)*y(x)/(x**2 - 1)**2 + (3*x**2 - 1)*Derivative(y(x), x)/(x*(x**2 - 1)),
        "kamke_2.399": Derivative(y(x), (x, 2)) - (3*x + 1)*Derivative(y(x), x)/((x - 1)*(x + 1)) + 36*(x + 1)**2*y(x)/((x - 1)**2*(3*x + 5)**2),
        "kamke_2.400": a*y(x)/x**6 + Derivative(y(x), (x, 2)) - Derivative(y(x), x)/x,
        "kamke_2.401": b*y(x)/x**6 + Derivative(y(x), (x, 2)) + (a + 3*x**2)*Derivative(y(x), x)/x**3,
        "kamke_2.402": Derivative(y(x), (x, 2)) + (x**2*(1 - 4*a) - 1)*Derivative(y(x), x)/(x*(x**2 - 1)) + (4*a*x**4*(a + 1) - 2*a*x**2*(x**2 - 1) + (-v**2 + x**2)*(x**2 - 1)**2)*y(x)/(x**2*(x**2 - 1)**2),
        "kamke_2.403": -(-(-a3 - b3 + 1)/(-c3 + x) - (-a2 - b2 + 1)/(-c2 + x) - (-a1 - b1 + 1)/(-c1 + x))*Derivative(y(x), x) + Derivative(y(x), (x, 2)) + (a1*b1*(c1 - c2)*(c1 - c3)/(-c1 + x) + a2*b2*(-c1 + c2)*(c2 - c3)/(-c2 + x) + a3*b3*(-c1 + c3)*(-c2 + c3)/(-c3 + x))*y(x)/((-c1 + x)*(-c2 + x)*(-c3 + x)),
        "kamke_2.404": Derivative(y(x), (x, 2)) - (-2*x**2 - 1)*Derivative(y(x), x)/x**3 + (1/4 - x**2/2)*y(x)/x**6,
        "kamke_2.405": Derivative(y(x), (x, 2)) - (2*x**2 + 1)*Derivative(y(x), x)/x**3 + (a*x**4/4 + 5*x**2/2 + 1/4)*y(x)/x**6,
        "kamke_2.406": 27*x*y(x)/(16*(x**3 - 1)**2) + Derivative(y(x), (x, 2)),
        "kamke_2.407": -(-b1*(-al1 - bl1 + 1)/(-a1 + b1*x) - b2*(-al2 - bl2 + 1)/(-a2 + b2*x) - b3*(-al3 - bl3 + 1)/(-a3 + b3*x))*Derivative(y(x), x) + Derivative(y(x), (x, 2)) + (al1*bl1*(a1*b2 - a2*b1)*(-a1*b3 + a3*b1)/(-a1 + b1*x) + al2*bl2*(a1*b2 - a2*b1)*(a2*b3 - a3*b2)/(-a2 + b2*x) + al3*bl3*(-a1*b3 + a3*b1)*(a2*b3 - a3*b2)/(-a3 + b3*x))*y(x)/((-a1 + b1*x)*(-a2 + b2*x)*(-a3 + b3*x)),
        "kamke_2.408": Derivative(y(x), (x, 2)) + (A*x**2 + B)*y(x)/(x*(-a1 + x**2)*(-a2 + x**2)*(-a3 + x**2)) - (-x**2*((-a1 + x**2)*(-a2 + x**2) + (-a1 + x**2)*(-a3 + x**2) + (-a2 + x**2)*(-a3 + x**2)) + (-a1 + x**2)*(-a2 + x**2)*(-a3 + x**2))*Derivative(y(x), x)/(x*(-a1 + x**2)*(-a2 + x**2)*(-a3 + x**2)),
        "kamke_2.409": a*x**(2*a - 1)*Derivative(y(x), x)/x**(2*a) + b**2*y(x)/x**(2*a) + Derivative(y(x), (x, 2)),
        "kamke_2.410": Derivative(y(x), (x, 2)) - (-a*ps*x**b - qs)*Derivative(y(x), x)/(x*(a*x**b - 1)) + (a*r*x**b + s)*y(x)/(x**2*(a*x**b - 1)),
        "kamke_2.411": Derivative(y(x), (x, 2)) - y(x)/(exp(x) + 1),
        "kamke_2.412": -y(x)*log(x)**2 + Derivative(y(x), (x, 2)) - Derivative(y(x), x)/(x*log(x)),
        "kamke_2.413": Derivative(y(x), (x, 2)) - Derivative(y(x), x)/(x*(log(x) - 1)) + y(x)/(x**2*(log(x) - 1)),
        "kamke_2.414": -(a**2*sinh(x)**2 + n*(n - 1))*y(x)/sinh(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.415": 2*n*cosh(x)*Derivative(y(x), x)/sinh(x) + (-a**2 + n**2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.416": -(-2*n - 1)*cos(x)*Derivative(y(x), x)/sin(x) + (-n + v)*(n + v + 1)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_2.417": -(-sin(x)**2 + cos(x))*Derivative(y(x), x)/sin(x) + y(x)*sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.418": x*sin(x)*Derivative(y(x), x)/(x*cos(x) - sin(x)) + Derivative(y(x), (x, 2)) - y(x)*sin(x)/(x*cos(x) - sin(x)),
        "kamke_2.419": Derivative(y(x), (x, 2)) + (-x*sin(x) + 2*cos(x))*y(x)/(x**2*cos(x)) - (-x**2*sin(x) + 2*x*cos(x))*Derivative(y(x), x)/(x**2*cos(x)),
        "kamke_2.420": -(a*cos(x)**2 + n*(n - 1))*y(x) + cos(x)**2*Derivative(y(x), (x, 2)),
        "kamke_2.421": a**2*n*((n - 1)*sin(a*x)**2 + cos(a*x)**2)*y(x)/cos(a*x)**2 + a*(n - 1)*sin(2*a*x)*Derivative(y(x), x)/cos(a*x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.422": -2*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.423": a*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.424": -(a*sin(x)**2 + n*(n - 1))*y(x) + sin(x)**2*Derivative(y(x), (x, 2)),
        "kamke_2.425": -(a**2*cos(x)**2 - 3*a + (3 - 2*a)*cos(x) + 3)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.426": -(a**2*cos(x)**2 + 3*a + b**2/(2*a - 3)**2 + b*cos(x) + 2)*y(x) + sin(x)**2*Derivative(y(x), (x, 2)),
        "kamke_2.427": -(a*b*(a + 1)*sin(2*x) + a*(a - 1) - (-a**2*b**2 + (a + 1)**2)*sin(x)**2)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.428": -(-a*cos(x)**2 - b*sin(x)**2 - c)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.429": -y(x)/sin(x)**2 + Derivative(y(x), (x, 2)) + cos(x)*Derivative(y(x), x)/sin(x),
        "kamke_2.430": (-n**2 + v*(v + 1)*sin(x)**2)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)) + cos(x)*Derivative(y(x), x)/sin(x),
        "kamke_2.431": 2*y(x) + Derivative(y(x), (x, 2)) - cos(2*x)*Derivative(y(x), x)/sin(2*x),
        "kamke_2.432": (-17*sin(x)**2/4 - 1/4)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)) + cos(x)*Derivative(y(x), x)/sin(x),
        "kamke_2.433": sin(x)*Derivative(y(x), x)/cos(x) - sqrt(cos(x)) + Derivative(y(x), (x, 2)) + (x**2*sin(x)**2/4 + x**2/2 - 6*cos(x)**2)*y(x)/(x**2*cos(x)**2),
        "kamke_2.434": Derivative(y(x), (x, 2)) + b*cos(x)*Derivative(y(x), x)/(a*sin(x)) + (c*cos(x)**2 + d*cos(x) + e)*y(x)/(a*sin(x)**2),
        "kamke_2.435": 4*y(x)*sin(3*x)/sin(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_2.436": -(n**2 - v*(v + 1)*sin(x)**2 + cos(x)**2/4 - 1/2)*y(x)/sin(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.437": -(3*sin(x)**2 + 1)*Derivative(y(x), x)/(sin(x)*cos(x)) - y(x)*sin(x)**2/cos(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_2.438": -(a*sin(x)**2*cos(x)**2 + m*(m - 1)*sin(x)**2 + n*(n - 1)*cos(x)**2)*y(x)/(sin(x)**2*cos(x)**2) + Derivative(y(x), (x, 2)),
        "kamke_2.440": Derivative(y(x), (x, 2)) - (phi(x)*Derivative(phi(x), x) - phi(x**3) + Derivative(phi(x), (x, 2)))*Derivative(y(x), x)/(phi(x)**2 + Derivative(phi(x), x)) + (-phi(x)**2*Derivative(phi(x), x) - phi(x)*Derivative(phi(x), (x, 2)) + Derivative(phi(x), x)**2)*y(x)/(phi(x)**2 + Derivative(phi(x), x)),
        "kamke_2.441": -(-(6*k**2*JacobiSN(a, k)**4 - 2*(2*k**2 + 2)*JacobiSN(a, k)**2 + 2)*y(x) + 2*JacobiCN(x, k)*JacobiDN(x, k)*JacobiSN(x, k)*Derivative(y(x), x))/(-JacobiSN(a, k) + JacobiSN(x, k)**2) + Derivative(y(x), (x, 2)),
        "kamke_2.442": x*Derivative(y(x), x)/f(x) + Derivative(y(x), (x, 2)) - y(x)/f(x),
        "kamke_2.443": Derivative(y(x), (x, 2)) + g(x)*y(x)/f(x) + Derivative(f(x), x)*Derivative(y(x), x)/(2*f(x)),
        "kamke_2.444": -a*Derivative(f(x), x)*Derivative(y(x), x)/f(x) + b*f(x)**(2*a + 1)*y(x)/f(x) + Derivative(y(x), (x, 2)),
        "kamke_2.445": (((f(x)*Derivative(g(x), (x, 2)) + 2*Derivative(f(x), x)*Derivative(g(x), x))*Derivative(f(x), x) - f(x)*Derivative(f(x), (x, 2))*Derivative(g(x), x))*(g(x)**2 - 1) - (v*(v + 1)*f(x)*Derivative(g(x), x) + 2*g(x)*Derivative(f(x), x))*f(x)*Derivative(g(x), x)**2)*y(x)/((g(x)**2 - 1)*f(x)**2*Derivative(g(x), x)) - ((f(x)*Derivative(g(x), (x, 2)) + 2*Derivative(f(x), x)*Derivative(g(x), x))*(g(x)**2 - 1) - 2*f(x)*g(x)*Derivative(g(x), x)**2)*Derivative(y(x), x)/((g(x)**2 - 1)*f(x)*Derivative(g(x), x)) + Derivative(y(x), (x, 2)),
        "kamke_2.446": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x + (x - 1)*y(x)/x**4,
        "kamke_2.447": Derivative(y(x), (x, 2)) + Derivative(y(x), x)/x + (-x - 1)*y(x)/x**4,
        "kamke_2.448": b**2*y(x)/(-a**2 + x**2)**2 + Derivative(y(x), (x, 2))
    }

    chapter_3 = {
        "kamke_3.1": -lambda_*y(x) + Derivative(y(x), (x, 3)),
        "kamke_3.2": a*x**3*y(x) - b*x + Derivative(y(x), (x, 3)),
        "kamke_3.3": -a*x**b*y(x) + Derivative(y(x), (x, 3)),
        "kamke_3.4": -4*y(x) + 3*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.5": -a**2*Derivative(y(x), x) - exp(2*a*x)*sin(x)**2 + Derivative(y(x), (x, 3)),
        "kamke_3.6": 2*a*x*Derivative(y(x), x) + a*y(x) + Derivative(y(x), (x, 3)),
        "kamke_3.7": -a*b*y(x) - x**2*Derivative(y(x), (x, 2)) + x*(a + b - 1)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.8": x**(2*c - 3)*(c - 1)*y(x) + x**(2*c - 2)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.9": b*y(x) - (3*a + 6*WeierstrassP(x, g2, g3))*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.10": (1 - n**2)*WeierstrassP(x, g2, g3)*Derivative(y(x), x) + (-a/2 + (1 - n**2)*WeierstrassPPrime(x, g2, g3)/2)*y(x) + Derivative(y(x), (x, 3)),
        "kamke_3.11": -2*n*(n + 1)*WeierstrassPPrime(x, g2, g3)*y(x) - (a + 4*n*(n + 1)*WeierstrassP(x, g2, g3))*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.12": B*WeierstrassPPrime(x, g2, g3)*y(x) + (A*WeierstrassP(x, g2, g3) + a)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.13": -(a + 3*k**2*JacobiSN(z, x)**2)*Derivative(y(x), x) + (b + c*JacobiSN(z, x)**2 - 3*k**2*JacobiCN(z, x)*JacobiDN(z, x)*JacobiSN(z, x))*y(x) + Derivative(y(x), (x, 3)),
        "kamke_3.14": b*y(x) - (a + 6*k**2*sin(x)**2)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.15": 2*f(x)*Derivative(y(x), x) + y(x)*Derivative(f(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.16": 10*y(x) - 3*Derivative(y(x), x) - 2*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_3.17": 2*a**2*y(x) - a**2*Derivative(y(x), x) - sinh(x) - 2*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_3.18": -a**3*y(x) + 3*a**2*Derivative(y(x), x) - 3*a*Derivative(y(x), (x, 2)) - exp(a*x) + Derivative(y(x), (x, 3)),
        "kamke_3.19": a0*y(x) + a1*Derivative(y(x), x) + a2*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_3.20": -8*a*x*y(x) - 6*x*Derivative(y(x), (x, 2)) + (4*a + 8*x**2 - 2)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.21": a**3*x**3*y(x) + 3*a**2*x**2*Derivative(y(x), x) + 3*a*x*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_3.22": y(x)*sin(x) - log(x) - sin(x)*Derivative(y(x), (x, 2)) - 2*cos(x)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.23": f(x)*y(x) + f(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.24": (x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + 2*y(x))*f(x) + Derivative(y(x), (x, 3)),
        "kamke_3.25": (f(x)*g(x) + Derivative(g(x), x))*y(x) + f(x)*Derivative(y(x), (x, 2)) + g(x)*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.26": (4*f(x)*g(x) + 2*Derivative(g(x), x))*y(x) + (2*f(x)**2 + 4*g(x) + Derivative(f(x), x))*Derivative(y(x), x) + 3*f(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_3.27": -3*y(x) + 18*exp(x) - 11*Derivative(y(x), x) - 8*Derivative(y(x), (x, 2)) + 4*Derivative(y(x), (x, 3)),
        "kamke_3.28": -36*n**2*WeierstrassP(x, g2, g3)*Derivative(y(x), x) - 2*n*(n + 3)*(4*n - 3)*WeierstrassPPrime(x, g2, g3)*y(x) + 27*Derivative(y(x), (x, 3)),
        "kamke_3.29": x*y(x) + x*Derivative(y(x), (x, 3)) + 3*Derivative(y(x), (x, 2)),
        "kamke_3.30": -a*x**2*y(x) + x*Derivative(y(x), (x, 3)) + 3*Derivative(y(x), (x, 2)),
        "kamke_3.31": -a*y(x) - x*Derivative(y(x), x) + x*Derivative(y(x), (x, 3)) + (a + b)*Derivative(y(x), (x, 2)),
        "kamke_3.32": x*Derivative(y(x), (x, 3)) - (2*v + x)*Derivative(y(x), (x, 2)) + (x - 1)*y(x) - (-2*v + x - 1)*Derivative(y(x), x),
        "kamke_3.33": 4*x*Derivative(y(x), x) + x*Derivative(y(x), (x, 3)) + (x**2 - 3)*Derivative(y(x), (x, 2)) - f(x) + 2*y(x),
        "kamke_3.34": a*x*y(x) - b + 2*x*Derivative(y(x), (x, 3)) + 3*Derivative(y(x), (x, 2)),
        "kamke_3.35": 2*x*Derivative(y(x), (x, 3)) + (1 - 2*nu)*y(x) - (4*nu + 4*x - 4)*Derivative(y(x), (x, 2)) + (6*nu + 2*x - 5)*Derivative(y(x), x),
        "kamke_3.36": 2*x*Derivative(y(x), (x, 3)) + (6*a*k + 6*b*x)*Derivative(y(x), x) + (6*a*x + 3*k)*Derivative(y(x), (x, 2)) + (3*b*k + 2*c*x)*y(x),
        "kamke_3.37": -x*(x - 2)*Derivative(y(x), (x, 2)) + x*(x - 2)*Derivative(y(x), (x, 3)) + 2*y(x) - 2*Derivative(y(x), x),
        "kamke_3.38": -8*x*Derivative(y(x), x) + (2*x - 1)*Derivative(y(x), (x, 3)) + 8*y(x),
        "kamke_3.39": (x + 4)*Derivative(y(x), (x, 2)) + (2*x - 1)*Derivative(y(x), (x, 3)) + 2*Derivative(y(x), x),
        "kamke_3.40": a*x**2*y(x) + x**2*Derivative(y(x), (x, 3)) - 6*Derivative(y(x), x),
        "kamke_3.41": x**2*Derivative(y(x), (x, 3)) + (x + 1)*Derivative(y(x), (x, 2)) - y(x),
        "kamke_3.42": x**2*Derivative(y(x), (x, 3)) - x*Derivative(y(x), (x, 2)) + (x**2 + 1)*Derivative(y(x), x),
        "kamke_3.43": -4*a**3*x**(2*a - 1)*y(x) + x**2*Derivative(y(x), (x, 3)) + 3*x*Derivative(y(x), (x, 2)) + (-4*a**2*nu**2 + 4*a**2*x**(2*a) + 1)*Derivative(y(x), x),
        "kamke_3.44": -2*n*(-2*m + 2*x + 1)*y(x) + x**2*Derivative(y(x), (x, 3)) - x*(-3*m + 3*x)*Derivative(y(x), (x, 2)) + (m*(2*m - 1) + 2*x**2 + x*(-4*m + 4*n))*Derivative(y(x), x),
        "kamke_3.45": x**2*Derivative(y(x), (x, 3)) + 3*x*y(x) + 4*x*Derivative(y(x), (x, 2)) + (x**2 + 2)*Derivative(y(x), x) - f(x),
        "kamke_3.46": x**2*Derivative(y(x), (x, 3)) + 5*x*Derivative(y(x), (x, 2)) - log(x) + 4*Derivative(y(x), x),
        "kamke_3.47": x**2*Derivative(y(x), (x, 3)) + 6*x*Derivative(y(x), (x, 2)) + 6*Derivative(y(x), x),
        "kamke_3.48": a*x**2*y(x) + x**2*Derivative(y(x), (x, 3)) + 6*x*Derivative(y(x), (x, 2)) + 6*Derivative(y(x), x),
        "kamke_3.49": 3*ps*(3*qs + 1)*Derivative(y(x), x) - x**2*y(x) + x**2*Derivative(y(x), (x, 3)) - x*(3*ps + 3*qs)*Derivative(y(x), (x, 2)),
        "kamke_3.50": -2*a*x*y(x) + x**2*Derivative(y(x), (x, 3)) - x*(2*n + 2)*Derivative(y(x), (x, 2)) + (a*x**2 + 6*n)*Derivative(y(x), x),
        "kamke_3.51": x**2*Derivative(y(x), (x, 3)) - (x**2 - 2*x)*Derivative(y(x), (x, 2)) - (nu**2 + x**2 - 1/4)*Derivative(y(x), x) + (nu**2 + x**2 - 2*x - 1/4)*y(x),
        "kamke_3.52": -nu*(x + 1)*y(x) + nu*(2*x + 1)*Derivative(y(x), x) + x**2*Derivative(y(x), (x, 3)) - x*(nu + x)*Derivative(y(x), (x, 2)),
        "kamke_3.53": x**2*Derivative(y(x), (x, 3)) + (nu**2 - 1/4)*y(x) - (2*x**2 - 2*x)*Derivative(y(x), (x, 2)) + (-nu**2 + x**2 - 2*x + 1/4)*Derivative(y(x), x),
        "kamke_3.54": 2*x**2*y(x) + x**2*Derivative(y(x), (x, 3)) - (2*x**3 - 6)*Derivative(y(x), x) - (x**4 - 6*x)*Derivative(y(x), (x, 2)),
        "kamke_3.55": 8*x*Derivative(y(x), (x, 2)) + (x**2 + 1)*Derivative(y(x), (x, 3)) - 2*log(x) + 10*Derivative(y(x), x) - 3 + x**(-2),
        "kamke_3.56": -2*x*y(x) - 2*x*Derivative(y(x), (x, 2)) + (x**2 + 2)*Derivative(y(x), x) + (x**2 + 2)*Derivative(y(x), (x, 3)),
        "kamke_3.57": a*y(x) + 2*x*(x - 1)*Derivative(y(x), (x, 3)) + (6*x - 3)*Derivative(y(x), (x, 2)) + (2*a*x + b)*Derivative(y(x), x),
        "kamke_3.58": 4*x**2*Derivative(y(x), (x, 3)) + (4*x + 4)*Derivative(y(x), x) + (x**2 + 14*x - 1)*Derivative(y(x), (x, 2)) + 2*y(x),
        "kamke_3.59": x*(a*x + b)*Derivative(y(x), (x, 3)) + x*Derivative(y(x), x) + (alpha*x + bbeta)*Derivative(y(x), (x, 2)) - f(x) + y(x),
        "kamke_3.60": x**3*Derivative(y(x), (x, 3)) + x*(1 - nu**2)*Derivative(y(x), x) + (a*x**3 + nu**2 - 1)*y(x),
        "kamke_3.61": x**3*Derivative(y(x), (x, 3)) + (4*nu**2 - 1)*y(x) + (4*x**3 + x*(1 - 4*nu**2))*Derivative(y(x), x),
        "kamke_3.62": x**3*Derivative(y(x), (x, 3)) + x*(a*x**(2*nu) - nu**2 + 1)*Derivative(y(x), x) + (a*x**(2*nu)*(nu - 1) + b*x**(3*nu) + nu**2 - 1)*y(x),
        "kamke_3.63": -6*x**3*(x - 1)*log(x) + x**3*(x + 8) + x**3*Derivative(y(x), (x, 3)) + 3*x**2*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x) + 2*y(x),
        "kamke_3.64": x**3*Derivative(y(x), (x, 3)) + 3*x**2*Derivative(y(x), (x, 2)) + x*(1 - a**2)*Derivative(y(x), x),
        "kamke_3.65": x**3*Derivative(y(x), (x, 3)) - 4*x**2*Derivative(y(x), (x, 2)) + x*(x**2 + 8)*Derivative(y(x), x) - (2*x**2 + 8)*y(x),
        "kamke_3.66": x**3*Derivative(y(x), (x, 3)) + 6*x**2*Derivative(y(x), (x, 2)) + (a*x**3 - 12)*y(x),
        "kamke_3.67": x**3*Derivative(y(x), (x, 3)) + x**2*(3 - 3*a)*Derivative(y(x), (x, 2)) + (a*(-a**2 + 4*c**2*nu**2) + 4*b**2*c**2*x**(2*c)*(-a + c))*y(x) + (3*a*x*(a - 1) + 4*b**2*c**2*x**(2*c + 1) - 4*c**2*nu**2 + 1)*Derivative(y(x), x),
        "kamke_3.68": x**3*Derivative(y(x), (x, 3)) + x**2*(x + 3)*Derivative(y(x), (x, 2)) + x*(5*x - 30)*Derivative(y(x), x) + (4*x + 30)*y(x),
        "kamke_3.69": x**3*Derivative(y(x), (x, 3)) - 2*x**3 + x**2*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x) - y(x) + log(x),
        "kamke_3.70": x*(x**2 + 1)*Derivative(y(x), (x, 3)) + (6*x**2 + 3)*Derivative(y(x), (x, 2)) - 12*y(x),
        "kamke_3.71": x**2*(x + 3)*Derivative(y(x), (x, 3)) - x*(3*x + 6)*Derivative(y(x), (x, 2)) + (6*x + 6)*Derivative(y(x), x) - 6*y(x),
        "kamke_3.72": -n*(n + 1)*y(x) + (-2*a1 + 2*x)*(-a2 + x)*(-a3 + x)*Derivative(y(x), (x, 3)) - (2*b + 2*x*(n**2 + n - 3))*Derivative(y(x), x) + (3*a1*a2 + 3*a1*a3 + 3*a2*a3 + 9*x**2 - x*(6*a1 + 6*a2 + 6*a3))*Derivative(y(x), (x, 2)),
        "kamke_3.73": x**3*(x + 1)*Derivative(y(x), (x, 3)) - x**2*(4*x + 2)*Derivative(y(x), (x, 2)) + x*(10*x + 4)*Derivative(y(x), x) - (12*x + 4)*y(x),
        "kamke_3.74": 4*x**4*Derivative(y(x), (x, 3)) - 4*x**3*Derivative(y(x), (x, 2)) + 4*x**2*Derivative(y(x), x) - 1,
        "kamke_3.75": x**3*(x**2 + 1)*Derivative(y(x), (x, 3)) - x**2*(4*x**2 + 2)*Derivative(y(x), (x, 2)) + x*(10*x**2 + 4)*Derivative(y(x), x) - (12*x**2 + 4)*y(x),
        "kamke_3.76": x**6*Derivative(y(x), (x, 3)) + x**2*Derivative(y(x), (x, 2)) - 2*y(x),
        "kamke_3.77": a*y(x) + x**6*Derivative(y(x), (x, 3)) + 6*x**5*Derivative(y(x), (x, 2)),
        "kamke_3.78": x**2*(x**4 + 2*x**2 + 2*x + 1)*Derivative(y(x), (x, 3)) + (x**4 + 4*x**3 + 8*x**2 + 6*x + 1)*y(x) + (x**6 - 6*x**3 - 15*x**2 - 12*x - 2)*Derivative(y(x), x) - (2*x**6 + 3*x**4 - 6*x**2 - 6*x - 1)*Derivative(y(x), (x, 2)),
        "kamke_3.79": -c*y(x) + (-a + x)**3*(-b + x)**3*Derivative(y(x), (x, 3)),
        "kamke_3.80": (2*cos(x) + 1)*Derivative(y(x), (x, 2)) - sin(x)*Derivative(y(x), x) + sin(x)*Derivative(y(x), (x, 3)) - cos(x),
        "kamke_3.81": (x + sin(x))*Derivative(y(x), (x, 3)) + (3*cos(x) + 3)*Derivative(y(x), (x, 2)) - y(x)*cos(x) - 3*sin(x)*Derivative(y(x), x) + sin(x),
        "kamke_3.82": 2*nu*(nu + 1)*y(x)*sin(2*x) + (4*nu*(nu + 1)*sin(x)**2 + cos(2*x))*Derivative(y(x), x) + sin(x)**2*Derivative(y(x), (x, 3)) + 3*sin(x)*cos(x)*Derivative(y(x), (x, 2)),
        "kamke_3.83": (f(x)*Derivative(y(x), (x, 2)) + g(x)*Derivative(y(x), x) + h(x)*y(x))*Af(x) + f(x)*Derivative(y(x), (x, 3)) + g(x)*Derivative(y(x), (x, 2)) + h(x)*Derivative(y(x), x) + y(x)*Derivative(h(x), x) + Derivative(f(x), x)*Derivative(y(x), (x, 2)) + Derivative(g(x), x)*Derivative(y(x), x),
        "kamke_3.84": n*y(x) + x*Derivative(y(x), x) + Derivative(y(x), (x, 3)),
        "kamke_3.85": -n*y(x) - x*Derivative(y(x), x) + Derivative(y(x), (x, 3))
    }

    chapter_4 = {
        "kamke_4.1": Derivative(y(x), (x, 4)),
        "kamke_4.2": -fs + 4*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.3": lambda_*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.4": -16*x**4*exp(x**2) + 12*y(x) - 12*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 4)),
        "kamke_4.5": a**4*y(x) + 2*a**2*Derivative(y(x), (x, 2)) - cosh(a*x) + Derivative(y(x), (x, 4)),
        "kamke_4.6": a**4*lambda_*y(x) + a**2*(lambda_ + 1)*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 4)),
        "kamke_4.7": a*b*Derivative(y(x), x) + a*(b*x - 1)*Derivative(y(x), (x, 2)) + lambda_*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.8": (a*x**2 + b*lambda_ + c)*Derivative(y(x), (x, 2)) + (a*x**2 + bbeta*lambda_ + ggamma)*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.9": a*WeierstrassP(x, g2, g3)*Derivative(y(x), (x, 2)) + b*WeierstrassPPrime(x, g2, g3)*Derivative(y(x), x) + (c*Derivative(WeierstrassP(x, g2, g3), (x, 2)) + d)*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.10": b*Derivative(y(x), x) - (a + 12*k**2*JacobiSN(z, x)**2)*Derivative(y(x), (x, 2)) + (alpha*JacobiSN(z, x)**2 + bbeta)*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.11": 10*df*Derivative(y(x), x) + 10*fs*Derivative(y(x), (x, 2)) + (3*ddf + 3*fs**2)*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.12": 4*y(x) - 32*sin(2*x) + 24*cos(2*x) - 4*Derivative(y(x), x) - 3*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), (x, 3)) + Derivative(y(x), (x, 4)),
        "kamke_4.13": a**4*x**4*y(x) + 4*a**3*x**3*Derivative(y(x), x) + 6*a**2*x**2*Derivative(y(x), (x, 2)) + 4*a*x*Derivative(y(x), (x, 3)) + Derivative(y(x), (x, 4)),
        "kamke_4.14": 6*fs*Derivative(y(x), (x, 3)) + (4*df + 11*fs**2 + 10*gs)*Derivative(y(x), (x, 2)) + (ddf + 7*df*fs + 10*dg + 6*fs**3 + 30*fs*gs)*Derivative(y(x), x) + (3*ddg + 6*df*gs + 15*dg*fs + 18*fs**2*gs + 9*gs**2)*y(x) + Derivative(y(x), (x, 4)),
        "kamke_4.15": -4*cos(x) - 3*Derivative(y(x), x) + 11*Derivative(y(x), (x, 2)) - 12*Derivative(y(x), (x, 3)) + 4*Derivative(y(x), (x, 4)),
        "kamke_4.16": x*Derivative(y(x), (x, 4)) + 5*Derivative(y(x), (x, 3)) - 24,
        "kamke_4.17": x**3*(2*x**2 - 6)*y(x) + 12*x**3*Derivative(y(x), (x, 2)) - x**2*(9*x**2 - 7)*Derivative(y(x), x) + x*Derivative(y(x), (x, 4)) - (6*x**2 + 1)*Derivative(y(x), (x, 3)),
        "kamke_4.18": nu**2*(nu**2*x**2 + 4)*y(x) + x**2*Derivative(y(x), (x, 4)) - (2*nu**2*x**2 + 12)*Derivative(y(x), (x, 2)),
        "kamke_4.19": a*y(x) - b*x**2 + x**2*Derivative(y(x), (x, 4)) + 2*x*Derivative(y(x), (x, 3)),
        "kamke_4.20": x**2*Derivative(y(x), (x, 4)) + 4*x*Derivative(y(x), (x, 3)) + 2*Derivative(y(x), (x, 2)),
        "kamke_4.21": x**2*Derivative(y(x), (x, 4)) + 6*x*Derivative(y(x), (x, 3)) + 6*Derivative(y(x), (x, 2)),
        "kamke_4.22": -lambda_**2*y(x) + x**2*Derivative(y(x), (x, 4)) + 6*x*Derivative(y(x), (x, 3)) + 6*Derivative(y(x), (x, 2)),
        "kamke_4.23": x**2*Derivative(y(x), (x, 4)) + 8*x*Derivative(y(x), (x, 3)) + 12*Derivative(y(x), (x, 2)),
        "kamke_4.24": -lambda_**2*y(x) + x**2*Derivative(y(x), (x, 4)) + 8*x*Derivative(y(x), (x, 3)) + 12*Derivative(y(x), (x, 2)),
        "kamke_4.25": -b**4*y(x)/16 + x**2*Derivative(y(x), (x, 4)) + x*(2*n - 2*nu + 4)*Derivative(y(x), (x, 3)) + (n - nu + 1)*(n - nu + 2)*Derivative(y(x), (x, 2)),
        "kamke_4.26": -a**4*x**3*y(x) + x**3*Derivative(y(x), (x, 4)) + 2*x**2*Derivative(y(x), (x, 3)) - x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_4.27": x**3*Derivative(y(x), (x, 4)) + 6*x**2*Derivative(y(x), (x, 3)) + 6*x*Derivative(y(x), (x, 2)),
        "kamke_4.28": -2*n*x**2*(n + 1)*Derivative(y(x), (x, 2)) + 4*n*x*(n + 1)*Derivative(y(x), x) + x**4*Derivative(y(x), (x, 4)) + (a*x**4 + n*(n - 2)*(n + 1)*(n + 3))*y(x),
        "kamke_4.29": -4*x**4*y(x) + x**4*Derivative(y(x), (x, 4)) + 4*x**3*Derivative(y(x), (x, 3)) - x**2*(4*n**2 - 1)*Derivative(y(x), (x, 2)) + x*(4*n**2 - 1)*Derivative(y(x), x),
        "kamke_4.30": x**4*Derivative(y(x), (x, 4)) + 4*x**3*Derivative(y(x), (x, 3)) - x**2*(4*n**2 - 1)*Derivative(y(x), (x, 2)) - x*(4*n**2 - 1)*Derivative(y(x), x) + (4*n**2 - 4*x**4 - 1)*y(x),
        "kamke_4.31": x**4*Derivative(y(x), (x, 4)) + 4*x**3*Derivative(y(x), (x, 3)) - x**2*(4*n**2 + 3)*Derivative(y(x), (x, 2)) + x*(12*n**2 - 3)*Derivative(y(x), x) - (12*n**2 + 4*x**4 - 3)*y(x),
        "kamke_4.32": x**4*Derivative(y(x), (x, 4)) + 6*x**3*Derivative(y(x), (x, 3)) + (16*x**3 + x*(-rho**2 - sigma**2 + 1))*Derivative(y(x), x) + (4*x**4 + x**2*(-rho**2 - sigma**2 + 7))*Derivative(y(x), (x, 2)) + (rho**2*sigma**2 + 8*x**2)*y(x),
        "kamke_4.33": x**4*Derivative(y(x), (x, 4)) + 6*x**3*Derivative(y(x), (x, 3)) + (8*x**2 + (mu**2 - nu**2)**2)*y(x) + (16*x**3 + x*(-2*mu**2 - 2*nu**2 + 1))*Derivative(y(x), x) + (4*x**4 + x**2*(-2*mu**2 - 2*nu**2 + 7))*Derivative(y(x), (x, 2)),
        "kamke_4.34": x**4*Derivative(y(x), (x, 4)) + 8*x**3*Derivative(y(x), (x, 3)) + 12*x**2*Derivative(y(x), (x, 2)),
        "kamke_4.35": a*y(x) + x**4*Derivative(y(x), (x, 4)) + 8*x**3*Derivative(y(x), (x, 3)) + 12*x**2*Derivative(y(x), (x, 2)),
        "kamke_4.36": x**4*Derivative(y(x), (x, 4)) + x**3*(6 - 4*a)*Derivative(y(x), (x, 3)) + x**2*(AAA + 4*b**2*c**2*x**(2*c))*Derivative(y(x), (x, 2)) + x*(4*BBB*b**2*c**2*x**(2*c) + CCC*(2*a - 1))*Derivative(y(x), x) + (4*DDD*b**2*c**2*x**(2*c) + EEE)*y(x),
        "kamke_4.37": x**4*Derivative(y(x), (x, 4)) + x**3*(-4*a - 4*c + 6)*Derivative(y(x), (x, 3)) + x**2*(2*a**2 - 2*c**2*nu**2 + (4*a - 4)*(c - 1) + 4*(a + c - 1)**2 - 1)*Derivative(y(x), (x, 2)) + x*(2*a + 2*c - 1)*(-2*a**2 + 2*c**2*nu**2 - (2*a - 1)*(2*c - 1))*Derivative(y(x), x) + (-b**4*c**4*x**(4*c) + (a**2 - c**2*nu**2)*(a**2 + 4*a*c - c**2*nu**2 + 4*c**2))*y(x),
        "kamke_4.38": -b**4*x**(2/nu)*y(x)/16 + nu**4*x**4*Derivative(y(x), (x, 4)) + nu**3*x**3*(4*nu - 2)*Derivative(y(x), (x, 3)) + nu**2*x**2*(nu - 1)*(2*nu - 1)*Derivative(y(x), (x, 2)),
        "kamke_4.39": 10*x*(x**2 - 1)*Derivative(y(x), (x, 3)) - 6*x*(mu*(mu + 1) + nu*(nu + 1) - 2)*Derivative(y(x), x) + (x**2 - 1)**2*Derivative(y(x), (x, 4)) + (24*x**2 - (x**2 - 1)*(2*mu*(mu + 1) + 2*nu*(nu + 1)) - 8)*Derivative(y(x), (x, 2)) + (-2*mu*(mu + 1) - 2*nu*(nu + 1) + (mu*(mu + 1) - nu*(nu + 1))**2)*y(x),
        "kamke_4.40": (2*x + exp(x))*Derivative(y(x), (x, 4)) + (4*exp(x) + 8)*Derivative(y(x), (x, 3)) + y(x)*exp(x) + 4*exp(x)*Derivative(y(x), x) + 6*exp(x)*Derivative(y(x), (x, 2)) - 1/x**5,
        "kamke_4.41": (a**4*sin(x)**4 - 3)*y(x) + (sin(x)**2 - 3)*sin(x)**2*Derivative(y(x), (x, 2)) + (2*sin(x)**2 + 3)*sin(x)*cos(x)*Derivative(y(x), x) + sin(x)**4*Derivative(y(x), (x, 4)) + 2*sin(x)**3*cos(x)*Derivative(y(x), (x, 3)),
        "kamke_4.42": -fs + y(x)*sin(x)**6 - 6*sin(x)**6*Derivative(y(x), (x, 2)) + sin(x)**6*Derivative(y(x), (x, 4)) - 4*sin(x)**5*cos(x)*Derivative(y(x), x) + 4*sin(x)**5*cos(x)*Derivative(y(x), (x, 3)),
        "kamke_4.43": 2*df*(-a**2*Derivative(y(x), x) + Derivative(y(x), (x, 3))) + fs*(a**4*y(x) - 2*a**2*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 4))),
        "kamke_4.44": fs*Derivative(y(x), (x, 4))
    }

    n = 5
    chapter_5 = {
        "kamke_5.1": a**4*y(x) - 2*a**2*Derivative(y(x), (x, 2)) - lambda_*(a*x - b)*(-a**2*y(x) + Derivative(y(x), (x, 2))) + Derivative(y(x), (x, 4)),
        "kamke_5.2": -a*x - b*sin(x) - c*cos(x) + Derivative(y(x), x) + 2*Derivative(y(x), (x, 3)) + Derivative(y(x), (x, 5)),
        "kamke_5.3": y(x) - sin(x/2)*sin(3*x/2) + Derivative(y(x), (x, 6)),
        "kamke_5.4": -a*x*y(x) - b + Derivative(y(x), (x, n)),
        "kamke_5.5": a*nu*x**(nu - 1)*y(x) + a*x**nu*Derivative(y(x), x) + Derivative(y(x), (x, n)),
        "kamke_5.6": a*Derivative(y(x), (x, n - 1)) - fs + Derivative(y(x), (x, n)),
        "kamke_5.7": a*x*y(x) - m*n*Derivative(y(x), (x, n - 1)) + x*Derivative(y(x), (x, n)),
        "kamke_5.8": x*(a*Derivative(y(x), x) + b*Derivative(y(x), (x, 2)) + c*Derivative(y(x), (x, 3)) + e*Derivative(y(x), (x, 4)))*y(x),
        "kamke_5.10": -a*y(x) + x**n*Derivative(y(x), (x, 2*n)),
        "kamke_5.11": -a*y(x) + x**(2*n)*Derivative(y(x), (x, n)),
        "kamke_5.12": -a*y(x) + x**(n + 1/2)*Derivative(y(x), (x, 2*n + 1)),
        "kamke_5.13": -c*y(x) + (-a + x)**n*(-b + x)**n*Derivative(y(x), (x, n))
    }

    n = symbols('n')
    chapter_6 = {
        "kamke_6.1": -y(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_6.2": -6*y(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_6.3": -x - 6*y(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_6.4": -6*y(x)**2 + 4*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.5": a*y(x)**2 + b*x + c + Derivative(y(x), (x, 2)),
        "kamke_6.6": a - x*y(x) - 2*y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.7": -a*y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.8": -2*a**2*y(x)**3 + 2*a*b*x*y(x) - b + Derivative(y(x), (x, 2)),
        "kamke_6.9": a*y(x)**3 + b*x*y(x) + c*y(x) + d + Derivative(y(x), (x, 2)),
        "kamke_6.10": a*y(x)**3 + b*y(x)**2 + c*y(x) + d + Derivative(y(x), (x, 2)),
        "kamke_6.11": a*x**r*y(x)**n + Derivative(y(x), (x, 2)),
        "kamke_6.12": a**(2*n)*(n + 1)*y(x)**(2*n + 1) - y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.13": Derivative(y(x), (x, 2)) - 1/(a*y(x)**2 + b*x*y(x) + c*x**2 + d*y(x) + e*x + k)**(3/2),
        "kamke_6.14": -exp(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.15": a*sqrt(y(x))*exp(x) + Derivative(y(x), (x, 2)),
        "kamke_6.16": exp(x)*sin(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.17": a*sin(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.18": a**2*sin(y(x)) - b*sin(x) + Derivative(y(x), (x, 2)),
        "kamke_6.19": a**2*sin(y(x)) - b*f(x) + Derivative(y(x), (x, 2)),
        "kamke_6.20": Derivative(y(x), (x, 2)) - h(y(x)/sqrt(x))/x**(3/2),
        "kamke_6.21": -y(x)**2 - 2*y(x) - 3*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.22": -y(x)**(3/2) + 12*y(x) - 7*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.23": 6*a**2*y(x) + 5*a*Derivative(y(x), x) - 6*y(x)**2 + Derivative(y(x), (x, 2)),
        "kamke_6.24": 2*a**2*y(x) + 3*a*Derivative(y(x), x) - 2*y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.25": Derivative(y(x), (x, 2)) - (3*n + 4)*Derivative(y(x), x)/n - (n + 2)*(2*n + 2)*(y(x)**(n/(n + 1)) - 1)*y(x)/n**2,
        "kamke_6.26": a*Derivative(y(x), x) + b*y(x)**n + (a**2/4 - 1/4)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.27": a*Derivative(y(x), x) + b*x**r*y(x)**n + Derivative(y(x), (x, 2)),
        "kamke_6.28": a*Derivative(y(x), x) - 2*a + b*exp(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.29": a*Derivative(y(x), x) + f(x)*sin(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.30": -y(x)**3 + y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.31": a*y(x) - y(x)**3 + y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.32": 2*a**2*y(x) + a*y(x)**2 + (3*a + y(x))*Derivative(y(x), x) - y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.33": (3*f(x) + y(x))*Derivative(y(x), x) + (2*f(x)**2 + Derivative(f(x), x))*y(x) + f(x)*y(x)**2 - y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.34": b*f(x)**3 - (f(x) + Derivative(f(x), x)/f(x))*(y(x)**2 + 3*Derivative(y(x), x)) + (a*f(x)**2 + 3*Derivative(f(x), x) - Derivative(f(x), (x, 2))/f(x) + 3*Derivative(f(x), x)**2/f(x)**2)*y(x) - y(x)**3 + y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.35": (y(x) - 3*Derivative(f(x), x)/(2*f(x)))*Derivative(y(x), x) + (f(x) - Derivative(f(x), (x, 2)) + Derivative(f(x), x)**2/f(x)**2)*y(x)/(2*f(x)) - y(x)**3 + Derivative(y(x), (x, 2)) - y(x)**2*Derivative(f(x), x)/(2*f(x)),
        "kamke_6.36": f(x)*Derivative(y(x), x) + y(x)*Derivative(f(x), x) + 2*y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.37": (y(x)**2 + Derivative(y(x), x))*f(x) - g(x) + 2*y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.38": f(x)*y(x) - g(x) + y(x)**3 + 3*y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.39": (f(x) + 3*y(x))*Derivative(y(x), x) + f(x)*y(x)**2 + y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.40": -4*a**2*y(x) - 3*a*y(x)**2 - b - 3*y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.41": -(f(x) + 3*y(x))*Derivative(y(x), x) + f(x)*y(x)**2 + y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.42": -2*a*y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.43": a*y(x)*Derivative(y(x), x) + b*y(x)**3 + Derivative(y(x), (x, 2)),
        "kamke_6.44": h(x, y(x))*Derivative(y(x), x) + j(x, y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.45": a*Derivative(y(x), x)**2 + b*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.46": a*Abs(Derivative(y(x), x))*Derivative(y(x), x) + b*Derivative(y(x), x) + c*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.47": a*Derivative(y(x), x)**2 + b*Derivative(y(x), x) + c*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.48": a*Derivative(y(x), x)**2 + b*sin(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.49": a*Abs(Derivative(y(x), x))*Derivative(y(x), x) + b*sin(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.50": a*y(x)*Derivative(y(x), x)**2 + b*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.51": g(x)*Derivative(y(x), x) + h(y(x))*Derivative(y(x), x)**2 + Derivative(y(x), (x, 2)),
        "kamke_6.52": f(x)*h(y(x)) + g(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2)) - j(y(x))*Derivative(y(x), x)**2/h(y(x)),
        "kamke_6.53": (1 - Derivative(j(y(x)), y(x)))*Derivative(y(x), x)**2/j(y(x)) + f(x)*Derivative(y(x), x) + g(x)*j(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.54": h(y(x))*Derivative(y(x), x)**2 + j(y(x))*Derivative(y(x), x) + kf(y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.55": (h(x, y(x))*Derivative(y(x), x) + j(x, y(x)))*(Derivative(y(x), x)**2 + 1) + Derivative(y(x), (x, 2)),
        "kamke_6.56": a*(Derivative(y(x), x)**2 + 1)**2*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.57": -a*(x*Derivative(y(x), x) - y(x))**r + Derivative(y(x), (x, 2)),
        "kamke_6.58": -k*x**a*y(x)**b*Derivative(y(x), x)**c + Derivative(y(x), (x, 2)),
        "kamke_6.59": (Derivative(y(x), x) - y(x)/x)**a*h(x, y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.60": -a*sqrt(Derivative(y(x), x)**2 + 1) + Derivative(y(x), (x, 2)),
        "kamke_6.61": -a*sqrt(Derivative(y(x), x)**2 + 1) - b + Derivative(y(x), (x, 2)),
        "kamke_6.62": -a*sqrt(b*y(x)**2 + Derivative(y(x), x)**2) + Derivative(y(x), (x, 2)),
        "kamke_6.63": -a*(Derivative(y(x), x)**2 + 1)**(3/2) + Derivative(y(x), (x, 2)),
        "kamke_6.64": -2*a*x*(Derivative(y(x), x)**2 + 1)**(3/2) + Derivative(y(x), (x, 2)),
        "kamke_6.65": -a*(Derivative(y(x), x)**2 + 1)**(3/2)*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.66": -a*(Derivative(y(x), x)**2 + 1)**(3/2)*(b*x + c + y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.67": -sqrt(y(x)**4 + 4*Derivative(y(x), x))*y(x)*Derivative(y(x), x) + y(x)**3*Derivative(y(x), x) + Derivative(y(x), (x, 2)),
        "kamke_6.68": -h(Derivative(y(x), x), a*x + b*y(x)) + Derivative(y(x), (x, 2)),
        "kamke_6.69": -h(x, Derivative(y(x), x)/y(x))*y(x) + Derivative(y(x), (x, 2)),
        "kamke_6.70": -x**(n - 2)*h(y(x)/x**n, x**(1 - n)*Derivative(y(x), x)) + Derivative(y(x), (x, 2)),
        "kamke_6.71": 9*Derivative(y(x), x)**4 + 8*Derivative(y(x), (x, 2)),
        "kamke_6.72": a*Derivative(y(x), (x, 2)) + c*y(x) + h(Derivative(y(x), x)),
        "kamke_6.73": -x*y(x)**n + x*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x),
        "kamke_6.74": a*x**m*y(x)**n + x*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x),
        "kamke_6.75": x*exp(y(x)) + x*Derivative(y(x), (x, 2)) + 2*Derivative(y(x), x),
        "kamke_6.76": a*Derivative(y(x), x) + b*x*exp(y(x)) + x*Derivative(y(x), (x, 2)),
        "kamke_6.77": a*Derivative(y(x), x) + b*x**(5 - 2*a)*exp(y(x)) + x*Derivative(y(x), (x, 2)),
        "kamke_6.78": x*Derivative(y(x), (x, 2)) - (1 - y(x))*Derivative(y(x), x),
        "kamke_6.79": -x**2*Derivative(y(x), x)**2 + x*Derivative(y(x), (x, 2)) + y(x)**2 + 2*Derivative(y(x), x),
        "kamke_6.80": a*(x*Derivative(y(x), x) - y(x))**2 - b + x*Derivative(y(x), (x, 2)),
        "kamke_6.81": 2*x*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**3 + Derivative(y(x), x),
        "kamke_6.82": -a*(-y(x) + y(x)**n) + x**2*Derivative(y(x), (x, 2)),
        "kamke_6.83": a*(exp(y(x)) - 1) + x**2*Derivative(y(x), (x, 2)),
        "kamke_6.84": x**2*Derivative(y(x), (x, 2)) - x*(2*a + b - 1)*Derivative(y(x), x) + (a*(a + b) + b**2*c**2*x**(2*b))*y(x),
        "kamke_6.85": x**2*Derivative(y(x), (x, 2)) + x*(a + 1)*Derivative(y(x), x) - x**k*h(x**k*y(x), k*y(x) + x*Derivative(y(x), x)),
        "kamke_6.86": a*(x*Derivative(y(x), x) - y(x))**2 - b*x**2 + x**2*Derivative(y(x), (x, 2)),
        "kamke_6.87": a*y(x)*Derivative(y(x), x)**2 + b*x + x**2*Derivative(y(x), (x, 2)),
        "kamke_6.88": x**2*Derivative(y(x), (x, 2)) - sqrt(a*x**2*Derivative(y(x), x)**2 + b*y(x)**2),
        "kamke_6.89": (x**2 + 1)*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 + 1,
        "kamke_6.90": -x**4*Derivative(y(x), x)**2 + 4*x**2*Derivative(y(x), (x, 2)) + 4*y(x),
        "kamke_6.91": a*y(x)**3 + 9*x**2*Derivative(y(x), (x, 2)) + 2*y(x),
        "kamke_6.92": x**3*(-y(x)**3 + y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2))) + 12*x*y(x) + 24,
        "kamke_6.93": -a*(x*Derivative(y(x), x) - y(x))**2 + x**3*Derivative(y(x), (x, 2)),
        "kamke_6.94": b + 2*x**3*Derivative(y(x), (x, 2)) + x**2*(2*x*y(x) + 9)*Derivative(y(x), x) + x*(a - 2*x**2*y(x)**2 + 3*x*y(x))*y(x),
        "kamke_6.95": a*x*y(x) + b + (8*x**3 - 2*x**k)*(-y(x)**3 + y(x)*Derivative(y(x), x) + Derivative(y(x), (x, 2))) - (k*x**(k - 1) - 12*x**2)*(y(x)**2 + 3*Derivative(y(x), x)),
        "kamke_6.96": a**2*y(x)**n + x**4*Derivative(y(x), (x, 2)),
        "kamke_6.97": x**4*Derivative(y(x), (x, 2)) - x*(x**2 + 2*y(x))*Derivative(y(x), x) + 4*y(x)**2,
        "kamke_6.98": x**4*Derivative(y(x), (x, 2)) - x**2*(x + Derivative(y(x), x))*Derivative(y(x), x) + 4*y(x)**2,
        "kamke_6.99": x**4*Derivative(y(x), (x, 2)) + (x*Derivative(y(x), x) - y(x))**3,
        "kamke_6.100": sqrt(x)*Derivative(y(x), (x, 2)) - y(x)**(3/2),
        "kamke_6.101": (a*x**2 + b*x + c)**(3/2)*Derivative(y(x), (x, 2)) - F(y(x)/sqrt(a*x**2 + b*x + c)),
        "kamke_6.102": x**(n/(n + 1))*Derivative(y(x), (x, 2)) - y(x)**((2*n + 1)/(n + 1)),
        "kamke_6.103": f(x)**2*Derivative(y(x), (x, 2)) + f(x)*Derivative(f(x), x)*Derivative(y(x), x) - h(y(x), f(x)*Derivative(y(x), x)),
        "kamke_6.104": -a + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.105": -a*x + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.106": -a*x**2 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.107": -a + y(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2,
        "kamke_6.108": -a*x - b + y(x)**2 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.109": y(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 - Derivative(y(x), x),
        "kamke_6.110": y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 + 1,
        "kamke_6.111": y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 - 1,
        "kamke_6.112": (a*y(x)**4 + b)*exp(2*x) + (c*y(x)**2 + d)*y(x)*exp(x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.113": -y(x)**2*log(y(x)) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.114": (Derivative(f(x), (x, 2))/f(x) - Derivative(f(x), x)**2/f(x)**2)*y(x)**2 + f(x)*y(x)**3 + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 - Derivative(y(x), x),
        "kamke_6.115": f(x)*Derivative(y(x), x) - y(x)**3 - y(x)*Derivative(f(x), x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.116": f(x)*y(x)**3 - y(x)**4 - y(x)*Derivative(f(x), (x, 2)) + y(x)*Derivative(y(x), (x, 2)) + Derivative(f(x), x)*Derivative(y(x), x) - Derivative(y(x), x)**2,
        "kamke_6.117": a*y(x)*Derivative(y(x), x) + b*y(x)**2 + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.118": -2*a*y(x)**2 + a*y(x)*Derivative(y(x), x) + b*y(x)**3 + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.119": 2*a**2*y(x)**2 + a*y(x) - 2*b**2*y(x)**3 - (a*y(x) - 1)*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.120": -(-a**2 + b**2*y(x)**2)*(y(x) + 1)*y(x) + (a*y(x) - 1)*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.121": (-n**2*cot(x)**2 + cos(x)**2)*y(x)**2*log(y(x)) + (tan(x) + cot(x))*y(x)*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.122": -f(x)*y(x)*Derivative(y(x), x) - g(x)*y(x)**2 + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.123": (f(x)*y(x)**2 + g(x))*Derivative(y(x), x) - (-y(x)**2*Derivative(f(x), x) + Derivative(g(x), x))*y(x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.124": -y(x)**2 + 3*y(x)*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.125": -a*Derivative(y(x), x)**2 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.126": a*(Derivative(y(x), x)**2 + 1) + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.127": a*Derivative(y(x), x)**2 + b*y(x)**3 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.128": a*Derivative(y(x), x)**2 + b*y(x)*Derivative(y(x), x) + c*y(x)**2 + d*y(x)**(1 - a) + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.129": a*Derivative(y(x), x)**2 + f(x)*y(x)*Derivative(y(x), x) + g(x)*y(x)**2 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.130": a*Derivative(y(x), x)**2 + b*y(x)**2*Derivative(y(x), x) + c*y(x)**4 + y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.131": -a*y(x)**3*Derivative(f(x), x)/(a + 2) + a*f(x)**2*y(x)**4/(a + 2)**2 - f(x)*y(x)**2*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 2)) - (a - 1)*Derivative(y(x), x)**2/a,
        "kamke_6.132": -2*a*(Derivative(y(x), x)**2 + 1)**(3/2)*y(x) + y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 - 1,
        "kamke_6.133": (x + y(x))*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 - Derivative(y(x), x),
        "kamke_6.134": (x - y(x))*Derivative(y(x), (x, 2)) + 2*(Derivative(y(x), x) + 1)*Derivative(y(x), x),
        "kamke_6.135": (x - y(x))*Derivative(y(x), (x, 2)) - (Derivative(y(x), x) + 1)*(Derivative(y(x), x)**2 + 1),
        "kamke_6.136": (x - y(x))*Derivative(y(x), (x, 2)) - h(Derivative(y(x), x)),
        "kamke_6.137": 2*y(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 + 1,
        "kamke_6.138": a + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.139": a + f(x)*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.140": -8*y(x)**3 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.141": -8*y(x)**3 - 4*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.142": -(4*x + 8*y(x))*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.143": (a*y(x) + b)*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.144": a*y(x)**3 + 2*x*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 + 1,
        "kamke_6.145": (a*y(x) + b*x)*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.146": -3*y(x)**4 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.147": b - 8*x*y(x)**3 - (4*a + 4*x**2)*y(x)**2 - 3*y(x)**4 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.148": (2*f(x)**2 + 2*Derivative(f(x), x))*y(x)**2 + 3*f(x)*y(x)*Derivative(y(x), x) - 8*y(x)**3 + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2,
        "kamke_6.149": f(x)*y(x)**2 + y(x)**4 + 4*y(x)**2*Derivative(y(x), x) + 2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 + 1,
        "kamke_6.150": 2*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.151": -4*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.152": f(x)*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.153": (a*y(x)**3 + 1)*y(x)**2 + 2*y(x)*Derivative(y(x), (x, 2)) - 6*Derivative(y(x), x)**2,
        "kamke_6.154": -(Derivative(y(x), x)**2 + 1)*Derivative(y(x), x)**2 + 2*y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.155": (-2*a + 2*y(x))*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 + 1,
        "kamke_6.156": -a*x**2 - b*x - c + 3*y(x)*Derivative(y(x), (x, 2)) - 2*Derivative(y(x), x)**2,
        "kamke_6.157": 3*y(x)*Derivative(y(x), (x, 2)) - 5*Derivative(y(x), x)**2,
        "kamke_6.158": 4*y(x)*Derivative(y(x), (x, 2)) + 4*y(x) - 3*Derivative(y(x), x)**2,
        "kamke_6.159": -12*y(x)**3 + 4*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.160": a*y(x)**3 + b*y(x)**2 + c*y(x) + 4*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.161": (6*y(x)**2 - 2*y(x)*Derivative(f(x), x)/f(x))*Derivative(y(x), x) + f(x)*y(x) + g(x)*y(x)**2 + y(x)**4 - 2*y(x)**2*Derivative(y(x), x) + 4*y(x)*Derivative(y(x), (x, 2)) - 3*Derivative(y(x), x)**2,
        "kamke_6.162": a*y(x)**2 + 4*y(x)*Derivative(y(x), (x, 2)) - 5*Derivative(y(x), x)**2,
        "kamke_6.163": 8*y(x)**3 + 12*y(x)*Derivative(y(x), (x, 2)) - 15*Derivative(y(x), x)**2,
        "kamke_6.164": n*y(x)*Derivative(y(x), (x, 2)) - (n - 1)*Derivative(y(x), x)**2,
        "kamke_6.165": a*y(x)*Derivative(y(x), (x, 2)) + b*Derivative(y(x), x)**2 + c0 + c1*y(x) + c2*y(x)**2 + c3*y(x)**3 + c4*y(x)**4,
        "kamke_6.166": a*y(x)*Derivative(y(x), (x, 2)) + b*Derivative(y(x), x)**2 - y(x)*Derivative(y(x), x)/sqrt(c**2 + x**2),
        "kamke_6.167": a*y(x)**3*Derivative(f(x), x) + a*y(x)*Derivative(y(x), (x, 2)) - (a - 1)*Derivative(y(x), x)**2 + (a + 2)*f(x)*y(x)**2*Derivative(y(x), x) + f(x)**2*y(x)**4,
        "kamke_6.168": c*Derivative(y(x), x)**2 + (a*y(x) + b)*Derivative(y(x), (x, 2)),
        "kamke_6.169": x*y(x)*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x)**2 - y(x)*Derivative(y(x), x),
        "kamke_6.170": a*y(x)*Derivative(y(x), x) + x*y(x)*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x)**2 + f(x),
        "kamke_6.171": x*(a*y(x)**4 + d) + x*y(x)*Derivative(y(x), (x, 2)) - x*Derivative(y(x), x)**2 + (b*y(x)**2 + c)*y(x) + y(x)*Derivative(y(x), x),
        "kamke_6.172": a*y(x)*Derivative(y(x), x) + b*x*y(x)**3 + x*y(x)*Derivative(y(x), (x, 2)) - x*Derivative(y(x), x)**2,
        "kamke_6.173": a*y(x)*Derivative(y(x), x) + x*y(x)*Derivative(y(x), (x, 2)) + 2*x*Derivative(y(x), x)**2,
        "kamke_6.174": x*y(x)*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x)**2 + (y(x) + 1)*Derivative(y(x), x),
        "kamke_6.175": a*y(x)*Derivative(y(x), x) + x*y(x)*Derivative(y(x), (x, 2)) - 2*x*Derivative(y(x), x)**2,
        "kamke_6.176": x*y(x)*Derivative(y(x), (x, 2)) - 4*x*Derivative(y(x), x)**2 + 4*y(x)*Derivative(y(x), x),
        "kamke_6.177": x*y(x)*Derivative(y(x), (x, 2)) + (a*x/sqrt(b**2 - x**2) - x)*Derivative(y(x), x)**2 - y(x)*Derivative(y(x), x),
        "kamke_6.178": x*(x + y(x))*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x)**2 + (x - y(x))*Derivative(y(x), x) - y(x),
        "kamke_6.179": 2*x*y(x)*Derivative(y(x), (x, 2)) - x*Derivative(y(x), x)**2 + y(x)*Derivative(y(x), x),
        "kamke_6.180": x**2*(y(x) - 1)*Derivative(y(x), (x, 2)) - 2*x**2*Derivative(y(x), x)**2 - 2*x*(y(x) - 1)*Derivative(y(x), x) - 2*(y(x) - 1)**2*y(x),
        "kamke_6.181": x**2*(x + y(x))*Derivative(y(x), (x, 2)) - (x*Derivative(y(x), x) - y(x))**2,
        "kamke_6.182": a*(x*Derivative(y(x), x) - y(x))**2 + x**2*(x - y(x))*Derivative(y(x), (x, 2)),
        "kamke_6.183": -x**2*(Derivative(y(x), x)**2 + 1) + 2*x**2*y(x)*Derivative(y(x), (x, 2)) + y(x)**2,
        "kamke_6.184": a*x**2*y(x)*Derivative(y(x), (x, 2)) + b*x**2*Derivative(y(x), x)**2 + c*x*y(x)*Derivative(y(x), x) + d*y(x)**2,
        "kamke_6.185": -a*(x + 2)*y(x)**2 + x*(x + 1)**2*y(x)*Derivative(y(x), (x, 2)) - x*(x + 1)**2*Derivative(y(x), x)**2 + 2*(x + 1)**2*y(x)*Derivative(y(x), x),
        "kamke_6.186": -12*x**2*y(x)*Derivative(y(x), x) + 3*x*y(x)**2 - (4 - 4*x**3)*Derivative(y(x), x)**2 + (8 - 8*x**3)*y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.187": f0(x)*y(x)*Derivative(y(x), (x, 2)) + f_1(x)*Derivative(y(x), x)**2 + f_2(x)*y(x)*Derivative(y(x), x) + f_3(x)*y(x)**2,
        "kamke_6.188": -a + y(x)**2*Derivative(y(x), (x, 2)),
        "kamke_6.189": a*x + y(x)**2*Derivative(y(x), (x, 2)) + y(x)*Derivative(y(x), x)**2,
        "kamke_6.190": -a*x - b + y(x)**2*Derivative(y(x), (x, 2)) + y(x)*Derivative(y(x), x)**2,
        "kamke_6.191": (1 - 2*y(x))*Derivative(y(x), x)**2 + (y(x)**2 + 1)*Derivative(y(x), (x, 2)),
        "kamke_6.192": (y(x)**2 + 1)*Derivative(y(x), (x, 2)) - 3*y(x)*Derivative(y(x), x)**2,
        "kamke_6.193": (x + y(x)**2)*Derivative(y(x), (x, 2)) - (2*x - 2*y(x)**2)*Derivative(y(x), x)**3 + (4*y(x)*Derivative(y(x), x) + 1)*Derivative(y(x), x),
        "kamke_6.194": (x**2 + y(x)**2)*Derivative(y(x), (x, 2)) - (x*Derivative(y(x), x) - y(x))*(Derivative(y(x), x)**2 + 1),
        "kamke_6.195": (x**2 + y(x)**2)*Derivative(y(x), (x, 2)) - (x*Derivative(y(x), x) - y(x))*(2*Derivative(y(x), x)**2 + 2),
        "kamke_6.196": -(1 - 2*y(x))*Derivative(y(x), x)**2 + (1 - y(x))*f(x)*y(x)*Derivative(y(x), x) + 2*(1 - y(x))*y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.197": -(1 - 3*y(x))*Derivative(y(x), x)**2 + 2*(1 - y(x))*y(x)*Derivative(y(x), (x, 2)) + h(y(x)),
        "kamke_6.198": 4*(f(x)*y(x) + g(x))*y(x)*Derivative(y(x), x) + 4*(y(x) - 1)*(-f(x)**2 + g(x)**2 - Derivative(f(x), x) - Derivative(g(x), x))*y(x)**2 + 2*(y(x) - 1)*y(x)*Derivative(y(x), (x, 2)) - (3*y(x) - 1)*Derivative(y(x), x)**2,
        "kamke_6.199": (1 - 3*y(x))*Derivative(y(x), x)**2 + (1 - y(x))**3*(f0(x)**2*y(x)**2 - f_1(x)**2) + 4*(1 - y(x))*(f(x)**2 - g(x)**2 - Derivative(f(x), x) - Derivative(g(x), x))*y(x)**2 - 2*(1 - y(x))*y(x)*Derivative(y(x), (x, 2)) - 4*(f(x)*y(x) + g(x))*y(x)*Derivative(y(x), x),
        "kamke_6.200": 3*(1 - y(x))*y(x)*Derivative(y(x), (x, 2)) - (2 - 4*y(x))*Derivative(y(x), x)**2 - h(y(x)),
        "kamke_6.201": (1 - y(x))*Derivative(y(x), (x, 2)) - (3 - 6*y(x))*Derivative(y(x), x)**2 - h(y(x)),
        "kamke_6.202": a*(y(x) - 1)*y(x)*Derivative(y(x), (x, 2)) + (b*y(x) + c)*Derivative(y(x), x)**2 + h(y(x)),
        "kamke_6.203": a*(y(x) - 1)*y(x)*Derivative(y(x), (x, 2)) + fs*(y(x) - 1)*y(x)*Derivative(y(x), x) - (a - 1)*(2*y(x) - 1)*Derivative(y(x), x)**2,
        "kamke_6.204": a*b*(y(x) - 1)*y(x)*Derivative(y(x), (x, 2)) + fs*(y(x) - 1)*y(x)*Derivative(y(x), x) - (b*(1 - a) + (2*a*b - a - b)*y(x))*Derivative(y(x), x)**2,
        "kamke_6.205": -a + x*y(x)**2*Derivative(y(x), (x, 2)),
        "kamke_6.206": -x*(a**2 - y(x)**2)*Derivative(y(x), x) + (a**2 - x**2)*(a**2 - y(x)**2)*Derivative(y(x), (x, 2)) + (a**2 - x**2)*y(x)*Derivative(y(x), x)**2,
        "kamke_6.207": c*x*(y(x) - 1)*y(x)**2 + d*x**2*(y(x) + 1)*y(x)**2 + 2*x**2*(y(x) - 1)*y(x)*Derivative(y(x), (x, 2)) - x**2*(3*y(x) - 1)*Derivative(y(x), x)**2 + 2*x*(y(x) - 1)*y(x)*Derivative(y(x), x) + (a*y(x)**2 + b)*(y(x) - 1)**3,
        "kamke_6.208": x**3*y(x)**2*Derivative(y(x), (x, 2)) + (x + y(x))*(x*Derivative(y(x), x) - y(x))**3,
        "kamke_6.209": -a + y(x)**3*Derivative(y(x), (x, 2)),
        "kamke_6.210": (1 - 3*y(x)**2)*Derivative(y(x), x)**2 + (y(x)**2 + 1)*y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.211": -a**2*x*y(x)**2 + y(x)**4 + 2*y(x)**3*Derivative(y(x), (x, 2)) - 1,
        "kamke_6.212": -a*x**2 - b*x - c + 2*y(x)**3*Derivative(y(x), (x, 2)) + y(x)**2*Derivative(y(x), x)**2,
        "kamke_6.213": -a0*(a - y(x))**2*(b - y(x))**2*(c - y(x))**2 - a1*(b - y(x))**2*(c - y(x))**2 - a2*(a - y(x))**2*(c - y(x))**2 - a3*(a - y(x))**2*(b - y(x))**2 + (a - y(x))*(b - y(x))*(2*c - 2*y(x))*Derivative(y(x), (x, 2)) + ((a - y(x))*(b - y(x)) + (a - y(x))*(c - y(x)) + (b - y(x))*(c - y(x)))*Derivative(y(x), x)**2,
        "kamke_6.214": -(-a/2 + 6*y(x)**2)*Derivative(y(x), x)**2 + (-a*y(x) - b + 4*y(x)**3)*Derivative(y(x), (x, 2)),
        "kamke_6.215": -(-a/2 + 6*y(x)**2)*Derivative(y(x), x)**2 + (fs*Derivative(y(x), x) + Derivative(y(x), (x, 2)))*(-a*y(x) - b + 4*y(x)**3),
        "kamke_6.216": -fs*((-x + y(x))*(y(x) - 1)*y(x))**(3/2) - 2*x*(1 - x)*(1 - y(x))*(x - y(x))*y(x)*Derivative(y(x), (x, 2)) + x*(1 - x)*(-2*x*y(x) + x + 3*y(x)**2 - 2*y(x))*Derivative(y(x), x)**2 - (1 - y(x))**2*y(x)**2 + 2*(1 - y(x))*(x**2 - 2*x*y(x) + y(x))*y(x)*Derivative(y(x), x),
        "kamke_6.217": a*(1 - y(x))**2*(x - y(x))**2*y(x)**2 + b*x*(1 - y(x))**2*(x - y(x))**2 - c*(1 - x)*(x - y(x))**2*y(x)**2 - d*x*(1 - x)*(1 - y(x))**2*y(x)**2 + 2*x**2*(1 - x)**2*(1 - y(x))*(x - y(x))*y(x)*Derivative(y(x), (x, 2)) - x**2*(1 - x)**2*(-2*x*y(x) + x + 3*y(x)**2 - 2*y(x))*Derivative(y(x), x)**2 - 2*x*(1 - x)*(1 - y(x))*(x**2 - 2*x*y(x) + y(x))*y(x)*Derivative(y(x), x),
        "kamke_6.218": b*sqrt((1 - y(x)**2)*(-a**2*y(x)**2 + 1))*Derivative(y(x), x)**2 + (a**2*y(x)**2 - 1)*(y(x)**2 - 1)*Derivative(y(x), (x, 2)) + (-2*a**2*y(x)**2 + a**2 + 1)*y(x)*Derivative(y(x), x)**2,
        "kamke_6.219": d*y(x) + (a*x**2 + 2*b*x + c + y(x)**2)**2*Derivative(y(x), (x, 2)),
        "kamke_6.220": -a + sqrt(y(x))*Derivative(y(x), (x, 2)),
        "kamke_6.221": -a*(Derivative(y(x), x)**2 + 1)**(3/2) + sqrt(x**2 + y(x)**2)*Derivative(y(x), (x, 2)),
        "kamke_6.222": (1 - log(y(x)))*y(x)*Derivative(y(x), (x, 2)) + (log(y(x)) + 1)*Derivative(y(x), x)**2,
        "kamke_6.223": A*(a*sin(y(x))**2 + c)*y(x) + a*sin(y(x))*cos(y(x))*Derivative(y(x), x)**2 + (a*sin(y(x))**2 + b)*Derivative(y(x), (x, 2)),
        "kamke_6.224": a*Derivative(h(y(x)), y(x))*Derivative(y(x), x)**2 + h(y(x))*Derivative(y(x), (x, 2)) + j(y(x)),
        "kamke_6.225": -h(y(x))**2*j(x, Derivative(y(x), x)/h(y(x))) + h(y(x))*Derivative(y(x), (x, 2)) - Derivative(h(y(x)), y(x))*Derivative(y(x), x)**2,
        "kamke_6.226": -x**2*y(x)*Derivative(y(x), x) - x*y(x)**2 + Derivative(y(x), x)*Derivative(y(x), (x, 2)),
        "kamke_6.227": (x*Derivative(y(x), x) - y(x))*Derivative(y(x), (x, 2)) + 4*Derivative(y(x), x)**2,
        "kamke_6.228": (x*Derivative(y(x), x) - y(x))*Derivative(y(x), (x, 2)) - (Derivative(y(x), x)**2 + 1)**2,
        "kamke_6.229": a*x**3*Derivative(y(x), x)*Derivative(y(x), (x, 2)) + b*y(x)**2,
        "kamke_6.230": f3*Derivative(y(x), x)**2 + (f1*Derivative(y(x), x) + f2*y(x))*Derivative(y(x), (x, 2)) + f4(x)*y(x)*Derivative(y(x), x) + f5(x)*y(x)**2,
        "kamke_6.231": 3*x*Derivative(y(x), x) + (x**2 + 2*y(x)**2*Derivative(y(x), x))*Derivative(y(x), (x, 2)) + 2*y(x)*Derivative(y(x), x)**3 + y(x),
        "kamke_6.232": (y(x)**2 + Derivative(y(x), x)**2)*Derivative(y(x), (x, 2)) + y(x)**3,
        "kamke_6.233": -b + (a*(x*Derivative(y(x), x) - y(x)) + Derivative(y(x), x)**2)*Derivative(y(x), (x, 2)),
        "kamke_6.234": (a*sqrt(Derivative(y(x), x)**2 + 1) - x*Derivative(y(x), x))*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 - 1,
        "kamke_6.235": fs + h(Derivative(y(x), x))*Derivative(y(x), (x, 2)) + j(y(x))*Derivative(y(x), x),
        "kamke_6.236": -a*y(x) - b + Derivative(y(x), (x, 2))**2,
        "kamke_6.237": a**2*Derivative(y(x), (x, 2))**2 - 2*a*x*Derivative(y(x), (x, 2)) + Derivative(y(x), x),
        "kamke_6.238": -x*(x + 4*Derivative(y(x), x))*Derivative(y(x), (x, 2)) + (2*x + 2*Derivative(y(x), x))*Derivative(y(x), x) + (2*x**2 + 2)*Derivative(y(x), (x, 2))**2 - 2*y(x),
        "kamke_6.239": 3*x**2*Derivative(y(x), (x, 2))**2 - (6*x*Derivative(y(x), x) + 2*y(x))*Derivative(y(x), (x, 2)) + 4*Derivative(y(x), x)**2,
        "kamke_6.240": x**2*(2 - 9*x)*Derivative(y(x), (x, 2))**2 - 6*x*(1 - 6*x)*Derivative(y(x), x)*Derivative(y(x), (x, 2)) - 36*x*Derivative(y(x), x)**2 + 6*y(x)*Derivative(y(x), (x, 2)),
        "kamke_6.241": ((F01(x) + F10(x))*y(x) + (F12(x) + F21(x))*Derivative(y(x), (x, 2)))*Derivative(y(x), x) + (F02(x) + F20(x))*y(x)*Derivative(y(x), (x, 2)) + F00(x)*y(x)**2 + F11(x)*Derivative(y(x), x)**2 + F22(x)*Derivative(y(x), (x, 2))**2,
        "kamke_6.242": -a*exp(2*x) + y(x)*Derivative(y(x), (x, 2))**2,
        "kamke_6.243": -2*a**2*y(x)*Derivative(y(x), x)**2*Derivative(y(x), (x, 2)) + (a**2*y(x)**2 - b**2)*Derivative(y(x), (x, 2))**2 + (a**2*Derivative(y(x), x)**2 - 1)*Derivative(y(x), x)**2,
        "kamke_6.244": -4*x*(x*Derivative(y(x), x) - y(x))**3*y(x) + (x**2*y(x)*Derivative(y(x), (x, 2)) - x**2*Derivative(y(x), x)**2 + y(x)**2)**2,
        "kamke_6.245": 32*(x*Derivative(y(x), (x, 2)) - Derivative(y(x), x))**3*Derivative(y(x), (x, 2)) + (2*y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2)**3,
        "kamke_6.246": c*y(x)*Derivative(y(x), (x, 2)) + d*Derivative(y(x), x)**2 + sqrt(a*Derivative(y(x), (x, 2))**2 + b*Derivative(y(x), x)**2)
    }

    chapter_7 = {
        "kamke_7.1": -a**2*(Derivative(y(x), x)**5 + 2*Derivative(y(x), x)**3 + Derivative(y(x), x)) + Derivative(y(x), (x, 3)),
        "kamke_7.2": y(x)*Derivative(y(x), (x, 2)) - Derivative(y(x), x)**2 + Derivative(y(x), (x, 3)) + 1,
        "kamke_7.3": -y(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), x)**2 + Derivative(y(x), (x, 3)),
        "kamke_7.4": a*y(x)*Derivative(y(x), (x, 2)) + Derivative(y(x), (x, 3)),
        "kamke_7.5": x**2*Derivative(y(x), (x, 3)) + x*Derivative(y(x), (x, 2)) + (2*x*y(x) - 1)*Derivative(y(x), x) - f(x) + y(x)**2,
        "kamke_7.6": x**2*Derivative(y(x), (x, 3)) + x*(y(x) - 1)*Derivative(y(x), (x, 2)) + x*Derivative(y(x), x)**2 + (1 - y(x))*Derivative(y(x), x),
        "kamke_7.7": y(x)**3*Derivative(y(x), x) + y(x)*Derivative(y(x), (x, 3)) - Derivative(y(x), x)*Derivative(y(x), (x, 2)),
        "kamke_7.8": 4*y(x)**2*Derivative(y(x), (x, 3)) - 18*y(x)*Derivative(y(x), x)*Derivative(y(x), (x, 2)) + 15*Derivative(y(x), x)**3,
        "kamke_7.9": 9*y(x)**2*Derivative(y(x), (x, 3)) - 45*y(x)*Derivative(y(x), x)*Derivative(y(x), (x, 2)) + 40*Derivative(y(x), x)**3,
        "kamke_7.10": -3*Derivative(y(x), x)**2 + 2*Derivative(y(x), x)*Derivative(y(x), (x, 3)),
        "kamke_7.11": (Derivative(y(x), x)**2 + 1)*Derivative(y(x), (x, 3)) - 3*Derivative(y(x), x)*Derivative(y(x), (x, 2))**2,
        "kamke_7.12": -(a + 3*Derivative(y(x), x))*Derivative(y(x), (x, 2))**2 + (Derivative(y(x), x)**2 + 1)*Derivative(y(x), (x, 3)),
        "kamke_7.13": -a*sqrt(b**2*Derivative(y(x), (x, 2))**2 + 1) + Derivative(y(x), (x, 2))*Derivative(y(x), (x, 3)),
        "kamke_7.14": Derivative(y(x), x)**3*Derivative(y(x), (x, 3)) + Derivative(y(x), x)*Derivative(y(x), (x, 4)) - Derivative(y(x), (x, 2))*Derivative(y(x), (x, 3)),
        "kamke_7.15": -fs*Derivative(y(x), (x, 2))*Derivative(y(x), (x, 3)) + (f(x)*Derivative(y(x), (x, 2)) + Derivative(f(x), x)*Derivative(y(x), x))*Derivative(y(x), x)**3 + (q(x)*Derivative(y(x), (x, 2)) - Derivative(q(x), x)*Derivative(y(x), x))*cos(y(x)) + (f(x)*Derivative(y(x), (x, 4)) + 3*Derivative(f(x), x)*Derivative(y(x), (x, 3)) + 3*Derivative(f(x), (x, 2))*Derivative(y(x), (x, 2)) + Derivative(f(x), (x, 3))*Derivative(y(x), x))*Derivative(y(x), x) + 2*q(x)*sin(y(x))*Derivative(y(x), x)**2,
        "kamke_7.16": 3*Derivative(y(x), (x, 2))*Derivative(y(x), (x, 4)) - 5*Derivative(y(x), (x, 3))**2,
        "kamke_7.17": 9*Derivative(y(x), (x, 2))**2*Derivative(y(x), (x, 5)) - 45*Derivative(y(x), (x, 2))*Derivative(y(x), (x, 3))*Derivative(y(x), (x, 4)) + 40*Derivative(y(x), (x, 3)),
        "kamke_7.18": -f(Derivative(y(x), (x, n - 1))) + Derivative(y(x), (x, n)),
        "kamke_7.19": -f(Derivative(y(x), (x, n - 2))) + Derivative(y(x), (x, n))
    }

    chapter_8 = {
        "kamke_8.1": ([-a*xf(t) + Derivative(xf(t), t), -b + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.2": ([-a*y(t) + Derivative(xf(t), t), a*xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.3": ([-a*y(t) + Derivative(xf(t), t), -b*xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.4": ([-a*xf(t) + y(t) + Derivative(xf(t), t), -a*y(t) - xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.5": ([-a*xf(t) - b*y(t) + Derivative(xf(t), t), -b*y(t) - c*xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.6": ([a*Derivative(xf(t), t) - alpha*xf(t) + b*Derivative(y(t), t) - bbeta*y(t), -a*Derivative(y(t), t) + alpha*y(t) + b*Derivative(xf(t), t) - bbeta*xf(t)], [xf(t), y(t)]),
        "kamke_8.7": ([y(t) + Derivative(xf(t), t), -2*xf(t) - 2*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.8": ([3*xf(t) + 4*y(t) + Derivative(xf(t), t), 2*xf(t) + 5*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.9": ([5*xf(t) + 2*y(t) + Derivative(xf(t), t), -xf(t) + 7*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.10": ([-a1*xf(t) - b1*y(t) - c1 + Derivative(xf(t), t), -a2*xf(t) - b2*y(t) - c2 + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.11": ([-3*t + 2*y(t) + Derivative(xf(t), t), -2*xf(t) + Derivative(y(t), t) - 4], [xf(t), y(t)]),
        "kamke_8.12": ([-t**2 + 6*t + y(t) + Derivative(xf(t), t) + 1, 3*t**2 - 3*t - xf(t) + Derivative(y(t), t) - 1], [xf(t), y(t)]),
        "kamke_8.13": ([3*xf(t) - y(t) - exp(2*t) + Derivative(xf(t), t), xf(t) + 5*y(t) - exp(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.14": ([-t + 2*xf(t) + y(t) - exp(2*t) + Derivative(xf(t), t) + Derivative(y(t), t), -xf(t) + 3*y(t) - exp(t) + Derivative(xf(t), t) + Derivative(y(t), t) - 1], [xf(t), y(t)]),
        "kamke_8.15": ([-y(t) - exp(t) + Derivative(xf(t), t) + Derivative(y(t), t), 2*y(t) - cos(t) + 2*Derivative(xf(t), t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.16": ([2*xf(t) + 31*y(t) - exp(t) + 4*Derivative(xf(t), t) + 9*Derivative(y(t), t), xf(t) + 24*y(t) + 3*Derivative(xf(t), t) + 7*Derivative(y(t), t) - 3], [xf(t), y(t)]),
        "kamke_8.17": ([11*xf(t) + 31*y(t) - exp(t) + 4*Derivative(xf(t), t) + 9*Derivative(y(t), t), 8*xf(t) + 24*y(t) - exp(2*t) + 3*Derivative(xf(t), t) + 7*Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.18": ([-t + 44*xf(t) + 49*y(t) + 4*Derivative(xf(t), t) + 9*Derivative(y(t), t), 34*xf(t) + 38*y(t) - exp(t) + 3*Derivative(xf(t), t) + 7*Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.19": ([-f(t)*xf(t) - g(t)*y(t) + Derivative(xf(t), t), -f(t)*y(t) + g(t)*xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.20": ([(a*xf(t) + b*y(t))*f(t) - g(t) + Derivative(xf(t), t), (c*xf(t) + d*y(t))*f(t) - h(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.21": ([-xf(t)*cos(t) + Derivative(xf(t), t), -xf(t)*exp(-sin(t)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.22": ([t*Derivative(xf(t), t) + y(t), t*Derivative(y(t), t) + xf(t)], [xf(t), y(t)]),
        "kamke_8.23": ([t*Derivative(xf(t), t) - t + 2*xf(t), -t*y(t) + t*Derivative(y(t), t) + t - (t + 2)*xf(t)], [xf(t), y(t)]),
        "kamke_8.24": ([t*Derivative(xf(t), t) - t + 2*xf(t) - 2*y(t), -t**2 + t*Derivative(y(t), t) + xf(t) + 5*y(t)], [xf(t), y(t)]),
        "kamke_8.25": ([t**2*(1 - sint(t))*Derivative(xf(t), t) - t**2*y(t) - t*(1 - 2*sin(t))*xf(t), t**2*(1 - sint(t))*Derivative(y(t), t) - t*(-t*cos(t) + 1)*y(t) - (t*cos(t) - sin(t))*xf(t)], [xf(t), y(t)]),
        "kamke_8.26": ([-f(t) + y(t) + Derivative(xf(t), t) + Derivative(y(t), t), -g(t) + xf(t) + y(t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.27": ([-3*xf(t) + 2*Derivative(xf(t), t) + Derivative(y(t), t), -2*y(t) - exp(2*t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.28": ([-2*t + xf(t) + Derivative(xf(t), t) - Derivative(y(t), t), -9*xf(t) + 3*y(t) - sin(2*t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.29": ([-xf(t) + 2*y(t) + Derivative(xf(t), t), -2*t + cos(2*t) + Derivative(xf(t), (t, 2)) - 2*Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_8.30": ([t*Derivative(xf(t), t) - t*Derivative(y(t), t) - 2*y(t), t*xf(t) + t*Derivative(xf(t), (t, 2)) + 2*Derivative(xf(t), t)], [xf(t), y(t)]),
        "kamke_8.31": ([a*y(t) + Derivative(xf(t), (t, 2)), -a**2*y(t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.32": ([-a*xf(t) - b*y(t) + Derivative(xf(t), (t, 2)), -c*xf(t) - d*y(t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.33": ([-a1*xf(t) - b1*y(t) - c1 + Derivative(xf(t), (t, 2)), -a2*xf(t) - b2*y(t) - c2 + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.34": ([xf(t) + y(t) + Derivative(xf(t), (t, 2)) + 5, -4*xf(t) - 3*y(t) + Derivative(y(t), (t, 2)) + 3], [xf(t), y(t)]),
        "kamke_8.35": ([-c**2*(3*cos(a*t + b)**2 - 1)*xf(t) - 3*c**2*y(t)*sin(2*a*b*t)/2 + Derivative(xf(t), (t, 2)), -c**2*(3*sin(a*t + b)**2 - 1)*y(t) - 3*c**2*xf(t)*sin(2*a*b*t)/2 + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.36": ([6*xf(t) + 7*y(t) + Derivative(xf(t), (t, 2)), -2*t + 3*xf(t) + 2*y(t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.37": ([-a*Derivative(y(t), t) + b*xf(t) + Derivative(xf(t), (t, 2)), a*Derivative(xf(t), t) + b*y(t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.38": ([-A*Derivative(y(t), t) - B*exp(I*omega*t) + a1*Derivative(xf(t), (t, 2)) + b1*Derivative(xf(t), t) + c1*xf(t), A*Derivative(xf(t), t) + a2*Derivative(y(t), (t, 2)) + b2*Derivative(y(t), t) + c2*y(t)], [xf(t), y(t)]),
        "kamke_8.39": ([a*(Derivative(xf(t), t) - Derivative(y(t), t)) + b1*xf(t) - c1*exp(I*omega*t) + Derivative(xf(t), (t, 2)), a*(-Derivative(xf(t), t) + Derivative(y(t), t)) + b2*y(t) - c2*exp(I*omega*t) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.40": ([a11*Derivative(xf(t), (t, 2)) + a12*Derivative(y(t), (t, 2)) + b11*Derivative(xf(t), t) + b12*Derivative(y(t), t) + c11*xf(t) + c12*y(t), a21*Derivative(xf(t), (t, 2)) + a22*Derivative(y(t), (t, 2)) + b21*Derivative(xf(t), t) + b22*Derivative(y(t), t) + c21*xf(t) + c22*y(t)], [xf(t), y(t)]),
        "kamke_8.41": ([y(t) - 2*Derivative(xf(t), t) + Derivative(xf(t), (t, 2)) - Derivative(y(t), t), -t - xf(t) + 2*Derivative(xf(t), t) - Derivative(y(t), (t, 2)) + Derivative(y(t), (t, 3))], [xf(t), y(t)]),
        "kamke_8.42": ([-sinh(2*t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), t) + Derivative(y(t), (t, 2)), -2*t + 2*Derivative(xf(t), (t, 2)) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.43": ([-Derivative(xf(t), t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), t), -xf(t) + Derivative(xf(t), (t, 2)) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_8.44": ([-2*xf(t) + Derivative(xf(t), t), -3*xf(t) + 2*y(t) + Derivative(y(t), t), -2*y(t) - 3*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.45": ([-4*xf(t) + Derivative(xf(t), t), -xf(t) + 2*y(t) + Derivative(y(t), t), -xf(t) + 4*y(t) - zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.46": ([-y(t) + zf(t) + Derivative(xf(t), t), -xf(t) - y(t) + Derivative(y(t), t), -xf(t) - zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.47": ([-y(t) + zf(t) + Derivative(xf(t), t), -t - xf(t) - y(t) + Derivative(y(t), t), -t - xf(t) - zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.48": ([a*Derivative(xf(t), t) - b*c*(y(t) - zf(t)), -a*c*(-xf(t) + zf(t)) + b*Derivative(y(t), t), -a*b*(xf(t) - y(t)) + c*Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.49": ([b*zf(t) - c*y(t) + Derivative(xf(t), t), -a*zf(t) + c*xf(t) + Derivative(y(t), t), -a*y(t) - b*xf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.50": ([g(t)*zf(t) - h(t)*y(t) + Derivative(xf(t), t), -f(t)*zf(t) + h(t)*xf(t) + Derivative(y(t), t), f(t)*y(t) - g(t)*xf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.51": ([-xf(t) - y(t) + zf(t) + Derivative(xf(t), t), xf(t) - y(t) - zf(t) + Derivative(y(t), t), -xf(t) + y(t) - zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.52": ([3*xf(t) - 48*y(t) + 28*zf(t) + Derivative(xf(t), t), 4*xf(t) - 40*y(t) + 22*zf(t) + Derivative(y(t), t), 6*xf(t) - 57*y(t) + 31*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.53": ([-6*xf(t) + 72*y(t) - 44*zf(t) + Derivative(xf(t), t), -4*xf(t) + 4*y(t) - 26*zf(t) + Derivative(y(t), t), -6*xf(t) + 63*y(t) - 38*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.54": ([-a*xf(t) - bbeta*zf(t) - gs*y(t) + Derivative(xf(t), t), -alpha*zf(t) - b*y(t) - gs*xf(t) + Derivative(y(t), t), -alpha*y(t) - bbeta*xf(t) - c*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.55": ([t*Derivative(xf(t), t) + t - 2*xf(t), t**3*Derivative(y(t), t) - t**2*y(t) - t + xf(t), t**4*Derivative(zf(t), t) - t**3*zf(t) + t**2*y(t) - t + xf(t)], [xf(t), y(t), zf(t)]),
        "kamke_8.56": ([a*t*Derivative(xf(t), t) - b*c*(y(t) - zf(t)), -a*c*(-xf(t) + zf(t)) + b*t*Derivative(y(t), t), -a*b*(xf(t) - y(t)) + c*t*Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_8.57": ([-a*x2(t) - b*x4*sin(c*t) - b*x3(t)*cos(c*t) + Derivative(x1(t), t), a*x1(t) - b*x3(t)*sin(c*t) + b*xf4(t)*cos(c*t) + Derivative(x2(t), t), -a*xf4(t) + b*x1(t)*cos(c*t) + b*x2(t)*sin(c*t) + Derivative(x3(t), t), a*x3(t) + b*x1(t)*sin(c*t) - b*x2(t)*cos(c*t) + Derivative(xf4(t), t)], [x1(t), x2(t), x3(t), xf4(t)])
    }

    chapter_9 = {
        "kamke_9.1": ([(xf(t) + y(t))*xf(t) + Derivative(xf(t), t), (xf(t) + y(t))*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.2": ([-(a*y(t) + b)*xf(t) + Derivative(xf(t), t), -(c*xf(t) + d)*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.3": ([-(a*(ps*xf(t) + qs*y(t)) + alpha)*xf(t) + Derivative(xf(t), t), -(b*(ps*xf(t) + qs*y(t)) + bbeta)*y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.4": ([-hs*(a - xf(t))*(c - xf(t) - y(t)) + Derivative(xf(t), t), -k*(b - y(t))*(c - xf(t) - y(t)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.5": ([-y(t)**2 + cos(xf(t)) + Derivative(xf(t), t), y(t)*sin(xf(t)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.6": ([xf(t)*y(t)**2 - xf(t) - y(t) + Derivative(xf(t), t), -xf(t)**2*y(t) + xf(t) + y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.7": ([(xf(t)**2 + y(t)**2)*xf(t) - xf(t) - y(t) + Derivative(xf(t), t), (xf(t)**2 + y(t)**2)*y(t) + xf(t) - y(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.8": ([-(xf(t)**2 + y(t)**2 - 1)*xf(t) + y(t) + Derivative(xf(t), t), -(xf(t)**2 + y(t)**2 - 1)*y(t) - xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.9": ([(xf(t)**2 + y(t)**2)*y(t) + Derivative(xf(t), t), -Piecewise((xf(t)**2 + y(t)**2, 2*xf(t) <= xf(t)**2 + y(t)**2), ((xf(t)/2 - y(t)**2/(2*xf(t)))*(xf(t)**2 + y(t)**2), True)) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.10": ([-Piecewise(((xf(t)**2 + y(t)**2 - 1)*xf(t)*sin(1/(xf(t)**2 + y(t)**2)), Ne(xf(t)**2 + y(t)**2, 1))) + y(t) + Derivative(xf(t), t), -Piecewise(((xf(t)**2 + y(t)**2 - 1)*y(t)*sin(1/(xf(t)**2 + y(t)**2)), Ne(xf(t)**2 + y(t)**2, 1))) - xf(t) + Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.11": ([t*xf(t) + (t**2 + 1)*Derivative(xf(t), t) - y(t), t*y(t) + (t**2 + 1)*Derivative(y(t), t) + xf(t)], [xf(t), y(t)]),
        "kamke_9.12": ([2*t*xf(t) + (-t**2 + xf(t)**2 + y(t)**2)*Derivative(xf(t), t), 2*t*y(t) + (-t**2 + xf(t)**2 + y(t)**2)*Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.13": ([a*Derivative(y(t), t) + t*Derivative(xf(t), t) - xf(t) + Derivative(xf(t), t)**2, t*Derivative(y(t), t) - y(t) + Derivative(xf(t), t)*Derivative(y(t), t)], [xf(t), y(t)]),
        "kamke_9.14": ([-t*Derivative(xf(t), t) - f(Derivative(xf(t), t), Derivative(y(t), t)) + xf(t), -t*Derivative(y(t), t) - g(Derivative(xf(t), t), Derivative(y(t), t)) + y(t)], [xf(t), y(t)]),
        "kamke_9.15": ([-a*exp(2*xf(t)) + Derivative(xf(t), (t, 2)) + exp(-xf(t)) - exp(-2*xf(t))*cos(y(t))**2, sin(y(t))/cos(y(t))**3 + Derivative(y(t), (t, 2)) - exp(-2*xf(t))*sin(y(t))*cos(y(t))], [xf(t), y(t)]),
        "kamke_9.16": ([-k*xf(t)/(xf(t)**2 + y(t)**2)**(3/2) + Derivative(xf(t), (t, 2)), -k*y(t)/(xf(t)**2 + y(t)**2)**(3/2) + Derivative(y(t), (t, 2))], [xf(t), y(t)]),
        "kamke_9.17": ([Derivative(xf(t), (t, 2)) + Cf(y(t))*f(sqrt(Derivative(y(t), t)**2))*Derivative(xf(t), t)/sqrt(Derivative(y(t), t)**2), gs + Derivative(y(t), (t, 2)) + Cf(y(t))*f(sqrt(Derivative(y(t), t)**2))*Derivative(y(t), t)/sqrt(Derivative(y(t), t)**2)], [xf(t), y(t)]),
        "kamke_9.18": ([-y(t) + zf(t) + Derivative(xf(t), t), -xf(t)**2 - y(t) + Derivative(y(t), t), -xf(t)**2 - zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.19": ([a*Derivative(xf(t), t) + (-b + c)*y(t)*zf(t), b*Derivative(y(t), t) + (a - c)*xf(t)*zf(t), c*Derivative(zf(t), t) + (-a + b)*xf(t)*y(t)], [xf(t), y(t), zf(t)]),
        "kamke_9.20": ([-(y(t) - zf(t))*xf(t) + Derivative(xf(t), t), -(-xf(t) + zf(t))*y(t) + Derivative(y(t), t), -(xf(t) - y(t))*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.21": ([-xf(t)*y(t) + Derivative(xf(t), t) + Derivative(y(t), t), -y(t)*zf(t) + Derivative(y(t), t) + Derivative(zf(t), t), -xf(t)*zf(t) + Derivative(xf(t), t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.22": ([-xf(t)**2/2 + y(t)/24 + Derivative(xf(t), t), -2*xf(t)*y(t) + 3*zf(t) + Derivative(y(t), t), -3*xf(t)*zf(t) + y(t)**2/6 + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.23": ([-(y(t)**2 - zf(t)**2)*xf(t) + Derivative(xf(t), t), -(-xf(t)**2 + zf(t)**2)*y(t) + Derivative(y(t), t), -(xf(t)**2 - y(t)**2)*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.24": ([-(y(t)**2 - zf(t)**2)*xf(t) + Derivative(xf(t), t), (xf(t)**2 + zf(t)**2)*y(t) + Derivative(y(t), t), -(xf(t)**2 + y(t)**2)*zf(t) + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.25": ([xf(t)*y(t)**2 - xf(t) - y(t) + Derivative(xf(t), t), -xf(t)**2*y(t) + xf(t) + y(t) + Derivative(y(t), t), xf(t)**2 - y(t)**2 + Derivative(zf(t), t)], [xf(t), y(t), zf(t)]),
        "kamke_9.26": ([Derivative(xf(t), (t, 2)) - xf(t)*Derivative(F(r), r)/r, Derivative(y(t), (t, 2)) - y(t)*Derivative(F(r), r)/r, Derivative(zf(t), (t, 2)) - zf(t)*Derivative(F(r), r)/r], [xf(t), y(t), zf(t)]),
        "kamke_9.27": ([(xf(t) - y(t))*(xf(t) - zf(t))*Derivative(xf(t), t) - f(t), (-xf(t) + y(t))*(y(t) - zf(t))*Derivative(y(t), t) - f(t), (-xf(t) + zf(t))*(-y(t) + zf(t))*Derivative(zf(t), t) - f(t)], [xf(t), y(t), zf(t)]),
        "kamke_9.28": ([-xf4(t)*sin(x3(t)) - x5(t)*cos(x3(t)) + sin(x2(t))*Derivative(x1(t), t), -xf4(t)*cos(x3(t)) + x5(t)*sin(x3(t)) + Derivative(x2(t), t), -a + cos(x2(t))*Derivative(x1(t), t) + Derivative(x3(t), t), -a*(1 - lambda_)*x5(t) + m*sin(x2(t))*cos(x3(t)) + Derivative(xf4(t), t), a*(1 - lambda_)*xf4(t) - m*sin(x2(t))*sin(x3(t)) + Derivative(x5(t), t)], [x1(t), x2(t), x3(t), xf4(t), x5(t)])
    }

    failing_examples = {
        "kamke_1.14": a*x**m + y(x)**2 + Derivative(y(x), x),
        "kamke_1.24": a*y(x)**2 - b*x**nu + Derivative(y(x), x),
        "kamke_1.31": -a*x**n*(y(x)**2 + 1) + Derivative(y(x), x),
        "kamke_1.47": -a*(-x + x**n)*y(x)**3 - y(x)**2 + Derivative(y(x), x),
        "kamke_1.48": -c*y(x)**2 - (a*x**n + b*x)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.94": a*y(x) + b*x**n + x*Derivative(y(x), x),
        "kamke_1.113": a*sqrt(x**2 + y(x)**2) + x*Derivative(y(x), x) - y(x),
        "kamke_1.205": a*y(x) + b*x**n + x*(a**2/4 - 1/4) + y(x)*Derivative(y(x), x),
        "kamke_5.9": x*diff(y(x),(x,n))-((a*AA[1]-AA[0])*x+AA[1])-Sum(((a*AA[v+1]-AA[v])*x+AA[v+1])*diff(y(x),(x,v)),(v,1,n-1)),
    }

    all_chapters = [chapter_1, chapter_2, chapter_3, chapter_4, chapter_5, chapter_6, chapter_7, chapter_8, chapter_9]

    css = """
    .container {
        max-width: 1000px;
        margin-left: auto;
        margin-right: auto;
        padding-left: 10px;
        padding-right: 10px;
        text-align: center;
    }

    h1, h2, h3 {
        margin: 20px 0;
        text-align: center;
    }

    li {
        border-radius: 3px;
        padding: 25px 30px;
        display: flex;
        justify-content: space-between;
        margin-bottom: 5px;
    }

    .table-header {
        background-color: rgba(0,0,0,0.75);
        color: white;
        font-size: 14px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.03em;
    }

    .table-row {
        background-color: #ffffff;
        box-shadow: 0px 0px 9px 0px rgba(0,0,0,0.1);
    }
    """


    def create_example_page(self, example, eq, status, sol, time, classify_output, checkodesol_output, all_hints):
        exno = int(example.split('.')[1])
        chno = int(example[6])
        prev = ""
        nxt = ""
        backslash = "\\"
        hints_section = f"""<h4>Matching Hints</h4>\n<ul>\n{classify_output}\n</ul>\n<br>\n"""
        if all_hints:
            hints_section = f"<h4>All Solutions</h4>\n{classify_output}"
        if exno != 1:
            prev = f"<a href='kamke_{chno}.{exno-1}.html'>&laquo; Previous</a>"
        if exno != len(self.all_chapters[chno-1]):
            nxt = f"<a href='kamke_{chno}.{exno+1}.html'>Next &raquo;</a>"

        example_page = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width">
        <title>{example.capitalize().replace("_", " ")}</title>
        <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
        <script id="MathJax-script" async
                src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
        </script>
        <style>
        body {{
            font-size: 1.2em;
        }}
        .MathJax {{
            font-size: 1.6em !important;
        }}
        a {{
            text-decoration: none;
            display: inline-block;
            padding: 8px 16px;
            background-color: rgba(0,0,0,0.75);
            color: white;
            border-radius: 7px;
            font-size: 0.8em;
        }}
        .collapsible {{
            background-color: #555;
            color: white;
            cursor: pointer;
            padding: 18px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 15px;
        }}

        .active, .collapsible:hover {{
            background-color: rgba(0,0,0,0.75);
        }}

        .content {{
            padding: 0 18px;
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.2s ease-out;
            background-color: #f1f1f1;
            overflow-x: scroll;
        }}
        </style>
        </head>
        <body>
        <p>
        <h2>{example.capitalize().replace("_", " ")}</h2>
        <div>
        {prev}
        <a href='chapter_summary.html'>Chapter Home</a>
        <a href='../summary.html'>Home</a>
        {nxt}
        </div>
        <h4>Equation</h4>
            \\({eq}\\) <br>
        <h4>Solution</h4>
            <!-- Render latex if solution was found by dsolve, else render error message -->
            {f"{backslash}({sol}{backslash})" if status != 1 else sol} <br>
            <h4>Verification (using checkodesol)</h4> {checkodesol_output} <br>
        <h4>Time Taken</h4>
            {time} seconds <br>
        </p>
        {hints_section}
        <script>
        var coll = document.getElementsByClassName("collapsible");
        var i;

        for (i = 0; i < coll.length; i++)
        coll[i].addEventListener("click", function() {{
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.maxHeight)
            content.style.maxHeight = null;
            else
            content.style.maxHeight = content.scrollHeight + "px";
        }});
        </script>
        </body>
        </html>
        """
        os.makedirs(f"kamke/chapter_{chno}/", exist_ok=True)
        file = open(f"kamke/chapter_{chno}/{example}.html", "w")
        file.write(example_page)
        file.close()


    def create_chapter_page(self, chno, rows):
        chapter_page = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width">
        <title>Chapter {chno}</title>
        <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
        <script id="MathJax-script" async
                src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
        </script>
        <style>
        {self.css}
        .col-1, .col-2, .col-3, .col-4 {{
            flex-basis: 10%;
        }}
        .col-5 {{
            flex-basis: 50%;
        }}
        .col-6 {{
            flex-basis: 10%;
        }}
        #home {{
            text-decoration: none;
            display: inline-block;
            padding: 8px 16px;
            background-color: rgba(0,0,0,0.75);
            color: white;
            border-radius: 7px;
        }}
        </style>
        </head>
        <body>
        <h1>Kamke Test Suite</h1>
        <h2>Chapter {chno}</h2>
        <div style="text-align: center">
            <a id="home" href="../summary.html">Home</a>
        </div>
        <div class="container">
            <ul>
                <li class="table-header">
                    <div class="col-1">Name</div>
                    <div class="col-2">Order</div>
                    <div class="col-3">Linearity</div>
                    <div class="col-4">Status</div>
                    <div class="col-5">Hint</div>
                    <div class="col-6">Time</div>
                </li>
                {rows}
            </ul>
        </div>
        </body>
        </html>
        """
        file = open(f"kamke/chapter_{chno}/chapter_summary.html", "w")
        file.write(chapter_page)
        file.close()


    def create_main_page(self, rows):
        main_page = f"""<!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width">
        <title>Kamke Test Suite</title>
        <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
        <script id="MathJax-script" async
                src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
        </script>
        <style>
        {self.css}
        .col-1 {{
            flex-basis: 15%;
        }}
        .col-2, .col-3, .col-4 {{
            flex-basis: 20%;
        }}
        .col-5 {{
            flex-basis: 25%;
        }}
        </style>
        </head>
        <body>
        <h1>Kamke Test Suite</h1>
        <h2>Chapter Summary</h2>
        <div class="container">
            <ul>
                <li class="table-header">
                    <div class="col-1">Chapter</div>
                    <div class="col-2">Solved</div>
                    <div class="col-3">Failed</div>
                    <div class="col-4">Solved and Checked</div>
                    <div class="col-5">Time</div>
                </li>
                {rows}
            </ul>
        </div>
        </body>
        </html>
        """
        file = open(f"kamke/summary.html", "w")
        file.write(main_page)
        file.close()


    def get_linearity(self, eq):
        if isinstance(eq, tuple):
            powers = set()
            for e in eq:
                powers |= e.atoms(Pow)
        else:
            powers = eq.atoms(Pow)
        for term in powers:
            if isinstance(term.base, AppliedUndef):
                if term.base.func == y:
                    return "Non-Linear"
        return "Linear"


    def get_order(self, eq):
        if isinstance(eq, tuple):
            eq = eq[0][0]
        order = ode_order(eq, y(x))
        if order == 1:
            return "1st"
        if order == 2:
            return "2nd"
        if order == 3:
            return "3rd"
        return f"{order}th"


    def get_example(self, example):
        return self.all_chapters[int(example[6])-1][example]


    def test_example(self, example, hint="default", verify=False, dsolve_time=10, checkodesol_time=10, single=False, all_hints=False):
        start = time.time()
        eq = self.get_example(example)
        final_sol = None
        # Status index
        final_status = 1
        status_messages = ["Solved", "Failed", "Solved and Checked"]
        final_error_message = ""
        final_checkodesol_output = "Skipped"
        classify_output = ""
        classify_hints = classify_ode(eq, y(x))

        if not all_hints:
            hints = [hint]
            for hnt in classify_hints:
                classify_output += f"<li>{hnt}</li>\n"
        else:
            hints = classify_hints

        for hnt in hints:
            sol = None
            status = 1
            error_message = ""
            checkodesol_output = "Skipped"
            hint_start = time.time()
            try:
                # Try to find the solution to the equation
                with time_limit(dsolve_time, 'dsolve'):
                    if isinstance(eq, tuple):
                        dsolve(*eq)
                    else:
                        sol = dsolve(eq, y(x), hnt)
                    if isinstance(sol, Eq):
                        sol = [sol]
                    status = 0

            except (TimeOutError, ValueError, NotImplementedError, TypeError) as e:
                # Solution not found / timeout
                error_message += str(e) + "\n"

            # If a solution is found
            if sol is not None and verify:
                try:
                    # Try to verify if the solution is correct
                    assert len(sol)
                    with time_limit(checkodesol_time, 'checkodesol'):
                        if isinstance(eq, tuple):
                            checkodesol_output = checkodesol(eq[0], sol)
                        else:
                            checkodesol_output = checkodesol(eq, sol, y(x))
                    assert any([x[0] for x in checkodesol_output])
                    status = 2

                except AssertionError as e:
                    # Wrong solution
                    error_message += str(e) + "\n"

                except (TimeOutError, ValueError, NotImplementedError, TypeError) as e:
                    # Checkodesol unable to verify / timeout
                    error_message += str(e) + "\n"

            if sol is not None:
                if final_sol is None:
                    final_sol = sol
                    final_status = status
                    final_error_message = error_message
                    final_checkodesol_output = checkodesol_output
                hint_sol = f"\\({latex(sol)}\\)"
            else:
                hint_sol = error_message

            if hnt != "default":
                classify_output += f"""<button class="collapsible">{hnt}</button>
                        <div class="content">
                            <h4>Solution</h4>
                            {hint_sol} <br>
                            <h4>Verification (using checkodesol)</h4> {checkodesol_output} <br>
                            <h4>Time Taken</h4>
                            {time.time() - hint_start} seconds <br> <br>
                        </div>"""

        elapsed = time.time() - start
        log = f"{example} {status_messages[final_status]} in {elapsed} seconds\n"

        if final_sol is None:
            if len(classify_hints) == 0:
                final_error_message += "Not Implemented\n"
            log += final_error_message
            self.create_example_page(example, latex(eq), final_status, final_error_message, elapsed, classify_output, final_checkodesol_output, all_hints)
        else:
            log += f"Equation: {eq}\nSolution: {final_sol}\n"
            self.create_example_page(example, latex(eq), final_status, latex(final_sol), elapsed, classify_output, final_checkodesol_output, all_hints)
        if single:
            print(log)
        return [log, final_status, classify_hints, elapsed]


    def test_chapter(self, chno, hint="default", verify=False, dsolve_time=10, checkodesol_time=10, single=False, all_hints=False):
        # Time elapsed
        total_time = 0
        # Counts for result
        counts = [0, 0, 0]
        status_messages = ["Solved", "Failed", "Solved and Checked"]
        rows = ""
        os.makedirs(f"kamke/chapter_{chno}/", exist_ok=True)

        for example in list(self.all_chapters[chno-1])[:2]:
            eq = self.get_example(example)
            if isinstance(eq, tuple):
                eq = eq[0]
            rows += f"""<li class="table-row">
            <div class="col-1">
                <a href='{example}.html'>{example.capitalize().replace("_", " ")}</a>
            </div>
            <div class="col-2">
                {self.get_order(eq)}
            </div>
            <div class="col-3">
                {self.get_linearity(eq)}
            </div>\n"""
            output, status, hints, elapsed = self.test_example(example, hint, verify, dsolve_time, checkodesol_time, all_hints=all_hints)
            print(output)
            counts[status] += 1
            total_time += elapsed
            hint = hints[0] if len(hints) else "Not Implemented"
            rows += f"""<div class="col-4">
            {status_messages[status]}
            </div>
            <div class="col-5">
                {hint}
            </div>
            <div class="col-6">
                {round(elapsed, 3)}
            </div>
            </li>\n"""

        self.create_chapter_page(chno, rows)

        # Summary
        if single:
            print("Total time taken:", total_time)
            print("No. of ODEs solved:", counts[0])
            print("No. of ODEs failed:", counts[1])
            print("No. of ODEs solved and checked:", counts[2])
        return [total_time, counts]


    def test_all_examples(self, hint="default", verify=False, dsolve_time=10, checkodesol_time=10, all_hints=False):
        # Time elapsed
        total_time = 0
        # Counts for result
        counts = [0, 0, 0]
        rows = ""

        for chapter in range(1, 8):
            elapsed, cts = self.test_chapter(chapter, hint, verify, dsolve_time, checkodesol_time, all_hints=all_hints)
            total_time += elapsed
            rows += f"""<li class="table-row">
            <div class="col-1">
                <a href="chapter_{chapter}/chapter_summary.html">Chapter {chapter}</a>
            </div>
            <div class="col-2">
                {cts[0]}
            </div>
            <div class="col-3">
                {cts[1]}
            </div>
            <div class="col-4">
                {cts[2]}
            </div>
            <div class="col-5">
                {round(elapsed, 3)}
            </div>
            </li>
            """
            for i in range(3):
                counts[i] += cts[i]

        self.create_main_page(rows)

        # Summary
        print("Total time taken:", total_time)
        print("No. of ODEs solved:", counts[0])
        print("No. of ODEs failed:", counts[1])
        print("No. of ODEs solved and checked:", counts[2])

if __name__ == '__main__':
    os.makedirs("kamke", exist_ok=True)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--example", help="Name of the example in the format kamke_{chapter_no}.{problem_no}\nSpecify all to test all examples", default="all")
    parser.add_argument("-ch", "--chapter", help="Chapter no. Tests all examples of a chapter", type=int)
    parser.add_argument("--hint", help="Hint to be used to solve the ODEs", default="default")
    parser.add_argument("--all_hints", help="Solve the ODE with all matching hints", action="store_true")
    parser.add_argument("--verify", help="Verify the solution from dsolve using checkodesol", action="store_true")
    parser.add_argument("--dsolve_time", help="Timeout duration (in seconds) for dsolve", type=int, default=10)
    parser.add_argument("--checkodesol_time", help="Timeout duration (in seconds) for checkodesol", type=int, default=10)

    args = parser.parse_args()
    kamke = Kamke()
    if args.chapter is not None:
        kamke.test_chapter(args.chapter, args.hint, args.verify, args.dsolve_time, args.checkodesol_time, True, args.all_hints)
    elif args.example == "all":
        kamke.test_all_examples(args.hint, args.verify, args.dsolve_time, args.checkodesol_time, all_hints=args.all_hints)
    else:
        kamke.test_example(args.example, args.hint, args.verify, args.dsolve_time, args.checkodesol_time, True, args.all_hints)
