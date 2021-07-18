"""
File to test Kamke ODEs. To run test suite, run

    python test_kamke.py

Information about failing ODEs is logged in test_report.txt.
"""


from contextlib import contextmanager
import threading
import _thread
import time

from sympy import (symbols, S, Function, Eq, dsolve, checkodesol,
    Abs, sqrt, cos, exp, log, sin, tan, Derivative)


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
        raise TimeOutError("Timed out for operation {}".format(msg))
    finally:
        # Cancel timer if process finishes in time
        timer.cancel()

A, B, a, a0, a1, a2, a3, a4, alpha, b, b0, b1, b2, b3, b4, bbeta, c, d, ggamma, k, m, n, nu, x = symbols('A, B, a, a0, a1, a2, a3, a4, alpha, b, b0, b1, b2, b3, b4, bbeta, c, d, ggamma, k, m, n, nu, x')
R1, R2, f, f0, f1, f2, f3, g, g0, g1, h, phi, tg, y = symbols('R1, R2, f, f0, f1, f2, f3, g, g0, g1, h, phi, tg, y', cls=Function)

class Kamke:

    eqs = {
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
        "kamke_1.14": a*x**m + y(x)**2 + Derivative(y(x), x),
        "kamke_1.15": x**4 - 2*x**2*y(x) - 2*x + y(x)**2 + Derivative(y(x), x) - 1,
        "kamke_1.16": (x*y(x) - 1)*f(x) + y(x)**2 + Derivative(y(x), x),
        "kamke_1.17": -y(x)**2 - 3*y(x) + Derivative(y(x), x) + 4,
        "kamke_1.18": -x*y(x) - x - y(x)**2 + Derivative(y(x), x) + 1,
        "kamke_1.19": -(x + y(x))**2 + Derivative(y(x), x),
        "kamke_1.20": -2*x + (x**2 + 1)*y(x) - y(x)**2 + Derivative(y(x), x),
        "kamke_1.21": -y(x)**2 + y(x)*sin(x) - cos(x) + Derivative(y(x), x),
        "kamke_1.22": -y(x)**2 - y(x)*sin(2*x) - cos(2*x) + Derivative(y(x), x),
        "kamke_1.23": a*y(x)**2 - b + Derivative(y(x), x),
        "kamke_1.24": a*y(x)**2 - b*x**nu + Derivative(y(x), x),
        "kamke_1.25": a*y(x)**2 - b*x**(2*nu) - c*x**(nu - 1) + Derivative(y(x), x),
        "kamke_1.26": -(A*y(x) - a)*(B*y(x) - b) + Derivative(y(x), x),
        "kamke_1.27": a*(-x + y(x))*y(x) + Derivative(y(x), x) - 1,
        "kamke_1.28": -x**3*y(x) + x*y(x)**2 - 2*x + Derivative(y(x), x),
        "kamke_1.29": -x*y(x)**2 - 3*x*y(x) + Derivative(y(x), x),
        "kamke_1.30": -x**a + x**(-a - 1)*y(x)**2 + Derivative(y(x), x),
        "kamke_1.31": -a*x**n*(y(x)**2 + 1) + Derivative(y(x), x),
        "kamke_1.32": y(x)**2*sin(x) - 2*sin(x)/cos(x)**2 + Derivative(y(x), x),
        "kamke_1.33": Derivative(y(x), x) - y(x)**2*Derivative(f(x), x)/g(x) + Derivative(g(x), x)/f(x),
        "kamke_1.34": f(x)*y(x)**2 + g(x)*y(x) + Derivative(y(x), x),
        "kamke_1.35": (2*a*y(x) + b + y(x)**2)*f(x) + Derivative(y(x), x),
        "kamke_1.36": a*x*y(x)**2 + y(x)**3 + Derivative(y(x), x),
        "kamke_1.37": -a*y(x)**2*exp(x) - y(x)**3 + Derivative(y(x), x),
        "kamke_1.38": -a*y(x)**3 - b/x**(S(3)/2) + Derivative(y(x), x),
        "kamke_1.39": -a0 - a1*y(x) - a2*y(x)**2 - a3*y(x)**3 + Derivative(y(x), x),
        "kamke_1.40": 6*a*x*y(x)**2 + 3*a*y(x)**3 + Derivative(y(x), x),
        "kamke_1.41": a*x*y(x)**3 + b*y(x)**2 + Derivative(y(x), x),
        "kamke_1.42": -x*(x + 2)*y(x)**3 - (x + 3)*y(x)**2 + Derivative(y(x), x),
        "kamke_1.43": 3*x*y(x)**2 + (4*a**2*x + 3*a*x**2 + b)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.44": 2*a*x**3*y(x)**3 + 2*x*y(x) + Derivative(y(x), x),
        "kamke_1.45": 3*b*y(x)**2 + (2*a**2*x**3 - 2*b**2*x)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.46": a*x**(-a - 1) - x**a*y(x)**3 + 3*y(x)**2 + Derivative(y(x), x) - y(x)/x**a - 1/x**(2*a),
        "kamke_1.47": -a*(-x + x**n)*y(x)**3 - y(x)**2 + Derivative(y(x), x),
        "kamke_1.48": -c*y(x)**2 - (a*x**n + b*x)*y(x)**3 + Derivative(y(x), x),
        "kamke_1.49": 6*a*phi(x)*y(x)**2 + a*y(x)**3*Derivative(phi(x), x) + 2*a + (2*a + 1)*y(x)*Derivative(phi(x), (x, 2))/Derivative(phi(x), x) + Derivative(y(x), x) + 2,
        "kamke_1.50": -f0(x) - f1(x)*y(x) - f2(x)*y(x)**2 - f3(x)*y(x)**3 + Derivative(y(x), x),
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
        "kamke_1.63": Derivative(y(x), x) - (y(x)**2 + 1)/((x + 1)**(S(3)/2)*Abs(sqrt(y(x) + 1) + y(x))),
        "kamke_1.64": -sqrt((a*y(x)**2 + b*y(x) + c)/(a*x**2 + b*x + c)) + Derivative(y(x), x),
        "kamke_1.65": -sqrt((y(x)**3 + 1)/(x**3 + 1)) + Derivative(y(x), x),
        "kamke_1.66": Derivative(y(x), x) - sqrt(Abs((1 - y(x))*(-a*y(x) + 1)*y(x)))/sqrt(Abs(x*(1 - x)*(-a*x + 1))),
        "kamke_1.67": Derivative(y(x), x) - sqrt(1 - y(x)**4)/sqrt(1 - x**4),
        "kamke_1.68": -sqrt((a*y(x)**4 + b*y(x)**2 + 1)/(a*x**4 + b*x**2 + 1)) + Derivative(y(x), x),
        "kamke_1.69": -sqrt((a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)*(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.70": -sqrt((a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)/(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.71": -sqrt((b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)/(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)) + Derivative(y(x), x),
        "kamke_1.72": -R1(x, sqrt(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4))*R2(y(x), sqrt(b0 + b1*y(x) + b2*y(x)**2 + b3*y(x)**3 + b4*y(x)**4)) + Derivative(y(x), x),
        "kamke_1.73": -((a0 + a1*x + a2*x**2 + a3*x**3)/(a0 + a1*y(x) + a2*y(x)**2 + a3*y(x)**3))**(S(2)/3) + Derivative(y(x), x),
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
        "kamke_1.94": a*y(x) + b*x**n + x*Derivative(y(x), x),
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
        "kamke_1.113": a*sqrt(x**2 + y(x)**2) + x*Derivative(y(x), x) - y(x),
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
        "kamke_1.205": a*y(x) + b*x**n + x*(a**2/4 - S(1)/4) + y(x)*Derivative(y(x), x),
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
        "kamke_1.219": (g(x) + y(x))*Derivative(y(x), x) - f0(x) - f1(x)*y(x) - f2(x)*y(x)**2,
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
        "kamke_1.269": (g0(x) + g1(x)*y(x))*Derivative(y(x), x) - f0(x) - f1(x)*y(x) - f2(x)*y(x)**2 - f3(x)*y(x)**3,
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
        "kamke_1.300": 6*x*y(x)**2*Derivative(y(x), x) + x + 2*y(x)**3
    }

    def get_example(self, example):
        return self.eqs[example]

    def test_example(self, example, hint="default", verify=True, dsolve_timelim=10, checkodesol_timelim=10, log=False):
        start = time.time()
        eq = self.get_example(example)
        sol = None
        output = example + " "

        try:
            # Try to find the solution to the equation
            with time_limit(dsolve_timelim):
                sol = dsolve(eq, y(x), hint)
                if isinstance(sol, Eq):
                    sol = [sol]
        except (TimeOutError, ValueError, NotImplementedError, TypeError):
            # Solution not found / timeout
            output += "dsolve failed"

        # If a solution is found
        if sol is not None and verify:
            try:
                # Try to verify if the solution is correct
                assert len(sol)
                with time_limit(checkodesol_timelim):
                    checks = checkodesol(eq, sol, y(x))
                assert any([x[0] for x in checks])
                output += "solved"
            except AssertionError:
                # Wrong solution
                output += "dsolve failed"
            except (TimeOutError, ValueError, NotImplementedError, TypeError):
                # Checkodesol unable to verify / timeout
                output += "checkodesol failed"
        output += f" in {time.time() - start} seconds\n"
        returns = sol
        if log:
            returns = [sol, output]
        return returns


    def test_all_examples(self, hint="default", verify=True, dsolve_timelim=10, checkodesol_timelim=10, report_path=None):
        if report_path is None:
            report_path = "test_report.txt"

        start = time.time()
        # Counts for result
        solved = 0
        check_slow = 0
        fail = 0

        with open(report_path, "w") as output_file:
            for example in self.eqs:
                # Process is not terminating despite using
                # time limit. Checkodesol is too slow for
                # these cases. This should be fixed.
                if example in ["kamke_1.14", "kamke_1.24", "kamke_1.31", "kamke_1.47", \
                "kamke_1.48", "kamke_1.94", "kamke_1.113", "kamke_1.205"]:
                    continue

                sol, output = self.test_example(example, hint, verify, dsolve_timelim, checkodesol_timelim, log=True)
                output_file.write(output)
                output_file.flush()

            # Summary
            output_file.write("Total time taken" + str(time.time() - start))
            output_file.write("No. of ODEs solved:" + str(solved))
            output_file.write("No. of ODEs failed:" + str(fail))
            output_file.write("No. of ODEs solved, but not checked:" + str(check_slow))
            output_file.close()
