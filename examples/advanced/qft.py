#!/usr/bin/env python

"""Quantum field theory example

* http://en.wikipedia.org/wiki/Quantum_field_theory

This particular example is a work in progress. Currently it calculates the
scattering amplitude of the process:

    electron + positron -> photon -> electron + positron

in QED (http://en.wikipedia.org/wiki/Quantum_electrodynamics). The aim
is to be able to do any kind of calculations in QED or standard model in
SymPy, but that's a long journey.

"""

from __future__ import division, print_function

from sympy import Basic, exp, Symbol, sin, Rational, I, Mul, Matrix, \
    ones, sqrt, pprint, simplify, Eq, sympify

from sympy.physics import msigma, mgamma

# gamma^mu
gamma0 = mgamma(0)
gamma1 = mgamma(1)
gamma2 = mgamma(2)
gamma3 = mgamma(3)
gamma5 = mgamma(5)

# sigma_i
sigma1 = msigma(1)
sigma2 = msigma(2)
sigma3 = msigma(3)

E = Symbol("E", real=True)
m = Symbol("m", real=True)


def u(p, r):
    """ p = (p1, p2, p3); r = 0,1 """
    if r not in [1, 2]:
        raise ValueError("Value of r should lie between 1 and 2")
    p1, p2, p3 = p
    if r == 1:
        ksi = Matrix([[1], [0]])
    else:
        ksi = Matrix([[0], [1]])
    a = (sigma1*p1 + sigma2*p2 + sigma3*p3) / (E + m)*ksi
    if a == 0:
        a = zeros(2, 1)
    return sqrt(E + m) *\
        Matrix([[ksi[0, 0]], [ksi[1, 0]], [a[0, 0]], [a[1, 0]]])


def v(p, r):
    """ p = (p1, p2, p3); r = 0,1 """
    if r not in [1, 2]:
        raise ValueError("Value of r should lie between 1 and 2")
    p1, p2, p3 = p
    if r == 1:
        ksi = Matrix([[1], [0]])
    else:
        ksi = -Matrix([[0], [1]])
    a = (sigma1*p1 + sigma2*p2 + sigma3*p3) / (E + m)*ksi
    if a == 0:
        a = zeros(2, 1)
    return sqrt(E + m) *\
        Matrix([[a[0, 0]], [a[1, 0]], [ksi[0, 0]], [ksi[1, 0]]])


def pslash(p):
    p1, p2, p3 = p
    p0 = sqrt(m**2 + p1**2 + p2**2 + p3**2)
    return gamma0*p0 - gamma1*p1 - gamma2*p2 - gamma3*p3


def Tr(M):
    return M.trace()


def xprint(lhs, rhs):
    pprint(Eq(sympify(lhs), rhs))


def main():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)
    c = Symbol("c", real=True)

    p = (a, b, c)

    assert u(p, 1).D*u(p, 2) == Matrix(1, 1, [0])
    assert u(p, 2).D*u(p, 1) == Matrix(1, 1, [0])

    p1, p2, p3 = [Symbol(x, real=True) for x in ["p1", "p2", "p3"]]
    pp1, pp2, pp3 = [Symbol(x, real=True) for x in ["pp1", "pp2", "pp3"]]
    k1, k2, k3 = [Symbol(x, real=True) for x in ["k1", "k2", "k3"]]
    kp1, kp2, kp3 = [Symbol(x, real=True) for x in ["kp1", "kp2", "kp3"]]

    p = (p1, p2, p3)
    pp = (pp1, pp2, pp3)

    k = (k1, k2, k3)
    kp = (kp1, kp2, kp3)

    mu = Symbol("mu")

    e = (pslash(p) + m*ones(4))*(pslash(k) - m*ones(4))
    f = pslash(p) + m*ones(4)
    g = pslash(p) - m*ones(4)

    xprint('Tr(f*g)', Tr(f*g))

    M0 = [(v(pp, 1).D*mgamma(mu)*u(p, 1))*(u(k, 1).D*mgamma(mu, True) *
                                                 v(kp, 1)) for mu in range(4)]
    M = M0[0] + M0[1] + M0[2] + M0[3]
    M = M[0]
    if not isinstance(M, Basic):
        raise TypeError("Invalid type of variable")

    d = Symbol("d", real=True)  # d=E+m

    xprint('M', M)
    print("-"*40)
    M = ((M.subs(E, d - m)).expand()*d**2).expand()
    xprint('M2', 1 / (E + m)**2*M)
    print("-"*40)
    x, y = M.as_real_imag()
    xprint('Re(M)', x)
    xprint('Im(M)', y)
    e = x**2 + y**2
    xprint('abs(M)**2', e)
    print("-"*40)
    xprint('Expand(abs(M)**2)', e.expand())

if __name__ == "__main__":
    main()
