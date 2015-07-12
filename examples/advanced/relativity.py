#!/usr/bin/env python

"""
This example calculates the Ricci tensor from the metric and does this
on the example of Schwarzschild solution.

If you want to derive this by hand, follow the wiki page here:

http://en.wikipedia.org/wiki/Deriving_the_Schwarzschild_solution

Also read the above wiki and follow the references from there if
something is not clear, like what the Ricci tensor is, etc.

"""

from __future__ import division, print_function

from sympy import (exp, Symbol, sin, Rational, Derivative, dsolve, Function,
                  Matrix, Eq, pprint, Pow, classify_ode, solve)


def grad(f, X):
    a = []
    for x in X:
        a.append(f.diff(x))
    return a


def d(m, x):
    return grad(m[0, 0], x)


class MT(object):
    def __init__(self, m):
        self.gdd = m
        self.guu = m.inv()

    def __str__(self):
        return "g_dd =\n" + str(self.gdd)

    def dd(self, i, j):
        return self.gdd[i, j]

    def uu(self, i, j):
        return self.guu[i, j]


class G(object):
    def __init__(self, g, x):
        self.g = g
        self.x = x

    def udd(self, i, k, l):
        g = self.g
        x = self.x
        r = 0
        for m in [0, 1, 2, 3]:
            r += g.uu(i, m)/2 * (g.dd(m, k).diff(x[l]) + g.dd(m, l).diff(x[k])
                    - g.dd(k, l).diff(x[m]))
        return r


class Riemann(object):
    def __init__(self, G, x):
        self.G = G
        self.x = x

    def uddd(self, rho, sigma, mu, nu):
        G = self.G
        x = self.x
        r = G.udd(rho, nu, sigma).diff(x[mu]) - G.udd(rho, mu, sigma).diff(x[nu])
        for lam in [0, 1, 2, 3]:
            r += G.udd(rho, mu, lam)*G.udd(lam, nu, sigma) \
                - G.udd(rho, nu, lam)*G.udd(lam, mu, sigma)
        return r


class Ricci(object):
    def __init__(self, R, x):
        self.R = R
        self.x = x
        self.g = R.G.g

    def dd(self, mu, nu):
        R = self.R
        x = self.x
        r = 0
        for lam in [0, 1, 2, 3]:
            r += R.uddd(lam, mu, lam, nu)
        return r

    def ud(self, mu, nu):
        r = 0
        for lam in [0, 1, 2, 3]:
            r += self.g.uu(mu, lam)*self.dd(lam, nu)
        return r.expand()


def curvature(Rmn):
    return Rmn.ud(0, 0) + Rmn.ud(1, 1) + Rmn.ud(2, 2) + Rmn.ud(3, 3)

nu = Function("nu")
lam = Function("lambda")

t = Symbol("t")
r = Symbol("r")
theta = Symbol(r"theta")
phi = Symbol(r"phi")

# general, spherically symmetric metric
gdd = Matrix((
    (-exp(nu(r)), 0, 0, 0),
    (0, exp(lam(r)), 0, 0),
    (0, 0, r**2, 0),
    (0, 0, 0, r**2*sin(theta)**2)
))
g = MT(gdd)
X = (t, r, theta, phi)
Gamma = G(g, X)
Rmn = Ricci(Riemann(Gamma, X), X)


def pprint_Gamma_udd(i, k, l):
    pprint(Eq(Symbol('Gamma^%i_%i%i' % (i, k, l)), Gamma.udd(i, k, l)))


def pprint_Rmn_dd(i, j):
    pprint(Eq(Symbol('R_%i%i' % (i, j)), Rmn.dd(i, j)))


# from Differential Equations example
def eq1():
    r = Symbol("r")
    e = Rmn.dd(0, 0)
    e = e.subs(nu(r), -lam(r))
    pprint(dsolve(e, lam(r)))


def eq2():
    r = Symbol("r")
    e = Rmn.dd(1, 1)
    C = Symbol("CC")
    e = e.subs(nu(r), -lam(r))
    pprint(dsolve(e, lam(r)))


def eq3():
    r = Symbol("r")
    e = Rmn.dd(2, 2)
    e = e.subs(nu(r), -lam(r))
    pprint(dsolve(e, lam(r)))


def eq4():
    r = Symbol("r")
    e = Rmn.dd(3, 3)
    e = e.subs(nu(r), -lam(r))
    pprint(dsolve(e, lam(r)))
    pprint(dsolve(e, lam(r), 'best'))


def main():

    print("Initial metric:")
    pprint(gdd)
    print("-"*40)
    print("Christoffel symbols:")
    pprint_Gamma_udd(0, 1, 0)
    pprint_Gamma_udd(0, 0, 1)
    print()
    pprint_Gamma_udd(1, 0, 0)
    pprint_Gamma_udd(1, 1, 1)
    pprint_Gamma_udd(1, 2, 2)
    pprint_Gamma_udd(1, 3, 3)
    print()
    pprint_Gamma_udd(2, 2, 1)
    pprint_Gamma_udd(2, 1, 2)
    pprint_Gamma_udd(2, 3, 3)
    print()
    pprint_Gamma_udd(3, 2, 3)
    pprint_Gamma_udd(3, 3, 2)
    pprint_Gamma_udd(3, 1, 3)
    pprint_Gamma_udd(3, 3, 1)
    print("-"*40)
    print("Ricci tensor:")
    pprint_Rmn_dd(0, 0)
    e = Rmn.dd(1, 1)
    pprint_Rmn_dd(1, 1)
    pprint_Rmn_dd(2, 2)
    pprint_Rmn_dd(3, 3)
    print("-"*40)
    print("Solve Einstein's equations:")
    e = e.subs(nu(r), -lam(r)).doit()
    l = dsolve(e, lam(r))
    pprint(l)
    lamsol = solve(l, lam(r))[0]
    metric = gdd.subs(lam(r), lamsol).subs(nu(r), -lamsol)  # .combine()
    print("metric:")
    pprint(metric)

if __name__ == "__main__":
    main()
