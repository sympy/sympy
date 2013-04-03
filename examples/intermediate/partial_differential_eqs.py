#!/usr/bin/env python

"""Partial Differential Equations example

Demonstrates various ways to solve partial differential equations
"""

from sympy import symbols, Eq, Function, pde_separate, pprint, sin, cos, latex
from sympy import Derivative as D


def main():
    r, phi, theta = symbols("r,phi,theta")
    Xi = Function('Xi')
    R, Phi, Theta, u = map(Function, ['R', 'Phi', 'Theta', 'u'])
    C1, C2 = symbols('C1,C2')

    pprint("Separation of variables in Laplace equation in spherical coordinates")
    pprint("Laplace equation in spherical coordinates:")
    eq = Eq(D(Xi(r, phi, theta), r, 2) + 2/r * D(Xi(r, phi, theta), r) +
            1/(r**2 * sin(phi)**2) * D(Xi(r, phi, theta), theta, 2) +
            cos(phi)/(r**2 * sin(phi)) * D(Xi(r, phi, theta), phi) +
            1/r**2 * D(Xi(r, phi, theta), phi, 2))
    pprint(eq)

    pprint("We can either separate this equation in regards with variable r:")
    res_r = pde_separate(eq, Xi(r, phi, theta), [R(r), u(phi, theta)])
    pprint(res_r)

    pprint("Or separate it in regards of theta:")
    res_theta = pde_separate(eq, Xi(r, phi, theta), [Theta(theta), u(r, phi)])
    pprint(res_theta)

    res_phi = pde_separate(eq, Xi(r, phi, theta), [Phi(phi), u(r, theta)])
    pprint("But we cannot separate it in regards of variable phi: ")
    pprint("Result: %s" % res_phi)

    pprint("\n\nSo let's make theta dependent part equal with -C1:")
    eq_theta = Eq(res_theta[0], -C1)

    pprint(eq_theta)
    pprint("\nThis also means that second part is also equal to -C1:")
    eq_left = Eq(res_theta[1], -C1)
    pprint(eq_left)

    pprint("\nLets try to separate phi again :)")
    res_theta = pde_separate(eq_left, u(r, phi), [Phi(phi), R(r)])
    pprint("\nThis time it is successful:")
    pprint(res_theta)

    pprint("\n\nSo our final equations with separated variables are:")
    pprint(eq_theta)
    pprint(Eq(res_theta[0], C2))
    pprint(Eq(res_theta[1], C2))


if __name__ == "__main__":
    main()
