#!/usr/bin/env python

"""
This example shows how to work with coordinate transformations, curvilinear
coordinates and a little bit with differential geometry.

It takes polar, cylindrical, spherical, rotating disk coordinates and others
and calculates all kinds of interesting properties, like Jacobian, metric
tensor, Laplace operator, ...
"""

from sympy import var, sin, cos, pprint, Matrix, eye, trigsimp, Eq, \
    Function, simplify, sinh, cosh, expand, symbols


def laplace(f, g_inv, g_det, X):
    """
    Calculates Laplace(f), using the inverse metric g_inv, the determinant of
    the metric g_det, all in variables X.
    """
    r = 0
    for i in range(len(X)):
        for j in range(len(X)):
            r += g_inv[i, j]*f.diff(X[i]).diff(X[j])
    for sigma in range(len(X)):
        for alpha in range(len(X)):
            r += g_det.diff(X[sigma]) * g_inv[sigma, alpha] * \
                f.diff(X[alpha]) / (2*g_det)
    return r


def transform(name, X, Y, g_correct=None, recursive=False):
    """
    Transforms from cartesian coordinates X to any curvilinear coordinates Y.

    It printing useful information, like Jacobian, metric tensor, determinant
    of metric, Laplace operator in the new coordinates, ...

    g_correct ... if not None, it will be taken as the metric --- this is
                  useful if sympy's trigsimp() is not powerful enough to
                  simplify the metric so that it is usable for later
                  calculation. Leave it as None, only if the metric that
                  transform() prints is not simplified, you can help it by
                  specifying the correct one.

    recursive ... apply recursive trigonometric simplification (use only when
                  needed, as it is an expensive operation)
    """
    print("_"*80)
    print("Transformation:", name)
    for x, y in zip(X, Y):
        pprint(Eq(y, x))
    J = X.jacobian(Y)
    print("Jacobian:")
    pprint(J)
    g = J.T*eye(J.shape[0])*J

    g = g.applyfunc(expand)
    print("metric tensor g_{ij}:")
    pprint(g)
    if g_correct is not None:
        g = g_correct
        print("metric tensor g_{ij} specified by hand:")
        pprint(g)
    print("inverse metric tensor g^{ij}:")
    g_inv = g.inv(method="ADJ")
    g_inv = g_inv.applyfunc(simplify)
    pprint(g_inv)
    print("det g_{ij}:")
    g_det = g.det()
    pprint(g_det)
    f = Function("f")(*list(Y))
    print("Laplace:")
    pprint(laplace(f, g_inv, g_det, Y))


def main():
    mu, nu, rho, theta, phi, sigma, tau, a, t, x, y, z, w = symbols(
        "mu, nu, rho, theta, phi, sigma, tau, a, t, x, y, z, w")

    transform("polar", Matrix([rho*cos(phi), rho*sin(phi)]), [rho, phi])

    transform("cylindrical", Matrix([rho*cos(phi), rho*sin(phi), z]),
              [rho, phi, z])

    transform("spherical",
              Matrix([rho*sin(theta)*cos(phi), rho*sin(theta)*sin(phi),
                      rho*cos(theta)]),
              [rho, theta, phi],
              recursive=True
              )

    transform("rotating disk",
              Matrix([t,
                      x*cos(w*t) - y*sin(w*t),
                      x*sin(w*t) + y*cos(w*t),
                      z]),
              [t, x, y, z])

    transform("parabolic",
              Matrix([sigma*tau, (tau**2 - sigma**2) / 2]),
              [sigma, tau])

    transform("bipolar",
            Matrix([a*sinh(tau)/(cosh(tau)-cos(sigma)),
                a*sin(sigma)/(cosh(tau)-cos(sigma))]),
            [sigma, tau]
            )

    transform("elliptic",
              Matrix([a*cosh(mu)*cos(nu), a*sinh(mu)*sin(nu)]),
              [mu, nu]
              )

if __name__ == "__main__":
    main()
