from sympy.tensor.components import *
from sympy.physics.gr.ricci import *
from sympy import symbols, sin, cos

r, theta, phi = symbols("r theta phi")
g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, [theta,phi])
rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
r = ricci_scalar(g, rab)

print(r) # gives 2.0/r**2, which is the curvature on the surface of a sphere

t, r, rs, theta, phi = symbols("t r rs theta phi")
coordinates = [t, r, theta, phi]
g = Metric([Index("_", "i"), Index("_", "j")], [[(1-rs/r), 0, 0, 0],
                                                [0, -1/(1-rs/r), 0, 0],
                                                [0, 0, -r**2, 0],
                                                [0, 0, 0, - r**2 * sin(theta)**2]])
gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, coordinates)
rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, coordinates)
rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
r = ricci_scalar(g, rab)

print(r) # gives 0, which is the ricci scalar of the Schwarzschild metric
