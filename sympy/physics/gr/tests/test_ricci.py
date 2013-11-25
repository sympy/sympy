from sympy.tensor.components import Index, Metric
from sympy.physics.gr.ricci import christoffel, riemann_tensor, ricci_curvature_tensor, ricci_scalar
from sympy import symbols, sin, cos

def test_ricci():

    # the curvature on the surface of a sphere
    r, theta, phi = symbols("r theta phi")
    g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
    gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, [theta,phi])
    rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
    ricci = ricci_scalar(g, rab)
    assert ricci == 2.0/r**2

    # the ricci scalar of the Schwarzschild metric
    t, r, rs, theta, phi = symbols("t r rs theta phi")
    coordinates = [t, r, theta, phi]
    g = Metric([Index("_", "i"), Index("_", "j")], [[(1-rs/r), 0, 0, 0],
                                                    [0, -1/(1-rs/r), 0, 0],
                                                    [0, 0, -r**2, 0],
                                                    [0, 0, 0, - r**2 * sin(theta)**2]])
    gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, coordinates)
    rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, coordinates)
    rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
    ricci = ricci_scalar(g, rab)
    assert ricci == 0
