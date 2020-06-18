from sympy import sin, cos, pi
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.integrals import ParametricIntegral
from sympy.vector.parametricregion import ParametricRegion
from sympy.abc import x, y, z, u, v, r, t, theta, phi, pi

C = CoordSys3D('C')

def test_lineintegrals():
    halfcircle = ParametricRegion((4*cos(theta), 4*sin(theta)), (theta, -pi/2, pi/2))
    assert ParametricIntegral(C.x*C.y**4, halfcircle) == 1638.4

    curve = ParametricRegion((t, t**2, t**3), (t, 0, 1))
    field1 = 8*C.x**2*C.y*C.z*C.i + 5*C.z*C.j - 4*C.x*C.y*C.k
    assert ParametricIntegral(field, curve) == 1
    line = ParametricRegion((4*t - 1, 2 - 2*t, t), (t, 0, 1))
    assert ParametricIntegral(C.x*C.z*C.i - C.y*C.z*C.k, line) == 3

    assert ParametricIntegral(4*R.x**3, ParametricRegion((1, t), (t, 0, 2))) == 8

    helix = ParametricRegion((cos(t), sin(t), 3*t), (t, 0, 4))
    assert ParametricIntegral(C.x*C.y*C.z, helix) == -3*sqrt(10)*pi

    field2 = C.y*C.i + C.z*C.j + C.z*C.k
    assert ParametricIntegral(field2, ParametricRegion((cos(t), sin(t), t**2), (t, 0, pi))) == 8*pi**4

def test_surfaceintegrals():

    semisphere = ParametricRegion((2*sin(phi)*cos(theta), 2*sin(phi)*sin(theta), 2*cos(phi)),\
                            (theta, 0, 2*pi), (phi, 0, pi/2))
    assert ParametricIntegral(C.z, semisphere) == 8*pi

    cylinder = ParametricRegion((sqrt(3)*cos(theta), sqrt(3)*sin(theta), z), (z, 0, 6), (theta, 0, 2*pi))
    assert ParametricIntegral(C.y, cylinder) == 0
    cone = ParametricRegion((v*cosu, v*sinu, v), (u, 0, 2*pi), (v, 0, 1))

    assert ParametricIntegral(C.x*C.i + C.y*C.z + C.z**4*C.k, cone) == pi/3
    
    triangle1 = ParametricRegion((x, y), (x, 0, 2), (y, 0, 10 - 5*x))
    triangle2 = ParametricRegion((x, y), (y, 0, 10 - 5*x), (x, 0, 2))
    assert ParametricIntegral(-15.6*C.y*C.k, traingle1) == ParametricIntegral(-15.6*C.y*C.k, traingle2)
    assert ParametricIntegral(C.z, triangle) == 10*C.z

def test_volumeintegrals():

    cube = ParametricRegion((x, y, z), (x, 0, 1), (y, 0, 1), (z, 0, 1))
    assert ParametricIntegral(1, cube) == 1
    
    solidsphere = ParametricRegion((r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)),\
                            (r, 0, 2), (theta, 0, 2*pi), (phi, 0, pi))
    assert ParametricIntegral(C.x**2 + C.y**2, solidsphere) == -256*pi/15

    region_under_plane1 = ParametricRegion((x, y, z), (x, 0, 3/2), (y, 0, -2*3/x + 2),\
                                    (z, 0, 6 - 2*C.x - 3*C,y))
    region_under_plane2 = ParametricRegion((x, y, z), (x, 0, 3/2), (z, 0, 6 - 2*C.x - 3*C,y),\
                                    (y, 0, -2*3/x + 2))
    
    assert ParametricIntegral(C.i*C.k + C.j - 100*C.k, region-under_plane1) == \
        ParametricIntegral(C.i*C.k + C.j - 100*C.k, region-under_plane2)    
    assert ParametricIntegral(2*C.x, region_under_plane) == 9
