from sympy.vector import CoordSys3D, Gradient, Divergence, Curl, VectorZero, Laplacian, gradient, directional_derivative, Vector
from sympy.printing.repr import srepr
from sympy import Symbol, S

R = CoordSys3D('R')
s1 = R.x*R.y*R.z  # type: ignore
s2 = R.x + 3*R.y**2  # type: ignore
s3 = R.x**2 + R.y**2 + R.z**2  # type: ignore
v1 = R.x*R.i + R.z*R.z*R.j  # type: ignore
v2 = R.x*R.i + R.y*R.j + R.z*R.k  # type: ignore
v3 = R.x**2*R.i + R.y**2*R.j + R.z**2*R.k  # type: ignore


def test_Gradient():
    assert Gradient(s1) == Gradient(R.x*R.y*R.z)
    assert Gradient(s2) == Gradient(R.x + 3*R.y**2)
    assert Gradient(s1).doit() == R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
    assert Gradient(s2).doit() == R.i + 6*R.y*R.j


def test_Divergence():
    assert Divergence(v1) == Divergence(R.x*R.i + R.z*R.z*R.j)
    assert Divergence(v2) == Divergence(R.x*R.i + R.y*R.j + R.z*R.k)
    assert Divergence(v1).doit() == 1
    assert Divergence(v2).doit() == 3
    # issue 22384
    Rc = CoordSys3D('R', transformation='cylindrical')
    assert Divergence(Rc.i).doit() == 1/Rc.r


def test_Curl():
    assert Curl(v1) == Curl(R.x*R.i + R.z*R.z*R.j)
    assert Curl(v2) == Curl(R.x*R.i + R.y*R.j + R.z*R.k)
    assert Curl(v1).doit() == (-2*R.z)*R.i
    assert Curl(v2).doit() == VectorZero()


def test_Laplacian():
    assert Laplacian(s3) == Laplacian(R.x**2 + R.y**2 + R.z**2)
    assert Laplacian(v3) == Laplacian(R.x**2*R.i + R.y**2*R.j + R.z**2*R.k)
    assert Laplacian(s3).doit() == 6
    assert Laplacian(v3).doit() == 2*R.i + 2*R.j + 2*R.k
    assert srepr(Laplacian(s3)) == \
            'Laplacian(Add(Pow(R.x, Integer(2)), Pow(R.y, Integer(2)), Pow(R.z, Integer(2))))'
    
def test_gradient_curvilinear():
    """Test gradient in different coordinate systems"""
    Rc = CoordSys3D('Rc', transformation='cylindrical')
    r, theta, z = Rc.base_scalars()
    
    f1 = r**2
    grad_f1 = gradient(f1)
    assert grad_f1.dot(Rc.i) == 2*r  
    assert grad_f1.dot(Rc.j) == 0    
    assert grad_f1.dot(Rc.k) == 0    
    
    f2 = r*theta
    grad_f2 = gradient(f2)
    assert grad_f2.dot(Rc.i) == theta  
    assert grad_f2.dot(Rc.j) == r      
    assert grad_f2.dot(Rc.k) == 0      
    
    Rs = CoordSys3D('Rs', transformation='spherical')
    r, phi, theta = Rs.base_scalars()
    
    f3 = r
    grad_f3 = gradient(f3)
    assert grad_f3.dot(Rs.i) == 1    
    assert grad_f3.dot(Rs.j) == 0    
    assert grad_f3.dot(Rs.k) == 0    

def test_directional_derivative_curvilinear():
    """Test directional derivative in curvilinear coordinates"""
    Rc = CoordSys3D('Rc', transformation='cylindrical')
    r, theta, z = Rc.base_scalars()
    Omega = Symbol('Omega')
    
    v_field = Omega*r*Rc.j
    dd = directional_derivative(v_field, v_field)
    
    assert dd.dot(Rc.i) == -Omega**2*r
    assert dd.dot(Rc.j) == 0
    assert dd.dot(Rc.k) == 0
    
    v_radial = r*Rc.i
    dd_radial = directional_derivative(v_radial, v_radial)
    assert dd_radial.dot(Rc.i) == r  
    assert dd_radial.dot(Rc.j) == 0
    assert dd_radial.dot(Rc.k) == 0
    
    v_mixed = r*Rc.i + theta*Rc.j
    dd_mixed = directional_derivative(v_mixed, v_mixed)

def test_gradient_zero():
    """Test gradient of zero and constants"""
    assert gradient(0) == Vector.zero
    assert gradient(S.One) == Vector.zero

def test_directional_derivative_errors():
    """Test error handling in directional derivative"""
    R1 = CoordSys3D('R1')
    R2 = CoordSys3D('R2')
    v1 = R1.i
    v2 = R2.i
    
    with raises(ValueError):
        directional_derivative(v1 + v2, v1)
