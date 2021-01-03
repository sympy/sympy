from sympy import Q, symbols

x, y = symbols('x y')
p = symbols('p', positive=True)
n = symbols('n', negative=True)
nz = symbols('nz', nonzero=True)
z = symbols('z', zero=True)

def test_rearrange_eq():
    assert Q.eq(2*x, 2*y).rearrange() == Q.eq(x,y)
    assert Q.eq(x+1, y+1).rearrange() == Q.eq(x,y)
    assert Q.eq(2**x, 2**y).rearrange() == Q.eq(x,y)
    assert Q.eq(3*2**x+3, 3*2**y+3).rearrange() == Q.eq(x,y)

    assert Q.eq(p*x, p*y).rearrange() == Q.eq(x,y)
    assert Q.eq(n*x, n*y).rearrange() == Q.eq(x,y)
    assert Q.eq(nz*x, nz*y).rearrange() == Q.eq(x,y)
    assert Q.eq(z*x, z*y).rearrange() == Q.eq(z*x,z*y)

def test_rearrange_ineq():
    assert Q.gt(-2*x, -2*y).rearrange() == Q.lt(x,y)
    assert Q.gt(2**(-x), 2**(-y)).rearrnage() == Q.lt(x,y)

    assert Q.gt(p*x, p*y).rearrange() == Q.gt(x,y)
    assert Q.gt(n*x, n*y).rearrange() == Q.lt(x,y)
    assert Q.gt(nz*x, nz*y).rearrange() == Q.gt(nz*x,nz*y)
    assert Q.gt(z*x, z*y).rearrange() == Q.gt(z*x,z*y)

def test_solveeqn_eq():
    assert Q.eq(z, x+y-3).solve(x) == Q.eq(x, -y+z+3)
