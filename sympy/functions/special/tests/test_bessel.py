from sympy import jn, yn, symbols, sin, cos, pi, S, jn_zeros
from sympy.functions.special.bessel import fn
from sympy.utilities.pytest import raises

def test_fn():
    x, z = symbols("x z")
    assert fn(1, z) == 1/z**2
    assert fn(2, z) == -1/z + 3/z**3
    assert fn(3, z) == -6/z**2 + 15/z**4
    assert fn(4, z) == 1/z - 45/z**3 + 105/z**5

    raises(TypeError, "fn(x, z)")
    raises(TypeError, "fn(1.5, z)")
    raises(TypeError, "fn(S(1)/2, z)")

def test_jn():
    z = symbols("z")
    assert jn(0, z) == sin(z)/z
    assert jn(1, z) == sin(z)/z**2 - cos(z)/z
    assert jn(2, z) == (3/z**3-1/z)*sin(z) - (3/z**2) * cos(z)
    assert jn(3, z) == (15/z**4 - 6/z**2)*sin(z) + (1/z - 15/z**3)*cos(z)
    assert jn(4, z) == (1/z + 105/z**5 - 45/z**3)*sin(z) + \
                (-105/z**4 + 10/z**2)*cos(z)
    assert jn(5, z) == (945/z**6 - 420/z**4 + 15/z**2)*sin(z) + \
                (-1/z - 945/z**5 + 105/z**3)*cos(z)
    assert jn(6, z) == (-1/z + 10395/z**7 - 4725/z**5 + 210/z**3)*sin(z) + \
                (-10395/z**6 + 1260/z**4 - 21/z**2)*cos(z)

def test_yn():
    z = symbols("z")
    assert yn(0, z) == -cos(z)/z
    assert yn(1, z) == -cos(z)/z**2-sin(z)/z
    assert yn(2, z) == -((3/z**3-1/z)*cos(z)+(3/z**2)*sin(z))

def test_sympify_yn():
    assert S(15) in yn(3, pi).atoms()
    assert yn(3, pi) == 15/pi**4 - 6/pi**2

def eq(a, b, tol=1e-6):
    for x, y in zip(a, b):
        if not (abs(x-y) < tol):
            return False
    return True

def test_jn_zeros():
    assert eq(jn_zeros(0, 4), [3.141592, 6.283185, 9.424777, 12.566370])
    assert eq(jn_zeros(1, 4), [4.493409, 7.725251, 10.904121, 14.066193])
    assert eq(jn_zeros(2, 4), [5.763459, 9.095011, 12.322940, 15.514603])
    assert eq(jn_zeros(3, 4), [6.987932, 10.417118, 13.698023, 16.923621])
    assert eq(jn_zeros(4, 4), [8.182561, 11.704907, 15.039664, 18.301255])

def test_jn_zeros_mpmath():
    try:
        from mpmath import besseljzero
    except ImportError:
        return
    zeros = lambda n, k: jn_zeros(n, k, method='mpmath')
    assert eq(zeros(0, 4), [3.141592, 6.283185, 9.424777, 12.566370])
    assert eq(zeros(1, 4), [4.493409, 7.725251, 10.904121, 14.066193])
    assert eq(zeros(2, 4), [5.763459, 9.095011, 12.322940, 15.514603])
    assert eq(zeros(3, 4), [6.987932, 10.417118, 13.698023, 16.923621])
    assert eq(zeros(4, 4), [8.182561, 11.704907, 15.039664, 18.301255])

