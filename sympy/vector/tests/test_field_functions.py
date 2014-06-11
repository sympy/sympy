from sympy.vector.vector import i, j, k, Vector
from sympy.vector.scalar import x, y, z
from sympy.vector.deloperator import delop
from sympy.simplify import simplify
from sympy.core.symbol import symbols
from sympy.core import S
from sympy import sin, cos

a, b, c = symbols('a b c')

def test_del_operator():

    #Tests for curl
    assert delop ^ Vector.Zero == Vector.Zero
    assert delop.cross(Vector.Zero) == Vector.Zero
    assert delop ^ i == Vector.Zero
    assert delop.cross(2*y**2*j) == Vector.Zero
    v = x*y*z * (i + j + k)
    assert delop ^ v == \
           (-x*y + x*z)*i + (x*y - y*z)*j + (-x*z + y*z)*k
    assert delop ^ v == delop.cross(v)
    assert delop.cross(2*x**2*j) == 4*x*k

    #Tests for divergence
    assert delop & Vector.Zero == S(0)
    assert delop.dot(Vector.Zero) == S(0)
    assert delop & i == S(0)
    assert delop & x**2*i == 2*x
    assert delop & v == x*y + y*z + z*x
    assert delop & v == delop.dot(v)
    assert delop.dot(1/(x*y*z) * (i + j + k)) == \
           - 1 / (x*y*z**2) - 1 / (x*y**2*z) - 1 / (x**2*y*z)
    v = x*i + y*j + z*k
    assert delop & v == 3
    assert simplify(delop & v) == 3

    #Tests for gradient
    assert delop(0) == Vector.Zero
    assert delop(S(0)) == Vector.Zero
    assert delop(x) == i
    assert delop(x*y*z) == y*z*i + z*x*j + x*y*k
    assert delop(2*x**2) == 4*x*i
    assert delop(a*sin(y) / x) == -a*sin(y)/x**2 * i + a*cos(y)/x * j

    #Tests for directional derivative
    assert (Vector.Zero & delop)(a) == S(0)
    assert ((v & delop)(Vector.Zero)) == Vector.Zero
    assert (v & delop)(S(0)) == S(0)
    assert (i & delop)(x) == 1
    assert (j & delop)(y) == 1
    assert (k & delop)(z) == 1
    assert (i & delop)(x*y*z) == y*z
    assert (v & delop)(x) == x
    assert (v & delop)(x*y*z) == 3*x*y*z
    assert (v & delop)(x + y + z) == x + y + z
    assert (v & delop)(v) == v
    assert (i & delop)(v) == i
    assert (j & delop)(v) == j
    assert (k & delop)(v) == k
    assert (v & delop)(Vector.Zero) == Vector.Zero


def test_product_rules():
    """
    Tests the six product rules defined with respect to the Del
    operator

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Del

    """

    #Define the scalar and vector functions
    f = 2*x*y*z
    g = x*y + y*z + z*x
    u = x**2*i + 4*j - y**2*z*k
    v = 4*i + x*y*z*k

    #First product rule
    lhs = delop(f * g)
    rhs = f * delop(g) + g * delop(f)
    assert simplify(lhs) == simplify(rhs)

    #Second product rule
    lhs = delop(u & v)
    rhs = (u ^ (delop ^ v)) + (v ^ (delop ^ u)) + \
          ((u & delop)(v)) + ((v & delop)(u))
    assert lhs.simplify() == rhs.simplify()

    #Third product rule
    lhs = delop & (f*v)
    rhs = (f * (delop & v)) + (v & (delop(f)))
    assert simplify(lhs) == simplify(rhs)

    #Fourth product rule
    lhs = delop & (u ^ v)
    rhs = (v & (delop ^ u)) - (u & (delop ^ v))
    assert simplify(lhs) == simplify(rhs)

    #Fifth product rule
    lhs = delop ^ (f * v)
    rhs = ((delop(f)) ^ v) + (f * (delop ^ v))
    assert lhs.simplify() == rhs.simplify()

    #Sixth product rule
    lhs = delop ^ (u ^ v)
    rhs = u * (delop & v) - v * (delop & u) + \
          (v & delop)(u) - (u & delop)(v)
    assert lhs.simplify() == rhs.simplify()
