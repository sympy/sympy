from sympy import symbols
from sympy.tensor.array import Array

from sympy.tensor.array.arrayop import tensorproduct, tensorcontraction


def test_tensorproduct():
    x,y,z,t = symbols('x y z t')
    from sympy.abc import a,b,c,d
    assert tensorproduct() == 1
    assert tensorproduct([x]) == Array([x])
    assert tensorproduct([x], [y]) == Array([[x*y]])
    assert tensorproduct([x], [y], [z]) == Array([[[x*y*z]]])
    assert tensorproduct([x], [y], [z], [t]) == Array([[[[x*y*z*t]]]])

    assert tensorproduct(x) == x
    assert tensorproduct(x, y) == x*y
    assert tensorproduct(x, y, z) == x*y*z
    assert tensorproduct(x, y, z, t) == x*y*z*t

    A = Array([x, y])
    B = Array([1, 2, 3])
    C = Array([a, b, c, d])

    assert tensorproduct(A, B, C) == Array([[[a*x, b*x, c*x, d*x], [2*a*x, 2*b*x, 2*c*x, 2*d*x], [3*a*x, 3*b*x, 3*c*x, 3*d*x]],
                                            [[a*y, b*y, c*y, d*y], [2*a*y, 2*b*y, 2*c*y, 2*d*y], [3*a*y, 3*b*y, 3*c*y, 3*d*y]]])

    assert tensorproduct([x, y], [1, 2, 3]) == tensorproduct(A, B)

    assert tensorproduct(A, 2) == Array([2*x, 2*y])
    assert tensorproduct(A, [2]) == Array([[2*x], [2*y]])
    assert tensorproduct([2], A) == Array([[2*x, 2*y]])
    assert tensorproduct(a, A) == Array([a*x, a*y])
    assert tensorproduct(a, A, B) == Array([[a*x, 2*a*x, 3*a*x], [a*y, 2*a*y, 3*a*y]])
    assert tensorproduct(A, B, a) == Array([[a*x, 2*a*x, 3*a*x], [a*y, 2*a*y, 3*a*y]])
    assert tensorproduct(B, a, A) == Array([[a*x, a*y], [2*a*x, 2*a*y], [3*a*x, 3*a*y]])

def test_tensorcontraction():
    from sympy.abc import a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x
    B = Array(range(18), (2, 3, 3))
    assert tensorcontraction(B, (1, 2)) == Array([12, 39])
    C1 = Array([a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x], (2, 3, 2, 2))

    assert tensorcontraction(C1, (0, 2)) == Array([[a + o, b + p], [e + s, f + t], [i + w, j + x]])
    assert tensorcontraction(C1, (0, 2, 3)) == Array([a + p, e + t, i + x])
    assert tensorcontraction(C1, (2, 3)) == Array([[a + d, e + h, i + l], [m + p, q + t, u + x]])
