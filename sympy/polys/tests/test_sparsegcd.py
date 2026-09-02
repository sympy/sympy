from sympy.polys.sparsegcd import _gcd_preprocess_polys, gcd_extract, gcd_prs, gcd_terms


def test_gcd_preprocess_polys():
    R, x, y = ring("x, y", ZZ)

    f = x**2*y + x*y
    g = x**3*y**2 + x*y**2
    polynomials = [f, g]
    result = _gcd_preprocess_polys(polynomials)
    expected = ([x**2*y + x*y, x**3*y**2 + x*y**2], {0, 1})
    assert list(ordered(result[0])), result[1] == expected

    f = x**2 - y**2
    g = x**2 - 2*x*y + y**2
    polynomials = [f, g]
    result = _gcd_preprocess_polys(polynomials)
    expected = ([x**2 - y**2, x**2 - 2*x*y + y**2], {0, 1})
    assert list(ordered(result[0])), result[1] == expected

def test_gcd_terms():
    from sympy import ring
    R, x, y = ring("x, y", ZZ)

    polynomials = [x**2 - y**2, x - y]
    ring = polynomials[0].ring
    domain = ring.domain
    result = gcd_terms(polynomials, ring, domain)
    assert result == 1

def test_gcd_extract():
    R, x, y, z = ring("x, y, z", ZZ)

    f = x - y
    g = x**2 - y**2
    assert gcd_extract([f, g], "heugcd") == x - y

    f = x*y + y*z
    g = x*z - y*z
    assert gcd_extract([f, g], "heugcd") == 1

    f = x**2 - x
    g = x**3 - x**2
    assert gcd_extract([f, g], "heugcd") == x**2 - x

    # Polynomials with rational coefficients
    Q, a, b = ring("a, b", QQ)
    f = a**2 - 2*a*b + b**2
    g = a**3 - b**3
    assert gcd_extract([f, g], "prs") == a - b

    f = x*y*z
    g = x**2*y - x*y**2
    h = x**3 - y**3
    assert gcd_extract([f, g, h], "prs") == 1

    f = x*y - y*z
    g = x*z - x*y
    h = x**2 - y**2
    assert gcd_extract([f, g, h], "prs") == 1

def test_gcd_prs():
    R, x, y, z = ring("x, y, z", ZZ)

    f = 4*x**2 + 8*x + 10
    g = 2*x**2 + 6*x + 10
    result = gcd_prs(f, g)
    assert result == 2

    h = 3*x**2 + 9*x + 15
    result = gcd_prs(f, h)

    i = 2*x**3 + 4*x**2 + 2*x
    result = gcd_prs(f, i)
    assert result == 2

    f = 4*x**2 + 8*x + 10*y
    g = 2*x**2 + 6*x + 10*y
    result = gcd_prs(f, g)
    assert result == 2

    h = 3*x**2 + 9*x + 15*y
    result = gcd_prs(f, h)
    assert result == 1

    i = 2*x**3 + 4*x**2*y + 2*x*y**2
    result = gcd_prs(f, i)
    assert result == 2


def test_PolyElement_gcd_ring():
    R, x, y = ring("x, y", QQ)

    f = 0.1*x + 0.1*y + 0.1*y + 0.3
    g = y + 1.0

    expected_gcd = (1, 1/10*x + 1/5*y + 3/10, y + 1)

    assert  f._gcd_ring(g) == expected_gcd

    f = -0.0001*x**4 + 0.0098*x**3 + 0.0299999999999999*x**2 + \
        0.0199999999999998*x
    g = -1.0e-7*x**6 + 3.0e-5*x**5 - 0.00297*x**4 + 0.094*x**3 + 0.297*x**2 \
        + 0.3*x + 0.1

    expected_gcd = (1, -1/10000*x**4 + 49/5000*x**3 + 3/100*x**2 + 1/50*x,
                -944473296573929/9444732965739290427392*x**6 + 3/100000*x**5 -
                297/100000*x**4 + 47/500*x**3 + 297/1000*x**2 + 3/10*x + 1/10)

    assert  f._gcd_ring(g) == expected_gcd

    R, x, y = ring("x, y", ZZ_I)
    f = (-1 + 0*I)*x*y + (0 + -1*I)*y**2
    g = x**3 + (0 + 1*I)*x**2*y + x*y**2 + (0 + 1*I)*y**3

    expected_gcd = (x + I*y, -y, x**2 + y**2)

    assert f._gcd_ring(g) == expected_gcd

    R, y, _t0, _t1 = ring("y, _t0, _t1", EX)
    f = (2 + 0*I)*y**2*_t0 + (0 + 2*I)*y**2*(EX(pi)) + (-2 + 0*I)*y*_t1 + \
                        (0 + -2*I)*y*(EX(pi)) + (-2 + 0*I)*y*_t0

    g = (4 + 0*I)*y**3 + (-4 + 0*I)*y**2

    expected_gcd = (
        y,
        EX(2)*y*_t0 + EX(2*I*pi)*y - EX(2)*_t0 - EX(2)*_t1 + EX(-2*I*pi),
        EX(4)*y**2 - EX(4)*y
    )

    assert f._gcd_ring(g) == expected_gcd

    R, y = ring("y", EX)
    f = (0 + 2*I)*y*(EX(pi)) + (-1 + 0*I)
    g = (0 + 2*I)*y*(EX(pi)) + (-1 + 0*I)

    expected_gcd = (y + EX(I/(2*pi)), EX(2*I*pi), EX(2*I*pi))

    assert f._gcd_ring(g) == expected_gcd

    f = (0 + 2*I)*y*pi + (-1 + 0*I)
    g = (0 + 2*I)*y*pi + (-1 + 0*I)

    expected_gcd = (y + EX(I/(2*pi)), EX(2*I*pi), EX(2*I*pi))

    assert f._gcd_ring(g) == expected_gcd

    R, x, _C3 = ring("x, _C3", ZZ_I)
    f = (0 + 4*I)*x**5 + (0 + -4*I)*x*_C3
    g = (4 + 0*I)*x**5 + (4 + 0*I)*x*_C3

    expected_gcd = (4*x, I*x**4 - I*_C3, x**4 + _C3)

    assert f._gcd_ring(g) == expected_gcd


def test_PolyElement_primitive_wrt():
    R, x = ring("x", ZZ)

    p = 4*x**3 + 8*x**2 + 12*x
    assert p.primitive_wrt(x) == (4, x**3 + 2*x**2 + 3*x)

    p = x**2 + 3*x
    assert p.primitive_wrt(x) == (1, x**2 + 3*x)
