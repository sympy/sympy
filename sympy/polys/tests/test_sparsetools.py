def test_PolyElement__coeff_split_syms():
    R, x, y, z = ring("x, y, z", ZZ)

    f = 2*x**4 + 3*y**4 + 10*z**2 + 10*x*z**2
    syms = {z}
    result = f._coeff_split_syms(syms)
    assert result == {(0, 0, 0): {(4, 0, 0): 2, (0, 4, 0): 3, (0, 0, 2): 10,
                    (1, 0, 2): 10}}

    g = 3*x**2 + 2*y**2 + 5*z**3 + 7*x**2*y
    syms = {x, y}
    result = g._coeff_split_syms(syms)
    assert result == {(0, 0, 0): {(2, 0, 0): 3, (0, 2, 0): 2, (0, 0, 3): 5,
                    (2, 1, 0): 7}}

def test_PolyElement_coeff_split():
    R, x, y, z = ring("x, y, z", ZZ)

    f = 2*x**4 + 3*y**4 + 10*z**2 + 10*x*z**2
    syms = {z}
    result = f.coeff_split(syms)
    assert result ==[2*x**4 + 3*y**4, 10*x + 10]

    g = 3*x**2 + 2*y**2 + 5*z**3 + 7*x**2*y
    syms = {x, y}
    result = g.coeff_split(syms)
    assert result == [3, 2, 5*z**3, 7]


def test_PolyElement_main_variable():
    R, x, y, z = ring("x, y, z", ZZ)

    p = x**2 + y**3 + z - 2*x*z**2
    assert p.main_variable() == 0

def test_monomial_extract():
    R, x, y, z = ring("x, y, z", ZZ)

    f = x**2*y + z**3
    g = x**3*y + z**3
    h = x**2*y**2 + z**3
    polynomials = [f, g, h]
    result = monomial_extract(polynomials)

    assert result == (polynomials, 1)

    f = x**2*y
    g = x**2*y + x*y
    h = x*y + y**2
    polynomials = [f, g, h]
    result = monomial_extract(polynomials)
    assert result == ([x**2, x**2 + x, x + y], y)