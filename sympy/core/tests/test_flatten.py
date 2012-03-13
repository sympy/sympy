from sympy.utilities.pytest import XFAIL

@XFAIL
def test_as_coeff_add():
    from sympy.abc import x
    r = (7 + 3*x + 4*x**2).as_coeff_add()
    assert r == (7, (3*x, 4*x**2))

def test_flatten_sort():
    from sympy.core.add import Add
    from sympy.abc import x, y, z
    e = Add(x, z, y)
    assert e.args == (x, z, y)
    assert (x + y + z).args == (x, y, z)
    e = Add(x, z)
    assert e.args == (x, z)
    e =  3 + x
    assert e.args == (3, x)

def test_add_flatten_nessted_add():
    from sympy.core.add import Add
    from sympy.core.mul import Mul
    from sympy.abc import x, y, z, a, b, c, d, e, f

    args = Add(z, e, Add(y, x, z), d, z**2, f).args
    assert args == (2*z, e, y, x, d, z**2, f)

    args = (z + e + Add(y, x, z) + d + z**2 + f).args
    assert args == (2*z, e, y, x, d, z**2, f)

def test_add_flatten_nessted_mul():
    from sympy.core.add import Add
    from sympy.core.mul import Mul
    from sympy.abc import x, y, z, a, b, c, d, e, f

    args = Add(z, e, Mul(2, Add(y, x, z)), d, z**2, f).args
    assert args == (3*z, e, 2*y, 2*x, d, z**2, f)

    # as Add(x, y) is hashed, and its hash is rather equal to Add(y, x),
    # (others tests have used it probebly)
    # then we clear it to do not influence on this test
    # Note that for the above assertrion (with manually Add construction)
    # it is not needed.
    from sympy.core.cache import clear_cache
    clear_cache()

    args = (z + e + 2*(y + x + z) + d + z**2 + f).args
    assert args == (3*z, e, 2*y, 2*x, d, z**2, f)

def test_mul_flatten():
    from sympy.abc import x, y, z
    from sympy.core.add import Add
    from sympy.core.mul import Mul

    args = Mul(2, y + x + z).args
    assert args == (2*y, 2*x, 2*z)   # it is related with Add.flatten only

    # as Add(x, y) is hashed, and its hash is rather equal to Add(y, x),
    # (others tests have used it probebly)
    # then we clear it to do not influence on this test
    # Note that for the above assertrion (with manually Add construction)
    # it is not needed.
    from sympy.core.cache import clear_cache
    clear_cache()

    args = (2*(y + x + z)).args
    assert args == (2*y, 2*x, 2*z)
