from sympy.utilities.pytest import XFAIL

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

def test_add_flatten_nested_add():
    from sympy.core.add import Add
    from sympy.core.mul import Mul
    from sympy.abc import x, y, z, a, b, c, d, e, f

    args = Add(z, e, Add(y, x, z), d, z**2, f).args
    assert args == (2*z, e, y, x, d, z**2, f)

    args = (z + e + Add(y, x, z) + d + z**2 + f).args
    assert args == (2*z, e, y, x, d, z**2, f)

def test_add_flatten_nested_mul():
    from sympy.core.add import Add
    from sympy.core.mul import Mul
    from sympy.abc import x, y, z, a, b, c, d, e, f

    args = Add(z, e, Mul(2, Add(y, x, z)), d, z**2, f).args
    assert args == (3*z, e, 2*y, 2*x, d, z**2, f)

    # in case Add(y, x), whose cache is equal to Add(x, y), is already hashed
    # we clear the cache to make sure that this test is independent of cache.
    # Note: the string `(y + x + z)` is parsing with sequential creation of
    # cached Add(y, x) and then Add(Add(y, x), z) and then `flatten` function
    # calls. But the cached Add(y, x) is just a link Add(x, y), so `flatten`
    # function treats `args` as the old ordered tuple (x, y) and insert it
    # incorrectly.
    # Note: with the manually constructed Add(y, x, z) this was not necessary.
    # because it is flatten already.

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

    # See the comment above
    from sympy.core.cache import clear_cache
    clear_cache()

    args = (2*(y + x + z)).args
    assert args == (2*y, 2*x, 2*z)
