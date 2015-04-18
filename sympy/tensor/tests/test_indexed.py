from sympy.core import symbols, Symbol, Tuple, oo, Wild
from sympy.core.compatibility import iterable, range
from sympy.tensor.indexed import IndexException
from sympy.utilities.pytest import raises

# Import test:
from sympy import IndexedBase, Idx, Indexed, DeltaIndexedBase


def assert_equal(expr1, expr2):
    if expr1 != expr2:
        raise AssertionError('{!s} != {!s}'.format(expr1, expr2))


def assert_not_equal(expr1, expr2):
    if expr1 == expr2:
        raise AssertionError('{!s} == {!s}'.format(expr1, expr2))


def assert_equal_modulo_dummies(ein_sum1, ein_sum2):
    index_structure = ein_sum1.index_structure
    all_inner = set()
    for inner in index_structure['inner_list']:
        all_inner |= set(inner)

    substitutions = [[index, Wild(str(index))] for index in all_inner]
    pattern = ein_sum1.subs(substitutions)
    if not (ein_sum2.match(pattern) or ein_sum2 == pattern):
        raise AssertionError('{!s} != {!s}'.format(ein_sum1, ein_sum2))


def test_Idx_construction():
    i, a, b = symbols('i a b', integer=True)
    assert_not_equal(Idx(i), Idx(i, 1))
    assert_equal(Idx(i, a), Idx(i, (0, a - 1)))
    assert_equal(Idx(i, oo), Idx(i, (0, oo)))

    x = symbols('x')
    raises(TypeError, lambda: Idx(x))
    raises(TypeError, lambda: Idx(0.5))
    raises(TypeError, lambda: Idx(i, x))
    raises(TypeError, lambda: Idx(i, 0.5))
    raises(TypeError, lambda: Idx(i, (x, 5)))
    raises(TypeError, lambda: Idx(i, (2, x)))
    raises(TypeError, lambda: Idx(i, (2, 3.5)))


def test_Idx_properties():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).is_integer


def test_Idx_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).lower is None
    assert Idx(i).upper is None
    assert_equal(Idx(i, a).lower, 0)
    assert_equal(Idx(i, a).upper, a - 1)
    assert_equal(Idx(i, 5).lower, 0)
    assert_equal(Idx(i, 5).upper, 4)
    assert_equal(Idx(i, oo).lower, 0)
    assert_equal(Idx(i, oo).upper, oo)
    assert_equal(Idx(i, (a, b)).lower, a)
    assert_equal(Idx(i, (a, b)).upper, b)
    assert_equal(Idx(i, (1, 5)).lower, 1)
    assert_equal(Idx(i, (1, 5)).upper, 5)
    assert_equal(Idx(i, (-oo, oo)).lower, -oo)
    assert_equal(Idx(i, (-oo, oo)).upper, oo)


def test_Idx_fixed_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(2).lower is None
    assert Idx(2).upper is None
    assert_equal(Idx(2, a).lower, 0)
    assert_equal(Idx(2, a).upper, a - 1)
    assert_equal(Idx(2, 5).lower, 0)
    assert_equal(Idx(2, 5).upper, 4)
    assert_equal(Idx(2, oo).lower, 0)
    assert_equal(Idx(2, oo).upper, oo)
    assert_equal(Idx(2, (a, b)).lower, a)
    assert_equal(Idx(2, (a, b)).upper, b)
    assert_equal(Idx(2, (1, 5)).lower, 1)
    assert_equal(Idx(2, (1, 5)).upper, 5)
    assert_equal(Idx(2, (-oo, oo)).lower, -oo)
    assert_equal(Idx(2, (-oo, oo)).upper, oo)


def test_Idx_func_args():
    i, a, b = symbols('i a b', integer=True)
    ii = Idx(i)
    assert_equal(ii.func(*ii.args), ii)
    ii = Idx(i, a)
    assert_equal(ii.func(*ii.args), ii)
    ii = Idx(i, (a, b))
    assert_equal(ii.func(*ii.args), ii)


def test_Idx_subs():
    i, a, b = symbols('i a b', integer=True)
    assert_equal(Idx(i, a).subs(a, b), Idx(i, b))
    assert_equal(Idx(i, a).subs(i, b), Idx(b, a))

    assert_equal(Idx(i).subs(i, 2), Idx(2))
    assert_equal(Idx(i, a).subs(a, 2), Idx(i, 2))
    assert_equal(Idx(i, (a, b)).subs(i, 2), Idx(2, (a, b)))


def test_IndexedBase_sugar():
    i, j = symbols('i j', integer=True)
    a = symbols('a')
    a_indexed = Indexed(a, i, j)
    a_base = IndexedBase(a)
    assert_equal(a_indexed, a_base[i, j])
    assert_equal(a_indexed, a_base[(i, j)])
    assert_equal(a_indexed, a_base[[i, j]])
    assert_equal(a_indexed, a_base[Tuple(i, j)])
    assert all(a.is_Integer for a in a_base[1, 0].args[1:])


def test_IndexedBase_subs():
    i, j, k = symbols('i j k', integer=True)
    a, b = symbols('a b')
    A = IndexedBase(a)
    B = IndexedBase(b)
    assert_equal(A[i], B[i].subs(b, a))


def test_IndexedBase_shape():
    i, j, m, n = symbols('i j m n', integer=True)
    a = IndexedBase('a', shape=(m, m))
    b = IndexedBase('a', shape=(m, n))
    assert_equal(b.shape, Tuple(m, n))
    assert_not_equal(a[i, j], b[i, j])
    assert_equal(a[i, j], b[i, j].subs(n, m))
    assert_equal(b.func(*b.args), b)
    assert_equal(b[i, j].func(*b[i, j].args), b[i, j])
    raises(IndexException, lambda: b[i])
    raises(IndexException, lambda: b[i, i, j])


def test_Indexed_constructor():
    i, j = symbols('i j', integer=True)
    a = Indexed('a', i, j)
    assert_equal(a, Indexed(Symbol('a'), i, j))
    assert_equal(a, Indexed(IndexedBase('a'), i, j))
    raises(TypeError, lambda: Indexed(a, i, j))
    raises(IndexException, lambda: Indexed('a'))


def test_Indexed_func_args():
    i, j = symbols('i j', integer=True)
    a = symbols('a')
    A = Indexed(a, i, j)
    assert_equal(A, A.func(*A.args))


def test_Indexed_subs():
    i, j, k = symbols('i j k', integer=True)
    a, b = symbols('a b')
    A = IndexedBase(a)
    B = IndexedBase(b)
    assert_equal(A[i, j], B[i, j].subs(b, a))
    assert_equal(A[i, j], A[i, k].subs(k, j))


def test_Indexed_properties():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, j)
    assert_equal(A.rank, 2)
    assert_equal(A.indices, (i, j))
    assert_equal(A.base, IndexedBase('A'))
    assert_equal(A.ranges, [None, None])
    raises(IndexException, lambda: A.shape)

    n, m = symbols('n m', integer=True)
    assert_equal(Indexed('A', Idx(i, m), Idx(j, n)).ranges,
                 [Tuple(0, m - 1), Tuple(0, n - 1)])
    assert_equal(Indexed('A', Idx(i, m), Idx(j, n)).shape, Tuple(m, n))
    raises(IndexException, lambda: Indexed("A", Idx(i, m), Idx(j)).shape)


def test_Indexed_shape_precedence():
    i, j = symbols('i j', integer=True)
    o, p = symbols('o p', integer=True)
    n, m = symbols('n m', integer=True)
    a = IndexedBase('a', shape=(o, p))
    assert_equal(a.shape, Tuple(o, p))
    assert_equal(Indexed(a, Idx(i, m), Idx(j, n)).ranges,
                 [Tuple(0, m - 1), Tuple(0, n - 1)])
    assert_equal(Indexed(a, Idx(i, m), Idx(j, n)).shape, Tuple(o, p))
    assert_equal(Indexed(a, Idx(i, m), Idx(j)).ranges,
                 [Tuple(0, m - 1), Tuple(None, None)])
    assert_equal(Indexed(a, Idx(i, m), Idx(j)).shape, Tuple(o, p))


def test_complex_indices():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, i + j)
    assert_equal(A.rank, 2)
    assert_equal(A.indices, (i, i + j))


def test_not_interable():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, i + j)
    assert not iterable(A)


def test_Indexed_coeff():
    N = Symbol('N', integer=True)
    len_y = N
    i = Idx('i', len_y - 1)
    y = IndexedBase('y', shape=(len_y,))
    a = (1 / y[i + 1] * y[i]).coeff(y[i])
    b = (y[i] / y[i + 1]).coeff(y[i])
    assert_equal(a, b)


def test_equality():
    i, j = symbols('i j', cls=Idx)
    a = symbols('a')
    h1 = IndexedBase('h')
    h2 = IndexedBase('h')
    h1i, h2i = h1[i], h1[i]
    delta1, delta2 = symbols('delta1 delta87', cls=DeltaIndexedBase)

    assert_equal(h1, h2)
    assert_equal(h1i, h2i)
    assert_equal(a * h1i, a * h2i)
    assert_equal(delta1[i, j], delta1[i, j])
    assert_not_equal(delta1[i, j], delta2[i, j])
