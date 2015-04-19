from sympy.functions import sin, cos
from sympy.core import symbols, S, Wild, Pow
from sympy.core.compatibility import range
from sympy.simplify import simplify
from sympy.tensor.indexed import IndexedBase, Idx, DeltaIndexedBase
from sympy.tensor.indexed_sums import IndexConformanceException
from sympy.utilities.pytest import raises

# Import test:
from sympy import EinsteinSum, get_indices


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


def test_differentiation():
    i, j, k, l = symbols('i j k l', cls=Idx)
    a = symbols('a')
    h, L = symbols('h L', cls=IndexedBase)
    hi, hj = h[i], h[j]
    delta = DeltaIndexedBase()

    expr = hi
    assert_equal(expr.diff(hj), delta[i, j])
    assert_equal(expr.diff(hi), delta[i, i])

    expr = S(2) * hi
    assert_equal(expr.diff(hj), S(2) * delta[i, j])
    assert_equal(expr.diff(hi), S(2) * delta[i, i])
    assert_equal(expr.diff(a), S.Zero)
    assert_equal(EinsteinSum(expr).diff(hj), EinsteinSum(S(2) * delta[i, j]))
    assert_equal(EinsteinSum(expr).diff(hi), EinsteinSum(S(2) * delta[i, i]))

    expr = a * hj * hj / S(2)
    assert_equal(expr.diff(hi), a * h[j] * delta[i, j])
    assert_equal(expr.diff(a), hj * hj / S(2))
    assert_equal(expr.diff(a, 2), S.Zero)
    assert_equal(EinsteinSum(expr).diff(hi), EinsteinSum(a * h[i]))

    expr = a * sin(hj * hj)
    assert_equal(expr.diff(hi), a * cos(hj * hj) * S(2) * hj * delta[i, j])

    expr = a * EinsteinSum(L[i, j] * h[j])
    assert_equal(expr.diff(L[k, l]), a * EinsteinSum(delta[i, k] * h[l]))


def test_simplify_contractions():
    i, j, k = symbols('i j k', cls=Idx)
    a = symbols('a')
    h = IndexedBase('h')
    hi, hj, hk = h[i], h[j], h[k]
    delta = DeltaIndexedBase()

    def run_test(expr1, expr2):
        assert_equal_modulo_dummies(expr1.simplify_deltas(), expr2)

    expr1 = EinsteinSum(a * hi * delta[i, j] * hj)
    expr2 = EinsteinSum(a * hi * hi)
    run_test(expr1, expr2)

    expr1 = EinsteinSum(a * hi * delta[i, k] * delta[k, j] * hj)
    expr2 = EinsteinSum(a * hj * hj)
    run_test(expr1, expr2)

    expr1 = EinsteinSum(S(2) * hj * delta[i, j] * hk)
    expr2 = EinsteinSum(S(2) * hi * hk)
    run_test(expr1, expr2)

    expr1 = EinsteinSum(delta[i, j] * delta[i, j])
    expr2 = EinsteinSum(delta[i, i])
    run_test(expr1, expr2)

    expr1 = EinsteinSum(delta[i, j] * delta[j, k] * delta[k, i])
    expr2 = EinsteinSum(delta[i, i])
    run_test(expr1, expr2)

    expr1 = EinsteinSum(S(2) * delta[i, i])
    expr2 = EinsteinSum(S(2) * delta[i, i])
    run_test(expr1, expr2)


def test_canonicalize_inner():
    i, j, k, l, m, n = symbols('i j k l m n')
    h, s, L, Q, phi = symbols('h s L Q phi', cls=IndexedBase)

    expr = EinsteinSum(s[i] + L[i, j] * h[j] + L[i, k] * h[k])
    expected_result = EinsteinSum(s[i] + 2 * L[i, j] * h[j])
    assert_equal(expr.canonicalize_inner(), expected_result)

    expr = EinsteinSum(L[i, j] * h[i] * h[j] + L[j, i] * h[i] * h[j])
    expected_result = EinsteinSum(2 * L[i, j] * h[i] * h[j])
    assert_equal(expr.canonicalize_inner(), expected_result)

    expr = EinsteinSum(s[i] + L[i, j] * h[j] + Q[i, l, m] * h[l] * h[m]
                       + Q[i, m, n] * h[m] * h[n])
    expected_result = EinsteinSum(s[i] + L[i, j] * h[j] +
                                  2 * Q[i, j, l] * h[j] * h[l])
    assert_equal(expr.canonicalize_inner(), expected_result)


def test_simplify():
    i, j = symbols('i j', cls=Idx)
    a = symbols('a')
    h = IndexedBase('h')
    hi, hj = h[i], h[j]
    delta = DeltaIndexedBase()

    expr_list = []
    expr_list.append(EinsteinSum(a * hi * delta[i, j] * hj))
    expr_list.append(EinsteinSum(a * hi * delta[i, j] * hj + S(7)))

    expected_list = []
    expected_list.append(EinsteinSum(a * hi * hi))
    expected_list.append(EinsteinSum(a * hi * hi + S(7)))

    for i in range(len(expr_list)):
        assert_equal_modulo_dummies(simplify(expr_list[i]), expected_list[i])


def test_index_structure():
    i, j, k = symbols('i j k', cls=Idx)
    a = symbols('a')
    h, L = IndexedBase('h'), IndexedBase('L')
    hi, hj, hk, Ljk = h[i], h[j], h[k], L[j, k]
    delta = DeltaIndexedBase()

    def run_test(expr, expected_structure):
        structure = EinsteinSum(expr).index_structure
        assert_equal(expected_structure, structure)

    expr = a * hi * delta[i, j] * hj
    expected_structure = {
        'monomial_list': [expr], 'outer': set(), 'inner_list': [set([i, j])]}
    run_test(expr, expected_structure)

    expr = a * hi * hi
    expected_structure = {
        'monomial_list': [expr], 'outer': set(), 'inner_list': [set([i])]}
    run_test(expr, expected_structure)

    expr = a * hi * delta[i, k] * delta[k, j] * hj
    expected_structure = {
        'monomial_list': [expr], 'outer': set(), 'inner_list': [set([i, j, k])]}
    run_test(expr, expected_structure)

    expr = S(2) * hj * delta[i, j] * hk
    expected_structure = {
        'monomial_list': [expr], 'outer': set([i, k]), 'inner_list': [set([j])]}
    run_test(expr, expected_structure)

    expr = delta[i, j] * delta[i, j]
    expected_structure = {
        'monomial_list': [expr], 'outer': set(), 'inner_list': [set([i, j])]}
    run_test(expr, expected_structure)

    expr = delta[i, j] * delta[j, k] * delta[k, i]
    expected_structure = {
        'monomial_list': [expr], 'outer': set(), 'inner_list': [set([i, j, k])]}
    run_test(expr, expected_structure)

    expr = Ljk * (delta[i, j] * hi + a * hj)
    expected_structure = {
        'monomial_list': [a * Ljk * hj, Ljk * delta[i, j] * hi],
        'outer': set([k]),
        'inner_list': [set([j]), set([i, j])]}
    run_test(expr, expected_structure)

    expr = a * hi * delta[i, j] * hj + a * hi * delta[i, k] * delta[k, j] * hj
    expected_structure = {
        'monomial_list': [a * hi * delta[i, j] * hj,
                          a * hi * delta[i, k] * delta[k, j] * hj],
        'outer': set(), 'inner_list': [set([i, j]), set([i, j, k])]}
    run_test(expr, expected_structure)

    expr = (a * hi * delta[i, k] * delta[k, j] * hj
            + delta[i, j] * delta[j, k] * delta[k, i])
    expected_structure = {
        'monomial_list': [delta[i, j] * delta[j, k] * delta[k, i],
                          a * hi * delta[i, k] * delta[k, j] * hj],
        'outer': set(), 'inner_list': [set([i, j, k]), set([i, j, k])]}
    run_test(expr, expected_structure)


def test_index_conformance():
    i, j, k = symbols('i j k', cls=Idx)
    h = IndexedBase('h')
    hi, hj, hk = h[i], h[j], h[k]
    delta = DeltaIndexedBase()

    expr_list = []
    expr_list.append(hi + hj)
    expr_list.append(S(2) * hj * delta[i, j] * hk + delta[i, j] * delta[i, j])

    for expr in expr_list:
        raises(IndexConformanceException, lambda: EinsteinSum(expr))


def test_numpify_norm():
    """Ensure numpification can compute vector norms."""
    try:
        import numpy as np

        i = symbols('i', cls=Idx)
        h = symbols('h', cls=IndexedBase)

        m = 3
        h_arr = np.random.rand(m)

        expr = EinsteinSum(S(-3) * h[i] * h[i])
        func = expr.numpify([h], [])
        val = func(h_arr)

        expected_val = 0.
        for ii in range(m):
            expected_val += -3. * h_arr[ii] * h_arr[ii]
        assert np.allclose(val, expected_val)
    except ImportError:
        pass


def test_numpify_matrix_mult():
    """Ensure numpification can compute matrix multiplication."""
    try:
        import numpy as np

        i, j, k = symbols('i j k', cls=Idx)
        A, B = symbols('A, B', cls=IndexedBase)

        m = 7
        A_arr = np.random.rand(m, m)
        B_arr = np.random.rand(m, m)

        expr = EinsteinSum(A[i, j] * B[j, k])
        func = expr.numpify([A, B], [i, k])
        mult = func(A_arr, B_arr)

        expected_mult = A_arr.dot(B_arr)
        assert np.allclose(mult, expected_mult)
    except ImportError:
        pass


def test_numpify_vec_1():
    try:
        import numpy as np

        i, j, k = symbols('i j k', cls=Idx)
        s, h, L, Q = symbols('s, h, L, Q', cls=IndexedBase)

        m, n = 3, 4
        arr_func = np.random.rand
        s_arr = arr_func(m)
        h_arr = arr_func(n)
        L_arr = arr_func(m * n).reshape(m, n)

        expr = EinsteinSum(2 * s[i] - 0.5 * L[i, j] * h[j])
        func = expr.numpify([L, h, s], [i])
        vec = func(L_arr, h_arr, s_arr)

        expected_vec = np.zeros(m)
        for ii in range(m):
            for jj in range(n):
                expected_vec[ii] += -0.5 * L_arr[ii, jj] * h_arr[jj]
            expected_vec[ii] += 2 * s_arr[ii]
        assert np.allclose(vec, expected_vec)
    except ImportError:
        pass


def test_numpify_vec_2():
    try:
        import numpy as np

        i, j, k = symbols('i j k', cls=Idx)
        s, h, L, Q = symbols('s, h, L, Q', cls=IndexedBase)

        m, n = 3, 4
        arr_func = np.random.rand
        s_arr = arr_func(m)
        h_arr = arr_func(n)
        L_arr = arr_func(m * n).reshape(m, n)
        Q_arr = arr_func(m * n * n).reshape(m, n, n)

        expr = EinsteinSum(S(-3) * s[i] + S(7) * L[i, j] * h[j]
                           + Q[i, j, k] * h[j] * h[k])
        func = expr.numpify([s, L, h, Q], [i])
        raises(IndexConformanceException,
               lambda: expr.numpify([s, L, h, Q], [i, j]))
        vec = func(s_arr, L_arr, h_arr, Q_arr)

        expected_vec = np.zeros(m)
        for ii in range(m):
            for jj in range(n):
                for kk in range(n):
                    expected_vec[ii] += (Q_arr[ii, jj, kk] * h_arr[jj]
                                         * h_arr[kk])
                expected_vec[ii] += 7. * L_arr[ii, jj] * h_arr[jj]
            expected_vec[ii] += -3. * s_arr[ii]
        assert np.allclose(vec, expected_vec)
    except ImportError:
        pass


def test_numpify_internal_contraction():
    try:
        import numpy as np

        i, j, k = symbols('i j k')
        A = symbols('A', cls=IndexedBase)
        delta = DeltaIndexedBase()

        A_array = np.array([[[1, 2], [-1, 7]], [[-3, 2], [-1, 2]]])
        delta_array = np.eye(2)
        expected_result = np.array([8., -1.])

        expr = EinsteinSum(A[i, j, j])
        func = expr.numpify([A], [i])
        assert np.allclose(func(A_array), expected_result)

        expr = EinsteinSum(A[i, j, k] * delta[j, k])
        func = expr.numpify([A, delta], [i])
        assert np.allclose(func(A_array, delta_array), expected_result)
    except ImportError:
        pass


def test_numpify_delta():
    try:
        import numpy as np

        i, j, k, l = symbols('i j k l', cls=Idx)
        A = symbols('A', cls=IndexedBase)
        delta1, delta2 = symbols('delta1 delta87', cls=DeltaIndexedBase)

        m, n = 4, 3
        A_array = np.random.rand(m, n)
        delta1_array = np.eye(m)
        delta2_array = np.eye(n)

        expr = EinsteinSum(delta1[i, j] * A[j, k] * delta2[k, l])
        func = expr.numpify([delta2, A, delta1], [i, l])
        assert np.allclose(func(delta2_array, A_array, delta1_array), A_array)
    except ImportError:
        pass


def test_numpify_complex_trace():
    try:
        import numpy as np

        i = symbols('i', cls=Idx)
        A = symbols('A', cls=IndexedBase)

        m = 4
        A_array = np.random.rand(m, m) + np.random.rand(m, m) * 1j

        expr = EinsteinSum(A[i, i])
        func = expr.numpify([A], [], dtype=complex)
        assert np.allclose(func(A_array), A_array.trace())
    except ImportError:
        pass


def test_numpify_exceptions():
    try:
        import numpy as np

        i, j, z = symbols('i j z', cls=Idx)
        s = symbols('s', cls=IndexedBase)

        # No exception.
        expr = EinsteinSum(S(-7.8) * s[i, j, j])
        func = expr.numpify([s], [i])
        func(np.random.rand(2, 3, 3))

        # Exceptions:
        expr = EinsteinSum(S(-7.8) * s[i, j, j])
        raises(IndexConformanceException, lambda: expr.numpify([s], [i, z]))
        raises(TypeError, lambda: expr.numpify([z], [i]))
        raises(ValueError, lambda: expr.numpify([], [i]))

        expr = EinsteinSum(z * s[i, j, j])
        raises(ValueError, lambda: expr.numpify([s], [i]))

        expr = EinsteinSum(s[1])
        raises(TypeError, lambda: expr.numpify([s], [1]))
    except ImportError:
        pass


def test_get_indices_trivial():
    x, y = symbols('x y')
    for expr in [x, x * y, x + y, x ** y]:
        assert_equal(get_indices(expr), set())


def test_get_indices_mul():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert_equal(get_indices(x[j] * y[i]), set([i, j]))
    assert_equal(get_indices(x[i] * y[j]), set([i, j]))


def test_get_indices_add():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j = Idx('i'), Idx('j')
    assert_equal(get_indices(x[i] + 2 * y[i]), set([i]))
    assert_equal(get_indices(y[i] + 2 * A[i, j] * x[j]), set([i, j]))
    assert_equal(get_indices(y[i] + 2 * x[j]), set([i, j]))


def test_get_indices_Pow():
    x = IndexedBase('x')
    i = Idx('i')
    assert_equal(get_indices(Pow(x[i], 2)), set([i]))
