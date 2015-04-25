# -*- coding: utf-8 -*-

from sympy.matrices import Matrix
from sympy.tensor.arraypy import Arraypy, Tensor
from sympy.tensor.riemannian_geometry import scal_prod, christoffel_1,\
    christoffel_2, covar_der, covar_der_XY, riemann, ricci, scal_curv, k_sigma
from sympy import symbols, cos, sin


def test_scal_prod_gxy_list():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    X = [1, 2]
    Y = [3, 4]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])

    ten_g = Arraypy((2, 2)).to_tensor((-1, -1))
    ten_g[0, 0] = cos(x2)**2
    ten_g[0, 1] = 0
    ten_g[1, 0] = 0
    ten_g[1, 1] = 1

    ten_g1 = Arraypy([2, 2, 1]).to_tensor((-1, -1))
    ten_g1[1, 1] = cos(x2)**2
    ten_g1[1, 2] = 0
    ten_g1[2, 1] = 0
    ten_g1[2, 2] = 1

    assert scal_prod(X, Y, g) == 3 * cos(x2)**2 + 8
    assert scal_prod(X, Y, ten_g) == 3 * cos(x2)**2 + 8
    assert scal_prod(X, Y, ten_g1) == 3 * cos(x2)**2 + 8


def test_christoffel_1_gtnsr():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Arraypy([2, 2, 1])
    g = g.to_tensor((-1, -1))
    g[1, 1] = cos(x2)**2
    g[1, 2] = 0
    g[2, 1] = 0
    g[2, 2] = 1

    res_arr = Arraypy([3, 2, 1])
    res_arr[1, 1, 1] = 0
    res_arr[1, 1, 2] = sin(x2) * cos(x2)
    res_arr[1, 2, 1] = -sin(x2) * cos(x2)
    res_arr[2, 1, 1] = -sin(x2) * cos(x2)
    res_arr[1, 2, 2] = 0
    res_arr[2, 2, 1] = 0
    res_arr[2, 1, 2] = 0
    res_arr[2, 2, 2] = 0
    res_ten = res_arr.to_tensor((-1, -1, -1))

    print('test_christoffel_1_gtnsr_t  <=== actual test code')
    assert christoffel_1(g, var_list) == res_ten
    assert isinstance(christoffel_1(g, var_list), Tensor)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    assert christoffel_1(g, var_list, 't') == res_ten
    assert isinstance(christoffel_1(g, var_list, 't'), Tensor)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    print('test_christoffel_1_gtnsr_a  <=== actual test code')
    assert christoffel_1(g, var_list, 'a') == res_arr
    assert isinstance(christoffel_1(g, var_list, 'a'), Arraypy)


def test_christoffel_1_gm():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])

    res_arr = Arraypy([3, 2, 0])
    res_arr[0, 0, 0] = 0
    res_arr[0, 0, 1] = sin(x2) * cos(x2)
    res_arr[0, 1, 0] = -sin(x2) * cos(x2)
    res_arr[1, 0, 0] = -sin(x2) * cos(x2)
    res_arr[0, 1, 1] = 0
    res_arr[1, 1, 0] = 0
    res_arr[1, 0, 1] = 0
    res_arr[1, 1, 1] = 0
    res_ten = res_arr.to_tensor((-1, -1, -1))

    print('test_christoffel_1_gm_t  <=== actual test code')
    assert christoffel_1(g, var_list) == res_ten
    assert isinstance(christoffel_1(g, var_list), Tensor)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    assert christoffel_1(g, var_list, 't') == res_ten
    assert isinstance(christoffel_1(g, var_list, 't'), Tensor)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    print('test_christoffel_1_gm_a  <=== actual test code')
    assert christoffel_1(g, var_list, 'a') == res_arr
    assert isinstance(christoffel_1(g, var_list, 'a'), Arraypy)


def test_christoffel_2_gtnsr():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Arraypy([2, 2, 1])
    g = g.to_tensor((-1, -1))
    g[1, 1] = cos(x2)**2
    g[1, 2] = 0
    g[2, 1] = 0
    g[2, 2] = 1

    res_arr = Arraypy([3, 2, 1])
    res_arr[1, 1, 1] = 0
    res_arr[1, 1, 2] = sin(x2) * cos(x2)
    res_arr[1, 2, 1] = -sin(x2) / cos(x2)
    res_arr[2, 1, 1] = -sin(x2) / cos(x2)
    res_arr[1, 2, 2] = 0
    res_arr[2, 2, 1] = 0
    res_arr[2, 1, 2] = 0
    res_arr[2, 2, 2] = 0
    res_ten = res_arr.to_tensor((-1, -1, -1))

    print('test_christoffel_2_gtnsr_t  <=== actual test code')
    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list), Tensor)
    assert christoffel_2(g, var_list).type_pq == (0, 3)

    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list, 't'), Tensor)
    assert christoffel_2(g, var_list).type_pq == (0, 3)

    print('test_christoffel_2_gtnsr_a  <=== actual test code')
    assert christoffel_2(g, var_list, 'a') == res_arr
    assert isinstance(christoffel_2(g, var_list, 'a'), Arraypy)


def test_christoffel_2_gm():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])

    res_arr = Arraypy([3, 2, 0])
    res_arr[0, 0, 0] = 0
    res_arr[0, 0, 1] = sin(x2) * cos(x2)
    res_arr[0, 1, 0] = -sin(x2) / cos(x2)
    res_arr[1, 0, 0] = -sin(x2) / cos(x2)
    res_arr[0, 1, 1] = 0
    res_arr[1, 1, 0] = 0
    res_arr[1, 0, 1] = 0
    res_arr[1, 1, 1] = 0
    res_ten = res_arr.to_tensor((-1, -1, -1))

    print('test_christoffel_2_gm_t  <=== actual test code')
    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list), Tensor)
    assert christoffel_2(g, var_list).type_pq == (0, 3)

    assert christoffel_2(g, var_list, 't') == res_ten
    assert isinstance(christoffel_2(g, var_list, 't'), Tensor)
    assert christoffel_2(g, var_list).type_pq == (0, 3)

    print('test_christoffel_2_gm_a  <=== actual test code')
    assert christoffel_2(g, var_list, 'a') == res_arr
    assert isinstance(christoffel_2(g, var_list, 'a'), Arraypy)


def test_covar_der():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])
    X = [x1 * x2**3, x1 - cos(x2)]

    res_arr = Arraypy([2, 2, 0])
    res_arr[0, 0] = x2**3 - (x1 - cos(x2)) * sin(x2) / cos(x2)
    res_arr[0, 1] = x1 * x2**3 * sin(x2) * cos(x2) + 1
    res_arr[1, 0] = -x1 * x2**3 * sin(x2) / cos(x2) + 3 * x1 * x2**2
    res_arr[1, 1] = sin(x2)
    res_ten = res_arr.to_tensor((-1, 1))

    print('test_covar_der_t  <=== actual test code')
    assert covar_der(X, g, var_list) == res_ten
    assert isinstance(covar_der(X, g, var_list), Tensor)
    assert covar_der(X, g, var_list).type_pq == (1, 1)

    assert covar_der(X, g, var_list, 't') == res_ten
    assert isinstance(covar_der(X, g, var_list, 't'), Tensor)
    assert covar_der(X, g, var_list, 't').type_pq == (1, 1)

    print('test_covar_der_a  <=== actual test code')
    assert covar_der(X, g, var_list, 'a') == res_arr
    assert isinstance(covar_der(X, g, var_list, 'a'), Arraypy)


def test_covar_der_XY():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])
    X = [x1 * x2**3, x1 - cos(x2)]
    Y = [1, 2]

    res_arr = Arraypy([1, 2, 0])
    res_arr[0] = -2 * x1 * x2**3 * \
        sin(x2) / cos(x2) + 6 * x1 * x2**2 + \
        x2**3 - (x1 - cos(x2)) * sin(x2) / cos(x2)
    res_arr[1] = x1 * x2**3 * sin(x2) * cos(x2) + 2 * sin(x2) + 1
    res_ten = res_arr.to_tensor(1)

    print('test_covar_der_XY_t  <=== actual test code')
    assert covar_der_XY(X, Y, g, var_list) == res_ten
    assert isinstance(covar_der_XY(X, Y, g, var_list), Tensor)
    assert covar_der_XY(X, Y, g, var_list).type_pq == (1, 0)

    assert covar_der_XY(X, Y, g, var_list, 't') == res_ten
    assert isinstance(covar_der_XY(X, Y, g, var_list, 't'), Tensor)
    assert covar_der_XY(X, Y, g, var_list).type_pq == (1, 0)

    print('test_covar_der_XY_a  <=== actual test code')
    assert covar_der_XY(X, Y, g, var_list, 'a') == res_arr
    assert isinstance(covar_der_XY(X, Y, g, var_list, 'a'), Arraypy)


def test_riemann_gtnsr():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Arraypy([2, 2, 1])
    g = g.to_tensor((-1, -1))
    g[1, 1] = cos(x2)**2
    g[1, 2] = 0
    g[2, 1] = 0
    g[2, 2] = 1

    res_arr = Arraypy([4, 2, 1])
    res_arr[1, 1, 1, 1] = 0
    res_arr[1, 1, 1, 2] = 0
    res_arr[1, 1, 2, 1] = 0
    res_arr[1, 1, 2, 2] = 0
    res_arr[1, 2, 1, 1] = 0
    res_arr[1, 2, 2, 1] = 1
    res_arr[1, 2, 2, 2] = 0
    res_arr[1, 2, 1, 2] = -cos(x2)**2
    res_arr[2, 1, 1, 1] = 0
    res_arr[2, 1, 1, 2] = cos(x2)**2
    res_arr[2, 2, 1, 1] = 0
    res_arr[2, 2, 2, 1] = 0
    res_arr[2, 1, 2, 2] = 0
    res_arr[2, 1, 2, 1] = -1
    res_arr[2, 2, 1, 2] = 0
    res_arr[2, 2, 2, 2] = 0
    res_ten = res_arr.to_tensor((-1, -1, -1, 1))

    print('test_riemann_gtnsr_t  <=== actual test code')
    assert riemann(g, var_list) == res_ten
    assert isinstance(riemann(g, var_list), Tensor)
    assert riemann(g, var_list).type_pq == (1, 3)

    assert riemann(g, var_list, 't') == res_ten
    assert isinstance(riemann(g, var_list, 't'), Tensor)
    assert riemann(g, var_list).type_pq == (1, 3)

    print('test_riemann_gtnsr_a  <=== actual test code')
    assert riemann(g, var_list, 'a') == res_arr
    assert isinstance(riemann(g, var_list, 'a'), Arraypy)


def test_riemann_gm():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])

    res_arr = Arraypy([4, 2, 0])
    res_arr[0, 0, 0, 0] = 0
    res_arr[0, 0, 0, 1] = 0
    res_arr[0, 0, 1, 1] = 0
    res_arr[0, 1, 1, 1] = 0
    res_arr[1, 1, 1, 1] = 0
    res_arr[1, 1, 1, 0] = 0
    res_arr[1, 1, 0, 0] = 0
    res_arr[1, 0, 0, 0] = 0
    res_arr[1, 0, 1, 0] = -1
    res_arr[0, 1, 0, 1] = -cos(x2)**2
    res_arr[1, 1, 0, 1] = 0
    res_arr[0, 0, 1, 0] = 0
    res_arr[1, 0, 1, 1] = 0
    res_arr[0, 1, 0, 0] = 0
    res_arr[0, 1, 1, 0] = 1
    res_arr[1, 0, 0, 1] = cos(x2)**2
    res_ten = res_arr.to_tensor((-1, -1, -1, 1))

    print('test_riemann_gm_t  <=== actual test code')
    assert riemann(g, var_list) == res_ten
    assert isinstance(riemann(g, var_list), Tensor)
    assert riemann(g, var_list).type_pq == (1, 3)

    assert riemann(g, var_list, 't') == res_ten
    assert isinstance(riemann(g, var_list, 't'), Tensor)
    assert riemann(g, var_list).type_pq == (1, 3)

    print('test_riemann_gm_a  <=== actual test code')
    assert riemann(g, var_list, 'a') == res_arr
    assert isinstance(riemann(g, var_list, 'a'), Arraypy)


def test_ricci_riemtnsr1():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    riemann_arr = Arraypy([4, 2, 1])
    riemann_arr[1, 1, 1, 1] = 0
    riemann_arr[1, 1, 1, 2] = 0
    riemann_arr[1, 1, 2, 1] = 0
    riemann_arr[1, 1, 2, 2] = 0
    riemann_arr[1, 2, 1, 1] = 0
    riemann_arr[1, 2, 2, 1] = 1
    riemann_arr[1, 2, 2, 2] = 0
    riemann_arr[1, 2, 1, 2] = -cos(x2)**2
    riemann_arr[2, 1, 1, 1] = 0
    riemann_arr[2, 1, 1, 2] = cos(x2)**2
    riemann_arr[2, 2, 1, 1] = 0
    riemann_arr[2, 2, 2, 1] = 0
    riemann_arr[2, 1, 2, 2] = 0
    riemann_arr[2, 1, 2, 1] = -1
    riemann_arr[2, 2, 1, 2] = 0
    riemann_arr[2, 2, 2, 2] = 0
    riemann_ten = riemann_arr.to_tensor((-1, -1, -1, 1))

    res_arr = Arraypy([2, 2, 1])
    res_arr[1, 1] = cos(x2)**2
    res_arr[1, 2] = 0
    res_arr[2, 2] = 0
    res_arr[2, 2] = 1
    res_ten = res_arr.to_tensor((-1, -1))

    print('test_ricci_riemtnsr1_t  <=== actual test code')
    assert ricci(riemann_ten, var_list) == res_ten
    assert isinstance(ricci(riemann_ten, var_list), Tensor)
    assert ricci(riemann_ten, var_list).type_pq == (0, 2)

    assert ricci(riemann_ten, var_list, 't') == res_ten
    assert isinstance(ricci(riemann_ten, var_list, 't'), Tensor)
    assert ricci(riemann_ten, var_list, 't').type_pq == (0, 2)

    print('test_ricci_riemtnsr1_a  <=== actual test code')
    assert ricci(riemann_ten, var_list, 'a') == res_arr
    assert isinstance(ricci(riemann_ten, var_list, 'a'), Arraypy)


def test_ricci_riemtnsr0():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    riemann_arr = Arraypy([4, 2, 0])
    riemann_arr[0, 0, 0, 0] = 0
    riemann_arr[0, 0, 0, 1] = 0
    riemann_arr[0, 0, 1, 1] = 0
    riemann_arr[0, 1, 1, 1] = 0
    riemann_arr[1, 1, 1, 1] = 0
    riemann_arr[1, 1, 1, 0] = 0
    riemann_arr[1, 1, 0, 0] = 0
    riemann_arr[1, 0, 0, 0] = 0
    riemann_arr[1, 0, 1, 0] = -1
    riemann_arr[0, 1, 0, 1] = -cos(x2)**2
    riemann_arr[1, 1, 0, 1] = 0
    riemann_arr[0, 0, 1, 0] = 0
    riemann_arr[1, 0, 1, 1] = 0
    riemann_arr[0, 1, 0, 0] = 0
    riemann_arr[0, 1, 1, 0] = 1
    riemann_arr[1, 0, 0, 1] = cos(x2)**2
    riemann_ten = riemann_arr.to_tensor((-1, -1, -1, 1))

    res_arr = Arraypy([2, 2, 0])
    res_arr[0, 0] = cos(x2)**2
    res_arr[0, 1] = 0
    res_arr[1, 1] = 0
    res_arr[1, 1] = 1
    res_ten = res_arr.to_tensor((-1, -1))

    print('test_ricci_riemtnsr0_t  <=== actual test code')
    assert ricci(riemann_ten, var_list) == res_ten
    assert isinstance(ricci(riemann_ten, var_list), Tensor)
    assert ricci(riemann_ten, var_list).type_pq == (0, 2)

    assert ricci(riemann_ten, var_list, 't') == res_ten
    assert isinstance(ricci(riemann_ten, var_list, 't'), Tensor)
    assert ricci(riemann_ten, var_list, 't').type_pq == (0, 2)

    print('test_ricci_riemtnsr0_a  <=== actual test code')
    assert ricci(riemann_ten, var_list, 'a') == res_arr
    assert isinstance(ricci(riemann_ten, var_list, 'a'), Arraypy)


def test_scal_curv():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    g = Matrix([[cos(x2)**2, 0], [0, 1]])

    g_ten = Arraypy([2, 2, 1]).to_tensor((-1, -1))
    g_ten[1, 1] = cos(x2)**2
    g_ten[1, 2] = 0
    g_ten[2, 1] = 0
    g_ten[2, 2] = 1

    g_ten0 = Arraypy([2, 2, 0]).to_tensor((-1, -1))
    g_ten0[0, 0] = cos(x2)**2
    g_ten0[0, 1] = 0
    g_ten0[1, 0] = 0
    g_ten0[1, 1] = 1

    ricci = Matrix([[cos(x2)**2, 0], [0, 1]])
    ricci_ten = Arraypy([2, 2, 1]).to_tensor((-1, -1))
    ricci_ten[1, 1] = cos(x2)**2
    ricci_ten[1, 2] = 0
    ricci_ten[2, 1] = 0
    ricci_ten[2, 2] = 1

    ricci_ten0 = Arraypy([2, 2, 0]).to_tensor((-1, -1))
    ricci_ten0[0, 0] = cos(x2)**2
    ricci_ten0[0, 1] = 0
    ricci_ten0[1, 0] = 0
    ricci_ten0[1, 1] = 1

    assert scal_curv(g, ricci, var_list) == 1
    assert scal_curv(g_ten, ricci, var_list) == 1
    assert scal_curv(g, ricci_ten, var_list) == 1
    assert scal_curv(g_ten, ricci_ten, var_list) == 1
    assert scal_curv(g_ten0, ricci, var_list) == 1
    assert scal_curv(g_ten0, ricci_ten, var_list) == 1
    assert scal_curv(g, ricci_ten0, var_list) == 1
    assert scal_curv(g_ten, ricci_ten0, var_list) == 1
    assert scal_curv(g_ten0, ricci_ten0, var_list) == 1


def test_k_sigma():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Matrix([[cos(x2)**2, 0], [0, 1]])
    X = [1, 2]
    Y = [3, 4]

    g_ten = Arraypy([2, 2, 1]).to_tensor((-1, -1))
    g_ten[1, 1] = cos(x2)**2
    g_ten[1, 2] = 0
    g_ten[2, 1] = 0
    g_ten[2, 2] = 1

    g_ten0 = Arraypy([2, 2, 0]).to_tensor((-1, -1))
    g_ten0[0, 0] = cos(x2)**2
    g_ten0[0, 1] = 0
    g_ten0[1, 0] = 0
    g_ten0[1, 1] = 1

    riemann_arr = Arraypy([4, 2, 0])
    riemann_arr[0, 0, 0, 0] = 0
    riemann_arr[0, 0, 0, 1] = 0
    riemann_arr[0, 0, 1, 1] = 0
    riemann_arr[0, 1, 1, 1] = 0
    riemann_arr[1, 1, 1, 1] = 0
    riemann_arr[1, 1, 1, 0] = 0
    riemann_arr[1, 1, 0, 0] = 0
    riemann_arr[1, 0, 0, 0] = 0
    riemann_arr[1, 0, 1, 0] = -1
    riemann_arr[0, 1, 0, 1] = -cos(x2)**2
    riemann_arr[1, 1, 0, 1] = 0
    riemann_arr[0, 0, 1, 0] = 0
    riemann_arr[1, 0, 1, 1] = 0
    riemann_arr[0, 1, 0, 0] = 0
    riemann_arr[0, 1, 1, 0] = 1
    riemann_arr[1, 0, 0, 1] = cos(x2)**2
    riemann_ten = riemann_arr.to_tensor((-1, -1, -1, 1))

    assert k_sigma(X, Y, riemann_ten, g, var_list) == 1
    assert k_sigma(X, Y, riemann_ten, g_ten, var_list) == 1
    assert k_sigma(X, Y, riemann_ten, g_ten0, var_list) == 1
