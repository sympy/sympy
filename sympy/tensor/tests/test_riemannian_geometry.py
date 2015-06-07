# -*- coding: utf-8 -*-

from sympy.matrices import Matrix
from sympy.tensor.arraypy import Arraypy, TensorArray
from sympy.tensor.riemannian_geometry import scal_prod, christoffel_1,\
    christoffel_2, covar_der, covar_der_xy, riemann, ricci, scal_curv, k_sigma, \
    kulkarni_nomizu, k_sigma_li, riemann_li, delta, nabla_x, nabla
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
    assert isinstance(christoffel_1(g, var_list), TensorArray)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    assert christoffel_1(g, var_list, 't') == res_ten
    assert isinstance(christoffel_1(g, var_list, 't'), TensorArray)
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
    assert isinstance(christoffel_1(g, var_list), TensorArray)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    assert christoffel_1(g, var_list, 't') == res_ten
    assert isinstance(christoffel_1(g, var_list, 't'), TensorArray)
    assert christoffel_1(g, var_list).type_pq == (0, 3)

    print('test_christoffel_1_gm_a  <=== actual test code')
    assert christoffel_1(g, var_list, 'a') == res_arr
    assert isinstance(christoffel_1(g, var_list, 'a'), Arraypy)


def test_christoffel_2_gtnsr():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]
    g = Arraypy([2, 2, 1]).to_tensor((-1, -1))
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
    res_ten = res_arr.to_tensor((1, -1, -1))

    print('test_christoffel_2_gtnsr_t  <=== actual test code')
    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list), TensorArray)
    assert christoffel_2(g, var_list).type_pq == (1, 2)

    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list, 't'), TensorArray)
    assert christoffel_2(g, var_list).type_pq == (1, 2)

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
    res_ten = res_arr.to_tensor((1, -1, -1))

    print('test_christoffel_2_gm_t  <=== actual test code')
    assert christoffel_2(g, var_list) == res_ten
    assert isinstance(christoffel_2(g, var_list), TensorArray)
    assert christoffel_2(g, var_list).type_pq == (1, 2)

    assert christoffel_2(g, var_list, 't') == res_ten
    assert isinstance(christoffel_2(g, var_list, 't'), TensorArray)
    assert christoffel_2(g, var_list, 't').type_pq == (1, 2)

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
    res_ten = res_arr.to_tensor(( 1, -1))

    print('test_covar_der_t  <=== actual test code')
    assert covar_der(X, g, var_list) == res_ten
    assert isinstance(covar_der(X, g, var_list), TensorArray)
    assert covar_der(X, g, var_list).type_pq == (1, 1)

    assert covar_der(X, g, var_list, 't') == res_ten
    assert isinstance(covar_der(X, g, var_list, 't'), TensorArray)
    assert covar_der(X, g, var_list, 't').type_pq == (1, 1)

    print('test_covar_der_a  <=== actual test code')
    assert covar_der(X, g, var_list, 'a') == res_arr
    assert isinstance(covar_der(X, g, var_list, 'a'), Arraypy)


def test_covar_der_xy():
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

    print('test_covar_der_xy_t  <=== actual test code')
    assert covar_der_xy(X, Y, g, var_list) == res_ten
    assert isinstance(covar_der_xy(X, Y, g, var_list), TensorArray)
    assert covar_der_xy(X, Y, g, var_list).type_pq == (1, 0)

    assert covar_der_xy(X, Y, g, var_list, 't') == res_ten
    assert isinstance(covar_der_xy(X, Y, g, var_list, 't'), TensorArray)
    assert covar_der_xy(X, Y, g, var_list, 't').type_pq == (1, 0)

    print('test_covar_der_xy_a  <=== actual test code')
    assert covar_der_xy(X, Y, g, var_list, 'a') == res_arr
    assert isinstance(covar_der_xy(X, Y, g, var_list, 'a'), Arraypy)


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
    res_ten = res_arr.to_tensor((1, -1, -1, -1))

    print('test_riemann_gtnsr_t  <=== actual test code')
    assert riemann(g, var_list) == res_ten
    assert isinstance(riemann(g, var_list), TensorArray)
    assert riemann(g, var_list).type_pq == (1, 3)

    assert riemann(g, var_list, 't') == res_ten
    assert isinstance(riemann(g, var_list, 't'), TensorArray)
    assert riemann(g, var_list, 't').type_pq == (1, 3)

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
    res_ten = res_arr.to_tensor((1, -1, -1, -1))

    print('test_riemann_gm_t  <=== actual test code')
    assert riemann(g, var_list) == res_ten
    assert isinstance(riemann(g, var_list), TensorArray)
    assert riemann(g, var_list).type_pq == (1, 3)

    assert riemann(g, var_list, 't') == res_ten
    assert isinstance(riemann(g, var_list, 't'), TensorArray)
    assert riemann(g, var_list, 't').type_pq == (1, 3)

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
    riemann_ten = riemann_arr.to_tensor((1, -1, -1, -1))

    res_arr = Arraypy([2, 2, 1])
    res_arr[1, 1] = cos(x2)**2
    res_arr[1, 2] = 0
    res_arr[2, 2] = 0
    res_arr[2, 2] = 1
    res_ten = res_arr.to_tensor((-1, -1))

    print('test_ricci_riemtnsr1_t  <=== actual test code')
    assert ricci(riemann_ten, var_list) == res_ten
    assert isinstance(ricci(riemann_ten, var_list), TensorArray)
    assert ricci(riemann_ten, var_list).type_pq == (0, 2)

    assert ricci(riemann_ten, var_list, 't') == res_ten
    assert isinstance(ricci(riemann_ten, var_list, 't'), TensorArray)
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
    riemann_ten = riemann_arr.to_tensor((1, -1, -1, -1))

    res_arr = Arraypy([2, 2, 0])
    res_arr[0, 0] = cos(x2)**2
    res_arr[0, 1] = 0
    res_arr[1, 1] = 0
    res_arr[1, 1] = 1
    res_ten = res_arr.to_tensor((-1, -1))

    print('test_ricci_riemtnsr0_t  <=== actual test code')
    assert ricci(riemann_ten, var_list) == res_ten
    assert isinstance(ricci(riemann_ten, var_list), TensorArray)
    assert ricci(riemann_ten, var_list).type_pq == (0, 2)

    assert ricci(riemann_ten, var_list, 't') == res_ten
    assert isinstance(ricci(riemann_ten, var_list, 't'), TensorArray)
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


def test_nabla():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    T = Arraypy([2, 2, 0]).to_tensor((1, -1))
    T[0, 0] = x2
    T[0, 1] = -x2
    T[1, 0] = -x1
    T[1, 1] = x1

    ch_2 = Arraypy([3, 2, 0])
    ch_2[0, 0, 0] = 0
    ch_2[0, 0, 1] = sin(x2) * cos(x2)
    ch_2[0, 1, 1] = 0
    ch_2[1, 1, 1] = 0
    ch_2[1, 0, 1] = 0
    ch_2[1, 1, 0] = 0
    ch_2[1, 0, 0] = -sin(x2) * cos(x2)
    ch_2[0, 1, 0] = -sin(x2) * cos(x2)

    res_ten = Arraypy([2, 2, 0]).to_tensor((1, -1, -1))
    res_ten[0, 0, 0] = -x1 * sin(x2) * cos(x2) + x2 * sin(x2) * cos(x2)
    res_ten[0, 0, 1] = 0
    res_ten[0, 1, 1] = x2 * sin(x2) / cos(x2) - 1
    res_ten[1, 1, 1] = 0
    res_ten[1, 0, 1] = -x1 * sin(x2) / cos(x2) - 1
    res_ten[1, 1, 0] = -x1 * sin(x2) / cos(x2) + x2 * sin(x2) / cos(x2)
    res_ten[1, 0, 0] = -x1 * sin(x2) * cos(x2) - x2 * sin(x2) / cos(x2)
    res_ten[0, 1, 0] = x1 * sin(x2) * cos(x2) + x2 * sin(x2) / cos(x2)

    assert nabla(T, ch_2, var_list) == res_ten
    assert isinstance(nabla(T, ch_2, var_list), TensorArray)


def test_nabla_x():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    X = [x1 * x2**3, x1 - cos(x2)]
    
    T = Arraypy([2, 2, 0]).to_tensor((-1, -1))
    T[0, 0] = x2
    T[0, 1] = -x2
    T[1, 0] = -x1
    T[1, 1] = x1

    ch_2 = Arraypy([3, 2, 0])
    ch_2[0, 0, 0] = 0
    ch_2[0, 0, 1] = sin(x2) * cos(x2)
    ch_2[0, 1, 1] = 0
    ch_2[1, 1, 1] = 0
    ch_2[1, 0, 1] = 0
    ch_2[1, 1, 0] = 0
    ch_2[1, 0, 0] = -sin(x2) * cos(x2)
    ch_2[0, 1, 0] = -sin(x2) * cos(x2)

    res_ten = Arraypy((2, 2)).to_tensor((-1, -1))
    res_ten[0, 0] = x1 * x2**3 * (x1 * sin(x2) * cos(x2) + x2 *
                                  sin(x2) * cos(x2)) + 2 * x2 * (x1 - cos(x2)) * sin(x2) / cos(x2)
    res_ten[0, 1] = x1 * x2**3 * (-x1 * sin(x2) * cos(x2) + x2 * sin(
        x2) / cos(x2)) + (x1 - cos(x2)) * (-x2 * sin(x2) / cos(x2) - 1)
    res_ten[1, 0] = x1 * x2**3 * (-x1 * sin(x2) * cos(x2) + x2 * sin(
        x2) / cos(x2)) + (x1 - cos(x2)) * (-x1 * sin(x2) / cos(x2) - 1)
    res_ten[1, 1] = x1 * x2**3 * \
        (-x1 * sin(x2) / cos(x2) - x2 * sin(x2) / cos(x2))

    assert nabla_x(T, ch_2, X, var_list) == res_ten
    assert isinstance(nabla_x(T, ch_2, X, var_list), TensorArray)


def test_delta():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    T = Arraypy([2, 2, 0]).to_tensor((1, -1))
    T[0, 0] = x2
    T[0, 1] = -x2
    T[1, 0] = -x1
    T[1, 1] = x1

    g = Arraypy((2, 2)).to_tensor((-1, -1))
    g[0, 0] = cos(x2)**2
    g[0, 1] = 0
    g[1, 0] = 0
    g[1, 1] = 1

    ch_2 = Arraypy([3, 2, 0])
    ch_2[0, 0, 0] = 0
    ch_2[0, 0, 1] = sin(x2) * cos(x2)
    ch_2[0, 1, 1] = 0
    ch_2[1, 1, 1] = 0
    ch_2[1, 0, 1] = 0
    ch_2[1, 1, 0] = 0
    ch_2[1, 0, 0] = -sin(x2) * cos(x2)
    ch_2[0, 1, 0] = -sin(x2) * cos(x2)

    res_ten = Arraypy((2)).to_tensor((-1))
    res_ten[0] = x1 * sin(x2) * cos(x2) + 1
    res_ten[1] = 0

    assert delta(T, g, ch_2, var_list) == res_ten
    assert isinstance(delta(T, g, ch_2, var_list), TensorArray)


def test_riemann_Li():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    C = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    C[0, 0, 0] = 0
    C[0, 0, 1] = sin(x2) * cos(x2)
    C[0, 1, 1] = 0
    C[1, 1, 1] = 0
    C[1, 0, 1] = 0
    C[1, 1, 0] = 0
    C[1, 0, 0] = -sin(x2) * cos(x2)
    C[0, 1, 0] = -sin(x2) * cos(x2)

    g = Arraypy((2, 2)).to_tensor((-1, -1))
    g[0, 0] = cos(x2)**2
    g[0, 1] = 0
    g[1, 0] = 0
    g[1, 1] = 1

    res_arr = Arraypy([4, 2, 0])
    res_arr[0, 0, 0, 0] = -0.25 * sin(x2)**2 * cos(x2)**2
    res_arr[0, 0, 0, 1] = 0
    res_arr[0, 0, 1, 1] = 0
    res_arr[0, 1, 1, 1] = 0
    res_arr[1, 1, 1, 1] = 0
    res_arr[1, 1, 1, 0] = 0
    res_arr[1, 1, 0, 0] = 0
    res_arr[1, 0, 0, 0] = 0
    res_arr[1, 0, 1, 0] = 0
    res_arr[0, 1, 0, 1] = 0
    res_arr[1, 1, 0, 1] = 0
    res_arr[0, 0, 1, 0] = 0
    res_arr[1, 0, 1, 1] = 0
    res_arr[0, 1, 0, 0] = 0
    res_arr[0, 1, 1, 0] = 0
    res_arr[1, 0, 0, 1] = 0
    res_ten = res_arr.to_tensor((1, -1, -1, -1))
    print('test_riemann_li_t  <=== actual test code')
    assert riemann_li(C, g, var_list) == res_ten
    assert isinstance(riemann_li(C, g, var_list), TensorArray)
    assert riemann_li(C, g, var_list).type_pq == (1, 3)

    assert riemann_li(C, g, var_list, 't') == res_ten
    assert isinstance(riemann_li(C, g, var_list, 't'), TensorArray)
    assert riemann_li(C, g, var_list, 't').type_pq == (1, 3)

    print('test_riemann_li_a  <=== actual test code')
    assert riemann_li(C, g, var_list, 'a') == res_arr
    assert isinstance(riemann_li(C, g, var_list, 'a'), Arraypy)


def test_kulkarni_nomizu():
    x1, x2 = symbols('x1, x2')
    var_list = [x1, x2]

    h = Arraypy((2, 2)).to_tensor((-1, -1))
    h[0, 0] = x1
    h[0, 1] = 0
    h[1, 0] = 0
    h[1, 1] = x2

    k = Arraypy((2, 2)).to_tensor((-1, -1))
    k[0, 0] = x2
    k[0, 1] = 0
    k[1, 0] = 0
    k[1, 1] = x1

    res_arr = Arraypy([4, 2, 0])
    res_arr[0, 0, 0, 0] = 0
    res_arr[0, 0, 0, 1] = 0
    res_arr[0, 0, 1, 1] = 0
    res_arr[0, 1, 1, 1] = 0
    res_arr[1, 1, 1, 1] = 0
    res_arr[1, 1, 1, 0] = 0
    res_arr[1, 1, 0, 0] = 0
    res_arr[1, 0, 0, 0] = 0
    res_arr[1, 0, 1, 0] = x1**2 + x2**2
    res_arr[0, 1, 0, 1] = x1**2 + x2**2
    res_arr[1, 1, 0, 1] = 0
    res_arr[0, 0, 1, 0] = 0
    res_arr[1, 0, 1, 1] = 0
    res_arr[0, 1, 0, 0] = 0
    res_arr[0, 1, 1, 0] = -x1**2 - x2**2
    res_arr[1, 0, 0, 1] = -x1**2 - x2**2
    res_ten = res_arr.to_tensor((-1, -1, -1, -1))
    
    print('test_kulkarni_nomizu_t  <=== actual test code')
    assert kulkarni_nomizu(h, k, var_list) == res_ten
    assert isinstance(kulkarni_nomizu(h, k, var_list), TensorArray)
    assert kulkarni_nomizu(h, k, var_list).type_pq == (0, 4)

    assert kulkarni_nomizu(h, k, var_list, 't') == res_ten
    assert isinstance(kulkarni_nomizu(h, k, var_list, 't'), TensorArray)
    assert kulkarni_nomizu(h, k, var_list, 't').type_pq == (0, 4)

    print('test_kulkarni_nomizu_a  <=== actual test code')
    assert kulkarni_nomizu(h, k, var_list, 'a') == res_arr
    assert isinstance(kulkarni_nomizu(h, k, var_list, 'a'), Arraypy)
