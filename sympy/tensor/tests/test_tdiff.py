from sympy.core.function import diff, Function
from sympy.core.symbol import symbols
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import sin
from sympy.tensor.tdiff import tdiff
from sympy.tensor.tensor import tensor_indices
from sympy.tensor.vtensor import VTensorIndexType, vtensorhead


def test_tdiff_generic():
    L = VTensorIndexType('L', [1, -1])
    i0, i1, i2 = tensor_indices('i0:3', L)
    x, y, z, t = symbols('x y z t')

    t22_values = [[x * y, sin(x) * y], [exp(x * y), log(x) / y]]
    A = vtensorhead('A', [L] * 2, [[1], [1]], values=t22_values)

    XY = vtensorhead('XY', [L], [[1]], values=[x, y])

    tdiff_result_rank_3 = tdiff(A(i1, i2), XY(i0))

    assert tdiff_result_rank_3.rank == 3

    # check that all derivatives in the 3x3 tensor are correct:
    for ki0, dvari in enumerate([x, y]):
        for ki1, el0 in enumerate(t22_values):
            for ki2, el1 in enumerate(el0):
                assert diff(el1, dvari) == tdiff_result_rank_3[ki0, ki1, ki2]

    # test contraction of indices
    L4 = VTensorIndexType('L4', [1, -1, -1, -1])
    mu0, mu1, mu2 = tensor_indices('mu0:3', L4)
    B = vtensorhead('B', [L4], [[1]], values=[x * y * z * t, x / y, x * y * log(z), x ** 2 * sin(y) * exp(-t)])
    D4 = vtensorhead('D4', [L4], [[1]], values=[x, y, z, t])
    scalar_diff_result = tdiff(B(-mu0), D4(mu0))

    # Assert diff contraction result
    assert scalar_diff_result == t * y * z + x ** 2 * exp(-t) * sin(y) - x * y / z + x / y ** 2

    scalar_diff_result2 = tdiff(B(mu0), D4(-mu0))
    assert scalar_diff_result == scalar_diff_result2

    # test contractions with many indices (the differential index contracted, the result is still a tensor expression )
    tdiff_result_rank_1 = tdiff(A(i1, i2), XY(-i2))

    for i in range(2):
        assert tdiff_result_rank_1[i] == t22_values[i][0].diff(x) - t22_values[i][1].diff(y)

    # TODO: test contractions on a function of tensors...
#     f = Function('f')
#
#     func_diff_result = tdiff(f(B(mu0)), D4(-mu0))
#
#     print func_diff_result
