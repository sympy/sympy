import itertools

import collections

from sympy import S, Tuple, MatrixBase
from sympy import S, Tuple, diff, MatrixBase

from sympy.tensor.array import ImmutableDenseNDimArray
from sympy.tensor.array.ndim_array import NDimArray


def _arrayfy(a):
    if isinstance(a, NDimArray):
        return a
    if isinstance(a, (MatrixBase, list, tuple, Tuple)):
        return ImmutableDenseNDimArray(a)
    return a


def tensorproduct(*args):
    """
    Tensor product among scalars or array-like objects.

    Examples
    ========

    >>> from sympy.tensor.array import tensorproduct, Array
    >>> from sympy.abc import x, y, z, t
    >>> A = Array([[1, 2], [3, 4]])
    >>> B = Array([x, y])
    >>> tensorproduct(A, B)
    [[[x, y], [2*x, 2*y]], [[3*x, 3*y], [4*x, 4*y]]]
    >>> tensorproduct(A, x)
    [[x, 2*x], [3*x, 4*x]]
    >>> tensorproduct(A, B, B)
    [[[[x**2, x*y], [x*y, y**2]], [[2*x**2, 2*x*y], [2*x*y, 2*y**2]]], [[[3*x**2, 3*x*y], [3*x*y, 3*y**2]], [[4*x**2, 4*x*y], [4*x*y, 4*y**2]]]]

    Applying this function on two matrices will result in a rank 4 array.

    >>> from sympy import Matrix, eye
    >>> m = Matrix([[x, y], [z, t]])
    >>> p = tensorproduct(eye(3), m)
    >>> p
    [[[[x, y], [z, t]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[x, y], [z, t]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [0, 0]], [[x, y], [z, t]]]]
    """
    if len(args) == 0:
        return S.One
    if len(args) == 1:
        return _arrayfy(args[0])
    if len(args) > 2:
        return tensorproduct(tensorproduct(args[0], args[1]), *args[2:])

    # length of args is 2:
    a, b = map(_arrayfy, args)

    if not isinstance(a, NDimArray) or not isinstance(b, NDimArray):
        return a*b

    al = list(a)
    bl = list(b)

    product_list = [i*j for i in al for j in bl]
    return ImmutableDenseNDimArray(product_list, a.shape + b.shape)


def tensorcontraction(array, *contraction_axes):
    """
    Contraction of an array-like object on the specified axes.

    Examples
    ========

    >>> from sympy.tensor.array import Array, tensorcontraction
    >>> from sympy import Matrix, eye
    >>> tensorcontraction(eye(3), (0, 1))
    3
    >>> A = Array(range(18), (3, 2, 3))
    >>> A
    [[[0, 1, 2], [3, 4, 5]], [[6, 7, 8], [9, 10, 11]], [[12, 13, 14], [15, 16, 17]]]
    >>> tensorcontraction(A, (0, 2))
    [21, 30]

    Matrix multiplication may be emulated with a proper combination of
    ``tensorcontraction`` and ``tensorproduct``

    >>> from sympy.tensor.array import tensorproduct
    >>> from sympy.abc import a,b,c,d,e,f,g,h
    >>> m1 = Matrix([[a, b], [c, d]])
    >>> m2 = Matrix([[e, f], [g, h]])
    >>> p = tensorproduct(m1, m2)
    >>> p
    [[[[a*e, a*f], [a*g, a*h]], [[b*e, b*f], [b*g, b*h]]], [[[c*e, c*f], [c*g, c*h]], [[d*e, d*f], [d*g, d*h]]]]
    >>> tensorcontraction(p, (1, 2))
    [[a*e + b*g, a*f + b*h], [c*e + d*g, c*f + d*h]]
    >>> m1*m2
    Matrix([
    [a*e + b*g, a*f + b*h],
    [c*e + d*g, c*f + d*h]])
    """
    array = _arrayfy(array)

    # Verify contraction_axes:
    taken_dims = set([])
    for axes_group in contraction_axes:
        if not isinstance(axes_group, collections.Iterable):
            raise ValueError("collections of contraction axes expected")

        dim = array.shape[axes_group[0]]

        for d in axes_group:
            if d in taken_dims:
                raise ValueError("dimension specified more than once")
            if dim != array.shape[d]:
                raise ValueError("cannot contract between axes of different dimension")
            taken_dims.add(d)

    rank = array.rank()

    remaining_shape = [dim for i, dim in enumerate(array.shape) if i not in taken_dims]
    cum_shape = [0]*rank
    _cumul = 1
    for i in range(rank):
        cum_shape[rank - i - 1] = _cumul
        _cumul *= int(array.shape[rank - i - 1])

    remaining_indices = [[cum_shape[i]*j for j in range(array.shape[i])]
                         for i in range(rank) if i not in taken_dims]

    contracted_array = []
    for icontrib in itertools.product(*remaining_indices):
        i = sum(icontrib)
        isum = S.Zero
        for axes_group in contraction_axes:
            for js in range(array.shape[axes_group[0]]):
                isum += array[i + sum([cum_shape[ig]*js for ig in axes_group])]
        contracted_array.append(isum)

    if len(remaining_indices) == 0:
        assert len(contracted_array) == 1
        return contracted_array[0]

    return type(array)(contracted_array, remaining_shape)


def derive_by_array(expr, dx):
    r"""
    Derivative by arrays. Supports both arrays and scalars.

    Given the array `A_{i_1, \ldots, i_N}` and the array `X_{j_1, \ldots, j_M}`
    this function will return a new array `B` defined by

    `B_{j_1,\ldots,j_M,i_1,\ldots,i_N} := \frac{\partial A_{i_1,\ldots,i_N}}{\partial X_{j_1,\ldots,j_M}}`

    Examples
    ========

    >>> from sympy.tensor.array import derive_by_array
    >>> from sympy.abc import x, y, z, t
    >>> from sympy import cos
    >>> derive_by_array(cos(x*t), x)
    -t*sin(t*x)
    >>> derive_by_array(cos(x*t), [x, y, z, t])
    [-t*sin(t*x), 0, 0, -x*sin(t*x)]
    >>> derive_by_array([x, y**2*z], [[x, y], [z, t]])
    [[[1, 0], [0, 2*y*z]], [[0, y**2], [0, 0]]]

    """
    array_types = (collections.Iterable, MatrixBase, NDimArray)

    if isinstance(dx, array_types):
        dx = ImmutableDenseNDimArray(dx)
        for i in dx:
            if not i._diff_wrt:
                raise ValueError("cannot derive by this array")

    if isinstance(expr, array_types):
        expr = ImmutableDenseNDimArray(expr)
        if isinstance(dx, array_types):
            new_array = [[y.diff(x) for y in expr] for x in dx]
            return type(expr)(new_array, dx.shape + expr.shape)
        else:
            return expr.diff(dx)
    else:
        if isinstance(dx, array_types):
            return ImmutableDenseNDimArray([expr.diff(i) for i in dx], dx.shape)
        else:
            return diff(expr, dx)
