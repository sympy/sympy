"""
N-dim array module for SymPy.

Four classes are provided to handle N-dim arrays, given by the combinations
dense/sparse (i.e. whether to store all elements or only the non-zero ones in
memory) and mutable/immutable (immutable classes are SymPy objects, but cannot
change after they have been created).

Examples
========

The following examples show the usage of ``ImmutableDenseNDimArray``, the other
classes are analogous. For mutable classes it is also possible to change
element values after the object has been constructed.

Array construction can detect the shape of nested lists and tuples:

>>> from sympy.tensor.array import ImmutableDenseNDimArray
>>> a1 = ImmutableDenseNDimArray([[1, 2], [3, 4], [5, 6]])
>>> a1
[[1, 2], [3, 4], [5, 6]]
>>> a1.shape
(3, 2)
>>> a1.rank()
2
>>> from sympy.abc import x, y, z
>>> a2 = ImmutableDenseNDimArray([[[x, y], [z, x*z]], [[1, x*y], [1/x, x/y]]])
>>> a2
[[[x, y], [z, x*z]], [[1, x*y], [1/x, x/y]]]
>>> a2.shape
(2, 2, 2)
>>> a2.rank()
3

Otherwise one could pass a 1-dim array followed by a shape tuple:

>>> m1 = ImmutableDenseNDimArray(range(12), (3, 4))
>>> m1
[[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]
>>> m2 = ImmutableDenseNDimArray(range(12), (3, 2, 2))
>>> m2
[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]]
>>> m2[1,1,1]
7
>>> m2.reshape(4, 3)
[[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]

Slice support:

>>> m2[:, 1, 1]
[3, 7, 11]

Elementwise derivative:

>>> from sympy.abc import x, y, z
>>> m3 = ImmutableDenseNDimArray([x**3, x*y, z])
>>> m3.diff(x)
[3*x**2, y, 0]
>>> m3.diff(z)
[0, 0, 1]

Multiplication with other SymPy expressions is applied elementwisely:

>>> (1+x)*m3
[x**3*(x + 1), x*y*(x + 1), z*(x + 1)]

To apply a function to each element of the N-dim array, use ``applyfunc``:

>>> m3.applyfunc(lambda x: x/2)
[x**3/2, x*y/2, z/2]

N-dim arrays can be converted to nested lists by the ``tolist()`` method:

>>> m2.tolist()
[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]]
>>> isinstance(m2.tolist(), list)
True

If the rank is 2, it is possible to convert them to matrices with ``tomatrix()``:

>>> m1.tomatrix()
Matrix([
[0, 1,  2,  3],
[4, 5,  6,  7],
[8, 9, 10, 11]])

"""

from .dense_ndim_array import MutableDenseNDimArray, ImmutableDenseNDimArray
from .sparse_ndim_array import MutableSparseNDimArray, ImmutableSparseNDimArray
from .arrayop import tensorproduct, tensorcontraction

Array = ImmutableDenseNDimArray
NDimArray = ImmutableDenseNDimArray
DenseNDimArray = ImmutableDenseNDimArray
SparseNDimArray = ImmutableSparseNDimArray
