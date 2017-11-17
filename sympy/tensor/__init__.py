"""A module to manipulate symbolic objects with indices including tensors

"""

__all__ = []

from .indexed import IndexedBase, Idx, Indexed
__all__ += ["IndexedBase", "Idx", "Indexed"]

from .index_methods import get_contraction_structure, get_indices
__all__ += ["get_contraction_structure", "get_indices"]

from .array import (
    MutableDenseNDimArray, ImmutableDenseNDimArray,
    MutableSparseNDimArray, ImmutableSparseNDimArray,
    NDimArray, tensorproduct, tensorcontraction, derive_by_array,
    permutedims, Array, DenseNDimArray, SparseNDimArray
)
__all__ += [
    "MutableDenseNDimArray", "ImmutableDenseNDimArray",
    "MutableSparseNDimArray", "ImmutableSparseNDimArray",
    "NDimArray", "tensorproduct", "tensorcontraction", "derive_by_array",
    "permutedims", "Array", "DenseNDimArray", "SparseNDimArray"
]
