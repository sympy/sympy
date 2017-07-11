"""A module to manipulate symbolic objects with indices including tensors

"""
from .array import Array, DenseNDimArray, ImmutableDenseNDimArray, \
    ImmutableSparseNDimArray, MutableDenseNDimArray, MutableSparseNDimArray, \
    NDimArray, SparseNDimArray, derive_by_array, permutedims, \
    tensorcontraction, tensorproduct
from .index_methods import get_contraction_structure, get_indices
from .indexed import Idx, Indexed, IndexedBase
